#!/usr/bin/env python3
"""Capture the frozen 20-block study or reconstruct its complete archive offline."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import io
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shutil
import statistics
import subprocess
import tarfile
import tempfile
import time

import immune_2d_measurement_report as measurement
import immune_measurement_report as common
import immune_seed_blocks as seeds
from immune_measurement_report import require, sha

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "analysis/immune-2d-replication"
OUT_JSON = ROOT / "analysis/immune-2d-replication.json"
OUT_MD = ROOT / "analysis/immune-2d-replication.md"
PLAN = "scripts/immune_2d_replication_plan.json"
PROTOCOL = "docs/IMMUNE_2D_REPLICATION_PROTOCOL.md"
CONTRACT = "scripts/immune_2d_measurement_config_v1.json"
CANONICAL_MANIFEST = "analysis/immune-2d-measurements/manifest.json"
ARMS = ("Control", "RSL3", "SDT")
ENDPOINT = "activation_sum_per_initial_tumor_cell"
BLOCKS = tuple(range(1, 21))
# Names of all 86 tracked Rust/Cargo inputs at the v1 source freeze, sorted
# and joined with a final newline. Pin the inventory independently of the
# archive's own manifest; contents remain governed by its per-file hashes.
RUST_INVENTORY_SHA256 = "3c01c856a56fe537554b9fc4d93f1a6f4c56777929bec659edcecd3e78069216"
REQUIRED_SOURCES = {
    PLAN, PROTOCOL, CONTRACT, CANONICAL_MANIFEST, measurement.GOLDEN, measurement.HISTORICAL,
    measurement.PROTOCOL, "scripts/immune_2d_measurement_report.py", "scripts/immune_measurement_report.py",
    "scripts/immune_2d_replication_report.py", "scripts/immune_seed_blocks.py",
    "tests/test_immune_seed_blocks.py", "tests/test_immune_2d_replication_report.py",
    "tests/test_immune_2d_replication_capture.py",
}
FIXED_PLAN = {
    "schema_version": 1, "study": "immune-2d-activation-replication",
    "simulator": "sim-tme", "dimension": 2, "block_ids": list(BLOCKS),
    "base_seed": 42, "block_seed_stride": 2**32, "arms": list(ARMS),
    "primary_endpoint": ENDPOINT, "contrast": ["SDT", "RSL3"],
    "bootstrap": {"unit": "whole_three_arm_block",
                  "statistic": "arithmetic_mean_of_paired_differences",
                  "resamples": 10000, "seed": 20260922,
                  "quantiles": [0.025, 0.975], "interpolation": "linear"},
}


def validate_plan(plan: dict) -> None:
    require(json.dumps(plan, sort_keys=True, allow_nan=False)
            == json.dumps(FIXED_PLAN, sort_keys=True), "frozen replication plan changed")


def block_name(block: int) -> str:
    seeds.block_seed(block)
    return f"block-{block:02d}.json.gz"


ARTIFACTS = {block_name(b) for b in BLOCKS} | {
    "baseline-summary.json", "canonical-observations.json.gz", "sources.tar.gz", "capture-logs.tar.gz"}


def validate_block(data: dict, block: int, contract: dict) -> dict:
    root = seeds.block_seed(block)
    require(type(data["replicate_block"]) is int and data["replicate_block"] == block,
            "replicate block identity")
    require(type(data["schema_version"]) is int and data["schema_version"] == 2
            and data["simulator"] == "sim-tme"
            and type(data["dimension"]) is int and data["dimension"] == 2,
            "replicate observation schema")
    cfg = data["config"]
    require(type(contract["schema_version"]) is int and contract["schema_version"] == 2, "frozen configuration schema")
    require(json.dumps(cfg, sort_keys=True, allow_nan=False)
            == json.dumps(dict(contract["config"], seed=root), sort_keys=True, allow_nan=False),
            "complete replicate configuration changed")
    require(type(cfg["seed"]) is int, "integer root seed")
    require([r["condition_name"] for r in data["conditions"]] == list(measurement.NAMES),
            "ordered replicate arms")
    n_tumor = sum(measurement.in_tumor_circle(i, cfg)
                  for i in range(cfg["grid_rows"] * cfg["grid_cols"]))
    require(n_tumor > 0, "nonempty initial tumor denominator")
    summaries = {}
    for arm_index, row in enumerate(data["conditions"]):
        arm, result = ARMS[arm_index], row["result"]
        require(type(row["seed"]) is int
                and row["seed"] == root + arm_index * cfg["rng"]["treatment_seed_stride"],
                "replicate treatment seed")
        require(result["treatment"] == arm and result["o2_condition"] == "gradient_120um"
                and result["o2_lambda_um"] == 120. and result["immune_mode"] == "immune_on"
                and result["stromal_mode"] == "off" and "ph_mode" not in result,
                "replicate condition identity")
        require(result["total_tumor"] == n_tumor, "replicate initial tumor census")
        summary = measurement.reconcile_condition(row, cfg)
        for population in ("ferroptotic_events", "immune_kill_events", "eligible_cells"):
            require(all(measurement.in_tumor_circle(c["cell_index"], cfg)
                        for c in row["measurements"][population]), "observed cell outside tumor circle")
        total_activation = math.fsum(c["activation_sum"] for c in row["measurements"]["eligible_cells"])
        summary.update(activation_sum=total_activation,
                       activation_sum_per_initial_tumor_cell=total_activation / n_tumor)
        summaries[arm] = summary
    return {"block": block, "root_seed": root, "arms": summaries,
            "sdt_rsl3_immune_kill_ratio": common.ratio(summaries["SDT"]["immune_kills"],
                                                      summaries["RSL3"]["immune_kills"])}


def quantile(values: list[float], p: float) -> float:
    ordered = sorted(values)
    require(bool(ordered) and 0 <= p <= 1, "nonempty quantile input")
    position = (len(ordered) - 1) * p
    lo, hi = math.floor(position), math.ceil(position)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (position - lo)


def distribution(values: list[tuple[int, float | None]]) -> dict:
    defined = [v for _, v in values if v is not None]
    return {"n_defined": len(defined), "undefined_blocks": [b for b, v in values if v is None],
            "median": statistics.median(defined) if defined else None,
            "minimum": min(defined) if defined else None,
            "maximum": max(defined) if defined else None}


def summarize(blocks: list[dict], plan: dict) -> dict:
    validate_plan(plan)
    require([b["block"] for b in blocks] == list(BLOCKS), "complete ordered block coverage")
    for block in blocks:
        require(block["root_seed"] == seeds.block_seed(block["block"]), "summary block seed")
    differences = [b["arms"]["SDT"][ENDPOINT] - b["arms"]["RSL3"][ENDPOINT] for b in blocks]
    require(all(math.isfinite(d) for d in differences), "finite primary contrast")
    rng = random.Random(plan["bootstrap"]["seed"])
    bootstrap = [math.fsum(differences[rng.randrange(len(blocks))] for _ in blocks) / len(blocks)
                 for _ in range(plan["bootstrap"]["resamples"])]
    lower, upper = [quantile(bootstrap, p) for p in plan["bootstrap"]["quantiles"]]
    fields = sorted(k for k in blocks[0]["arms"]["Control"] if k not in ("condition", "condition_seed"))
    return {
        "primary": {"endpoint": ENDPOINT, "contrast": "SDT minus RSL3",
                    "units": "accumulated activation steps per initial tumor cell",
                    "mean_paired_difference": math.fsum(differences) / len(blocks),
                    "percentile_interval": [lower, upper], "paired_differences": differences,
                    "interpretation": "positive" if lower > 0 else "negative" if upper < 0 else "unresolved",
                    "bootstrap": plan["bootstrap"]},
        "arm_distributions": {arm: {key: distribution([(b["block"], b["arms"][arm][key])
                                                       for b in blocks]) for key in fields} for arm in ARMS},
        "threshold_positive_blocks": {
            arm: {key: [b["block"] for b in blocks if b["arms"][arm][key] > 0]
                  for key in measurement.COUNTS} for arm in ARMS},
        "sdt_rsl3_immune_kill_ratio": distribution([
            (b["block"], b["sdt_rsl3_immune_kill_ratio"]) for b in blocks]),
    }


def pack_files(contents: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as bundle:
        for name, blob in sorted(contents.items()):
            entry = tarfile.TarInfo(name)
            entry.size, entry.mode = len(blob), 0o644
            bundle.addfile(entry, io.BytesIO(blob))
    return gzip.compress(buffer.getvalue(), mtime=0)


def unpack_sources(blob: bytes, hashes: dict) -> dict[str, bytes]:
    require(REQUIRED_SOURCES <= set(hashes), "required frozen source inventory")
    rust_paths = sorted(p for p in hashes if p.startswith("simulations/") and p.endswith((".rs", ".toml", ".lock")))
    require(sha(("\n".join(rust_paths) + "\n").encode()) == RUST_INVENTORY_SHA256,
            "complete frozen Rust source inventory")
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as bundle:
        members = bundle.getmembers()
        require(len(members) == len(hashes) and {m.name for m in members} == set(hashes),
                "frozen source inventory")
        contents = {}
        for member in members:
            require(member.isfile() and not Path(member.name).is_absolute()
                    and ".." not in Path(member.name).parts, "invalid source archive entry")
            contents[member.name] = bundle.extractfile(member).read()
            require(sha(contents[member.name]) == hashes[member.name], f"source hash: {member.name}")
    return contents


def load_archive(path: Path) -> dict:
    require(path.is_dir() and not path.is_symlink(), "regular archive directory")
    require({p.name for p in path.iterdir()} == ARTIFACTS | {"manifest.json"}, "archive file inventory")
    require(all(p.is_file() and not p.is_symlink() for p in path.iterdir()), "regular archive files")
    manifest = json.loads((path / "manifest.json").read_text())
    require(type(manifest["schema_version"]) is int and manifest["schema_version"] == 1
            and manifest["study"] == FIXED_PLAN["study"]
            and manifest["execution"] == "serial" and manifest["complete"] is True, "complete archive schema")
    require(isinstance(manifest["source_commit"], str)
            and re.fullmatch(r"[0-9a-f]{40}", manifest["source_commit"])
            and isinstance(manifest["binary_sha256"], str)
            and re.fullmatch(r"[0-9a-f]{64}", manifest["binary_sha256"]), "capture source/binary provenance")
    require(set(manifest["artifacts"]) == ARTIFACTS, "manifest artifact inventory")
    for name in sorted(ARTIFACTS):
        require(sha((path / name).read_bytes()) == manifest["artifacts"][name], f"artifact hash: {name}")
    frozen = unpack_sources((path / "sources.tar.gz").read_bytes(), manifest["sources"])
    plan = json.loads(frozen[PLAN])
    validate_plan(plan)
    contract = json.loads(frozen[CONTRACT])
    baseline = (path / "baseline-summary.json").read_bytes()
    require(sha(baseline) == measurement.golden_sha(frozen[measurement.GOLDEN].decode()), "full baseline parity")
    require(manifest["baseline_reference"] == measurement.baseline_reference(frozen[measurement.GOLDEN].decode()),
            "baseline provenance")
    canonical = (path / "canonical-observations.json.gz").read_bytes()
    old_manifest = json.loads(frozen[CANONICAL_MANIFEST])
    require(sha(canonical) == old_manifest["artifacts"]["observations.json.gz"], "canonical observation parity")
    measurement.validate_observations(json.loads(gzip.decompress(canonical)), json.loads(baseline),
                                      contract=contract, historical=json.loads(frozen[measurement.HISTORICAL]))
    require([r["block"] for r in manifest["runs"]] == list(BLOCKS), "capture run coverage")
    blocks = []
    for block, run in zip(BLOCKS, manifest["runs"]):
        require(type(run["block"]) is int and type(run["seed"]) is int
                and run["seed"] == seeds.block_seed(block) and run["artifact"] == block_name(block), "capture run identity")
        common.finite(run["elapsed_seconds"], "capture duration")
        data = json.loads(gzip.decompress((path / block_name(block)).read_bytes()))
        blocks.append(validate_block(data, block, contract))
        del data
    return {"schema_version": 1, "study": plan["study"], "plan": plan,
            "provenance": manifest, "seed_audit": seeds.audit(contract["config"]),
            "blocks": blocks, **summarize(blocks, plan)}


def source_paths() -> list[str]:
    return sorted(set(measurement.source_paths()) | {
        PLAN, PROTOCOL, CANONICAL_MANIFEST, "scripts/immune_seed_blocks.py",
        "scripts/immune_2d_replication_report.py", "tests/test_immune_seed_blocks.py",
        "tests/test_immune_2d_replication_report.py", "tests/test_immune_2d_replication_capture.py",
        # Both replication test modules import fixtures from this older suite.
        # The original frozen archive is preserved; its recovery recipe is in
        # docs/IMMUNE_2D_REPLICATION_RESULTS.md.
        "tests/test_immune_2d_measurement_report.py"})


def capture(destination: Path) -> None:
    require(not destination.exists() and not destination.is_symlink(), "immutable capture destination already exists")
    require(not any(k.startswith("FERRO_") for k in os.environ), "unset all FERRO_* overrides")
    require(not subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip(),
            "commit the protocol, plan, implementation and tests before capture")
    sources = {p: (ROOT / p).read_bytes() for p in source_paths()}
    validate_plan(json.loads(sources[PLAN]))
    contract = json.loads(sources[CONTRACT])
    seeds.audit(contract["config"])
    canonical = (measurement.ARCHIVE / "observations.json.gz").read_bytes()
    require(sha(canonical) == json.loads(sources[CANONICAL_MANIFEST])["artifacts"]["observations.json.gz"],
            "historical canonical archive hash")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    rustc = subprocess.check_output(["rustup", "run", common.TOOLCHAIN, "rustc", "--version"], text=True).strip()
    destination.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".immune-replication-capture-", dir=destination.parent))
    print(f"Capture evidence (retained on failure): {work}", flush=True)
    staged = work / "archive"
    staged.mkdir()
    try:
        binary = common.build_binary(["rustup", "run", common.TOOLCHAIN, "cargo"], ROOT / "simulations", "sim-tme")
        private = work / "sim-tme"
        shutil.copy2(binary, private)
        binary_sha = sha(private.read_bytes())
        logs = {}

        def execute(name: str, args: list[str], output: str) -> tuple[bytes, float]:
            cwd = work / name
            cwd.mkdir()
            start = time.monotonic()
            log = work / f"{name}.log"
            with log.open("wb") as stream:
                subprocess.run([str(private), *args], cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, check=True)
            logs[log.name] = log.read_bytes()
            raw = (cwd / "output/tme" / output).read_bytes()
            elapsed = time.monotonic() - start
            print(f"Completed {name} ({elapsed:.1f}s)", flush=True)
            return raw, elapsed

        baseline, _ = execute("baseline", [], "tme_summary.json")
        require(sha(baseline) == measurement.golden_sha(sources[measurement.GOLDEN].decode()), "full baseline parity")
        raw, _ = execute("canonical", ["--immune-measurements"], "immune_measurements.json")
        require(raw == gzip.decompress(canonical), "canonical observer byte parity")
        measurement.validate_observations(json.loads(raw), json.loads(baseline), contract=contract,
                                          historical=json.loads(sources[measurement.HISTORICAL]))
        (staged / "baseline-summary.json").write_bytes(baseline)
        (staged / "canonical-observations.json.gz").write_bytes(canonical)
        del raw
        runs = []
        for block in BLOCKS:
            raw, elapsed = execute(f"block-{block:02d}", ["--immune-replicate", str(block)], "immune_replicate.json")
            # Retain raw/compressed observations even if their validation fails.
            (staged / block_name(block)).write_bytes(gzip.compress(raw, mtime=0))
            validate_block(json.loads(raw), block, contract)
            runs.append({"block": block, "seed": seeds.block_seed(block),
                         "artifact": block_name(block), "elapsed_seconds": elapsed})
            del raw
        require(sources == {p: (ROOT / p).read_bytes() for p in source_paths()}, "frozen sources changed during capture")
        require(commit == subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "source commit changed during capture")
        (staged / "sources.tar.gz").write_bytes(pack_files(sources))
        (staged / "capture-logs.tar.gz").write_bytes(pack_files(logs))
        manifest = {"schema_version": 1, "study": FIXED_PLAN["study"], "complete": True,
                    "source_commit": commit, "captured_at_utc": datetime.now(timezone.utc).isoformat(),
                    "platform": platform.platform(), "rustc": rustc, "execution": "serial",
                    "python": {"implementation": platform.python_implementation(), "version": platform.python_version()},
                    "binary_sha256": binary_sha, "runs": runs,
                    "baseline_reference": measurement.baseline_reference(sources[measurement.GOLDEN].decode()),
                    "sources": {p: sha(blob) for p, blob in sources.items()},
                    "artifacts": {p: sha((staged / p).read_bytes()) for p in sorted(ARTIFACTS)}}
        (staged / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        render(load_archive(staged))
        require(not destination.exists(), "capture destination appeared during run")
        staged.rename(destination)
    except BaseException as exc:
        (work / "failure.json").write_text(json.dumps({"source_commit": commit, "error": str(exc),
                                                     "complete": False}, indent=2) + "\n")
        raise
    else:
        shutil.rmtree(work)


def render(report: dict) -> str:
    def f(value):
        return "undefined" if value is None else f"{value:,}" if isinstance(value, int) else f"{value:.6g}"

    primary = report["primary"]
    low, high = primary["percentile_interval"]
    lines = ["# Twenty-block 2D immune activation replication", "",
             "A prospectively frozen simulation study; these results are not independent biological validation.", "",
             f"The mean paired SDT − RSL3 activation contrast is **{f(primary['mean_paired_difference'])}** "
             f"accumulated activation steps per initial tumor cell. The 95% whole-block percentile "
             f"bootstrap interval is **[{f(low)}, {f(high)}]** ({primary['interpretation']} under the frozen design).", "",
             "Activation is DAMP/(DAMP+50), summed over eligible cell-step opportunities and divided by "
             "the initial tumor count. It is a composite exposure measure, not a probability, immunogenicity "
             "per death, or an assay-calibrated biological quantity. Zero opportunities give zero exposure.", "",
             "## All blocks", "", "| Block | Root seed | Control exposure | RSL3 exposure | SDT exposure | SDT − RSL3 |",
             "|---:|---:|---:|---:|---:|---:|"]
    for block in sorted(report["blocks"], key=lambda b: b["block"]):
        values = [block["arms"][arm][ENDPOINT] for arm in ARMS]
        lines.append(f"| {block['block']} | {block['root_seed']} | " + " | ".join(map(f, values))
                     + f" | {f(values[2] - values[1])} |")
    lines += ["", "## Distributions across blocks", "",
              "Each cell gives median [minimum, maximum]; n is the number of defined blocks. "
              "The JSON names every undefined block and retains every arm's full event-accounting summary. "
              "Empty conditional means, maxima and zero-denominator ratios are undefined, not zero.", "",
              "| Measurement | Control | RSL3 | SDT |", "|---|---|---|---|"]
    for key in sorted(report["arm_distributions"]["Control"]):
        cells = []
        for arm in ARMS:
            d = report["arm_distributions"][arm][key]
            cells.append(f"{f(d['median'])} [{f(d['minimum'])}, {f(d['maximum'])}]; n={d['n_defined']}")
        lines.append("| " + key.replace("_", " ") + " | " + " | ".join(cells) + " |")
    lines += ["", "## Threshold observations", "",
              "Counts refer to blocks containing any eligible opportunity at each activation threshold.", "",
              "| Arm | Activation ≥0.5 | Activation ≥0.9 |", "|---|---:|---:|"]
    for arm in ARMS:
        lines.append(f"| {arm} | " + " | ".join(f"{len(report['threshold_positive_blocks'][arm][k])}/20"
                                               for k in measurement.COUNTS) + " |")
    kill_ratio = report["sdt_rsl3_immune_kill_ratio"]
    lines += ["", f"The descriptive SDT/RSL3 immune-kill ratio has median {f(kill_ratio['median'])} "
              f"and range [{f(kill_ratio['minimum'])}, {f(kill_ratio['maximum'])}] among "
              f"{kill_ratio['n_defined']} defined blocks. Undefined blocks: "
              + (", ".join(map(str, kill_ratio["undefined_blocks"])) or "none") + "."]
    lines += ["", "Zero observed threshold-positive blocks would not establish impossibility. "
              "Secondary outcomes are descriptive; no confirmatory secondary significance tests were planned.", "",
              "## Design and limits", "",
              "Twenty fixed roots, 42 + block×2^32, retain the original three-arm model and within-block RNG "
              "structure. Complete blocks, never cells or arms, are resampled 10,000 times using Python "
              "random.Random(20260922). The interval uses linearly interpolated 2.5th and 97.5th percentiles "
              "of the mean paired difference. It is conditional on this model, configuration and "
              "independent-block interpretation; it excludes parameter, structural and experimental uncertainty.", "",
              "A source audit verifies disjoint additive seed-address envelopes without u64 wrapping. "
              "This removes known cross-block additive aliases; it does not prove statistical independence "
              "of pseudorandom streams. The earlier seeds 42–61 reuse seed addresses and their historical "
              "bootstrap interval should not be interpreted as validated independent-run uncertainty.", "",
              "The new binary reproduced both the full historical 33-condition summary and canonical "
              "observer bytes before any new block. All 20 planned blocks are included; there was no "
              "outcome-dependent replacement, tuning or stopping. The observation window, fixed initial "
              "tumor count and O2-independent SDT assumption limit interpretation. This does not establish "
              "a causal saturation mechanism or clinical efficacy.", "",
              "## Reproduction", "", "```bash", "python3 scripts/immune_2d_replication_report.py", "```", "",
              "This validates frozen sources and payload hashes, canonical parity, all 20 event ledgers, "
              "activation summaries and seed coverage, then recomputes the endpoint, interval and both reports "
              "without Rust. `--render-only` refreshes only Markdown from the saved derived JSON; it does "
              "not revalidate the archive. JSON and Markdown publication is sequential, not an atomic pair.", "",
              "- [Frozen protocol](../docs/IMMUNE_2D_REPLICATION_PROTOCOL.md).",
              "- [Plan](../scripts/immune_2d_replication_plan.json).",
              "- [Archive manifest](immune-2d-replication/manifest.json).",
              "- [Complete derived data](immune-2d-replication.json).",
              f"- Source freeze: `{report['provenance']['source_commit']}`.",
              f"- Binary SHA-256: `{report['provenance']['binary_sha256']}`.", ""]
    return "\n".join(lines)


def validate_paths(archive: Path) -> None:
    for output in (OUT_JSON, OUT_MD):
        common.validate_output_path(archive, output)
        require(not output.is_symlink(), "derived output must not be a symlink")
        for protected in (archive, measurement.ARCHIVE):
            if protected.exists():
                require(not output.resolve().is_relative_to(protected.resolve()), "protected archive output")
                if output.exists():
                    require(not any(p.is_file() and output.samefile(p) for p in protected.iterdir()),
                            "derived output aliases an archive payload")
    require(OUT_JSON.resolve() != OUT_MD.resolve()
            and not (OUT_JSON.exists() and OUT_MD.exists() and OUT_JSON.samefile(OUT_MD)), "distinct report outputs")
    require(not archive.resolve().is_relative_to(measurement.ARCHIVE.resolve())
            and not measurement.ARCHIVE.resolve().is_relative_to(archive.resolve()), "protect canonical archive")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--capture", action="store_true")
    mode.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    validate_paths(ARCHIVE)
    if args.capture:
        capture(ARCHIVE)
    report = json.loads(OUT_JSON.read_text()) if args.render_only else load_archive(ARCHIVE)
    markdown = render(report)
    json_text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    validate_paths(ARCHIVE)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    if not args.render_only:
        OUT_JSON.write_text(json_text)
    OUT_MD.write_text(markdown)
    print(OUT_MD)


if __name__ == "__main__":
    main()
