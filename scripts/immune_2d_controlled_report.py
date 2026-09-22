#!/usr/bin/env python3
"""Capture or reconstruct the frozen deterministic source/recipient comparison."""

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
import re
import shutil
import subprocess
import tarfile
import tempfile
import time

import immune_2d_measurement_report as measurement
import immune_measurement_report as common
from immune_measurement_report import finite, integer, require, sha

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "analysis/immune-2d-controlled"
OUT_JSON = ROOT / "analysis/immune-2d-controlled.json"
OUT_MD = ROOT / "analysis/immune-2d-controlled.md"
PLAN = "scripts/immune_2d_controlled_plan.json"
PROTOCOL = "docs/IMMUNE_2D_CONTROLLED_PROTOCOL.md"
CONTRACT = "scripts/immune_2d_measurement_config_v1.json"
CANONICAL_MANIFEST = "analysis/immune-2d-measurements/manifest.json"
STUDY = "immune-2d-controlled-source-recipient"
ENDPOINT = "activation_sum_per_reference_recipient"
PLAN_SHA256 = "0d3e0e701fbdde8f19e1a740c93d588e01c156686ff9ae8b9cc4cb7ca964ac0a"
RUST_INVENTORY_SHA256 = "3c01c856a56fe537554b9fc4d93f1a6f4c56777929bec659edcecd3e78069216"
UNITS = {"damp": "model field units", "time": "model steps; no physical-time conversion",
         "cell_size": "micrometers", "activation": "dimensionless",
         "activation_sum_per_reference_recipient": "activation-steps per fixed reference recipient"}
FINGERPRINT_DEFINITION = {
    "step": "FNV-1a-64 of row-major finite nonnegative f64 bit patterns, each encoded little-endian",
    "condition": "FNV-1a-64 of the ordered step fingerprint u64 values, each encoded little-endian",
    "cryptographic": False,
}
COUNTS = ("damp_ge_kd_opportunities", "damp_ge_9kd_opportunities")
REQUIRED_SOURCES = {
    PLAN, PROTOCOL, CONTRACT, CANONICAL_MANIFEST, measurement.GOLDEN, measurement.HISTORICAL,
    measurement.PROTOCOL, "scripts/immune_2d_measurement_report.py", "scripts/immune_measurement_report.py",
    "scripts/immune_2d_controlled_report.py", "tests/test_immune_2d_controlled_report.py",
}
ARTIFACTS = {"observations.json.gz", "baseline-summary.json", "canonical-observations.json.gz",
             "sources.tar.gz", "capture-logs.tar.gz"}


def exact(left, right, label: str) -> None:
    require(json.dumps(left, sort_keys=True, allow_nan=False)
            == json.dumps(right, sort_keys=True, allow_nan=False), label)


def validate_plan(plan: dict) -> None:
    require(sha(json.dumps(plan, sort_keys=True, allow_nan=False).encode()) == PLAN_SHA256,
            "frozen controlled plan changed")


def close(left, right, label: str, plan: dict) -> None:
    finite(left, label)
    finite(right, label)
    require(abs(left - right) <= plan["tolerances"]["absolute"]
            + plan["tolerances"]["relative"] * abs(right), label)


def at_most(left, right, label: str, plan: dict) -> None:
    finite(left, label)
    finite(right, label)
    require(left <= right + plan["tolerances"]["absolute"]
            + plan["tolerances"]["relative"] * abs(right), label)


def expected_masks(plan: dict) -> dict:
    cfg, mask = plan["config"], plan["masks"]
    start, stop, tile = mask["region_start"], mask["region_stop"], mask["tile_size"]
    sparse, dense, region, reduced = [], [], [], []
    for r in range(start, stop):
        for c in range(start, stop):
            idx = r * cfg["grid_cols"] + c
            rr, cc = r - start, c - start
            region.append(idx)
            if rr % tile in mask["source_offsets"] and cc % tile in mask["source_offsets"]:
                dense.append(idx)
                if (rr % tile == mask["source_offsets"][(rr // tile) % 2]
                        and cc % tile == mask["source_offsets"][(cc // tile) % 2]):
                    sparse.append(idx)
            elif rr % 2 == 0 and cc % 2 == 0:
                reduced.append(idx)
    reserved = set(dense)
    return {"reserved_source_indices": dense, "source_indices_64": sparse,
            "source_indices_256": dense, "reference_recipient_indices": [i for i in region if i not in reserved],
            "available_recipient_indices_960": reduced}


def checked_response(row: dict, n: int, plan: dict) -> None:
    """Bounds available from aggregate observations, not a transport replay."""
    cfg = plan["config"]
    activation = finite(row["activation_sum"], "finite activation sum")
    damp_sum = finite(row["local_damp_sum"], "local DAMP sum")
    lower, upper = (integer(row[k], k) for k in COUNTS)
    require(upper <= lower <= n, "nested activation threshold counts")
    maximum = row["max_eligible_damp"]
    if n == 0:
        require(activation == 0 and maximum is None and row["local_damp_sum"] == 0,
                "empty response must have zero sums and null maximum")
        return
    finite(maximum, "finite eligible maximum")
    require(maximum >= cfg["damp_kill_threshold"], "maximum below eligibility floor")
    floor = cfg["damp_kill_threshold"] / (cfg["damp_kill_threshold"] + cfg["dc_activation_kd"])
    ceiling = maximum / (maximum + cfg["dc_activation_kd"])
    eps = plan["tolerances"]["absolute"] + plan["tolerances"]["relative"] * n
    require(n * floor - eps <= activation <= n * ceiling + eps, "activation bounds from eligibility and maximum")
    require(activation + eps >= 0.5 * lower + 0.4 * upper, "threshold counts exceed activation sum")
    for count, threshold in zip((lower, upper), (cfg["dc_activation_kd"], 9 * cfg["dc_activation_kd"])):
        require((count > 0) == (maximum >= threshold), "threshold maximum consistency")
    eps = plan["tolerances"]["absolute"] + plan["tolerances"]["relative"] * damp_sum
    require(maximum + (n - 1) * cfg["damp_kill_threshold"] - eps <= damp_sum <= n * maximum + eps,
            "DAMP sum and maximum")
    require(damp_sum + eps >= n * cfg["damp_kill_threshold"], "DAMP eligibility sum")
    require(activation <= damp_sum / cfg["dc_activation_kd"] + eps, "activation exceeds DAMP sum bound")
    kd, threshold = cfg["dc_activation_kd"], cfg["damp_kill_threshold"]
    at_most(ceiling + (n - 1) * floor, activation, "activation below attained maximum", plan)
    at_most(upper * 9 * kd + (lower - upper) * kd + (n - lower) * threshold,
            damp_sum, "threshold counts exceed DAMP sum", plan)
    at_most(damp_sum, upper * maximum + (lower - upper) * min(maximum, 9 * kd)
            + (n - lower) * min(maximum, kd), "DAMP sum exceeds threshold bins", plan)
    at_most(damp_sum / (maximum + kd), activation, "activation below DAMP bound", plan)
    at_most(activation, n * damp_sum / (damp_sum + n * kd), "activation exceeds concave response bound", plan)
    if n > 1:
        remaining = max(0.0, damp_sum - maximum)
        at_most(activation, ceiling + (n - 1) * remaining / (remaining + (n - 1) * kd),
                "activation exceeds attained-maximum bound", plan)
        if n == 2:
            close(activation, ceiling + remaining / (remaining + kd), "two-opportunity activation", plan)
    else:
        close(activation, ceiling, "single-opportunity activation", plan)


def fingerprint(value) -> None:
    require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{16}", value) is not None,
            "field fingerprint format")


def trajectory_fingerprint(steps: list[dict]) -> str:
    value = 14695981039346656037
    for step in steps:
        fingerprint(step["field_fingerprint"])
        for byte in int(step["field_fingerprint"], 16).to_bytes(8, "little"):
            value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return f"{value:016x}"


def validate_observations(data: dict, plan: dict) -> list[dict]:
    validate_plan(plan)
    for key in ("schema_version", "study", "simulator", "dimension", "config"):
        exact(data[key], plan[key], f"observation {key}")
    exact(data["mask_config"], plan["masks"], "mask configuration")
    exact(data["units"], UNITS, "observation units")
    exact(data["field_fingerprint_definition"], FINGERPRINT_DEFINITION, "fingerprint definition")
    masks = expected_masks(plan)
    exact(data["masks"], masks, "exact source and recipient masks")
    exact([r["condition_id"] for r in data["conditions"]],
          [r["condition_id"] for r in plan["conditions"]], "ordered complete condition coverage")
    cfg = plan["config"]
    window = cfg["n_steps"] - cfg["immune_start_step"]
    reference = masks["reference_recipient_indices"]
    reduced = set(masks["available_recipient_indices_960"])
    summaries = []
    for condition, intended in zip(data["conditions"], plan["conditions"]):
        exact({k: condition[k] for k in intended}, intended, "condition identity")
        steps, cells = condition["steps"], condition["recipients"]
        exact([s["step"] for s in steps], list(range(cfg["n_steps"])), "complete ordered steps")
        exact([c["cell_index"] for c in cells], reference, "complete ordered recipient identities")
        fingerprint(condition["field_fingerprint"])
        require(condition["field_fingerprint"] == trajectory_fingerprint(steps), "trajectory fingerprint chain")
        expected_mass = 0.0
        for s in steps:
            t = s["step"]
            injected = intended["total_damp"] if t == intended["release_step"] else 0.0
            exact(s["injected_damp"], injected, "release amount and phase")
            expected_mass = (expected_mass + injected) * (1 - cfg["damp_clearance_rate"])
            close(s["field_mass"], expected_mass, "whole-field mass recurrence", plan)
            if expected_mass == 0:
                require(s["field_mass"] == 0, "zero-release field must remain exactly zero")
            available = integer(s["available_cells"], "available cells")
            require(available == intended["recipient_count"], "fixed availability")
            eligible = integer(s["eligible_cells"], "eligible cells")
            require(eligible <= available, "eligible cells exceed availability")
            if t < cfg["immune_start_step"] or expected_mass == 0:
                require(eligible == 0, "opportunity before immune window or without DAMP")
            checked_response(s, eligible, plan)
            at_most(s["local_damp_sum"], s["field_mass"], "recipient DAMP exceeds whole-field mass", plan)
            fingerprint(s["field_fingerprint"])
        for cell in cells:
            available = intended["recipient_count"] == len(reference) or cell["cell_index"] in reduced
            require(type(cell["available"]) is bool and cell["available"] == available, "recipient availability mask")
            n = integer(cell["opportunities"], "recipient opportunities")
            require(n <= (window if available else 0), "recipient opportunity bound")
            checked_response(cell, n, plan)
        activation = math.fsum(c["activation_sum"] for c in cells)
        close(activation, math.fsum(s["activation_sum"] for s in steps), "cell/step activation reconciliation", plan)
        damp_sum = math.fsum(c["local_damp_sum"] for c in cells)
        close(damp_sum, math.fsum(s["local_damp_sum"] for s in steps), "cell/step DAMP reconciliation", plan)
        opportunities = sum(c["opportunities"] for c in cells)
        require(opportunities == sum(s["eligible_cells"] for s in steps), "cell/step opportunity reconciliation")
        maximum = max((c["max_eligible_damp"] for c in cells if c["opportunities"]), default=None)
        step_max = max((s["max_eligible_damp"] for s in steps if s["eligible_cells"]), default=None)
        exact(maximum, step_max, "cell/step maximum reconciliation")
        summary = {"activation_sum": activation, ENDPOINT: activation / len(reference),
                   "opportunities": opportunities,
                   "unique_eligible_cells": sum(c["opportunities"] > 0 for c in cells),
                   "max_eligible_damp": maximum, "final_field_mass": steps[-1]["field_mass"]}
        for key in COUNTS:
            summary[key] = sum(c[key] for c in cells)
            require(summary[key] == sum(s[key] for s in steps), "cell/step threshold reconciliation")
        require(set(condition["summary"]) == set(summary), "summary fields")
        for key, value in summary.items():
            if value is None or type(value) is int:
                exact(condition["summary"][key], value, f"summary {key}")
            else:
                close(condition["summary"][key], value, f"summary {key}", plan)
        summary.update(intended)
        summary["conditional_mean_activation"] = activation / opportunities if opportunities else None
        summary["local_damp_sum"] = damp_sum
        summary["conditional_mean_damp"] = damp_sum / opportunities if opportunities else None
        summary["threshold_fractions"] = {key: summary[key] / opportunities if opportunities else None for key in COUNTS}
        summary["activation_sum_per_available_recipient"] = activation / intended["recipient_count"]
        summaries.append(summary)
    # Availability changes the receiver mask only, including the two zero controls.
    for full, subset in zip(data["conditions"][::2], data["conditions"][1::2]):
        require(full["field_fingerprint"] == subset["field_fingerprint"], "availability changed field trajectory")
        for left, right in zip(full["steps"], subset["steps"]):
            for key in ("field_mass", "field_fingerprint", "injected_damp"):
                exact(left[key], right[key], "availability field invariance")
        for left, right in zip(full["recipients"], subset["recipients"]):
            if right["available"]:
                exact(left, right, "available subset recipient invariance")
    return summaries


def summarize(summaries: list[dict], plan: dict) -> dict:
    validate_plan(plan)
    exact([s["condition_id"] for s in summaries], [s["condition_id"] for s in plan["conditions"]],
          "summary condition coverage")
    by_id = {s["condition_id"]: s for s in summaries}
    contrasts = []
    for comparison in plan["contrasts"]:
        left, right = (by_id[comparison[k]] for k in ("left", "right"))
        contrasts.append({**comparison, "endpoint": ENDPOINT,
                          "difference": left[ENDPOINT] - right[ENDPOINT]})
    return {"conditions": summaries, "contrasts": contrasts}


def pack_files(contents: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as bundle:
        for name, blob in sorted(contents.items()):
            entry = tarfile.TarInfo(name)
            entry.size, entry.mode = len(blob), 0o644
            bundle.addfile(entry, io.BytesIO(blob))
    return gzip.compress(buffer.getvalue(), mtime=0)


def source_paths() -> list[str]:
    return sorted(set(measurement.source_paths()) | REQUIRED_SOURCES)


def unpack_sources(blob: bytes, hashes: dict) -> dict[str, bytes]:
    require(REQUIRED_SOURCES <= set(hashes), "required frozen source inventory")
    rust_paths = sorted(p for p in hashes if p.startswith("simulations/") and p.endswith((".rs", ".toml", ".lock")))
    require(sha(("\n".join(rust_paths) + "\n").encode()) == RUST_INVENTORY_SHA256, "complete frozen Rust source inventory")
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as bundle:
        members = bundle.getmembers()
        require(len(members) == len(hashes) and {m.name for m in members} == set(hashes), "source bundle inventory")
        contents = {}
        for member in members:
            require(member.isfile() and not Path(member.name).is_absolute() and ".." not in Path(member.name).parts,
                    "invalid source archive entry")
            contents[member.name] = bundle.extractfile(member).read()
            require(sha(contents[member.name]) == hashes[member.name], f"source hash: {member.name}")
    return contents


def load_archive(path: Path) -> dict:
    require(path.is_dir() and not path.is_symlink(), "regular archive directory")
    require({p.name for p in path.iterdir()} == ARTIFACTS | {"manifest.json"}, "archive file inventory")
    require(all(p.is_file() and not p.is_symlink() for p in path.iterdir()), "regular archive files")
    manifest = json.loads((path / "manifest.json").read_text())
    require(type(manifest["schema_version"]) is int and manifest["schema_version"] == 1
            and manifest["study"] == STUDY and manifest["execution"] == "serial"
            and manifest["complete"] is True, "complete archive schema")
    for key, width in (("source_commit", 40), ("binary_sha256", 64)):
        require(isinstance(manifest[key], str) and re.fullmatch(r"[0-9a-f]{%d}" % width, manifest[key]),
                "capture source/binary provenance")
    require(set(manifest["artifacts"]) == ARTIFACTS, "manifest artifact inventory")
    blobs = {name: (path / name).read_bytes() for name in ARTIFACTS}
    for name, blob in blobs.items():
        require(sha(blob) == manifest["artifacts"][name], f"artifact hash: {name}")
    frozen = unpack_sources(blobs["sources.tar.gz"], manifest["sources"])
    plan = json.loads(frozen[PLAN])
    validate_plan(plan)
    require(sha(blobs["baseline-summary.json"]) == measurement.golden_sha(frozen[measurement.GOLDEN].decode()),
            "full baseline parity")
    exact(manifest["baseline_reference"], measurement.baseline_reference(frozen[measurement.GOLDEN].decode()),
          "baseline provenance")
    require(sha(blobs["canonical-observations.json.gz"])
            == json.loads(frozen[CANONICAL_MANIFEST])["artifacts"]["observations.json.gz"], "canonical observation parity")
    measurement.validate_observations(json.loads(gzip.decompress(blobs["canonical-observations.json.gz"])),
                                      json.loads(blobs["baseline-summary.json"]), contract=json.loads(frozen[CONTRACT]),
                                      historical=json.loads(frozen[measurement.HISTORICAL]))
    exact([r["name"] for r in manifest["runs"]], ["baseline", "canonical", "controlled"], "capture run coverage")
    for run in manifest["runs"]:
        finite(run["elapsed_seconds"], "capture duration")
    data = json.loads(gzip.decompress(blobs["observations.json.gz"]))
    summaries = validate_observations(data, plan)
    return {"schema_version": 1, "study": STUDY, "plan": plan, "provenance": manifest,
            **summarize(summaries, plan)}


def capture(destination: Path) -> None:
    require(not destination.exists() and not destination.is_symlink(), "immutable capture destination already exists")
    require(not any(k.startswith("FERRO_") for k in os.environ), "unset all FERRO_* overrides")
    require(not subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip(),
            "commit the protocol, plan, implementation and tests before capture")
    sources = {p: (ROOT / p).read_bytes() for p in source_paths()}
    plan = json.loads(sources[PLAN])
    validate_plan(plan)
    canonical = (measurement.ARCHIVE / "observations.json.gz").read_bytes()
    require(sha(canonical) == json.loads(sources[CANONICAL_MANIFEST])["artifacts"]["observations.json.gz"],
            "historical canonical archive hash")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    rustc = subprocess.check_output(["rustup", "run", common.TOOLCHAIN, "rustc", "--version"], text=True).strip()
    destination.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".immune-controlled-capture-", dir=destination.parent))
    print(f"Capture evidence (retained on failure): {work}", flush=True)
    staged = work / "archive"
    staged.mkdir()
    try:
        binary = common.build_binary(["rustup", "run", common.TOOLCHAIN, "cargo"], ROOT / "simulations", "sim-tme")
        private = work / "sim-tme"
        shutil.copy2(binary, private)
        binary_sha = sha(private.read_bytes())
        logs, runs = {}, []

        def execute(name: str, args: list[str], output: str) -> bytes:
            cwd = work / name
            cwd.mkdir()
            start = time.monotonic()
            log = work / f"{name}.log"
            with log.open("wb") as stream:
                subprocess.run([str(private), *args], cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, check=True)
            logs[log.name] = log.read_bytes()
            raw = (cwd / "output/tme" / output).read_bytes()
            elapsed = time.monotonic() - start
            runs.append({"name": name, "elapsed_seconds": elapsed})
            print(f"Completed {name} ({elapsed:.1f}s)", flush=True)
            return raw

        baseline = execute("baseline", [], "tme_summary.json")
        (staged / "baseline-summary.json").write_bytes(baseline)
        require(sha(baseline) == measurement.golden_sha(sources[measurement.GOLDEN].decode()), "full baseline parity")
        raw = execute("canonical", ["--immune-measurements"], "immune_measurements.json")
        require(raw == gzip.decompress(canonical), "canonical observer byte parity")
        measurement.validate_observations(json.loads(raw), json.loads(baseline), contract=json.loads(sources[CONTRACT]),
                                          historical=json.loads(sources[measurement.HISTORICAL]))
        (staged / "canonical-observations.json.gz").write_bytes(canonical)
        raw = execute("controlled", ["--immune-controlled"], "immune_controlled.json")
        (staged / "observations.json.gz").write_bytes(gzip.compress(raw, mtime=0))
        validate_observations(json.loads(raw), plan)
        require(sources == {p: (ROOT / p).read_bytes() for p in source_paths()}, "frozen sources changed during capture")
        require(commit == subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "source commit changed during capture")
        (staged / "sources.tar.gz").write_bytes(pack_files(sources))
        (staged / "capture-logs.tar.gz").write_bytes(pack_files(logs))
        manifest = {"schema_version": 1, "study": STUDY, "complete": True,
                    "source_commit": commit, "captured_at_utc": datetime.now(timezone.utc).isoformat(),
                    "platform": platform.platform(), "rustc": rustc, "execution": "serial",
                    "python": {"implementation": platform.python_implementation(), "version": platform.python_version()},
                    "binary_sha256": binary_sha, "runs": runs,
                    "baseline_reference": measurement.baseline_reference(sources[measurement.GOLDEN].decode()),
                    "sources": {p: sha(blob) for p, blob in sources.items()},
                    "artifacts": {p: sha((staged / p).read_bytes()) for p in sorted(ARTIFACTS)}}
        (staged / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        render(load_archive(staged))
        require(not destination.exists() and not destination.is_symlink(), "capture destination appeared during run")
        staged.rename(destination)
    except BaseException as exc:
        (work / "failure.json").write_text(json.dumps({"source_commit": commit, "error": str(exc),
                                                     "complete": False}, indent=2) + "\n")
        raise
    else:
        shutil.rmtree(work)


def render(report: dict) -> str:
    def f(value):
        return "undefined" if value is None else f"{value:.7g}"

    lines = ["# Controlled 2D DAMP source/recipient comparison", "",
             "A deterministic experiment on the existing transport and activation rules. "
             "These artificial source and recipient masks are not biological validation or an attribution "
             "of the historical SDT–RSL3 immune-kill contrast.", "",
             "The primary response sums DAMP/(DAMP+50) at eligible, available recipient sites during "
             "steps 60–179 and divides by the same 3,840 reference sites in every condition. "
             "Units are activation steps per reference recipient. Zero opportunities give zero exposure; "
             "conditional means remain undefined. No stochastic kills, source feedback or inference intervals are used.", "",
             "## Every planned condition", "",
             "| Condition | DAMP/source | Reference exposure | Opportunities | Conditional activation | Max eligible DAMP | ≥50 | ≥450 |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for c in report["conditions"]:
        fields = ("per_source_damp", ENDPOINT, "opportunities", "conditional_mean_activation",
                  "max_eligible_damp", *COUNTS)
        lines.append(f"| {c['condition_id']} | " + " | ".join(f(c[k]) for k in fields) + " |")
    lines += ["", "Condition IDs name source count (n), total injection (m), release step (t), "
              "and available recipients (r). Both zero controls use a zero injection at step 60.", "",
              "## All prespecified contrasts", "",
              "Differences below are left minus right on the primary response. Equal-total-mass comparisons "
              "change source positions/count and the amount per source together. Equal-per-source-amount "
              "comparisons change source count and total injection together. Neither is a pure density effect "
              "holding both quantities fixed. Availability comparisons use an unchanged field and a fixed subset.", "",
              "| Comparison | Left | Right | Difference |", "|---|---|---|---:|"]
    for c in report["contrasts"]:
        lines.append(f"| {c['kind']} | {c['left']} | {c['right']} | {f(c['difference'])} |")
    lines += ["", "## Accounting and limits", "",
              "All planned conditions and contrasts are published. Release occurs before transport and "
              "clearance, observation after both. Whole-field mass follows M[t]=0.97×(M[t−1]+injection[t]); "
              "the rectangular boundary loses no mass to unrepresented neighbors. Recipient availability "
              "does not remove DAMP or alter transport. Cell and step ledgers reconcile; selected recipients "
              "retain identical observations in both availability arms.", "",
              "The archived observations contain per-step and per-recipient aggregates, not every field value "
              "at every step. Offline reconstruction checks those ledgers and analytic invariants; it cannot "
              "independently recompute transport or recover each activation evaluation. FNV fingerprints "
              "are noncryptographic consistency checks; SHA-256 identifies the archived payloads. Replay the "
              "frozen Rust driver for full transport reproduction. Neither check establishes biological accuracy.", "",
              "DAMP is in model field units. Model steps have no new conversion to hours. The 0.001 diffusion "
              "cutoff and 0.01 eligibility floor remain active; general linear-superposition claims do not follow. "
              "The masks and timing are artificial, and the result is conditional on this finite design. "
              "It supplies no per-death immunogenicity estimate, calibrated assay, clinical conclusion or P5 replacement.", "",
              "## Reproduction", "", "```bash", "python3 scripts/immune_2d_controlled_report.py", "```", "",
              "The default command validates the archive and regenerates both reports without running Rust. "
              "`--render-only` also validates and reconstructs the archive, then writes only Markdown. "
              "Report files are written sequentially, not as an atomic pair. A new capture requires a clean "
              "committed checkout, an absent archive destination, and unchanged historical production and observer outputs.", "",
              "- [Frozen protocol](../docs/IMMUNE_2D_CONTROLLED_PROTOCOL.md).",
              "- [Plan](../scripts/immune_2d_controlled_plan.json).",
              "- [Archive manifest](immune-2d-controlled/manifest.json).",
              "- [Complete derived data](immune-2d-controlled.json).",
              f"- Source freeze: `{report['provenance']['source_commit']}`.",
              f"- Binary SHA-256: `{report['provenance']['binary_sha256']}`.", ""]
    return "\n".join(lines)


def validate_paths(archive: Path) -> None:
    for output in (OUT_JSON, OUT_MD):
        common.validate_output_path(archive, output)
        require(not output.is_symlink(), "derived output must not be a symlink")
        for protected in (archive, measurement.ARCHIVE):
            require(not output.resolve().is_relative_to(protected.resolve()), "protected archive output")
            if protected.exists() and output.exists():
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
    report = load_archive(ARCHIVE)
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
