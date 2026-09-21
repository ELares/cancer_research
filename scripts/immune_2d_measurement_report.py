#!/usr/bin/env python3
"""Capture or reconcile passive observations of the historical 2D immune comparison.

The default reconstructs a frozen archive without Rust or simulation execution.
Capture requires committed sources and an unchanged full 33-condition baseline.
See docs/IMMUNE_2D_MEASUREMENT_PROTOCOL.md for the prospective interpretation.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import heapq
import io
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tarfile
import tempfile

import immune_measurement_report as common
from immune_measurement_report import (at_most, close, finite, integer, ratio,
                                       report_link, require, sha, validate_output_path)

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "analysis/immune-2d-measurements"
REPORT = ROOT / "analysis/immune-2d-measurement-report.md"
CONFIG = ROOT / "scripts/immune_2d_measurement_config_v1.json"
GOLDEN = "simulations/sim-tme/expected_summary.sha256"
HISTORICAL = "tests/fixtures/flagship_tme_rows.json"
PROTOCOL = "docs/IMMUNE_2D_MEASUREMENT_PROTOCOL.md"
NAMES = ("immune_Control", "immune_RSL3", "immune_SDT")
ARTIFACTS = common.ARTIFACTS
TOOLCHAIN = common.TOOLCHAIN
SCHEMA = 2
COUNTS = ("damp_ge_kd_opportunities", "damp_ge_9kd_opportunities")


def golden_sha(text: str) -> str:
    return next(line.split()[0] for line in text.splitlines() if line and not line.startswith("#"))


def baseline_reference(text: str) -> dict:
    labels = {"source_commit": "# Source commit: ", "binary_sha256": "# Binary SHA-256: ",
              "toolchain": "# Toolchain: "}
    reference = {key: next(line.removeprefix(prefix) for line in text.splitlines()
                           if line.startswith(prefix)) for key, prefix in labels.items()}
    return dict(reference, summary_sha256=golden_sha(text))


def baseline_rows(baseline: dict) -> dict[str, dict]:
    rows = [r for r in baseline["conditions"] if r.get("o2_condition") == "gradient_120um"
            and r.get("o2_lambda_um") == 120.0 and r.get("immune_mode") == "immune_on"
            and r.get("stromal_mode") == "off" and "ph_mode" not in r]
    require(len(rows) == 3 and {r["treatment"] for r in rows} == {"Control", "RSL3", "SDT"},
            "exactly one historical baseline immune row per treatment")
    return {r["treatment"]: r for r in rows}


def validate_baseline(baseline: dict, historical: dict) -> dict[str, dict]:
    require(baseline["schema_version"] == 1 and len(baseline["conditions"]) == 33,
            "complete 33-condition baseline matrix")
    expected = baseline_rows(baseline)
    require(expected == baseline_rows(historical), "historical 2D fixture changed")
    return expected


def check_activation(population: dict, count: int, damp_sum: float, cfg: dict) -> None:
    """Bound aggregate activation using its recorded DAMP population, not a fitted hazard."""
    kd = cfg["immune_config"]["dc_activation_kd"]
    threshold = cfg["damp_kill_threshold"]
    activation = finite(population["activation_sum"], "activation sum")
    ge_kd, ge_9kd = (integer(population[k], k) for k in COUNTS)
    maximum = population["max_local_damp"]
    require(ge_9kd <= ge_kd <= count, "nested activation thresholds")
    if count == 0:
        require(activation == 0 and maximum is None and ge_kd == ge_9kd == 0,
                "empty activation population")
        return
    finite(maximum, "maximum eligible DAMP")
    require(maximum >= threshold, "maximum below eligibility threshold")
    at_most(maximum + (count - 1) * threshold, damp_sum, "maximum exceeds DAMP sum")
    at_most(damp_sum, count * maximum, "DAMP sum exceeds recorded maximum")
    require((ge_kd > 0) == (maximum >= kd) and (ge_9kd > 0) == (maximum >= 9 * kd),
            "maximum disagrees with threshold counts")
    lower_damp = ge_9kd * 9 * kd + (ge_kd - ge_9kd) * kd + (count - ge_kd) * threshold
    upper_damp = (ge_9kd * maximum + (ge_kd - ge_9kd) * min(maximum, 9 * kd)
                  + (count - ge_kd) * min(maximum, kd))
    at_most(lower_damp, damp_sum, "threshold counts exceed DAMP sum")
    at_most(damp_sum, upper_damp, "DAMP sum exceeds threshold bins")
    low_activation = (ge_9kd * .9 + (ge_kd - ge_9kd) * .5
                      + (count - ge_kd) * threshold / (threshold + kd))
    at_most(low_activation, activation, "activation below threshold populations")
    at_most(damp_sum / (maximum + kd), activation, "activation below DAMP bound")
    at_most(activation, count * damp_sum / (damp_sum + count * kd),
            "activation exceeds concave-response bound")
    if count > 1:
        remaining = max(0., damp_sum - maximum)
        attained_activation = maximum / (maximum + kd)
        at_most(activation, attained_activation
                + (count - 1) * remaining / (remaining + (count - 1) * kd),
                "activation exceeds attained-maximum bound")
        if count == 2:
            close(activation, attained_activation + remaining / (remaining + kd),
                  "two-opportunity activation")
    if count == 1:
        close(maximum, damp_sum, "single-opportunity maximum")
        close(activation, damp_sum / (damp_sum + kd), "single-opportunity activation")


def reconcile_condition(row: dict, cfg: dict) -> dict:
    n_cells = cfg["grid_rows"] * cfg["grid_cols"]
    total = finite(row["final_damp"]["total"], "final field total")
    peak = finite(row["final_damp"]["peak"], "final field peak")
    at_most(peak, total, "final field peak exceeds total")
    at_most(total, n_cells * peak, "final field total exceeds peak bound")
    # Explicit adapter for shared event physics; raw 2D configuration/results
    # retain their own schema, including the separate final field observation.
    accounting_cfg = dict(cfg, immune=cfg["immune_config"],
                          post_death_steps=cfg["params"]["post_death_steps"])
    summary = common.reconcile_condition(row, accounting_cfg, n_cells=n_cells,
                                         final_damp=total, runtime_seed=row["seed"])
    obs = row["measurements"]
    at_most(max((event["terminal_damp"] for event in obs["ferroptotic_events"]
                 if event["terminal_damp"] is not None), default=0.), peak,
            "terminal addition exceeds final field peak")
    cells, steps = obs["eligible_cells"], obs["steps"]
    for cell in cells:
        check_activation(cell, cell["opportunities"], cell["local_damp_sum"], cfg)
    for step in steps:
        check_activation(step, step["eligible_cells"], step["eligible_local_damp_sum"], cfg)
    close(math.fsum(c["activation_sum"] for c in cells),
          math.fsum(s["activation_sum"] for s in steps), "activation population reconciliation")
    for key in COUNTS:
        require(sum(c[key] for c in cells) == sum(s[key] for s in steps),
                "activation threshold population reconciliation")
    max_cell = max((c["max_local_damp"] for c in cells), default=None)
    max_step = max((s["max_local_damp"] for s in steps if s["max_local_damp"] is not None), default=None)
    require(max_cell == max_step, "maximum eligible DAMP reconciliation")
    # Maxima must be possible within each cell's observed eligibility interval.
    # Range maxima and an active-interval heap avoid a cell-by-step Cartesian scan.
    interval_maxima = [[0.] * len(steps) for _ in steps]
    for first in range(len(steps)):
        maximum = 0.
        for last in range(first, len(steps)):
            maximum = max(maximum, steps[last]["max_local_damp"] or 0.)
            interval_maxima[first][last] = maximum
    for cell in cells:
        at_most(cell["max_local_damp"], interval_maxima[cell["first_step"]][cell["last_step"]],
                "cell DAMP maximum exceeds eligible step maxima")
    ordered = iter(sorted(cells, key=lambda c: c["first_step"]))
    pending, active = next(ordered, None), []
    for step in steps:
        while pending is not None and pending["first_step"] <= step["step"]:
            heapq.heappush(active, (-pending["max_local_damp"], pending["last_step"]))
            pending = next(ordered, None)
        while active and active[0][1] < step["step"]:
            heapq.heappop(active)
        at_most(step["max_local_damp"] or 0., -active[0][0] if active else 0.,
                "step DAMP maximum exceeds eligible cell maxima")
    by_cell = {c["cell_index"]: c for c in cells}
    kd = cfg["immune_config"]["dc_activation_kd"]
    floor_activation = cfg["damp_kill_threshold"] / (cfg["damp_kill_threshold"] + kd)
    kills_by_step = [[] for _ in steps]
    for kill in obs["immune_kill_events"]:
        cell = by_cell[kill["cell_index"]]
        local = kill["local_damp"]
        activation = local / (local + kd)
        at_most(local, cell["max_local_damp"], "kill exceeds cell DAMP maximum")
        at_most(activation + (cell["opportunities"] - 1) * floor_activation,
                cell["activation_sum"], "kill exceeds cell activation sum")
        for key, cut in zip(COUNTS, (kd, 9 * kd)):
            require(cell[key] >= int(local >= cut), "kill exceeds cell threshold count")
        kills_by_step[kill["step"]].append(kill)
    for step, kills in zip(steps, kills_by_step):
        if not kills:
            continue
        kill_activation = math.fsum(k["local_damp"] / (k["local_damp"] + kd) for k in kills)
        at_most(max(k["local_damp"] for k in kills), step["max_local_damp"],
                "kill exceeds step DAMP maximum")
        at_most(kill_activation + (step["eligible_cells"] - len(kills)) * floor_activation,
                step["activation_sum"], "kill exceeds step activation sum")
        for key, cut in zip(COUNTS, (kd, 9 * kd)):
            require(step[key] >= sum(k["local_damp"] >= cut for k in kills),
                    "kill exceeds step threshold count")
        if len(kills) == step["eligible_cells"]:
            close(kill_activation, step["activation_sum"], "all-eligible killed activation")
    opportunities = summary["eligible_cell_steps"]
    summary.update(final_damp_peak=peak, max_eligible_damp=max_cell,
                   mean_activation=ratio(math.fsum(c["activation_sum"] for c in cells), opportunities))
    for key in COUNTS:
        summary[key] = sum(c[key] for c in cells)
        summary[key + "_fraction"] = ratio(summary[key], opportunities)
    return summary


def in_tumor_circle(index: int, cfg: dict) -> bool:
    r, c = divmod(index, cfg["grid_cols"])
    return ((r - cfg["grid_rows"] / 2) ** 2 + (c - cfg["grid_cols"] / 2) ** 2
            <= (cfg["tumor_radius_um"] / cfg["cell_size_um"]) ** 2)


def validate_observations(data: dict, baseline: dict, *, contract: dict | None = None,
                          historical: dict | None = None) -> list[dict]:
    require(data["schema_version"] == SCHEMA and data["simulator"] == "sim-tme"
            and data["dimension"] == 2, "2D measurement schema")
    contract = json.loads(CONFIG.read_text()) if contract is None else contract
    historical = json.loads((ROOT / HISTORICAL).read_text()) if historical is None else historical
    require(contract["schema_version"] == SCHEMA, "2D configuration contract schema")
    cfg = data["config"]
    require(cfg == contract["config"], "complete canonical 2D configuration changed")
    expected = validate_baseline(baseline, historical)
    require([r["condition_name"] for r in data["conditions"]] == list(NAMES), "canonical 2D conditions")
    summaries = []
    tumor_count = sum(in_tumor_circle(i, cfg) for i in range(cfg["grid_rows"] * cfg["grid_cols"]))
    for treatment_index, row in enumerate(data["conditions"]):
        treatment = row["result"]["treatment"]
        require(row["condition_name"] == "immune_" + treatment, "2D condition identity")
        require(row["seed"] == cfg["seed"] + treatment_index * cfg["rng"]["treatment_seed_stride"],
                "canonical 2D runtime seed")
        require(row["result"] == expected[treatment], "2D baseline row changed")
        require(row["result"]["total_tumor"] == tumor_count, "2D tumor census")
        summary = reconcile_condition(row, cfg)
        for population in ("ferroptotic_events", "immune_kill_events", "eligible_cells"):
            for cell in row["measurements"][population]:
                require(in_tumor_circle(cell["cell_index"], cfg), "observed cell outside 2D tumor circle")
        summaries.append(summary)
    return summaries


def source_paths() -> list[str]:
    tracked = subprocess.check_output(["git", "ls-files", "simulations"], cwd=ROOT, text=True).splitlines()
    return sorted([p for p in tracked if p.endswith((".rs", ".toml", ".lock"))] + [
        "scripts/immune_2d_measurement_report.py", "scripts/immune_measurement_report.py",
        str(CONFIG.relative_to(ROOT)), PROTOCOL, GOLDEN, HISTORICAL,
    ])


def load_archive(path: Path) -> tuple[dict, list[dict]]:
    manifest = json.loads((path / "manifest.json").read_text())
    require(manifest["schema_version"] == SCHEMA and manifest["simulator"] == "sim-tme"
            and manifest["dimension"] == 2, "2D archive schema")
    require(manifest["execution"] == "serial", "2D execution mode")
    require(set(manifest["artifacts"]) == ARTIFACTS, "archive artifact set")
    require({p.name for p in path.iterdir()} == ARTIFACTS | {"manifest.json"}, "unexpected archive files")
    blobs = {name: (path / name).read_bytes() for name in ARTIFACTS}
    for name, blob in blobs.items():
        require(sha(blob) == manifest["artifacts"][name], f"artifact hash: {name}")
    with tarfile.open(fileobj=io.BytesIO(blobs["sources.tar.gz"]), mode="r:gz") as bundle:
        members = bundle.getmembers()
        require(len(members) == len(manifest["sources"]), "source inventory size")
        require({m.name for m in members} == set(manifest["sources"]), "source inventory")
        frozen = {}
        for member in members:
            require(member.isfile() and not Path(member.name).is_absolute()
                    and ".." not in Path(member.name).parts, "invalid source archive entry")
            frozen[member.name] = bundle.extractfile(member).read()
            require(sha(frozen[member.name]) == manifest["sources"][member.name],
                    f"frozen source hash: {member.name}")
    require(sha(blobs["baseline-summary.json"]) == golden_sha(frozen[GOLDEN].decode()),
            "frozen 2D production SHA changed")
    require(manifest["baseline_reference"] == baseline_reference(frozen[GOLDEN].decode()),
            "frozen baseline provenance changed")
    contract = json.loads(frozen[str(CONFIG.relative_to(ROOT))])
    historical = json.loads(frozen[HISTORICAL])
    data = json.loads(gzip.decompress(blobs["observations.json.gz"]))
    return manifest, validate_observations(data, json.loads(blobs["baseline-summary.json"]),
                                          contract=contract, historical=historical)


def render(manifest: dict, summaries: list[dict], archive: Path = ARCHIVE,
           report: Path = REPORT) -> str:
    def f(value):
        return "undefined (n=0)" if value is None else f"{value:,}" if isinstance(value, int) else f"{value:.6g}"

    labels = [("Tumor cells", "total_tumor"), ("Ferroptotic deaths", "ferroptotic_deaths"),
              ("Completed releases", "completed_releases"), ("Horizon-censored releases", "censored_deaths"),
              ("Mean death LP, all deaths", "death_lp_all_deaths"),
              ("Mean death LP, completed cohort", "death_lp_completed_cohort"),
              ("Mean release LP, completed cohort", "release_lp_completed_cohort"),
              ("Mean horizon LP, censored cohort", "horizon_lp_censored_cohort"),
              ("DAMP injected at completed releases", "completed_injected_damp"),
              ("DAMP injected after last immune update", "terminal_injected_damp"),
              ("Field DAMP before terminal additions", "damp_before_terminal"),
              ("Field DAMP after terminal additions", "damp_after_terminal"),
              ("Peak field DAMP after terminal additions", "final_damp_peak"),
              ("Unique eligible cells", "unique_eligible_cells"),
              ("Eligible cell-step opportunities", "eligible_cell_steps"), ("Immune kills", "immune_kills"),
              ("Kills / unique eligible cells", "kills_per_unique_eligible_cell"),
              ("Kills / eligible cell-step opportunities", "kills_per_eligible_cell_step"),
              ("Mean local DAMP per eligible opportunity", "mean_damp_per_eligible_cell_step"),
              ("Maximum local DAMP among eligible opportunities", "max_eligible_damp"),
              ("Mean activation per eligible opportunity", "mean_activation"),
              ("Eligible opportunities with DAMP ≥ Kd", COUNTS[0]),
              ("Fraction of eligible opportunities with DAMP ≥ Kd", COUNTS[0] + "_fraction"),
              ("Eligible opportunities with DAMP ≥ 9Kd", COUNTS[1]),
              ("Fraction of eligible opportunities with DAMP ≥ 9Kd", COUNTS[1] + "_fraction")]
    command = "python3 scripts/immune_2d_measurement_report.py"
    if archive.resolve() != ARCHIVE.resolve() or report.resolve() != REPORT.resolve():
        command = shlex.join(["python3", str(ROOT / "scripts/immune_2d_measurement_report.py"),
                              "--archive", str(archive.resolve()), "--report", str(report.resolve())])
    lines = ["# Passive 2D immune measurements", "",
             "These observations reconstruct the historical model comparison; they are not independent biological validation.", "",
             "| Measurement | Control | RSL3 | SDT |", "|---|---:|---:|---:|"]
    lines += ["| " + label + " | " + " | ".join(f(s[key]) for s in summaries) + " |" for label, key in labels]
    lines += ["", "## What was observed", "",
              "The historical baseline uses a 500×500 grid, 20 µm cells, a 4,500 µm tumor radius, "
              "geometry seed 42, 180 steps, λ=120 µm oxygen gradient, and the existing O2-independent "
              "SDT setting. Immune coupling uses the 2D defaults; anti-PD-1, stromal shielding and pH "
              "modulation are off. Treatment-specific runtime seeds and all parameters are archived.", "",
              "Living tumor cells are eligible immediately before immune killing at step ≥60 and "
              "local DAMP ≥0.01. Cells killed by that update are included. Repeated eligibility forms "
              "cell-step opportunities; unique eligibility counts each cell once. Both kill rates are "
              "descriptive, not causal treatment effects or competing-risk corrections.", "",
              "Activation is DAMP / (DAMP + Kd), with Kd=50 model units. DAMP ≥Kd corresponds to "
              "activation ≥0.5; DAMP ≥9Kd corresponds to activation ≥0.9. The fractions use eligible "
              "cell-step opportunities, not all cells or terminal heatmap pixels. Empty means, maxima "
              "and fractions are undefined. Aggregates bound and reconcile these observations; they "
              "do not reconstruct every unarchived cell-step exposure.", "",
              "Completed DAMP releases occur five steps after death, before diffusion and immune "
              "killing. Releases scheduled at step 180 or later are censored; their horizon LP and "
              "terminal additions are separate. Terminal additions occur after all immune updates "
              "and cannot explain earlier immune kills. Field totals differ from cumulative injections "
              "because of clearance. LP and DAMP remain uncalibrated model quantities.", "",
              "## Interpretation and limits", "",
              f"The observed SDT:RSL3 total immune-kill ratio is {f(ratio(summaries[2]['immune_kills'], summaries[1]['immune_kills']))}:1. "
              "This reconstructs the counts behind the historical approximately 104:1 comparison. "
              "It does not measure DAMP potency, DC maturation or immunogenicity per dead cell. "
              "Release means refer to the completed-release cohort; deaths censored at the horizon "
              "have separate means and denominators.", "",
              "This is one geometry seed and one runtime seed per arm. No replicate uncertainty or "
              "experimental calibration is added. The eligible-exposure measurements allow the "
              "saturation explanation to be inspected, without establishing it as a causal explanation. "
              "The older final non-ferroptotic-pool normalized statistic is a different denominator. "
              "The frozen 3D report remains a separate model comparison with different geometry, "
              "RNG streams and transport settings.", "",
              "## Reproduction", "", "```bash", command, "```", "",
              "Reconstruction checks event accounting, eligibility and activation summaries, archive "
              "hashes, frozen configuration, historical rows and the unchanged full 33-condition "
              "production summary without building or running a simulation.", "",
              f"- [Protocol]({report_link(ROOT / PROTOCOL, report)}).",
              f"- [Archive and provenance]({report_link(archive / 'manifest.json', report)}).",
              f"- [Unchanged full baseline]({report_link(archive / 'baseline-summary.json', report)}).",
              f"- [Separate frozen 3D report]({report_link(common.REPORT, report)}).",
              f"- Source commit: `{manifest['source_commit']}`.",
              f"- Captured: `{manifest['captured_at_utc']}`.",
              f"- Platform: `{manifest['platform']}`; `{manifest['rustc']}`; execution: serial.",
              f"- Built binary SHA-256: `{manifest['binary_sha256']}`.",
              f"- Full baseline SHA-256: `{manifest['artifacts']['baseline-summary.json']}`.", ""]
    reference = manifest["baseline_reference"]
    lines += [f"- Uninstrumented baseline source commit: `{reference['source_commit']}`.",
              f"- Uninstrumented baseline binary SHA-256: `{reference['binary_sha256']}`.",
              f"- Uninstrumented baseline toolchain: `{reference['toolchain']}`.", ""]
    lines += [f"- `{s['condition']}` runtime seed: `{s['condition_seed']}`." for s in summaries]
    return "\n".join(lines) + "\n"


def write_archive(path: Path, raw: bytes, baseline: bytes, sources: dict[str, bytes], metadata: dict) -> None:
    path.mkdir()
    (path / "observations.json.gz").write_bytes(gzip.compress(raw, mtime=0))
    (path / "baseline-summary.json").write_bytes(baseline)
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as bundle:
        for name, content in sorted(sources.items()):
            info = tarfile.TarInfo(name)
            info.size, info.mode = len(content), 0o644
            bundle.addfile(info, io.BytesIO(content))
    (path / "sources.tar.gz").write_bytes(gzip.compress(buffer.getvalue(), mtime=0))
    manifest = dict(metadata, schema_version=SCHEMA, simulator="sim-tme", dimension=2,
                    baseline_reference=baseline_reference(sources[GOLDEN].decode()),
                    sources={name: sha(content) for name, content in sources.items()},
                    artifacts={name: sha((path / name).read_bytes()) for name in sorted(ARTIFACTS)})
    (path / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def capture(destination: Path) -> None:
    require(not destination.exists(), "capture directory already exists; archives are immutable")
    require(not any(k.startswith("FERRO_") for k in os.environ), "unset FERRO_* overrides before capture")
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True)
    require(not status.strip(), "commit the protocol and implementation before capture")
    sources = {name: (ROOT / name).read_bytes() for name in source_paths()}
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    sim = ROOT / "simulations"
    rustc = subprocess.check_output(["rustup", "run", TOOLCHAIN, "rustc", "--version"],
                                    cwd=sim, text=True).strip()
    binary = common.build_binary(["rustup", "run", TOOLCHAIN, "cargo"], sim, "sim-tme")
    with tempfile.TemporaryDirectory(prefix="immune-2d-observations-") as tmp:
        work = Path(tmp)
        private_binary = work / "sim-tme"
        shutil.copy2(binary, private_binary)
        binary_digest = sha(private_binary.read_bytes())
        baseline = b""
        for folder, args in (("baseline", []), ("measurement", ["--immune-measurements"])):
            cwd = work / folder
            cwd.mkdir()
            subprocess.run([str(private_binary), *args], cwd=cwd, check=True)
            if folder == "baseline":
                baseline = (cwd / "output/tme/tme_summary.json").read_bytes()
                require(sha(baseline) == golden_sha(sources[GOLDEN].decode()), "2D production SHA changed")
                validate_baseline(json.loads(baseline), json.loads(sources[HISTORICAL]))
        raw = (work / "measurement/output/tme/immune_measurements.json").read_bytes()
        validate_observations(json.loads(raw), json.loads(baseline))
        require(sources == {name: (ROOT / name).read_bytes() for name in source_paths()},
                "sources changed during capture")
        staged = work / "archive"
        write_archive(staged, raw, baseline, sources, {
            "source_commit": commit, "captured_at_utc": datetime.now(timezone.utc).isoformat(),
            "platform": platform.platform(), "rustc": rustc, "execution": "serial",
            "binary_sha256": binary_digest,
        })
        load_archive(staged)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(staged, destination)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--capture", action="store_true")
    args = parser.parse_args()
    validate_output_path(args.archive, args.report)
    if args.capture:
        capture(args.archive)
    manifest, summaries = load_archive(args.archive)
    text = render(manifest, summaries, args.archive, args.report)
    validate_output_path(args.archive, args.report)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(text)
    print(args.report)


if __name__ == "__main__":
    main()
