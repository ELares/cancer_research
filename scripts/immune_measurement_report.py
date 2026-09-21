#!/usr/bin/env python3
"""Capture or independently reconcile the canonical 3D immune observations.

Default: render the committed archive, without Rust or simulation execution.
Capture builds pinned, committed sources, checks the complete production SHA,
then observes the three existing immune conditions in a separate directory.
See docs/IMMUNE_MEASUREMENT_PROTOCOL.md for populations and interpretation.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
import tarfile

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "analysis" / "immune-measurements"
REPORT = ROOT / "analysis" / "immune-measurement-report.md"
NAMES = ("immune_Control", "immune_RSL3", "immune_SDT")
ARTIFACTS = {"observations.json.gz", "baseline-summary.json", "sources.tar.gz"}
TOOLCHAIN = "1.96.0"


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def finite(value, label: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool), label)
    require(math.isfinite(value) and value >= 0, label)
    return value


def integer(value, label: str) -> int:
    require(type(value) is int and value >= 0, label)
    return value


def close(a: float, b: float, label: str) -> None:
    finite(a, label)
    finite(b, label)
    require(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-8), label)


def mean(values: list[float]) -> float | None:
    return math.fsum(values) / len(values) if values else None


def ratio(a: float, b: float) -> float | None:
    return a / b if b else None


def condition_seed(name: str) -> int:
    h = 0xCBF29CE484222325
    for b in name.encode():
        h = ((h ^ b) * 0x100000001B3) & ((1 << 64) - 1)
    return (h + 42) & ((1 << 64) - 1)


def reconcile_condition(row: dict, cfg: dict) -> dict:
    """Recompute totals from events; reject internally inconsistent archives."""
    result, obs = row["result"], row["measurements"]
    n_steps, delay = cfg["n_steps"], cfg["immune_start_step"]
    grace, threshold = cfg["post_death_steps"], cfg["damp_kill_threshold"]
    factor = cfg["immune"]["damp_per_lp"]
    events, kills = obs["ferroptotic_events"], obs["immune_kill_events"]
    cells, steps, terminal = obs["eligible_cells"], obs["steps"], obs["terminal"]
    require([s["step"] for s in steps] == list(range(n_steps)), "complete ordered steps")

    def cell_id(value):
        integer(value, "cell index")
        require(value < cfg["grid_dim"] ** 3, "cell index outside grid")
        return value

    def index(rows):
        ids = [cell_id(r["cell_index"]) for r in rows]
        require(len(ids) == len(set(ids)), "duplicate cell events")
        return dict(zip(ids, rows))

    death_by_id, kill_by_id, eligible_by_id = index(events), index(kills), index(cells)
    require(not death_by_id.keys() & kill_by_id.keys(), "death causes overlap")
    completed, censored = [], []
    death_steps, release_steps, kill_steps = Counter(), Counter(), Counter()
    released = [[] for _ in steps]
    for e in events:
        ds = integer(e["death_step"], "death step")
        require(ds < n_steps, "death outside simulation")
        finite(e["death_lp"], "death LP")
        require(e["death_lp"] > cfg["params"]["death_threshold"], "death below LP threshold")
        require(e["scheduled_release_step"] == ds + grace, "release schedule")
        death_steps[ds] += 1
        if ds + grace < n_steps:
            require(e["release_step"] == ds + grace, "completed release step")
            require(e["horizon_lp"] is None and e["terminal_damp"] is None,
                    "completed event mislabeled censored")
            finite(e["release_lp"], "release LP")
            require(e["release_lp"] >= e["death_lp"], "LP fell during release grace")
            close(e["release_damp"], e["release_lp"] * factor, "release DAMP")
            release_steps[e["release_step"]] += 1
            released[e["release_step"]].append(e["release_damp"])
            completed.append(e)
        else:
            require(all(e[k] is None for k in ("release_step", "release_lp", "release_damp")),
                    "censored event labeled completed")
            finite(e["horizon_lp"], "horizon LP")
            require(e["horizon_lp"] >= e["death_lp"], "LP fell before horizon")
            close(e["terminal_damp"], e["horizon_lp"] * factor, "terminal event DAMP")
            censored.append(e)

    eligible_interval_changes = Counter()
    eligible_endpoints = Counter()
    for c in cells:
        first = integer(c["first_step"], "first eligible step")
        last = integer(c["last_step"], "last eligible step")
        count = integer(c["opportunities"], "cell opportunities")
        require(delay <= first <= last < n_steps, "eligibility timing")
        require(1 <= count <= last - first + 1, "cell opportunity count")
        require(first == last or count >= 2, "distinct eligibility endpoints need two opportunities")
        eligible_interval_changes[first] += 1
        eligible_interval_changes[last + 1] -= 1
        eligible_endpoints[first] += 1
        if last != first:
            eligible_endpoints[last] += 1
        require(finite(c["local_damp_sum"], "eligible DAMP") >= threshold * count - 1e-8,
                "eligible DAMP below threshold")
        if c["cell_index"] in death_by_id:
            require(last < death_by_id[c["cell_index"]]["death_step"],
                    "ferroptotic cell eligible after biochemistry death")
    for k in kills:
        step = integer(k["step"], "immune kill step")
        require(delay <= step < n_steps, "immune kill timing")
        require(finite(k["local_damp"], "kill DAMP") >= threshold, "kill below threshold")
        require(k["cell_index"] in eligible_by_id, "kill without eligible cell")
        c = eligible_by_id[k["cell_index"]]
        require(c["last_step"] == step, "immune cell eligible after kill")
        kill_steps[step] += 1

    dead_before_immune = possible_eligible = 0
    for s in steps:
        step = s["step"]
        possible_eligible += eligible_interval_changes[step]
        for field in ("ferroptotic_deaths", "completed_releases", "eligible_cells", "immune_kills"):
            integer(s[field], field)
        require(s["ferroptotic_deaths"] == death_steps[step], "step death reconciliation")
        require(s["completed_releases"] == release_steps[step], "step release reconciliation")
        require(s["immune_kills"] == kill_steps[step], "step kill reconciliation")
        close(s["released_damp"], math.fsum(released[step]), "step release DAMP")
        dead_before_immune += s["ferroptotic_deaths"]
        require(s["immune_kills"] <= s["eligible_cells"] <= result["total_tumor"] - dead_before_immune,
                "step eligibility count")
        require(eligible_endpoints[step] <= s["eligible_cells"] <= possible_eligible,
                "step eligibility outside cell observation intervals")
        dead_before_immune += s["immune_kills"]
        if step < delay:
            require(s["eligible_cells"] == 0, "eligibility before activation")
        finite(s["eligible_local_damp_sum"], "step eligible DAMP")
        require(s["eligible_cells"] > 0 or s["eligible_local_damp_sum"] == 0,
                "eligible DAMP without eligible cells")
        require(s["eligible_local_damp_sum"] >= threshold * s["eligible_cells"] - 1e-8,
                "step DAMP below threshold")

    opportunities = sum(c["opportunities"] for c in cells)
    require(opportunities == sum(s["eligible_cells"] for s in steps), "opportunity reconciliation")
    close(math.fsum(c["local_damp_sum"] for c in cells),
          math.fsum(s["eligible_local_damp_sum"] for s in steps), "eligible DAMP reconciliation")
    require(len(events) == result["ferroptosis_kills"], "ferroptotic census")
    require(len(kills) == result["immune_kills"], "immune census")
    require(len(events) + len(kills) == result["total_dead"] <= result["total_tumor"], "dead census")
    require(len(kills) <= len(cells) <= opportunities, "unique/opportunity denominators")
    require(len(cells) <= result["total_tumor"], "unique eligible cells exceed tumor population")
    require(len(censored) == terminal["censored_deaths"] == terminal["terminal_additions"],
            "terminal event reconciliation")
    terminal_damp = math.fsum(e["terminal_damp"] for e in censored)
    close(terminal_damp, terminal["terminal_damp"], "terminal DAMP sum")
    close(terminal["damp_before_terminal"] + terminal_damp, terminal["damp_after_terminal"],
          "terminal DAMP balance")
    close(terminal["damp_after_terminal"], result["total_damp"], "final DAMP census")
    return {
        "condition": row["condition_name"], "condition_seed": row["condition_seed"],
        "total_tumor": result["total_tumor"], "ferroptotic_deaths": len(events),
        "completed_releases": len(completed), "censored_deaths": len(censored),
        "death_lp_all_deaths": mean([e["death_lp"] for e in events]),
        "death_lp_completed_cohort": mean([e["death_lp"] for e in completed]),
        "release_lp_completed_cohort": mean([e["release_lp"] for e in completed]),
        "horizon_lp_censored_cohort": mean([e["horizon_lp"] for e in censored]),
        "completed_injected_damp": math.fsum(e["release_damp"] for e in completed),
        "terminal_injected_damp": terminal_damp,
        "damp_before_terminal": terminal["damp_before_terminal"],
        "damp_after_terminal": terminal["damp_after_terminal"],
        "unique_eligible_cells": len(cells), "eligible_cell_steps": opportunities,
        "immune_kills": len(kills), "kills_per_unique_eligible_cell": ratio(len(kills), len(cells)),
        "kills_per_eligible_cell_step": ratio(len(kills), opportunities),
        "mean_damp_per_eligible_cell_step": ratio(math.fsum(c["local_damp_sum"] for c in cells), opportunities),
    }


def expected_sha() -> str:
    text = (ROOT / "simulations/sim-tme-3d/expected_summary.sha256").read_text()
    return next(line.split()[0] for line in text.splitlines() if line and not line.startswith("#"))


def validate_observations(data: dict, baseline: dict) -> list[dict]:
    require(data["schema_version"] == 1, "measurement schema")
    cfg = data["config"]
    for k, v in {"grid_dim": 60, "cell_size_um": 20.0, "tumor_radius_um": 540.0,
                 "n_steps": 180, "grid_seed": 42, "post_death_steps": 5,
                 "immune_start_step": 60, "damp_kill_threshold": 0.01}.items():
        require(cfg[k] == v, f"noncanonical configuration: {k}")
    require([r["condition_name"] for r in data["conditions"]] == list(NAMES), "canonical conditions")
    require(cfg["immune"] == {
        "damp_per_lp": 1.0, "damp_diffusion_fraction": .025, "damp_clearance_rate": .03,
        "dc_activation_kd": 50.0, "immune_kill_rate": .02, "pd1_brake": .7,
        "anti_pd1_efficacy": 0.0, "exhaustion_rate": 0.0, "ferro_immunosuppression_strength": 0.0,
    }, "canonical immune configuration")
    require(len(baseline["conditions"]) == 24, "complete baseline matrix")
    summaries = []
    for row in data["conditions"]:
        require(row["condition_seed"] == condition_seed(row["condition_name"]), "canonical seed")
        require(row["configuration"] == {
            "treatment": row["result"]["treatment"], "o2_lambda_um": 120.0,
            "immune_on": True, "stromal_on": False, "ph_on": False, "dose_schedule": "Constant",
        }, "canonical condition settings")
        candidates = [r for r in baseline["conditions"]
                      if r["treatment"] == row["result"]["treatment"]
                      and r["immune_mode"] == "immune_on"
                      and "stromal_mode" not in r and "ph_mode" not in r]
        require(len(candidates) == 1 and candidates[0] == row["result"], "baseline row changed")
        require(row["condition_name"] == "immune_" + row["result"]["treatment"], "condition identity")
        summaries.append(reconcile_condition(row, cfg))
    return summaries


def source_paths() -> list[str]:
    tracked = subprocess.check_output(["git", "ls-files", "simulations"], cwd=ROOT, text=True).splitlines()
    return sorted([p for p in tracked if p.endswith((".rs", ".toml", ".lock"))] + [
        "scripts/immune_measurement_report.py", "docs/IMMUNE_MEASUREMENT_PROTOCOL.md",
        "simulations/sim-tme-3d/expected_summary.sha256",
    ])


def load_archive(path: Path) -> tuple[dict, list[dict]]:
    manifest = json.loads((path / "manifest.json").read_text())
    require(manifest["schema_version"] == 1, "archive schema")
    require(set(manifest["artifacts"]) == ARTIFACTS, "archive artifact set")
    require({p.name for p in path.iterdir()} == ARTIFACTS | {"manifest.json"}, "unexpected archive files")
    blobs = {name: (path / name).read_bytes() for name in ARTIFACTS}
    for name, blob in blobs.items():
        require(sha(blob) == manifest["artifacts"][name], f"artifact hash: {name}")
    # Verify the frozen implementation, not today's evolving Rust checkout.
    # Never extract archive paths; reconstruction also works in shallow clones.
    with tarfile.open(fileobj=io.BytesIO(blobs["sources.tar.gz"]), mode="r:gz") as bundle:
        members = bundle.getmembers()
        require(len(members) == len(manifest["sources"]), "source inventory size")
        require({m.name for m in members} == set(manifest["sources"]), "source inventory")
        for member in members:
            require(member.isfile() and not Path(member.name).is_absolute()
                    and ".." not in Path(member.name).parts, "invalid source archive entry")
            require(sha(bundle.extractfile(member).read()) == manifest["sources"][member.name],
                    f"frozen source hash: {member.name}")
        expected = bundle.extractfile("simulations/sim-tme-3d/expected_summary.sha256").read().decode()
        frozen_sha = next(line.split()[0] for line in expected.splitlines() if line and not line.startswith("#"))
        require(sha(blobs["baseline-summary.json"]) == frozen_sha, "frozen production SHA changed")
    data = json.loads(gzip.decompress(blobs["observations.json.gz"]))
    summaries = validate_observations(data, json.loads(blobs["baseline-summary.json"]))
    return manifest, summaries


def render(manifest: dict, summaries: list[dict]) -> str:
    def f(value):
        if value is None:
            return "undefined (n=0)"
        if isinstance(value, int):
            return f"{value:,}"
        return f"{value:.6g}"

    labels = [
        ("Tumor cells", "total_tumor"),
        ("Ferroptotic death events", "ferroptotic_deaths"),
        ("Completed releases", "completed_releases"),
        ("Deaths with release censored at horizon", "censored_deaths"),
        ("Mean death-time LP, all ferroptotic deaths", "death_lp_all_deaths"),
        ("Mean death-time LP, completed-release cohort", "death_lp_completed_cohort"),
        ("Mean release-time LP, completed-release cohort", "release_lp_completed_cohort"),
        ("Mean horizon LP, censored-death cohort", "horizon_lp_censored_cohort"),
        ("DAMP injected at completed releases", "completed_injected_damp"),
        ("DAMP injected after last immune update", "terminal_injected_damp"),
        ("DAMP field sum before terminal additions", "damp_before_terminal"),
        ("DAMP field sum after terminal additions", "damp_after_terminal"),
        ("Unique eligible living tumor cells", "unique_eligible_cells"),
        ("Eligible cell-step opportunities", "eligible_cell_steps"),
        ("Immune kills", "immune_kills"),
        ("Kills / unique eligible cells", "kills_per_unique_eligible_cell"),
        ("Kills / eligible cell-step opportunities", "kills_per_eligible_cell_step"),
        ("Mean local DAMP per eligible cell-step", "mean_damp_per_eligible_cell_step"),
    ]
    rsl3, sdt = summaries[1:]
    lines = ["# Passive 3D immune measurements", "",
             "These are observations of the existing model, not independent biological validation.", "",
             "| Measurement | Control | RSL3 | SDT |", "|---|---:|---:|---:|"]
    lines += ["| " + label + " | " + " | ".join(f(s[key]) for s in summaries) + " |" for label, key in labels]
    lines += ["", "## What this measures", "",
              "The three canonical conditions share a 60³ grid, 20 µm cells, a 540 µm tumor radius, "
              "an oxygen gradient with λ=120 µm, 180 simulation steps, and the default 3D immune configuration. "
              "Geometry seed 42 is shared; treatment-specific runtime seeds are preserved. "
              "Stromal shielding, pH modulation, repopulation and additional killing modalities are off.", "",
              "Eligibility is counted immediately before each immune update: a living tumor cell, "
              "step ≥60, and local DAMP ≥0.01. Unique cells count the union of those cells; "
              "cell-step opportunities count repeated eligibility. A cell killed in that update is included. "
              "Both rates are descriptive conditional summaries, not causal treatment effects or a "
              "correction for competing ferroptotic deaths. Different cells and exposure durations contribute across arms.", "",
              "LP is in model units. DAMP is an uncalibrated model field with one injected unit per LP unit. "
              "Completed release occurs five steps after ferroptotic death, before diffusion and immune killing. "
              "A release scheduled at or beyond step 180 is right-censored; its horizon LP and terminal "
              "DAMP addition are reported separately. Terminal additions occur after every immune update "
              "and cannot explain kills already counted. Field sums differ from cumulative injections because of clearance.", "",
              "## Interpretation and limits", "",
              f"The observed 3D SDT:RSL3 total immune-kill ratio is {f(ratio(sdt['immune_kills'], rsl3['immune_kills']))}. "
              "This does not measure DAMP potency per dead cell. The completed-release means compare a "
              "named event cohort; they exclude horizon-censored releases and should not be substituted "
              "for all-death means.", "",
              "This is one existing geometry seed and one existing runtime seed per condition. There are "
              "no replicate uncertainty intervals, no calibrated conversion from cell-steps to hours, "
              "and no fitted per-cell experimental release threshold. The historical 104:1 figure belongs "
              "to the 2D model; this 3D observation neither reconstructs that experiment nor validates its mechanism. "
              "The older `immune_kills / max(total_tumor - ferroptosis_kills, 1)` statistic uses a final "
              "non-ferroptotic pool, not the observed eligible population.", "",
              "All death, release, censoring, eligibility and kill counts reconcile from the archived records. "
              "Every measured legacy result equals its corresponding default-matrix row; the complete "
              f"24-condition summary retains SHA-256 `{manifest['artifacts']['baseline-summary.json']}`.", "",
              "## Reproduction", "",
              "```bash", "python3 scripts/immune_measurement_report.py", "```", "",
              "This validates compressed raw event records, per-cell eligibility counts, every step, "
              "artifact hashes and frozen source hashes, then regenerates this report without running a simulation. "
              "See the [protocol](../docs/IMMUNE_MEASUREMENT_PROTOCOL.md) for capture and validation commands.", "",
              f"- Frozen source commit: `{manifest['source_commit']}`.",
              f"- Captured: `{manifest['captured_at_utc']}`.",
              f"- Platform: `{manifest['platform']}`; `{manifest['rustc']}`; Rayon threads: {manifest['rayon_threads']}.",
              f"- Built binary SHA-256: `{manifest['binary_sha256']}`.",
              "- [Archive and provenance](immune-measurements/manifest.json).",
              "- [Unchanged production summary](immune-measurements/baseline-summary.json).", ""]
    lines += [f"- `{s['condition']}` runtime seed: `{s['condition_seed']}`." for s in summaries]
    return "\n".join(lines) + "\n"


def capture(destination: Path, threads: int) -> None:
    require(threads > 0, "threads must be positive")
    require(not destination.exists(), "capture directory already exists; archives are immutable")
    require(not any(k.startswith("FERRO_") for k in os.environ), "unset FERRO_* overrides before capture")
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True)
    require(not status.strip(), "commit the protocol and implementation before capture")
    sources = {p: sha((ROOT / p).read_bytes()) for p in source_paths()}
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    sim = ROOT / "simulations"
    cargo = ["rustup", "run", TOOLCHAIN, "cargo"]
    rustc = subprocess.check_output(["rustup", "run", TOOLCHAIN, "rustc", "--version"], cwd=sim, text=True).strip()
    subprocess.run([*cargo, "build", "--locked", "--release", "-p", "sim-tme-3d"], cwd=sim, check=True)
    meta = json.loads(subprocess.check_output([*cargo, "metadata", "--format-version=1", "--no-deps"], cwd=sim))
    binary = Path(meta["target_directory"]) / "release" / "sim-tme-3d"
    env = dict(os.environ, RAYON_NUM_THREADS=str(threads))
    with tempfile.TemporaryDirectory(prefix="immune-observations-") as tmp:
        work = Path(tmp)
        # A shared target directory may be used by other worktrees. Run and
        # identify one private copy, not a path a concurrent build can replace.
        private_binary = work / "sim-tme-3d"
        shutil.copy2(binary, private_binary)
        binary_digest = sha(private_binary.read_bytes())
        for folder, args in (("baseline", []), ("measurement", ["--immune-measurements"])):
            cwd = work / folder
            cwd.mkdir()
            subprocess.run([str(private_binary), *args], cwd=cwd, env=env, check=True)
        baseline = (work / "baseline/output/tme-3d/summary.json").read_bytes()
        raw = (work / "measurement/output/tme-3d/immune_measurements.json").read_bytes()
        require(sha(baseline) == expected_sha(), "production SHA changed")
        validate_observations(json.loads(raw), json.loads(baseline))
        require(sources == {p: sha((ROOT / p).read_bytes()) for p in source_paths()}, "sources changed during capture")
        staged = work / "archive"
        staged.mkdir()
        (staged / "observations.json.gz").write_bytes(gzip.compress(raw, mtime=0))
        (staged / "baseline-summary.json").write_bytes(baseline)
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w") as bundle:
            for name in sources:
                content = (ROOT / name).read_bytes()
                info = tarfile.TarInfo(name)
                info.size = len(content)
                info.mode = 0o644
                bundle.addfile(info, io.BytesIO(content))
        (staged / "sources.tar.gz").write_bytes(gzip.compress(buffer.getvalue(), mtime=0))
        manifest = {
            "schema_version": 1, "source_commit": commit,
            "captured_at_utc": datetime.now(timezone.utc).isoformat(),
            "platform": platform.platform(), "rayon_threads": threads,
            "rustc": rustc,
            "binary_sha256": binary_digest, "sources": sources,
            "artifacts": {p: sha((staged / p).read_bytes()) for p in sorted(ARTIFACTS)},
        }
        (staged / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        load_archive(staged)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(staged, destination)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--capture", action="store_true")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.capture:
        capture(args.archive, args.threads)
    manifest, summaries = load_archive(args.archive)
    text = render(manifest, summaries)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(text)
    print(args.report)


if __name__ == "__main__":
    main()
