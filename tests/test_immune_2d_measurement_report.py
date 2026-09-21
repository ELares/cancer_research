"""Independent 2D accounting, frozen provenance and passive capture checks."""

import copy
import json
from pathlib import Path
import re
import shlex
import sys
from urllib.parse import unquote

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import immune_2d_measurement_report as report


def population(values, kd=50):
    return {"activation_sum": sum(v / (v + kd) for v in values),
            "max_local_damp": max(values, default=None),
            "damp_ge_kd_opportunities": sum(v >= kd for v in values),
            "damp_ge_9kd_opportunities": sum(v >= 9 * kd for v in values)}


def fixture():
    cfg = json.loads(report.CONFIG.read_text())["config"]
    cfg.update(grid_rows=10, grid_cols=10, tumor_radius_um=80., n_steps=8, immune_start_step=2)
    cfg["params"]["post_death_steps"] = 2
    events = [
        {"cell_index": 56, "death_step": 5, "death_lp": 11., "scheduled_release_step": 7,
         "release_step": 7, "release_lp": 14., "release_damp": 14.,
         "horizon_lp": None, "terminal_damp": None},
        {"cell_index": 57, "death_step": 6, "death_lp": 12., "scheduled_release_step": 8,
         "release_step": None, "release_lp": None, "release_damp": None,
         "horizon_lp": 13., "terminal_damp": 13.},
        {"cell_index": 54, "death_step": 0, "death_lp": 11., "scheduled_release_step": 2,
         "release_step": 2, "release_lp": 1200., "release_damp": 1200.,
         "horizon_lp": None, "terminal_damp": None},
    ]
    cells = [dict(cell_index=55, first_step=2, last_step=2, opportunities=1,
                  local_damp_sum=450., **population([450.])),
             dict(cell_index=56, first_step=2, last_step=4, opportunities=3,
                  local_damp_sum=150.02, **population([50., .02, 100.]))]
    exposure = {2: [450., 50.], 3: [.02], 4: [100.]}
    steps = [dict(step=s, ferroptotic_deaths=int(s in (0, 5, 6)),
                  completed_releases=int(s in (2, 7)), released_damp=1200. if s == 2 else 14. if s == 7 else 0.,
                  eligible_cells=len(exposure.get(s, [])),
                  eligible_local_damp_sum=sum(exposure.get(s, [])), immune_kills=int(s == 2),
                  **population(exposure.get(s, []))) for s in range(8)]
    mass = 1200 * .97 ** 6 + 14 * .97
    row = {"condition_name": "immune_RSL3", "seed": 10000042,
           "result": {"treatment": "RSL3", "o2_condition": "gradient_120um", "o2_lambda_um": 120.,
                      "immune_mode": "immune_on", "stromal_mode": "off", "total_tumor": 49,
                      "total_dead": 4, "ferroptosis_kills": 3, "immune_kills": 1},
           "final_damp": {"total": mass + 13., "peak": 500.},
           "measurements": {"ferroptotic_events": events, "eligible_cells": cells, "steps": steps,
                            "immune_kill_events": [{"cell_index": 55, "step": 2, "local_damp": 450.}],
                            "terminal": {"censored_deaths": 1, "terminal_additions": 1,
                                         "terminal_damp": 13., "damp_before_terminal": mass,
                                         "damp_after_terminal": mass + 13.}}}
    return row, cfg


def dataset():
    row, cfg = fixture()
    rows = []
    for index, name in enumerate(report.NAMES):
        current = copy.deepcopy(row)
        current.update(condition_name=name, seed=42 + index * 10000000)
        current["result"]["treatment"] = name.removeprefix("immune_")
        rows.append(current)
    data = dict(schema_version=2, simulator="sim-tme", dimension=2, config=cfg, conditions=rows)
    historic = dict(schema_version=1, conditions=[copy.deepcopy(r["result"]) for r in rows])
    baseline = copy.deepcopy(historic)
    baseline["conditions"] += [{"treatment": "Control", "o2_condition": f"other_{i}",
                                "immune_mode": "off"} for i in range(30)]
    contract = dict(schema_version=2, config=cfg)
    return data, baseline, contract, historic


def source_fixture():
    data, baseline, contract, historical = dataset()
    raw = json.dumps(data).encode()
    baseline_bytes = json.dumps(baseline).encode()
    golden = ("# Source commit: " + "1" * 40 + "\n# Binary SHA-256: " + "2" * 64
              + "\n# Toolchain: rustc 1.96.0; serial execution.\n"
              + report.sha(baseline_bytes) + "  tme_summary.json\n")
    sources = {report.GOLDEN: golden.encode(), report.HISTORICAL: json.dumps(historical).encode(),
               str(report.CONFIG.relative_to(report.ROOT)): json.dumps(contract).encode(),
               report.PROTOCOL: b"synthetic protocol\n",
               "scripts/immune_2d_measurement_report.py": b"synthetic reader\n",
               "scripts/immune_measurement_report.py": b"synthetic shared reader\n"}
    return raw, baseline_bytes, sources


def archived(tmp_path):
    path = tmp_path / "archive #1 (frozen)"
    report.write_archive(path, *source_fixture(), {
        "source_commit": "3" * 40, "captured_at_utc": "2026-01-01T00:00:00+00:00",
        "platform": "synthetic", "rustc": "rustc 1.96.0", "execution": "serial", "binary_sha256": "4" * 64})
    return path


def test_release_cohorts_and_exposure_populations_remain_distinct():
    row, cfg = fixture()
    before = copy.deepcopy(row)
    result = report.reconcile_condition(row, cfg)
    assert row == before and "total_damp" not in row["result"] and "condition_seed" not in row
    assert result["ferroptotic_deaths"] == 3
    assert result["death_lp_all_deaths"] == 34 / 3
    assert result["release_lp_completed_cohort"] == 607
    assert result["horizon_lp_censored_cohort"] == 13
    assert result["kills_per_unique_eligible_cell"] == .5
    assert result["kills_per_eligible_cell_step"] == .25
    assert result["max_eligible_damp"] == 450
    assert result["damp_ge_kd_opportunities"] == 3
    assert result["damp_ge_9kd_opportunities"] == 1
    assert result["damp_ge_kd_opportunities_fraction"] == .75
    assert result["damp_ge_9kd_opportunities_fraction"] == .25
    assert result["mean_activation"] == pytest.approx((.9 + .5 + .02 / 50.02 + 2/3) / 4)


def test_empty_populations_have_undefined_means_maxima_and_fractions():
    row, cfg = fixture()
    obs = row["measurements"]
    for key in ("ferroptotic_events", "eligible_cells", "immune_kill_events"):
        obs[key] = []
    for step in obs["steps"]:
        step.update({key: 0 for key in step if key != "step"})
        step.update(population([]))
    obs["terminal"] = {key: 0 for key in obs["terminal"]}
    row["result"].update(total_dead=0, ferroptosis_kills=0, immune_kills=0)
    row["final_damp"] = dict(total=0, peak=0)
    result = report.reconcile_condition(row, cfg)
    for key in ("mean_activation", "max_eligible_damp", "damp_ge_kd_opportunities_fraction",
                "damp_ge_9kd_opportunities_fraction", "death_lp_all_deaths",
                "release_lp_completed_cohort", "horizon_lp_censored_cohort"):
        assert result[key] is None


@pytest.mark.parametrize("values", [[.01], [50.], [450.], [.01, 50., 450.], []])
def test_activation_threshold_boundaries_are_inclusive(values):
    _, cfg = fixture()
    report.check_activation(population(values), len(values), sum(values), cfg)


def test_known_maximum_tightens_activation_beyond_unconditional_jensen_bound():
    _, cfg = fixture()
    values = [99.99, .01]
    stats = population(values)
    report.check_activation(stats, 2, 100., cfg)
    stats["activation_sum"] = .95
    with pytest.raises(ValueError, match="attained-maximum"):
        report.check_activation(stats, 2, 100., cfg)


@pytest.mark.parametrize("mutation", [
    "activation_nan", "activation_high", "activation_low", "activation_cell_step",
    "threshold_nested", "threshold_maximum", "threshold_cell_step", "maximum_sum",
    "maximum_cell_step", "interval_maximum", "empty_maximum", "final_peak", "terminal_peak", "final_mass", "outside_grid",
    "after_death", "censoring", "missing_kill", "kill_maximum",
])
def test_corrupt_physical_or_activation_records_are_rejected(mutation):
    row, cfg = fixture()
    obs = row["measurements"]
    cell, step = obs["eligible_cells"][1], obs["steps"][2]
    if mutation == "activation_nan":
        cell["activation_sum"] = float("nan")
    elif mutation == "activation_high":
        cell["activation_sum"] = 3.
    elif mutation == "activation_low":
        cell["activation_sum"] = 0.
    elif mutation == "activation_cell_step":
        cell["activation_sum"] += .001
    elif mutation == "threshold_nested":
        cell[report.COUNTS[1]] = 3
    elif mutation == "threshold_maximum":
        cell[report.COUNTS[1]] = 1
    elif mutation == "threshold_cell_step":
        cell[report.COUNTS[0]] = 1
    elif mutation == "maximum_sum":
        cell["max_local_damp"] = 1000.
    elif mutation == "maximum_cell_step":
        step["max_local_damp"] = 451.
    elif mutation == "interval_maximum":
        cell["max_local_damp"] = 80.
    elif mutation == "empty_maximum":
        obs["steps"][0]["max_local_damp"] = 0.
    elif mutation == "final_peak":
        row["final_damp"]["peak"] = 0.
    elif mutation == "terminal_peak":
        row["final_damp"]["peak"] = row["final_damp"]["total"] / (cfg["grid_rows"] * cfg["grid_cols"])
    elif mutation == "final_mass":
        row["final_damp"]["total"] += 1.
    elif mutation == "outside_grid":
        obs["ferroptotic_events"][1]["cell_index"] = 100
    elif mutation == "after_death":
        cell["last_step"] = 5
    elif mutation == "censoring":
        obs["ferroptotic_events"][1]["release_step"] = 8
    elif mutation == "missing_kill":
        obs["immune_kill_events"] = []
    elif mutation == "kill_maximum":
        obs["immune_kill_events"][0]["local_damp"] = 451.
    with pytest.raises(ValueError):
        report.reconcile_condition(row, cfg)


@pytest.mark.parametrize("mutation", ["schema", "seed", "parameter", "rng_stride", "diffusion",
                                     "scenario", "duplicate_scenario", "historical", "circle"])
def test_contract_scenario_seed_and_tumor_identity_are_enforced(mutation):
    data, baseline, contract, historical = dataset()
    # Keep an independent pre-mutation contract.
    contract = copy.deepcopy(contract)
    if mutation == "schema":
        data["dimension"] = 3
    elif mutation == "seed":
        data["conditions"][1]["seed"] += 1
    elif mutation == "parameter":
        data["config"]["params"]["sdt_ros"] += 1
    elif mutation == "rng_stride":
        data["config"]["rng"]["immune_step_stride"] += 1
    elif mutation == "diffusion":
        data["config"]["immune_config"]["damp_diffusion_fraction"] = .025
    elif mutation == "scenario":
        baseline["conditions"][1]["ph_mode"] = "off"
    elif mutation == "duplicate_scenario":
        baseline["conditions"][-1] = copy.deepcopy(baseline["conditions"][0])
    elif mutation == "historical":
        historical["conditions"][0]["immune_kills"] += 1
    elif mutation == "circle":
        data["conditions"][0]["measurements"]["ferroptotic_events"][1]["cell_index"] = 0
    with pytest.raises(ValueError):
        report.validate_observations(data, baseline, contract=contract, historical=historical)


def test_complete_synthetic_dataset_reconciles():
    data, baseline, contract, historical = dataset()
    summaries = report.validate_observations(data, baseline, contract=contract, historical=historical)
    assert [s["condition_seed"] for s in summaries] == [42, 10000042, 20000042]


def test_frozen_archive_reconstructs_without_live_golden_or_contract(tmp_path, monkeypatch):
    archive = archived(tmp_path)
    manifest, summaries = report.load_archive(archive)
    # The same relative source names survive a different/no-longer-present checkout.
    original_root = report.ROOT
    new_root = tmp_path / "absent-checkout"
    monkeypatch.setattr(report, "ROOT", new_root)
    monkeypatch.setattr(report, "CONFIG", new_root / report.CONFIG.relative_to(original_root))
    assert report.load_archive(archive) == (manifest, summaries)


@pytest.mark.parametrize("mutation", ["artifact", "source_hash", "baseline_hash", "provenance", "extra"])
def test_invalid_archive_preserves_published_report(tmp_path, monkeypatch, mutation):
    archive = archived(tmp_path)
    manifest_path = archive / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if mutation == "artifact":
        with (archive / "observations.json.gz").open("ab") as stream:
            stream.write(b"corrupt")
    elif mutation == "source_hash":
        manifest["sources"][report.GOLDEN] = "0" * 64
    elif mutation == "baseline_hash":
        blob = (archive / "baseline-summary.json").read_bytes() + b"\n"
        (archive / "baseline-summary.json").write_bytes(blob)
        manifest["artifacts"]["baseline-summary.json"] = report.sha(blob)
    elif mutation == "provenance":
        manifest["baseline_reference"]["source_commit"] = "0" * 40
    elif mutation == "extra":
        (archive / "report.md").write_text("unrecorded")
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / "published.md"
    output.write_text("preserve this report\n")
    monkeypatch.setattr(sys, "argv", ["report", "--archive", str(archive), "--report", str(output)])
    with pytest.raises(ValueError):
        report.main()
    assert output.read_text() == "preserve this report\n"


@pytest.mark.parametrize("alias", ["child", "symlink", "hardlink"])
def test_output_aliases_cannot_modify_frozen_archive(tmp_path, monkeypatch, alias):
    archive = archived(tmp_path)
    before = {p.name: p.read_bytes() for p in archive.iterdir()}
    output = archive / "report.md" if alias == "child" else tmp_path / "report.md"
    if alias == "symlink":
        output.symlink_to(archive / "manifest.json")
    elif alias == "hardlink":
        output.hardlink_to(archive / "baseline-summary.json")
    monkeypatch.setattr(sys, "argv", ["report", "--capture", "--archive", str(archive), "--report", str(output)])
    monkeypatch.setattr(report, "capture", lambda *args: pytest.fail("capture started before path validation"))
    with pytest.raises(ValueError, match="immutable archive"):
        report.main()
    assert {p.name: p.read_bytes() for p in archive.iterdir()} == before


def test_custom_paths_links_and_reproduction_command(tmp_path, monkeypatch):
    archive = archived(tmp_path)
    output = tmp_path / "reports" / "one.md"
    nested = tmp_path / "nested"
    nested.mkdir()
    alias = tmp_path / "links" / "alias"
    alias.parent.mkdir()
    alias.symlink_to(nested, target_is_directory=True)
    monkeypatch.setattr(sys, "argv", ["report", "--archive", str(alias / ".." / archive.name),
                                     "--report", str(alias / ".." / "reports" / "one.md")])
    report.main()
    text = output.read_text()
    links = dict(re.findall(r"\[([^\]]+)\]\(([^)]+)\)", text))
    assert (output.parent / unquote(links["Archive and provenance"])).resolve() == archive.resolve() / "manifest.json"
    assert (output.parent / unquote(links["Protocol"])).resolve() == report.ROOT / report.PROTOCOL
    command = shlex.split(re.search(r"```bash\n([^\n]+)\n```", text).group(1))
    assert command[2:] == ["--archive", str(archive.resolve()), "--report", str(output.resolve())]
    monkeypatch.chdir(report.ROOT)
    monkeypatch.setattr(sys, "argv", command[1:])
    report.main()
    assert output.read_text() == text
    assert "not causal treatment effects" in text and "experimental calibration" in text
    report.load_archive(archive)


def test_committed_archive_and_report_are_required_and_reproducible():
    assert report.ARCHIVE.is_dir(), "the frozen 2D production archive must be committed"
    manifest, summaries = report.load_archive(report.ARCHIVE)
    assert report.REPORT.read_text() == report.render(manifest, summaries)
    assert [s["immune_kills"] for s in summaries] == [0, 5, 521]
    assert [s["ferroptotic_deaths"] for s in summaries] == [1, 163, 139640]
    assert "104.2:1" in report.REPORT.read_text()


@pytest.mark.parametrize("failure", [None, "baseline", "historical", "observations", "source_change"])
def test_capture_uses_private_binary_and_publishes_only_validated_data(tmp_path, monkeypatch, failure):
    raw, baseline, sources = source_fixture()
    if failure == "historical":
        historic = json.loads(sources[report.HISTORICAL])
        historic["conditions"][0]["immune_kills"] += 1
        sources[report.HISTORICAL] = json.dumps(historic).encode()
    original_root = report.ROOT
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    for name, content in sources.items():
        path = checkout / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    config_relative = report.CONFIG.relative_to(original_root)
    monkeypatch.setattr(report, "ROOT", checkout)
    monkeypatch.setattr(report, "CONFIG", checkout / config_relative)
    monkeypatch.setattr(report, "source_paths", lambda: sorted(sources))
    binary = tmp_path / "shared-cargo-target" / "sim-tme"
    binary.parent.mkdir()
    binary.write_bytes(b"compiled executable")
    monkeypatch.setattr(report.common, "build_binary", lambda cargo, sim, name: binary)

    def command_output(command, **kwargs):
        if command[:2] == ["git", "status"]:
            return ""
        return "3" * 40 if command[0] == "git" else "rustc 1.96.0"

    monkeypatch.setattr(report.subprocess, "check_output", command_output)
    monkeypatch.setattr(report.platform, "platform", lambda: "synthetic test platform")
    calls = []

    def run(command, *, cwd, check):
        calls.append(command)
        assert Path(command[0]) != binary
        assert Path(command[0]).read_bytes() == b"compiled executable"
        directory = cwd / "output/tme"
        directory.mkdir(parents=True)
        if len(command) == 1:
            (directory / "tme_summary.json").write_bytes(baseline + (b"\n" if failure == "baseline" else b""))
            # A later build changing the shared path cannot affect either execution.
            binary.write_bytes(b"unrelated later build")
        else:
            assert command[1:] == ["--immune-measurements"]
            value = json.loads(raw)
            if failure == "observations":
                value["conditions"][1]["seed"] += 1
            (directory / "immune_measurements.json").write_text(json.dumps(value))
            if failure == "source_change":
                (checkout / report.PROTOCOL).write_text("changed during capture")

    monkeypatch.setattr(report.subprocess, "run", run)
    for key in list(report.os.environ):
        if key.startswith("FERRO_"):
            monkeypatch.delenv(key)
    destination = tmp_path / "published-archive"
    if failure:
        with pytest.raises(ValueError):
            report.capture(destination)
        assert not destination.exists()
        assert len(calls) == (1 if failure in ("baseline", "historical") else 2)
    else:
        report.capture(destination)
        manifest, summaries = report.load_archive(destination)
        assert manifest["binary_sha256"] == report.sha(b"compiled executable")
        assert manifest["execution"] == "serial"
        assert len(summaries) == 3 and len(calls) == 2


def test_capture_refuses_existing_archive_and_overrides(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="immutable"):
        report.capture(tmp_path)
    monkeypatch.setenv("FERRO_SEED", "42")
    with pytest.raises(ValueError, match="FERRO_"):
        report.capture(tmp_path / "new")
