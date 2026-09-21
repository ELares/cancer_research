"""Independent event accounting for passive immune observations (no simulator)."""

import copy
import json
from pathlib import Path
import shutil
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import immune_measurement_report as report


def fixture():
    cfg = {"grid_dim": 3, "n_steps": 8, "immune_start_step": 2,
           "post_death_steps": 2, "damp_kill_threshold": .01,
           "immune": {"damp_per_lp": 1}, "params": {"death_threshold": 10}}
    events = [
        {"cell_index": 1, "death_step": 5, "death_lp": 11., "scheduled_release_step": 7,
         "release_step": 7, "release_lp": 14., "release_damp": 14.,
         "horizon_lp": None, "terminal_damp": None},
        {"cell_index": 2, "death_step": 6, "death_lp": 12., "scheduled_release_step": 8,
         "release_step": None, "release_lp": None, "release_damp": None,
         "horizon_lp": 13., "terminal_damp": 13.},
    ]
    cells = [
        {"cell_index": 0, "first_step": 2, "last_step": 2, "opportunities": 1, "local_damp_sum": .02},
        {"cell_index": 1, "first_step": 2, "last_step": 4, "opportunities": 3, "local_damp_sum": .06},
    ]
    steps = [{"step": s, "ferroptotic_deaths": int(s in (5, 6)),
              "completed_releases": int(s == 7), "released_damp": 14. if s == 7 else 0.,
              "eligible_cells": 2 if s == 2 else int(s in (3, 4)),
              "eligible_local_damp_sum": .04 if s == 2 else .02 if s in (3, 4) else 0.,
              "immune_kills": int(s == 2)} for s in range(8)]
    row = {"condition_name": "immune_RSL3", "condition_seed": report.condition_seed("immune_RSL3"),
           "result": {"total_tumor": 5, "total_dead": 3, "ferroptosis_kills": 2,
                      "immune_kills": 1, "total_damp": 20.},
           "measurements": {"ferroptotic_events": events, "eligible_cells": cells, "steps": steps,
                            "immune_kill_events": [{"cell_index": 0, "step": 2, "local_damp": .02}],
                            "terminal": {"censored_deaths": 1, "terminal_additions": 1,
                                         "terminal_damp": 13., "damp_before_terminal": 7.,
                                         "damp_after_terminal": 20.}}}
    return row, cfg


def test_completed_and_censored_cohorts_have_distinct_means_and_denominators():
    row, cfg = fixture()
    result = report.reconcile_condition(row, cfg)
    assert result["death_lp_all_deaths"] == 11.5
    assert result["death_lp_completed_cohort"] == 11
    assert result["release_lp_completed_cohort"] == 14
    assert result["horizon_lp_censored_cohort"] == 13
    assert result["unique_eligible_cells"] == 2
    assert result["eligible_cell_steps"] == 4
    assert result["kills_per_unique_eligible_cell"] == .5
    assert result["kills_per_eligible_cell_step"] == .25
    assert result["terminal_injected_damp"] == 13
    assert result["completed_injected_damp"] == 14


def test_empty_populations_produce_undefined_means_and_rates():
    row, cfg = fixture()
    obs = row["measurements"]
    for key in ("ferroptotic_events", "eligible_cells", "immune_kill_events"):
        obs[key] = []
    for step in obs["steps"]:
        step.update({k: 0 for k in step if k != "step"})
    obs["terminal"] = {k: 0 for k in obs["terminal"]}
    row["result"].update(total_dead=0, ferroptosis_kills=0, immune_kills=0, total_damp=0)
    result = report.reconcile_condition(row, cfg)
    for key in ("death_lp_all_deaths", "release_lp_completed_cohort", "horizon_lp_censored_cohort",
                "kills_per_unique_eligible_cell", "kills_per_eligible_cell_step"):
        assert result[key] is None
    assert result["ferroptotic_deaths"] == 0


@pytest.mark.parametrize("mutation", [
    "duplicate_death", "censor_as_complete", "wrong_release", "kill_before_delay",
    "dead_eligible", "unobserved_kill", "opportunity_total", "terminal_sum",
    "step_count", "nonfinite_lp", "threshold", "lp_falls", "endpoint_count",
])
def test_corrupt_event_or_denominator_is_rejected(mutation):
    row, cfg = fixture()
    obs = row["measurements"]
    if mutation == "duplicate_death":
        obs["ferroptotic_events"].append(copy.deepcopy(obs["ferroptotic_events"][0]))
    elif mutation == "censor_as_complete":
        obs["ferroptotic_events"][1]["release_step"] = 8
    elif mutation == "wrong_release":
        obs["ferroptotic_events"][0]["release_step"] = 6
    elif mutation == "kill_before_delay":
        obs["immune_kill_events"][0]["step"] = 1
    elif mutation == "dead_eligible":
        obs["eligible_cells"][1]["last_step"] = 5
    elif mutation == "unobserved_kill":
        obs["immune_kill_events"][0]["cell_index"] = 4
    elif mutation == "opportunity_total":
        obs["eligible_cells"][1]["opportunities"] = 2
    elif mutation == "terminal_sum":
        obs["terminal"]["terminal_damp"] = 14
    elif mutation == "step_count":
        obs["steps"][7]["completed_releases"] = 0
    elif mutation == "nonfinite_lp":
        obs["ferroptotic_events"][0]["death_lp"] = float("nan")
    elif mutation == "threshold":
        obs["immune_kill_events"][0]["local_damp"] = .009
    elif mutation == "lp_falls":
        obs["ferroptotic_events"][0]["release_lp"] = 9
    elif mutation == "endpoint_count":
        obs["eligible_cells"][1]["opportunities"] = 1
    with pytest.raises(ValueError):
        report.reconcile_condition(row, cfg)


def archive():
    if not report.ARCHIVE.exists():
        pytest.skip("production archive is captured after the implementation freeze")
    return report.ARCHIVE


def test_committed_archive_reconciles_and_report_is_reproducible():
    manifest, summaries = report.load_archive(archive())
    assert report.REPORT.read_text() == report.render(manifest, summaries)
    assert len(summaries) == 3
    text = report.render(manifest, summaries)
    assert "undefined (n=0)" in text
    assert "not causal treatment effects" in text
    assert "after every immune update" in text
    assert "104:1" in text and "2D" in text


@pytest.mark.parametrize("mutation", ["missing", "extra", "corrupt", "source"])
def test_archive_fails_closed_and_preserves_report(tmp_path, monkeypatch, mutation):
    target = tmp_path / "archive"
    shutil.copytree(archive(), target)
    if mutation == "missing":
        (target / "observations.json.gz").unlink()
    elif mutation == "extra":
        (target / "unrecorded.json").write_text("{}")
    elif mutation == "corrupt":
        with (target / "observations.json.gz").open("ab") as f:
            f.write(b"corruption")
    else:
        manifest = json.loads((target / "manifest.json").read_text())
        first = next(iter(manifest["sources"]))
        manifest["sources"][first] = "0" * 64
        (target / "manifest.json").write_text(json.dumps(manifest))
    output = tmp_path / "report.md"
    output.write_text("published finding\n")
    monkeypatch.setattr(sys, "argv", ["report", "--archive", str(target), "--report", str(output)])
    with pytest.raises((ValueError, FileNotFoundError)):
        report.main()
    assert output.read_text() == "published finding\n"


def test_capture_refuses_replacing_any_archive(tmp_path):
    with pytest.raises(ValueError, match="immutable"):
        report.capture(tmp_path, 8)


def test_unique_population_cannot_exceed_tumor_even_at_disjoint_steps():
    row, cfg = fixture()
    obs = row["measurements"]
    obs["ferroptotic_events"] = []
    obs["immune_kill_events"] = []
    obs["terminal"] = {k: 0 for k in obs["terminal"]}
    row["result"].update(total_tumor=1, total_dead=0, ferroptosis_kills=0, immune_kills=0, total_damp=0)
    for i, cell in enumerate(obs["eligible_cells"]):
        cell.update(first_step=2+i, last_step=2+i, opportunities=1, local_damp_sum=.02)
    for s in obs["steps"]:
        s.update({k: 0 for k in s if k != "step"})
        if s["step"] in (2, 3):
            s.update(eligible_cells=1, eligible_local_damp_sum=.02)
    with pytest.raises(ValueError, match="unique eligible cells exceed"):
        report.reconcile_condition(row, cfg)


def test_historical_archive_does_not_require_current_production_golden(monkeypatch):
    target = archive()
    monkeypatch.setattr(report, "expected_sha", lambda: "0" * 64)
    manifest, summaries = report.load_archive(target)
    assert "0" * 64 not in report.render(manifest, summaries)
