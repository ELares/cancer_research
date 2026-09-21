"""Independent event accounting for passive immune observations (no simulator)."""

import copy
import gzip
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
           "immune": {"damp_per_lp": 1, "damp_clearance_rate": .03},
           "params": {"death_threshold": 10}}
    events = [
        {"cell_index": 1, "death_step": 5, "death_lp": 11., "scheduled_release_step": 7,
         "release_step": 7, "release_lp": 14., "release_damp": 14.,
         "horizon_lp": None, "terminal_damp": None},
        {"cell_index": 2, "death_step": 6, "death_lp": 12., "scheduled_release_step": 8,
         "release_step": None, "release_lp": None, "release_damp": None,
         "horizon_lp": 13., "terminal_damp": 13.},
        # An early source makes the immune opportunities at steps 2--4
        # physically possible; the canonical field starts with no DAMP.
        {"cell_index": 3, "death_step": 0, "death_lp": 11., "scheduled_release_step": 2,
         "release_step": 2, "release_lp": 14., "release_damp": 14.,
         "horizon_lp": None, "terminal_damp": None},
    ]
    cells = [
        {"cell_index": 0, "first_step": 2, "last_step": 2, "opportunities": 1, "local_damp_sum": .02},
        {"cell_index": 1, "first_step": 2, "last_step": 4, "opportunities": 3, "local_damp_sum": .06},
    ]
    steps = [{"step": s, "ferroptotic_deaths": int(s in (0, 5, 6)),
              "completed_releases": int(s in (2, 7)), "released_damp": 14. if s in (2, 7) else 0.,
              "eligible_cells": 2 if s == 2 else int(s in (3, 4)),
              "eligible_local_damp_sum": .04 if s == 2 else .02 if s in (3, 4) else 0.,
              "immune_kills": int(s == 2)} for s in range(8)]
    preterminal = 14 * .97 ** 6 + 14 * .97
    row = {"condition_name": "immune_RSL3", "condition_seed": report.condition_seed("immune_RSL3"),
           "result": {"total_tumor": 5, "total_dead": 4, "ferroptosis_kills": 3,
                      "immune_kills": 1, "total_damp": preterminal + 13.},
           "measurements": {"ferroptotic_events": events, "eligible_cells": cells, "steps": steps,
                            "immune_kill_events": [{"cell_index": 0, "step": 2, "local_damp": .02}],
                            "terminal": {"censored_deaths": 1, "terminal_additions": 1,
                                         "terminal_damp": 13., "damp_before_terminal": preterminal,
                                         "damp_after_terminal": preterminal + 13.}}}
    return row, cfg


def test_completed_and_censored_cohorts_have_distinct_means_and_denominators():
    row, cfg = fixture()
    result = report.reconcile_condition(row, cfg)
    assert result["death_lp_all_deaths"] == 34 / 3
    assert result["death_lp_completed_cohort"] == 11
    assert result["release_lp_completed_cohort"] == 14
    assert result["horizon_lp_censored_cohort"] == 13
    assert result["unique_eligible_cells"] == 2
    assert result["eligible_cell_steps"] == 4
    assert result["kills_per_unique_eligible_cell"] == .5
    assert result["kills_per_eligible_cell_step"] == .25
    assert result["terminal_injected_damp"] == 13
    assert result["completed_injected_damp"] == 28


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
    "step_count", "nonfinite_lp", "threshold", "lp_falls", "endpoint_count", "shift_opportunity",
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
    elif mutation == "shift_opportunity":
        obs["steps"][3].update(eligible_cells=0, eligible_local_damp_sum=0)
        obs["steps"][6].update(eligible_cells=1, eligible_local_damp_sum=.02)
    with pytest.raises(ValueError):
        report.reconcile_condition(row, cfg)


def archive():
    assert report.ARCHIVE.is_dir(), "the frozen production archive must be committed"
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
    obs["ferroptotic_events"] = [obs["ferroptotic_events"][2]]
    obs["immune_kill_events"] = []
    obs["terminal"] = {k: 0 for k in obs["terminal"]}
    mass = 14 * .97 ** 6
    obs["terminal"].update(damp_before_terminal=mass, damp_after_terminal=mass)
    row["result"].update(total_tumor=2, total_dead=1, ferroptosis_kills=1,
                         immune_kills=0, total_damp=mass)
    obs["eligible_cells"].append(copy.deepcopy(obs["eligible_cells"][0]))
    for i, cell in enumerate(obs["eligible_cells"]):
        cell.update(cell_index=i, first_step=2+i, last_step=2+i, opportunities=1, local_damp_sum=.02)
    for s in obs["steps"]:
        s.update({k: 0 for k in s if k != "step"})
        if s["step"] in (2, 3, 4):
            s.update(eligible_cells=1, eligible_local_damp_sum=.02)
        if s["step"] == 0:
            s["ferroptotic_deaths"] = 1
        if s["step"] == 2:
            s.update(completed_releases=1, released_damp=14.)
    with pytest.raises(ValueError, match="unique eligible cells exceed"):
        report.reconcile_condition(row, cfg)


def test_historical_archive_does_not_require_current_production_golden(monkeypatch):
    target = archive()
    monkeypatch.setattr(report, "expected_sha", lambda: "0" * 64)
    manifest, summaries = report.load_archive(target)
    assert "0" * 64 not in report.render(manifest, summaries)


@pytest.mark.parametrize("mutation, message", [
    ("release_mass", "release/clearance DAMP mass balance"),
    ("cell_kill_damp", "kill DAMP exceeds cell opportunity sum"),
    ("step_kill_damp", "kill DAMP exceeds step opportunity sum"),
    ("eligible_above_field", "eligible DAMP exceeds whole field mass"),
])
def test_physically_inconsistent_damp_ledgers_are_rejected(mutation, message):
    row, cfg = fixture()
    obs = row["measurements"]
    if mutation == "release_mass":
        # Internal release totals still match; the unchanged final field
        # cannot be produced by these injections and the declared clearance.
        event = obs["ferroptotic_events"][0]
        event["release_lp"] += 100
        event["release_damp"] += 100
        obs["steps"][event["release_step"]]["released_damp"] += 100
    elif mutation == "cell_kill_damp":
        # Preserve the overall eligible sum but contradict the killed cell's
        # one observed opportunity (its kill event still records .02).
        obs["eligible_cells"][0]["local_damp_sum"] -= .01
        obs["eligible_cells"][1]["local_damp_sum"] += .01
    elif mutation == "step_kill_damp":
        # At step 2, the kill contributes .02 and the other eligible cell
        # must contribute at least .01. Both cells still reconcile globally.
        obs["steps"][2]["eligible_local_damp_sum"] -= .02
        obs["steps"][3]["eligible_local_damp_sum"] += .02
    else:
        obs["eligible_cells"][1]["local_damp_sum"] += 100
        obs["steps"][3]["eligible_local_damp_sum"] += 100
    with pytest.raises(ValueError, match=message):
        report.reconcile_condition(row, cfg)


@pytest.mark.parametrize("clearance", [0., 1.])
def test_release_mass_balance_handles_clearance_boundaries(clearance):
    row, cfg = fixture()
    cfg["immune"]["damp_clearance_rate"] = clearance
    obs = row["measurements"]
    if clearance == 1:
        # Immediate full clearance leaves no immune opportunities, even
        # though ferroptotic sources and a terminal addition still exist.
        obs["eligible_cells"] = []
        obs["immune_kill_events"] = []
        row["result"].update(immune_kills=0, total_dead=3)
        for step in obs["steps"]:
            step.update(eligible_cells=0, eligible_local_damp_sum=0., immune_kills=0)
    mass = 0.
    for step in obs["steps"]:
        mass = (mass + step["released_damp"]) * (1 - clearance)
    obs["terminal"].update(damp_before_terminal=mass, damp_after_terminal=mass + 13.)
    row["result"]["total_damp"] = mass + 13.
    summary = report.reconcile_condition(row, cfg)
    assert summary["damp_before_terminal"] == (28. if clearance == 0 else 0.)
    assert summary["terminal_injected_damp"] == 13.


def test_eligible_subset_may_equal_whole_field_mass():
    row, cfg = fixture()
    obs = row["measurements"]
    # One eligible cell at step 3 contains the whole current field; equality
    # is a valid subset boundary, not an overflow of the field mass.
    mass = 14 * .97 ** 2
    obs["steps"][3]["eligible_local_damp_sum"] = mass
    obs["eligible_cells"][1]["local_damp_sum"] = mass + .04
    report.reconcile_condition(row, cfg)


@pytest.fixture(scope="module")
def production_data():
    return (
        json.loads(gzip.decompress((archive() / "observations.json.gz").read_bytes())),
        json.loads((archive() / "baseline-summary.json").read_text()),
    )


@pytest.mark.parametrize("population", ["ferroptotic_events", "immune_kill_events", "eligible_cells"])
def test_stromal_id_cannot_enter_any_tumor_event_population(production_data, population):
    raw, baseline = production_data
    data = copy.deepcopy(raw)
    obs = data["conditions"][2]["measurements"]
    if population == "immune_kill_events":
        kill = obs[population][0]
        cell = next(c for c in obs["eligible_cells"] if c["cell_index"] == kill["cell_index"])
        cell["cell_index"] = kill["cell_index"] = 0
    elif population == "eligible_cells":
        dead = {e["cell_index"] for key in ("ferroptotic_events", "immune_kill_events") for e in obs[key]}
        next(c for c in obs[population] if c["cell_index"] not in dead)["cell_index"] = 0
    else:
        # SDT's first ferroptotic death precedes all immune opportunities.
        obs[population][0]["cell_index"] = 0
    with pytest.raises(ValueError, match="outside canonical tumor sphere"):
        report.validate_observations(data, baseline)


def test_canonical_tumor_sphere_includes_exact_radius(production_data):
    raw, _ = production_data
    cfg = raw["config"]
    # The center is (30, 30, 30) and radius is 27 cells. This also pins the
    # flat-index axis order used to recover the geometric membership.
    assert report.in_tumor_sphere(3 * 60 * 60 + 30 * 60 + 30, cfg)
    assert not report.in_tumor_sphere(2 * 60 * 60 + 30 * 60 + 30, cfg)
    assert report.in_tumor_sphere(30 * 60 * 60 + 30 * 60 + 30, cfg)
    assert not report.in_tumor_sphere(0, cfg)


@pytest.mark.parametrize("path, value", [
    (("biochem_seed_salt",), 0),
    (("immune_seed_salt",), 0),
    (("params", "sdt_ros"), 999999.),
    (("spatial_params", "neighbor_iron_fraction"), 0.),
])
def test_complete_configuration_contract_rejects_mislabeled_inputs(production_data, path, value):
    raw, baseline = production_data
    data = copy.deepcopy(raw)
    parent = data["config"]
    for key in path[:-1]:
        parent = parent[key]
    parent[path[-1]] = value
    with pytest.raises(ValueError, match="complete canonical configuration changed"):
        report.validate_observations(data, baseline)
