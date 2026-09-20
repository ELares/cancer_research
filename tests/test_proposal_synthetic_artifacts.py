"""Keep the failed coverage study visible when later methods are added."""

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import proposal_synthetic_validation as historical


def study():
    return json.loads(historical.OUT_JSON.read_text())


def test_failed_prerequisite_retains_every_declared_run():
    result = study()
    assert not result["passed"]
    assert {(r["fixture"], r["seed"]) for r in result["positive_runs"]} == {
        (fixture, seed) for fixture in historical.FIXTURES for seed in historical.POSITIVE_SEEDS
    }
    for run in result["positive_runs"]:
        assert run["assessment"]["passed"] == (run["fixture"] == "boundary_slab")
        assert run["pilot"]["all_islands_reached_epsilon"]


def test_apparently_adequate_weights_can_hide_a_missing_region():
    run = next(r for r in study()["positive_runs"]
               if r["fixture"] == "separated_boxes" and r["seed"] == 2026092102)
    result = run["assessment"]
    assert all(result["usual_checks"].values())
    assert not any(result["truth_checks"].values())
    assert result["mode_counts"]["A"] == 0
    assert result["importance"]["ess"] > 1000
    assert result["relative_mass_error"] < -0.28


def test_minority_neighborhoods_in_failed_runs_cross_the_gap():
    for run in study()["positive_runs"]:
        if run["fixture"] != "separated_boxes":
            continue
        points = np.unique(run["pilot"]["final_points"], axis=0)
        in_a = historical.membership(points, "separated_boxes")["A"]
        assert 0 < in_a.sum() < historical.PILOT_PLAN["kernel_neighbors"] < len(points)


def test_missing_region_controls_expose_weight_diagnostic_limit():
    controls = study()["negative_controls"]
    assert len(controls) == len(historical.NEGATIVE_SEEDS)
    assert all(not r["assessment"]["passed"] for r in controls)
    assert sum(all(r["assessment"]["usual_checks"].values()) for r in controls) == 2
    assert all(not r["assessment"]["truth_checks"]["mode_masses"] for r in controls)
