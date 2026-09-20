"""Fixed-target decisions, archive integrity and independent-run reporting gates.

The simulator is replaced with controlled curves. These tests exercise the
driver's accounting and decisions without requiring the compiled extension.
"""

import copy
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import abc_joint_importance as driver  # noqa: E402


@pytest.fixture(autouse=True)
def no_compiled_extension(monkeypatch):
    def forbidden():
        pytest.fail("driver unit tests must not load the compiled extension")

    monkeypatch.setattr(driver.abc.ck, "_fc", forbidden)


@pytest.fixture
def target():
    rd, ed = list(driver.abc.ck.DOSE_GRID_UM), list(driver.abc.ce.DOSE_GRID_UM)
    return {"priors": [list(p) for p in driver.abc.PRIORS],
            "rsl3_doses_um": rd, "erastin_doses_um": ed,
            "empirical_rsl3": [0.0] * len(rd),
            "empirical_erastin": [0.0] * len(ed),
            "empirical_heldout": [0.0] * len(rd),
            "epsilon": 0.125,
            "reference_simulator_dose_calls": len(rd) + len(ed)}


def install_models(monkeypatch, rsl3_offset=0.0, erastin_offset=0.0):
    calls = []

    def model(name, offset):
        def run(doses, params):
            calls.append((name, list(doses), dict(params)))
            return [offset] * len(doses)
        return run

    monkeypatch.setattr(driver.abc, "model_rsl3", model("rsl3", rsl3_offset))
    monkeypatch.setattr(driver.abc, "model_erastin", model("erastin", erastin_offset))
    return calls


@pytest.mark.parametrize("rsl3,erastin", [(0.0, 0.125), (0.0625, 0.0625), (0.125, 0.0)])
@pytest.mark.parametrize("direction", [-1, 0, 1])
def test_early_rejection_preserves_exact_and_neighboring_threshold_decisions(
        monkeypatch, target, rsl3, erastin, direction):
    calls = install_models(monkeypatch, rsl3, erastin)
    distance = rsl3 + erastin
    target["epsilon"] = (math.nextafter(distance, -math.inf) if direction < 0 else
                         math.nextafter(distance, math.inf) if direction > 0 else distance)
    full = driver.evaluate_unit([0.5] * 7, target)
    calls.clear()
    screened = driver.evaluate_unit([0.5] * 7, target, early_reject=True)
    assert full["distance"] == distance
    assert screened["accepted"] == full["accepted"] == (direction >= 0)
    if erastin > target["epsilon"]:
        assert screened["status"] == "early_rejected"
        assert screened["distance"] is None and screened["rsl3"] is None
        assert screened["distance_lower_bound"] == erastin
        assert screened["simulator_dose_calls"] == len(target["erastin_doses_um"])
        assert [c[0] for c in calls] == ["erastin"]
    else:
        assert screened == full
        assert [c[0] for c in calls] == ["erastin", "rsl3"]


def test_pilot_default_completes_scores_even_when_erastin_alone_exceeds_epsilon(monkeypatch, target):
    calls = install_models(monkeypatch, 0.0625, 0.25)
    result = driver.evaluate_unit([0.5] * 7, target)
    assert not result["accepted"] and result["status"] == "complete"
    assert result["distance"] == 0.3125
    assert result["simulator_dose_calls"] == 13
    assert [c[0] for c in calls] == ["erastin", "rsl3"]


@pytest.mark.parametrize("early_reject", [False, True])
def test_heldout_targets_do_not_change_training_distance_or_acceptance(monkeypatch, target, early_reject):
    install_models(monkeypatch, 0.03125, 0.0625)
    before = driver.evaluate_unit([0.5] * 7, target, early_reject=early_reject)
    changed = copy.deepcopy(target)
    changed["empirical_heldout"] = [1e6] * len(target["empirical_heldout"])
    after = driver.evaluate_unit([0.5] * 7, changed, early_reject=early_reject)
    assert after == before


@pytest.mark.parametrize("coordinate", range(7))
@pytest.mark.parametrize("value", [math.nextafter(0.0, -math.inf), math.nextafter(1.0, math.inf)])
def test_every_outside_coordinate_is_rejected_without_simulator_calls(monkeypatch, target, coordinate, value):
    calls = install_models(monkeypatch)
    point = [0.5] * 7
    point[coordinate] = value
    result = driver.evaluate_unit(point, target, early_reject=True)
    assert result["status"] == "outside_prior" and not result["accepted"]
    assert result["simulator_dose_calls"] == 0 and calls == []
    assert result["distance"] is None and result["distance_lower_bound"] is None
    assert result["rsl3"] is None and result["erastin"] is None


@pytest.mark.parametrize("endpoint", [0.0, 1.0])
def test_cube_endpoints_are_included_and_map_to_the_documented_prior(monkeypatch, target, endpoint):
    calls = install_models(monkeypatch)
    result = driver.evaluate_unit([endpoint] * 7, target)
    assert result["accepted"] and result["status"] == "complete"
    expected = {name: lo if endpoint == 0 else hi for name, lo, hi in target["priors"]}
    for _, _, params in calls:
        assert params == pytest.approx(expected, rel=1e-15, abs=1e-15)


@pytest.mark.parametrize("point", [[0.5] * 6, [0.5] * 8, [math.nan] + [0.5] * 6,
                                   [math.inf] + [0.5] * 6])
def test_malformed_unit_points_fail_before_simulation(monkeypatch, target, point):
    calls = install_models(monkeypatch)
    with pytest.raises(ValueError):
        driver.evaluate_unit(point, target, early_reject=True)
    assert calls == []


@pytest.mark.parametrize("compound", ["rsl3", "erastin"])
@pytest.mark.parametrize("bad", [math.nan, math.inf])
@pytest.mark.parametrize("early_reject", [False, True])
def test_nonfinite_model_curves_are_errors_not_valid_rejections(monkeypatch, target, compound, bad, early_reject):
    install_models(monkeypatch)
    monkeypatch.setattr(driver.abc, f"model_{compound}", lambda doses, params: [bad] * len(doses))
    with pytest.raises(ValueError):
        driver.evaluate_unit([0.5] * 7, target, early_reject=early_reject)


@pytest.fixture
def archive(target):
    points = [[u] * 7 for u in (0.2, 0.7, 0.4, 0.8, 0.5, 0.3)]
    points[4][0] = -0.1
    proposal = driver.GaussianMixture(7, 0.2, [[0.5] * 7], [np.eye(7) * 0.1])

    def complete(r, e):
        rsl3, erastin = [r] * len(target["empirical_rsl3"]), [e] * len(target["empirical_erastin"])
        distance = driver.abc.ck.rmse(rsl3, target["empirical_rsl3"]) + driver.abc.ck.rmse(erastin, target["empirical_erastin"])
        accepted = distance <= target["epsilon"]
        return {"status": "complete", "distance": distance, "distance_lower_bound": distance,
                "rsl3": rsl3 if accepted else None, "erastin": erastin if accepted else None,
                "simulator_dose_calls": 13, "accepted": accepted}

    records = [complete(0.03125, 0.03125), complete(0.0625, 0.03125), complete(0.125, 0.125),
               {"status": "early_rejected", "distance": None, "distance_lower_bound": 0.25,
                "rsl3": None, "erastin": None, "simulator_dose_calls": 6, "accepted": False},
               {"status": "outside_prior", "distance": None, "distance_lower_bound": None,
                "rsl3": None, "erastin": None, "simulator_dose_calls": 0, "accepted": False},
               complete(0.0, 0.0)]
    return {"schema_version": 1, "target": target, "seed": driver.SEEDS[0], "proposal": proposal.to_dict(),
            "plan": dict(driver.PLAN, production_attempts=len(points)),
            "production": {"unit_points": points, "log_q": proposal.log_density(points).tolist(),
                           "component_ids": [1] * len(points), "records": records},
            "pilot": [{"attempts": 4}], "pilot_simulator_dose_calls": 12, "wall_seconds": 1.0}


def test_summary_reconstructs_weights_counts_and_costs_from_all_attempts(archive):
    result = driver.summarize_run(archive)
    assert result["importance"]["n_attempts"] == 6
    assert result["importance"]["n_accepted"] == 3
    assert result["production_outside_prior"] == 1
    assert result["production_early_rejections"] == 1
    assert result["simulator_dose_calls_including_pilot_reference"] == 83
    raw = [math.exp(-q) if r["accepted"] else 0.0
           for q, r in zip(archive["production"]["log_q"], archive["production"]["records"])]
    assert result["importance"]["normalizer_estimate"] == pytest.approx(sum(raw) / 6)
    assert result["importance"]["ess"] == pytest.approx(sum(raw) ** 2 / sum(w * w for w in raw))
    assert len(result["heldout"]["band_median"]) == len(archive["target"]["rsl3_doses_um"])
    assert not result["adequacy"]["passed"]


def test_summary_uses_archived_prior_bounds_not_current_module_globals(archive, monkeypatch):
    before = driver.summarize_run(archive)
    monkeypatch.setattr(driver, "LOWS", driver.LOWS + 100)
    monkeypatch.setattr(driver, "WIDTHS", driver.WIDTHS * 3)
    after = driver.summarize_run(archive)
    assert after["diagnostic_quantiles"] == before["diagnostic_quantiles"]
    assert after["prior_widths"] == before["prior_widths"]


def test_zero_acceptance_has_no_quantiles_or_heldout_inference(archive):
    for record in archive["production"]["records"]:
        if record["accepted"]:
            record.update(accepted=False, distance=0.5, distance_lower_bound=0.5,
                          rsl3=None, erastin=None)
    result = driver.summarize_run(archive)
    assert result["importance"]["n_accepted"] == 0
    assert result["importance"]["normalizer_relative_mcse"] is None
    assert result["diagnostic_quantiles"] is None and result["heldout"] is None
    assert not result["adequacy"]["passed"]


def test_summary_rejects_tampered_proposal_density(archive):
    archive["production"]["log_q"][0] += 0.1
    with pytest.raises(ValueError, match="densit"):
        driver.summarize_run(archive)


def test_summary_rejects_attempt_count_drift(archive):
    archive["plan"]["production_attempts"] += 1
    with pytest.raises(ValueError, match="count"):
        driver.summarize_run(archive)


@pytest.mark.parametrize("index,changes", [
    (2, {"accepted": True}),
    (4, {"accepted": True}),
    (4, {"simulator_dose_calls": 1}),
    (3, {"distance_lower_bound": 0.125}),
    (3, {"distance_lower_bound": math.inf}),
    (3, {"distance_lower_bound": math.nan}),
    (3, {"simulator_dose_calls": 13}),
    (0, {"simulator_dose_calls": 6}),
    (0, {"distance_lower_bound": -0.1}),
    (0, {"status": "unknown"}),
])
def test_summary_rejects_inconsistent_decisions_or_accounting(archive, index, changes):
    archive["production"]["records"][index].update(changes)
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


@pytest.mark.parametrize("compound", ["rsl3", "erastin"])
@pytest.mark.parametrize("malformation", ["truncated", "empty", "nonfinite", "wrong_score"])
def test_accepted_curve_integrity_is_checked_before_bands_or_broadcasting(archive, compound, malformation):
    for record in archive["production"]["records"]:
        if not record["accepted"]:
            continue
        if malformation == "truncated":
            record[compound] = record[compound][:1]
        elif malformation == "empty":
            record[compound] = []
        elif malformation == "nonfinite":
            record[compound][0] = math.inf
        else:
            record[compound][0] += 0.25
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


@pytest.mark.parametrize("key", ["empirical_rsl3", "empirical_erastin", "empirical_heldout"])
def test_summary_rejects_mismatched_target_lengths(archive, key):
    archive["target"][key] = archive["target"][key][:1]
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


@pytest.fixture
def adequate_summary():
    importance = driver.importance_diagnostics(np.zeros(512), np.arange(512) < 256)
    result = {"seed": driver.SEEDS[0], "importance": importance,
              "diagnostic_quantiles": {
                  name: {"q2_5": lo + 0.2 * (hi - lo), "median": lo + 0.5 * (hi - lo),
                         "q97_5": lo + 0.8 * (hi - lo)} for name, lo, hi in driver.abc.PRIORS},
              "prior_widths": {name: hi - lo for name, lo, hi in driver.abc.PRIORS},
              "heldout": {"rmse_quantiles": [0.01, 0.02, 0.03]}}
    result["adequacy"] = driver.run_adequacy(result)
    return result


@pytest.mark.parametrize("field,threshold,direction,gate", [
    ("n_accepted", "minimum_accepted", -1, "accepted_count"),
    ("ess", "minimum_ess", -1, "ess"),
    ("max_normalized_weight", "maximum_normalized_weight", 1, "maximum_weight"),
    ("normalizer_relative_mcse", "maximum_normalizer_relative_mcse", 1, "normalizer_mcse"),
])
def test_each_per_run_gate_is_inclusive_at_threshold_and_rejects_just_outside(
        adequate_summary, field, threshold, direction, gate):
    value = driver.GATES[threshold]
    adequate_summary["importance"][field] = value
    assert driver.run_adequacy(adequate_summary)["passed"]
    adequate_summary["importance"][field] = (value - 1 if field == "n_accepted" else
        math.nextafter(value, math.inf if direction > 0 else -math.inf))
    result = driver.run_adequacy(adequate_summary)
    assert not result["passed"] and not result["checks"][gate]


@pytest.mark.parametrize("field", ["max_normalized_weight", "normalizer_relative_mcse"])
def test_unavailable_weight_diagnostics_fail_adequacy(adequate_summary, field):
    adequate_summary["importance"][field] = None
    assert not driver.run_adequacy(adequate_summary)["passed"]


def planned_runs(summary):
    return [dict(copy.deepcopy(summary), seed=seed) for seed in driver.SEEDS]


def test_three_stable_prespecified_runs_pass_regardless_of_list_order(adequate_summary):
    runs = planned_runs(adequate_summary)
    assert driver.compare_runs(runs)["passed"]
    assert driver.compare_runs(list(reversed(runs)))["passed"]


@pytest.mark.parametrize("case", ["empty", "missing", "duplicate", "unplanned", "extra"])
def test_missing_duplicate_or_unplanned_seeds_cannot_establish_independence(adequate_summary, case):
    runs = planned_runs(adequate_summary)
    if case == "empty":
        runs = []
    elif case == "missing":
        runs.pop()
    elif case == "duplicate":
        runs[-1]["seed"] = runs[0]["seed"]
    elif case == "unplanned":
        runs[-1]["seed"] = max(driver.SEEDS) + 1
    else:
        runs.append(copy.deepcopy(runs[0]))
    result = driver.compare_runs(runs)
    assert not result["passed"]
    assert not result["checks"]["three_prespecified_independent_runs"]


@pytest.mark.parametrize("gate", ["all_runs_adequate", "median_stability", "tail_stability",
                                   "heldout_rmse_stability", "normalizer_agreement"])
def test_each_independent_run_screen_can_withhold_inference(adequate_summary, gate):
    runs = planned_runs(adequate_summary)
    changed = runs[-1]
    name = driver.NAMES[0]
    width = changed["prior_widths"][name]
    if gate == "all_runs_adequate":
        changed["importance"]["ess"] = driver.GATES["minimum_ess"] - 1
        changed["adequacy"] = driver.run_adequacy(changed)
    elif gate == "median_stability":
        changed["diagnostic_quantiles"][name]["median"] += 2 * driver.GATES["maximum_median_span_prior_fraction"] * width
    elif gate == "tail_stability":
        changed["diagnostic_quantiles"][name]["q97_5"] += 2 * driver.GATES["maximum_tail_span_prior_fraction"] * width
    elif gate == "heldout_rmse_stability":
        changed["heldout"]["rmse_quantiles"][1] += 2 * driver.GATES["maximum_heldout_rmse_median_span"]
    else:
        changed["importance"]["normalizer_estimate"] += 0.5
    result = driver.compare_runs(runs)
    assert not result["passed"] and not result["checks"][gate]


def test_incompatible_archived_prior_widths_cannot_be_compared(adequate_summary):
    runs = planned_runs(adequate_summary)
    runs[-1]["prior_widths"][driver.NAMES[0]] *= 2
    result = driver.compare_runs(runs)
    assert not result["passed"]


@pytest.mark.parametrize("field", ["diagnostic_quantiles", "heldout"])
def test_missing_run_summaries_do_not_pass_stability(adequate_summary, field):
    runs = planned_runs(adequate_summary)
    runs[0][field] = None
    result = driver.compare_runs(runs)
    assert not result["passed"]
    assert not result["checks"]["median_stability"]


def test_recomputed_accepted_curves_must_meet_the_exact_threshold(archive):
    record = archive["production"]["records"][0]
    epsilon = archive["target"]["epsilon"]
    record["rsl3"] = [0.0] * len(archive["target"]["rsl3_doses_um"])
    record["erastin"] = [math.nextafter(epsilon, math.inf)] * len(archive["target"]["erastin_doses_um"])
    record["distance"] = record["distance_lower_bound"] = epsilon
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


@pytest.fixture
def small_experiment(monkeypatch, target):
    plan = dict(driver.PLAN, pilot_rounds=2, pilot_attempts=32, production_attempts=64,
                min_elite=8, elite_fraction=0.25)
    target["epsilon"] = 0.5

    def model(coordinate):
        name, lo, hi = driver.abc.PRIORS[coordinate]

        def run(doses, params):
            return [(params[name] - lo) / (4 * (hi - lo))] * len(doses)
        return run

    monkeypatch.setattr(driver.abc, "model_rsl3", model(0))
    monkeypatch.setattr(driver.abc, "model_erastin", model(1))
    calls = []

    def evaluator(point, target, early_reject=False):
        calls.append((point.copy(), early_reject))
        return driver.evaluate_unit(point, target, early_reject=early_reject)

    result = driver.run_experiment(driver.SEEDS[0], target, plan=plan, evaluator=evaluator)
    return result, calls


def test_frozen_production_draws_reproduce_the_independent_spawned_rng(small_experiment):
    archive, calls = small_experiment
    pilot_stream, production_stream = np.random.SeedSequence(archive["seed"]).spawn(2)
    proposal = driver.GaussianMixture.from_dict(archive["proposal"])
    assert len(proposal.means) == archive["plan"]["pilot_rounds"]
    points, ids = proposal.sample(np.random.default_rng(production_stream),
                                  archive["plan"]["production_attempts"])
    assert np.array_equal(points, archive["production"]["unit_points"])
    assert np.array_equal(ids, archive["production"]["component_ids"])
    assert np.array_equal(proposal.log_density(points), archive["production"]["log_q"])
    first_pilot, _ = driver.GaussianMixture(7).sample(np.random.default_rng(pilot_stream), 32)
    assert np.array_equal(first_pilot, [p for p, _ in calls[:32]])
    assert archive["rng"]["pilot_spawn_key"] != archive["rng"]["production_spawn_key"]


def test_pilot_attempts_are_excluded_from_production_estimates(small_experiment):
    archive, calls = small_experiment
    summary = driver.summarize_run(archive)
    assert len(calls) == 128
    assert all(not early for _, early in calls[:64])
    assert all(early for _, early in calls[64:])
    assert summary["pilot_attempts"] == 64
    assert summary["importance"]["n_attempts"] == 64
    assert len(archive["production"]["records"]) == 64


def test_heldout_changes_cannot_train_proposals_or_change_production_decisions(small_experiment):
    archive, _ = small_experiment
    changed = copy.deepcopy(archive["target"])
    changed["empirical_heldout"] = [1.0] * len(changed["empirical_heldout"])
    rerun = driver.run_experiment(archive["seed"], changed, plan=archive["plan"])
    assert rerun["proposal"] == archive["proposal"]
    assert rerun["pilot"] == archive["pilot"]
    assert rerun["production"] == archive["production"]
    assert driver.summarize_run(rerun)["heldout"] != driver.summarize_run(archive)["heldout"]


def test_gzip_archive_roundtrip_preserves_rederived_summaries(small_experiment, tmp_path):
    archive, _ = small_experiment
    path = tmp_path / "run.json.gz"
    driver.write_archive(path, archive)
    restored = driver.read_archive(path)
    assert restored == archive
    before = driver.summarize_run(archive)
    assert before["importance"]["n_accepted"] > 0
    assert driver.summarize_run(restored) == before
    first_bytes = path.read_bytes()
    driver.write_archive(path, restored)
    assert path.read_bytes() == first_bytes


@pytest.fixture
def report_files(archive, tmp_path, monkeypatch):
    archive["target"].update(source={"sha256": "fixture-input"}, support={"fixture": "cohort"},
                             tolerance_factor=1.1, reference_distance=archive["target"]["epsilon"] / 1.1)
    archive.update(gates=copy.deepcopy(driver.GATES), source_hashes={"driver.py": "fixture-source"},
                   runtime={"extension_binaries": {"core.so": "fixture-binary"},
                            "python": "fixture-python", "numpy": "fixture-numpy", "scipy": "fixture-scipy"})
    monkeypatch.setattr(driver, "PLAN", archive["plan"])
    monkeypatch.setattr(driver, "OUT", tmp_path)
    for seed in driver.SEEDS:
        item = copy.deepcopy(archive)
        item["seed"] = seed
        driver.write_archive(driver.archive_path(tmp_path, seed), item)
    target = archive["target"]
    baseline = {"target_source": target["source"], "target_support": target["support"],
                "curves": {"rsl3_doses_um": target["rsl3_doses_um"],
                           "erastin_doses_um": target["erastin_doses_um"]},
                "tolerance_factor": target["tolerance_factor"],
                "epsilon_joint_distance": round(target["epsilon"], 4),
                "n_draws": 40000, "n_accepted": 4, "underpowered": True}
    (tmp_path / "joint-posterior.json").write_text(json.dumps(baseline))
    return tmp_path


def test_report_rebuild_records_the_matching_baseline_and_archive_hashes(report_files):
    report = driver.build_report(report_files)
    assert not report["stability"]["passed"]
    assert report["baseline"]["sha256"] == hashlib.sha256(
        (report_files / "joint-posterior.json").read_bytes()).hexdigest()
    assert len(report["archive_sha256"]) == len(driver.SEEDS)
    for name, digest in report["archive_sha256"].items():
        assert hashlib.sha256((report_files / name).read_bytes()).hexdigest() == digest


@pytest.mark.parametrize("field", ["target", "source_hashes", "extension_binaries", "python", "numpy", "scipy"])
def test_report_rejects_mixed_targets_source_or_numerical_runtimes(report_files, field):
    path = driver.archive_path(report_files, driver.SEEDS[-1])
    archive = driver.read_archive(path)
    if field == "target":
        archive["target"]["epsilon"] += 0.01
    elif field == "source_hashes":
        archive["source_hashes"]["driver.py"] = "different-source"
    elif field == "extension_binaries":
        archive["runtime"][field]["core.so"] = "different-binary"
    else:
        archive["runtime"][field] = "different-runtime"
    driver.write_archive(path, archive)
    with pytest.raises(ValueError):
        driver.build_report(report_files)


@pytest.mark.parametrize("field", ["target_source", "target_support", "rsl3_doses_um",
                                   "erastin_doses_um", "tolerance_factor", "epsilon_joint_distance"])
def test_report_rejects_baseline_metadata_from_a_different_target(report_files, field):
    path = report_files / "joint-posterior.json"
    baseline = json.loads(path.read_text())
    if field in {"target_source", "target_support"}:
        baseline[field] = {"different": "input"}
    elif field.endswith("doses_um"):
        baseline["curves"][field][0] *= 2
    else:
        baseline[field] += 0.01
    path.write_text(json.dumps(baseline))
    with pytest.raises(ValueError, match="baseline"):
        driver.build_report(report_files)
