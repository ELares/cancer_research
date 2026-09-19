"""Fixed-target driver, archive replay and provenance checks on small fake runs."""

import copy
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import abc_joint_resample as driver  # noqa: E402


@pytest.fixture(autouse=True)
def no_compiled_extension(monkeypatch):
    def forbidden():
        pytest.fail("driver tests must not load the biological simulator")

    monkeypatch.setattr(driver.abc.ck, "_fc", forbidden)


@pytest.fixture
def settings():
    return {"pilot": {**driver.PLAN["pilot"], "particles_per_island": 16,
                       "levels": 3, "moves_per_level": 2, "kernel_neighbors": 4},
            "production_attempts": 96}


@pytest.fixture
def target():
    rd, ed = list(driver.abc.ck.DOSE_GRID_UM), list(driver.abc.ce.DOSE_GRID_UM)
    return {"priors": [list(p) for p in driver.abc.PRIORS],
            "rsl3_doses_um": rd, "erastin_doses_um": ed,
            "empirical_rsl3": [0.0] * len(rd), "empirical_erastin": [0.0] * len(ed),
            "empirical_heldout": [0.1] * len(rd), "epsilon": 0.2,
            "reference_distance": 0.2 / 1.1, "tolerance_factor": 1.1,
            "reference_simulator_dose_calls": len(rd) + len(ed),
            "source": {"sha256": "fixture-source"}, "support": {"cohort": "fixture"}}


@pytest.fixture
def models(monkeypatch):
    calls = []
    priors = {name: (lo, hi) for name, lo, hi in driver.abc.PRIORS}

    def model(compound, parameter, factor):
        def evaluate(doses, params):
            calls.append(compound)
            lo, hi = priors[parameter]
            value = float(factor * (params[parameter] - lo) / (hi - lo))
            return [value] * len(doses)
        return evaluate

    monkeypatch.setattr(driver.abc, "model_rsl3", model("rsl3", "lp_propagation", 0.25))
    monkeypatch.setattr(driver.abc, "model_erastin", model("erastin", "k_erastin", 1.0))
    return calls


@pytest.fixture
def archive(target, settings, models):
    return driver.run_experiment(driver.SEEDS[0], target, settings)


def test_small_run_reconciles_pilot_production_and_actual_dose_call_cost(archive, models):
    summary = driver.summarize_run(archive)
    records = archive["production"]["records"]
    assert any(r["accepted"] for r in records)
    assert any(r["status"] == "early_rejected" for r in records)
    assert any(r["status"] == "complete" and not r["accepted"] for r in records)
    assert all(r["simulator_dose_calls"] == (6 if r["status"] == "early_rejected" else 13)
               for r in records)
    expected = 7 * models.count("rsl3") + 6 * models.count("erastin") + 13
    assert summary["simulator_dose_calls_including_pilot_reference"] == expected
    assert summary["pilot_attempts"] == 2 * 16 * (1 + 3 * 2)
    assert summary["importance"]["n_attempts"] == 96
    assert summary["production_outside_prior"] == 0
    assert summary["importance"]["n_accepted"] == sum(r["accepted"] for r in records)
    raw = np.array([np.exp(-q) if r["accepted"] else 0 for q, r in
                    zip(archive["production"]["log_q"], records)])
    assert summary["importance"]["normalizer_estimate"] == pytest.approx(raw.mean())
    assert summary["importance"]["ess"] == pytest.approx(raw.sum()**2 / np.dot(raw, raw))
    assert not summary["adequacy"]["passed"]  # Tiny test budgets do not relax screens.
    driver.replay_production(archive)


def test_heldout_only_changes_neither_pilot_nor_production(archive, settings, models):
    changed = copy.deepcopy(archive["target"])
    changed["empirical_heldout"] = [0.9] * len(changed["empirical_heldout"])
    rerun = driver.run_experiment(archive["seed"], changed, settings)
    assert rerun["pilot"] == archive["pilot"]
    assert rerun["proposal"] == archive["proposal"]
    assert rerun["production"] == archive["production"]
    assert driver.summarize_run(rerun)["heldout"] != driver.summarize_run(archive)["heldout"]


def test_zero_production_acceptance_withholds_all_quantiles(archive):
    for r in archive["production"]["records"]:
        if r["accepted"]:
            r.update(accepted=False, distance=1.0, distance_lower_bound=1.0,
                     rsl3=[0.0] * len(archive["target"]["rsl3_doses_um"]),
                     erastin=[1.0] * len(archive["target"]["erastin_doses_um"]))
    summary = driver.summarize_run(archive)
    assert summary["importance"]["ess"] == 0
    assert summary["diagnostic_quantiles"] is None and summary["heldout"] is None
    assert not summary["adequacy"]["passed"]


def test_archive_roundtrip_reconstructs_the_same_statistics_and_stream(archive, tmp_path):
    path = driver.archive_path(tmp_path, archive["seed"])
    driver.write_archive(path, archive)
    restored = driver.read_archive(path)
    assert restored == archive
    assert driver.summarize_run(restored) == driver.summarize_run(archive)
    driver.replay_production(restored)
    original_bytes = path.read_bytes()
    driver.write_archive(path, restored)
    assert path.read_bytes() == original_bytes


@pytest.mark.parametrize("field", ["seed", "component", "point", "pilot_stream", "production_stream", "generator"])
def test_independent_production_replay_rejects_stream_tampering(archive, field):
    if field == "seed":
        archive["seed"] += 1
    elif field == "component":
        archive["production"]["component_ids"][0] += 1
    elif field == "point":
        archive["production"]["unit_points"][0][0] += 0.001
    elif field == "pilot_stream":
        archive["rng"]["pilot_spawn_key"] = [1]
    elif field == "production_stream":
        archive["rng"]["production_spawn_key"] = [0]
    else:
        archive["rng"]["bit_generator"] = "wrong-generator"
    with pytest.raises(ValueError):
        driver.replay_production(archive)


@pytest.mark.parametrize("change", ["attempts", "records", "density", "accepted_type", "early_cost", "early_bound"])
def test_production_count_density_and_decision_tampering_is_rejected(archive, change):
    p = archive["production"]
    if change == "attempts":
        archive["plan"]["production_attempts"] += 1
    elif change == "records":
        p["records"].pop()
    elif change == "density":
        p["log_q"][0] += 0.01
    elif change == "accepted_type":
        p["records"][0]["accepted"] = 0
    else:
        early = next(r for r in p["records"] if r["status"] == "early_rejected")
        if change == "early_cost":
            early["simulator_dose_calls"] = 13
        else:
            early["distance_lower_bound"] = archive["target"]["epsilon"]
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


@pytest.mark.parametrize("compound", ["rsl3", "erastin"])
@pytest.mark.parametrize("malformation", ["short", "long", "nan", "different_score"])
def test_accepted_curves_must_reconstruct_exact_fixed_target(archive, compound, malformation):
    row = next(r for r in archive["production"]["records"] if r["accepted"])
    if malformation == "short":
        row[compound].pop()
    elif malformation == "long":
        row[compound].append(row[compound][-1])
    elif malformation == "nan":
        row[compound][0] = float("nan")
    else:
        row[compound][0] += 0.01
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


def test_an_accepted_record_cannot_exceed_epsilon_even_by_one_float(archive):
    row = next(r for r in archive["production"]["records"] if r["accepted"])
    row["distance"] = row["distance_lower_bound"] = np.nextafter(archive["target"]["epsilon"], np.inf)
    with pytest.raises(ValueError, match="criterion"):
        driver.summarize_run(archive)


def test_every_evaluated_production_curve_is_retained_even_after_rejection(archive):
    for row in archive["production"]["records"]:
        assert row["erastin"] is not None
        if row["status"] == "early_rejected":
            assert row["rsl3"] is None
        else:
            assert row["rsl3"] is not None


@pytest.mark.parametrize("compound", ["rsl3", "erastin"])
@pytest.mark.parametrize("malformation", ["missing", "short", "nan", "changed_curve"])
def test_rejected_full_curves_are_required_and_must_reconstruct_score(archive, compound, malformation):
    row = next(r for r in archive["production"]["records"]
               if r["status"] == "complete" and not r["accepted"])
    if malformation == "missing":
        row[compound] = None
    elif malformation == "short":
        row[compound].pop()
    elif malformation == "nan":
        row[compound][0] = float("nan")
    else:
        row[compound][0] += 0.01
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


def test_rejected_full_score_cannot_change_while_staying_above_epsilon(archive):
    row = next(r for r in archive["production"]["records"]
               if r["status"] == "complete" and not r["accepted"])
    row["distance"] += 0.01
    row["distance_lower_bound"] = row["distance"]
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


@pytest.mark.parametrize("malformation", ["missing", "short", "nan", "changed_curve", "changed_bound", "uncomputed_rsl3"])
def test_early_rejection_requires_its_actual_erastin_curve_and_bound(archive, malformation):
    row = next(r for r in archive["production"]["records"] if r["status"] == "early_rejected")
    if malformation == "missing":
        row["erastin"] = None
    elif malformation == "short":
        row["erastin"].pop()
    elif malformation == "nan":
        row["erastin"][0] = float("nan")
    elif malformation == "changed_curve":
        row["erastin"][0] += 0.01
    elif malformation == "changed_bound":
        row["distance_lower_bound"] += 0.01  # Still strict, but no longer its curve's RMSE.
    else:
        row["rsl3"] = [0.0] * len(archive["target"]["rsl3_doses_um"])
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


def test_rejected_curves_are_validated_even_when_the_run_has_no_accepted_points(archive):
    for row in archive["production"]["records"]:
        if row["accepted"]:
            row.update(accepted=False, distance=1.0, distance_lower_bound=1.0,
                       rsl3=[0.0] * len(archive["target"]["rsl3_doses_um"]),
                       erastin=[1.0] * len(archive["target"]["erastin_doses_um"]))
    rejected = next(r for r in archive["production"]["records"] if r["status"] == "complete")
    rejected["erastin"] = None
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


@pytest.mark.parametrize("change", [
    "pilot_attempts", "pilot_calls", "island_count", "island_stream", "initial_cost", "initial_cloud", "initial_seed",
    "history_length", "move_attempts", "evaluated_count", "early_cost", "negative_count",
    "threshold_below_epsilon", "threshold_increases", "resampling_flag", "final_score",
    "final_cloud", "ancestor", "reached_flag", "kernel_mean", "kernel_scale", "kernel_weight",
])
def test_pilot_and_kernel_integrity_rejects_inconsistent_archive(archive, change):
    pilot = archive["pilot"]
    island = pilot["islands"][0]
    history = island["history"]
    if change == "pilot_attempts":
        pilot["attempts"] += 1
    elif change == "pilot_calls":
        pilot["simulator_dose_calls"] += 1
    elif change == "island_count":
        pilot["islands"].pop()
    elif change == "island_stream":
        island["spawn_key"] = [0, 1]
    elif change == "initial_cost":
        island["initial_simulator_dose_calls"] -= 1
    elif change == "initial_cloud":
        island["initial_points"][0][0] = 0.5
    elif change == "initial_seed":
        archive["seed"] += 1
    elif change == "history_length":
        history.pop()
    elif change == "move_attempts":
        history[0]["attempts"] += 1
    elif change == "evaluated_count":
        history[0]["evaluated"] += 1
    elif change == "early_cost":
        history[0]["simulator_dose_calls"] += 1
    elif change == "negative_count":
        history[0]["early_rejected"] = -1
    elif change == "threshold_below_epsilon":
        history[0]["threshold"] = 0
    elif change == "threshold_increases":
        history[1]["threshold"] = history[0]["threshold"] + 1
    elif change == "resampling_flag":
        history[0]["resampled"] = not history[0]["resampled"]
    elif change == "final_score":
        island["final_distances"][0] = history[-1]["threshold"] + 1
    elif change == "final_cloud":
        pilot["final_points"][0][0] += 0.01
    elif change == "ancestor":
        island["final_ancestors"][0] = archive["plan"]["pilot"]["particles_per_island"]
    elif change == "reached_flag":
        pilot["all_islands_reached_epsilon"] = not pilot["all_islands_reached_epsilon"]
    elif change == "kernel_mean":
        archive["proposal"]["means"][0][0] += 0.01
    elif change == "kernel_scale":
        archive["proposal"]["scales"][0][0] += 0.01
    else:
        archive["proposal"]["uniform_weight"] = 0.3
    with pytest.raises(ValueError):
        driver.summarize_run(archive)


def _adequate_runs():
    return [{"seed": seed, "adequacy": {"passed": True}, "all_islands_reached_epsilon": True,
             "prior_widths": dict.fromkeys(driver.NAMES, 1.0),
             "diagnostic_quantiles": {n: {"q2_5": 0.2, "median": 0.5, "q97_5": 0.8}
                                      for n in driver.NAMES},
             "heldout": {"rmse_quantiles": [0.1, 0.2, 0.3]},
             "importance": {"normalizer_estimate": 0.01, "normalizer_mcse": 0.001}}
            for seed in driver.SEEDS]


@pytest.mark.parametrize("case", ["missing", "duplicate", "unplanned", "pilot_failed", "run_failed"])
def test_every_planned_independent_run_and_pilot_must_pass(case):
    runs = _adequate_runs()
    assert driver.compare_runs(runs)["passed"]
    if case == "missing":
        runs.pop()
    elif case == "duplicate":
        runs[-1]["seed"] = runs[0]["seed"]
    elif case == "unplanned":
        runs[-1]["seed"] = 123
    elif case == "pilot_failed":
        runs[-1]["all_islands_reached_epsilon"] = False
    else:
        runs[-1]["adequacy"]["passed"] = False
    assert not driver.compare_runs(runs)["passed"]


@pytest.fixture
def report_files(tmp_path, monkeypatch, target, settings, models):
    monkeypatch.setattr(driver, "PLAN", settings)
    monkeypatch.setattr(driver, "OUT", tmp_path)
    prerequisite = {"artifact": "proposal-synthetic-validation-v2.json", "sha256": "fixture-study",
                    "source_hashes": {"validator": "fixture-source"}, "passed": True}
    monkeypatch.setattr(driver, "verify_synthetic_study", lambda: copy.deepcopy(prerequisite))
    for seed in driver.SEEDS:
        archive = driver.run_experiment(seed, target, settings)
        archive.update(source_hashes={"driver": "fixture-frozen-code"}, synthetic_validation=prerequisite,
                       runtime={"extension_binaries": {"core.so": "fixture-binary"},
                                "python": "fixture-python", "numpy": "fixture-numpy", "scipy": "fixture-scipy"})
        driver.write_archive(driver.archive_path(tmp_path, seed), archive)
    uniform_bytes = b'{"fixture": "historical uniform result"}\n'
    (tmp_path / "joint-posterior.json").write_bytes(uniform_bytes)
    previous = {"target": target,
                "baseline": {"sha256": hashlib.sha256(uniform_bytes).hexdigest(), "accepted": 4, "draws": 40000},
                "runs": [{"seed": 11, "importance": {"ess": 1}, "ess_per_1000_simulator_dose_calls": 0.01}]}
    (tmp_path / "joint-importance-sampling.json").write_text(json.dumps(previous))
    return tmp_path


def test_report_uses_replayed_archives_and_hashed_matching_baselines(report_files):
    report = driver.build_report(report_files)
    assert len(report["runs"]) == 3 and not report["stability"]["passed"]
    for seed in driver.SEEDS:
        path = driver.archive_path(report_files, seed)
        assert report["archive_sha256"][path.name] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert report["baselines"]["first_importance"]["sha256"] == hashlib.sha256(
        (report_files / "joint-importance-sampling.json").read_bytes()).hexdigest()
    text = driver.render(report)
    assert "pooled inference is withheld" in text and "No pooled posterior is published" in text


@pytest.mark.parametrize("field", ["seed", "plan", "gates", "target", "source_hashes", "synthetic_validation",
                                    "extension_binaries", "python", "numpy", "scipy"])
def test_report_rejects_incompatible_run_metadata(report_files, field):
    path = driver.archive_path(report_files, driver.SEEDS[1])
    archive = driver.read_archive(path)
    if field == "seed":
        archive[field] = driver.SEEDS[0]
    elif field == "plan":
        archive[field]["production_attempts"] += 1
    elif field == "gates":
        archive[field]["minimum_ess"] = 1
    elif field == "target":
        archive[field]["epsilon"] += 0.01
    elif field == "source_hashes":
        archive[field]["driver"] = "different-source"
    elif field == "synthetic_validation":
        archive[field]["sha256"] = "different-study"
    else:
        archive["runtime"][field] = "different-runtime"
    driver.write_archive(path, archive)
    with pytest.raises(ValueError):
        driver.build_report(report_files)


@pytest.mark.parametrize("field", ["target", "uniform_bytes"])
def test_report_rejects_changed_historical_target_or_uniform_baseline(report_files, field):
    if field == "target":
        path = report_files / "joint-importance-sampling.json"
        previous = json.loads(path.read_text())
        previous["target"]["epsilon"] += 0.01
        path.write_text(json.dumps(previous))
    else:
        (report_files / "joint-posterior.json").write_text("changed baseline\n")
    with pytest.raises(ValueError):
        driver.build_report(report_files)


@pytest.mark.parametrize("condition", ["missing", "failed", "different_plan"])
def test_prerequisite_blocks_cli_before_biological_runtime_or_target_loading(
        tmp_path, monkeypatch, condition):
    import proposal_synthetic_validation_v2 as synthetic

    monkeypatch.setattr(driver, "OUT", tmp_path)
    monkeypatch.setattr(sys, "argv", ["abc_joint_resample.py", "--seed", str(driver.SEEDS[0]),
                                      "--output-dir", str(tmp_path / "out")])
    if condition != "missing":
        (tmp_path / "proposal-synthetic-validation-v2.json").write_text('{"fixture": true}')
        pilot_plan = copy.deepcopy(driver.PLAN["pilot"])
        if condition == "different_plan":
            pilot_plan["levels"] += 1
        monkeypatch.setattr(synthetic, "assemble", lambda payload: {
            "passed": condition != "failed", "pilot_plan": pilot_plan}, raising=False)

    def forbidden(*args):
        pytest.fail("missing or invalid prerequisite reached biological work")

    monkeypatch.setattr(driver.base, "runtime_provenance", forbidden)
    monkeypatch.setattr(driver.base, "prepare_target", forbidden)
    monkeypatch.setattr(driver, "run_experiment", forbidden)
    with pytest.raises(FileNotFoundError if condition == "missing" else ValueError):
        driver.main()
    assert not list(tmp_path.glob("out/*.gz"))


def test_verified_prerequisite_records_hash_of_the_exact_validated_bytes(tmp_path, monkeypatch):
    import proposal_synthetic_validation_v2 as synthetic

    monkeypatch.setattr(driver, "OUT", tmp_path)
    raw = b'{"fixture": "study"}\n'
    (tmp_path / "proposal-synthetic-validation-v2.json").write_bytes(raw)
    seen = []

    def assemble(payload):
        seen.append(payload)
        return {"passed": True, "pilot_plan": driver.PLAN["pilot"],
                "source_hashes": {"validator": "reviewed"}}

    monkeypatch.setattr(synthetic, "assemble", assemble, raising=False)
    result = driver.verify_synthetic_study()
    assert seen == [{"fixture": "study"}]
    assert result["sha256"] == hashlib.sha256(raw).hexdigest()
    assert result["passed"] and result["source_hashes"] == {"validator": "reviewed"}


def test_source_change_during_run_prevents_archive_publication(
        tmp_path, monkeypatch, target, settings, models):
    monkeypatch.setattr(driver, "PLAN", settings)
    monkeypatch.setattr(sys, "argv", ["abc_joint_resample.py", "--seed", str(driver.SEEDS[0]),
                                      "--output-dir", str(tmp_path)])
    monkeypatch.setattr(driver, "verify_synthetic_study", lambda: {"passed": True})
    versions = iter([{"driver": "before"}, {"driver": "after"}])
    monkeypatch.setattr(driver, "source_provenance", lambda: next(versions))
    monkeypatch.setattr(driver.base, "runtime_provenance", lambda: {})
    monkeypatch.setattr(driver.base, "prepare_target", lambda: target)
    with pytest.raises(SystemExit, match="source changed"):
        driver.main()
    assert not driver.archive_path(tmp_path, driver.SEEDS[0]).exists()
