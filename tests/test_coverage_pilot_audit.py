"""Artificial archive checks; no prospective fixture or sampler is executed."""

import copy
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import coverage_pilot_audit as audit  # noqa: E402

TEST_SEED = 41017


def evaluate(point, fixture, threshold=None):
    if fixture == "test-unreached":
        distance = 2.0
    elif fixture == "test-resampled":
        distance = 4.0 * float(point[0])
    else:
        distance = 0.25 + 0.25 * float(point[0])
    return {"status": "complete", "distance": distance,
            "distance_lower_bound": distance, "simulator_dose_calls": 0}


def artificial_pilot(fixture="test-reached"):
    """Supply trace doubles, not draws from a competing sampler design.

    All move attempts are declared global and rejected before evaluation. A
    resampled trace declares every selected parent to be its minimum-score
    initial particle. These trace facts are possible but not a replay of the
    island streams past initialization; the validator does not claim otherwise.
    """
    plan = copy.deepcopy(audit.strategy.PILOT_PLAN)
    n, levels, moves = (plan[key] for key in ("particles_per_island", "levels", "moves_per_level"))
    streams = np.random.SeedSequence(TEST_SEED).spawn(2)[0].spawn(plan["islands"])
    islands = []
    for index, stream in enumerate(streams):
        initial = np.random.default_rng(stream).uniform(size=(n, audit.DIMENSION))
        scores = [evaluate(point, fixture)["distance"] for point in initial]
        first_threshold = max(1.0, float(np.quantile(scores, plan["retained_fraction"], method="higher")))
        retained = sum(score <= first_threshold for score in scores)
        if fixture == "test-resampled":
            ancestor = int(np.argmin(scores))
            ancestors = [ancestor] * n
            final = np.repeat(initial[ancestor][None, :], n, axis=0)
        else:
            ancestors = list(range(n))
            final = initial.copy()
        final_scores = [evaluate(point, fixture)["distance"] for point in final]
        later_threshold = max(1.0, float(np.quantile(final_scores, plan["retained_fraction"], method="higher")))
        history = []
        for level in range(1, levels + 1):
            threshold = first_threshold if level == 1 else later_threshold
            keep = retained if level == 1 else n
            history.append({"level": level, "threshold": threshold,
                            "retained_before_resampling": keep, "resampled": keep < n,
                            "unique_ancestors": len(set(ancestors)),
                            "unique_particles": len(np.unique(final, axis=0)),
                            "minimum_distance": min(final_scores),
                            "median_distance": float(np.median(final_scores)),
                            "attempts": n * moves, "evaluated": 0,
                            "mh_rejected_without_evaluation": n * moves,
                            "early_rejected": 0, "accepted_moves": 0, "simulator_dose_calls": 0,
                            "local_attempts": 0, "global_attempts": n * moves,
                            "local_accepted_moves": 0})
        first_reached = next((step["level"] for step in history if step["threshold"] == 1.0), None)
        islands.append({"island": index, "spawn_key": list(stream.spawn_key),
                        "initial_points": initial.tolist(), "initial_distances": scores,
                        "initial_simulator_dose_calls": 0, "final_points": final.tolist(),
                        "final_distances": final_scores, "final_ancestors": ancestors,
                        "first_epsilon_level": first_reached, "history": history})
    return {"plan": plan, "islands": islands,
            "final_points": [point for island in islands for point in island["final_points"]],
            "final_distances": [score for island in islands for score in island["final_distances"]],
            "all_islands_reached_epsilon": all(island["first_epsilon_level"] is not None for island in islands),
            "attempts": plan["islands"] * (n + levels * n * moves), "simulator_dose_calls": 0}


@pytest.fixture
def pilot():
    return artificial_pilot()


@pytest.mark.parametrize("fixture, reached", [("test-reached", True), ("test-unreached", False),
                                              ("test-resampled", True)])
def test_valid_trace_recomputes_reach_without_mutation_or_training(fixture, reached, monkeypatch):
    pilot = artificial_pilot(fixture)
    original = copy.deepcopy(pilot)

    def forbidden(*args, **kwargs):
        pytest.fail("archive auditing must not train a pilot")

    monkeypatch.setattr(audit.strategy, "train_pilot", forbidden)
    assert audit.validate_pilot(pilot, fixture, TEST_SEED, evaluate) is reached
    assert pilot == original


def test_exact_initial_randomness_rejects_a_single_ulp_change(pilot):
    old = pilot["islands"][0]["initial_points"][0][0]
    pilot["islands"][0]["initial_points"][0][0] = float(np.nextafter(old, 1.0))
    with pytest.raises(ValueError, match="exact RNG"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


def test_each_endpoint_is_independently_scored(pilot):
    calls = []

    def recording(point, fixture, threshold=None):
        calls.append((fixture, threshold))
        return evaluate(point, fixture, threshold)

    audit.validate_pilot(pilot, "test-reached", TEST_SEED, recording)
    assert calls == [("test-reached", None)] * (2 * 2 * 256)


def test_strategy_counter_hook_is_required(pilot, monkeypatch):
    def rejecting_hook(stored, fixture, seed):
        assert stored is pilot and fixture == "test-reached" and seed == TEST_SEED
        raise ValueError("strategy-specific rejection")

    monkeypatch.setattr(audit.strategy, "validate_pilot_archive", rejecting_hook)
    with pytest.raises(ValueError, match="strategy-specific"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("seed", [True, -1, 1.0, "41017", None])
def test_invalid_seed_types_fail(pilot, seed):
    with pytest.raises(ValueError, match="seed"):
        audit.validate_pilot(pilot, "test-reached", seed, evaluate)


@pytest.mark.parametrize("epsilon", [True, False, 0, -1, "1", None, float("nan"), float("inf"), 10**1000])
def test_invalid_epsilon_fails(pilot, epsilon):
    with pytest.raises(ValueError, match="epsilon"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate, epsilon)


@pytest.mark.parametrize("key, value", [("islands", 2.0), ("particles_per_island", 255),
    ("levels", 13), ("moves_per_level", True), ("uniform_weight", "0.2"),
    ("local_probability", 0.7), ("neighbor_gap_ratio", 3), ("kernel_scale", float("nan"))])
def test_fixed_plan_checks_values_and_scalar_types(pilot, key, value):
    pilot["plan"][key] = value
    with pytest.raises(ValueError, match="pilot plan"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("key, value", [("status", "early_rejected"), ("status", True),
    ("distance", "0.5"), ("distance", True), ("distance", None), ("distance", -1),
    ("distance", float("inf")), ("distance", float("nan")),
    ("simulator_dose_calls", False), ("simulator_dose_calls", 1),
    ("distance_lower_bound", None), ("distance_lower_bound", "0.3"),
    ("distance_lower_bound", 999)])
def test_analytic_evaluator_must_provide_complete_finite_zero_cost_records(pilot, key, value):
    def malformed(point, fixture, threshold=None):
        return {**evaluate(point, fixture), key: value}

    with pytest.raises(ValueError, match="analytic"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, malformed)


@pytest.mark.parametrize("record", [None, [], 1, {"status": "complete"}])
def test_missing_analytic_record_fields_fail(pilot, record):
    with pytest.raises(ValueError, match="analytic"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, lambda *args: record)


@pytest.mark.parametrize("field", ["initial_points", "final_points", "initial_distances", "final_distances"])
@pytest.mark.parametrize("value", [True, "0.2", None, float("nan"), float("inf")])
def test_array_scalar_types_and_finiteness(pilot, field, value):
    points = pilot["islands"][0][field]
    if field.endswith("points"):
        points[0][0] = value
    else:
        points[0] = value
    with pytest.raises(ValueError, match="finite real"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("field", ["initial_points", "final_points", "initial_distances", "final_distances"])
def test_endpoint_array_shape_is_fixed(pilot, field):
    pilot["islands"][0][field].pop()
    with pytest.raises(ValueError, match="wrong shape"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("field", ["initial_distances", "final_distances"])
def test_altered_endpoint_scores_fail_even_if_summaries_are_left_unchanged(pilot, field):
    pilot["islands"][0][field][0] += 0.01
    with pytest.raises(ValueError, match="scores"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_final_points_cannot_escape_unit_cube(pilot, value):
    pilot["islands"][0]["final_points"][0][0] = value
    with pytest.raises(ValueError, match="unit cube"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("value", [True, 0.0, "0", -1, 256])
def test_ancestry_indices_have_exact_integer_types_and_range(pilot, value):
    pilot["islands"][0]["final_ancestors"][0] = value
    with pytest.raises(ValueError, match="ancestry"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


def test_final_ancestry_cannot_contain_an_initially_discarded_parent():
    pilot = artificial_pilot("test-resampled")
    island = pilot["islands"][0]
    island["final_ancestors"][0] = int(np.argmax(island["initial_distances"]))
    with pytest.raises(ValueError, match="discarded initial"):
        audit.validate_pilot(pilot, "test-resampled", TEST_SEED, evaluate)


@pytest.mark.parametrize("key, value", [("level", True), ("level", 2), ("threshold", "1"),
    ("threshold", 0.99), ("threshold", 1.01), ("minimum_distance", -1),
    ("minimum_distance", True), ("median_distance", float("nan")),
    ("median_distance", 1.1), ("resampled", 0), ("resampled", True),
    ("attempts", 1023), ("attempts", 1024.0), ("evaluated", 1),
    ("evaluated", False), ("accepted_moves", 1), ("early_rejected", 1),
    ("simulator_dose_calls", False), ("simulator_dose_calls", 1),
    ("retained_before_resampling", 255), ("unique_particles", 0),
    ("unique_particles", 257), ("unique_ancestors", 255), ("local_attempts", -1),
    ("global_attempts", True), ("local_accepted_moves", 1)])
def test_history_facts_are_reconciled(pilot, key, value):
    pilot["islands"][0]["history"][0][key] = value
    with pytest.raises(ValueError):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


def test_thresholds_cannot_increase():
    pilot = artificial_pilot("test-unreached")
    pilot["islands"][0]["history"][1]["threshold"] = 2.1
    with pytest.raises(ValueError, match="nonincreasing"):
        audit.validate_pilot(pilot, "test-unreached", TEST_SEED, evaluate)


def test_first_threshold_is_recomputed_from_fixed_quantile():
    pilot = artificial_pilot("test-resampled")
    pilot["islands"][0]["history"][0]["threshold"] -= 0.01
    with pytest.raises(ValueError, match="first threshold"):
        audit.validate_pilot(pilot, "test-resampled", TEST_SEED, evaluate)


def test_first_quantile_decision_is_exact_even_for_one_ulp_changes():
    pilot = artificial_pilot("test-resampled")
    step = pilot["islands"][0]["history"][0]
    step["threshold"] = float(np.nextafter(step["threshold"], float("inf")))
    with pytest.raises(ValueError, match="first threshold"):
        audit.validate_pilot(pilot, "test-resampled", TEST_SEED, evaluate)


def test_analytic_endpoint_rounding_does_not_change_stored_quantile_decisions():
    pilot = artificial_pilot("test-resampled")

    def rounding_variant(point, fixture, threshold=None):
        record = evaluate(point, fixture, threshold)
        distance = float(np.nextafter(record["distance"], float("inf")))
        return {**record, "distance": distance, "distance_lower_bound": distance}

    assert audit.validate_pilot(pilot, "test-resampled", TEST_SEED, rounding_variant)


@pytest.mark.parametrize("ulps", [1, 8])
def test_recomputed_endpoint_on_selected_threshold_allows_narrow_roundoff(ulps):
    pilot = artificial_pilot("test-unreached")

    def rounding_variant(point, fixture, threshold=None):
        record = evaluate(point, fixture, threshold)
        distance = record["distance"] + ulps * float(np.spacing(record["distance"]))
        return {**record, "distance": distance, "distance_lower_bound": distance}

    assert audit.validate_pilot(pilot, "test-unreached", TEST_SEED, rounding_variant) is False


def test_recomputed_endpoint_threshold_roundoff_is_narrower_than_score_tolerance():
    pilot = artificial_pilot("test-unreached")

    def rounding_variant(point, fixture, threshold=None):
        record = evaluate(point, fixture, threshold)
        distance = record["distance"] + 9 * float(np.spacing(record["distance"]))
        return {**record, "distance": distance, "distance_lower_bound": distance}

    with pytest.raises(ValueError, match="escaped its conditioning"):
        audit.validate_pilot(pilot, "test-unreached", TEST_SEED, rounding_variant)


def test_stored_endpoint_cannot_exceed_threshold_even_by_one_ulp():
    pilot = artificial_pilot("test-unreached")
    pilot["islands"][0]["final_distances"][0] = float(np.nextafter(2.0, float("inf")))
    with pytest.raises(ValueError, match="escaped its conditioning"):
        audit.validate_pilot(pilot, "test-unreached", TEST_SEED, evaluate)


def test_ancestry_diversity_cannot_reappear_after_resampling():
    pilot = artificial_pilot("test-resampled")
    pilot["islands"][0]["history"][1]["unique_ancestors"] = 2
    with pytest.raises(ValueError, match="cannot increase"):
        audit.validate_pilot(pilot, "test-resampled", TEST_SEED, evaluate)


@pytest.mark.parametrize("key", ["minimum_distance", "median_distance", "unique_particles", "unique_ancestors"])
def test_last_level_summaries_are_recomputed_from_endpoints(pilot, key):
    history = pilot["islands"][0]["history"]
    if key == "unique_ancestors":
        pilot["islands"][0]["final_ancestors"][-1] = 0
    else:
        history[-1][key] -= 0.01 if key.endswith("distance") else 1
    with pytest.raises(ValueError, match="final"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("field", ["final_points", "final_distances"])
def test_combined_cloud_is_exact_concatenation(pilot, field):
    # The builder shares nested lists; isolate the aggregate before corrupting it.
    pilot[field] = copy.deepcopy(pilot[field])
    if field == "final_points":
        old = pilot[field][0][0]
        pilot[field][0][0] = float(np.nextafter(old, 1.0))
    else:
        pilot[field][0] = float(np.nextafter(pilot[field][0], 1.0))
    with pytest.raises(ValueError, match="combined"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("key, value", [("attempts", 29183), ("attempts", 29184.0),
    ("simulator_dose_calls", False), ("simulator_dose_calls", 1),
    ("all_islands_reached_epsilon", False), ("all_islands_reached_epsilon", 1)])
def test_aggregate_budget_cost_and_reach_flag_are_recomputed(pilot, key, value):
    pilot[key] = value
    with pytest.raises(ValueError):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("value", [None, 2, True, 1.0])
def test_first_target_level_is_recomputed(pilot, value):
    pilot["islands"][0]["first_epsilon_level"] = value
    with pytest.raises(ValueError, match="first epsilon"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("key, value", [("island", False), ("island", 1),
    ("spawn_key", [0, False]), ("spawn_key", [1, 0]), ("spawn_key", "0,0"),
    ("initial_simulator_dose_calls", False), ("initial_simulator_dose_calls", 1)])
def test_island_identity_stream_and_initial_cost_are_strict(pilot, key, value):
    pilot["islands"][0][key] = value
    with pytest.raises(ValueError):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("location", ["pilot", "plan", "island", "step"])
def test_missing_and_undeclared_fields_fail(pilot, location):
    target = {"pilot": pilot, "plan": pilot["plan"], "island": pilot["islands"][0],
              "step": pilot["islands"][0]["history"][0]}[location]
    target["undeclared"] = "not part of the frozen trace"
    with pytest.raises(ValueError, match="declared fields"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)
    del target["undeclared"]
    del target[next(iter(target))]
    with pytest.raises(ValueError, match="declared fields"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


@pytest.mark.parametrize("location", ["islands", "history"])
def test_fixed_island_and_level_counts_cannot_be_shortened(pilot, location):
    target = pilot["islands"] if location == "islands" else pilot["islands"][0]["history"]
    target.pop()
    with pytest.raises(ValueError, match="wrong number"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


def test_local_and_global_counter_constraints_are_applied(pilot):
    step = pilot["islands"][0]["history"][0]
    step["local_attempts"], step["global_attempts"] = 1, 1023
    # Common attempts still reconcile, but a local attempt must be evaluated.
    with pytest.raises(ValueError, match="local and global"):
        audit.validate_pilot(pilot, "test-reached", TEST_SEED, evaluate)


def test_audit_scope_explicitly_excludes_independent_intermediate_replay():
    assert "not independently replayed" in audit.__doc__
    assert "not the intervening MH trajectory" in audit.validate_pilot.__doc__
