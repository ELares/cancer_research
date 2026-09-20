"""Audit fixed-plan analytic pilot archives without replaying their MH paths.

Initial randomness and both endpoint scores are independently reproduced. The
intermediate particles, resampling choices, and Metropolis proposals are absent
from these archives: their trajectories and the recorded intermediate summaries
are therefore not independently replayed. Counter consistency is a weaker check.
"""

import math

import numpy as np

import resample_move_local as strategy

DIMENSION = 7
# A selected order statistic can be a particle's exact archived score. Replaying
# its analytic norm on another platform may place it a few ULPs above that same
# threshold; this allowance applies only to recomputed endpoint scores.
ENDPOINT_THRESHOLD_ULPS = 8


def _object(value, keys, label):
    if type(value) is not dict or set(value) != set(keys):
        raise ValueError(f"{label} must contain exactly its declared fields")


def _integer(value, label):
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def _real(value, label):
    if type(value) not in (int, float):
        raise ValueError(f"{label} must be a finite real number")
    try:
        finite = math.isfinite(value)
    except OverflowError:
        finite = False
    if not finite:
        raise ValueError(f"{label} must be a finite real number")
    return value


def _array(value, shape, label):
    def check(part, remaining):
        if not remaining:
            _real(part, label)
        elif type(part) is not list or len(part) != remaining[0]:
            raise ValueError(f"{label} has the wrong shape")
        else:
            for child in part:
                check(child, remaining[1:])

    check(value, shape)
    return np.asarray(value, dtype=float)


def _same(actual, expected, label):
    if type(actual) is not type(expected) or actual != expected:
        raise ValueError(f"stored {label} differs from the fixed pilot")


def _close(actual, expected, label):
    if not np.allclose(actual, expected, rtol=1e-12, atol=1e-14):
        raise ValueError(f"stored {label} differs from its recomputation")


def _scores(points, fixture, evaluate):
    scores = []
    for point in points:
        record = evaluate(point, fixture)
        if type(record) is not dict:
            raise ValueError("analytic evaluation must return a record")
        distance = _real(record.get("distance"), "analytic distance")
        if distance < 0:
            raise ValueError("analytic distance must be nonnegative")
        _same(record.get("status"), "complete", "analytic evaluation status")
        _same(record.get("simulator_dose_calls"), 0, "analytic evaluation cost")
        bound = _real(record.get("distance_lower_bound"), "analytic distance lower bound")
        _same(bound, distance, "complete analytic distance lower bound")
        scores.append(distance)
    return np.asarray(scores)


def validate_pilot(pilot, fixture, seed, evaluate, epsilon=1.0):
    """Validate a seven-dimensional revision-2 archive and return target reach.

    This checks stored facts, not the intervening MH trajectory. Only derived
    endpoint scores/summaries allow numerical tolerance; initial RNG coordinates,
    combined clouds, counts, plan, and reach flags are checked exactly. Stored
    endpoint scores must satisfy the conditioning threshold exactly; recomputed
    scores allow at most eight ULPs of threshold rounding in addition to passing
    the independent endpoint-score comparison. The input is never mutated.
    Malformed JSON scalar types (including bool-as-int) fail.
    """
    _integer(seed, "pilot seed")
    if _real(epsilon, "epsilon") <= 0:
        raise ValueError("epsilon must be positive")
    _object(pilot, ("plan", "islands", "final_points", "final_distances",
                    "all_islands_reached_epsilon", "attempts", "simulator_dose_calls"), "pilot")
    plan = strategy.PILOT_PLAN
    _object(pilot["plan"], plan, "pilot plan")
    for key, expected in plan.items():
        _same(pilot["plan"][key], expected, f"pilot plan {key}")
    islands = pilot["islands"]
    if type(islands) is not list or len(islands) != plan["islands"]:
        raise ValueError("pilot has the wrong number of islands")
    streams = np.random.SeedSequence(seed).spawn(2)[0].spawn(plan["islands"])
    n, levels, moves = (plan[key] for key in ("particles_per_island", "levels", "moves_per_level"))
    reached, final_clouds, final_scores = [], [], []
    for index, (island, stream) in enumerate(zip(islands, streams)):
        _object(island, ("island", "spawn_key", "initial_points", "initial_distances",
                         "initial_simulator_dose_calls", "final_points", "final_distances",
                         "final_ancestors", "first_epsilon_level", "history"), "island")
        _same(island["island"], index, "island index")
        key = island["spawn_key"]
        if type(key) is not list or any(type(v) is not int for v in key):
            raise ValueError("island spawn key must be a list of integers")
        _same(key, list(stream.spawn_key), "island RNG spawn key")
        initial = _array(island["initial_points"], (n, DIMENSION), "initial points")
        expected_initial = np.random.default_rng(stream).uniform(size=(n, DIMENSION))
        if not np.array_equal(initial, expected_initial):
            raise ValueError("initial pilot points differ from exact RNG replay")
        final = _array(island["final_points"], (n, DIMENSION), "final points")
        if np.any((final < 0) | (final > 1)):
            raise ValueError("final pilot points must lie in the unit cube")
        scores, stored_scores = {}, {}
        for name, points in (("initial", initial), ("final", final)):
            stored = _array(island[f"{name}_distances"], (n,), f"{name} distances")
            if np.any(stored < 0):
                raise ValueError("pilot distances must be nonnegative")
            scores[name] = _scores(points, fixture, evaluate)
            stored_scores[name] = stored
            _close(stored, scores[name], f"{name} scores")
        ancestor_list = island["final_ancestors"]
        if (type(ancestor_list) is not list or len(ancestor_list) != n or
                any(type(a) is not int or not 0 <= a < n for a in ancestor_list)):
            raise ValueError("final pilot ancestry must contain valid integer indices")
        ancestors = np.asarray(ancestor_list)
        history = island["history"]
        if type(history) is not list or len(history) != levels:
            raise ValueError("pilot has the wrong number of levels")
        # The quantile selects a stored endpoint score, so this rule can be
        # checked exactly without turning harmless platform-dependent score
        # rounding into a different archived resampling decision.
        first_threshold = max(epsilon, float(np.quantile(stored_scores["initial"], plan["retained_fraction"], method="higher")))
        initial_keep = np.flatnonzero(stored_scores["initial"] <= first_threshold)
        if not np.isin(ancestors, initial_keep).all():
            raise ValueError("final ancestors include a discarded initial particle")
        previous_threshold = max(epsilon, float(np.max(stored_scores["initial"])))
        previous_ancestors, previous_median = n, None
        first_reached = None
        minimum_retained = math.ceil(plan["retained_fraction"] * (n - 1)) + 1
        for level, step in enumerate(history, 1):
            _object(step, ("level", "threshold", "retained_before_resampling", "resampled",
                           "unique_ancestors", "unique_particles", "minimum_distance", "median_distance",
                           "attempts", "evaluated", "mh_rejected_without_evaluation", "early_rejected",
                           "accepted_moves", "simulator_dose_calls", "local_attempts", "global_attempts",
                           "local_accepted_moves"), "pilot history step")
            _same(step["level"], level, "pilot level")
            threshold = _real(step["threshold"], "pilot threshold")
            if not epsilon <= threshold <= previous_threshold:
                raise ValueError("pilot thresholds must be nonincreasing and at least epsilon")
            if level == 1:
                if threshold != first_threshold:
                    raise ValueError("stored first threshold differs from the fixed quantile rule")
                _same(step["retained_before_resampling"], len(initial_keep), "initial retained count")
            elif threshold == previous_threshold:
                _same(step["retained_before_resampling"], n, "unrestricted retained count")
            if previous_median is not None and threshold < max(epsilon, previous_median):
                raise ValueError("pilot threshold is below the preceding median")
            for key in ("attempts", "evaluated", "mh_rejected_without_evaluation", "early_rejected",
                        "accepted_moves", "simulator_dose_calls", "retained_before_resampling",
                        "unique_ancestors", "unique_particles", "local_attempts", "global_attempts",
                        "local_accepted_moves"):
                _integer(step[key], key)
            _same(step["attempts"], n * moves, "move attempts")
            _same(step["evaluated"] + step["mh_rejected_without_evaluation"], step["attempts"], "move accounting")
            _same(step["early_rejected"], 0, "complete analytic scoring")
            _same(step["simulator_dose_calls"], 0, "analytic move cost")
            if step["accepted_moves"] > step["evaluated"]:
                raise ValueError("accepted moves exceed evaluated moves")
            retained = step["retained_before_resampling"]
            if not minimum_retained <= retained <= n:
                raise ValueError("retained count violates the fixed quantile rule")
            _same(step["resampled"], retained < n, "resampling flag")
            diversity = step["unique_ancestors"]
            if not 1 <= diversity <= min(previous_ancestors, retained):
                raise ValueError("pilot ancestry diversity cannot increase")
            if not step["resampled"] and diversity != previous_ancestors:
                raise ValueError("ancestry diversity cannot change without resampling")
            if not 1 <= step["unique_particles"] <= n:
                raise ValueError("invalid pilot particle diversity")
            minimum = _real(step["minimum_distance"], "minimum distance")
            median = _real(step["median_distance"], "median distance")
            if not 0 <= minimum <= median <= threshold:
                raise ValueError("pilot distance summaries violate their threshold")
            if threshold == epsilon and first_reached is None:
                first_reached = level
            previous_threshold, previous_ancestors, previous_median = threshold, diversity, median
        _same(island["first_epsilon_level"], first_reached, "first epsilon level")
        _same(island["initial_simulator_dose_calls"], 0, "analytic initial cost")
        stored_final_scores = np.asarray(island["final_distances"], dtype=float)
        replay_threshold = previous_threshold + ENDPOINT_THRESHOLD_ULPS * math.ulp(previous_threshold)
        if (np.any(stored_final_scores > previous_threshold) or
                np.any(scores["final"] > replay_threshold)):
            raise ValueError("final pilot escaped its conditioning region")
        _close(history[-1]["minimum_distance"], np.min(scores["final"]), "final minimum distance")
        _close(history[-1]["median_distance"], np.median(scores["final"]), "final median distance")
        _same(history[-1]["unique_particles"], len(np.unique(final, axis=0)), "final particle diversity")
        _same(history[-1]["unique_ancestors"], len(np.unique(ancestors)), "final ancestry diversity")
        final_clouds.append(final)
        final_scores.append(stored_final_scores)
        reached.append(first_reached is not None)
    for name, shape, parts in (("final_points", (len(islands) * n, DIMENSION), final_clouds),
                               ("final_distances", (len(islands) * n,), final_scores)):
        combined = _array(pilot[name], shape, f"combined {name}")
        if not np.array_equal(combined, np.concatenate(parts)):
            raise ValueError(f"combined {name} differs from its islands")
    _same(pilot["attempts"], len(islands) * (n + levels * n * moves), "total pilot attempts")
    _same(pilot["simulator_dose_calls"], 0, "analytic total cost")
    _same(pilot["all_islands_reached_epsilon"], all(reached), "all-islands target flag")
    strategy.validate_pilot_archive(pilot, fixture, seed)
    return all(reached)
