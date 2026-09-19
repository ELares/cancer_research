"""Bounded resample/move pilots used only to train a frozen proposal.

No particle or resampling multiplicity is an independent posterior observation.
Inference uses a separate stream of independent, importance-weighted attempts.
"""

import math

import numpy as np

from bounded_proposal import BoundedGaussianMixture


PILOT_PLAN = {
    "islands": 2, "particles_per_island": 256, "levels": 14,
    "moves_per_level": 4, "retained_fraction": 0.5,
    "uniform_weight": 0.2, "kernel_neighbors": 32,
    "kernel_scale": 1.5, "kernel_floor": 0.005, "kernel_ceiling": 0.5,
}


def validate_plan(plan):
    for name in ("islands", "particles_per_island", "levels", "moves_per_level", "kernel_neighbors"):
        value = plan[name]
        if type(value) is not int or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if plan["particles_per_island"] < 2:
        raise ValueError("each island needs at least two particles")
    for name in ("retained_fraction", "uniform_weight"):
        if not math.isfinite(plan[name]) or not 0 < plan[name] < 1:
            raise ValueError(f"{name} must be strictly between zero and one")
    for name in ("kernel_scale", "kernel_floor", "kernel_ceiling"):
        if not math.isfinite(plan[name]) or plan[name] <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if plan["kernel_floor"] > plan["kernel_ceiling"]:
        raise ValueError("kernel floor cannot exceed its ceiling")


def local_scales(points, plan=None):
    """Diagonal kernel widths from unique neighbors, retaining particle weights.

    Duplicate particles get the same width; their multiplicity remains in the
    mixture. Duplicates cannot create a falsely zero neighborhood variance.
    """
    plan = PILOT_PLAN if plan is None else plan
    validate_plan(plan)
    points = np.asarray(points, dtype=float)
    if (points.ndim != 2 or min(points.shape) < 1 or
            not np.all(np.isfinite(points)) or np.any((points < 0) | (points > 1))):
        raise ValueError("kernel centers must be finite points in the unit cube")
    unique, inverse = np.unique(points, axis=0, return_inverse=True)
    standardizer = np.maximum(np.std(unique, axis=0), plan["kernel_floor"])
    scaled = unique / standardizer
    squared = np.sum((scaled[:, None, :] - scaled[None, :, :])**2, axis=2)
    n_neighbors = min(plan["kernel_neighbors"], len(unique))
    neighbors = np.argsort(squared, axis=1, kind="stable")[:, :n_neighbors]
    widths = plan["kernel_scale"] * np.std(unique[neighbors], axis=1)
    widths = np.clip(widths, plan["kernel_floor"], plan["kernel_ceiling"])
    return widths[inverse]


def fit_proposal(points, plan=None):
    plan = PILOT_PLAN if plan is None else plan
    points = np.asarray(points, dtype=float)
    scales = local_scales(points, plan)
    return BoundedGaussianMixture(points.shape[1], uniform_weight=plan["uniform_weight"],
                                  means=points, scales=scales)


def _evaluation(value, threshold):
    """Validate a full score or a proved early-rejection bound."""
    distance = value["distance"]
    calls = value["simulator_dose_calls"]
    if type(calls) is not int or calls < 0:
        raise ValueError("evaluation cost must be a nonnegative integer")
    if distance is None:
        bound = value["distance_lower_bound"]
        if (threshold is None or bound is None or not math.isfinite(bound) or
                bound <= threshold):
            raise ValueError("a missing score requires a strict early-rejection bound")
    elif not math.isfinite(distance) or distance < 0:
        raise ValueError("distance must be finite and nonnegative")
    return distance, calls


def independent_move(points, distances, ancestors, proposal, threshold, rng, evaluate):
    """One fixed-proposal independence-MH sweep under a uniform hard target.

    Accept with I[d(y)<=threshold] * min(1, q(x)/q(y)). The proposal is fixed
    throughout the sweep. A failed MH test can skip the expensive score exactly.
    """
    points = np.asarray(points, dtype=float).copy()
    distances = np.asarray(distances, dtype=float).copy()
    ancestors = np.asarray(ancestors)
    if (points.ndim != 2 or len(points) < 1 or distances.shape != (len(points),) or
            ancestors.shape != (len(points),) or not np.all(np.isfinite(distances)) or
            not math.isfinite(threshold) or threshold <= 0 or
            np.any((distances < 0) | (distances > threshold))):
        raise ValueError("MH states must have valid scores inside the conditioning region")
    proposed, _ = proposal.sample(rng, len(points))
    log_q_old, log_q_new = proposal.log_density(points), proposal.log_density(proposed)
    if not np.all(np.isfinite(log_q_old)) or not np.all(np.isfinite(log_q_new)):
        raise ValueError("MH states need finite proposal densities")
    # U=0 is a valid always-pass event in the floating-point uniform generator.
    with np.errstate(divide="ignore"):
        log_u = np.log(rng.random(len(points)))
    eligible = log_u <= np.minimum(0.0, log_q_old - log_q_new)
    accepted = 0
    calls = 0
    early = 0
    for i in np.flatnonzero(eligible):
        distance, cost = _evaluation(evaluate(proposed[i], threshold), threshold)
        calls += cost
        if distance is None:
            early += 1
        elif distance <= threshold:
            points[i], distances[i] = proposed[i], distance
            accepted += 1
    return points, distances, np.asarray(ancestors).copy(), {
        "attempts": len(points), "evaluated": int(eligible.sum()),
        "mh_rejected_without_evaluation": int((~eligible).sum()),
        "early_rejected": early, "accepted_moves": accepted,
        "simulator_dose_calls": calls,
    }


def train_pilot(evaluate, dimension, epsilon, rng, plan=None, progress=None):
    """Train independent islands without using production randomness.

    ``rng`` is a SeedSequence, whose children define the island streams.
    ``evaluate(point, threshold=None)`` returns a distance, a lower bound for a
    proved early rejection (if distance is None), and simulator_dose_calls.
    The initial population requires full scores; later evaluations may use
    exact rejection at that level's intermediate threshold.
    """
    plan = dict(PILOT_PLAN if plan is None else plan)
    validate_plan(plan)
    if type(dimension) is not int or dimension < 1:
        raise ValueError("dimension must be a positive integer")
    if not math.isfinite(epsilon) or epsilon <= 0:
        raise ValueError("epsilon must be positive and finite")
    if not isinstance(rng, np.random.SeedSequence):
        raise ValueError("pilot randomness must be a SeedSequence")
    streams = rng.spawn(plan["islands"])
    islands = []
    for island_id, stream in enumerate(streams):
        island_rng = np.random.default_rng(stream)
        n = plan["particles_per_island"]
        points = island_rng.uniform(size=(n, dimension))
        initial_points = points.tolist()
        distances = []
        initial_calls = 0
        for point in points:
            distance, calls = _evaluation(evaluate(point, None), None)
            distances.append(distance)
            initial_calls += calls
        distances = np.asarray(distances)
        initial_distances = distances.tolist()
        ancestors = np.arange(n)
        history = []
        previous = max(epsilon, float(np.max(distances)))
        first_reached = None
        for level in range(plan["levels"]):
            candidate = float(np.quantile(distances, plan["retained_fraction"], method="higher"))
            threshold = max(epsilon, min(previous, candidate))
            keep = np.flatnonzero(distances <= threshold)
            if len(keep) == 0:
                raise ValueError("pilot threshold removed every particle")
            # Resample only when restriction actually discards particles. At
            # epsilon, repeated needless resampling would erase ancestry.
            parents = (island_rng.choice(keep, size=n, replace=True)
                       if len(keep) < n else np.arange(n))
            points, distances, ancestors = points[parents], distances[parents], ancestors[parents]
            proposal = fit_proposal(points, plan)
            counters = {k: 0 for k in ("attempts", "evaluated", "mh_rejected_without_evaluation",
                                      "early_rejected", "accepted_moves", "simulator_dose_calls")}
            for _ in range(plan["moves_per_level"]):
                points, distances, ancestors, stats = independent_move(
                    points, distances, ancestors, proposal, threshold, island_rng, evaluate)
                for key in counters:
                    counters[key] += stats[key]
            if np.any(distances > threshold):
                raise ValueError("pilot move escaped its conditioning region")
            if threshold == epsilon and first_reached is None:
                first_reached = level + 1
            history.append({"level": level + 1, "threshold": threshold,
                            "retained_before_resampling": len(keep),
                            "resampled": len(keep) < n, "unique_ancestors": len(np.unique(ancestors)),
                            "unique_particles": len(np.unique(points, axis=0)),
                            "minimum_distance": float(np.min(distances)),
                            "median_distance": float(np.median(distances)), **counters})
            previous = threshold
            if progress:
                progress(island_id, history[-1])
        islands.append({"island": island_id, "spawn_key": list(stream.spawn_key),
                        "initial_points": initial_points, "initial_distances": initial_distances,
                        "initial_simulator_dose_calls": initial_calls,
                        "final_points": points.tolist(), "final_distances": distances.tolist(),
                        "final_ancestors": ancestors.tolist(), "first_epsilon_level": first_reached,
                        "history": history})
    final_points = np.concatenate([island["final_points"] for island in islands])
    final_scores = np.concatenate([island["final_distances"] for island in islands])
    return {"plan": plan, "islands": islands, "final_points": final_points.tolist(),
            "final_distances": final_scores.tolist(),
            "all_islands_reached_epsilon": all(i["first_epsilon_level"] is not None for i in islands),
            "attempts": sum(len(i["initial_points"]) + sum(h["attempts"] for h in i["history"]) for i in islands),
            "simulator_dose_calls": sum(i["initial_simulator_dose_calls"] +
                                        sum(h["simulator_dose_calls"] for h in i["history"]) for i in islands)}
