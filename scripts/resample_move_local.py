"""Revision 2: local reflected moves plus frozen independence proposals.

The historical resample_move module remains unchanged for failed-study replay.
Only separate, independent production attempts contribute inference weights.
"""

import math

import numpy as np

import resample_move as historical
from resample_move import _evaluation, independent_move

PILOT_PLAN = {
    **historical.PILOT_PLAN, "kernel_neighbors": 32,
    "neighbor_gap_ratio": 3.0, "singleton_scale": 0.04,
    "local_probability": 0.5, "local_scale_multiplier": 0.5,
    "local_scale_floor": 0.0025, "local_scale_ceiling": 0.25,
}


def validate_plan(plan):
    historical.validate_plan(plan)
    if not math.isfinite(plan["local_probability"]) or not 0 < plan["local_probability"] < 1:
        raise ValueError("local probability must be strictly between zero and one")
    for key in ("local_scale_multiplier", "local_scale_floor", "local_scale_ceiling",
                "neighbor_gap_ratio", "singleton_scale"):
        if not math.isfinite(plan[key]) or plan[key] <= 0:
            raise ValueError("local scales must be positive and finite")
    if plan["neighbor_gap_ratio"] <= 1:
        raise ValueError("neighbor gap ratio must exceed one")
    if plan["local_scale_floor"] > plan["local_scale_ceiling"]:
        raise ValueError("local scale floor cannot exceed its ceiling")


def fit_proposal(points, plan=None):
    plan = PILOT_PLAN if plan is None else plan
    validate_plan(plan)
    points = np.asarray(points, dtype=float)
    scales = local_scales(points, plan)
    return historical.BoundedGaussianMixture(points.shape[1],
        uniform_weight=plan["uniform_weight"], means=points, scales=scales)


def local_scales(points, plan=None):
    """Stop a unique-center neighborhood at a large positive distance jump.

    This geometry heuristic does not certify modes or discover lost regions.
    Particle multiplicities remain in the density, but cannot erase variance.
    """
    plan = PILOT_PLAN if plan is None else plan
    validate_plan(plan)
    points = np.asarray(points, dtype=float)
    if (points.ndim != 2 or min(points.shape) < 1 or not np.isfinite(points).all()
            or np.any((points < 0) | (points > 1))):
        raise ValueError("kernel centers must be finite points in the unit cube")
    unique, inverse = np.unique(points, axis=0, return_inverse=True)
    standardizer = np.maximum(np.std(unique, axis=0), plan["kernel_floor"])
    scaled = unique / standardizer
    distances = np.sqrt(np.sum((scaled[:, None, :] - scaled[None, :, :])**2, axis=2))
    neighbors = np.argsort(distances, axis=1, kind="stable")[:, :min(plan["kernel_neighbors"], len(unique))]
    widths = []
    for index, nearest in enumerate(neighbors):
        radii = distances[index, nearest]
        for j in range(1, len(nearest) - 1):
            if radii[j] > 0 and radii[j + 1] > plan["neighbor_gap_ratio"] * radii[j]:
                nearest = nearest[:j + 1]
                break
        width = (np.full(points.shape[1], plan["singleton_scale"]) if len(nearest) == 1
                 else plan["kernel_scale"] * np.std(unique[nearest], axis=0))
        widths.append(np.clip(width, plan["kernel_floor"], plan["kernel_ceiling"]))
    return np.asarray(widths)[inverse]


def local_step(proposal, plan=None):
    """One shared diagonal scale, frozen for a level, counting particle copies."""
    plan = PILOT_PLAN if plan is None else plan
    validate_plan(plan)
    return np.clip(plan["local_scale_multiplier"] * np.median(proposal.scales, axis=0),
                   plan["local_scale_floor"], plan["local_scale_ceiling"])


def reflect_cube(values):
    """Repeated independent coordinate reflections; never clip at a boundary."""
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("reflection requires finite values")
    folded = np.remainder(values, 2.0)
    return np.where(folded <= 1.0, folded, 2.0 - folded)


def mixed_move(points, distances, ancestors, proposal, threshold, rng, evaluate,
               *, step, local_probability):
    """Mixture of symmetric reflected random walk and independence MH.

    The diagonal step is shared by all states and fixed for the level. Each
    reflected Gaussian coordinate has a symmetric image-sum transition density;
    a local proposal therefore needs only the hard-target acceptance test.
    Per-state bandwidths or correlated reflection would require another ratio.
    """
    points = np.asarray(points, dtype=float).copy()
    distances = np.asarray(distances, dtype=float).copy()
    ancestors = np.asarray(ancestors).copy()
    step = np.asarray(step, dtype=float)
    if (points.ndim != 2 or len(points) < 1 or distances.shape != (len(points),) or
            ancestors.shape != (len(points),) or not np.isfinite(points).all() or
            np.any((points < 0) | (points > 1)) or not np.isfinite(distances).all() or
            not math.isfinite(threshold) or threshold <= 0 or
            np.any((distances < 0) | (distances > threshold))):
        raise ValueError("MH states must be valid and inside the conditioning region")
    if step.shape != (points.shape[1],) or not np.isfinite(step).all() or np.any(step <= 0):
        raise ValueError("local step must be a shared positive finite diagonal vector")
    if not math.isfinite(local_probability) or not 0 <= local_probability <= 1:
        raise ValueError("local probability must be between zero and one")
    # Fixed stream order includes global candidates and MH uniforms for every
    # particle, even when its local move does not use them.
    proposed, _ = proposal.sample(rng, len(points))
    local = rng.random(len(points)) < local_probability
    proposed[local] = reflect_cube(points[local] + rng.normal(size=(int(local.sum()), points.shape[1])) * step)
    with np.errstate(divide="ignore"):
        log_u = np.log(rng.random(len(points)))
    eligible = local.copy()
    if np.any(~local):
        old_q = proposal.log_density(points[~local])
        new_q = proposal.log_density(proposed[~local])
        if not np.isfinite(old_q).all() or not np.isfinite(new_q).all():
            raise ValueError("global MH states need finite proposal densities")
        eligible[~local] = log_u[~local] <= np.minimum(0.0, old_q - new_q)
    accepted = calls = early = local_accepted = 0
    for index in np.flatnonzero(eligible):
        distance, cost = _evaluation(evaluate(proposed[index], threshold), threshold)
        calls += cost
        if distance is None:
            early += 1
        elif distance <= threshold:
            points[index], distances[index] = proposed[index], distance
            accepted += 1
            local_accepted += int(local[index])
    return points, distances, ancestors, {
        "attempts": len(points), "evaluated": int(eligible.sum()),
        "mh_rejected_without_evaluation": int((~eligible).sum()),
        "early_rejected": early, "accepted_moves": accepted,
        "simulator_dose_calls": calls, "local_attempts": int(local.sum()),
        "global_attempts": int((~local).sum()), "local_accepted_moves": local_accepted,
    }



def validate_pilot_archive(pilot, fixture=None, seed=None):
    """Reconcile revision-specific counters; the caller checks common facts."""
    validate_plan(pilot["plan"])
    for island in pilot["islands"]:
        for step in island["history"]:
            for key in ("local_attempts", "global_attempts", "local_accepted_moves"):
                if type(step[key]) is not int or step[key] < 0:
                    raise ValueError("local move counters must be nonnegative integers")
            local, global_ = step["local_attempts"], step["global_attempts"]
            accepted_local = step["local_accepted_moves"]
            if (local + global_ != step["attempts"] or local > step["evaluated"] or
                    step["mh_rejected_without_evaluation"] > global_ or
                    accepted_local > min(local, step["accepted_moves"]) or
                    step["accepted_moves"] - accepted_local > step["evaluated"] - local):
                raise ValueError("local and global move counters do not reconcile")


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
            step = local_step(proposal, plan)
            counters = {k: 0 for k in ("attempts", "evaluated", "mh_rejected_without_evaluation",
                                      "early_rejected", "accepted_moves", "simulator_dose_calls",
                                      "local_attempts", "global_attempts", "local_accepted_moves") }
            for _ in range(plan["moves_per_level"]):
                points, distances, ancestors, stats = mixed_move(
                    points, distances, ancestors, proposal, threshold, island_rng, evaluate,
                    step=step, local_probability=plan["local_probability"])
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
