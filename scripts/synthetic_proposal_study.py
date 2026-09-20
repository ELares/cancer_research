"""Reusable replayable study engine for explicit pilot-design revisions.

The frozen first study supplies the unchanged analytic fixtures and assessments.
Pilot implementations are passed explicitly; archived data never choose code.
"""

import copy
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import scipy

import proposal_synthetic_validation as frozen

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class StudySpec:
    revision_id: str
    strategy_module: str
    positive_seeds: tuple
    source_paths: tuple
    prior_study_path: str
    prior_study_sha256: str


def _equal(actual, expected, label):
    if actual != expected:
        raise ValueError(f"stored {label} differs from this fixed study revision")


def _count(value, label):
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def _close(actual, expected, label):
    frozen._require_close(actual, expected, label)


def source_hashes(spec):
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for path in spec.source_paths}


def _specification(spec, strategy):
    if not spec.revision_id or len(spec.positive_seeds) != 3 or len(set(spec.positive_seeds)) != 3:
        raise ValueError("a revision requires an ID and three distinct positive seeds")
    if any(type(seed) is not int or seed < 0 for seed in spec.positive_seeds):
        raise ValueError("positive seeds must be nonnegative integers")
    if set(spec.positive_seeds) & set(frozen.POSITIVE_SEEDS + frozen.NEGATIVE_SEEDS):
        raise ValueError("revised positive seeds must differ from the prior study's seeds")
    _equal(getattr(strategy, "__name__", None), spec.strategy_module, "strategy module")
    if len(set(spec.source_paths)) != len(spec.source_paths):
        raise ValueError("numerical source paths must be distinct")
    prior_hash = hashlib.sha256((ROOT / spec.prior_study_path).read_bytes()).hexdigest()
    _equal(prior_hash, spec.prior_study_sha256, "prior study SHA256")
    return {"schema_version": 1, "revision_id": spec.revision_id,
            "strategy_module": spec.strategy_module, "pilot_plan": copy.deepcopy(strategy.PILOT_PLAN),
            "dimension": frozen.DIMENSION, "epsilon": frozen.EPSILON,
            "fixtures": list(frozen.FIXTURES), "production_attempts": frozen.PRODUCTION_ATTEMPTS,
            "positive_seeds": list(spec.positive_seeds), "negative_seeds": list(frozen.NEGATIVE_SEEDS),
            "gates": dict(frozen.GATES), "source_hashes": source_hashes(spec),
            "prior_study": {"path": spec.prior_study_path, "sha256": prior_hash}}


def check_pilot_archive(pilot, fixture, seed, strategy):
    """Validate common trace facts, with an optional strategy-specific hook.

    This replays initial randomness and analytic scores, not the intervening MH
    trajectory. Extra local-move counters belong to the strategy's optional
    validate_pilot_archive(pilot, fixture, seed) function.
    """
    plan = strategy.PILOT_PLAN
    _equal(pilot["plan"], plan, "pilot plan")
    islands = pilot["islands"]
    _equal([i["island"] for i in islands], list(range(plan["islands"])), "pilot islands")
    streams = np.random.SeedSequence(seed).spawn(2)[0].spawn(plan["islands"])
    n, levels, moves = (plan[k] for k in ("particles_per_island", "levels", "moves_per_level"))
    reached = []
    for island, stream in zip(islands, streams):
        _equal(island["spawn_key"], list(stream.spawn_key), "island RNG")
        initial = np.random.default_rng(stream).uniform(size=(n, frozen.DIMENSION))
        _close(island["initial_points"], initial, "initial pilot points")
        final = frozen._points(island["final_points"])
        if final.shape != (n, frozen.DIMENSION):
            raise ValueError("final pilot cloud has the wrong shape")
        for name, points in (("initial", initial), ("final", final)):
            scores = [frozen.evaluate(point, fixture)["distance"] for point in points]
            _close(island[f"{name}_distances"], scores, f"{name} pilot scores")
        ancestors = np.asarray(island["final_ancestors"])
        if ancestors.shape != (n,) or ancestors.dtype.kind not in "iu" or np.any((ancestors < 0) | (ancestors >= n)):
            raise ValueError("final pilot ancestry is invalid")
        history = island["history"]
        _equal([h["level"] for h in history], list(range(1, levels + 1)), "pilot levels")
        thresholds = np.asarray([h["threshold"] for h in history], dtype=float)
        if not np.isfinite(thresholds).all() or np.any(thresholds < frozen.EPSILON) or np.any(np.diff(thresholds) > 0):
            raise ValueError("pilot thresholds violate the fixed final threshold")
        final_scores = np.asarray(island["final_distances"])
        if np.any(final_scores > thresholds[-1]):
            raise ValueError("final pilot escaped its conditioning region")
        hits = np.flatnonzero(thresholds == frozen.EPSILON)
        first = int(hits[0] + 1) if len(hits) else None
        _equal(island["first_epsilon_level"], first, "first target level")
        reached.append(first is not None)
        _equal(island["initial_simulator_dose_calls"], 0, "analytic initial cost")
        for step in history:
            for key in ("attempts", "evaluated", "mh_rejected_without_evaluation", "early_rejected",
                        "accepted_moves", "simulator_dose_calls", "retained_before_resampling",
                        "unique_ancestors", "unique_particles"):
                _count(step[key], key)
            _equal(step["attempts"], n * moves, "move attempts")
            _equal(step["evaluated"] + step["mh_rejected_without_evaluation"], step["attempts"], "move accounting")
            _equal(step["early_rejected"], 0, "complete analytic scoring")
            _equal(step["simulator_dose_calls"], 0, "analytic move cost")
            if step["accepted_moves"] > step["evaluated"]:
                raise ValueError("accepted moves exceed evaluated moves")
            if not 1 <= step["retained_before_resampling"] <= n:
                raise ValueError("invalid retained particle count")
            if type(step["resampled"]) is not bool:
                raise ValueError("resampling flag must be a boolean")
            _equal(step["resampled"], step["retained_before_resampling"] < n, "resampling flag")
            if not 1 <= step["unique_particles"] <= n or not 1 <= step["unique_ancestors"] <= n:
                raise ValueError("invalid pilot diversity count")
            minimum, median = step["minimum_distance"], step["median_distance"]
            if (not np.isfinite([minimum, median]).all() or
                    not 0 <= minimum <= median <= step["threshold"]):
                raise ValueError("pilot level distance summaries violate their threshold")
        _close(history[-1]["minimum_distance"], np.min(final_scores), "final minimum distance")
        _close(history[-1]["median_distance"], np.median(final_scores), "final median distance")
        _equal(history[-1]["unique_particles"], len(np.unique(final, axis=0)), "final particle diversity")
        _equal(history[-1]["unique_ancestors"], len(np.unique(ancestors)), "final ancestry diversity")
    for name in ("final_points", "final_distances"):
        _close(pilot[name], np.concatenate([island[name] for island in islands]), f"combined {name}")
    _equal(pilot["attempts"], len(islands) * (n + levels * n * moves), "total pilot attempts")
    _equal(pilot["simulator_dose_calls"], 0, "analytic total cost")
    hook = getattr(strategy, "validate_pilot_archive", None)
    if hook is not None:
        hook(pilot, fixture, seed)
    return all(reached)


def _production(pilot, proposal, fixture, seed):
    stream = np.random.SeedSequence(seed).spawn(2)[1]
    points, components = proposal.sample(np.random.default_rng(stream), frozen.PRODUCTION_ATTEMPTS)
    assessment = frozen.summarize(points, proposal.log_density(points), fixture)
    assessment["pilot_reached_epsilon"] = pilot["all_islands_reached_epsilon"]
    assessment["passed"] &= assessment["pilot_reached_epsilon"]
    counts = np.bincount(components, minlength=len(proposal.means) + 1).tolist()
    return assessment, counts


def run_study(spec, strategy):
    result = _specification(spec, strategy)
    initial_specification = copy.deepcopy(result)
    positive = []
    for fixture in frozen.FIXTURES:
        for seed in spec.positive_seeds:
            started = time.perf_counter()
            pilot_stream = np.random.SeedSequence(seed).spawn(2)[0]
            pilot = strategy.train_pilot(lambda point, threshold=None: frozen.evaluate(point, fixture, threshold),
                                         frozen.DIMENSION, frozen.EPSILON, pilot_stream,
                                         plan=copy.deepcopy(strategy.PILOT_PLAN))
            pilot["all_islands_reached_epsilon"] = check_pilot_archive(pilot, fixture, seed, strategy)
            proposal = strategy.fit_proposal(pilot["final_points"], strategy.PILOT_PLAN)
            assessment, counts = _production(pilot, proposal, fixture, seed)
            positive.append({"fixture": fixture, "seed": seed, "kind": "learned_proposal",
                             "pilot": pilot, "proposal": proposal.to_dict(), "assessment": assessment,
                             "production_attempts": frozen.PRODUCTION_ATTEMPTS,
                             "production_component_counts": counts,
                             "rng": {"pilot_spawn_key": [0], "production_spawn_key": [1]},
                             "wall_seconds": time.perf_counter() - started})
            print(f"{spec.revision_id}: {fixture} seed {seed}: ESS={assessment['importance']['ess']:.1f}, "
                  f"relative mass error={assessment['relative_mass_error']:.3f}, pass={assessment['passed']}", flush=True)
    result.update(runtime={"numpy": np.__version__, "scipy": scipy.__version__}, positive_runs=positive,
                  negative_controls=[frozen.run_negative(seed) for seed in frozen.NEGATIVE_SEEDS])
    _equal(_specification(spec, strategy), initial_specification, "specification after study")
    return assemble(result, spec, strategy)


def assemble(stored, spec, strategy):
    """Replay the frozen production proposal and recompute all derived fields."""
    result = copy.deepcopy(stored)
    for key, expected in _specification(spec, strategy).items():
        _equal(result[key], expected, key)
    positives, negatives = result["positive_runs"], result["negative_controls"]
    _equal(sorted((r["fixture"], r["seed"]) for r in positives),
           sorted((fixture, seed) for fixture in frozen.FIXTURES for seed in spec.positive_seeds), "positive run identities")
    _equal(sorted(r["seed"] for r in negatives), sorted(frozen.NEGATIVE_SEEDS), "negative run identities")
    for run in positives:
        _equal(run["kind"], "learned_proposal", "positive kind")
        _equal(run["production_attempts"], frozen.PRODUCTION_ATTEMPTS, "production count")
        _equal(run["rng"], {"pilot_spawn_key": [0], "production_spawn_key": [1]}, "production RNG")
        pilot = run["pilot"]
        pilot["all_islands_reached_epsilon"] = check_pilot_archive(pilot, run["fixture"], run["seed"], strategy)
        proposal = frozen._checked_proposal(run["proposal"], strategy.fit_proposal(pilot["final_points"], strategy.PILOT_PLAN))
        run["assessment"], run["production_component_counts"] = _production(pilot, proposal, run["fixture"], run["seed"])
    for run in negatives:
        replay = frozen.run_negative(run["seed"])
        for key in ("fixture", "kind", "production_attempts", "proposal"):
            _equal(run[key], replay[key], f"negative {key}")
        run["assessment"] = replay["assessment"]
    result["truth"] = {fixture: frozen.truth(fixture) for fixture in frozen.FIXTURES}
    result["oracle_precision_bound"] = frozen.oracle_precision_bound()
    result["checks"] = frozen.study_checks(positives, negatives)
    result["passed"] = all(result["checks"].values())
    return result


def render(result):
    lines = ["# Revised proposal validation on known synthetic targets", "",
             f"Revision `{result['revision_id']}` uses `{result['strategy_module']}` and analytic scores only.", "",
             "**All fixed synthetic checks passed.**" if result["passed"] else
             "**Fixed synthetic checks failed. All outcomes are retained.**", "",
             f"The [prior study](../../{result['prior_study']['path']}) remains separate; its JSON SHA256 is",
             f"`{result['prior_study']['sha256']}`.", "",
             "The unchanged seven-dimensional fixtures are A = [0, 1/4]^7 and",
             "B = [1/4, 1] × [5/8, 7/8]^6 (mass 1/4096; mode masses 1/4 and 3/4),",
             "plus the boundary slab [0, 1/8] × [0, 1]^6 (mass 1/8). Exact moments and",
             "the same fixed error gates are stored in JSON. Fixture derivations are in the prior report.", "",
             f"Each learned run uses {result['production_attempts']:,} independent production attempts.",
             "Pilot particles do not enter the estimates. The fixed pilot plan is:", "",
             "```json", json.dumps(result["pilot_plan"], indent=2, sort_keys=True), "```", "",
             "## Learned proposals", "",
             "| Fixture | Seed | Accepted | ESS | Relative mass error | Largest moment error | Result |",
             "|---|---:|---:|---:|---:|---:|---|"]
    failures = []
    for run in sorted(result["positive_runs"], key=lambda r: (r["fixture"], r["seed"])):
        a = run["assessment"]
        error = max(a["absolute_moment_errors"].values()) if a["absolute_moment_errors"] else None
        error_text = "unavailable" if error is None else f"{error:.4f}"
        lines.append(f"| {run['fixture']} | {run['seed']} | {a['importance']['n_accepted']} | "
                     f"{a['importance']['ess']:.1f} | {a['relative_mass_error']:.3f} | {error_text} | "
                     f"{'pass' if a['passed'] else 'fail'} |")
        failed = [name for checks in (a["usual_checks"], a["truth_checks"])
                  for name, passed in sorted(checks.items()) if not passed]
        if not a["pilot_reached_epsilon"]:
            failed.append("pilot_reached_epsilon")
        if failed:
            failures.append(f"- {run['fixture']}, seed {run['seed']}: {', '.join(failed)}.")
    if failures:
        lines += ["", "Failed learned-run checks:", "", *failures]
    lines += ["", "## Reused negative controls", "",
              "The three previously inspected proposals centered only on B are retained with their",
              "20% uniform component. These are known failure controls, not fresh validation runs.", "",
              "| Seed | A hits | B hits | ESS | Usual screens | Known-truth screens |",
              "|---|---:|---:|---:|---|---|"]
    for run in sorted(result["negative_controls"], key=lambda r: r["seed"]):
        a = run["assessment"]
        lines.append(f"| {run['seed']} | {a['mode_counts']['A']} | {a['mode_counts']['B']} | "
                     f"{a['importance']['ess']:.1f} | {'pass' if all(a['usual_checks'].values()) else 'fail'} | "
                     f"{'pass' if all(a['truth_checks'].values()) else 'fail'} |")
    lines += ["", "The reused fixtures informed this revision. Fresh random seeds do not make these",
              "independent validation of general mode discovery. The unchanged oracle precision bound",
              "applies only to its ideal exact-target proposal, not to these learned kernels.",
              "High ESS and agreement cannot exclude an unseen region. These checks establish neither",
              "biological validity nor coverage of unknown biological acceptance regions.", ""]
    return "\n".join(lines)
