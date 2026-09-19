#!/usr/bin/env python3
"""Known-truth validation of bounded resample/move proposal training.

The fixtures and error screens are fixed before the learned positive runs.
This script uses analytic scores only and never imports the biological simulator.
"""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import scipy

from bounded_proposal import BoundedGaussianMixture
from importance_sampling import importance_diagnostics, normalized_weights
from resample_move import PILOT_PLAN, fit_proposal, train_pilot

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "analysis" / "calibration" / "proposal-synthetic-validation.json"
OUT_MD = ROOT / "analysis" / "calibration" / "proposal-synthetic-validation.md"
DIMENSION = 7
EPSILON = 1.0
PRODUCTION_ATTEMPTS = 8192
POSITIVE_SEEDS = (2026092101, 2026092102, 2026092103)
NEGATIVE_SEEDS = (2026092001, 2026092002, 2026092003)
FIXTURES = ("separated_boxes", "boundary_slab")
GATES = {"minimum_accepted": 20, "minimum_ess": 200.0,
         "maximum_normalized_weight": 0.02, "maximum_relative_mcse": 0.1,
         "maximum_relative_mass_error": 0.15, "maximum_mode_mass_error": 0.075,
         "maximum_moment_error": 0.05}
BOX_A_LOW = np.zeros(DIMENSION)
BOX_A_HIGH = np.full(DIMENSION, 0.25)
BOX_B_LOW = np.array([0.25] + [0.625] * 6)
BOX_B_HIGH = np.array([1.0] + [0.875] * 6)


def _points(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != DIMENSION or not np.isfinite(points).all():
        raise ValueError("fixture points must be a finite (n, 7) array")
    return points


def membership(points, fixture):
    points = _points(points)
    inside = ((points >= 0) & (points <= 1)).all(axis=1)
    if fixture == "separated_boxes":
        return {"A": inside & ((points >= BOX_A_LOW) & (points <= BOX_A_HIGH)).all(axis=1),
                "B": inside & ((points >= BOX_B_LOW) & (points <= BOX_B_HIGH)).all(axis=1)}
    if fixture == "boundary_slab":
        return {"slab": inside & (points[:, 0] <= 0.125)}
    raise ValueError(f"unknown fixture: {fixture}")


def evaluate(point, fixture, threshold=None):
    """Complete analytic score; threshold never changes the fixture's target."""
    point = _points([point])[0]
    if np.any(point < 0) or np.any(point > 1):
        raise ValueError("analytic scoring expects a point in the unit cube")
    if fixture == "separated_boxes":
        da = np.max(np.abs(point - (BOX_A_LOW + BOX_A_HIGH) / 2) / ((BOX_A_HIGH - BOX_A_LOW) / 2))
        db = np.max(np.abs(point - (BOX_B_LOW + BOX_B_HIGH) / 2) / ((BOX_B_HIGH - BOX_B_LOW) / 2))
        distance = float(min(da, db))
    elif fixture == "boundary_slab":
        distance = float(8 * point[0])
    else:
        raise ValueError(f"unknown fixture: {fixture}")
    return {"distance": distance, "distance_lower_bound": distance, "simulator_dose_calls": 0}


def moment_values(points):
    points = _points(points)
    return {**{f"mean_x{j+1}": points[:, j] for j in range(DIMENSION)},
            "second_x1": points[:, 0] ** 2, "second_x2": points[:, 1] ** 2,
            "cross_x1_x2": points[:, 0] * points[:, 1],
            "cross_x2_x3": points[:, 1] * points[:, 2]}


def truth(fixture):
    if fixture == "separated_boxes":
        return {"mass": 1 / 4096, "mode_masses": {"A": 0.25, "B": 0.75},
                "moments": {"mean_x1": 0.5, **{f"mean_x{j}": 19 / 32 for j in range(2, 8)},
                            "second_x1": 1 / 3, "second_x2": 331 / 768,
                            "cross_x1_x2": 91 / 256, "cross_x2_x3": 109 / 256}}
    if fixture == "boundary_slab":
        return {"mass": 1 / 8, "mode_masses": {"slab": 1.0},
                "moments": {"mean_x1": 1 / 16, **{f"mean_x{j}": 0.5 for j in range(2, 8)},
                            "second_x1": 1 / 192, "second_x2": 1 / 3,
                            "cross_x1_x2": 1 / 32, "cross_x2_x3": 1 / 4}}
    raise ValueError(f"unknown fixture: {fixture}")


def summarize(points, log_q, fixture):
    """Assess independent production attempts against external analytic truth."""
    points = _points(points)
    log_q = np.asarray(log_q, dtype=float)
    modes = membership(points, fixture)
    accepted = np.logical_or.reduce(list(modes.values()))
    diagnostics = importance_diagnostics(log_q, accepted)
    gold = truth(fixture)
    relative_mass_error = diagnostics["normalizer_estimate"] / gold["mass"] - 1
    usual = {
        "accepted_count": diagnostics["n_accepted"] >= GATES["minimum_accepted"],
        "ess": diagnostics["ess"] >= GATES["minimum_ess"],
        "maximum_weight": diagnostics["max_normalized_weight"] is not None and
            diagnostics["max_normalized_weight"] <= GATES["maximum_normalized_weight"],
        "relative_mcse": diagnostics["normalizer_relative_mcse"] is not None and
            diagnostics["normalizer_relative_mcse"] <= GATES["maximum_relative_mcse"],
    }
    counts = {name: int(mask.sum()) for name, mask in modes.items()}
    if np.any(accepted):
        weights = normalized_weights(np.where(accepted, -log_q, -np.inf))
        moments = {name: float(weights @ values) for name, values in moment_values(points).items()}
        mode_masses = {name: float(weights[mask].sum()) for name, mask in modes.items()}
        moment_errors = {name: abs(value - gold["moments"][name]) for name, value in moments.items()}
        mode_errors = {name: abs(value - gold["mode_masses"][name]) for name, value in mode_masses.items()}
    else:
        moments = mode_masses = moment_errors = mode_errors = None
    gold_checks = {
        "normalizing_mass": abs(relative_mass_error) <= GATES["maximum_relative_mass_error"],
        "every_region_observed": all(n > 0 for n in counts.values()),
        "mode_masses": mode_errors is not None and max(mode_errors.values()) <= GATES["maximum_mode_mass_error"],
        "moments": moment_errors is not None and max(moment_errors.values()) <= GATES["maximum_moment_error"],
    }
    return {"importance": diagnostics, "mode_counts": counts, "mode_masses": mode_masses,
            "moments": moments, "absolute_moment_errors": moment_errors,
            "absolute_mode_mass_errors": mode_errors, "relative_mass_error": relative_mass_error,
            "usual_checks": usual, "truth_checks": gold_checks,
            "passed": all(usual.values()) and all(gold_checks.values())}


def oracle_precision_bound():
    """Hoeffding/union bound for an ideal exact-box proposal, not learned kernels."""
    n, lower_k = PRODUCTION_ATTEMPTS, 6000
    beta = 0.8 + 0.2 / 4096
    terms = {
        "fewer_than_6000_hits": math.exp(-2 * n * (beta - lower_k / n) ** 2),
        "normalizer_error": 2 * math.exp(-2 * n * (GATES["maximum_relative_mass_error"] * beta) ** 2),
        "eleven_moment_errors": 22 * math.exp(-2 * lower_k * GATES["maximum_moment_error"] ** 2),
        "mode_mass_error": 2 * math.exp(-2 * lower_k * GATES["maximum_mode_mass_error"] ** 2),
    }
    return {"proposal": "0.2 Uniform(cube) + 0.8 Uniform(A union B)",
            "accepted_density": 3277.0, "acceptance_probability": beta,
            "expected_accepted": n * beta, "conditioning_hit_count": lower_k,
            "bound_terms": terms, "truth_error_probability_bound_per_run": sum(terms.values()),
            "scope": "analytic precision benchmark only; not a guarantee for learned proposals or ESS screens"}


def run_positive(fixture, seed, plan=None, production_attempts=PRODUCTION_ATTEMPTS):
    pilot_stream, production_stream = np.random.SeedSequence(seed).spawn(2)
    started = time.perf_counter()
    pilot = train_pilot(lambda point, threshold=None: evaluate(point, fixture, threshold),
                        DIMENSION, EPSILON, pilot_stream, plan=plan)
    proposal = fit_proposal(pilot["final_points"], pilot["plan"])
    points, components = proposal.sample(np.random.default_rng(production_stream), production_attempts)
    assessment = summarize(points, proposal.log_density(points), fixture)
    assessment["pilot_reached_epsilon"] = pilot["all_islands_reached_epsilon"]
    assessment["passed"] = assessment["passed"] and assessment["pilot_reached_epsilon"]
    return {"fixture": fixture, "seed": seed, "kind": "learned_proposal",
            "pilot": pilot, "proposal": proposal.to_dict(), "assessment": assessment,
            "production_attempts": production_attempts,
            "production_component_counts": np.bincount(components, minlength=len(proposal.means) + 1).tolist(),
            "rng": {"pilot_spawn_key": list(pilot_stream.spawn_key),
                    "production_spawn_key": list(production_stream.spawn_key)},
            "wall_seconds": time.perf_counter() - started}


def run_negative(seed, production_attempts=PRODUCTION_ATTEMPTS):
    # Deliberately omit box A from the fitted component. The defensive uniform
    # retains support, but at this budget it usually supplies no A observations.
    center = (BOX_B_LOW + BOX_B_HIGH) / 2
    scale = (BOX_B_HIGH - BOX_B_LOW) / math.sqrt(12)
    proposal = BoundedGaussianMixture(DIMENSION, 0.2, [center], [scale])
    points, _ = proposal.sample(np.random.default_rng(seed), production_attempts)
    return {"fixture": "separated_boxes", "seed": seed, "kind": "known_missing_region",
            "proposal": proposal.to_dict(), "production_attempts": production_attempts,
            "assessment": summarize(points, proposal.log_density(points), "separated_boxes")}


def source_hashes():
    return {f"scripts/{name}": hashlib.sha256((ROOT / "scripts" / name).read_bytes()).hexdigest()
            for name in ("proposal_synthetic_validation.py", "resample_move.py",
                         "bounded_proposal.py", "importance_sampling.py")}


def study_checks(positive, negative):
    return {
        "all_learned_runs_passed": all(r["assessment"]["passed"] for r in positive),
        "all_negative_controls_rejected_by_truth": all(not all(r["assessment"]["truth_checks"].values()) for r in negative),
        "false_assurance_negative_control_observed": any(all(r["assessment"]["usual_checks"].values()) and
            not all(r["assessment"]["truth_checks"].values()) for r in negative),
    }


def _require_equal(actual, expected, label):
    if actual != expected:
        raise ValueError(f"stored {label} differs from the fixed study specification")


def _require_close(actual, expected, label):
    actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=1e-12, atol=1e-14):
        raise ValueError(f"stored {label} does not reproduce from its archived inputs")


def _checked_proposal(stored, reconstructed):
    expected = reconstructed.to_dict()
    _require_equal(set(stored), set(expected), "proposal fields")
    for key in ("dimension", "uniform_weight"):
        _require_equal(stored[key], expected[key], f"proposal {key}")
    for key in ("means", "scales"):
        _require_close(stored[key], expected[key], f"proposal {key}")
    # Replay the exact archived density after checking its independently fitted
    # reconstruction. Tolerance accommodates harmless arithmetic differences.
    return BoundedGaussianMixture.from_dict(stored)


def _check_pilot(pilot, fixture, seed):
    """Check final conditioning and trace consistency without rerunning MH."""
    _require_equal(pilot["plan"], PILOT_PLAN, "pilot plan")
    islands = pilot["islands"]
    _require_equal([i["island"] for i in islands], list(range(PILOT_PLAN["islands"])), "pilot islands")
    streams = np.random.SeedSequence(seed).spawn(2)[0].spawn(PILOT_PLAN["islands"])
    n, levels, moves = (PILOT_PLAN[k] for k in ("particles_per_island", "levels", "moves_per_level"))
    reached = []
    for island, stream in zip(islands, streams):
        _require_equal(island["spawn_key"], list(stream.spawn_key), "island RNG")
        initial = np.random.default_rng(stream).uniform(size=(n, DIMENSION))
        _require_close(island["initial_points"], initial, "initial pilot points")
        final = _points(island["final_points"])
        if final.shape != (n, DIMENSION):
            raise ValueError("stored final pilot cloud has the wrong size")
        for name, points in (("initial", initial), ("final", final)):
            scores = [evaluate(point, fixture)["distance"] for point in points]
            _require_close(island[f"{name}_distances"], scores, f"{name} pilot scores")
        history = island["history"]
        _require_equal([h["level"] for h in history], list(range(1, levels + 1)), "pilot levels")
        thresholds = np.array([h["threshold"] for h in history], dtype=float)
        if not np.isfinite(thresholds).all() or np.any(thresholds < EPSILON) or np.any(np.diff(thresholds) > 0):
            raise ValueError("stored pilot thresholds violate the fixed final threshold")
        if np.any(np.asarray(island["final_distances"]) > thresholds[-1]):
            raise ValueError("stored final pilot points escaped their conditioning region")
        hits = np.flatnonzero(thresholds == EPSILON)
        first = int(hits[0] + 1) if len(hits) else None
        _require_equal(island["first_epsilon_level"], first, "first target level")
        reached.append(first is not None)
        _require_equal(island["initial_simulator_dose_calls"], 0, "analytic initial cost")
        for step in history:
            _require_equal(step["attempts"], n * moves, "pilot move attempts")
            _require_equal(step["simulator_dose_calls"], 0, "analytic move cost")
    for name in ("final_points", "final_distances"):
        _require_close(pilot[name], np.concatenate([i[name] for i in islands]), f"combined {name}")
    _require_equal(pilot["attempts"], len(islands) * (n + levels * n * moves), "total pilot attempts")
    _require_equal(pilot["simulator_dose_calls"], 0, "analytic total cost")
    pilot["all_islands_reached_epsilon"] = all(reached)


def assemble(stored):
    """Replay frozen production; never trust cached assessment fields.

    Pilot histories are checked for consistency, not certified by rerunning
    adaptation. Runtime versions and wall times remain historical measurements.
    """
    result = copy.deepcopy(stored)
    for key, expected in (("schema_version", 1), ("dimension", DIMENSION), ("epsilon", EPSILON),
                          ("source_hashes", source_hashes()),
                          ("pilot_plan", PILOT_PLAN), ("gates", GATES),
                          ("production_attempts", PRODUCTION_ATTEMPTS),
                          ("positive_seeds", list(POSITIVE_SEEDS)),
                          ("negative_seeds", list(NEGATIVE_SEEDS))):
        _require_equal(result[key], expected, key)
    positive, negative = result["positive_runs"], result["negative_controls"]
    _require_equal(sorted((r["fixture"], r["seed"]) for r in positive),
                   sorted((fixture, seed) for fixture in FIXTURES for seed in POSITIVE_SEEDS), "positive run identities")
    _require_equal(sorted(r["seed"] for r in negative), sorted(NEGATIVE_SEEDS), "negative run identities")
    for run in positive:
        _require_equal(run["kind"], "learned_proposal", "positive run kind")
        _require_equal(run["production_attempts"], PRODUCTION_ATTEMPTS, "production count")
        _require_equal(run["rng"], {"pilot_spawn_key": [0], "production_spawn_key": [1]}, "production RNG")
        _check_pilot(run["pilot"], run["fixture"], run["seed"])
        proposal = _checked_proposal(run["proposal"], fit_proposal(run["pilot"]["final_points"], PILOT_PLAN))
        stream = np.random.SeedSequence(run["seed"]).spawn(2)[1]
        points, components = proposal.sample(np.random.default_rng(stream), PRODUCTION_ATTEMPTS)
        run["assessment"] = summarize(points, proposal.log_density(points), run["fixture"])
        run["assessment"]["pilot_reached_epsilon"] = run["pilot"]["all_islands_reached_epsilon"]
        run["assessment"]["passed"] &= run["assessment"]["pilot_reached_epsilon"]
        run["production_component_counts"] = np.bincount(components, minlength=len(proposal.means) + 1).tolist()
    for run in negative:
        _require_equal(run["fixture"], "separated_boxes", "negative fixture")
        _require_equal(run["kind"], "known_missing_region", "negative run kind")
        _require_equal(run["production_attempts"], PRODUCTION_ATTEMPTS, "negative production count")
        replay = run_negative(run["seed"], PRODUCTION_ATTEMPTS)
        _checked_proposal(run["proposal"], BoundedGaussianMixture.from_dict(replay["proposal"]))
        run["assessment"] = replay["assessment"]
    result["truth"] = {name: truth(name) for name in FIXTURES}
    result["oracle_precision_bound"] = oracle_precision_bound()
    result["checks"] = study_checks(positive, negative)
    result["passed"] = all(result["checks"].values())
    return result


def run_study():
    hashes = source_hashes()
    positive = []
    for fixture in FIXTURES:
        for seed in POSITIVE_SEEDS:
            result = run_positive(fixture, seed)
            positive.append(result)
            d = result["assessment"]
            print(f"{fixture} seed {seed}: ESS={d['importance']['ess']:.1f}, "
                  f"relative mass error={d['relative_mass_error']:.3f}, pass={d['passed']}", flush=True)
    negative = [run_negative(seed) for seed in NEGATIVE_SEEDS]
    if source_hashes() != hashes:
        raise RuntimeError("validation sources changed during the experiment")
    checks = study_checks(positive, negative)
    return {"schema_version": 1, "dimension": DIMENSION, "epsilon": EPSILON, "source_hashes": hashes,
            "runtime": {"numpy": np.__version__, "scipy": scipy.__version__},
            "pilot_plan": dict(PILOT_PLAN), "production_attempts": PRODUCTION_ATTEMPTS,
            "positive_seeds": list(POSITIVE_SEEDS), "negative_seeds": list(NEGATIVE_SEEDS),
            "gates": dict(GATES), "truth": {name: truth(name) for name in FIXTURES},
            "oracle_precision_bound": oracle_precision_bound(),
            "positive_runs": positive, "negative_controls": negative,
            "checks": checks, "passed": all(checks.values())}


def render(result):
    lines = ["# Proposal validation on known synthetic targets", "",
             "Generated by `scripts/proposal_synthetic_validation.py`; no biological simulator is used.", "",
             "**All prespecified synthetic checks passed.**" if result["passed"] else
             "**Prespecified synthetic checks failed. Failures are retained.**", "",
             "## Fixtures and fixed criteria", "",
             "The seven-dimensional prior is uniform on the unit cube. The separated target is",
             "A = [0, 1/4]^7 and B = [1/4, 1] × [5/8, 7/8]^6. Its mass is 1/4096,",
             "with conditional mode masses 1/4 and 3/4. E[X1] = 1/2 and E[Xj] = 19/32",
             "for j ≥ 2. Second and cross moments are stored as independent analytic truth in JSON.", "",
             "The boundary slab [0, 1/8] × [0, 1]^6 has mass 1/8; its other six coordinates",
             "remain uniform. It checks boundary behavior and unnecessary contraction of uninformed directions.", "",
             f"Each learned run uses the same two-island pilot plan and {result['production_attempts']:,}",
             "independent production attempts from a frozen exact-density proposal. Pilot particles",
             "and resampling copies never enter the estimates. Each of three fixed seeds is run on both fixtures.", "",
             "Required: at least 20 accepted draws, ESS ≥ 200, maximum normalized weight ≤ 0.02,",
             "reported relative MCSE ≤ 0.10, absolute relative mass error ≤ 0.15, mode-mass error ≤ 0.075,",
             "and all eleven mean/second/cross-moment errors ≤ 0.05. Every accepted region must be observed",
             "and both pilot islands must reach the unchanged final threshold.", "",
             "## Learned proposals", "",
             "| Fixture | Seed | Accepted | ESS | Relative mass error | Largest moment error | Pilot reached target | Result |",
             "|---|---:|---:|---:|---:|---:|---|---|"]
    for run in sorted(result["positive_runs"], key=lambda r: (r["fixture"], r["seed"])):
        a = run["assessment"]
        error = max(a["absolute_moment_errors"].values()) if a["absolute_moment_errors"] else None
        error_text = f"{error:.4f}" if error is not None else "unavailable"
        lines.append(f"| {run['fixture']} | {run['seed']} | {a['importance']['n_accepted']} | "
                     f"{a['importance']['ess']:.1f} | {a['relative_mass_error']:.3f} | {error_text} | "
                     f"{a['pilot_reached_epsilon']} | {'pass' if a['passed'] else 'fail'} |")
    failures = []
    for run in sorted(result["positive_runs"], key=lambda r: (r["fixture"], r["seed"])):
        a = run["assessment"]
        failed = [name for checks in (a["usual_checks"], a["truth_checks"])
                  for name, passed in sorted(checks.items()) if not passed]
        if not a["pilot_reached_epsilon"]:
            failed.append("pilot_reached_epsilon")
        if failed:
            failures.append(f"- {run['fixture']}, seed {run['seed']}: {', '.join(failed)}.")
    if failures:
        lines += ["", "Failed learned-run checks:", "", *failures]
    lines += ["", "## Required negative control: an omitted region", "",
              "These proposals contain a 20% uniform component and a product truncated-normal component",
              "centered only on B. Its coordinate scales are B's widths divided by sqrt(12). The uniform",
              "component gives support to A but does not guarantee a hit. With 8,192 attempts, the chance",
              "of no uniform-component observation in A is about 90.5%; Gaussian leakage into A is negligible.", "",
              "| Seed | A hits | B hits | ESS | Relative mass error | Usual weight screens | Known-truth screens |",
              "|---|---:|---:|---:|---:|---|---|"]
    for run in sorted(result["negative_controls"], key=lambda r: r["seed"]):
        a = run["assessment"]
        lines.append(f"| {run['seed']} | {a['mode_counts']['A']} | {a['mode_counts']['B']} | "
                     f"{a['importance']['ess']:.1f} | {a['relative_mass_error']:.3f} | "
                     f"{'pass' if all(a['usual_checks'].values()) else 'fail'} | "
                     f"{'pass' if all(a['truth_checks'].values()) else 'fail'} |")
    bound = result["oracle_precision_bound"]
    lines += ["", "## Analytic calibration and limits", "",
              "An ideal comparison proposal is 0.2 Uniform(cube) + 0.8 Uniform(A union B).",
              "Its accepted density is exactly 3277, and its expected accepted count is 6,554 of 8,192.",
              "Conditional on at least 6,000 hits, accepted points are independent draws from the exact",
              "target. Hoeffding bounds for the bounded moments and a union bound, plus the binomial",
              f"normalizer bound, give a truth-error-screen failure bound of {bound['truth_error_probability_bound_per_run']:.2e}",
              "per ideal run. This calibrates an attainable precision benchmark; it is not a guarantee",
              "for the learned kernels, their ESS, or mode discovery.", "",
              "The truth-error thresholds and seeds were fixed before the learned positive experiments.",
              "Negative-control seeds were inspected during fixture design; all three are retained.",
              "Three successful runs on a fixture would still give weak evidence about its general failure rate.",
              "This finite validation cannot establish coverage of unknown biological acceptance regions.",
              "ESS and agreement among runs cannot replace a known-truth mode check. No biological posterior",
              "or spatial headline is validated by these synthetic results.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    if args.render_only:
        result = assemble(json.loads(OUT_JSON.read_text()))
    else:
        result = assemble(run_study())
    json_text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    md_text = render(result)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"wrote {OUT_MD}; all checks passed: {result['passed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
