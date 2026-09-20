#!/usr/bin/env python3
"""Prospective geometry challenges for the unchanged local-move sampler.

Learned proposals use no oracle information. Direct-target controls separately
test achievable precision and the inability of ESS to certify full support.
"""
import argparse
import copy
import gzip
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from bounded_proposal import BoundedGaussianMixture
import coverage_challenge_geometry as geometry
from coverage_pilot_audit import validate_pilot
from importance_sampling import importance_diagnostics, normalized_weights
import resample_move_local as strategy

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "analysis/calibration/proposal-coverage-challenges.json"
OUT_MD = ROOT / "analysis/calibration/proposal-coverage-challenges.md"
ARCHIVES = ROOT / "analysis/calibration/coverage-challenges-20260920"
PROTOCOL = "docs/COVERAGE_CHALLENGE_PLAN.md"
REVISION = "unseen-geometry-20260920"
ATTEMPTS = 8192
EPSILON = 1.0
LEARNED_SEEDS = {
    "rotated_box": (2026092501, 2026092502, 2026092503),
    "annular_cylinder": (2026092511, 2026092512, 2026092513),
    "unequal_balls": (2026092521, 2026092522, 2026092523),
}
ORACLE_SEEDS = {name: tuple(seed + 100 for seed in seeds)
                for name, seeds in LEARNED_SEEDS.items()}
NEGATIVE_SEEDS = {name: tuple(seed + 200 for seed in seeds)
                  for name, seeds in LEARNED_SEEDS.items()}
GATES = {
    "minimum_accepted": 20, "minimum_ess": 200.0,
    "maximum_normalized_weight": 0.02, "maximum_relative_mcse": 0.1,
    "maximum_relative_mass_error": 0.15, "maximum_region_mass_error": 0.075,
    "maximum_moment_error": 0.05,
}
SOURCE_PATHS = (
    "scripts/proposal_coverage_challenges.py",
    "scripts/coverage_challenge_geometry.py", "scripts/coverage_pilot_audit.py",
    "scripts/resample_move_local.py", "scripts/resample_move.py",
    "scripts/bounded_proposal.py", "scripts/importance_sampling.py", PROTOCOL,
)
PRIOR_PATH = "analysis/calibration/proposal-synthetic-validation-v2.json"
PRIOR_SHA256 = "cc801aef90938127a15d7a671a88953f081d19b1f11c2d4a59b73c0da7eef46e"


def _equal(actual, expected, label):
    """JSON-level equality also distinguishes bool/int and int/float."""
    if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(
            expected, sort_keys=True, allow_nan=False):
        raise ValueError(f"{label} differs from the fixed study")


def _close(actual, expected, label):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if (actual.shape != expected.shape or not np.isfinite(actual).all() or
            not np.isfinite(expected).all() or
            not np.allclose(actual, expected, rtol=1e-12, atol=1e-14)):
        raise ValueError(f"{label} does not reproduce")


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def specification():
    """The checked-in protocol chooses all inputs, never an archive."""
    prior_path = ROOT / PRIOR_PATH
    if _hash(prior_path) != PRIOR_SHA256:
        raise ValueError("the earlier frozen synthetic study changed")
    prior = json.loads(prior_path.read_text())
    for path in SOURCE_PATHS[3:7]:
        if _hash(ROOT / path) != prior["source_hashes"][path]:
            raise ValueError(f"the frozen sampler source changed: {path}")
    return {
        "revision": REVISION, "dimension": geometry.DIMENSION,
        "epsilon": EPSILON, "production_attempts": ATTEMPTS,
        "pilot_plan": copy.deepcopy(strategy.PILOT_PLAN),
        "learned_seeds": {name: list(seeds) for name, seeds in LEARNED_SEEDS.items()},
        "oracle_seeds": {name: list(seeds) for name, seeds in ORACLE_SEEDS.items()},
        "negative_seeds": {name: list(seeds) for name, seeds in NEGATIVE_SEEDS.items()},
        "gates": dict(GATES),
        "source_hashes": {path: _hash(ROOT / path) for path in SOURCE_PATHS},
        "prior_study": {"path": PRIOR_PATH, "sha256": PRIOR_SHA256},
    }


def summarize(points, log_q, fixture):
    points = np.asarray(points, dtype=float)
    log_q = np.asarray(log_q, dtype=float)
    if (points.shape != (ATTEMPTS, geometry.DIMENSION) or
            not np.isfinite(points).all() or np.any((points < 0) | (points > 1)) or
            log_q.shape != (ATTEMPTS,) or not np.isfinite(log_q).all()):
        raise ValueError("production must contain the fixed number of finite cube points and densities")
    regions = geometry.membership(points, fixture)
    region_hits = np.stack(list(regions.values()))
    if np.any(region_hits.sum(axis=0) > 1):
        raise ValueError("fixture regions overlap")
    accepted = region_hits.any(axis=0)
    distance_acceptance = geometry.scores(points, fixture) <= EPSILON
    if not np.array_equal(accepted, distance_acceptance):
        raise ValueError("analytic scores and region membership disagree")
    diagnostics = importance_diagnostics(log_q, accepted)
    gold = geometry.truth(fixture)
    counts = {name: int(mask.sum()) for name, mask in regions.items()}
    relative_error = diagnostics["normalizer_estimate"] / gold["mass"] - 1
    usual = {
        "accepted_count": diagnostics["n_accepted"] >= GATES["minimum_accepted"],
        "ess": diagnostics["ess"] >= GATES["minimum_ess"],
        "maximum_weight": diagnostics["max_normalized_weight"] is not None and
            diagnostics["max_normalized_weight"] <= GATES["maximum_normalized_weight"],
        "relative_mcse": diagnostics["normalizer_relative_mcse"] is not None and
            diagnostics["normalizer_relative_mcse"] <= GATES["maximum_relative_mcse"],
    }
    if np.any(accepted):
        weights = normalized_weights(-log_q[accepted])
        features = geometry.features(points[accepted], fixture)
        if set(features) != set(gold["moments"]):
            raise ValueError("geometry features and analytic moments differ")
        for values in features.values():
            if (np.shape(values) != weights.shape or not np.isfinite(values).all() or
                    np.any((values < -1e-12) | (values > 1 + 1e-12))):
                raise ValueError("geometry features must be bounded on the target")
        moments = {name: float(weights @ values) for name, values in features.items()}
        masses = {name: float(weights[mask[accepted]].sum()) for name, mask in regions.items()}
        moment_errors = {name: abs(value - gold["moments"][name]) for name, value in moments.items()}
        region_errors = {name: abs(value - gold["region_masses"][name]) for name, value in masses.items()}
    else:
        moments = masses = moment_errors = region_errors = None
    truth_checks = {
        "normalizing_mass": abs(relative_error) <= GATES["maximum_relative_mass_error"],
        "every_region_observed": all(counts.values()),
        "region_masses": region_errors is not None and
            max(region_errors.values()) <= GATES["maximum_region_mass_error"],
        "moments": moment_errors is not None and
            max(moment_errors.values()) <= GATES["maximum_moment_error"],
    }
    return {"importance": diagnostics, "region_counts": counts, "region_masses": masses,
            "moments": moments, "absolute_moment_errors": moment_errors,
            "absolute_region_mass_errors": region_errors, "relative_mass_error": relative_error,
            "usual_checks": usual, "truth_checks": truth_checks,
            "passed": all(usual.values()) and all(truth_checks.values())}


def archive_path(directory, fixture, seed):
    if fixture not in LEARNED_SEEDS or seed not in LEARNED_SEEDS[fixture]:
        raise ValueError("unplanned fixture or seed")
    return directory / f"{fixture}-{seed}.json.gz"


def run_learned(fixture, seed):
    archive_path(ARCHIVES, fixture, seed)
    spec = specification()
    started = time.perf_counter()
    pilot_stream, production_stream = np.random.SeedSequence(seed).spawn(2)
    pilot = strategy.train_pilot(
        lambda point, threshold=None: geometry.evaluate(point, fixture, threshold),
        geometry.DIMENSION, EPSILON, pilot_stream, plan=copy.deepcopy(strategy.PILOT_PLAN))
    reached = validate_pilot(pilot, fixture, seed, geometry.evaluate, EPSILON)
    proposal = strategy.fit_proposal(pilot["final_points"], strategy.PILOT_PLAN)
    rng = np.random.default_rng(production_stream)
    points, ids = proposal.sample(rng, ATTEMPTS)
    scores = geometry.scores(points, fixture)
    result = {
        "schema_version": 1, "fixture": fixture, "seed": seed, "specification": spec,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "scipy": scipy.__version__, "bit_generator": type(rng.bit_generator).__name__},
        "rng": {"pilot_spawn_key": list(pilot_stream.spawn_key),
                "production_spawn_key": list(production_stream.spawn_key)},
        "pilot": pilot, "proposal": proposal.to_dict(),
        "production": {"unit_points": points.tolist(), "component_ids": ids.tolist(),
                       "log_q": proposal.log_density(points).tolist(),
                       "scores": scores.tolist(), "accepted": (scores <= EPSILON).tolist()},
        "wall_seconds": time.perf_counter() - started,
    }
    _equal(specification(), spec, "sources after learned run")
    assessment = assess_archive(result, spec)
    print(f"{fixture}, seed {seed}: ESS={assessment['importance']['ess']:.1f}, "
          f"mass error={assessment['relative_mass_error']:+.3f}, "
          f"pilot reached={reached}, pass={assessment['passed']}", flush=True)
    return result


def assess_archive(raw, spec=None):
    spec = specification() if spec is None else spec
    _equal(raw["schema_version"], 1, "archive schema")
    fixture, seed = raw["fixture"], raw["seed"]
    if type(seed) is not int:
        raise ValueError("archive seed must be an integer")
    archive_path(ARCHIVES, fixture, seed)
    _equal(raw["specification"], spec, "archive specification")
    _equal(raw["rng"], {"pilot_spawn_key": [0], "production_spawn_key": [1]}, "streams")
    runtime = raw["runtime"]
    if (not isinstance(runtime, dict) or
            set(runtime) != {"python", "numpy", "scipy", "bit_generator"} or
            any(type(value) is not str or not value.strip() for value in runtime.values())):
        raise ValueError("runtime must contain nonempty Python, NumPy, SciPy and bit-generator labels")
    if (isinstance(raw["wall_seconds"], bool) or not isinstance(raw["wall_seconds"], (int, float)) or
            not math.isfinite(raw["wall_seconds"]) or raw["wall_seconds"] < 0):
        raise ValueError("invalid elapsed time")
    reached = validate_pilot(raw["pilot"], fixture, seed, geometry.evaluate, EPSILON)
    expected = strategy.fit_proposal(raw["pilot"]["final_points"], strategy.PILOT_PLAN)
    actual = raw["proposal"]
    _equal(sorted(actual), sorted(expected.to_dict()), "proposal fields")
    _equal(actual["dimension"], geometry.DIMENSION, "proposal dimension")
    _equal(actual["uniform_weight"], strategy.PILOT_PLAN["uniform_weight"], "uniform mass")
    for key in ("means", "scales"):
        _close(actual[key], expected.to_dict()[key], f"proposal {key}")
    proposal = BoundedGaussianMixture.from_dict(actual)
    stream = np.random.SeedSequence(seed).spawn(2)[1]
    rng = np.random.default_rng(stream)
    _equal(raw["runtime"]["bit_generator"], type(rng.bit_generator).__name__, "bit generator")
    points, ids = proposal.sample(rng, ATTEMPTS)
    production = raw["production"]
    _close(production["unit_points"], points, "production points")
    _equal(production["component_ids"], ids.tolist(), "component IDs")
    log_q = proposal.log_density(np.asarray(production["unit_points"]))
    _close(production["log_q"], log_q, "production log densities")
    scores = geometry.scores(np.asarray(production["unit_points"]), fixture)
    _close(production["scores"], scores, "analytic production scores")
    _equal(production["accepted"], (scores <= EPSILON).tolist(), "acceptance decisions")
    result = summarize(np.asarray(production["unit_points"]), log_q, fixture)
    result["pilot_reached_epsilon"] = reached
    result["passed"] = result["passed"] and reached
    return result


def run_control(fixture, seed, restricted):
    rng = np.random.default_rng(seed)
    points, log_q = geometry.oracle_sample(rng, ATTEMPTS, fixture, restricted=restricted)
    return {
        "fixture": fixture, "seed": seed,
        "kind": "support_hole_oracle" if restricted else "full_target_oracle",
        "retained_mass_fraction": geometry.oracle_retained_fraction(fixture, restricted),
        "assessment": summarize(points, log_q, fixture),
    }


def _archive_hashes(directory):
    expected = [archive_path(directory, name, seed)
                for name, seeds in LEARNED_SEEDS.items() for seed in seeds]
    if set(directory.glob("*.json.gz")) != set(expected):
        raise ValueError("the archive set must contain every planned run and no extras")
    return {path.name: _hash(path) for path in sorted(expected)}


def build_report(directory=ARCHIVES):
    spec = specification()
    archive_hashes = _archive_hashes(directory)
    runs, runtime = [], None
    for fixture, seeds in LEARNED_SEEDS.items():
        for seed in seeds:
            path = archive_path(directory, fixture, seed)
            raw = json.loads(gzip.decompress(path.read_bytes()))
            _equal([raw["fixture"], raw["seed"]], [fixture, seed], "archive identity")
            if runtime is None:
                runtime = raw["runtime"]
            _equal(raw["runtime"], runtime, "run runtimes")
            assessment = assess_archive(raw, spec)
            runs.append({"fixture": fixture, "seed": seed, "archive": path.name,
                         "wall_seconds": raw["wall_seconds"], "assessment": assessment,
                         "pilot_attempts": raw["pilot"]["attempts"],
                         "islands": [
                             {key: island[key] for key in ("island", "first_epsilon_level", "history")}
                             for island in raw["pilot"]["islands"]]})
    oracles = [run_control(name, seed, False)
               for name, seeds in ORACLE_SEEDS.items() for seed in seeds]
    negatives = [run_control(name, seed, True)
                 for name, seeds in NEGATIVE_SEEDS.items() for seed in seeds]
    checks = {
        "all_nine_learned_runs_pass": all(run["assessment"]["passed"] for run in runs),
        "all_nine_full_target_oracles_pass": all(run["assessment"]["passed"] for run in oracles),
        "all_nine_support_holes_pass_usual_screens":
            all(all(run["assessment"]["usual_checks"].values()) for run in negatives),
        "all_nine_support_holes_fail_truth":
            all(not all(run["assessment"]["truth_checks"].values()) for run in negatives),
    }
    _equal(specification(), spec, "sources after report reconstruction")
    _equal(_archive_hashes(directory), archive_hashes, "archives after report reconstruction")
    return {"schema_version": 1, "specification": spec, "runtime": runtime,
            "runtime_scope": "Labels describe archived learned runs; controls are regenerated in the reader's environment.",
            "archive_sha256": archive_hashes,
            "truth": {name: geometry.truth(name) for name in LEARNED_SEEDS},
            "positive_runs": runs, "oracle_controls": oracles, "negative_controls": negatives,
            "checks": checks, "passed": all(checks.values())}


def assemble(stored):
    _equal(stored["specification"], specification(), "report specification")
    _equal(stored["archive_sha256"], _archive_hashes(ARCHIVES), "report archive hashes")
    return build_report(ARCHIVES)


def render(result):
    lines = [
        "# Coverage challenges on new synthetic geometries", "",
        "**All prespecified checks passed.**" if result["passed"] else
        "**Prespecified checks failed. All planned outcomes are retained.**", "",
        "The local-move sampler is unchanged. These geometries and all seeds, budgets",
        "and criteria were committed before evaluation. They were not used to develop",
        "the sampler, but were deliberately designed with knowledge of its earlier",
        "failure modes. This is a bounded challenge study, not an unbiased sample of",
        "possible targets or independent biological validation.", "",
        "See the [frozen protocol](../../docs/COVERAGE_CHALLENGE_PLAN.md) for exact",
        "geometry, analytic masses, moments, partition rules and complete seed lists.",
        "Each learned run uses 29,184 pilot attempts and 8,192 independent production",
        "attempts. No pilot particles contribute to estimates; no model dose calls",
        "are made. All proposal weights use the full normalized mixture.", "",
        "## Learned proposals", "",
        "| Geometry | Seed | Accepted | ESS | Relative mass error | Maximum region error | Maximum moment error | Result |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    failures = []
    for run in sorted(result["positive_runs"], key=lambda r: (r["fixture"], r["seed"])):
        a = run["assessment"]
        def maximum(values):
            return f"{max(values.values()):.4f}" if values is not None else "unavailable"
        lines.append(
            f"| {run['fixture']} | {run['seed']} | {a['importance']['n_accepted']} | "
            f"{a['importance']['ess']:.1f} | {a['relative_mass_error']:+.4f} | "
            f"{maximum(a['absolute_region_mass_errors'])} | "
            f"{maximum(a['absolute_moment_errors'])} | {'pass' if a['passed'] else 'fail'} |")
        failed = [name for group in ("usual_checks", "truth_checks")
                  for name, passed in sorted(a[group].items()) if not passed]
        if not a["pilot_reached_epsilon"]:
            failed.append("pilot_reached_epsilon")
        if failed:
            failures.append(f"- {run['fixture']}, seed {run['seed']}: {', '.join(failed)}.")
    lines += ["", "Failed checks:", "", *(failures or ["None."]), "",
              "## Analytic controls", "",
              "Full-target oracles draw directly from known analytic targets. Support-hole",
              "oracles deliberately omit half of the rotated or annular target, or the",
              "smaller ball. Their densities are correct only on their restricted support.",
              "They are intentionally invalid importance proposals for the full target;",
              "they show that ordinary weight diagnostics cannot establish support.",
              "They are distinct from the earlier defensive-mixture failure controls.",
              "No oracle points or identities enter learned pilots or production.", "",
              "| Geometry | Seed | Control | ESS | Usual screens | Truth screens |",
              "|---|---:|---|---:|---|---|"]
    for run in sorted(result["oracle_controls"] + result["negative_controls"],
                      key=lambda r: (r["fixture"], r["kind"], r["seed"])):
        a = run["assessment"]
        lines.append(
            f"| {run['fixture']} | {run['seed']} | {run['kind']} | {a['importance']['ess']:.1f} | "
            f"{'pass' if all(a['usual_checks'].values()) else 'fail'} | "
            f"{'pass' if all(a['truth_checks'].values()) else 'fail'} |")
    lines += ["", "## Scope and reproducibility", ""]
    for key, value in sorted(result["checks"].items()):
        lines.append(f"- {key}: {'pass' if value else 'fail'}.")
    lines += [
        "", "Archives retain pilot endpoints and histories, fitted kernels, independent",
        "production points, component identities, densities, scores and decisions.",
        "Reassembly checks source and archive hashes, replays production randomness,",
        "recomputes analytic endpoint scores and reconciles pilot counts. It does not",
        "replay unarchived intermediate pilot trajectories. Numerical assessment",
        "roundoff can vary across platforms; counts and pass/fail decisions remain exact.",
        "Runtime labels describe the archived learned runs; reconstructed analytic",
        "controls execute in the reader's environment rather than the archived runtime.",
        "", "The controls establish only that the fixed checks can detect specified",
        "support defects and that their accuracy is achievable with analytic oracles.",
        "Passing learned runs cannot exclude a region outside these challenges.",
        "A failed challenge does not alter the historical biological scores or prove",
        "that the biological target shares that geometry. No sampler tuning, pooled",
        "posterior, assay validation, or clinical claim is inferred from this study.", "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--resume", action="store_true",
                        help="validate and retain existing planned archives; never overwrite")
    args = parser.parse_args()
    if args.render_only:
        result = assemble(json.loads(OUT_JSON.read_text()))
    else:
        ARCHIVES.mkdir(parents=True, exist_ok=True)
        for fixture, seeds in LEARNED_SEEDS.items():
            for seed in seeds:
                path = archive_path(ARCHIVES, fixture, seed)
                if path.exists():
                    if not args.resume:
                        raise FileExistsError(f"{path.name} already exists; use --resume to retain it")
                    raw = json.loads(gzip.decompress(path.read_bytes()))
                    _equal([raw["fixture"], raw["seed"]], [fixture, seed], "existing archive identity")
                    assess_archive(raw)
                    print(f"Retained verified archive: {path.name}", flush=True)
                    continue
                raw = run_learned(fixture, seed)
                payload = json.dumps(raw, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
                with path.open("xb") as stream:
                    stream.write(gzip.compress(payload, mtime=0))
        result = build_report()
    json_text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    md_text = render(result)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"All prespecified checks passed: {result['passed']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
