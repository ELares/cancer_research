#!/usr/bin/env python3
"""Pilot/frozen-proposal importance sampling of the fixed joint ABC target.

See docs/JOINT_SAMPLING_PLAN.md for the evaluation plan. This does not replace
joint-posterior.json, change its threshold, or validate experimental outcomes.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import abc_joint_posterior as abc
from importance_sampling import (GaussianMixture, fit_component,
                                 importance_diagnostics, normalized_weights,
                                 weighted_quantiles)

OUT = ROOT / "analysis" / "calibration"
OUT_MD = OUT / "joint-importance-sampling.md"
OUT_JSON = OUT / "joint-importance-sampling.json"
SEEDS = (2026091901, 2026091902, 2026091903)
PLAN = {"pilot_rounds": 4, "pilot_attempts": 1024, "production_attempts": 8192,
        "elite_fraction": 0.2, "min_elite": 32, "uniform_weight": 0.2,
        "covariance_scale": 2.0, "covariance_floor": 1e-6}
GATES = {"minimum_accepted": 20, "minimum_ess": 200.0,
         "maximum_normalized_weight": 0.02, "maximum_normalizer_relative_mcse": 0.1,
         "maximum_median_span_prior_fraction": 0.05,
         "maximum_tail_span_prior_fraction": 0.1,
         "maximum_heldout_rmse_median_span": 0.02,
         "normalizer_combined_mcse_multiplier": 4.0}
PROBS = (0.025, 0.5, 0.975)
NAMES = [p[0] for p in abc.PRIORS]
LOWS = np.array([p[1] for p in abc.PRIORS])
WIDTHS = np.array([p[2] - p[1] for p in abc.PRIORS])


def _curve(values, length, label):
    array = np.asarray(values, dtype=float)
    if array.shape != (length,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{label} must have {length} finite values")
    return array


def validate_target(target):
    for key, empirical in (("rsl3_doses_um", "empirical_rsl3"),
                           ("erastin_doses_um", "empirical_erastin")):
        doses = np.asarray(target[key], dtype=float)
        if doses.ndim != 1 or not len(doses) or not np.all(np.isfinite(doses)) or np.any(doses <= 0):
            raise ValueError("target doses must be nonempty, positive and finite")
        _curve(target[empirical], len(doses), empirical)
    _curve(target["empirical_heldout"], len(target["rsl3_doses_um"]), "heldout target")
    if not math.isfinite(target["epsilon"]) or target["epsilon"] <= 0:
        raise ValueError("target epsilon must be positive and finite")


def prepare_target():
    curves, source = abc.ck.load_target_data()
    rd, ed = list(abc.ck.DOSE_GRID_UM), list(abc.ce.DOSE_GRID_UM)
    er, sr = abc.ck.empirical_target(curves[abc.ck.FIT_COMPOUND], rd)
    ee, se = abc.ck.empirical_target(curves[abc.ce.COMPOUND], ed)
    eh, sh = abc.ck.empirical_target(curves[abc.ck.HELDOUT_GPX4I], rd)
    reference = (abc.ck.rmse(abc.model_rsl3(rd, abc.REFERENCE_VECTOR), er)
                 + abc.ck.rmse(abc.model_erastin(ed, abc.REFERENCE_VECTOR), ee))
    return {"source": source, "support": {"ML162": sr, "ERASTIN": se, "ML210": sh},
            "priors": [list(p) for p in abc.PRIORS], "rsl3_doses_um": rd,
            "erastin_doses_um": ed, "empirical_rsl3": er,
            "empirical_erastin": ee, "empirical_heldout": eh,
            "reference_vector": abc.REFERENCE_VECTOR, "reference_distance": reference,
            "tolerance_factor": abc.TOLERANCE_FACTOR,
            "epsilon": reference * abc.TOLERANCE_FACTOR,
            "simulation_n": abc.SIM_N, "simulation_seed": abc.SIM_SEED,
            "reference_simulator_dose_calls": len(rd) + len(ed)}


def evaluate_unit(point, target, early_reject=False):
    """Score one proposal; out-of-prior attempts and exact rejections are explicit."""
    validate_target(target)
    u = np.asarray(point, dtype=float)
    if u.shape != LOWS.shape or not np.all(np.isfinite(u)):
        raise ValueError("unit point must contain seven finite coordinates")
    empty = {"distance": None, "distance_lower_bound": None, "rsl3": None,
             "erastin": None, "simulator_dose_calls": 0, "accepted": False}
    if np.any(u < 0) or np.any(u > 1):
        return {**empty, "status": "outside_prior"}
    params = dict(zip(NAMES, LOWS + WIDTHS * u))
    erastin = abc.model_erastin(target["erastin_doses_um"], params)
    _curve(erastin, len(target["erastin_doses_um"]), "simulated erastin curve")
    ee = abc.ck.rmse(erastin, target["empirical_erastin"])
    if not math.isfinite(ee):
        raise ValueError("simulator returned a nonfinite erastin distance")
    if early_reject and ee > target["epsilon"]:
        return {**empty, "status": "early_rejected", "erastin": erastin,
                "distance_lower_bound": ee,
                "simulator_dose_calls": len(target["erastin_doses_um"])}
    rsl3 = abc.model_rsl3(target["rsl3_doses_um"], params)
    _curve(rsl3, len(target["rsl3_doses_um"]), "simulated RSL3 curve")
    distance = abc.ck.rmse(rsl3, target["empirical_rsl3"]) + ee
    if not math.isfinite(distance):
        raise ValueError("simulator returned a nonfinite joint distance")
    return {"distance": distance, "distance_lower_bound": distance,
            "rsl3": rsl3, "erastin": erastin, "status": "complete",
            "simulator_dose_calls": len(rsl3) + len(erastin),
            "accepted": distance <= target["epsilon"]}


def source_provenance():
    paths = [ROOT / "scripts" / name for name in
             ("abc_joint_importance.py", "importance_sampling.py", "abc_joint_posterior.py",
              "calibrate_kill_switch.py", "calibrate_erastin.py", "ctrp_dose_support.py")]
    paths += sorted((ROOT / "simulations" / "ferroptosis-core" / "src").rglob("*.rs"))
    paths += [ROOT / "simulations" / "ferroptosis-python" / "src" / "lib.rs",
              ROOT / "simulations" / "Cargo.lock",
              ROOT / "simulations" / "rust-toolchain.toml"]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def runtime_provenance():
    import scipy
    fc = abc.ck._fc()
    binaries = sorted(Path(fc.__file__).parent.glob("*.so"))
    return {"python": platform.python_version(), "platform": platform.platform(),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "logical_cpus": os.cpu_count(), "rayon_num_threads_env": os.getenv("RAYON_NUM_THREADS"),
            "extension_binaries": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in binaries}}


def run_experiment(seed, target, plan=None, evaluator=evaluate_unit):
    plan = dict(PLAN if plan is None else plan)
    pilot_stream, production_stream = np.random.SeedSequence(seed).spawn(2)
    pilot_rng = np.random.default_rng(pilot_stream)
    production_rng = np.random.default_rng(production_stream)
    proposal = GaussianMixture(len(NAMES))
    means, covariances, history = [], [], []
    pilot_calls = 0
    started = time.perf_counter()
    for step in range(plan["pilot_rounds"]):
        points, _ = proposal.sample(pilot_rng, plan["pilot_attempts"])
        log_q = proposal.log_density(points)
        scores = np.full(len(points), np.inf)
        calls = 0
        for i, point in enumerate(points):
            record = evaluator(point, target, early_reject=False)
            calls += record["simulator_dose_calls"]
            if record["distance"] is not None:
                scores[i] = record["distance"]
        usable = np.flatnonzero(np.isfinite(scores))
        n_elite = min(len(usable), max(plan["min_elite"],
                                      math.ceil(plan["elite_fraction"] * len(usable))))
        elite = usable[np.argsort(scores[usable], kind="stable")[:n_elite]]
        fitted = len(elite) >= len(NAMES) + 1
        if fitted:
            mean, covariance = fit_component(
                points[elite], -log_q[elite],
                covariance_scale=plan["covariance_scale"],
                covariance_floor=plan["covariance_floor"])
            means.append(mean)
            covariances.append(covariance)
            proposal = GaussianMixture(len(NAMES), uniform_weight=plan["uniform_weight"],
                                       means=means, covariances=covariances)
        pilot_calls += calls
        history.append({"round": step + 1, "attempts": len(points), "in_prior": len(usable),
                        "simulator_dose_calls": calls, "n_elite": len(elite),
                        "selection_distance": float(scores[elite[-1]]) if len(elite) else None,
                        "minimum_distance": float(scores[elite[0]]) if len(elite) else None,
                        "component_fitted": fitted})
        print(f"seed {seed}: pilot {step+1}/{plan['pilot_rounds']}, "
              f"in-prior={len(usable)}, elite={len(elite)}, "
              f"cutoff={history[-1]['selection_distance']}", flush=True)
    # This proposal is now frozen. No production observation can alter it.
    points, component_ids = proposal.sample(production_rng, plan["production_attempts"])
    log_q = proposal.log_density(points)
    records = []
    for i, point in enumerate(points):
        record = evaluator(point, target, early_reject=True)
        # Only accepted curves are needed for weighted diagnostics; retain every
        # attempt's point, density, decision, and distance/bound in the archive.
        if not record["accepted"]:
            record["rsl3"] = record["erastin"] = None
        records.append(record)
        if (i + 1) % 2048 == 0:
            print(f"seed {seed}: production {i+1}/{len(points)}", flush=True)
    return {"schema_version": 1, "seed": int(seed), "plan": plan, "gates": GATES,
            "target": target, "pilot": history, "proposal": proposal.to_dict(),
            "pilot_simulator_dose_calls": pilot_calls,
            "rng": {"bit_generator": type(pilot_rng.bit_generator).__name__,
                    "pilot_spawn_key": list(pilot_stream.spawn_key),
                    "production_spawn_key": list(production_stream.spawn_key)},
            "production": {"unit_points": points.tolist(), "log_q": log_q.tolist(),
                           "component_ids": component_ids.tolist(), "records": records},
            "wall_seconds": time.perf_counter() - started}


def run_adequacy(summary):
    d = summary["importance"]
    checks = {"accepted_count": d["n_accepted"] >= GATES["minimum_accepted"],
              "ess": d["ess"] >= GATES["minimum_ess"],
              "maximum_weight": d["max_normalized_weight"] is not None and
                  d["max_normalized_weight"] <= GATES["maximum_normalized_weight"],
              "normalizer_mcse": d["normalizer_relative_mcse"] is not None and
                  d["normalizer_relative_mcse"] <= GATES["maximum_normalizer_relative_mcse"]}
    return {"passed": all(checks.values()), "checks": checks}


def cdf_mcse(values, weights, cut, attempts):
    indicator = np.asarray(values) <= cut
    estimate = float(np.sum(weights * indicator))
    mcse = math.sqrt(attempts / (attempts - 1) *
                     float(np.sum(weights**2 * (indicator - estimate)**2)))
    return {"cdf": estimate, "mcse": mcse}


def summarize_run(archive):
    target, production = archive["target"], archive["production"]
    if archive["schema_version"] != 1:
        raise ValueError("unsupported importance archive schema")
    validate_target(target)
    priors = target["priors"]
    names = [p[0] for p in priors]
    lows = np.array([p[1] for p in priors], dtype=float)
    widths = np.array([p[2] - p[1] for p in priors], dtype=float)
    if names != NAMES or not np.all(np.isfinite(lows)) or not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("invalid archived prior definition")
    records = production["records"]
    points = np.asarray(production["unit_points"], dtype=float)
    log_q = np.asarray(production["log_q"], dtype=float)
    proposal = GaussianMixture.from_dict(archive["proposal"])
    if (points.shape != (archive["plan"]["production_attempts"], len(names)) or
            not np.all(np.isfinite(points)) or len(records) != len(points) or
            log_q.shape != (len(points),) or not np.all(np.isfinite(log_q))):
        raise ValueError("archive production count differs from its plan")
    if not np.allclose(log_q, proposal.log_density(points), rtol=1e-12, atol=1e-12):
        raise ValueError("archived proposal densities do not reproduce")
    if any(type(r["accepted"]) is not bool for r in records):
        raise ValueError("archive acceptance flags must be booleans")
    accepted = np.array([r["accepted"] for r in records], dtype=bool)
    nr, ne = len(target["rsl3_doses_um"]), len(target["erastin_doses_um"])
    for point, record in zip(points, records):
        inside = bool(np.all((point >= 0) & (point <= 1)))
        if record["status"] == "outside_prior":
            valid = (not inside and not record["accepted"] and record["simulator_dose_calls"] == 0 and
                     record["distance"] is None and record["distance_lower_bound"] is None)
        elif record["status"] == "early_rejected":
            bound = record["distance_lower_bound"]
            valid = (inside and not record["accepted"] and record["distance"] is None and
                     bound is not None and math.isfinite(bound) and bound > target["epsilon"] and
                     record["simulator_dose_calls"] == ne)
        else:
            valid = (record["status"] == "complete" and inside and
                     record["distance"] is not None and math.isfinite(record["distance"]) and
                     record["distance"] >= 0 and record["simulator_dose_calls"] == nr + ne and
                     record["distance_lower_bound"] == record["distance"] and
                     record["accepted"] == (record["distance"] <= target["epsilon"]))
        if not valid:
            raise ValueError("archive decision is inconsistent with the fixed criterion")
    diagnostics = importance_diagnostics(log_q, accepted)
    total_calls = (archive["pilot_simulator_dose_calls"] +
                   sum(r["simulator_dose_calls"] for r in records) +
                   target["reference_simulator_dose_calls"])
    result = {"seed": archive["seed"], "importance": diagnostics,
              "prior_widths": dict(zip(names, widths.tolist())),
              "pilot_attempts": sum(h["attempts"] for h in archive["pilot"]),
              "production_outside_prior": sum(r["status"] == "outside_prior" for r in records),
              "production_early_rejections": sum(r["status"] == "early_rejected" for r in records),
              "simulator_dose_calls_including_pilot_reference": total_calls,
              "ess_per_1000_simulator_dose_calls": 1000 * diagnostics["ess"] / total_calls,
              "wall_seconds": archive["wall_seconds"], "pilot": archive["pilot"]}
    result["adequacy"] = run_adequacy(result)
    if not np.any(accepted):
        result["diagnostic_quantiles"] = None
        result["heldout"] = None
        return result
    weights = normalized_weights(-log_q[accepted])
    params = lows + points[accepted] * widths
    quantiles = {}
    for j, name in enumerate(names):
        q = weighted_quantiles(params[:, j], weights, PROBS)
        quantiles[name] = {"q2_5": float(q[0]), "median": float(q[1]), "q97_5": float(q[2]),
                           "cdf_mcse_at_quantiles": [cdf_mcse(params[:, j], weights, x, len(points)) for x in q]}
    selected = [r for r in records if r["accepted"]]
    for record in selected:
        _curve(record["rsl3"], nr, "archived accepted RSL3 curve")
        _curve(record["erastin"], ne, "archived accepted erastin curve")
    rsl3 = np.array([r["rsl3"] for r in selected])
    erastin = np.array([r["erastin"] for r in selected])
    # Recheck every accepted score from its archived curves, at full precision.
    distances = [abc.ck.rmse(r["rsl3"], target["empirical_rsl3"]) +
                 abc.ck.rmse(r["erastin"], target["empirical_erastin"]) for r in selected]
    if (any(d > target["epsilon"] for d in distances) or
            not np.allclose(distances, [r["distance"] for r in selected], rtol=0, atol=1e-14)):
        raise ValueError("accepted curves do not reproduce their distances")
    rmse = np.array([abc.ck.rmse(r, target["empirical_heldout"]) for r in rsl3])
    bands = np.array([weighted_quantiles(rsl3[:, j], weights, PROBS) for j in range(rsl3.shape[1])]).T
    inside = (np.asarray(target["empirical_heldout"]) >= bands[0]) & (np.asarray(target["empirical_heldout"]) <= bands[2])
    result["diagnostic_quantiles"] = quantiles
    result["heldout"] = {"rmse_quantiles": weighted_quantiles(rmse, weights, PROBS).tolist(),
                         "band_q2_5": bands[0].tolist(), "band_median": bands[1].tolist(),
                         "band_q97_5": bands[2].tolist(), "coverage": f"{int(inside.sum())}/{len(inside)}"}
    return result


def compare_runs(runs):
    unique_seeds = {r["seed"] for r in runs}
    checks = {"three_prespecified_independent_runs": len(runs) == len(SEEDS) and unique_seeds == set(SEEDS),
              "all_runs_adequate": all(r["adequacy"]["passed"] for r in runs)}
    prior_widths = runs[0]["prior_widths"] if runs else {}
    checks["same_prior_widths"] = bool(runs) and all(r["prior_widths"] == prior_widths for r in runs)
    spans = {}
    if runs and all(r["diagnostic_quantiles"] is not None and r["heldout"] is not None for r in runs):
        for name in NAMES:
            width = prior_widths[name]
            spans[name] = {key: float(np.ptp([r["diagnostic_quantiles"][name][key] for r in runs]) / width)
                           for key in ("q2_5", "median", "q97_5")}
        checks["median_stability"] = all(v["median"] <= GATES["maximum_median_span_prior_fraction"] for v in spans.values())
        checks["tail_stability"] = all(max(v["q2_5"], v["q97_5"]) <= GATES["maximum_tail_span_prior_fraction"] for v in spans.values())
        heldout_span = float(np.ptp([r["heldout"]["rmse_quantiles"][1] for r in runs]))
        checks["heldout_rmse_stability"] = heldout_span <= GATES["maximum_heldout_rmse_median_span"]
    else:
        heldout_span = None
        checks.update(median_stability=False, tail_stability=False, heldout_rmse_stability=False)
    comparisons = []
    for i, a in enumerate(runs):
        for b in runs[i+1:]:
            da, db = a["importance"], b["importance"]
            difference = abs(da["normalizer_estimate"] - db["normalizer_estimate"])
            sa, sb = da["normalizer_mcse"], db["normalizer_mcse"]
            combined = math.hypot(sa, sb) if sa is not None and sb is not None else None
            comparisons.append({"seeds": [a["seed"], b["seed"]], "difference": difference,
                                "combined_mcse": combined,
                                "passed": combined is not None and difference <= GATES["normalizer_combined_mcse_multiplier"] * combined})
    checks["normalizer_agreement"] = bool(comparisons) and all(c["passed"] for c in comparisons)
    return {"passed": all(checks.values()), "checks": checks, "parameter_spans_prior_fraction": spans,
            "heldout_rmse_median_span": heldout_span, "normalizer_comparisons": comparisons}


def archive_path(directory, seed):
    return directory / f"joint-importance-run-{seed}.json.gz"


def write_archive(path, archive):
    # Stable gzip header; archive timings themselves are explicitly runtime data.
    payload = json.dumps(archive, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    path.write_bytes(gzip.compress(payload, mtime=0))


def read_archive(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def build_report(directory=OUT):
    archives = [read_archive(archive_path(directory, seed)) for seed in SEEDS]
    first = archives[0]
    for item in archives:
        if item["plan"] != PLAN or item["gates"] != GATES:
            raise ValueError("archive does not use the prespecified plan and screens")
        if item["target"] != first["target"] or item["source_hashes"] != first["source_hashes"]:
            raise ValueError("independent runs must share the same target and code")
        for key in ("extension_binaries", "python", "numpy", "scipy"):
            if item["runtime"][key] != first["runtime"][key]:
                raise ValueError("independent runs must share the simulator build and numerical runtime")
    runs = [summarize_run(item) for item in archives]
    baseline_bytes = (OUT / "joint-posterior.json").read_bytes()
    baseline = json.loads(baseline_bytes)
    target = first["target"]
    if (baseline["target_source"] != target["source"] or
            baseline["target_support"] != target["support"] or
            baseline["curves"]["rsl3_doses_um"] != target["rsl3_doses_um"] or
            baseline["curves"]["erastin_doses_um"] != target["erastin_doses_um"] or
            baseline["tolerance_factor"] != target["tolerance_factor"] or
            abs(baseline["epsilon_joint_distance"] - target["epsilon"]) > 0.00005):
        raise ValueError("uniform baseline no longer matches the archived sampling target")
    dose_calls = baseline["n_draws"] * (len(baseline["curves"]["rsl3_doses_um"]) + len(baseline["curves"]["erastin_doses_um"]))
    return {"schema_version": 1, "plan": PLAN, "gates": GATES, "target": first["target"],
            "source_hashes": first["source_hashes"], "runtime": first["runtime"],
            "archive_sha256": {archive_path(directory, seed).name: hashlib.sha256(archive_path(directory, seed).read_bytes()).hexdigest() for seed in SEEDS},
            "runs": runs, "stability": compare_runs(runs),
            "baseline": {"artifact": "joint-posterior.json", "draws": baseline["n_draws"],
                         "sha256": hashlib.sha256(baseline_bytes).hexdigest(),
                         "accepted": baseline["n_accepted"], "underpowered": baseline["underpowered"],
                         "sampling_simulator_dose_calls_excluding_overhead": dose_calls,
                         "ess_per_1000_sampling_dose_calls": 1000 * baseline["n_accepted"] / dose_calls}}


def render(report):
    passed = report["stability"]["passed"]
    lines = ["# Joint calibration: frozen-proposal importance sampling", "",
             "Generated by `scripts/abc_joint_importance.py` from three archived runs.", "",
             "**Computational screens passed.**" if passed else "**Computational screens failed; pooled inference is withheld.**", "",
             "These screens assess finite-sample behavior under the fixed tolerance target.",
             "They do not establish biological validity, precise tails, or complete mode discovery.", "",
             "## Fixed target and method", "",
             f"Reference distance: {report['target']['reference_distance']:.10f}; unchanged factor "
             f"{report['target']['tolerance_factor']}; final epsilon: {report['target']['epsilon']:.10f}.",
             "Priors, supported cohorts, simulation seed and cell count match the corrected rejection run.",
             "Pilot draws train a proposal and are excluded from the estimates. Each production run uses",
             "8,192 independent attempts from its frozen mixture, including zero-weight outside-prior attempts.",
             "Weights use the entire proposal density; there is no resampling to inflate the sample count.",
             "The three experiments have independent pilot and production streams. See",
             "[the prespecified plan](../../docs/JOINT_SAMPLING_PLAN.md) for budgets and screens.", "",
             "## Per-run diagnostics", "",
             "| Seed | Accepted | ESS | Maximum weight | Normalizer | Relative MCSE | Dose calls incl. pilot/reference | Screens |",
             "|---|---:|---:|---:|---:|---:|---:|---|"]
    for run in sorted(report["runs"], key=lambda r: r["seed"]):
        d = run["importance"]
        mx = f"{d['max_normalized_weight']:.4f}" if d["max_normalized_weight"] is not None else "unavailable"
        se = f"{d['normalizer_relative_mcse']:.4f}" if d["normalizer_relative_mcse"] is not None else "unavailable"
        lines.append(f"| {run['seed']} | {d['n_accepted']} | {d['ess']:.1f} | {mx} | {d['normalizer_estimate']:.6g} | {se} | {run['simulator_dose_calls_including_pilot_reference']} | {'pass' if run['adequacy']['passed'] else 'fail'} |")
    lines += ["", "ESS is an importance-weight diagnostic, not a count of experimental replicates.",
              "The normalizer estimates the fraction of the original uniform prior meeting the tolerance.",
              "Its MCSE is Monte Carlo error conditional on the frozen proposal, not model or assay uncertainty.", "",
              "The 20% uniform component preserves support but does not guarantee discovery of missed modes.",
              "At the historical rejection rate, each production run expects only about 0.16 accepted",
              "uniform-component draws. Agreement across pilots can still miss a common blind spot.", "",
              "## Independent-run screens", ""]
    for name, value in sorted(report["stability"]["checks"].items()):
        lines.append(f"- `{name}`: **{'pass' if value else 'fail'}**.")
    lines += ["", "| Parameter | Lower-quantile span / prior width | Median span / prior width | Upper-quantile span / prior width |",
              "|---|---:|---:|---:|"]
    for name in NAMES:
        span = report["stability"]["parameter_spans_prior_fraction"].get(name)
        if span:
            lines.append(f"| `{name}` | {span['q2_5']:.4f} | {span['median']:.4f} | {span['q97_5']:.4f} |")
    lines += ["", "These spreads are diagnostic comparisons, not confidence intervals. The JSON retains",
              "per-run weighted quantiles and CDF Monte Carlo errors for audit. No pooled posterior is published.", "",
              "## Held-out compound and baseline", "",
              "| Seed | Weighted median ML210 RMSE | Fitted targets within parameter-draw band |",
              "|---|---:|---|"]
    for run in sorted(report["runs"], key=lambda r: r["seed"]):
        h = run["heldout"]
        if h:
            lines.append(f"| {run['seed']} | {h['rmse_quantiles'][1]:.4f} | {h['coverage']} |")
    b = report["baseline"]
    lines += ["", "These held-out diagnostics do not enter pilot selection or acceptance. ML210 is held out",
              "by compound within the same screen, with overlapping cell lines. Bands reflect parameter",
              "draws conditional on fixed simulation settings, not experimental-outcome uncertainty.", "",
              f"The original uniform run accepted {b['accepted']} of {b['draws']:,} draws and remains underpowered.",
              f"Its ESS per 1,000 sampling dose calls is {b['ess_per_1000_sampling_dose_calls']:.4f}, excluding overhead.",
              "The new per-run efficiencies include pilot and reference costs:"]
    for run in sorted(report["runs"], key=lambda r: r["seed"]):
        lines.append(f"- Seed {run['seed']}: {run['ess_per_1000_simulator_dose_calls']:.4f} ESS per 1,000 dose calls.")
    lines += ["", "This compares observed computational efficiency across different budgets, not a controlled",
              "wall-clock benchmark. The four-draw rejection result is not a posterior oracle.", "",
              "The historical rejection artifact and equal-weight information analysis remain separate.",
              "Do not apply their accepted-count contraction null to these weighted samples. Production",
              "simulation parameters and spatial headline claims are unchanged; independent assay validation",
              "remains pending. A failed screen is reported without changing the tolerance or selecting a",
              "favorable seed. Further runs require a separate plan.", ""]
    return "\n".join(lines)


def write_report(directory=OUT):
    report = build_report(directory)
    OUT_JSON = directory / "joint-importance-sampling.json"
    OUT_MD = directory / "joint-importance-sampling.md"
    json_text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    md_text = render(report)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, choices=SEEDS, help="run one planned experiment")
    parser.add_argument("--output-dir", type=Path, default=OUT)
    parser.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.render_only:
        report = write_report(args.output_dir)
    else:
        hashes = source_provenance()
        runtime = runtime_provenance()
        for seed in ([args.seed] if args.seed is not None else SEEDS):
            target = prepare_target()
            archive = run_experiment(seed, target)
            archive.update(source_hashes=hashes, runtime=runtime)
            if source_provenance() != hashes:
                raise SystemExit("source changed during sampling; refusing to publish mixed provenance")
            write_archive(archive_path(args.output_dir, seed), archive)
            d = summarize_run(archive)
            print(f"seed {seed}: accepted={d['importance']['n_accepted']}, ESS={d['importance']['ess']:.1f}, adequate={d['adequacy']['passed']}", flush=True)
        if not all(archive_path(args.output_dir, seed).exists() for seed in SEEDS):
            print("Individual run archived; summary awaits all three planned runs.")
            return 0
        report = write_report(args.output_dir)
    print(f"Independent-run screens passed: {report['stability']['passed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
