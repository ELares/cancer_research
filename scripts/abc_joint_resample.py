#!/usr/bin/env python3
"""Resample/move pilots and bounded frozen-proposal joint calibration.

The scientific target is unchanged. See docs/JOINT_RESAMPLE_LOCAL_PLAN.md for the
prespecified design and its limits. Historical experiments are not overwritten.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import abc_joint_importance as base
from bounded_proposal import BoundedGaussianMixture
from importance_sampling import importance_diagnostics, normalized_weights, weighted_quantiles
from resample_move_local import PILOT_PLAN, fit_proposal, train_pilot, validate_pilot_archive

OUT = ROOT / "analysis" / "calibration"
OUT_MD = OUT / "joint-resample-sampling.md"
OUT_JSON = OUT / "joint-resample-sampling.json"
SEEDS = (2026091911, 2026091912, 2026091913)
PLAN = {"pilot": dict(PILOT_PLAN), "production_attempts": 8192}
GATES = dict(base.GATES)
NAMES = list(base.NAMES)
PROBS = base.PROBS
abc = base.abc
_curve = base._curve
validate_target = base.validate_target
run_adequacy = base.run_adequacy
cdf_mcse = base.cdf_mcse

def summarize_run(archive):
    target, production = archive["target"], archive["production"]
    if archive["schema_version"] != 2:
        raise ValueError("unsupported importance archive schema")
    validate_target(target)
    validate_pilot(archive)
    priors = target["priors"]
    names = [p[0] for p in priors]
    lows = np.array([p[1] for p in priors], dtype=float)
    widths = np.array([p[2] - p[1] for p in priors], dtype=float)
    if names != NAMES or not np.all(np.isfinite(lows)) or not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("invalid archived prior definition")
    records = production["records"]
    points = np.asarray(production["unit_points"], dtype=float)
    log_q = np.asarray(production["log_q"], dtype=float)
    proposal = BoundedGaussianMixture.from_dict(archive["proposal"])
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
        if record["status"] != "outside_prior":
            _curve(record["erastin"], ne, "archived erastin curve")
            erastin_error = abc.ck.rmse(record["erastin"], target["empirical_erastin"])
            if record["status"] == "early_rejected":
                if (record["rsl3"] is not None or erastin_error <= target["epsilon"] or
                        not math.isclose(erastin_error, record["distance_lower_bound"], rel_tol=0, abs_tol=1e-14)):
                    raise ValueError("early-rejected curve does not reproduce its strict bound")
            else:
                _curve(record["rsl3"], nr, "archived RSL3 curve")
                distance = abc.ck.rmse(record["rsl3"], target["empirical_rsl3"]) + erastin_error
                if (record["accepted"] != (distance <= target["epsilon"]) or
                        not math.isclose(distance, record["distance"], rel_tol=0, abs_tol=1e-14)):
                    raise ValueError("complete curves do not reproduce the archived decision")
    diagnostics = importance_diagnostics(log_q, accepted)
    total_calls = (archive["pilot"]["simulator_dose_calls"] +
                   sum(r["simulator_dose_calls"] for r in records) +
                   target["reference_simulator_dose_calls"])
    result = {"seed": archive["seed"], "importance": diagnostics,
              "prior_widths": dict(zip(names, widths.tolist())),
              "pilot_attempts": archive["pilot"]["attempts"],
              "production_outside_prior": sum(r["status"] == "outside_prior" for r in records),
              "production_early_rejections": sum(r["status"] == "early_rejected" for r in records),
              "simulator_dose_calls_including_pilot_reference": total_calls,
              "ess_per_1000_simulator_dose_calls": 1000 * diagnostics["ess"] / total_calls,
              "wall_seconds": archive["wall_seconds"],
              "all_islands_reached_epsilon": archive["pilot"]["all_islands_reached_epsilon"],
              "pilot": [{key: island[key] for key in ("island", "first_epsilon_level", "history")}
                        for island in archive["pilot"]["islands"]]}
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
    checks["all_pilots_reached_epsilon"] = all(r["all_islands_reached_epsilon"] for r in runs)
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
    return directory / f"joint-resample-run-{seed}.json.gz"

def write_archive(path, archive):
    # Stable gzip header; archive timings themselves are explicitly runtime data.
    payload = json.dumps(archive, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    path.write_bytes(gzip.compress(payload, mtime=0))

def read_archive(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def validate_pilot(archive):
    """Reconcile the fixed pilot budget and final cloud before using its cost."""
    pilot, plan = archive['pilot'], archive['plan']['pilot']
    validate_pilot_archive(pilot)
    target = archive['target']
    if pilot['plan'] != plan or len(pilot['islands']) != plan['islands']:
        raise ValueError('pilot plan or island count differs')
    n, d = plan['particles_per_island'], len(NAMES)
    nr, ne = len(target['rsl3_doses_um']), len(target['erastin_doses_um'])
    pilot_stream = np.random.SeedSequence(archive['seed']).spawn(2)[0]
    island_streams = pilot_stream.spawn(plan['islands'])
    total_attempts = total_calls = 0
    all_final, all_scores = [], []
    for number, island in enumerate(pilot['islands']):
        if island['island'] != number or island['spawn_key'] != [0, number]:
            raise ValueError('pilot island does not match its independent stream')
        for key in ('initial_points', 'final_points'):
            points = np.asarray(island[key], dtype=float)
            if points.shape != (n, d) or not np.all(np.isfinite(points)) or np.any((points < 0) | (points > 1)):
                raise ValueError('invalid pilot cloud')
        initial = np.random.default_rng(island_streams[number]).uniform(size=(n, d))
        if not np.array_equal(initial, island['initial_points']):
            raise ValueError('initial pilot cloud does not reproduce its independent stream')
        for key in ('initial_distances', 'final_distances'):
            scores = np.asarray(island[key], dtype=float)
            if scores.shape != (n,) or not np.all(np.isfinite(scores)) or np.any(scores < 0):
                raise ValueError('invalid pilot scores')
        ancestors = np.asarray(island['final_ancestors'])
        if ancestors.shape != (n,) or ancestors.dtype.kind not in 'iu' or np.any((ancestors < 0) | (ancestors >= n)):
            raise ValueError('invalid pilot ancestry')
        if island['initial_simulator_dose_calls'] != n * (nr + ne):
            raise ValueError('initial pilot cost is inconsistent')
        history = island['history']
        if len(history) != plan['levels']:
            raise ValueError('pilot level budget differs from plan')
        previous = math.inf
        first_reached = None
        for level, entry in enumerate(history, start=1):
            threshold = entry['threshold']
            attempts = n * plan['moves_per_level']
            if (entry['level'] != level or not math.isfinite(threshold) or
                    not target['epsilon'] <= threshold <= previous):
                raise ValueError('pilot threshold sequence is inconsistent')
            for key in ('attempts', 'evaluated', 'mh_rejected_without_evaluation', 'early_rejected',
                        'accepted_moves', 'simulator_dose_calls', 'unique_ancestors', 'unique_particles',
                        'retained_before_resampling'):
                if type(entry[key]) is not int or entry[key] < 0:
                    raise ValueError('pilot counters must be nonnegative integers')
            if (entry['attempts'] != attempts or
                    entry['evaluated'] + entry['mh_rejected_without_evaluation'] != attempts or
                    not entry['accepted_moves'] <= entry['evaluated'] - entry['early_rejected'] or
                    not 1 <= entry['unique_ancestors'] <= n or
                    not 1 <= entry['unique_particles'] <= n or
                    not 1 <= entry['retained_before_resampling'] <= n or
                    type(entry['resampled']) is not bool or
                    entry['resampled'] != (entry['retained_before_resampling'] < n)):
                raise ValueError('pilot counters do not reconcile')
            expected_calls = ne * entry['early_rejected'] + (nr + ne) * (entry['evaluated'] - entry['early_rejected'])
            if entry['simulator_dose_calls'] != expected_calls:
                raise ValueError('pilot model calls do not reconcile')
            if threshold == target['epsilon'] and first_reached is None:
                first_reached = level
            if not 0 <= entry['minimum_distance'] <= entry['median_distance'] <= threshold:
                raise ValueError('pilot score summaries escape their threshold')
            total_attempts += attempts
            total_calls += expected_calls
            previous = threshold
        final_scores = np.asarray(island['final_distances'])
        last = history[-1]
        if (np.any(final_scores > last['threshold']) or
                float(np.min(final_scores)) != last['minimum_distance'] or
                float(np.median(final_scores)) != last['median_distance'] or
                len(np.unique(ancestors)) != last['unique_ancestors'] or
                len(np.unique(island['final_points'], axis=0)) != last['unique_particles'] or
                island['first_epsilon_level'] != first_reached):
            raise ValueError('final pilot state does not match its history')
        total_attempts += n
        total_calls += island['initial_simulator_dose_calls']
        all_final.extend(island['final_points'])
        all_scores.extend(island['final_distances'])
    if (pilot['final_points'] != all_final or pilot['final_distances'] != all_scores or
            pilot['attempts'] != total_attempts or pilot['simulator_dose_calls'] != total_calls or
            pilot['all_islands_reached_epsilon'] != all(i['first_epsilon_level'] is not None for i in pilot['islands'])):
        raise ValueError('combined pilot state or cost differs from its islands')
    expected_proposal = fit_proposal(all_final, plan).to_dict()
    actual = archive['proposal']
    if actual['dimension'] != d or actual['uniform_weight'] != plan['uniform_weight']:
        raise ValueError('production proposal does not follow pilot plan')
    for key in ('means', 'scales'):
        if np.shape(actual[key]) != np.shape(expected_proposal[key]) or not np.allclose(actual[key], expected_proposal[key], rtol=1e-12, atol=1e-12):
            raise ValueError('production proposal does not reproduce from final pilot state')


def source_provenance():
    hashes = base.source_provenance()
    for name in ('abc_joint_resample.py', 'resample_move.py', 'bounded_proposal.py',
                 'proposal_synthetic_validation.py', 'resample_move_local.py',
                 'proposal_synthetic_validation_v2.py', 'synthetic_proposal_study.py'):
        relative = 'scripts/' + name
        hashes[relative] = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
    return hashes


def verify_synthetic_study():
    """Require the committed known-target study before any biological run."""
    import proposal_synthetic_validation_v2 as synthetic
    path = OUT / 'proposal-synthetic-validation-v2.json'
    raw = path.read_bytes()
    study = synthetic.assemble(json.loads(raw))
    if not study['passed'] or study['pilot_plan'] != PLAN['pilot']:
        raise ValueError('known-target checks must pass under the current pilot plan before biological production')
    return {'artifact': path.name, 'sha256': hashlib.sha256(raw).hexdigest(),
            'source_hashes': study['source_hashes'], 'passed': study['passed']}


def run_experiment(seed, target, plan=None, evaluator=base.evaluate_unit):
    plan = dict(PLAN if plan is None else plan)
    pilot_stream, production_stream = np.random.SeedSequence(seed).spawn(2)
    started = time.perf_counter()

    def evaluate(point, threshold=None):
        temporary = target if threshold is None else dict(target, epsilon=threshold)
        return evaluator(point, temporary, early_reject=threshold is not None)

    def progress(island, entry):
        print(f"seed {seed}: island {island+1}/{plan['pilot']['islands']}, "
              f"level {entry['level']}/{plan['pilot']['levels']}, "
              f"threshold={entry['threshold']:.6g}, moves={entry['accepted_moves']}/{entry['attempts']}, "
              f"ancestors={entry['unique_ancestors']}", flush=True)

    pilot = train_pilot(evaluate, len(NAMES), target['epsilon'], pilot_stream, plan['pilot'], progress)
    proposal = fit_proposal(pilot['final_points'], plan['pilot'])
    production_rng = np.random.default_rng(production_stream)
    points, ids = proposal.sample(production_rng, plan['production_attempts'])
    log_q = proposal.log_density(points)
    records = []
    for i, point in enumerate(points):
        value = evaluator(point, target, early_reject=True)
        # Retain every evaluated production curve, including rejects, so an
        # offline audit can reconstruct every distance or strict early bound.
        records.append(value)
        if (i + 1) % 2048 == 0:
            print(f"seed {seed}: production {i+1}/{len(points)}", flush=True)
    return {'schema_version': 2, 'seed': int(seed), 'plan': plan, 'gates': GATES,
            'target': target, 'pilot': pilot, 'proposal': proposal.to_dict(),
            'rng': {'bit_generator': type(production_rng.bit_generator).__name__,
                    'pilot_spawn_key': list(pilot_stream.spawn_key),
                    'production_spawn_key': list(production_stream.spawn_key)},
            'production': {'unit_points': points.tolist(), 'log_q': log_q.tolist(),
                           'component_ids': ids.tolist(), 'records': records},
            'wall_seconds': time.perf_counter() - started}


def replay_production(archive):
    pilot_stream, stream = np.random.SeedSequence(archive['seed']).spawn(2)
    if (archive['rng']['pilot_spawn_key'] != list(pilot_stream.spawn_key) or
            archive['rng']['production_spawn_key'] != list(stream.spawn_key)):
        raise ValueError('archived random streams do not match seed')
    rng = np.random.default_rng(stream)
    if archive['rng']['bit_generator'] != type(rng.bit_generator).__name__:
        raise ValueError('archived bit generator differs')
    proposal = BoundedGaussianMixture.from_dict(archive['proposal'])
    points, ids = proposal.sample(rng, archive['plan']['production_attempts'])
    production = archive['production']
    if (not np.array_equal(ids, production['component_ids']) or
            np.shape(production['unit_points']) != points.shape or
            not np.allclose(points, production['unit_points'], rtol=1e-12, atol=1e-12)):
        raise ValueError('production draws do not reproduce their independent stream')


def build_report(directory=OUT):
    archives = [read_archive(archive_path(directory, seed)) for seed in SEEDS]
    first = archives[0]
    prerequisite = verify_synthetic_study()
    if first['source_hashes'] != source_provenance():
        raise ValueError('archived numerical sources differ from the current report implementation')
    for seed, item in zip(SEEDS, archives):
        if item['seed'] != seed or item['plan'] != PLAN or item['gates'] != GATES:
            raise ValueError('archive does not use the prespecified seed, plan and screens')
        if item['target'] != first['target'] or item['source_hashes'] != first['source_hashes']:
            raise ValueError('independent runs must share the same target and code')
        if item['synthetic_validation'] != prerequisite:
            raise ValueError('archive does not match the required known-target validation')
        for key in ('extension_binaries', 'python', 'numpy', 'scipy'):
            if item['runtime'][key] != first['runtime'][key]:
                raise ValueError('independent runs must share the build and numerical runtime')
        replay_production(item)
    runs = [summarize_run(item) for item in archives]
    previous_path = OUT / 'joint-importance-sampling.json'
    previous_bytes = previous_path.read_bytes()
    previous = json.loads(previous_bytes)
    if previous['target'] != first['target']:
        raise ValueError('historical experiment no longer matches the fixed target')
    uniform_bytes = (OUT / 'joint-posterior.json').read_bytes()
    if hashlib.sha256(uniform_bytes).hexdigest() != previous['baseline']['sha256']:
        raise ValueError('historical uniform baseline changed')
    return {'schema_version': 2, 'plan': PLAN, 'gates': GATES, 'target': first['target'],
            'source_hashes': first['source_hashes'], 'runtime': first['runtime'],
            'synthetic_validation': prerequisite,
            'archive_sha256': {archive_path(directory, seed).name: hashlib.sha256(archive_path(directory, seed).read_bytes()).hexdigest() for seed in SEEDS},
            'runs': runs, 'stability': compare_runs(runs),
            'baselines': {'uniform': previous['baseline'],
                          'first_importance': {'artifact': previous_path.name,
                              'sha256': hashlib.sha256(previous_bytes).hexdigest(),
                              'runs': [{key: run[key] for key in ('seed', 'importance', 'ess_per_1000_simulator_dose_calls')} for run in previous['runs']]}}}


def render(report):
    passed = report['stability']['passed']
    lines = ['# Joint calibration: resample/move pilots and bounded importance sampling', '',
             'Generated by `scripts/abc_joint_resample.py` from the archived attempts.', '',
             '**Computational screens passed.**' if passed else '**Computational screens failed; pooled inference is withheld.**', '',
             'These are operational sampling checks under a deterministic tolerance target. Passing',
             'them does not establish experimental validity, precise tails, or discovery of every region.', '',
             '## Fixed target and method', '',
             f"Reference distance {report['target']['reference_distance']:.10f}; factor {report['target']['tolerance_factor']}; final epsilon {report['target']['epsilon']:.10f}.",
             'Priors, dose support, empirical curves, simulation seed and cell count are unchanged.',
             'Each run trains two independent islands using local reflected and global independence',
             'moves at intermediate thresholds, with gap-limited neighborhoods. It then freezes',
             'a mixture of bounded Gaussian product kernels and 20% uniform prior. A separate stream',
             'draws 8,192 independent production attempts; weights use the entire normalized mixture.',
             'Pilots and resampling multiplicities are excluded from inference. See the',
             '[prespecified plan](../../docs/JOINT_RESAMPLE_LOCAL_PLAN.md) and',
             '[known-target validation](proposal-synthetic-validation-v2.md).', '',
             '## Per-run diagnostics', '',
             '| Seed | Accepted | ESS | Maximum weight | Relative normalizer MCSE | Dose calls incl. pilot/reference | Adequacy |',
             '|---|---:|---:|---:|---:|---:|---|']
    for run in sorted(report['runs'], key=lambda x: x['seed']):
        d = run['importance']
        mx = 'unavailable' if d['max_normalized_weight'] is None else f"{d['max_normalized_weight']:.4f}"
        se = 'unavailable' if d['normalizer_relative_mcse'] is None else f"{d['normalizer_relative_mcse']:.4f}"
        lines.append(f"| {run['seed']} | {d['n_accepted']} | {d['ess']:.1f} | {mx} | {se} | {run['simulator_dose_calls_including_pilot_reference']} | {'pass' if run['adequacy']['passed'] else 'fail'} |")
    lines += ['', 'ESS is a weight diagnostic, not experimental replication or proof of independent-region coverage.',
              'All production points are inside the prior. Exact MH pre-rejections and erastin bounds save',
              'pilot model calls; the archived counters include every attempted move.', '',
              '## Pilot progression', '',
              '| Seed | Island | First level at epsilon | Final unique particles | Surviving initial ancestors | Final threshold |',
              '|---|---:|---:|---:|---:|---:|']
    for run in sorted(report['runs'], key=lambda x: x['seed']):
        for island in sorted(run['pilot'], key=lambda x: x['island']):
            last = island['history'][-1]
            reached = island['first_epsilon_level'] if island['first_epsilon_level'] is not None else 'not reached'
            lines.append(f"| {run['seed']} | {island['island']+1} | {reached} | {last['unique_particles']} | {last['unique_ancestors']} | {last['threshold']:.6f} |")
    lines += ['', 'Ancestry counts expose resampling loss; unique particles are not independent posterior draws.', '',
              '## Independent-run screens', '']
    for key, value in sorted(report['stability']['checks'].items()):
        lines.append(f"- `{key}`: **{'pass' if value else 'fail'}**.")
    lines += ['', '| Parameter | Lower-quantile span / prior width | Median span / prior width | Upper-quantile span / prior width |',
              '|---|---:|---:|---:|']
    for name in NAMES:
        value = report['stability']['parameter_spans_prior_fraction'].get(name)
        if value:
            lines.append(f"| `{name}` | {value['q2_5']:.4f} | {value['median']:.4f} | {value['q97_5']:.4f} |")
    lines += ['', 'These are run-to-run diagnostic spreads, not confidence intervals. JSON retains per-run',
              'weighted quantiles and CDF Monte Carlo errors. No pooled posterior is published.', '',
              '## Held-out compound and normalizing mass', '',
              '| Seed | Weighted median ML210 RMSE | Fitted targets within parameter-draw band | Normalizer | MCSE |',
              '|---|---:|---|---:|---:|']
    for run in sorted(report['runs'], key=lambda x: x['seed']):
        h, d = run['heldout'], run['importance']
        if h:
            lines.append(f"| {run['seed']} | {h['rmse_quantiles'][1]:.4f} | {h['coverage']} | {d['normalizer_estimate']:.6g} | {d['normalizer_mcse']:.6g} |")
    lines += ['', 'ML210 never enters pilot adaptation or acceptance. It is held out by compound within the',
              'same screen, with overlapping cell lines. Bands describe parameter draws, not experimental',
              'outcomes. The normalizer estimates prior mass meeting the fixed tolerance; its MCSE is',
              'conditional on the frozen proposal. Neither it nor CDF MCSE measures unobserved modes.', '',
              '## Historical comparison and limits', '',
              f"The uniform baseline accepted {report['baselines']['uniform']['accepted']} of {report['baselines']['uniform']['draws']:,} draws.",
              'The first frozen-Gaussian experiment accepted 1, 4 and 4 draws and failed its screens.',
              'The new ESS per 1,000 simulator dose calls includes pilot and reference costs:', '']
    for run in sorted(report['runs'], key=lambda x: x['seed']):
        lines.append(f"- Seed {run['seed']}: {run['ess_per_1000_simulator_dose_calls']:.4f}.")
    lines += ['', 'This is observed sampling efficiency across different designs and budgets, not a controlled',
              'wall-clock speedup. Historical underpowered draws are not an oracle or extra production samples.', '',
              'The known-target negative control shows that usual ESS screens can pass while a separated',
              'region is missed. The uniform component preserves support but may visit a rare region too',
              'infrequently to reveal it. Three agreeing runs do not rule out shared blind spots.',
              'Independent assay validation and transfer to spatial headline claims remain unresolved.',
              'A failed screen is recorded without changing the criterion or selecting a favorable seed.', '']
    return '\n'.join(lines)


def write_report(directory=OUT):
    report = build_report(directory)
    OUT_JSON = directory / 'joint-resample-sampling.json'
    OUT_MD = directory / 'joint-resample-sampling.md'
    json_text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n'
    md_text = render(report)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, choices=SEEDS)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    parser.add_argument('--render-only', action='store_true')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.render_only:
        report = write_report(args.output_dir)
    else:
        prerequisite = verify_synthetic_study()
        hashes = source_provenance()
        runtime = base.runtime_provenance()
        for seed in ([args.seed] if args.seed is not None else SEEDS):
            target = base.prepare_target()
            archive = run_experiment(seed, target)
            archive.update(source_hashes=hashes, runtime=runtime, synthetic_validation=prerequisite)
            if source_provenance() != hashes:
                raise SystemExit('source changed during sampling; refusing mixed provenance')
            summary = summarize_run(archive)
            replay_production(archive)
            write_archive(archive_path(args.output_dir, seed), archive)
            print(f"seed {seed}: accepted={summary['importance']['n_accepted']}, ESS={summary['importance']['ess']:.1f}, adequate={summary['adequacy']['passed']}", flush=True)
        if not all(archive_path(args.output_dir, seed).exists() for seed in SEEDS):
            print('Individual run archived; summary awaits all three planned runs.')
            return 0
        report = write_report(args.output_dir)
    print(f"Independent-run screens passed: {report['stability']['passed']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
