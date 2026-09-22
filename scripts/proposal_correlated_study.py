#!/usr/bin/env python3
"""Frozen comparison of final covariance proposals on shared, unchanged pilots."""
import argparse
import copy
import gzip
import json
import math
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import scipy

import correlated_challenge_geometry as geometry
import correlated_proposal as correlated
import proposal_coverage_challenges as historical
from proposal_coverage_challenges import _equal, _close, _hash
from coverage_pilot_audit import validate_pilot
from bounded_proposal import BoundedGaussianMixture
from importance_sampling import GaussianMixture, importance_diagnostics, normalized_weights
import resample_move_local as strategy

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / 'analysis/calibration/proposal-correlated-study.json'
OUT_MD = ROOT / 'analysis/calibration/proposal-correlated-study.md'
ARCHIVES = ROOT / 'analysis/calibration/correlated-study-20260921'
PROTOCOL = 'docs/CORRELATED_PROPOSAL_PLAN.md'
ATTEMPTS = 8192
EPSILON = 1.0
ARMS = ('bounded_diagonal', 'unbounded_diagonal', 'correlated')
LEARNED_SEEDS = {
    'rotated_box': (2026092801, 2026092802, 2026092803),
    'annular_cylinder': (2026092811, 2026092812, 2026092813),
    'unequal_balls': (2026092821, 2026092822, 2026092823),
    'shifted_rotated_box': (2026092831, 2026092832, 2026092833),
}
ORACLE_SEEDS = {name: tuple(seed + 100 for seed in seeds) for name, seeds in LEARNED_SEEDS.items()}
NEGATIVE_SEEDS = {name: tuple(seed + 200 for seed in seeds) for name, seeds in LEARNED_SEEDS.items()}
GATES = dict(historical.GATES)
PRIOR_PATH = 'analysis/calibration/proposal-coverage-challenges.json'
PRIOR_SHA256 = '6252872268d765da98485a4f36dac099ef1f4f2fd5883e846780bf334af04940'
SOURCE_PATHS = historical.SOURCE_PATHS + (
    'scripts/correlated_proposal.py', 'scripts/correlated_challenge_geometry.py',
    'scripts/proposal_correlated_study.py', PROTOCOL,
)


def specification():
    if _hash(ROOT / PRIOR_PATH) != PRIOR_SHA256:
        raise ValueError('the completed coverage study changed')
    prior = json.loads((ROOT / PRIOR_PATH).read_text())
    _equal(prior['specification'], historical.specification(), 'historical specification')
    return {
        'revision': 'final-covariance-20260921', 'dimension': geometry.DIMENSION,
        'epsilon': EPSILON, 'production_attempts_per_arm': ATTEMPTS,
        'pilot_plan': copy.deepcopy(strategy.PILOT_PLAN),
        'fit_plan': copy.deepcopy(correlated.FIT_PLAN), 'arms': list(ARMS),
        'learned_seeds': {name: list(seeds) for name, seeds in LEARNED_SEEDS.items()},
        'oracle_seeds': {name: list(seeds) for name, seeds in ORACLE_SEEDS.items()},
        'negative_seeds': {name: list(seeds) for name, seeds in NEGATIVE_SEEDS.items()},
        'gates': dict(GATES), 'source_hashes': {p: _hash(ROOT / p) for p in SOURCE_PATHS},
        'prior_study': {'path': PRIOR_PATH, 'sha256': PRIOR_SHA256},
    }


def committed_sources():
    """Refuse an experiment before every numerical source/protocol is committed."""
    for path in SOURCE_PATHS:
        committed = subprocess.run(['git', 'show', f'HEAD:{path}'], cwd=ROOT,
                                   capture_output=True, check=True).stdout
        if committed != (ROOT / path).read_bytes():
            raise ValueError(f'commit the protocol and numerical sources before evaluation: {path}')
    return subprocess.check_output(
        ['git', 'log', '-1', '--format=%H', '--', *SOURCE_PATHS], cwd=ROOT, text=True).strip()


def proposals(points):
    full = correlated.fit_proposal(points)
    diagonal = [np.diag(np.diag(cov)) for cov in full.covariances]
    return {
        'bounded_diagonal': strategy.fit_proposal(points, strategy.PILOT_PLAN),
        'unbounded_diagonal': GaussianMixture(full.dimension, full.uniform_weight,
                                             full.means, diagonal),
        'correlated': full,
    }


def summarize(points, log_q, fixture):
    points, log_q = np.asarray(points, dtype=float), np.asarray(log_q, dtype=float)
    if (points.shape != (ATTEMPTS, geometry.DIMENSION) or not np.isfinite(points).all()
            or log_q.shape != (ATTEMPTS,) or not np.isfinite(log_q).all()):
        raise ValueError('production requires the fixed number of finite points and densities')
    inside = ((points >= 0) & (points <= 1)).all(axis=1)
    regions = geometry.membership(points, fixture)
    hits = np.stack(list(regions.values()))
    if np.any(hits.sum(axis=0) > 1):
        raise ValueError('fixture regions overlap')
    accepted = hits.any(axis=0)
    if not np.array_equal(accepted, inside & (geometry.scores(points, fixture) <= EPSILON)):
        raise ValueError('analytic scores and region membership disagree')
    diagnostics = importance_diagnostics(log_q, accepted)
    gold = geometry.truth(fixture)
    counts = {name: int(mask.sum()) for name, mask in regions.items()}
    relative_error = diagnostics['normalizer_estimate'] / gold['mass'] - 1
    usual = {
        'accepted_count': diagnostics['n_accepted'] >= GATES['minimum_accepted'],
        'ess': diagnostics['ess'] >= GATES['minimum_ess'],
        'maximum_weight': diagnostics['max_normalized_weight'] is not None and
            diagnostics['max_normalized_weight'] <= GATES['maximum_normalized_weight'],
        'relative_mcse': diagnostics['normalizer_relative_mcse'] is not None and
            diagnostics['normalizer_relative_mcse'] <= GATES['maximum_relative_mcse'],
    }
    if np.any(accepted):
        weights = normalized_weights(-log_q[accepted])
        features = geometry.features(points[accepted], fixture)
        if set(features) != set(gold['moments']):
            raise ValueError('feature names differ from analytic truth')
        for values in features.values():
            if (np.shape(values) != weights.shape or not np.isfinite(values).all()
                    or np.any((values < -1e-12) | (values > 1 + 1e-12))):
                raise ValueError('features must be bounded on the target')
        moments = {name: float(weights @ value) for name, value in features.items()}
        masses = {name: float(weights[mask[accepted]].sum()) for name, mask in regions.items()}
        moment_errors = {name: abs(value - gold['moments'][name]) for name, value in moments.items()}
        region_errors = {name: abs(value - gold['region_masses'][name]) for name, value in masses.items()}
    else:
        moments = masses = moment_errors = region_errors = None
    truth_checks = {
        'normalizing_mass': abs(relative_error) <= GATES['maximum_relative_mass_error'],
        'every_region_observed': all(counts.values()),
        'region_masses': region_errors is not None and
            max(region_errors.values()) <= GATES['maximum_region_mass_error'],
        'moments': moment_errors is not None and
            max(moment_errors.values()) <= GATES['maximum_moment_error'],
    }
    return {'importance': diagnostics, 'outside_cube': int((~inside).sum()),
            'region_counts': counts, 'region_masses': masses, 'moments': moments,
            'absolute_moment_errors': moment_errors, 'absolute_region_mass_errors': region_errors,
            'relative_mass_error': relative_error, 'usual_checks': usual, 'truth_checks': truth_checks,
            'passed': all(usual.values()) and all(truth_checks.values())}


def archive_path(directory, fixture, seed):
    if type(seed) is not int or fixture not in LEARNED_SEEDS or seed not in LEARNED_SEEDS[fixture]:
        raise ValueError('unplanned fixture or seed')
    return directory / f'{fixture}-{seed}.json.gz'


def run_learned(fixture, seed):
    archive_path(ARCHIVES, fixture, seed)
    commit, spec = committed_sources(), specification()
    started = time.perf_counter()
    streams = np.random.SeedSequence(seed).spawn(4)
    pilot = strategy.train_pilot(
        lambda point, threshold=None: geometry.evaluate(point, fixture, threshold),
        geometry.DIMENSION, EPSILON, streams[0], plan=copy.deepcopy(strategy.PILOT_PLAN))
    validate_pilot(pilot, fixture, seed, geometry.evaluate, EPSILON)
    arms = {}
    for index, (name, proposal) in enumerate(proposals(pilot['final_points']).items(), 1):
        rng = np.random.default_rng(streams[index])
        points, ids = proposal.sample(rng, ATTEMPTS)
        scores = geometry.scores(points, fixture)
        inside = ((points >= 0) & (points <= 1)).all(axis=1)
        arms[name] = {'spawn_key': list(streams[index].spawn_key), 'proposal': proposal.to_dict(),
                      'production': {'points': points.tolist(), 'component_ids': ids.tolist(),
                                     'log_q': proposal.log_density(points).tolist(),
                                     'scores': scores.tolist(), 'accepted': (inside & (scores <= EPSILON)).tolist()}}
    result = {'schema_version': 1, 'fixture': fixture, 'seed': seed, 'specification': spec,
              'implementation_commit': commit,
              'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                          'scipy': scipy.__version__, 'bit_generator': type(rng.bit_generator).__name__},
              'pilot': pilot, 'arms': arms, 'wall_seconds': time.perf_counter() - started}
    _equal(specification(), spec, 'sources after run')
    for name, a in assess_archive(result, spec).items():
        print(f"{fixture} {seed} {name}: ESS={a['importance']['ess']:.1f}, "
              f"mass error={a['relative_mass_error']:+.4f}, pass={a['passed']}", flush=True)
    return result


def assess_archive(raw, spec=None):
    spec = specification() if spec is None else spec
    _equal(raw['schema_version'], 1, 'archive schema')
    fixture, seed = raw['fixture'], raw['seed']
    archive_path(ARCHIVES, fixture, seed)
    _equal(raw['specification'], spec, 'archive specification')
    commit = raw['implementation_commit']
    if type(commit) is not str or len(commit) != 40 or any(c not in '0123456789abcdef' for c in commit):
        raise ValueError('invalid implementation commit')
    runtime = raw['runtime']
    if (type(runtime) is not dict or set(runtime) != {'python', 'numpy', 'scipy', 'bit_generator'}
            or any(type(v) is not str or not v.strip() for v in runtime.values())):
        raise ValueError('invalid runtime labels')
    if (type(raw['wall_seconds']) not in (int, float) or
            not math.isfinite(raw['wall_seconds']) or raw['wall_seconds'] < 0):
        raise ValueError('invalid elapsed time')
    reached = validate_pilot(raw['pilot'], fixture, seed, geometry.evaluate, EPSILON)
    expected = proposals(raw['pilot']['final_points'])
    _equal(sorted(raw['arms']), sorted(ARMS), 'proposal arms')
    streams = np.random.SeedSequence(seed).spawn(4)
    results = {}
    for index, name in enumerate(ARMS, 1):
        arm = raw['arms'][name]
        _equal(arm['spawn_key'], list(streams[index].spawn_key), 'production stream')
        actual, reference = arm['proposal'], expected[name].to_dict()
        _equal(sorted(actual), sorted(reference), 'proposal fields')
        for key in ('dimension', 'uniform_weight'):
            _equal(actual[key], reference[key], f'proposal {key}')
        for key in reference.keys() - {'dimension', 'uniform_weight'}:
            _close(actual[key], reference[key], f'proposal {key}')
        cls = BoundedGaussianMixture if name == 'bounded_diagonal' else GaussianMixture
        proposal = cls.from_dict(actual)
        rng = np.random.default_rng(streams[index])
        _equal(runtime['bit_generator'], type(rng.bit_generator).__name__, 'bit generator')
        points, ids = proposal.sample(rng, ATTEMPTS)
        stored = arm['production']
        _close(stored['points'], points, 'production points')
        _equal(stored['component_ids'], ids.tolist(), 'component IDs')
        points = np.asarray(stored['points'])
        log_q = proposal.log_density(points)
        _close(stored['log_q'], log_q, 'full mixture density')
        scores = geometry.scores(points, fixture)
        _close(stored['scores'], scores, 'production scores')
        accepted = ((points >= 0) & (points <= 1)).all(axis=1) & (scores <= EPSILON)
        _equal(stored['accepted'], accepted.tolist(), 'acceptance decisions')
        a = summarize(points, log_q, fixture)
        a['pilot_reached_epsilon'] = reached
        a['passed'] = a['passed'] and reached
        results[name] = a
    return results


def run_control(fixture, seed, restricted):
    points, log_q = geometry.oracle_sample(np.random.default_rng(seed), ATTEMPTS, fixture, restricted)
    return {'fixture': fixture, 'seed': seed,
            'kind': 'support_hole_oracle' if restricted else 'full_target_oracle',
            'assessment': summarize(points, log_q, fixture)}


def archive_hashes(directory):
    expected = {archive_path(directory, name, seed) for name, seeds in LEARNED_SEEDS.items() for seed in seeds}
    if set(directory.glob('*.json.gz')) != expected:
        raise ValueError('archive set must contain every planned run and no extras')
    return {p.name: _hash(p) for p in sorted(expected)}


def build_report(directory=None):
    directory = ARCHIVES if directory is None else directory
    spec, hashes = specification(), archive_hashes(directory)
    runs, runtime, commit = [], None, None
    for fixture, seeds in LEARNED_SEEDS.items():
        for seed in seeds:
            path = archive_path(directory, fixture, seed)
            raw = json.loads(gzip.decompress(path.read_bytes()))
            _equal([raw['fixture'], raw['seed']], [fixture, seed], 'archive identity')
            if runtime is None:
                runtime, commit = raw['runtime'], raw['implementation_commit']
            _equal(raw['runtime'], runtime, 'run runtimes')
            _equal(raw['implementation_commit'], commit, 'implementation commit')
            assessments = assess_archive(raw, spec)
            for arm in ARMS:
                runs.append({'fixture': fixture, 'seed': seed, 'arm': arm, 'archive': path.name,
                             'assessment': assessments[arm], 'shared_pilot_attempts': raw['pilot']['attempts']})
    oracles = [run_control(f, seed, False) for f, seeds in ORACLE_SEEDS.items() for seed in seeds]
    negatives = [run_control(f, seed, True) for f, seeds in NEGATIVE_SEEDS.items() for seed in seeds]
    checks = {
        'all_correlated_runs_pass': all(r['assessment']['passed'] for r in runs if r['arm'] == 'correlated'),
        'all_full_target_oracles_pass': all(r['assessment']['passed'] for r in oracles),
        'all_support_holes_pass_usual': all(all(r['assessment']['usual_checks'].values()) for r in negatives),
        'all_support_holes_fail_truth': all(not all(r['assessment']['truth_checks'].values()) for r in negatives),
    }
    _equal(specification(), spec, 'sources after reconstruction')
    _equal(archive_hashes(directory), hashes, 'archives after reconstruction')
    return {'schema_version': 1, 'specification': spec, 'runtime': runtime, 'implementation_commit': commit,
            'archive_sha256': hashes, 'truth': {f: geometry.truth(f) for f in LEARNED_SEEDS},
            'positive_runs': runs, 'oracle_controls': oracles, 'negative_controls': negatives,
            'checks': checks, 'passed': all(checks.values())}


def assemble(stored):
    _equal(stored['specification'], specification(), 'report specification')
    _equal(stored['archive_sha256'], archive_hashes(ARCHIVES), 'report archive hashes')
    return build_report()


def render(result):
    lines = ['# Final-proposal covariance comparison', '',
             '**All prespecified checks passed.**' if result['passed'] else
             '**Prespecified checks failed. Every planned outcome is retained.**', '',
             'Only the final independent production proposal changes. Each trio shares an unchanged',
             '29,184-attempt pilot; each arm then makes 8,192 independent production attempts.',
             'Out-of-cube Gaussian draws retain their full mixture densities, have zero target',
             'weight, and remain in the denominator. No boundary clipping or retry is used.', '',
             'The old three geometries are regression cases already seen during earlier work.',
             'Only shifted_rotated_box is a new prospective geometry, selected with knowledge of',
             'earlier failures. None of these results is independent biological validation.', '',
             'The unbounded diagonal and correlated arms share centers, multiplicities and marginal',
             'variances; they differ only in off-diagonal covariance. Comparison with the historical',
             'bounded arm also changes boundary handling and kernel marginal scales.', '',
             '## All proposal outcomes', '',
             '| Geometry | Seed | Arm | Outside | Accepted | ESS | Max weight | Mass error | Result |',
             '|---|---:|---|---:|---:|---:|---:|---:|---|']
    failures = []
    for r in sorted(result['positive_runs'], key=lambda r: (r['fixture'], r['seed'], ARMS.index(r['arm']))):
        a, d = r['assessment'], r['assessment']['importance']
        weight = 'unavailable' if d['max_normalized_weight'] is None else f"{d['max_normalized_weight']:.5f}"
        lines.append(f"| {r['fixture']} | {r['seed']} | {r['arm']} | {a['outside_cube']} | "
                     f"{d['n_accepted']} | {d['ess']:.1f} | {weight} | {a['relative_mass_error']:+.4f} | "
                     f"{'pass' if a['passed'] else 'fail'} |")
        failed = [k for group in ('usual_checks', 'truth_checks') for k, v in sorted(a[group].items()) if not v]
        if not a['pilot_reached_epsilon']:
            failed.append('pilot_reached_epsilon')
        if failed:
            failures.append(f"- {r['fixture']}, {r['seed']}, {r['arm']}: {', '.join(failed)}.")
    lines += ['', 'Failed checks:', '', *(failures or ['None.']), '', '## Analytic controls', '',
              '| Geometry | Seed | Kind | Usual screens | Truth screens |', '|---|---:|---|---|---|']
    for r in sorted(result['oracle_controls'] + result['negative_controls'], key=lambda r: (r['fixture'], r['kind'], r['seed'])):
        a = r['assessment']
        lines.append(f"| {r['fixture']} | {r['seed']} | {r['kind']} | "
                     f"{'pass' if all(a['usual_checks'].values()) else 'fail'} | "
                     f"{'pass' if all(a['truth_checks'].values()) else 'fail'} |")
    lines += ['', '## Interpretation and replay', '']
    lines += [f"- {key}: {'pass' if value else 'fail'}." for key, value in sorted(result['checks'].items())]
    lines += ['', 'The all-run rule is a fixed benchmark criterion, not a calibrated confidence statement.',
              'Arm comparisons are descriptive across three pilots per geometry; no superiority test or',
              'pooled posterior is inferred. Normalizer estimates include every attempted draw;',
              'self-normalized region masses and moments retain finite-sample bias.',
              'Finite moment and CDF checks cannot establish a complete distribution or precise tails.', '',
              'Archives retain shared pilot endpoints/history and every arm’s proposal, points, component',
              'IDs, densities, scores and decisions. Replay verifies source/archive hashes, reconstructs',
              'fits and production randomness, and audits pilot endpoints and counters. Intermediate',
              'pilot paths are not archived or independently replayed. Runtime labels describe archived',
              'learned runs; analytic controls are regenerated in the reader’s environment.',
              'Offline replay verifies source bytes against their recorded hashes, not the Git',
              'history behind the informational commit label. No Git database is needed for replay.', '',
              f"Recorded protocol/source commit: `{result['implementation_commit']}`.",
              'See [the protocol](../../docs/CORRELATED_PROPOSAL_PLAN.md) for all fixed inputs and criteria.', '',
              'Rebuild without training: `python scripts/proposal_correlated_study.py --render-only`.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--render-only', action='store_true')
    mode.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if args.render_only:
        result = assemble(json.loads(OUT_JSON.read_text()))
    else:
        committed_sources()
        ARCHIVES.mkdir(parents=True, exist_ok=True)
        for fixture, seeds in LEARNED_SEEDS.items():
            for seed in seeds:
                path = archive_path(ARCHIVES, fixture, seed)
                if path.exists():
                    if not args.resume:
                        raise FileExistsError(f'{path.name} exists; --resume validates and retains it')
                    raw = json.loads(gzip.decompress(path.read_bytes()))
                    _equal([raw['fixture'], raw['seed']], [fixture, seed], 'existing archive identity')
                    assess_archive(raw)
                    continue
                raw = run_learned(fixture, seed)
                payload = json.dumps(raw, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
                with path.open('xb') as stream:
                    stream.write(gzip.compress(payload, mtime=0))
        result = build_report()
    json_text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n'
    md_text = render(result)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"All prespecified checks passed: {result['passed']}", flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
