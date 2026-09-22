#!/usr/bin/env python3
"""Post hoc feature diagnostics computed only from the completed covariance archives.

No sampler, pilot, oracle, or random-number stream runs here. The completed
study is pinned by hash; its criteria, estimates, and outcomes are not revised.
"""
import argparse
import copy
import gzip
import json
import math
from pathlib import Path

import numpy as np

import correlated_challenge_geometry as geometry
import proposal_correlated_study as study
from bounded_proposal import BoundedGaussianMixture
from importance_sampling import GaussianMixture
from proposal_coverage_challenges import _close, _equal, _hash


ROOT = Path(__file__).resolve().parents[1]
STUDY_PATH = 'analysis/calibration/proposal-correlated-study.json'
STUDY_JSON = ROOT / STUDY_PATH
STUDY_SHA256 = 'df4bbab24034f9bea804ce7425a5c7a8cf6849ab1b60f498d33cf04e80a5d153'
ARCHIVES = ROOT / 'analysis/calibration/correlated-study-20260921'
SCRIPT_PATH = 'scripts/proposal_correlated_diagnostics.py'
OUT_JSON = ROOT / 'analysis/calibration/proposal-correlated-diagnostics.json'
OUT_MD = ROOT / 'analysis/calibration/proposal-correlated-diagnostics.md'
SELECTED_FEATURES = {
    'annular_cylinder': {
        'feature': 'angular_sin_2',
        'conditional_region_truth': {
            f'sector_{j}': 0.5 + sign / math.pi
            for j, sign in enumerate((1, 1, -1, -1, 1, 1, -1, -1))},
    },
    'shifted_rotated_box': {
        'feature': 'latent_cdf_6_50',
        'conditional_region_truth': {f'octant_{j}': 0.5 for j in range(8)},
    },
}


def _array(value, name, shape=None):
    array = np.asarray(value)
    if (array.dtype.kind not in 'fiu' or not np.isfinite(array).all()
            or (shape is not None and array.shape != shape)):
        raise ValueError(f'{name} must have finite real values and the expected shape')
    return array.astype(float)


def _weights(raw_weights):
    raw = _array(raw_weights, 'raw weights')
    if raw.ndim != 1 or len(raw) < 2 or np.any(raw < 0) or not np.any(raw > 0):
        raise ValueError('raw weights require at least two attempts and positive total weight')
    # Scaling first also makes these diagnostics invariant to the units of the
    # unnormalized weights without overflowing their sum.
    scaled = raw / raw.max()
    return raw, scaled / scaled.sum()


def weight_concentration(raw_weights):
    """Describe all attempted weights, with rejected attempts remaining zero."""
    raw, weights = _weights(raw_weights)
    ordered = np.sort(weights)[::-1]
    return {
        'n_attempts': len(weights), 'n_positive': int(np.count_nonzero(raw)),
        'ess': float(1 / (weights @ weights)),
        'max_normalized_weight': float(ordered[0]),
        **{f'top_{k}_weight': float(ordered[:k].sum()) for k in (1, 5, 20)},
    }


def feature_diagnostics(values, raw_weights, truth):
    """Self-normalized estimate, conditional delta-method MCSE, and sensitivity.

    Both vectors retain the original production order and all n attempts. The
    estimated variance is n/(n-1) * sum(w_i**2 * (f_i - estimate)**2), where
    w_i are normalized weights. This is conditional on the frozen proposal;
    it neither measures pilot uncertainty nor yields a calibrated post-selection
    test. A deletion changes only this descriptive sensitivity calculation.
    """
    raw, weights = _weights(raw_weights)
    values = _array(values, 'feature values', weights.shape)
    if type(truth) not in (int, float) or not math.isfinite(truth):
        raise ValueError('feature truth must be a finite real number')
    positive_values = values[raw > 0]
    constant = np.all(positive_values == positive_values[0])
    estimate = float(positive_values[0] if constant else weights @ values)
    residual = values - estimate
    mcse = (float(np.sqrt(len(weights) / (len(weights) - 1)
                          * np.sum((weights * residual) ** 2)))
            if len(positive_values) >= 2 else None)
    deletion = None
    if weights.max() < 1:
        changes = -weights * residual / (1 - weights)
        # np.argmax returns the first occurrence: an exact tie selects the
        # lowest original attempt index, including a zero-weight attempt.
        index = int(np.argmax(np.abs(changes)))
        deletion = {'attempt_index': index, 'signed_change': float(changes[index]),
                    'absolute_change': float(abs(changes[index]))}
    error = estimate - truth
    return {'estimate': estimate, 'signed_error': error, 'conditional_mcse': mcse,
            'error_over_mcse': float(error / mcse) if mcse is not None and mcse > 0 else None,
            'single_deletion': deletion}


def decompose_feature(values, raw_weights, regions, region_truth, conditional_truth):
    """Split signed error into between-region mass and within-region terms.

    mu_hat - mu = sum_r (p_hat_r-p_r)*mu_r
                  + sum_r sum_i_in_r w_i*(f_i-mu_r).
    The direct second sum remains defined when a region has no observations.
    """
    raw, weights = _weights(raw_weights)
    values = _array(values, 'feature values', weights.shape)
    if not regions or set(regions) != set(region_truth) or set(regions) != set(conditional_truth):
        raise ValueError('region masks and both truth mappings must have identical names')
    probabilities = _array(list(region_truth.values()), 'region probabilities')
    if np.any(probabilities < 0) or not np.isclose(probabilities.sum(), 1, rtol=0, atol=1e-14):
        raise ValueError('region probabilities must be nonnegative and sum to one')
    _array(list(conditional_truth.values()), 'conditional feature truth')
    coverage = np.zeros(len(weights), dtype=int)
    rows = {}
    for name, region in regions.items():
        mask = np.asarray(region)
        if mask.dtype.kind != 'b' or mask.shape != weights.shape:
            raise ValueError('region masks must be boolean arrays for every attempt')
        coverage += mask
        mass = float(weights[mask].sum())
        numerator = float(weights[mask] @ values[mask])
        rows[name] = {
            'n_accepted': int(np.count_nonzero(mask & (raw > 0))),
            'estimated_mass': mass,
            'conditional_estimate': numerator / mass if mass > 0 else None,
            'between_contribution': (mass - region_truth[name]) * conditional_truth[name],
            'within_contribution': float(weights[mask] @ (values[mask] - conditional_truth[name])),
        }
    if np.any(coverage > 1) or np.any(coverage[raw > 0] != 1):
        raise ValueError('regions must partition all positive-weight attempts without overlap')
    truth = sum(region_truth[name] * conditional_truth[name] for name in regions)
    error = float(weights @ values - truth)
    # Since both region distributions sum to one, subtracting the global
    # truth leaves the same between-region term. This form is exactly zero
    # when all conditional truths coincide, without cancellation roundoff.
    between = sum((rows[name]['estimated_mass'] - region_truth[name])
                  * (conditional_truth[name] - truth) for name in regions)
    within = sum(row['within_contribution'] for row in rows.values())
    _close(between + within, error, 'feature error decomposition')
    return {'total_signed_error': error, 'between_regions': between,
            'within_regions': within, 'regions': rows}


def _load_inputs():
    if _hash(STUDY_JSON) != STUDY_SHA256:
        raise ValueError('the completed covariance study hash changed')
    completed = json.loads(STUDY_JSON.read_text())
    _equal(completed['schema_version'], 1, 'completed study schema')
    _equal(completed['specification'], study.specification(), 'frozen covariance specification')
    _equal(completed['archive_sha256'], study.archive_hashes(ARCHIVES), 'completed archive hashes')
    _equal(completed['truth'], {f: geometry.truth(f) for f in study.LEARNED_SEEDS}, 'analytic truth')
    provenance = {
        'completed_study': {'path': STUDY_PATH, 'sha256': STUDY_SHA256},
        'archive_sha256': copy.deepcopy(completed['archive_sha256']),
        'source_hashes': {**completed['specification']['source_hashes'], SCRIPT_PATH: _hash(ROOT / SCRIPT_PATH)},
        'implementation_commit': completed['implementation_commit'],
        'archived_runtime': copy.deepcopy(completed['runtime']),
    }
    return completed, provenance


def _endpoint_summary(points, fixture):
    features = geometry.features(points, fixture)
    regions = geometry.membership(points, fixture)
    return {
        'n_endpoints': len(points), 'n_unique_endpoints': len(np.unique(points, axis=0)),
        'feature_means': {name: float(values.mean()) for name, values in features.items()},
        'region_counts': {name: int(mask.sum()) for name, mask in regions.items()},
        'region_fractions': {name: float(mask.mean()) for name, mask in regions.items()},
    }


def _pilot_endpoints(pilot, fixture):
    """Check archived endpoints and accounting without replaying any randomness."""
    plan = study.strategy.PILOT_PLAN
    _equal(pilot['plan'], plan, 'pilot plan')
    _equal(len(pilot['islands']), plan['islands'], 'pilot island count')
    n = plan['particles_per_island']
    clouds, scores, summaries = [], [], []
    for index, island in enumerate(pilot['islands']):
        _equal(island['island'], index, 'island index')
        _equal(island['spawn_key'], [0, index], 'archived island stream identity')
        points = _array(island['final_points'], 'pilot endpoints', (n, geometry.DIMENSION))
        stored_scores = _array(island['final_distances'], 'pilot endpoint scores', (n,))
        if np.any((points < 0) | (points > 1)) or np.any(stored_scores > study.EPSILON):
            raise ValueError('completed pilot endpoints escaped their final target')
        _close(stored_scores, geometry.scores(points, fixture), 'pilot endpoint scores')
        level = island['first_epsilon_level']
        if type(level) is not int or not 1 <= level <= plan['levels']:
            raise ValueError('completed pilot island did not reach the final threshold')
        _equal(len(island['history']), plan['levels'], 'pilot level count')
        _equal(island['history'][-1]['threshold'], study.EPSILON, 'final pilot threshold')
        clouds.append(points)
        scores.append(stored_scores)
        summaries.append({'island': index, 'summary': _endpoint_summary(points, fixture)})
    pooled = np.concatenate(clouds)
    _equal(pilot['final_points'], pooled.tolist(), 'combined pilot endpoints')
    _equal(pilot['final_distances'], np.concatenate(scores).tolist(), 'combined pilot scores')
    _equal(pilot['all_islands_reached_epsilon'], True, 'completed pilot reach')
    _equal(pilot['attempts'], plan['islands'] * n * (1 + plan['levels'] * plan['moves_per_level']),
           'pilot attempt accounting')
    _equal(pilot['simulator_dose_calls'], 0, 'analytic pilot cost')
    return pooled, {'pooled': _endpoint_summary(pooled, fixture), 'islands': summaries}


def _assessment_equal(actual, expected, label='original assessment'):
    """Only recomputed float values tolerate rounding; types and decisions do not."""
    if type(expected) is dict:
        if type(actual) is not dict or set(actual) != set(expected):
            raise ValueError(f'{label} fields differ')
        for key in expected:
            _assessment_equal(actual[key], expected[key], f'{label}.{key}')
    elif type(expected) is float:
        if type(actual) is not float:
            raise ValueError(f'{label} numeric type differs')
        _close(actual, expected, label)
    else:
        _equal(actual, expected, label)


def _archive_runs(raw, completed):
    fixture, seed = raw['fixture'], raw['seed']
    archive = study.archive_path(ARCHIVES, fixture, seed).name
    _equal(raw['schema_version'], 1, 'archive schema')
    _equal(raw['specification'], completed['specification'], 'archive specification')
    _equal(raw['implementation_commit'], completed['implementation_commit'], 'archive implementation')
    _equal(raw['runtime'], completed['runtime'], 'archived runtime')
    endpoints, pilot = _pilot_endpoints(raw['pilot'], fixture)
    expected_proposals = study.proposals(endpoints)
    _equal(sorted(raw['arms']), sorted(study.ARMS), 'archive arms')
    originals = {run['arm']: run['assessment'] for run in completed['positive_runs']
                 if (run['fixture'], run['seed']) == (fixture, seed)}
    _equal(sorted(originals), sorted(study.ARMS), 'completed outcome arms')
    runs = []
    for index, name in enumerate(study.ARMS, 1):
        arm = raw['arms'][name]
        _equal(arm['spawn_key'], [index], 'archived production stream identity')
        actual, expected = arm['proposal'], expected_proposals[name].to_dict()
        _equal(sorted(actual), sorted(expected), 'proposal fields')
        for key in ('dimension', 'uniform_weight'):
            _equal(actual[key], expected[key], f'proposal {key}')
        for key in expected.keys() - {'dimension', 'uniform_weight'}:
            _close(actual[key], expected[key], f'proposal {key}')
        cls = BoundedGaussianMixture if name == 'bounded_diagonal' else GaussianMixture
        proposal = cls.from_dict(actual)
        production = arm['production']
        points = _array(production['points'], 'production points', (study.ATTEMPTS, geometry.DIMENSION))
        log_q = _array(production['log_q'], 'archived log density', (study.ATTEMPTS,))
        _close(log_q, proposal.log_density(points), 'full mixture density')
        scores = geometry.scores(points, fixture)
        _close(production['scores'], scores, 'production scores')
        inside = ((points >= 0) & (points <= 1)).all(axis=1)
        accepted = inside & (scores <= study.EPSILON)
        _equal(production['accepted'], accepted.tolist(), 'archived acceptance decisions')
        ids = production['component_ids']
        if (type(ids) is not list or len(ids) != study.ATTEMPTS
                or any(type(k) is not int or not 0 <= k <= len(proposal.means) for k in ids)):
            raise ValueError('invalid archived component IDs')
        if np.any((np.asarray(ids) == 0) & ~inside):
            raise ValueError('uniform-component draw escaped the cube')
        # Recompute against the generating density stored with these exact draws.
        # No alternative proposal density substitutes for their actual q.
        assessment = study.summarize(points, log_q, fixture)
        assessment['pilot_reached_epsilon'] = True
        _assessment_equal(assessment, originals[name])
        raw_weights = np.zeros(study.ATTEMPTS)
        raw_weights[accepted] = np.exp(-log_q[accepted])
        if np.any(raw_weights > 5 + 1e-12):
            raise ValueError('defensive-mixture target weights exceed their bound')
        features = {}
        full_values = {}
        for feature, values in geometry.features(points[accepted], fixture).items():
            all_values = np.zeros(study.ATTEMPTS)
            all_values[accepted] = values
            full_values[feature] = all_values
            features[feature] = feature_diagnostics(
                all_values, raw_weights, completed['truth'][fixture]['moments'][feature])
        decomposition = None
        if fixture in SELECTED_FEATURES:
            selected = SELECTED_FEATURES[fixture]
            decomposition = decompose_feature(
                full_values[selected['feature']], raw_weights, geometry.membership(points, fixture),
                completed['truth'][fixture]['region_masses'], selected['conditional_region_truth'])
        failed = [key for group in ('usual_checks', 'truth_checks')
                  for key, passed in sorted(originals[name][group].items()) if not passed]
        runs.append({'fixture': fixture, 'seed': seed, 'arm': name, 'archive': archive,
                     'original_passed': originals[name]['passed'], 'original_failed_checks': failed,
                     'diagnostics': {'weights': weight_concentration(raw_weights),
                                     'features': features, 'decomposition': decomposition}})
    return runs, {'fixture': fixture, 'seed': seed, 'archive': archive, 'diagnostics': pilot}


def _rebuild(completed, provenance):
    runs, pilots = [], []
    for fixture, seeds in study.LEARNED_SEEDS.items():
        for seed in seeds:
            path = study.archive_path(ARCHIVES, fixture, seed)
            raw = json.loads(gzip.decompress(path.read_bytes()))
            _equal([raw['fixture'], raw['seed']], [fixture, seed], 'archive identity')
            trio, pilot = _archive_runs(raw, completed)
            runs.extend(trio)
            pilots.append(pilot)
    # A concurrent edit must not produce a report mixing different inputs.
    _, current = _load_inputs()
    _equal(current, provenance, 'inputs after diagnostic reconstruction')
    return {
        'schema_version': 1, 'provenance': provenance,
        'truth': copy.deepcopy(completed['truth']),
        'selected_features': copy.deepcopy(SELECTED_FEATURES),
        'completed_benchmark': {
            'passed': completed['passed'], 'checks': copy.deepcopy(completed['checks']),
            'correlated_runs_passed': sum(r['original_passed'] for r in runs if r['arm'] == 'correlated'),
            'correlated_run_count': sum(r['arm'] == 'correlated' for r in runs),
        },
        'runs': runs, 'pilots': pilots,
    }


def build_report():
    """Validate the pinned completed study and rebuild solely from its archives."""
    return _rebuild(*_load_inputs())


def assemble(stored=None):
    """Validate cached provenance, ignore cached diagnostics, and reread archives."""
    completed, provenance = _load_inputs()
    if stored is not None:
        if type(stored) is not dict:
            raise ValueError('stored diagnostic report must be an object')
        _equal(stored.get('schema_version'), 1, 'diagnostic report schema')
        _equal(stored.get('provenance'), provenance, 'diagnostic report provenance')
    return _rebuild(completed, provenance)


def render(result):
    benchmark = result['completed_benchmark']
    lines = [
        '# Archived covariance proposal feature diagnostics', '',
        f"The completed benchmark remains **failed: {benchmark['correlated_runs_passed']}/"
        f"{benchmark['correlated_run_count']} correlated runs passed**. Its estimates, thresholds and outcomes are unchanged.", '',
        'This is a post hoc descriptive analysis of all 12 shared pilots and all 36 production arms.',
        'The two highlighted features were selected after inspecting the completed failures.',
        'Every declared feature is retained in the numerical report. No new pilot, production draw,',
        'oracle draw, or random-number stream is run, and no biological sampler is changed.', '',
        '## Highlighted feature errors across every paired pilot', '',
        '| Fixture | Seed | Arm | Estimate | Signed error | Conditional MCSE | Largest deletion change | Attempt index |',
        '|---|---:|---|---:|---:|---:|---:|---:|',
    ]
    highlighted = sorted(
        (r for r in result['runs'] if r['fixture'] in result['selected_features']),
        key=lambda r: (r['fixture'], r['seed'], study.ARMS.index(r['arm'])))
    for run in highlighted:
        feature = result['selected_features'][run['fixture']]['feature']
        row = run['diagnostics']['features'][feature]
        deletion = row['single_deletion']
        change = 'unavailable' if deletion is None else f"{deletion['signed_change']:+.6f}"
        index = 'unavailable' if deletion is None else str(deletion['attempt_index'])
        mcse = 'unavailable' if row['conditional_mcse'] is None else f"{row['conditional_mcse']:.6f}"
        lines.append(f"| {run['fixture']} | {run['seed']} | {run['arm']} | {row['estimate']:.6f} | "
                     f"{row['signed_error']:+.6f} | {mcse} | {change} | {index} |")
    lines += [
        '', 'The annular feature is `(1 + sin(2 theta))/2`; the shifted-box feature is `P(T6 <= 0.5)`.',
        'Both full-target truths are 0.5. Attempt indices are zero-based original production indices.',
        'Deletion changes are the signed change in the self-normalized estimate after removing one',
        'attempt, selecting the largest absolute change and the lowest index for an exact tie.',
        'They are sensitivity summaries only; no point is removed from an original result.', '',
        '## Exact decomposition of the selected errors', '',
        '| Fixture | Seed | Arm | Between-region contribution | Within-region contribution | Total error |',
        '|---|---:|---|---:|---:|---:|',
    ]
    for run in highlighted:
        row = run['diagnostics']['decomposition']
        lines.append(f"| {run['fixture']} | {run['seed']} | {run['arm']} | "
                     f"{row['between_regions']:+.6f} | {row['within_regions']:+.6f} | {row['total_signed_error']:+.6f} |")
    lines += [
        '', 'For each region r, the error is `(p_hat_r - p_r) * mu_r + sum_i_in_r w_i * (f_i - mu_r)`,',
        'summed across regions. The eight annular sector conditional means are `0.5 +/- 1/pi`,',
        'with signs `++--++--`. The shifted feature has conditional mean 0.5 in every octant.',
        'Empty regions have an unavailable conditional estimate and a zero direct within-region term.',
        'Detailed per-region counts, masses and contributions remain in the JSON.', '',
        'The annular failing run has a substantial combined sector-mass contribution even though',
        'every individual region passed its original error limit. The shifted failure lies within',
        'the declared octants; their identities involve different latent coordinates from T6.', '',
        '## Descriptive pilot endpoint summaries', '',
        '| Fixture | Seed | Pooled selected-feature mean | Island 0 mean | Island 1 mean | Unique endpoints / total |',
        '|---|---:|---:|---:|---:|---:|',
    ]
    for pilot in sorted(result['pilots'], key=lambda p: (p['fixture'], p['seed'])):
        if pilot['fixture'] not in result['selected_features']:
            continue
        feature = result['selected_features'][pilot['fixture']]['feature']
        pooled = pilot['diagnostics']['pooled']
        islands = {row['island']: row['summary'] for row in pilot['diagnostics']['islands']}
        _equal(len(pilot['diagnostics']['islands']), 2, 'rendered pilot island count')
        _equal(sorted(islands), [0, 1], 'rendered pilot island identities')
        lines.append(f"| {pilot['fixture']} | {pilot['seed']} | {pooled['feature_means'][feature]:.6f} | "
                     f"{islands[0]['feature_means'][feature]:.6f} | "
                     f"{islands[1]['feature_means'][feature]:.6f} | "
                     f"{pooled['n_unique_endpoints']} / {pooled['n_endpoints']} |")
    lines += [
        '', 'Endpoint means and region fractions retain particle multiplicities. Unique endpoint counts',
        'describe repetition, not independent information. The JSON retains these summaries for every',
        'fixture, island and declared feature. Pilot endpoints are dependent training particles;',
        'no inferential ESS, MCSE, or posterior estimate is assigned to them.', '',
        '## Numerical definitions and limits', '',
        'Production weights use each archived draw’s actual full generating density, after independent',
        'density validation. Rejected and outside-cube attempts retain zero weight and stay in n.',
        'For normalized weights w, the feature variance estimate is',
        '`n/(n-1) * sum(w_i^2 * (f_i - estimate)^2)`. Its square root is the delta-method MCSE',
        'conditional on the frozen proposal. The JSON also records signed error divided by this MCSE',
        '(unavailable when MCSE is zero or unavailable), plus ESS, maximum weight and top-1/5/20 weight shares.',
        'MCSE and the largest deletion are unavailable when only one attempt has positive weight.',
        'These are descriptive diagnostics, not new gates or calibrated significance tests after selection.',
        'They do not measure pilot uncertainty or certify tails absent from the realized sample.', '',
        'The shared pilots permit descriptive arm comparisons, but production draws use different streams.',
        'Neither the decompositions, endpoint imbalances nor deletion sensitivities identify a causal',
        'failure mechanism. One realized production sample per arm cannot separate pilot effects,',
        'proposal geometry and production variation. No new calibration tier or coverage claim follows.', '',
        'Any covariance regularization or mixture revision requires a separate frozen protocol before',
        'evaluation, with fixed coefficients, budgets, seeds, controls and decision rules. All four',
        'current geometries are seen development cases; prospective claims require new declared cases.', '',
        '## Inputs and reproduction', '',
        f"Completed study SHA256: `{result['provenance']['completed_study']['sha256']}`.",
        f"Historical numerical/protocol commit: `{result['provenance']['implementation_commit']}`.",
        'Source bytes and the complete archive set are verified against the pinned completed report.',
        'Endpoint arrays, proposal fits, full mixture densities, scores, acceptance decisions and',
        'original assessments are recomputed without sampling. Component IDs and stream identities',
        'are validated as archived metadata; generating randomness and intermediate pilot paths are',
        'not replayed here. Archived runtime labels describe the original experiment only.',
        'Analytic truths, identities, counts, flags and selected attempt indices remain exact;',
        'derived numerical comparisons allow only the established narrow floating-point tolerance.', '',
        '    python scripts/proposal_correlated_diagnostics.py',
        '    python scripts/proposal_correlated_diagnostics.py --render-only', '',
        'Both commands reread and validate the archives. Render-only validates stored provenance and',
        'rebuilds all diagnostics rather than trusting cached derived values. Both output strings are',
        'prepared before either report is replaced. The completed study and every original archive',
        'remain unchanged.', '',
    ]
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--render-only', action='store_true')
    args = parser.parse_args(argv)
    if args.render_only:
        result = assemble(json.loads(OUT_JSON.read_text()))
    else:
        result = assemble()
    json_text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n'
    md_text = render(result)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"Rebuilt archive-only diagnostics for {len(result['runs'])} production arms.")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
