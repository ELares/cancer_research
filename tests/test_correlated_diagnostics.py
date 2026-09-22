"""Archive-only feature diagnostics, checked against independent hand arithmetic."""

import copy
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import proposal_correlated_diagnostics as diagnostics  # noqa: E402
import proposal_correlated_study as study  # noqa: E402
import resample_move  # noqa: E402


@pytest.fixture(scope='module', autouse=True)
def forbid_new_randomness():
    """Even archive replay must not draw again in this retrospective report."""
    def forbidden(*args, **kwargs):
        pytest.fail('retrospective diagnostics must never draw or train a pilot')

    with pytest.MonkeyPatch.context() as patch:
        for name in ('default_rng', 'SeedSequence', 'RandomState', 'Generator',
                     'random', 'uniform', 'normal', 'seed'):
            patch.setattr(np.random, name, forbidden)
        patch.setattr(study.GaussianMixture, 'sample', forbidden)
        patch.setattr(study.BoundedGaussianMixture, 'sample', forbidden)
        patch.setattr(study.strategy, 'train_pilot', forbidden)
        patch.setattr(resample_move, 'train_pilot', forbidden)
        patch.setattr(study.geometry, 'oracle_sample', forbidden)
        patch.setattr(study.geometry.baseline, 'oracle_sample', forbidden)
        for name in ('run_learned', 'run_control', 'assess_archive', 'assemble', 'build_report'):
            patch.setattr(study, name, forbidden)
        patch.setattr(study, 'validate_pilot', forbidden)
        yield


def _same_numbers(actual, expected):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected:
            _same_numbers(actual[key], expected[key])
    elif isinstance(expected, (list, tuple)):
        assert len(actual) == len(expected)
        for left, right in zip(actual, expected):
            _same_numbers(left, right)
    elif type(expected) is float:
        assert actual == pytest.approx(expected, rel=2e-13, abs=2e-15)
    else:
        assert actual == expected


def test_hand_calculated_weighted_feature_mcse_and_largest_deletion():
    result = diagnostics.feature_diagnostics([0., 1., 0., 0.], [1., 3., 0., 0.], .5)
    assert result['estimate'] == .75
    assert result['signed_error'] == .25
    # Four attempts, including two rejections: (4/3) * 2 * (3/16)^2.
    assert result['conditional_mcse'] == pytest.approx(math.sqrt(3 / 32))
    assert result['error_over_mcse'] == pytest.approx(math.sqrt(2 / 3))
    assert result['single_deletion'] == {
        'attempt_index': 1, 'signed_change': -.75, 'absolute_change': .75,
    }


def test_rejections_remain_in_the_mcse_attempt_count():
    full = diagnostics.feature_diagnostics([0., 1., 0., 0.], [1., 3., 0., 0.], .5)
    accepted_only = diagnostics.feature_diagnostics([0., 1.], [1., 3.], .5)
    assert full['estimate'] == accepted_only['estimate']
    assert full['conditional_mcse'] == pytest.approx(math.sqrt(3 / 32))
    assert accepted_only['conditional_mcse'] == .375


def test_hand_calculated_weight_concentration():
    result = diagnostics.weight_concentration([1., 3., 0., 0.])
    assert result['n_attempts'] == 4
    assert result['n_positive'] == 2
    assert result['ess'] == pytest.approx(1.6)
    assert result['max_normalized_weight'] == .75
    assert result['top_1_weight'] == .75
    assert result['top_5_weight'] == result['top_20_weight'] == 1.


def test_one_hit_has_no_defined_positive_weight_deletion():
    result = diagnostics.feature_diagnostics([0., .4, 0.], [0., 3., 0.], .5)
    assert result['estimate'] == .4
    assert result['signed_error'] == pytest.approx(-.1)
    assert result['conditional_mcse'] is None
    assert result['error_over_mcse'] is None
    assert result['single_deletion'] is None


def test_constant_feature_has_zero_mcse_and_lowest_original_index_tie():
    result = diagnostics.feature_diagnostics([.4, .4, .4], [0., 1., 1.], .5)
    assert result['conditional_mcse'] == 0.
    assert result['error_over_mcse'] is None
    assert result['single_deletion'] == {
        'attempt_index': 0, 'signed_change': 0., 'absolute_change': 0.,
    }


def test_equal_nonzero_deletion_influences_choose_lowest_original_index():
    result = diagnostics.feature_diagnostics([0., 0., 1., 0.], [0., 1., 1., 0.], .5)
    assert result['single_deletion'] == {
        'attempt_index': 1, 'signed_change': .5, 'absolute_change': .5,
    }


@pytest.mark.parametrize('weights', [[], [0., 0.], [-1., 2.], [np.nan, 1.], [np.inf, 1.], [[1., 2.]]])
def test_invalid_weight_vectors_are_rejected(weights):
    with pytest.raises(ValueError):
        diagnostics.weight_concentration(weights)


@pytest.mark.parametrize('values,weights,truth', [
    ([], [], .5), ([0., 1.], [0., 0.], .5), ([0.], [1., 2.], .5),
    ([np.nan, 1.], [1., 2.], .5), ([np.inf, 1.], [1., 2.], .5),
    ([0., 1.], [-1., 2.], .5), ([0., 1.], [1., 2.], np.nan),
    ([[0., 1.]], [1., 2.], .5),
])
def test_invalid_feature_inputs_are_rejected(values, weights, truth):
    with pytest.raises(ValueError):
        diagnostics.feature_diagnostics(values, weights, truth)


def test_scaling_and_permutation_preserve_metrics_except_attempt_indices():
    values = np.array([.1, .8, .6, .2, .9])
    weights = np.array([1., 0., 3., 2., 0.])
    order = np.array([3, 0, 4, 2, 1])
    reference = diagnostics.feature_diagnostics(values, weights, .5)
    scaled = diagnostics.feature_diagnostics(values, 37 * weights, .5)
    permuted = diagnostics.feature_diagnostics(values[order], weights[order], .5)
    _same_numbers(scaled, reference)
    assert order[permuted['single_deletion']['attempt_index']] == reference['single_deletion']['attempt_index']
    for result in (reference, permuted):
        result['single_deletion'].pop('attempt_index')
    _same_numbers(permuted, reference)
    _same_numbers(diagnostics.weight_concentration(37 * weights[order]),
                  diagnostics.weight_concentration(weights))


def _regional_fixture():
    values = np.array([0., 1., 0., 0.])
    weights = np.array([1., 3., 0., 0.])
    regions = {'left': np.array([True, False, False, False]),
               'right': np.array([False, True, False, False]),
               'empty': np.zeros(4, dtype=bool)}
    probabilities = {'left': .25, 'right': .5, 'empty': .25}
    conditional = {'left': .2, 'right': .6, 'empty': .8}
    return values, weights, regions, probabilities, conditional


def test_region_decomposition_closes_with_a_missing_region():
    result = diagnostics.decompose_feature(*_regional_fixture())
    # Truth .25*.2 + .5*.6 + .25*.8 = .55, sample mean .75.
    assert result['total_signed_error'] == pytest.approx(.2)
    assert result['between_regions'] == pytest.approx(-.05)
    assert result['within_regions'] == pytest.approx(.25)
    assert result['between_regions'] + result['within_regions'] == pytest.approx(result['total_signed_error'])
    empty = result['regions']['empty']
    assert empty['n_accepted'] == 0
    assert empty['estimated_mass'] == 0.
    assert empty['conditional_estimate'] is None
    assert empty['between_contribution'] == pytest.approx(-.2)
    assert empty['within_contribution'] == 0.
    assert result['regions']['left']['conditional_estimate'] == 0.
    assert result['regions']['right']['conditional_estimate'] == 1.


def test_equal_region_conditional_truth_has_exactly_zero_between_component():
    values, weights, regions, probabilities, _ = _regional_fixture()
    result = diagnostics.decompose_feature(values, weights, regions, probabilities,
                                           {name: .5 for name in regions})
    assert result['between_regions'] == 0.
    assert result['within_regions'] == result['total_signed_error'] == .25


def test_region_decomposition_is_invariant_to_weight_scale_and_attempt_order():
    values, weights, regions, probabilities, conditional = _regional_fixture()
    expected = diagnostics.decompose_feature(values, weights, regions, probabilities, conditional)
    order = np.array([3, 1, 0, 2])
    result = diagnostics.decompose_feature(values[order], 17 * weights[order],
                                           {name: mask[order] for name, mask in regions.items()},
                                           probabilities, conditional)
    _same_numbers(result, expected)


@pytest.mark.parametrize('fault', ['overlap', 'unassigned', 'missing_truth', 'bad_probability'])
def test_invalid_regional_partitions_or_truth_are_rejected(fault):
    values, weights, regions, probabilities, conditional = _regional_fixture()
    if fault == 'overlap':
        regions['left'][1] = True
    elif fault == 'unassigned':
        regions['right'][1] = False
    elif fault == 'missing_truth':
        del conditional['left']
    else:
        probabilities['left'] = -.25
    with pytest.raises(ValueError):
        diagnostics.decompose_feature(values, weights, regions, probabilities, conditional)


def _frozen_bytes():
    paths = [diagnostics.STUDY_JSON, *diagnostics.ARCHIVES.glob('*.json.gz'),
             *(study.ROOT / relative for relative in study.SOURCE_PATHS)]
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


@pytest.fixture(scope='module')
def report(forbid_new_randomness):
    before = _frozen_bytes()
    result = diagnostics.assemble()
    assert _frozen_bytes() == before
    return result


def test_every_frozen_record_and_arm_is_present_without_new_draws(report):
    expected = {(fixture, seed, arm) for fixture, seeds in study.LEARNED_SEEDS.items()
                for seed in seeds for arm in study.ARMS}
    assert {(run['fixture'], run['seed'], run['arm']) for run in report['runs']} == expected
    assert len(report['runs']) == 36
    assert len(report['pilots']) == 12
    assert len({run['archive'] for run in report['runs']}) == 12
    original = json.loads(diagnostics.STUDY_JSON.read_text())
    original_runs = {(run['fixture'], run['seed'], run['arm']): run
                     for run in original['positive_runs']}
    for run in report['runs']:
        key = (run['fixture'], run['seed'], run['arm'])
        assert run['original_passed'] is original_runs[key]['assessment']['passed']
        assert run['diagnostics']['weights']['n_attempts'] == 8192
        assert run['diagnostics']['weights']['n_positive'] == original_runs[key]['assessment']['importance']['n_accepted']
        assert run['diagnostics']['weights']['ess'] == pytest.approx(
            original_runs[key]['assessment']['importance']['ess'], rel=1e-12)
        for name, feature in run['diagnostics']['features'].items():
            assert feature['estimate'] == pytest.approx(
                original_runs[key]['assessment']['moments'][name], rel=1e-12)


def test_pilot_endpoint_summaries_describe_the_archived_clouds(report):
    for row in report['pilots']:
        raw = json.loads(gzip.decompress((diagnostics.ARCHIVES / row['archive']).read_bytes()))
        pooled = row['diagnostics']['pooled']
        points = raw['pilot']['final_points']
        assert pooled['n_endpoints'] == len(points) == 512
        assert pooled['n_unique_endpoints'] == len({tuple(point) for point in points})
        assert sum(pooled['region_counts'].values()) == len(points)
        assert len(row['diagnostics']['islands']) == 2
        for summary in row['diagnostics']['islands']:
            island = raw['pilot']['islands'][summary['island']]
            assert summary['summary']['n_endpoints'] == len(island['final_points']) == 256


def test_shifted_failure_decomposition_has_exactly_zero_between_region_component(report):
    shifted = [run for run in report['runs'] if run['fixture'] == 'shifted_rotated_box']
    assert len(shifted) == 9
    for run in shifted:
        decomposition = run['diagnostics']['decomposition']
        assert decomposition['between_regions'] == 0.
        assert decomposition['within_regions'] == pytest.approx(decomposition['total_signed_error'], abs=2e-15)


def test_stale_cached_diagnostics_are_rebuilt_deterministically(report):
    stored = copy.deepcopy(report)
    stored['runs'][0]['diagnostics']['weights']['ess'] = -123.
    stored['pilots'][0]['diagnostics']['pooled']['n_unique_endpoints'] = -1
    stored['runs'].reverse()
    assert diagnostics.assemble(stored) == report


def test_direct_render_is_stable_when_runs_pilots_and_islands_are_reversed(report):
    shuffled = copy.deepcopy(report)
    shuffled['runs'].reverse()
    shuffled['pilots'].reverse()
    for pilot in shuffled['pilots']:
        pilot['diagnostics']['islands'].reverse()
    rendered = diagnostics.render(shuffled)
    assert rendered == diagnostics.render(report)
    pilot_table = rendered.split('## Descriptive pilot endpoint summaries', 1)[1].split(
        '## Numerical definitions and limits', 1)[0]
    for pilot in report['pilots']:
        if pilot['fixture'] not in report['selected_features']:
            continue
        feature = report['selected_features'][pilot['fixture']]['feature']
        prefix = f"| {pilot['fixture']} | {pilot['seed']} |"
        line = next(line for line in pilot_table.splitlines() if line.startswith(prefix))
        cells = [cell.strip() for cell in line.strip('|').split('|')]
        islands = {row['island']: row['summary'] for row in pilot['diagnostics']['islands']}
        assert cells[3] == f"{islands[0]['feature_means'][feature]:.6f}"
        assert cells[4] == f"{islands[1]['feature_means'][feature]:.6f}"


def test_highlighted_feature_with_unavailable_mcse_renders_explicitly(report):
    altered = copy.deepcopy(report)
    run = next(run for run in altered['runs'] if run['fixture'] in altered['selected_features'])
    feature = altered['selected_features'][run['fixture']]['feature']
    run['diagnostics']['features'][feature]['conditional_mcse'] = None
    run['diagnostics']['features'][feature]['error_over_mcse'] = None
    rendered = diagnostics.render(altered)
    feature_table = rendered.split('## Highlighted feature errors across every paired pilot', 1)[1].split(
        '## Exact decomposition of the selected errors', 1)[0]
    prefix = f"| {run['fixture']} | {run['seed']} | {run['arm']} |"
    line = next(line for line in feature_table.splitlines() if line.startswith(prefix))
    cells = [cell.strip() for cell in line.strip('|').split('|')]
    assert cells[5] == 'unavailable'
    assert cells[3] == f"{run['diagnostics']['features'][feature]['estimate']:.6f}"


@pytest.fixture
def isolated_outputs(monkeypatch, tmp_path):
    out_json, out_md = tmp_path / 'diagnostics.json', tmp_path / 'diagnostics.md'
    out_json.write_text('original JSON')
    out_md.write_text('original Markdown')
    monkeypatch.setattr(diagnostics, 'OUT_JSON', out_json)
    monkeypatch.setattr(diagnostics, 'OUT_MD', out_md)
    return out_json, out_md


@pytest.mark.parametrize('fault', ['missing', 'extra', 'modified'])
def test_invalid_archive_inventory_preserves_existing_outputs(monkeypatch, tmp_path, isolated_outputs, fault):
    archive_dir = tmp_path / 'archives'
    archive_dir.mkdir()
    originals = list(diagnostics.ARCHIVES.glob('*.json.gz'))
    for original in originals:
        shutil.copyfile(original, archive_dir / original.name)
    first = archive_dir / originals[0].name
    if fault == 'missing':
        first.unlink()
    elif fault == 'extra':
        (archive_dir / 'extra.json.gz').write_bytes(b'extra')
    else:
        first.write_bytes(first.read_bytes() + b'modified')
    monkeypatch.setattr(diagnostics, 'ARCHIVES', archive_dir)
    with pytest.raises(ValueError):
        diagnostics.main([])
    assert tuple(path.read_text() for path in isolated_outputs) == ('original JSON', 'original Markdown')


def test_tampered_study_report_preserves_existing_outputs(monkeypatch, tmp_path, isolated_outputs):
    altered = json.loads(diagnostics.STUDY_JSON.read_text())
    altered['implementation_commit'] = 'b' * 40
    path = tmp_path / 'altered-study.json'
    path.write_text(json.dumps(altered))
    monkeypatch.setattr(diagnostics, 'STUDY_JSON', path)
    with pytest.raises(ValueError):
        diagnostics.main([])
    assert tuple(path.read_text() for path in isolated_outputs) == ('original JSON', 'original Markdown')


def test_tampered_cached_provenance_preserves_existing_outputs(report, isolated_outputs):
    stored = copy.deepcopy(report)
    stored['provenance']['unexpected'] = 'changed'
    isolated_outputs[0].write_text(json.dumps(stored))
    before = tuple(path.read_bytes() for path in isolated_outputs)
    with pytest.raises(ValueError):
        diagnostics.main(['--render-only'])
    assert tuple(path.read_bytes() for path in isolated_outputs) == before


def test_render_failure_preserves_both_existing_outputs(monkeypatch, report, isolated_outputs):
    monkeypatch.setattr(diagnostics, 'assemble', lambda *args: copy.deepcopy(report))

    def broken_render(result):
        raise ValueError('render failed')

    monkeypatch.setattr(diagnostics, 'render', broken_render)
    with pytest.raises(ValueError, match='render failed'):
        diagnostics.main([])
    assert tuple(path.read_text() for path in isolated_outputs) == ('original JSON', 'original Markdown')


def test_render_only_rebuilds_cached_values_without_touching_frozen_inputs(report, isolated_outputs):
    stale = copy.deepcopy(report)
    stale['runs'][0]['diagnostics']['weights']['ess'] = -123.
    isolated_outputs[0].write_text(json.dumps(stale))
    before = _frozen_bytes()
    assert diagnostics.main(['--render-only']) == 0
    assert json.loads(isolated_outputs[0].read_text()) == report
    assert isolated_outputs[1].read_text() == diagnostics.render(report)
    assert _frozen_bytes() == before
