"""Analytic geometry checks, with no pilot training or production experiments.

Deterministic quadrature checks the declared truth. Small fixed-seed oracle
probes check transforms and replay; they never apply experimental pass criteria.
"""

import itertools
import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import correlated_challenge_geometry as geometry
import coverage_challenge_geometry as baseline

FIXTURE = "shifted_rotated_box"


def _parameters():
    center = np.array([0.18, 0.18, 0.18, 0.5, 0.5, 0.5, 0.5])
    widths = np.array([0.25, 1 / 48, 1 / 32, 0.25, 0.25, 0.25, 0.25])
    basis = np.eye(7)
    basis[:3, :3] = [[1 / math.sqrt(3), 1 / math.sqrt(2), 1 / math.sqrt(6)],
                      [1 / math.sqrt(3), -1 / math.sqrt(2), 1 / math.sqrt(6)],
                      [1 / math.sqrt(3), 0, -2 / math.sqrt(6)]]
    return center, widths, basis


def _transform(normalized):
    center, widths, basis = _parameters()
    return center + (np.asarray(normalized) * widths) @ basis.T


def test_declared_rotation_volume_and_all_corner_bounds():
    center, widths, basis = _parameters()
    np.testing.assert_allclose(basis.T @ basis, np.eye(7), atol=3e-16)
    assert np.linalg.det(basis) == pytest.approx(1, abs=3e-16)
    assert geometry.truth(FIXTURE)["mass"] == pytest.approx(
        abs(np.linalg.det(basis)) * np.prod(2 * widths), rel=3e-16)
    assert geometry.truth(FIXTURE)["mass"] == 1 / 12288
    corners = _transform(list(itertools.product([-1, 1], repeat=7)))
    extents = np.abs(basis) @ widths
    np.testing.assert_allclose(corners.min(axis=0), center - extents, atol=6e-17)
    np.testing.assert_allclose(corners.max(axis=0), center + extents, atol=6e-17)
    assert np.all(corners > 0)
    assert np.all(corners < 1)
    assert 0 < corners[:, :3].min() < 0.01
    # The irrational rotation can round a mathematical boundary to either side
    # of one. Check boundary scores numerically, without broadening acceptance.
    np.testing.assert_allclose(geometry.scores(corners, FIXTURE), 1, atol=3e-15, rtol=0)


def test_polynomial_moments_and_regions_by_independent_tensor_quadrature():
    nodes = np.array([-1, 1]) / math.sqrt(3)
    points = _transform(list(itertools.product(nodes, repeat=7)))
    observed = geometry.features(points, FIXTURE)
    truth = geometry.truth(FIXTURE)
    assert observed.keys() == truth["moments"].keys()
    assert len(observed) == 69
    for name, values in observed.items():
        assert np.all((values >= 0) & (values <= 1)), name
        if "cdf" not in name:
            assert values.mean() == pytest.approx(truth["moments"][name], abs=2e-15), name
    regions = geometry.membership(points, FIXTURE)
    assert regions.keys() == truth["region_masses"].keys()
    np.testing.assert_array_equal(np.sum(list(regions.values()), axis=0), np.ones(len(points)))
    for name, mask in regions.items():
        assert mask.mean() == truth["region_masses"][name]


def test_raw_covariance_moments_have_independent_closed_forms():
    moments = geometry.truth(FIXTURE)["moments"]
    long_variance = (1 / 4) ** 2 / 3
    short_variance = (1 / 48) ** 2 / 3
    third_variance = (1 / 32) ** 2 / 3
    diagonal = long_variance / 3 + short_variance / 2 + third_variance / 6
    covariance12 = long_variance / 3 - short_variance / 2 + third_variance / 6
    covariance23 = long_variance / 3 - third_variance / 3
    for name in ("second_x1", "second_x2"):
        assert moments[name] == pytest.approx((9 / 50) ** 2 + diagonal, abs=1e-17)
    assert moments["cross_x1_x2"] == pytest.approx((9 / 50) ** 2 + covariance12, abs=1e-17)
    assert moments["cross_x2_x3"] == pytest.approx((9 / 50) ** 2 + covariance23, abs=1e-17)


def test_latent_cdf_indicators_by_exact_equal_volume_bin_quadrature():
    # Four equal-volume bins per independent coordinate exactly integrate the
    # three declared indicators. Midpoints avoid numerical cut-point ties.
    nodes = [-0.75, -0.25, 0.25, 0.75]
    points = _transform(list(itertools.product(nodes, repeat=7)))
    observed = geometry.features(points, FIXTURE)
    truth = geometry.truth(FIXTURE)["moments"]
    for j in range(1, 8):
        for suffix in ("25", "50", "75"):
            name = f"latent_cdf_{j}_{suffix}"
            assert set(np.unique(observed[name])) == {0, 1}
            assert observed[name].mean() == truth[name]


def test_box_radial_cdf_by_independent_shell_quadrature():
    # Volume inside normalized max-radius r is r**7, giving shell density
    # 7*r**6. Eleven Gauss-Legendre nodes integrate r**14 times this density.
    nodes, weights = np.polynomial.legendre.leggauss(11)
    radii = (nodes + 1) / 2
    weights = weights / 2 * 7 * radii ** 6
    latent = np.zeros((len(radii), 7))
    latent[:, 0] = radii
    observed = geometry.features(_transform(latent), FIXTURE)
    for name, expected in (("radial_cdf", 1 / 2), ("radial_cdf_second", 1 / 3)):
        assert weights @ observed[name] == pytest.approx(expected, abs=3e-15)
        assert geometry.truth(FIXTURE)["moments"][name] == expected


def test_octants_and_center_tie_assignment():
    latent = np.zeros((9, 7))
    latent[:8, :3] = list(itertools.product([-0.5, 0.5], repeat=3))
    regions = geometry.membership(_transform(latent), FIXTURE)
    for index in range(8):
        np.testing.assert_array_equal(regions[f"octant_{index}"][:8], np.arange(8) == index)
    assert regions["octant_7"][-1]
    center = geometry.features(_transform(np.zeros((1, 7))), FIXTURE)
    for j in range(1, 8):
        assert center[f"latent_cdf_{j}_25"][0] == 0
        assert center[f"latent_cdf_{j}_50"][0] == 1
        assert center[f"latent_cdf_{j}_75"][0] == 1


def test_cube_membership_and_scalar_score_contract():
    points = np.vstack([_transform(np.zeros((1, 7))), np.zeros((1, 7)),
                        np.ones((1, 7)), np.full((1, 7), -0.5), np.full((1, 7), 1.5)])
    scores = geometry.scores(points, FIXTURE)
    assert np.isfinite(scores).all()
    expected = ((points >= 0) & (points <= 1)).all(axis=1) & (scores <= 1)
    regions = geometry.membership(points, FIXTURE)
    np.testing.assert_array_equal(np.sum(list(regions.values()), axis=0), expected)
    for point, score in zip(points[:3], scores[:3]):
        expected_result = {"status": "complete", "distance": score,
                           "distance_lower_bound": score, "simulator_dose_calls": 0}
        for threshold in (None, 0.01, 100):
            assert geometry.evaluate(point, FIXTURE, threshold) == expected_result
    for point in points[3:]:
        with pytest.raises(ValueError, match="unit cube"):
            geometry.evaluate(point, FIXTURE)


@pytest.mark.parametrize("fixture", baseline.FIXTURES)
def test_original_geometry_and_oracles_delegate_without_changes(fixture):
    points = np.array([[0.25] * 7, [0.5] * 7, [0.75] * 7, [-0.5] * 7])
    np.testing.assert_array_equal(geometry.scores(points, fixture), baseline.scores(points, fixture))
    for function in ("membership", "features"):
        wrapped = getattr(geometry, function)(points, fixture)
        original = getattr(baseline, function)(points, fixture)
        assert wrapped.keys() == original.keys()
        for name in wrapped:
            np.testing.assert_array_equal(wrapped[name], original[name])
    assert geometry.truth(fixture) == baseline.truth(fixture)
    assert geometry.evaluate(points[0], fixture, 0.2) == baseline.evaluate(points[0], fixture, 0.2)
    for restricted in (False, True):
        assert geometry.oracle_retained_fraction(fixture, restricted) == baseline.oracle_retained_fraction(
            fixture, restricted)
        actual = geometry.oracle_sample(np.random.default_rng(61), 9, fixture, restricted)
        original = baseline.oracle_sample(np.random.default_rng(61), 9, fixture, restricted)
        for left, right in zip(actual, original):
            np.testing.assert_array_equal(left, right)


@pytest.mark.parametrize("restricted", [False, True])
@pytest.mark.parametrize("n", [0, 1, 17])
def test_oracle_transform_support_density_replay_and_retained_regions(restricted, n):
    # Fixed unit-test seeds are separate from any declared experiment seeds.
    uniforms = np.random.default_rng(67).random((n, 7))
    normalized = 2 * uniforms - 1
    if restricted:
        normalized[:, 0] = uniforms[:, 0]
    points, log_q = geometry.oracle_sample(np.random.default_rng(67), n, FIXTURE, restricted)
    replay, replay_log_q = geometry.oracle_sample(np.random.default_rng(67), n, FIXTURE, restricted)
    np.testing.assert_allclose(points, _transform(normalized), atol=6e-17, rtol=0)
    np.testing.assert_array_equal(points, replay)
    np.testing.assert_array_equal(log_q, replay_log_q)
    assert points.shape == (n, 7)
    assert log_q.shape == (n,)
    assert np.all((points >= 0) & (points <= 1))
    retained = geometry.oracle_retained_fraction(FIXTURE, restricted)
    assert retained == (0.5 if restricted else 1)
    np.testing.assert_allclose(np.exp(-log_q), np.full(n, retained / 12288), rtol=2e-15)
    regions = geometry.membership(points, FIXTURE)
    np.testing.assert_array_equal(np.sum(list(regions.values()), axis=0), np.ones(n))
    if restricted:
        assert not any(regions[f"octant_{j}"].any() for j in range(4))
        truth = geometry.truth(FIXTURE)
        assert sum(truth["region_masses"][f"octant_{j}"] for j in range(4, 8)) == retained


def test_empty_arrays_preserve_schemas_and_declared_interface():
    assert geometry.DIMENSION == 7
    assert geometry.EPSILON == 1
    assert geometry.FIXTURES == baseline.FIXTURES + (FIXTURE,)
    empty = np.empty((0, 7))
    assert geometry.scores(empty, FIXTURE).shape == (0,)
    for function, schema in ((geometry.features, "moments"), (geometry.membership, "region_masses")):
        values = function(empty, FIXTURE)
        assert values.keys() == geometry.truth(FIXTURE)[schema].keys()
        assert all(value.shape == (0,) for value in values.values())


@pytest.mark.parametrize("points", [np.zeros(7), np.zeros((2, 6)), [[float("nan")] * 7],
                                  [[float("inf")] * 7]])
def test_invalid_points_are_rejected(points):
    for function in (geometry.scores, geometry.membership, geometry.features):
        with pytest.raises(ValueError, match="finite.*array"):
            function(points, FIXTURE)


def test_unknown_fixtures_are_rejected():
    for function in (geometry.scores, geometry.membership, geometry.features):
        with pytest.raises(ValueError, match="unknown coverage fixture"):
            function(np.zeros((1, 7)), "unknown")
    for function in (geometry.truth, geometry.oracle_retained_fraction):
        with pytest.raises(ValueError, match="unknown coverage fixture"):
            function("unknown")
    with pytest.raises(ValueError, match="unknown coverage fixture"):
        geometry.evaluate(np.zeros(7), "unknown")
    with pytest.raises(ValueError, match="unknown coverage fixture"):
        geometry.oracle_sample(np.random.default_rng(67), 1, "unknown")


@pytest.mark.parametrize("restricted", [0, 1, None, "false"])
def test_restriction_requires_a_boolean(restricted):
    with pytest.raises(ValueError, match="boolean"):
        geometry.oracle_retained_fraction(FIXTURE, restricted)
    with pytest.raises(ValueError, match="boolean"):
        geometry.oracle_sample(np.random.default_rng(67), 1, FIXTURE, restricted)


@pytest.mark.parametrize("n", [-1, 1.5, True])
def test_invalid_count_is_rejected_before_sampling(n):
    with pytest.raises(ValueError, match="nonnegative integer"):
        geometry.oracle_sample(None, n, FIXTURE)


def test_invalid_randomness_is_rejected_before_sampling():
    with pytest.raises(ValueError, match="numpy Generator"):
        geometry.oracle_sample(None, 1, FIXTURE)
