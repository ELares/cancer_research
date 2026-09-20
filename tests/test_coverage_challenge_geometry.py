"""Independent, deterministic checks of the prospectively declared geometry.

These tests use polynomial quadrature, explicit boundary points, and small
nonstudy-seed probes of the exact sampling transforms. They do not run the
learned proposal, inspect its outcomes, or evaluate reserved oracle experiments.
"""

import itertools
import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import coverage_challenge_geometry as geometry


def _rotation_quadrature():
    nodes = np.array([-1, 1]) / math.sqrt(3)
    latent = np.array(list(itertools.product(nodes, repeat=7)))
    latent *= np.array([3 / 8, 1 / 32] + [1 / 6] * 5)
    points = latent + 0.5
    points[:, 0] = 0.5 + (latent[:, 0] + latent[:, 1]) / math.sqrt(2)
    points[:, 1] = 0.5 + (latent[:, 0] - latent[:, 1]) / math.sqrt(2)
    return points, np.full(len(points), 1 / len(points))


def _annulus_quadrature():
    nodes = (1 + np.array([-1, 1]) / math.sqrt(3)) / 2
    cartesian = np.array(list(itertools.product(nodes, repeat=6)))
    # Sixteen midpoint angles exactly integrate all declared angular harmonics
    # and quadratic raw moments, with two nodes in squared radius and each slab.
    angles = (np.arange(16) + 0.5) * 2 * math.pi / 16
    values = np.repeat(cartesian, len(angles), axis=0)
    theta = np.tile(angles, len(cartesian))
    radius = np.sqrt(1 / 16 + values[:, 0] * 5 / 64)
    points = np.empty((len(values), 7))
    points[:, 0] = 0.5 + radius * np.cos(theta)
    points[:, 1] = 0.5 + radius * np.sin(theta)
    points[:, 2:] = 3 / 8 + values[:, 1:] / 4
    return points, np.full(len(points), 1 / len(points))


def _ball_quadrature():
    # Radial density in a 7-ball is 7*r**6. Eleven Gauss-Legendre nodes integrate
    # r**14 (the squared radial CDF) times that density exactly. Axis directions
    # integrate the declared degree-two angular moments without Monte Carlo.
    nodes, weights = np.polynomial.legendre.leggauss(11)
    radial, radial_weights = (nodes + 1) / 2, weights / 2
    radial_weights *= 7 * radial**6
    directions = np.concatenate([np.eye(7), -np.eye(7)])
    unit = (radial[:, None, None] * directions[None, :, :]).reshape(-1, 7)
    weights = np.repeat(radial_weights, len(directions)) / len(directions)
    a, b = 1 / 5**7, 1 / 4**7
    return (np.concatenate([1 / 4 + unit / 5, 3 / 4 + unit / 4]),
            np.concatenate([weights * a / (a + b), weights * b / (a + b)]))


QUADRATURES = dict(zip(geometry.FIXTURES,
                       (_rotation_quadrature, _annulus_quadrature, _ball_quadrature)))


@pytest.mark.parametrize("fixture", geometry.FIXTURES)
def test_exact_moments_by_independent_deterministic_quadrature(fixture):
    points, weights = QUADRATURES[fixture]()
    assert weights.sum() == pytest.approx(1, abs=2e-14)
    observed = geometry.features(points, fixture)
    truth = geometry.truth(fixture)
    assert observed.keys() == truth["moments"].keys()
    for name, values in observed.items():
        assert np.all((values >= 0) & (values <= 1)), name
        assert weights @ values == pytest.approx(truth["moments"][name], abs=2e-14), name
    regions = geometry.membership(points, fixture)
    assert regions.keys() == truth["region_masses"].keys()
    assert np.all(np.sum(list(regions.values()), axis=0) == 1)
    for name, mask in regions.items():
        assert weights[mask].sum() == pytest.approx(truth["region_masses"][name], abs=2e-14)


def test_rotation_volume_and_strong_correlation_are_exact():
    determinant_one_volume = np.prod(2 * np.array([3 / 8, 1 / 32] + [1 / 6] * 5))
    truth = geometry.truth("rotated_box")
    assert truth["mass"] == pytest.approx(determinant_one_volume, rel=1e-14)
    variance = truth["moments"]["second_x1"] - 1 / 4
    covariance = truth["moments"]["cross_x1_x2"] - 1 / 4
    assert covariance / variance == pytest.approx(143 / 145, rel=1e-14)
    assert (3 / 8 + 1 / 32) / math.sqrt(2) < 1 / 2


def test_annulus_volume_from_area_times_slab_widths():
    expected = math.pi * ((3 / 8)**2 - (1 / 4)**2) * (1 / 4)**5
    assert geometry.truth("annular_cylinder")["mass"] == expected


def test_ball_volume_from_gamma_function_and_nonoverlap():
    unit_volume = math.pi**(7 / 2) / math.gamma(7 / 2 + 1)
    expected = unit_volume * ((1 / 5)**7 + (1 / 4)**7)
    truth = geometry.truth("unequal_balls")
    assert truth["mass"] == pytest.approx(expected, rel=1e-14)
    assert truth["region_masses"]["A"] == pytest.approx((4 / 5)**7 / (1 + (4 / 5)**7))
    assert math.sqrt(7) / 2 > 1 / 5 + 1 / 4
    assert 1 / 4 - 1 / 5 > 0
    assert 3 / 4 + 1 / 4 == 1


@pytest.mark.parametrize("fixture", geometry.FIXTURES)
def test_scalar_and_vector_scores_agree_without_threshold_dependence(fixture):
    points, _ = QUADRATURES[fixture]()
    points = np.vstack([points[:10], np.zeros(7), np.ones(7), np.full(7, 0.5)])
    vector = geometry.scores(points, fixture)
    for point, distance in zip(points, vector):
        result = geometry.evaluate(point, fixture)
        assert result == geometry.evaluate(point, fixture, threshold=0.01)
        assert result == geometry.evaluate(point, fixture, threshold=100)
        assert result == {"status": "complete", "distance": distance,
                          "distance_lower_bound": distance, "simulator_dose_calls": 0}


@pytest.mark.parametrize("fixture", geometry.FIXTURES)
def test_membership_matches_scores_and_is_disjoint(fixture):
    points = np.array(list(itertools.product([0, 1 / 4, 1 / 2, 3 / 4, 1], repeat=3)))
    points = np.column_stack([points, np.full((len(points), 4), 0.5)])
    outside = np.array([[2] * 7, [-1] * 7])
    points = np.vstack([points, outside])
    counts = np.sum(list(geometry.membership(points, fixture).values()), axis=0)
    expected = ((points >= 0) & (points <= 1)).all(axis=1) & (geometry.scores(points, fixture) <= 1)
    np.testing.assert_array_equal(counts, expected)


def test_annular_hole_inner_outer_boundaries_and_slab_boundary():
    points = np.full((6, 7), 0.5)
    points[:, 0] += [0, 1 / 4, 3 / 8, 1 / 8, 7 / 16, 5 / 16]
    points[-1, 2] = 5 / 8
    np.testing.assert_array_equal(geometry.scores(points, "annular_cylinder") <= 1,
                                  [False, True, True, False, False, True])
    points[-1, 2] = np.nextafter(5 / 8, 1)
    assert geometry.scores(points[-1:], "annular_cylinder")[0] > 1


def test_rotated_regions_use_both_signs_and_have_fixed_tie_assignment():
    latent = np.array([[-0.1, -0.01], [-0.1, 0.01], [0.1, -0.01], [0.1, 0.01], [0, 0]])
    points = np.full((5, 7), 0.5)
    points[:, 0] += (latent[:, 0] + latent[:, 1]) / math.sqrt(2)
    points[:, 1] += (latent[:, 0] - latent[:, 1]) / math.sqrt(2)
    regions = geometry.membership(points, "rotated_box")
    for index in range(4):
        assert regions[f"quadrant_{index}"][index]
    assert regions["quadrant_3"][-1]


def test_annular_sector_names_and_axis_ties_are_fixed():
    angles = np.arange(8) * math.pi / 4 + math.pi / 8
    points = np.full((8, 7), 0.5)
    points[:, 0] += 5 / 16 * np.cos(angles)
    points[:, 1] += 5 / 16 * np.sin(angles)
    regions = geometry.membership(points, "annular_cylinder")
    for index in range(8):
        np.testing.assert_array_equal(regions[f"sector_{index}"], np.arange(8) == index)
    on_axis = np.full((1, 7), 0.5)
    on_axis[0, 0] += 5 / 16
    assert geometry.membership(on_axis, "annular_cylinder")["sector_0"][0]


def test_ball_centers_and_binary_exact_boundary_have_correct_membership():
    points = np.array([[1 / 4] * 7, [3 / 4] * 7, [1 / 2] * 7, [3 / 4] * 7])
    points[-1, 0] = 1
    regions = geometry.membership(points, "unequal_balls")
    np.testing.assert_array_equal(regions["A"], [True, False, False, False])
    np.testing.assert_array_equal(regions["B"], [False, True, False, True])
    assert geometry.scores(points, "unequal_balls")[-1] == 1


def test_annular_point_immediately_below_the_positive_axis_is_assigned():
    point = np.full((1, 7), 0.5)
    point[0, 0] = 13 / 16
    point[0, 1] = np.nextafter(0.5, 0)
    regions = geometry.membership(point, "annular_cylinder")
    assert regions["sector_7"][0]
    assert sum(mask[0] for mask in regions.values()) == 1


@pytest.mark.parametrize("fixture", geometry.FIXTURES)
def test_oracle_retained_fractions_follow_the_declared_region_partition(fixture):
    truth = geometry.truth(fixture)
    retained = geometry.oracle_retained_fraction(fixture, restricted=True)
    assert geometry.oracle_retained_fraction(fixture) == 1
    if fixture == "rotated_box":
        expected = sum(truth["region_masses"][f"quadrant_{j}"] for j in (2, 3))
    elif fixture == "annular_cylinder":
        expected = sum(truth["region_masses"][f"sector_{j}"] for j in range(4))
    else:
        expected = truth["region_masses"]["B"]
    assert retained == expected
    assert retained < 0.85  # All support-hole controls fail the fixed mass gate.


@pytest.mark.parametrize("fixture", geometry.FIXTURES)
def test_empty_arrays_preserve_feature_and_region_schema(fixture):
    empty = np.empty((0, 7))
    assert geometry.scores(empty, fixture).shape == (0,)
    assert all(value.shape == (0,) for value in geometry.membership(empty, fixture).values())
    assert all(value.shape == (0,) for value in geometry.features(empty, fixture).values())


@pytest.mark.parametrize("points", [np.zeros(7), np.zeros((2, 6)), [[float("nan")] * 7],
                                  [[float("inf")] * 7]])
def test_invalid_points_are_rejected(points):
    for function in (geometry.scores, geometry.membership, geometry.features):
        with pytest.raises(ValueError, match="finite.*array"):
            function(points, "rotated_box")


def test_scalar_evaluation_requires_the_unit_cube():
    with pytest.raises(ValueError, match="unit cube"):
        geometry.evaluate([-0.1] * 7, "rotated_box")


def test_unknown_fixtures_are_rejected():
    for function in (geometry.scores, geometry.membership, geometry.features):
        with pytest.raises(ValueError, match="unknown coverage fixture"):
            function(np.zeros((1, 7)), "unknown")
    for function in (geometry.truth, geometry.oracle_retained_fraction):
        with pytest.raises(ValueError, match="unknown coverage fixture"):
            function("unknown")


@pytest.mark.parametrize("restricted", [0, 1, None, "false"])
def test_oracle_restriction_requires_a_boolean(restricted):
    with pytest.raises(ValueError, match="boolean"):
        geometry.oracle_retained_fraction("rotated_box", restricted)


@pytest.mark.parametrize("n", [-1, 1.5, True])
def test_invalid_oracle_count_is_rejected_before_sampling(n):
    with pytest.raises(ValueError, match="nonnegative integer"):
        geometry.oracle_sample(None, n, "rotated_box")


def test_invalid_oracle_randomness_is_rejected_before_sampling():
    with pytest.raises(ValueError, match="numpy Generator"):
        geometry.oracle_sample(None, 1, "rotated_box")


@pytest.mark.parametrize("fixture", geometry.FIXTURES)
@pytest.mark.parametrize("restricted", [False, True])
@pytest.mark.parametrize("n", [0, 1, 17])
def test_valid_oracle_transform_shape_support_density_and_replay(fixture, restricted, n):
    # This small API probe has no adequacy screen and uses no reserved seed.
    points, log_q = geometry.oracle_sample(np.random.default_rng(41), n, fixture, restricted)
    replay, replay_log_q = geometry.oracle_sample(np.random.default_rng(41), n, fixture, restricted)
    assert points.shape == (n, 7)
    assert log_q.shape == (n,)
    np.testing.assert_array_equal(points, replay)
    np.testing.assert_array_equal(log_q, replay_log_q)
    assert np.isfinite(points).all()
    assert np.all((points >= 0) & (points <= 1))
    regions = geometry.membership(points, fixture)
    np.testing.assert_array_equal(np.sum(list(regions.values()), axis=0), np.ones(n, dtype=int))
    expected_mass = geometry.truth(fixture)["mass"] * geometry.oracle_retained_fraction(fixture, restricted)
    np.testing.assert_allclose(np.exp(-log_q), np.full(n, expected_mass), rtol=2e-15)
    if restricted:
        if fixture == "rotated_box":
            assert not regions["quadrant_0"].any()
            assert not regions["quadrant_1"].any()
        elif fixture == "annular_cylinder":
            assert not any(regions[f"sector_{index}"].any() for index in range(4, 8))
        else:
            assert not regions["A"].any()


@pytest.mark.parametrize("restricted", [False, True])
def test_rotation_inverse_transform_recovers_independent_uniform_coordinates(restricted):
    uniforms = np.random.default_rng(43).random((9, 7))
    points, _ = geometry.oracle_sample(np.random.default_rng(43), 9, "rotated_box", restricted)
    latent = np.empty_like(points)
    latent[:, 0] = (points[:, 0] + points[:, 1] - 1) / math.sqrt(2)
    latent[:, 1] = (points[:, 0] - points[:, 1]) / math.sqrt(2)
    latent[:, 2:] = points[:, 2:] - 0.5
    recovered = (latent / np.array([3 / 8, 1 / 32] + [1 / 6] * 5) + 1) / 2
    if restricted:
        recovered[:, 0] = 2 * recovered[:, 0] - 1
    np.testing.assert_allclose(recovered, uniforms, rtol=0, atol=3e-15)


@pytest.mark.parametrize("restricted", [False, True])
def test_annular_inverse_transform_recovers_uniform_area_and_angle(restricted):
    uniforms = np.random.default_rng(47).random((9, 7))
    points, _ = geometry.oracle_sample(np.random.default_rng(47), 9, "annular_cylinder", restricted)
    centered = points[:, :2] - 0.5
    recovered = np.empty_like(points)
    recovered[:, 0] = (np.sum(centered**2, axis=1) - (1 / 4)**2) / ((3 / 8)**2 - (1 / 4)**2)
    angle = np.mod(np.arctan2(centered[:, 1], centered[:, 0]), 2 * math.pi)
    recovered[:, 1] = angle / (math.pi if restricted else 2 * math.pi)
    recovered[:, 2:] = 4 * (points[:, 2:] - 3 / 8)
    np.testing.assert_allclose(recovered, uniforms, rtol=0, atol=2e-15)


@pytest.mark.parametrize("restricted", [False, True])
def test_ball_inverse_transform_recovers_radial_cdf_and_component_draw(restricted):
    rng = np.random.default_rng(53)
    component_draws = rng.random(9)
    normals = rng.normal(size=(9, 7))
    radial_draws = rng.random(9)
    points, _ = geometry.oracle_sample(np.random.default_rng(53), 9, "unequal_balls", restricted)
    a_mass, b_mass = (1 / 5)**7, (1 / 4)**7
    expected_b = np.ones(9, dtype=bool) if restricted else component_draws >= a_mass / (a_mass + b_mass)
    centers = np.where(expected_b, 3 / 4, 1 / 4)
    scales = np.where(expected_b, 1 / 4, 1 / 5)
    local = (points - centers[:, None]) / scales[:, None]
    lengths = np.linalg.norm(local, axis=1)
    np.testing.assert_allclose(lengths**7, radial_draws, rtol=0, atol=3e-15)
    expected_directions = normals / np.linalg.norm(normals, axis=1)[:, None]
    np.testing.assert_allclose(local / lengths[:, None], expected_directions, rtol=0, atol=2e-15)
    np.testing.assert_array_equal(geometry.membership(points, "unequal_balls")["B"], expected_b)
