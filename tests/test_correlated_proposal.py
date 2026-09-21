"""Endpoint-only covariance fitting and exact unrestricted mixture checks."""

import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import multivariate_normal

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from correlated_proposal import (CorrelatedGaussianMixture, FIT_PLAN,
                                 fit_proposal)  # noqa: E402
from importance_sampling import GaussianMixture  # noqa: E402
from resample_move_local import local_scales  # noqa: E402


ENDPOINTS = np.array([[.2, .2], [.3, .6], [.7, .4], [.8, .8]])


def test_fit_reuses_established_unrestricted_density_implementation():
    assert CorrelatedGaussianMixture is GaussianMixture
    assert type(fit_proposal(ENDPOINTS)) is GaussianMixture
    assert FIT_PLAN == {
        "uniform_weight": .2, "kernel_neighbors": 32, "neighbor_gap_ratio": 3.,
        "kernel_scale": 1.5, "kernel_floor": .005, "kernel_ceiling": .5,
        "singleton_scale": .04,
    }


def test_local_covariances_are_inflated_population_covariances_with_off_diagonals():
    # All four points are within the gap rule; neither eigenvalue needs clipping.
    proposal = fit_proposal(ENDPOINTS)
    expected = np.array([[.14625, .07875], [.07875, .1125]])
    assert np.asarray(proposal.means) == pytest.approx(ENDPOINTS)
    assert np.asarray(proposal.covariances) == pytest.approx(
        np.repeat(expected[None], 4, axis=0), abs=1e-15)
    assert proposal.uniform_weight == .2


def test_fitted_mixture_density_matches_scipy_inside_and_outside_cube():
    proposal = fit_proposal(ENDPOINTS)
    points = np.array([[0, 0], [1, 1], [.2, .35], [.5, .5],
                       [-.2, .4], [.7, 1.3], [-.4, 1.2]])
    inside = ((points >= 0) & (points <= 1)).all(axis=1)
    expected = .2 * inside.astype(float)
    for mean, covariance in zip(proposal.means, proposal.covariances):
        expected += .2 * multivariate_normal.pdf(points, mean=mean, cov=covariance)
    density = np.exp(proposal.log_density(points))
    assert density == pytest.approx(expected, rel=3e-14, abs=1e-15)
    assert np.all(density[inside] > .2)
    assert np.all(density[~inside] > 0)


def test_fitted_density_integrates_over_real_line_including_outside_mass():
    proposal = fit_proposal([[0.0], [.35], [.6], [1.0]])

    def density(x):
        return float(np.exp(proposal.log_density([[x]])[0]))

    # Split at the uniform component's two discontinuities.
    left = quad(density, -np.inf, 0, epsabs=1e-10)[0]
    center = quad(density, 0, 1, epsabs=1e-10)[0]
    right = quad(density, 1, np.inf, epsabs=1e-10)[0]
    assert left + center + right == pytest.approx(1, abs=2e-10)
    assert left + right > .1
    assert center < .9


def test_samples_recover_fitted_mixture_mean_covariance_and_component_shares():
    proposal = fit_proposal(ENDPOINTS)
    points, ids = proposal.sample(np.random.default_rng(432187), 80000)
    means = np.asarray(proposal.means)
    covariances = np.asarray(proposal.covariances)
    expected_mean = .2 * np.full(2, .5) + .8 * means.mean(axis=0)
    uniform_second = np.array([[1 / 3, 1 / 4], [1 / 4, 1 / 3]])
    expected_second = (.2 * uniform_second
                       + .8 * np.mean(covariances + means[:, :, None] * means[:, None, :], axis=0))
    expected_covariance = expected_second - np.outer(expected_mean, expected_mean)
    assert points.mean(axis=0) == pytest.approx(expected_mean, abs=.005)
    assert np.cov(points, rowvar=False, bias=True) == pytest.approx(expected_covariance, abs=.004)
    assert np.bincount(ids, minlength=5) / len(ids) == pytest.approx(np.full(5, .2), abs=.005)
    assert len(points) == len(ids) == 80000
    assert np.sum(np.any((points < 0) | (points > 1), axis=1)) > 10000
    assert not np.any((points == 0) | (points == 1))


class OutsideRng:
    """Provide a fixed component schedule and tails that must stay attempts."""

    def __init__(self):
        self.normal_counts = []
        self.uniform_counts = []

    def choice(self, count, size, p):
        assert count == 3 and size == 5
        assert p == pytest.approx([.2, .4, .4])
        return np.array([1, 2, 0, 1, 2])

    def uniform(self, size):
        self.uniform_counts.append(size[0])
        return np.full(size, .5)

    def standard_normal(self, size):
        self.normal_counts.append(size[0])
        assert len(self.normal_counts) <= 2, "an outside attempt was redrawn"
        return np.full(size, 20.0 if len(self.normal_counts) == 1 else -20.0)


def test_outside_draws_are_retained_once_with_original_component_ids():
    proposal = fit_proposal([[0, 0], [1, 1]])
    rng = OutsideRng()
    points, ids = proposal.sample(rng, 5)
    outside = np.any((points < 0) | (points > 1), axis=1)
    assert ids.tolist() == [1, 2, 0, 1, 2]
    assert outside.tolist() == [True, True, False, True, True]
    assert points.shape == (5, 2)
    assert rng.normal_counts == [2, 2] and rng.uniform_counts == [1]
    assert np.isfinite(proposal.log_density(points)).all()


@pytest.mark.parametrize("points,expected_largest", [
    ([[.25, .25], [.5, .5], [.75, .75]], .1875),
    ([[0, 0], [1, 1]], .25),
])
def test_singular_local_clouds_receive_declared_spectral_floor_and_ceiling(points, expected_largest):
    proposal = fit_proposal(points)
    eigenvalues = np.linalg.eigvalsh(proposal.covariances)
    assert eigenvalues[:, 0] == pytest.approx(np.full(len(points), .005**2), abs=3e-17)
    assert eigenvalues[:, 1] == pytest.approx(np.full(len(points), expected_largest), abs=1e-15)
    np.linalg.cholesky(proposal.covariances)


@pytest.mark.parametrize("point", [[.2, .7], [0, 1]])
def test_singleton_and_all_duplicate_clouds_use_fixed_isotropic_fallback(point):
    points = np.repeat([point], 20, axis=0)
    proposal = fit_proposal(points)
    assert np.asarray(proposal.covariances) == pytest.approx(
        np.repeat((.04**2 * np.eye(2))[None], 20, axis=0), abs=1e-16)
    assert len(proposal.means) == 20


def test_duplicate_particles_preserve_mixture_multiplicity_without_shrinking_neighbors():
    unique = np.concatenate([[[.10], [.12]], np.linspace(.75, .95, 40)[:, None]])
    points = np.concatenate([unique, np.repeat([[.10], [.12]], 30, axis=0)])
    proposal = fit_proposal(points)
    base = fit_proposal(unique)
    covariance = np.asarray(proposal.covariances)
    assert len(proposal.means) == len(points)
    assert np.asarray(proposal.means) == pytest.approx(points)
    assert covariance[:len(unique)] == pytest.approx(np.asarray(base.covariances), abs=1e-16)
    assert covariance[:2, 0, 0] == pytest.approx([.015**2, .015**2], abs=1e-16)
    assert covariance[:, 0, 0] == pytest.approx(local_scales(points)[:, 0]**2, abs=1e-15)
    probe = np.array([[.11], [.85]])
    reference = .2 + sum(.8 / len(points) * multivariate_normal.pdf(probe, mean=mean, cov=cov)
                         for mean, cov in zip(proposal.means, proposal.covariances))
    assert np.exp(proposal.log_density(probe)) == pytest.approx(reference, rel=2e-14)
    assert proposal.log_density(probe)[0] > base.log_density(probe)[0]


def test_neighborhood_uses_32_nearest_unique_centers_with_stable_ties():
    center = np.array([.5, .5])
    inner = [center + direction * radius
             for radius in np.arange(1, 8) / 20
             for direction in np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])]
    points = np.array([center, *inner, [0, 0], [0, 1], [1, 0], [1, 1]])
    proposal = fit_proposal(points)
    # The four corners tie for last place. Stable unique-center ordering drops
    # only (1, 1), giving negative cross covariance; dropping (1, 0) instead
    # would flip its sign. The inner/corner distance gap stays below three.
    expected = 1.5**2 * np.cov(points[:32], rowvar=False, bias=True)
    assert proposal.covariances[0] == pytest.approx(expected, abs=1e-16)
    assert proposal.covariances[0][0, 1] < 0
    assert abs(np.cov(points, rowvar=False, bias=True)[0, 1]) < 1e-16


def test_unique_center_standardization_uses_coordinate_floor_before_neighbor_selection():
    # The tiny y-coordinate differences must not be expanded to unit variance.
    points = np.column_stack([np.arange(40) / 40,
                              np.where(np.arange(40) % 2, .50001, .49999)])
    proposal = fit_proposal(points)
    local = points[:32]
    raw = 1.5**2 * np.cov(local, rowvar=False, bias=True)
    values, vectors = np.linalg.eigh(raw)
    expected = (vectors * np.maximum(values, .005**2)) @ vectors.T
    assert proposal.covariances[0] == pytest.approx(expected, abs=1e-15)


def test_fit_copies_endpoints_and_serialization_replays_frozen_proposal():
    points = ENDPOINTS.copy()
    proposal = fit_proposal(points)
    points[:] = .5
    assert np.asarray(proposal.means) == pytest.approx(ENDPOINTS)
    with pytest.raises(ValueError):
        proposal.means[0][0] = .5
    with pytest.raises(ValueError):
        proposal.covariances[0][0, 0] = .5
    state = json.loads(json.dumps(proposal.to_dict(), allow_nan=False))
    restored = CorrelatedGaussianMixture.from_dict(state)
    state["means"][0][0] = .9
    state["covariances"][0][0][0] = 9
    first, ids = proposal.sample(np.random.default_rng(29871), 200)
    replay, replay_ids = restored.sample(np.random.default_rng(29871), 200)
    assert np.array_equal(first, replay) and np.array_equal(ids, replay_ids)
    assert np.array_equal(proposal.log_density(first), restored.log_density(first))
    assert np.asarray(restored.means) == pytest.approx(ENDPOINTS)


@pytest.mark.parametrize("points", [
    [], [0.2], [[.2], [.2, .5]], np.empty((0, 2)), np.empty((2, 0)),
    [[-.01]], [[1.01]], [[np.nan]], [[np.inf]], [[-np.inf]],
    [[True]], [[.5, False]], np.array([[False]], dtype=bool),
    [[".5"]], [[1j]], None,
])
def test_invalid_or_boolean_pilot_endpoints_are_rejected(points):
    with pytest.raises(ValueError):
        fit_proposal(points)
