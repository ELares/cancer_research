"""Density, inverse-CDF and known-moment checks for bounded product kernels."""

import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import truncnorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bounded_proposal import BoundedGaussianMixture  # noqa: E402


class QuantileRng:
    """Force one component and specified quantiles to check the inverse exactly."""

    def __init__(self, quantiles):
        self.quantiles = np.asarray(quantiles)

    def choice(self, count, size, p):
        return np.ones(size, dtype=int)

    def uniform(self, size):
        return np.zeros(size)

    def random(self, size):
        if not size[0]:
            return np.empty(size)
        return np.broadcast_to(self.quantiles[:, None], size).copy()


@pytest.mark.parametrize("mean,scale", [(0, 0.03), (1, 0.03), (0.3, 0.2), (0.8, 2)])
def test_full_density_matches_normalized_reference_and_integrates(mean, scale):
    proposal = BoundedGaussianMixture(1, means=[[mean]], scales=[[scale]])
    xs = np.linspace(0, 1, 101)
    reference = 0.2 + 0.8 * truncnorm.pdf(xs, -mean / scale, (1 - mean) / scale,
                                         loc=mean, scale=scale)
    assert np.exp(proposal.log_density(xs[:, None])) == pytest.approx(reference, rel=1e-12)
    mass = quad(lambda x: float(np.exp(proposal.log_density([[x]])[0])), 0, 1,
                epsabs=1e-11)[0]
    assert mass == pytest.approx(1, abs=1e-10)
    assert np.isneginf(proposal.log_density([[-1e-12], [1 + 1e-12]])).all()


def test_product_density_sums_all_components_with_each_truncation_factor():
    means = np.array([[0, 0.8], [0.7, 1]])
    scales = np.array([[0.2, 0.4], [0.1, 0.3]])
    proposal = BoundedGaussianMixture(2, uniform_weight=0.3, means=means, scales=scales)
    points = np.array([[0, 1], [0.1, 0.9], [0.7, 0.5], [1, 0]])
    expected = np.full(len(points), 0.3)
    for mu, sigma in zip(means, scales):
        expected += 0.35 * np.prod(truncnorm.pdf(points, -mu / sigma, (1 - mu) / sigma,
                                                loc=mu, scale=sigma), axis=1)
    assert np.exp(proposal.log_density(points)) == pytest.approx(expected, rel=1e-12)
    assert np.isneginf(proposal.log_density([[-1e-12, 0.5], [0.5, 1.00001]])).all()


@pytest.mark.parametrize("mean,scale", [(0, 1e-8), (1, 1e-8), (0.3, 0.2), (0.8, 3)])
def test_inverse_conditional_cdf_matches_reference_including_endpoints(mean, scale):
    probabilities = np.array([0, 0.001, 0.1, 0.5, 0.9, 0.999, 1])
    proposal = BoundedGaussianMixture(1, means=[[mean]], scales=[[scale]])
    points, ids = proposal.sample(QuantileRng(probabilities), len(probabilities))
    expected = truncnorm.ppf(probabilities, -mean / scale, (1 - mean) / scale,
                            loc=mean, scale=scale)
    assert points[:, 0] == pytest.approx(expected, rel=1e-9, abs=1e-15)
    assert points[0, 0] == 0 and points[-1, 0] == 1
    assert (ids == 1).all()


@pytest.mark.parametrize("mean", [0, 0.3, 1])
@pytest.mark.parametrize("scale", [1e8, 1e100, 1e308])
def test_extremely_wide_kernels_have_stable_uniform_limit(mean, scale):
    probabilities = np.array([0, 0.01, 0.25, 0.5, 0.75, 0.99, 1])
    proposal = BoundedGaussianMixture(1, means=[[mean]], scales=[[scale]])
    # Ordinary normal-CDF subtraction is zero for the widest cases. Both the
    # conditional density and inverse must still approach the uniform law.
    assert proposal.log_density(probabilities[:, None]) == pytest.approx(
        np.zeros(len(probabilities)), abs=2e-13)
    points, _ = proposal.sample(QuantileRng(probabilities), len(probabilities))
    assert points[:, 0] == pytest.approx(probabilities, abs=1e-13)


def test_samples_recover_known_mixture_moments_and_component_shares():
    means = np.array([[0, 0.75], [0.65, 1]])
    scales = np.array([[0.12, 0.15], [0.3, 0.2]])
    proposal = BoundedGaussianMixture(2, uniform_weight=0.2, means=means, scales=scales)
    points, ids = proposal.sample(np.random.default_rng(394), 60000)
    assert points.shape == (60000, 2)
    assert np.isfinite(points).all() and ((points >= 0) & (points <= 1)).all()
    for k, p in enumerate([0.2, 0.4, 0.4]):
        assert np.mean(ids == k) == pytest.approx(p, abs=0.01)
    km, kv = truncnorm.stats(-means / scales, (1 - means) / scales,
                             loc=means, scale=scales, moments="mv")
    expected_mean = 0.2 * 0.5 + 0.8 * np.mean(km, axis=0)
    expected_second = 0.2 / 3 + 0.8 * np.mean(kv + km**2, axis=0)
    expected_cross = 0.2 / 4 + 0.8 * np.mean(np.prod(km, axis=1))
    assert points.mean(axis=0) == pytest.approx(expected_mean, abs=0.004)
    assert np.mean(points**2, axis=0) == pytest.approx(expected_second, abs=0.004)
    assert np.mean(np.prod(points, axis=1)) == pytest.approx(expected_cross, abs=0.004)


def test_serialization_preserves_density_and_exact_rng_replay_without_input_aliasing():
    means, scales = np.array([[0.2, 0.9]]), np.array([[0.1, 0.15]])
    proposal = BoundedGaussianMixture(2, means=means, scales=scales)
    means[:] = 0.5
    scales[:] = 9
    assert proposal.means.tolist() == [[0.2, 0.9]]
    assert proposal.scales.tolist() == [[0.1, 0.15]]
    with pytest.raises(ValueError):
        proposal.means[0, 0] = 0.8
    with pytest.raises(ValueError):
        proposal.scales[0, 0] = 0.8
    restored = BoundedGaussianMixture.from_dict(json.loads(json.dumps(proposal.to_dict(), allow_nan=False)))
    first, ids = proposal.sample(np.random.default_rng(187), 500)
    second, ids_again = restored.sample(np.random.default_rng(187), 500)
    assert np.array_equal(first, second) and np.array_equal(ids, ids_again)
    assert np.array_equal(proposal.log_density(first), restored.log_density(first))


@pytest.mark.parametrize("explicit_component", [False, True])
def test_pure_uniform_and_empty_batches(explicit_component):
    kwargs = {"means": [[0.2, 0.8]], "scales": [[0.1, 0.1]]} if explicit_component else {}
    proposal = BoundedGaussianMixture(2, uniform_weight=1, **kwargs)
    points, ids = proposal.sample(np.random.default_rng(41), 50)
    assert proposal.log_density(points).tolist() == [0] * 50 and (ids == 0).all()
    assert proposal.log_density([[0, 1]]).tolist() == [0]
    assert np.isneginf(proposal.log_density([[0, 1.001]])).all()
    empty, empty_ids = proposal.sample(np.random.default_rng(1), 0)
    assert empty.shape == (0, 2) and empty_ids.shape == (0,)
    assert proposal.log_density(empty).shape == (0,)


def test_defensive_mass_bounds_importance_weights_everywhere_in_cube():
    proposal = BoundedGaussianMixture(2, means=[[0, 1]], scales=[[1e-8, 1e-8]])
    points = np.array([[0, 1], [1, 0], [0.5, 0.5]])
    log_q = proposal.log_density(points)
    assert np.isfinite(log_q).all()
    assert np.all(np.exp(-log_q) <= 5 * (1 + 1e-14))


@pytest.mark.parametrize("kwargs", [
    {"dimension": 0}, {"dimension": True}, {"dimension": 1.5},
    {"dimension": 1, "uniform_weight": 0}, {"dimension": 1, "uniform_weight": np.nan},
    {"dimension": 1, "means": [[0.2]], "scales": []},
    {"dimension": 2, "means": [[0.2]], "scales": [[0.1]]},
    {"dimension": 1, "means": [[-0.01]], "scales": [[0.1]]},
    {"dimension": 1, "means": [[1.01]], "scales": [[0.1]]},
    {"dimension": 1, "means": [[np.nan]], "scales": [[0.1]]},
    {"dimension": 1, "means": [[0.2]], "scales": [[0]]},
    {"dimension": 1, "means": [[0.2]], "scales": [[-1]]},
    {"dimension": 1, "means": [[0.2]], "scales": [[np.inf]]},
    {"dimension": 1, "means": [[0.2]], "scales": [[np.nan]]},
])
def test_invalid_kernel_parameters_are_rejected(kwargs):
    with pytest.raises(ValueError):
        BoundedGaussianMixture(**kwargs)


@pytest.mark.parametrize("points", [[[np.nan]], [[np.inf]], [[0.2, 0.4]], [0.2]])
def test_invalid_density_inputs_are_rejected(points):
    with pytest.raises(ValueError):
        BoundedGaussianMixture(1).log_density(points)


@pytest.mark.parametrize("n", [-1, 1.1, True])
def test_invalid_sample_count_is_rejected(n):
    with pytest.raises(ValueError):
        BoundedGaussianMixture(1).sample(np.random.default_rng(1), n)
