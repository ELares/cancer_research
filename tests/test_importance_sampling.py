"""Distributional and numerical checks for the frozen-proposal sampler core."""

import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import multivariate_normal, norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from importance_sampling import (GaussianMixture, fit_component,  # noqa: E402
                                 importance_diagnostics, normalized_weights,
                                 weighted_quantiles)


def test_full_mixture_density_includes_every_component_and_cube_indicator():
    mixture = GaussianMixture(2, 0.2, [[0.2, 0.4], [0.8, 0.7]],
                              [[[0.1, 0.02], [0.02, 0.2]], [[0.3, 0], [0, 0.1]]])
    points = np.array([[0.4, 0.6], [0, 1], [-0.1, 0.5], [1.1, 2]])
    expected = 0.2 * ((points >= 0) & (points <= 1)).all(axis=1)
    for mean, covariance in zip(mixture.means, mixture.covariances):
        expected += 0.4 * multivariate_normal.pdf(points, mean=mean, cov=covariance)
    assert np.exp(mixture.log_density(points)) == pytest.approx(expected, rel=1e-12)
    assert np.isfinite(mixture.log_density(points)).all()


def test_untruncated_density_integrates_to_one_on_real_line():
    mixture = GaussianMixture(1, 0.2, [[0.05]], [[[0.25]]])

    def density(x):
        return float(np.exp(mixture.log_density([[x]])[0]))

    segments = [(-np.inf, 0), (0, 1), (1, np.inf)]
    masses = [quad(density, lo, hi, epsabs=1e-10)[0] for lo, hi in segments]
    assert sum(masses) == pytest.approx(1, abs=1e-9)
    expected_inside = 0.2 + 0.8 * (norm.cdf(1, 0.05, 0.5) - norm.cdf(0, 0.05, 0.5))
    assert masses[1] == pytest.approx(expected_inside, abs=1e-10)
    assert masses[0] > 0.3  # Substantial outside mass was not silently truncated.


@pytest.mark.parametrize("with_component", [False, True])
def test_pure_uniform_has_unit_density_and_zero_weight_outside(with_component):
    kwargs = {"means": [[0.3]], "covariances": [[[0.1]]]} if with_component else {}
    mixture = GaussianMixture(1, uniform_weight=1, **kwargs)
    assert mixture.log_density([[0], [0.5], [1]]).tolist() == [0, 0, 0]
    assert np.isneginf(mixture.log_density([[-0.01], [1.01]])).all()
    points, components = mixture.sample(np.random.default_rng(12), 20)
    assert ((points >= 0) & (points <= 1)).all()
    assert (components == 0).all()
    result = importance_diagnostics(mixture.log_density(points), points[:, 0] < 0.5)
    assert result["ess"] == pytest.approx(result["n_accepted"])
    assert result["normalizer_estimate"] == pytest.approx(result["n_accepted"] / 20)


def test_sampling_retains_outside_attempts_and_round_trips_exactly():
    mixture = GaussianMixture(1, 0.2, [[-0.5], [1.5]], [[[0.01]], [[0.02]]])
    restored = GaussianMixture.from_dict(json.loads(json.dumps(mixture.to_dict())))
    points, components = mixture.sample(np.random.default_rng(61), 4000)
    points_again, components_again = restored.sample(np.random.default_rng(61), 4000)
    assert np.array_equal(points, points_again)
    assert np.array_equal(components, components_again)
    assert np.array_equal(mixture.log_density(points), restored.log_density(points))
    assert points.shape == (4000, 1)
    assert ((points < 0) | (points > 1)).sum() > 2500
    assert set(components) == {0, 1, 2}
    assert (components == 0).sum() == pytest.approx(800, abs=100)
    empty, empty_ids = mixture.sample(np.random.default_rng(1), 0)
    assert empty.shape == (0, 1) and empty_ids.shape == (0,)
    assert mixture.log_density(empty).shape == (0,)


def test_defensive_density_bounds_raw_weights_on_cube():
    mixture = GaussianMixture(2, 0.2, [[0.1, 0.9]], [[[0.002, 0], [0, 0.002]]])
    points = np.random.default_rng(9).uniform(size=(1000, 2))
    raw_weights = np.exp(-mixture.log_density(points))
    assert np.max(raw_weights) <= 5 * (1 + 1e-14)


def test_outside_attempts_enter_normalizer_and_mcse_denominators():
    # One unit-weight accepted attempt out of four, including outside attempts.
    result = importance_diagnostics(np.zeros(4), np.array([True, False, False, False]))
    assert result["n_attempts"] == 4 and result["n_accepted"] == 1
    assert result["normalizer_estimate"] == pytest.approx(0.25)
    assert result["normalizer_mcse"] == pytest.approx(0.25)
    assert result["normalizer_relative_mcse"] == pytest.approx(1)
    assert result["ess"] == pytest.approx(1)


def test_normalizer_mcse_matches_direct_independent_attempt_variance():
    raw = np.array([0.5, 1, 0, 4, 0, 2])
    accepted = raw > 0
    q = np.where(accepted, 1 / np.maximum(raw, 1e-100), 1)
    result = importance_diagnostics(np.log(q), accepted)
    assert result["normalizer_estimate"] == pytest.approx(raw.mean())
    assert result["normalizer_mcse"] == pytest.approx(raw.std(ddof=1) / np.sqrt(len(raw)))
    assert result["ess"] == pytest.approx(raw.sum() ** 2 / np.dot(raw, raw))


def test_many_accepted_attempts_do_not_imply_large_effective_sample_size():
    raw = np.ones(1000)
    raw[0] = 1e6
    result = importance_diagnostics(-np.log(raw), np.ones(1000, dtype=bool))
    assert result["n_accepted"] == 1000
    assert result["ess"] < 1.003
    assert result["max_normalized_weight"] > 0.999
    assert result["normalizer_relative_mcse"] > 0.99


def test_zero_acceptance_is_explicitly_unassessable():
    result = importance_diagnostics(np.zeros(200), np.zeros(200, dtype=bool))
    assert result == {"n_attempts": 200, "n_accepted": 0, "assessable": False,
                      "ess": 0.0, "max_normalized_weight": None,
                      "normalizer_estimate": 0.0, "normalizer_mcse": None,
                      "normalizer_relative_mcse": None}
    assert "NaN" not in json.dumps(result, allow_nan=False)


def test_one_attempt_cannot_estimate_mcse():
    result = importance_diagnostics([0], np.array([True]))
    assert result["ess"] == 1 and result["normalizer_estimate"] == 1
    assert not result["assessable"]
    assert result["normalizer_mcse"] is None


def test_known_multimodal_conditioning_recovers_mass_and_distribution():
    # The proposal deliberately concentrates on the *smaller* accepted interval.
    # A posterior made from unweighted accepted counts gets both answers wrong.
    mixture = GaussianMixture(1, 0.2, [[0.045]], [[[0.0004]]])
    points, _ = mixture.sample(np.random.default_rng(921), 60000)
    x = points[:, 0]
    left = (x >= 0.02) & (x <= 0.07)
    right = (x >= 0.6) & (x <= 0.95)
    accepted = left | right
    log_q = mixture.log_density(points)
    weights = normalized_weights(np.where(accepted, -log_q, -np.inf))
    result = importance_diagnostics(log_q, accepted)
    assert result["normalizer_estimate"] == pytest.approx(0.4, abs=5 * result["normalizer_mcse"])
    assert weights[left].sum() == pytest.approx(0.125, abs=0.01)
    assert np.dot(weights, x) == pytest.approx(0.68375, abs=0.01)
    assert abs(x[accepted].mean() - 0.68375) > 0.3
    # Exact quantiles of Uniform([.02,.07] union [.6,.95]) by interval length.
    assert weighted_quantiles(x, weights, [0.025, 0.5, 0.975]) == pytest.approx(
        [0.03, 0.75, 0.94], abs=0.01)


def test_stable_log_weights_and_inverse_cdf_with_zero_weight_rows():
    assert normalized_weights([-10000, -10000, -np.inf]).tolist() == [0.5, 0.5, 0]
    quantiles = weighted_quantiles([-100, 3, 1, 2, 100], [0, 7, 1, 2, 0],
                                  [0, 0.1, 0.2, 0.3, 0.5, 1])
    assert quantiles.tolist() == [1, 1, 2, 2, 3, 3]
    assert weighted_quantiles([0, 1], [1e308, 1e308], [0.5, 1]).tolist() == [0, 1]


def test_weighted_component_fit_uses_population_covariance_and_floor():
    mean, covariance = fit_component([[0, 0], [2, 4]], np.log([0.75, 0.25]),
                                     covariance_scale=2, covariance_floor=0.1)
    assert mean == pytest.approx([0.5, 1])
    assert covariance == pytest.approx(np.array([[1.6, 3], [3, 6.1]]))
    one_mean, one_covariance = fit_component([[0.2, 0.9]], [0])
    assert one_mean == pytest.approx([0.2, 0.9])
    assert one_covariance == pytest.approx(np.eye(2) * 1e-6)


@pytest.mark.parametrize("kwargs", [
    {"dimension": 0}, {"dimension": 1.5}, {"dimension": True},
    {"dimension": 1, "uniform_weight": 0},
    {"dimension": 1, "uniform_weight": float("nan")},
    {"dimension": 1, "means": [[0]], "covariances": []},
    {"dimension": 1, "means": [[float("inf")]], "covariances": [[[1]]]},
    {"dimension": 1, "means": [[0]], "covariances": [[[-1]]]},
    {"dimension": 2, "means": [[0, 0]], "covariances": [[[1, 0.2], [0, 1]]]},
])
def test_invalid_proposal_parameters_are_rejected(kwargs):
    with pytest.raises(ValueError):
        GaussianMixture(**kwargs)


@pytest.mark.parametrize("log_weights", [[-np.inf, -np.inf], [np.nan], [np.inf], [], [[0]]])
def test_invalid_or_all_zero_log_weights_are_rejected(log_weights):
    with pytest.raises(ValueError):
        normalized_weights(log_weights)


@pytest.mark.parametrize("log_q,accepted", [
    ([np.nan], [True]), ([np.inf], [True]), ([-np.inf], [True]),
    ([0, 0], [True]), ([0], [1]),
])
def test_invalid_diagnostic_inputs_are_rejected(log_q, accepted):
    with pytest.raises(ValueError):
        importance_diagnostics(log_q, np.asarray(accepted))


@pytest.mark.parametrize("values,weights,probs", [
    ([1, 2], [0, 0], [0.5]), ([1], [-1], [0.5]), ([np.inf], [1], [0.5]),
    ([1], [np.nan], [0.5]), ([1], [1], [1.1]), ([1], [1], [np.nan]),
])
def test_invalid_weighted_quantiles_are_rejected(values, weights, probs):
    with pytest.raises(ValueError):
        weighted_quantiles(values, weights, probs)
