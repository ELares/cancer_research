"""Numerical primitives for frozen-proposal importance sampling on a unit cube.

Gaussian components are *not* truncated. Outside-cube proposals remain attempts
with zero target weight; dropping them changes the estimated normalizing mass.
All component densities, including the defensive uniform, enter each weight.
This module supplies numerical diagnostics, not a posterior adequacy decision.
"""

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import logsumexp


def _points(values, dimension=None):
    result = np.asarray(values, dtype=float)
    if result.ndim != 2 or (dimension is not None and result.shape[1] != dimension):
        raise ValueError("points must be a two-dimensional array of the expected dimension")
    if result.shape[1] == 0 or not np.isfinite(result).all():
        raise ValueError("points must have positive dimension and finite coordinates")
    return result


def _log_weights(values):
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size == 0:
        raise ValueError("log weights must be a nonempty one-dimensional array")
    if np.isnan(result).any() or np.isposinf(result).any():
        raise ValueError("log weights must be finite or negative infinity for zero weight")
    return result


def normalized_weights(log_weights):
    """Normalize finite log weights, allowing -inf for zero-weight attempts."""
    values = _log_weights(log_weights)
    if not np.isfinite(values).any():
        raise ValueError("cannot normalize all-zero importance weights")
    weights = np.exp(values - np.max(values))
    return weights / weights.sum()


class GaussianMixture:
    """Equal-weight unrestricted Gaussians plus a uniform on [0, 1]**d.

    Component IDs are 0 for the uniform and 1..K for the Gaussians. Empty
    component lists mean a pure uniform proposal, irrespective of the supplied
    defensive weight. Input arrays are copied and made read-only.
    """

    def __init__(self, dimension, uniform_weight=0.2, means=None, covariances=None):
        if isinstance(dimension, (bool, np.bool_)) or not isinstance(dimension, (int, np.integer)) or dimension < 1:
            raise ValueError("dimension must be a positive integer")
        uniform_weight = float(uniform_weight)
        if not np.isfinite(uniform_weight) or not 0 < uniform_weight <= 1:
            raise ValueError("uniform_weight must be finite and in (0, 1]")
        means = [] if means is None else list(means)
        covariances = [] if covariances is None else list(covariances)
        if len(means) != len(covariances):
            raise ValueError("means and covariances must have matching lengths")
        self.dimension = int(dimension)
        self.uniform_weight = uniform_weight if means else 1.0
        self.means = []
        self.covariances = []
        self._cholesky = []
        self._log_normalizers = []
        for mean, covariance in zip(means, covariances):
            mean = np.array(mean, dtype=float, copy=True)
            covariance = np.array(covariance, dtype=float, copy=True)
            if mean.shape != (self.dimension,) or not np.isfinite(mean).all():
                raise ValueError("each mean must have the expected dimension and finite values")
            if covariance.shape != (self.dimension, self.dimension) or not np.isfinite(covariance).all():
                raise ValueError("each covariance must have the expected shape and finite values")
            if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=1e-14):
                raise ValueError("covariances must be symmetric")
            # Avoid drawing and evaluating from different triangles when an input
            # is symmetric only up to floating-point roundoff.
            covariance = 0.5 * covariance + 0.5 * covariance.T
            try:
                chol = np.linalg.cholesky(covariance)
            except np.linalg.LinAlgError as exc:
                raise ValueError("covariances must be positive definite") from exc
            mean.setflags(write=False)
            covariance.setflags(write=False)
            self.means.append(mean)
            self.covariances.append(covariance)
            self._cholesky.append(chol)
            self._log_normalizers.append(
                -0.5 * self.dimension * np.log(2 * np.pi) - np.log(np.diag(chol)).sum())

    def sample(self, rng, n):
        """Draw exactly n attempts, retaining Gaussian draws outside the cube."""
        if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 0:
            raise ValueError("n must be a nonnegative integer")
        if not self.means or self.uniform_weight == 1:
            return rng.uniform(size=(n, self.dimension)), np.zeros(n, dtype=int)
        probabilities = [self.uniform_weight] + [
            (1 - self.uniform_weight) / len(self.means)] * len(self.means)
        components = rng.choice(len(probabilities), size=n, p=probabilities)
        points = np.empty((n, self.dimension))
        uniform = components == 0
        points[uniform] = rng.uniform(size=(int(uniform.sum()), self.dimension))
        for k, (mean, chol) in enumerate(zip(self.means, self._cholesky), start=1):
            chosen = components == k
            points[chosen] = rng.standard_normal((int(chosen.sum()), self.dimension)) @ chol.T + mean
        if not np.isfinite(points).all():
            raise ValueError("proposal sampling produced nonfinite coordinates")
        return points, components

    def log_density(self, points):
        """Evaluate the full, untruncated mixture density at finite points."""
        points = _points(points, self.dimension)
        inside = ((points >= 0) & (points <= 1)).all(axis=1)
        log_uniform = np.where(inside, np.log(self.uniform_weight), -np.inf)
        if not self.means or self.uniform_weight == 1:
            return log_uniform
        components = [log_uniform]
        log_component_weight = np.log1p(-self.uniform_weight) - np.log(len(self.means))
        for mean, chol, normalizer in zip(self.means, self._cholesky, self._log_normalizers):
            with np.errstate(over="ignore", invalid="ignore"):
                centered = (points - mean).T
                standardized = solve_triangular(chol, centered, lower=True, check_finite=False)
                squared_distance = np.sum(standardized * standardized, axis=0)
                log_component = normalizer - 0.5 * squared_distance
            if np.isnan(log_component).any():
                raise ValueError("proposal density calculation produced nonfinite arithmetic")
            components.append(log_component_weight + log_component)
        return logsumexp(np.stack(components), axis=0)

    def to_dict(self):
        return {"dimension": self.dimension, "uniform_weight": self.uniform_weight,
                "means": [mean.tolist() for mean in self.means],
                "covariances": [cov.tolist() for cov in self.covariances]}

    @classmethod
    def from_dict(cls, value):
        return cls(dimension=value["dimension"], uniform_weight=value["uniform_weight"],
                   means=value["means"], covariances=value["covariances"])


def fit_component(points, log_weights, covariance_scale=2.0, covariance_floor=1e-6):
    """Fit a weighted mean and population covariance, then inflate and regularize.

    This fits a proposal only: the covariance is the weighted MLE, without a
    sample-variance correction, and is not a posterior covariance estimate.
    """
    points = _points(points)
    weights = normalized_weights(log_weights)
    if len(weights) != len(points):
        raise ValueError("points and log weights must have matching lengths")
    if not np.isfinite(covariance_scale) or covariance_scale <= 0:
        raise ValueError("covariance_scale must be finite and positive")
    if not np.isfinite(covariance_floor) or covariance_floor <= 0:
        raise ValueError("covariance_floor must be finite and positive")
    with np.errstate(over="ignore", invalid="ignore"):
        mean = np.sum(points * weights[:, None], axis=0)
        centered = points - mean
        covariance = covariance_scale * ((centered * weights[:, None]).T @ centered)
        covariance += covariance_floor * np.eye(points.shape[1])
    if not np.isfinite(mean).all() or not np.isfinite(covariance).all():
        raise ValueError("component fit produced nonfinite parameters")
    return mean, covariance


def weighted_quantiles(values, weights, probs):
    """Return inverse weighted empirical CDF quantiles, without interpolation."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    probs = np.asarray(probs, dtype=float)
    if values.ndim != 1 or weights.shape != values.shape or values.size == 0:
        raise ValueError("values and weights must be matching nonempty one-dimensional arrays")
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or (weights < 0).any():
        raise ValueError("values must be finite and weights must be finite and nonnegative")
    if probs.ndim != 1 or not np.isfinite(probs).all() or ((probs < 0) | (probs > 1)).any():
        raise ValueError("quantile probabilities must be a one-dimensional array in [0, 1]")
    positive = weights > 0
    if not positive.any():
        raise ValueError("cannot compute quantiles with all-zero weights")
    values, weights = values[positive], weights[positive]
    order = np.argsort(values, kind="stable")
    values, weights = values[order], weights[order]
    with np.errstate(over="ignore"):
        cumulative = np.cumsum(weights)
    if not np.isfinite(cumulative[-1]):
        cumulative = np.cumsum(weights / weights.max())
    cumulative = cumulative / cumulative[-1]
    cumulative[-1] = 1.0
    indices = np.minimum(np.searchsorted(cumulative, probs, side="left"), len(values) - 1)
    return values[indices]


def importance_diagnostics(log_q, accepted):
    """Describe weights A/q for a uniform-cube target across *all* attempts.

    The caller must mark every outside-cube proposal unaccepted. ESS is the
    inverse sum of squared normalized weights, not an independent sample count.
    MCSE uses the ordinary independent-attempt sample variance of A/q. No hits
    give an explicit unassessable result, not a misleading zero standard error.
    """
    log_q = _log_weights(log_q)
    accepted = np.asarray(accepted)
    if accepted.dtype.kind != "b" or accepted.shape != log_q.shape:
        raise ValueError("accepted must be a matching one-dimensional boolean array")
    if not np.isfinite(log_q[accepted]).all():
        raise ValueError("accepted points must have finite log proposal density")
    n = len(log_q)
    n_accepted = int(accepted.sum())
    result = {"n_attempts": n, "n_accepted": n_accepted, "assessable": False,
              "ess": 0.0, "max_normalized_weight": None,
              "normalizer_estimate": 0.0, "normalizer_mcse": None,
              "normalizer_relative_mcse": None}
    if not n_accepted:
        return result
    log_raw = np.where(accepted, -log_q, -np.inf)
    weights = normalized_weights(log_raw)
    ess = min(float(1 / np.dot(weights, weights)), float(n_accepted))
    with np.errstate(over="ignore", under="ignore"):
        normalizer = float(np.exp(logsumexp(log_raw) - np.log(n)))
    if not np.isfinite(normalizer) or normalizer == 0:
        raise ValueError("normalizer estimate is outside floating-point range")
    relative_mcse = (float(np.sqrt(max(n / ess - 1, 0) / (n - 1))) if n > 1 else None)
    result.update(assessable=n > 1, ess=ess, max_normalized_weight=float(weights.max()),
                  normalizer_estimate=normalizer,
                  normalizer_mcse=(normalizer * relative_mcse if n > 1 else None),
                  normalizer_relative_mcse=relative_mcse)
    return result
