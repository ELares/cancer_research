"""Exact-density defensive mixtures of Gaussian product kernels on a unit cube.

Each Gaussian coordinate is independently conditioned to [0, 1]. Its mass is
computed with erf terms of the same sign, avoiding cancellation for wide scales.
Sampling uses the inverse conditional CDF; no clipping or rejection is used.
"""

import numpy as np
from scipy.special import erf, erfinv, logsumexp


class BoundedGaussianMixture:
    """Uniform-cube component plus equally weighted truncated product normals.

    ``means`` and ``scales`` have shape (number of kernels, dimension). Centers
    must lie inside the closed unit cube; scales must be positive and finite.
    Input arrays are copied and made read-only. With no kernels, the proposal
    is pure uniform. Component IDs are zero for uniform and 1..K for kernels.
    """

    def __init__(self, dimension, uniform_weight=0.2, means=None, scales=None):
        if (isinstance(dimension, (bool, np.bool_)) or
                not isinstance(dimension, (int, np.integer)) or dimension < 1):
            raise ValueError("dimension must be a positive integer")
        uniform_weight = float(uniform_weight)
        if not np.isfinite(uniform_weight) or not 0 < uniform_weight <= 1:
            raise ValueError("uniform_weight must be finite and in (0, 1]")
        self.dimension = int(dimension)
        means = np.array([] if means is None else means, dtype=float, copy=True)
        scales = np.array([] if scales is None else scales, dtype=float, copy=True)
        if means.shape == (0,):
            means = means.reshape(0, self.dimension)
        if scales.shape == (0,):
            scales = scales.reshape(0, self.dimension)
        if (means.ndim != 2 or means.shape[1] != self.dimension or
                scales.shape != means.shape):
            raise ValueError("means and scales must have matching (kernels, dimension) shapes")
        if not np.isfinite(means).all() or np.any(means < 0) or np.any(means > 1):
            raise ValueError("kernel centers must be finite and inside the unit cube")
        if not np.isfinite(scales).all() or np.any(scales <= 0):
            raise ValueError("kernel scales must be finite and positive")
        self.uniform_weight = uniform_weight if len(means) else 1.0
        means.setflags(write=False)
        scales.setflags(write=False)
        self.means = means
        self.scales = scales
        # Both standardized interval endpoints lie on opposite sides of zero.
        # Phi(b)-Phi(a) = (erf(b/sqrt(2))+erf(-a/sqrt(2)))/2, with
        # addition instead of subtraction even when both CDFs round to 1/2.
        with np.errstate(over="ignore"):
            self._erf_left = erf((means / scales) / np.sqrt(2.0))
            self._erf_right = erf(((1.0 - means) / scales) / np.sqrt(2.0))
        mass = 0.5 * (self._erf_left + self._erf_right)
        if not np.isfinite(mass).all() or np.any(mass <= 0):
            raise ValueError("kernel truncation mass is not numerically representable")
        self._log_normalizer = -np.log(scales) - 0.5 * np.log(2 * np.pi) - np.log(mass)

    def sample(self, rng, n):
        """Return exactly n in-cube points and their sampled component IDs."""
        if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 0:
            raise ValueError("n must be a nonnegative integer")
        if not len(self.means) or self.uniform_weight == 1:
            return rng.uniform(size=(n, self.dimension)), np.zeros(n, dtype=int)
        probabilities = [self.uniform_weight] + [
            (1.0 - self.uniform_weight) / len(self.means)] * len(self.means)
        components = rng.choice(len(probabilities), size=n, p=probabilities)
        points = np.empty((n, self.dimension))
        uniform = components == 0
        points[uniform] = rng.uniform(size=(int(uniform.sum()), self.dimension))
        for k in range(len(self.means)):
            selected = components == k + 1
            u = rng.random(size=(int(selected.sum()), self.dimension))
            # Solve erf(z/sqrt(2)) = u*erf(b/sqrt(2))-(1-u)*erf(-a/sqrt(2)).
            # This remains stable for scales so wide that normal CDF endpoints
            # are indistinguishable in floating point.
            e = u * self._erf_right[k] - (1.0 - u) * self._erf_left[k]
            with np.errstate(invalid="ignore", over="ignore"):
                draw = self.means[k] + self.scales[k] * (np.sqrt(2.0) * erfinv(e))
            # Exact inverse-CDF endpoint handling, not clipping. Generator.random
            # can return zero; the conditional distribution's zero quantile is 0.
            draw[u == 0] = 0.0
            draw[u == 1] = 1.0
            points[selected] = draw
        if not np.isfinite(points).all() or np.any(points < 0) or np.any(points > 1):
            raise ValueError("bounded inverse CDF produced nonfinite or out-of-cube values")
        return points, components

    def log_density(self, points):
        """Evaluate the full normalized mixture; density is zero outside cube."""
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != self.dimension or not np.isfinite(points).all():
            raise ValueError("points must be a finite (n, dimension) array")
        inside = ((points >= 0) & (points <= 1)).all(axis=1)
        result = np.full(len(points), -np.inf)
        if not inside.any():
            return result
        if not len(self.means) or self.uniform_weight == 1:
            result[inside] = 0.0
            return result
        values = points[inside]
        terms = [np.full(len(values), np.log(self.uniform_weight))]
        log_share = np.log1p(-self.uniform_weight) - np.log(len(self.means))
        for mean, scale, normalizer in zip(self.means, self.scales, self._log_normalizer):
            with np.errstate(over="ignore", invalid="ignore"):
                z = (values - mean) / scale
                log_kernel = np.sum(normalizer - 0.5 * z * z, axis=1)
            if np.isnan(log_kernel).any():
                raise ValueError("bounded kernel density produced nonfinite arithmetic")
            terms.append(log_share + log_kernel)
        result[inside] = logsumexp(np.stack(terms), axis=0)
        return result

    def to_dict(self):
        return {"dimension": self.dimension, "uniform_weight": self.uniform_weight,
                "means": self.means.tolist(), "scales": self.scales.tolist()}

    @classmethod
    def from_dict(cls, value):
        return cls(dimension=value["dimension"], uniform_weight=value["uniform_weight"],
                   means=value["means"], scales=value["scales"])
