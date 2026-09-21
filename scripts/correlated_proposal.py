"""Frozen defensive mixtures with untruncated full-covariance Gaussian kernels.

The density is normalized on R^d. Draws outside the unit cube are ordinary
proposal attempts with zero target density; callers must count them, never
clip them, truncate the Gaussians, or redraw them. This module fits production
proposals from pilot endpoints and does not change any pilot or MH transition.
"""

from numbers import Real

import numpy as np

from importance_sampling import GaussianMixture


# Reuse the established density and sampling implementation without changing it.
CorrelatedGaussianMixture = GaussianMixture


FIT_PLAN = {
    "uniform_weight": 0.2,
    "kernel_neighbors": 32,
    "neighbor_gap_ratio": 3.0,
    "kernel_scale": 1.5,
    "kernel_floor": 0.005,
    "kernel_ceiling": 0.5,
    "singleton_scale": 0.04,
}


def _finite_array(value, name):
    """Reject booleans and non-real inputs before numeric coercion hides them."""
    try:
        if isinstance(value, np.ndarray) and value.dtype.kind in "fiu":
            result = np.asarray(value, dtype=float)
        else:
            objects = np.asarray(value, dtype=object)
            if any(isinstance(item, (bool, np.bool_)) or not isinstance(item, Real)
                   for item in objects.flat):
                raise ValueError(f"{name} must contain real numbers, not booleans")
            result = np.asarray(value, dtype=float)
    except (TypeError, OverflowError) as error:
        raise ValueError(f"{name} must contain finite real numbers") from error
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite real numbers")
    return result


def fit_proposal(points):
    """Fit fixed local covariances using endpoint geometry alone.

    Neighborhoods match resample_move_local.local_scales: distinct endpoints,
    coordinate standardization, stable nearest-neighbor ordering and the first
    positive-radius gap exceeding a factor of three. Duplicates do not change
    covariance geometry, but each particle keeps an equal mixture share.
    """
    points = _finite_array(points, "kernel centers")
    if (points.ndim != 2 or min(points.shape) < 1 or
            np.any((points < 0) | (points > 1))):
        raise ValueError("kernel centers must be nonempty points in the unit cube")
    dimension = points.shape[1]
    unique, inverse = np.unique(points, axis=0, return_inverse=True)
    standardizer = np.maximum(np.std(unique, axis=0), FIT_PLAN["kernel_floor"])
    scaled = unique / standardizer
    distances = np.sqrt(np.sum((scaled[:, None, :] - scaled[None, :, :])**2, axis=2))
    neighbors = np.argsort(distances, axis=1, kind="stable")[
        :, :min(FIT_PLAN["kernel_neighbors"], len(unique))]
    covariances = []
    for index, nearest in enumerate(neighbors):
        radii = distances[index, nearest]
        for j in range(1, len(nearest) - 1):
            if radii[j] > 0 and radii[j + 1] > FIT_PLAN["neighbor_gap_ratio"] * radii[j]:
                nearest = nearest[:j + 1]
                break
        if len(nearest) == 1:
            covariance = FIT_PLAN["singleton_scale"]**2 * np.eye(dimension)
        else:
            local = unique[nearest]
            centered = local - local.mean(axis=0)
            raw = FIT_PLAN["kernel_scale"]**2 * (centered.T @ centered) / len(local)
            values, vectors = np.linalg.eigh(raw)
            values = np.clip(values, FIT_PLAN["kernel_floor"]**2,
                             FIT_PLAN["kernel_ceiling"]**2)
            covariance = (vectors * values) @ vectors.T
            covariance = 0.5 * (covariance + covariance.T)
        covariances.append(covariance)
    return CorrelatedGaussianMixture(
        dimension, uniform_weight=FIT_PLAN["uniform_weight"], means=points,
        covariances=np.asarray(covariances)[inverse])
