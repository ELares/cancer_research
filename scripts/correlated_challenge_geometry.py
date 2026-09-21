"""Frozen analytic fixtures for a prospective correlated-proposal comparison.

The original three fixtures delegate to their unchanged geometry module. The
new target is a shifted, rotated seven-dimensional box, with exact uniform
oracles and analytic moments. Restricted oracles intentionally omit support;
they are diagnostic controls, never valid full-target importance proposals.
"""

import math
from fractions import Fraction

import numpy as np

import coverage_challenge_geometry as baseline

DIMENSION = baseline.DIMENSION
EPSILON = baseline.EPSILON
FIXTURES = baseline.FIXTURES + ("shifted_rotated_box",)

SHIFTED_CENTER = np.array([9 / 50] * 3 + [1 / 2] * 4)
SHIFTED_HALF_WIDTHS = np.array([1 / 4, 1 / 48, 1 / 32] + [1 / 4] * 4)
SHIFTED_BASIS = np.eye(DIMENSION)
SHIFTED_BASIS[:3, :3] = np.column_stack((
    np.array([1, 1, 1]) / math.sqrt(3),
    np.array([1, -1, 0]) / math.sqrt(2),
    np.array([1, 1, -2]) / math.sqrt(6),
))
CDF_THRESHOLDS = (("25", 1 / 4), ("50", 1 / 2), ("75", 3 / 4))


def _points(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != DIMENSION or not np.isfinite(points).all():
        raise ValueError("coverage points must be a finite (n, 7) array")
    return points


def _latent(points):
    return (points - SHIFTED_CENTER) @ SHIFTED_BASIS


def scores(points, fixture):
    """Complete scores, including finite points outside the prior unit cube."""
    if fixture != "shifted_rotated_box":
        return baseline.scores(points, fixture)
    points = _points(points)
    return np.max(np.abs(_latent(points)) / SHIFTED_HALF_WIDTHS, axis=1)


def evaluate(point, fixture, threshold=None):
    """Score a cube point fully; the caller's threshold does not alter target."""
    if fixture != "shifted_rotated_box":
        return baseline.evaluate(point, fixture, threshold)
    points = _points([point])
    if np.any((points < 0) | (points > 1)):
        raise ValueError("analytic scoring expects a point in the unit cube")
    distance = float(scores(points, fixture)[0])
    return {"status": "complete", "distance": distance,
            "distance_lower_bound": distance, "simulator_dose_calls": 0}


def membership(points, fixture):
    """Eight disjoint latent octants; zero coordinates take the positive side."""
    if fixture != "shifted_rotated_box":
        return baseline.membership(points, fixture)
    points = _points(points)
    latent = _latent(points)
    accepted = ((points >= 0) & (points <= 1)).all(axis=1) & (scores(points, fixture) <= EPSILON)
    signs = (latent[:, :3] >= 0).astype(int)
    region = signs @ np.array([4, 2, 1])
    return {f"octant_{index}": accepted & (region == index) for index in range(8)}


def features(points, fixture):
    """Moment and CDF functions bounded in [0, 1] on the accepted target.

    The extension to points outside the target is not used for target estimates
    and need not be bounded. CDF indicators include their named cut point.
    """
    if fixture != "shifted_rotated_box":
        return baseline.features(points, fixture)
    points = _points(points)
    normalized = _latent(points) / SHIFTED_HALF_WIDTHS
    latent = (normalized + 1) / 2
    values = {**{f"mean_x{j + 1}": points[:, j] for j in range(DIMENSION)},
              "second_x1": points[:, 0] ** 2, "second_x2": points[:, 1] ** 2,
              "cross_x1_x2": points[:, 0] * points[:, 1],
              "cross_x2_x3": points[:, 1] * points[:, 2]}
    values.update({f"latent_mean_{j + 1}": latent[:, j] for j in range(DIMENSION)})
    values.update({f"latent_second_{j + 1}": latent[:, j] ** 2 for j in range(DIMENSION)})
    values.update({f"latent_cross_{j + 1}_{k + 1}": latent[:, j] * latent[:, k]
                   for j in range(DIMENSION) for k in range(j + 1, DIMENSION)})
    values.update({f"latent_cdf_{j + 1}_{name}": (latent[:, j] <= threshold).astype(float)
                   for j in range(DIMENSION) for name, threshold in CDF_THRESHOLDS})
    radial = np.max(np.abs(normalized), axis=1) ** DIMENSION
    values.update(radial_cdf=radial, radial_cdf_second=radial ** 2)
    return values


def truth(fixture):
    """Exact target mass, octant probabilities, and bounded-feature expectations."""
    if fixture != "shifted_rotated_box":
        return baseline.truth(fixture)
    # Keep serialized analytic truth independent of BLAS reduction order.
    # For x = center + Q*y, Cov(x) = Q*diag(a**2 / 3)*Q.T; simplify the
    # requested entries with exact rationals before rounding each once.
    center = Fraction(9, 50)
    a1, a2, a3 = Fraction(1, 4), Fraction(1, 48), Fraction(1, 32)
    variance = a1 ** 2 / 9 + a2 ** 2 / 6 + a3 ** 2 / 18
    cross12 = a1 ** 2 / 9 - a2 ** 2 / 6 + a3 ** 2 / 18
    cross23 = a1 ** 2 / 9 - a3 ** 2 / 9
    moments = {f"mean_x{j + 1}": float(SHIFTED_CENTER[j]) for j in range(DIMENSION)}
    moments.update(second_x1=float(center ** 2 + variance),
                   second_x2=float(center ** 2 + variance),
                   cross_x1_x2=float(center ** 2 + cross12),
                   cross_x2_x3=float(center ** 2 + cross23))
    moments.update({f"latent_mean_{j + 1}": 0.5 for j in range(DIMENSION)})
    moments.update({f"latent_second_{j + 1}": 1 / 3 for j in range(DIMENSION)})
    moments.update({f"latent_cross_{j + 1}_{k + 1}": 0.25
                    for j in range(DIMENSION) for k in range(j + 1, DIMENSION)})
    moments.update({f"latent_cdf_{j + 1}_{name}": threshold
                    for j in range(DIMENSION) for name, threshold in CDF_THRESHOLDS})
    moments.update(radial_cdf=0.5, radial_cdf_second=1 / 3)
    return {"mass": 1 / 12288,
            "region_masses": {f"octant_{index}": 1 / 8 for index in range(8)},
            "moments": moments}


def oracle_retained_fraction(fixture, restricted=False):
    """Restricted controls retain the positive first-latent-coordinate half."""
    if fixture != "shifted_rotated_box":
        return baseline.oracle_retained_fraction(fixture, restricted)
    if type(restricted) is not bool:
        raise ValueError("restricted must be a boolean")
    return 0.5 if restricted else 1.0


def oracle_sample(rng, n, fixture, restricted=False):
    """Draw uniform latent coordinates and return their actual constant log q.

    The restriction retains y1 >= 0, hence octants 4 through 7 and half the
    target volume. Full and restricted calls consume the same number of draws.
    """
    if fixture != "shifted_rotated_box":
        return baseline.oracle_sample(rng, n, fixture, restricted)
    retained = oracle_retained_fraction(fixture, restricted)
    if type(n) is not int or n < 0:
        raise ValueError("oracle count must be a nonnegative integer")
    if not isinstance(rng, np.random.Generator):
        raise ValueError("oracle randomness must be a numpy Generator")
    latent = (2 * rng.random((n, DIMENSION)) - 1) * SHIFTED_HALF_WIDTHS
    if restricted:
        latent[:, 0] = (latent[:, 0] + SHIFTED_HALF_WIDTHS[0]) / 2
    points = SHIFTED_CENTER + latent @ SHIFTED_BASIS.T
    log_q = np.full(n, -math.log(truth(fixture)["mass"] * retained))
    return points, log_q
