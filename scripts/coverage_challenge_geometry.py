"""Analytic unseen-geometry challenges for the frozen local-move proposal.

These seven-dimensional targets, moments and exact samplers are specified
without consulting learned-proposal outcomes. All densities use Lebesgue
measure on the unit cube. Restricted oracles deliberately remove target support;
they are truth-screen controls, not defensive importance proposals.
"""

import math

import numpy as np

DIMENSION = 7
EPSILON = 1.0
FIXTURES = ("rotated_box", "annular_cylinder", "unequal_balls")

ROTATED_HALF_WIDTHS = np.array([3 / 8, 1 / 32] + [1 / 6] * 5)
ANNULUS_INNER_RADIUS = 1 / 4
ANNULUS_OUTER_RADIUS = 3 / 8
ANNULUS_HALF_WIDTH = 1 / 8
BALL_CENTERS = np.array([[1 / 4] * DIMENSION, [3 / 4] * DIMENSION])
BALL_RADII = np.array([1 / 5, 1 / 4])
BALL_A_FRACTION = 4**7 / (4**7 + 5**7)


def _fixture(fixture):
    if fixture not in FIXTURES:
        raise ValueError(f"unknown coverage fixture: {fixture}")


def _points(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != DIMENSION or not np.isfinite(points).all():
        raise ValueError("coverage points must be a finite (n, 7) array")
    return points


def _rotated(points):
    values = points - 0.5
    values[:, 0] = (points[:, 0] + points[:, 1] - 1) / math.sqrt(2)
    values[:, 1] = (points[:, 0] - points[:, 1]) / math.sqrt(2)
    return values


def _ball_coordinates(points):
    coordinates = (points[:, None, :] - BALL_CENTERS) / BALL_RADII[None, :, None]
    radii = np.linalg.norm(coordinates, axis=2)
    # The balls are disjoint. Outside the target this tie rule only provides a
    # deterministic extension of features that are assessed on accepted points.
    components = np.argmin(radii, axis=1)
    return coordinates[np.arange(len(points)), components], radii, components


def scores(points, fixture):
    """Complete vector scores; acceptance is exactly score <= 1 in the cube."""
    _fixture(fixture)
    points = _points(points)
    if fixture == "rotated_box":
        return np.max(np.abs(_rotated(points)) / ROTATED_HALF_WIDTHS, axis=1)
    if fixture == "annular_cylinder":
        radius = np.linalg.norm(points[:, :2] - 0.5, axis=1)
        radial = np.abs(radius - 5 / 16) / (1 / 16)
        slabs = np.max(np.abs(points[:, 2:] - 0.5) / ANNULUS_HALF_WIDTH, axis=1)
        return np.maximum(radial, slabs)
    return np.min(_ball_coordinates(points)[1], axis=1)


def evaluate(point, fixture, threshold=None):
    """Score a cube point fully; the caller's threshold never changes target."""
    points = _points([point])
    if np.any((points < 0) | (points > 1)):
        raise ValueError("analytic scoring expects a point in the unit cube")
    distance = float(scores(points, fixture)[0])
    return {"status": "complete", "distance": distance,
            "distance_lower_bound": distance, "simulator_dose_calls": 0}


def membership(points, fixture):
    """Disjoint declared regions, with boundaries assigned deterministically."""
    points = _points(points)
    accepted = ((points >= 0) & (points <= 1)).all(axis=1) & (scores(points, fixture) <= EPSILON)
    if fixture == "rotated_box":
        latent = _rotated(points)
        region = 2 * (latent[:, 0] >= 0).astype(int) + (latent[:, 1] >= 0)
        return {f"quadrant_{index}": accepted & (region == index) for index in range(4)}
    if fixture == "annular_cylinder":
        angle = np.mod(np.arctan2(points[:, 1] - 0.5, points[:, 0] - 0.5), 2 * np.pi)
        # A tiny negative angle can round to exactly 2*pi under remainder.
        # It still belongs to the last octant, never an unassigned ninth region.
        sector = np.minimum(7, np.floor(angle / (np.pi / 4)).astype(int))
        return {f"sector_{index}": accepted & (sector == index) for index in range(8)}
    components = _ball_coordinates(points)[2]
    return {name: accepted & (components == index) for index, name in enumerate(("A", "B"))}


def features(points, fixture):
    """Moment functions bounded in [0, 1] on the accepted target.

    Values outside that target are not estimates and need not remain bounded.
    A common set of raw coordinate moments accompanies geometry-specific moments.
    """
    _fixture(fixture)
    points = _points(points)
    values = {**{f"mean_x{j + 1}": points[:, j] for j in range(DIMENSION)},
              "second_x1": points[:, 0] ** 2, "second_x2": points[:, 1] ** 2,
              "cross_x1_x2": points[:, 0] * points[:, 1],
              "cross_x2_x3": points[:, 1] * points[:, 2]}
    if fixture == "rotated_box":
        latent = (_rotated(points) / ROTATED_HALF_WIDTHS + 1) / 2
        values.update({f"latent_mean_{j + 1}": latent[:, j] for j in range(DIMENSION)})
        values.update(latent_second_1=latent[:, 0] ** 2, latent_second_2=latent[:, 1] ** 2,
                      latent_cross_1_2=latent[:, 0] * latent[:, 1])
    elif fixture == "annular_cylinder":
        centered = points[:, :2] - 0.5
        radial = (np.sum(centered**2, axis=1) - 1 / 16) / (5 / 64)
        angle = np.arctan2(centered[:, 1], centered[:, 0])
        values.update(radial_cdf=radial, radial_cdf_second=radial**2)
        for frequency in (1, 2, 4):
            values[f"angular_cos_{frequency}"] = (1 + np.cos(frequency * angle)) / 2
            values[f"angular_sin_{frequency}"] = (1 + np.sin(frequency * angle)) / 2
    else:
        local, _, _ = _ball_coordinates(points)
        radial = np.linalg.norm(local, axis=1) ** DIMENSION
        values.update(radial_cdf=radial, radial_cdf_second=radial**2)
        values.update({f"local_mean_{j + 1}": (local[:, j] + 1) / 2 for j in range(DIMENSION)})
        values.update({f"local_second_{j + 1}": local[:, j] ** 2 for j in range(DIMENSION)})
    return values


def truth(fixture):
    """Exact masses and moments of the uniform distribution on each target."""
    _fixture(fixture)
    moments = {**{f"mean_x{j + 1}": 0.5 for j in range(DIMENSION)},
               "second_x1": 0.0, "second_x2": 0.0,
               "cross_x1_x2": 0.25, "cross_x2_x3": 0.25}
    if fixture == "rotated_box":
        mass = 1 / 5184
        regions = {f"quadrant_{index}": 0.25 for index in range(4)}
        moments.update(second_x1=1681 / 6144, second_x2=1681 / 6144, cross_x1_x2=1679 / 6144)
        moments.update({f"latent_mean_{j + 1}": 0.5 for j in range(DIMENSION)})
        moments.update(latent_second_1=1 / 3, latent_second_2=1 / 3, latent_cross_1_2=0.25)
    elif fixture == "annular_cylinder":
        mass = 5 * math.pi / 65536
        regions = {f"sector_{index}": 1 / 8 for index in range(8)}
        moments.update(second_x1=77 / 256, second_x2=77 / 256,
                       radial_cdf=0.5, radial_cdf_second=1 / 3)
        for frequency in (1, 2, 4):
            moments[f"angular_cos_{frequency}"] = 0.5
            moments[f"angular_sin_{frequency}"] = 0.5
    else:
        mass = 16 * math.pi**3 / 105 * (1 / 5**7 + 1 / 4**7)
        p = BALL_A_FRACTION
        regions = {"A": p, "B": 1 - p}
        moments.update({f"mean_x{j + 1}": 3 / 4 - p / 2 for j in range(DIMENSION)})
        second = p * (1 / 16 + 1 / 225) + (1 - p) * (9 / 16 + 1 / 144)
        moments.update(second_x1=second, second_x2=second,
                       cross_x1_x2=9 / 16 - p / 2, cross_x2_x3=9 / 16 - p / 2,
                       radial_cdf=0.5, radial_cdf_second=1 / 3)
        moments.update({f"local_mean_{j + 1}": 0.5 for j in range(DIMENSION)})
        moments.update({f"local_second_{j + 1}": 1 / 9 for j in range(DIMENSION)})
    return {"mass": mass, "region_masses": regions, "moments": moments}


def oracle_retained_fraction(fixture, restricted=False):
    """Support retained by the exact positive or deliberately restricted oracle."""
    _fixture(fixture)
    if type(restricted) is not bool:
        raise ValueError("restricted must be a boolean")
    if not restricted:
        return 1.0
    return 1 - BALL_A_FRACTION if fixture == "unequal_balls" else 0.5


def oracle_sample(rng, n, fixture, restricted=False):
    """Draw exact uniform target points and their actual (possibly restricted) q.

    Positive oracles retain the complete target. Restricted controls keep the
    rotated positive long-axis half, the annular upper semicircle, or ball B.
    Neither kind uses the learned proposal or measures its ability to discover
    target regions. The draw order is fixed by this implementation.
    """
    retained = oracle_retained_fraction(fixture, restricted)
    if type(n) is not int or n < 0:
        raise ValueError("oracle count must be a nonnegative integer")
    if not isinstance(rng, np.random.Generator):
        raise ValueError("oracle randomness must be a numpy Generator")
    if fixture == "rotated_box":
        latent = (2 * rng.random((n, DIMENSION)) - 1) * ROTATED_HALF_WIDTHS
        if restricted:
            latent[:, 0] = (latent[:, 0] + ROTATED_HALF_WIDTHS[0]) / 2
        points = latent + 0.5
        points[:, 0] = 0.5 + (latent[:, 0] + latent[:, 1]) / math.sqrt(2)
        points[:, 1] = 0.5 + (latent[:, 0] - latent[:, 1]) / math.sqrt(2)
    elif fixture == "annular_cylinder":
        uniforms = rng.random((n, DIMENSION))
        radius = np.sqrt(1 / 16 + uniforms[:, 0] * (5 / 64))
        angle = uniforms[:, 1] * (math.pi if restricted else 2 * math.pi)
        points = 3 / 8 + uniforms / 4
        points[:, 0] = 0.5 + radius * np.cos(angle)
        points[:, 1] = 0.5 + radius * np.sin(angle)
    else:
        # Consume the component draw even for B-only controls, retaining a
        # straightforward reproducible stream structure across both controls.
        component_uniforms = rng.random(n)
        components = (component_uniforms >= BALL_A_FRACTION).astype(int)
        if restricted:
            components[:] = 1
        directions = rng.normal(size=(n, DIMENSION))
        lengths = np.linalg.norm(directions, axis=1)
        # A zero normal vector is a probability-zero event in exact arithmetic;
        # fail explicitly instead of silently changing the claimed distribution.
        if np.any(lengths == 0):
            raise ValueError("oracle ball direction has zero norm")
        radii = rng.random(n) ** (1 / DIMENSION)
        directions = directions / lengths[:, None]
        points = BALL_CENTERS[components] + (BALL_RADII[components] * radii)[:, None] * directions
    log_q = np.full(n, -math.log(truth(fixture)["mass"] * retained))
    return points, log_q
