"""Target preservation and accounting for local proposal-training moves.

All targets are analytic. Stationarity checks start from exact target draws;
they do not treat pilot particles or clones as independent posterior samples.
"""

from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.stats import truncnorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import resample_move as historical  # noqa: E402
import resample_move_local as local  # noqa: E402
from bounded_proposal import BoundedGaussianMixture  # noqa: E402


def plan(**changes):
    return {**local.PILOT_PLAN, "islands": 2, "particles_per_island": 16,
            "levels": 3, "moves_per_level": 3, "kernel_neighbors": 4, **changes}


class ControlledRng:
    """Control component-choice and MH coins independently of fixed proposals."""

    def __init__(self, local_coins, mh_coins, normals=None):
        self.coins = iter([np.asarray(local_coins), np.asarray(mh_coins)])
        self.normals = None if normals is None else np.asarray(normals)

    def random(self, n):
        values = next(self.coins)
        assert values.shape == (n,)
        return values.copy()

    def normal(self, size):
        if size[0] == 0:
            return np.empty(size)
        assert self.normals is not None and self.normals.shape == size
        return self.normals.copy()


class FixedProposal:
    def __init__(self, points, density=None):
        self.points = np.asarray(points, dtype=float)
        self.density = density

    def sample(self, rng, n):
        assert len(self.points) == n
        # Component identity must not replace the full mixture density in MH.
        return self.points.copy(), np.zeros(n, dtype=int)

    def log_density(self, points):
        assert self.density is not None, "local moves must not use independence-MH densities"
        return self.density(points)


def test_reflection_handles_multiple_crossings_without_boundary_atoms():
    values = np.array([-6.25, -4, -3.2, -2, -1.75, -1, -.2,
                       0, .2, 1, 1.2, 2, 3.8, 4, 6.25])
    expected = [.25, 0, .8, 0, .25, 1, .2, 0, .2, 1, .8, 0, .2, 0, .25]
    assert local.reflect_cube(values) == pytest.approx(expected, abs=1e-14)
    assert local.reflect_cube(-values) == pytest.approx(expected, abs=1e-14)
    assert local.reflect_cube(values + 12) == pytest.approx(expected, abs=1e-14)
    assert local.reflect_cube(values.reshape(3, 5)).shape == (3, 5)


@pytest.mark.parametrize("fixture", ["cube", "slab", "disconnected"])
@pytest.mark.parametrize("local_probability", [1.0, 0.5])
def test_reflected_and_mixed_moves_preserve_exact_conditional_uniform_targets(
        fixture, local_probability):
    rng = np.random.default_rng(9206)
    n = 20000
    points = rng.uniform(size=(n, 2))
    if fixture == "slab":
        points[:, 0] *= .125
        mean, second, left_mass = .0625, 1 / 192, None
    elif fixture == "disconnected":
        u = points[:, 0] * .5
        points[:, 0] = np.where(u < .2, u, u + .45)
        mean, second, left_mass = .52, .4 * .2**2 / 3 + .6 * (.65**2 + .65 * .95 + .95**2) / 3, .4
    else:
        mean, second, left_mass = .5, 1 / 3, None

    def inside(p):
        if fixture == "slab":
            return p[0] <= .125
        if fixture == "disconnected":
            return p[0] <= .2 or .65 <= p[0] <= .95
        return True

    def evaluate(p, threshold):
        return {"distance": 0.0 if inside(p) else 1.0, "simulator_dose_calls": 1}

    # Deliberately far from the target's uniform distribution, especially in y.
    proposal = BoundedGaussianMixture(2, means=[[.04, .15]], scales=[[.04, .07]])
    scores, ancestors = np.zeros(n), np.arange(n)
    for _ in range(3):
        points, scores, ancestry, stats = local.mixed_move(
            points, scores, ancestors, proposal, .5, rng, evaluate,
            step=np.array([.35, .8]), local_probability=local_probability)
        assert np.array_equal(ancestry, ancestors)
        assert stats["local_attempts"] + stats["global_attempts"] == n
    assert np.all((points > 0) & (points < 1))  # Clipping would create atoms.
    assert all(inside(p) for p in points)
    assert np.all(scores == 0)
    assert points[:, 0].mean() == pytest.approx(mean, abs=.0015 if fixture == "slab" else .008)
    assert np.mean(points[:, 0]**2) == pytest.approx(second, abs=.0003 if fixture == "slab" else .01)
    assert points[:, 1].mean() == pytest.approx(.5, abs=.008)
    assert np.mean(points[:, 1]**2) == pytest.approx(1 / 3, abs=.01)
    assert np.mean(points[:, 0] * points[:, 1]) == pytest.approx(mean / 2, abs=.0015 if fixture == "slab" else .008)
    if left_mass is not None:
        assert np.mean(points[:, 0] <= .2) == pytest.approx(left_mass, abs=.012)


@pytest.mark.parametrize("coin_position,accepted", [(0.5, True), (1.5, False)])
def test_global_move_uses_full_mixture_density_not_chosen_uniform_component(coin_position, accepted):
    means, scales = np.array([.15, .85]), np.array([.07, .15])
    proposal = BoundedGaussianMixture(1, means=means[:, None], scales=scales[:, None])
    old, new = .4, .14

    def reference_density(x):
        return .2 + .4 * truncnorm.pdf(x, -means / scales, (1 - means) / scales,
                                       loc=means, scale=scales).sum()

    ratio = reference_density(old) / reference_density(new)
    assert 0 < ratio < .5
    coin = coin_position * ratio
    rng = ControlledRng([.9], [coin])
    calls = []

    def evaluate(p, threshold):
        calls.append(p.copy())
        return {"distance": .2, "simulator_dose_calls": 7}

    points, scores, ancestors, stats = local.mixed_move(
        [[old]], [.1], [17], FixedProposal([[new]], proposal.log_density), .2,
        rng, evaluate, step=[.04], local_probability=0)
    assert points[0, 0] == (new if accepted else old)
    assert scores[0] == (.2 if accepted else .1)
    assert ancestors.tolist() == [17]
    assert len(calls) == int(accepted)
    assert stats["simulator_dose_calls"] == 7 * int(accepted)
    assert stats["mh_rejected_without_evaluation"] == int(not accepted)


def test_mixed_sweep_reconciles_local_global_rejection_and_exact_costs():
    before = np.full((5, 1), .2)
    old_scores = np.full(5, .1)
    ancestor_ids = np.arange(5)
    densities = iter([np.ones(3), np.array([4., 4., 1.])])
    proposal = FixedProposal([[.6], [.7], [.9], [.8], [.75]],
                             lambda p: np.log(next(densities)))
    rng = ControlledRng([.1, .9, .9, .1, .9], [.99, .5, .1, .99, .3], [[1.], [-1.]])
    calls = []

    def evaluate(p, threshold):
        calls.append(float(p[0]))
        if p[0] > .85:
            return {"distance": None, "distance_lower_bound": np.nextafter(threshold, 1),
                    "simulator_dose_calls": 2}
        return {"distance": .7 if p[0] > .7 else (.4 if p[0] > .2 else .5),
                "simulator_dose_calls": 5}

    points, scores, ancestors, stats = local.mixed_move(
        before, old_scores, ancestor_ids, proposal, .5, rng, evaluate,
        step=[.1], local_probability=.5)
    assert points[:, 0] == pytest.approx([.3, .2, .2, .1, .2])
    assert scores == pytest.approx([.4, .1, .1, .5, .1])
    assert calls == pytest.approx([.3, .9, .1, .75])
    assert stats == {"attempts": 5, "evaluated": 4, "mh_rejected_without_evaluation": 1,
                     "early_rejected": 1, "accepted_moves": 2, "simulator_dose_calls": 17,
                     "local_attempts": 2, "global_attempts": 3, "local_accepted_moves": 2}
    assert np.array_equal(ancestors, ancestor_ids)
    assert np.array_equal(before, np.full((5, 1), .2)) and np.all(old_scores == .1)


def test_local_clones_receive_distinct_perturbations_and_no_independence_ratio():
    n = 64
    points = np.full((n, 2), .5)
    ancestors = np.full(n, 27)
    out, scores, ancestry, stats = local.mixed_move(
        points, np.zeros(n), ancestors, FixedProposal(np.zeros((n, 2))), .5,
        np.random.default_rng(9711),
        lambda p, t: {"distance": 0.0, "simulator_dose_calls": 1},
        step=[.04, .06], local_probability=1)
    assert len(np.unique(out, axis=0)) == n
    assert np.all(out != points)
    assert np.array_equal(ancestors, ancestry) and np.all(scores == 0)
    assert stats["local_accepted_moves"] == n and stats["global_attempts"] == 0


@pytest.mark.parametrize("strict", [True, False])
def test_local_early_rejection_requires_a_strict_bound_and_retains_its_cost(strict):
    def evaluate(p, threshold):
        return {"distance": None,
                "distance_lower_bound": np.nextafter(threshold, 1) if strict else threshold,
                "simulator_dose_calls": 2}

    def run():
        return local.mixed_move(
            [[.2]], [.1], [7], FixedProposal([[.8]]), .2,
            ControlledRng([0.0], [.9], [[1.0]]), evaluate,
            step=[.04], local_probability=1)

    if not strict:
        with pytest.raises(ValueError, match="strict early-rejection"):
            run()
    else:
        points, scores, ancestors, stats = run()
        assert points.tolist() == [[.2]] and scores.tolist() == [.1]
        assert ancestors.tolist() == [7]
        assert stats["early_rejected"] == 1 and stats["simulator_dose_calls"] == 2
        assert stats["local_attempts"] == 1 and stats["local_accepted_moves"] == 0


def test_gap_neighborhood_protects_two_minority_centers_from_thirty_majority_neighbors():
    settings = plan(kernel_neighbors=32)
    unique = np.concatenate([[[.10], [.12]], np.linspace(.75, .95, 40)[:, None]])
    points = np.concatenate([unique, np.repeat([[.10], [.12]], 30, axis=0)])
    scales = local.local_scales(points, settings)
    assert scales[:2, 0] == pytest.approx([.015, .015])
    assert historical.local_scales(points, settings)[:2, 0].min() > .1
    assert local.local_scales(unique, settings) == pytest.approx(scales[:len(unique)])
    proposal = local.fit_proposal(points, settings)
    assert len(proposal.means) == len(points) and np.array_equal(proposal.means, points)
    assert proposal.scales == pytest.approx(scales)


@pytest.mark.parametrize("points,neighbors", [([[.2, .7]] * 20, 32), ([[.2, .7], [.8, .3]], 1)])
def test_singleton_neighborhood_uses_declared_finite_fallback(points, neighbors):
    settings = plan(kernel_neighbors=neighbors, singleton_scale=.04)
    assert local.local_scales(points, settings) == pytest.approx(np.full((len(points), 2), .04))


def test_gap_rule_does_not_claim_to_identify_an_isolated_singleton_mode():
    # With no second point in a putative minority mode there is no short
    # positive radius before its distant neighbors. The gap rule cannot know
    # that these observations lie in separate acceptance regions.
    points = np.concatenate([[[.1]], np.linspace(.75, .95, 40)[:, None]])
    settings = plan(kernel_neighbors=32)
    assert local.local_scales(points, settings)[0, 0] > .1


def test_shared_local_step_counts_duplicates_and_respects_declared_bounds():
    scales = np.array([[.02, .2], [.2, .02], [.2, .02], [.2, .02]])
    proposal = BoundedGaussianMixture(2, means=[[.2, .7]] * 4, scales=scales)
    assert local.local_step(proposal, plan()) == pytest.approx([.1, .01])
    settings = plan(local_scale_floor=.02, local_scale_ceiling=.08)
    assert local.local_step(proposal, settings) == pytest.approx([.08, .02])


def _score(point, threshold):
    distance = float(point[0])
    if threshold is not None and distance > threshold:
        return {"distance": None, "distance_lower_bound": distance, "simulator_dose_calls": 2}
    return {"distance": distance, "simulator_dose_calls": 5}


def test_each_level_freezes_one_proposal_and_one_shared_scale_for_every_sweep(monkeypatch):
    real_fit, real_step, real_move = local.fit_proposal, local.local_step, local.mixed_move
    fitted, steps, used = [], [], []

    def fit(points, settings):
        result = real_fit(points, settings)
        fitted.append(result)
        return result

    def step(proposal, settings):
        result = real_step(proposal, settings)
        steps.append(result)
        return result

    def move(points, scores, ancestors, proposal, threshold, rng, evaluate, **kwargs):
        used.append((proposal, kwargs["step"], kwargs["step"].copy()))
        return real_move(points, scores, ancestors, proposal, threshold, rng, evaluate, **kwargs)

    monkeypatch.setattr(local, "fit_proposal", fit)
    monkeypatch.setattr(local, "local_step", step)
    monkeypatch.setattr(local, "mixed_move", move)
    settings = plan()
    local.train_pilot(_score, 2, .2, np.random.SeedSequence(1927), settings)
    assert len(fitted) == len(steps) == settings["islands"] * settings["levels"]
    assert len(used) == len(fitted) * settings["moves_per_level"]
    for i, (proposal, scale) in enumerate(zip(fitted, steps)):
        assert scale.shape == (2,)
        for used_proposal, used_scale, observed in used[i * 3:(i + 1) * 3]:
            assert used_proposal is proposal and used_scale is scale
            assert np.array_equal(observed, scale)


def test_pilot_replay_independent_islands_and_sweep_costs_are_reconciled():
    settings = plan(levels=5)
    calls = []

    def evaluate(point, threshold):
        result = _score(point, threshold)
        calls.append(result["simulator_dose_calls"])
        return result

    seed = 48210
    result = local.train_pilot(evaluate, 2, .15, np.random.SeedSequence(seed), settings)
    streams = np.random.SeedSequence(seed).spawn(settings["islands"])
    n = settings["particles_per_island"]
    for island, stream in zip(result["islands"], streams):
        expected = np.random.default_rng(stream).uniform(size=(n, 2))
        assert np.array_equal(island["initial_points"], expected)
        assert island["spawn_key"] == list(stream.spawn_key)
        previous = max(.15, max(island["initial_distances"]))
        for h in island["history"]:
            assert .15 <= h["threshold"] <= previous
            previous = h["threshold"]
            assert h["attempts"] == n * settings["moves_per_level"]
            assert h["local_attempts"] + h["global_attempts"] == h["attempts"]
            assert h["evaluated"] + h["mh_rejected_without_evaluation"] == h["attempts"]
            assert h["local_accepted_moves"] <= h["local_attempts"]
            assert h["local_accepted_moves"] <= h["accepted_moves"] <= h["evaluated"] - h["early_rejected"]
            assert h["simulator_dose_calls"] == 2 * h["early_rejected"] + 5 * (h["evaluated"] - h["early_rejected"])
        assert np.all(np.array(island["final_distances"]) <= previous)
    assert result["attempts"] == 2 * n * (1 + settings["levels"] * settings["moves_per_level"])
    assert result["simulator_dose_calls"] == sum(calls)
    assert result["islands"][0]["initial_points"] != result["islands"][1]["initial_points"]
    assert local.train_pilot(_score, 2, .15, np.random.SeedSequence(seed), settings) == result


def test_archive_counter_hook_accepts_actual_pilot_accounting():
    pilot = local.train_pilot(_score, 2, .15, np.random.SeedSequence(32181), plan())
    local.validate_pilot_archive(pilot, fixture="analytic_test", seed=32181)


def _counter_archive(**changes):
    # A valid sweep: two local evaluations accepted; among three global
    # attempts, one skips evaluation, one rejects early, and one rejects fully.
    step = {"attempts": 5, "evaluated": 4, "mh_rejected_without_evaluation": 1,
            "early_rejected": 1, "accepted_moves": 2, "simulator_dose_calls": 17,
            "local_attempts": 2, "global_attempts": 3, "local_accepted_moves": 2,
            **changes}
    return {"plan": plan(), "islands": [{"history": [step]}]}


@pytest.mark.parametrize("changes", [
    {"global_attempts": 4},  # Branch attempts no longer sum to all attempts.
    {"local_attempts": 5, "global_attempts": 0},  # Every local move must be evaluated.
    {"local_attempts": 0, "global_attempts": 5, "local_accepted_moves": 1},
    {"local_accepted_moves": 3},  # More local successes than total successes.
    {"local_attempts": 4, "global_attempts": 1, "local_accepted_moves": 0},
    # Last case invents global successes, although all evaluated moves were local.
])
def test_archive_counter_hook_rejects_impossible_branch_accounting(changes):
    with pytest.raises(ValueError, match="do not reconcile"):
        local.validate_pilot_archive(_counter_archive(**changes))


@pytest.mark.parametrize("key", ["local_attempts", "global_attempts", "local_accepted_moves"])
@pytest.mark.parametrize("value", [-1, True, 1.5])
def test_archive_branch_counters_require_actual_nonnegative_integers(key, value):
    with pytest.raises(ValueError, match="nonnegative integers"):
        local.validate_pilot_archive(_counter_archive(**{key: value}))


@pytest.mark.parametrize("step", [[[.04], [.04]], [0.0], [-.1], [np.nan], [np.inf]])
def test_state_dependent_or_invalid_local_scales_are_rejected_before_scoring(step):
    def forbidden(*args):
        pytest.fail("invalid local scale reached a target evaluation")

    with pytest.raises(ValueError, match="shared positive finite diagonal"):
        local.mixed_move([[.2]], [.1], [0], FixedProposal([[.3]]), .2,
                         np.random.default_rng(1), forbidden, step=step, local_probability=.5)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_nonfinite_reflection_is_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        local.reflect_cube([bad])


@pytest.mark.parametrize("points", [[], [.2]])
def test_malformed_proposal_cloud_raises_value_error(points):
    with pytest.raises(ValueError, match="kernel centers"):
        local.fit_proposal(points, plan())
