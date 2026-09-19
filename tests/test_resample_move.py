"""Scientific invariants of proposal-training resample/move pilots."""

from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import resample_move as rm  # noqa: E402
from bounded_proposal import BoundedGaussianMixture  # noqa: E402


def plan(**changes):
    return {**rm.PILOT_PLAN, "islands": 2, "particles_per_island": 16,
            "levels": 3, "moves_per_level": 2, "kernel_neighbors": 4, **changes}


class FixedProposal:
    def __init__(self, proposed, old_density=1, new_density=4):
        self.proposed = np.asarray(proposed, dtype=float)
        self.old_density, self.new_density = old_density, new_density
        self.density_calls = 0

    def sample(self, rng, n):
        assert n == len(self.proposed)
        return self.proposed.copy(), np.ones(n, dtype=int)

    def log_density(self, points):
        density = self.old_density if self.density_calls == 0 else self.new_density
        self.density_calls += 1
        return np.full(len(points), np.log(density))


class FixedUniform:
    def __init__(self, value):
        self.value = value

    def random(self, n):
        return np.full(n, self.value)


@pytest.mark.parametrize("old,new,u,moves,evaluations", [
    (1, 4, 0.5, 0, 0),  # Reverse ratio would incorrectly accept this proposal.
    (4, 1, 0.9, 1, 1),
    (1, 4, 0.25, 1, 1),  # Inclusive exact MH boundary.
    (1, 4, 0, 1, 1),
])
def test_independent_mh_uses_old_over_new_density_and_skips_failed_tests(
        old, new, u, moves, evaluations):
    calls = []

    def evaluate(point, threshold):
        calls.append((point.copy(), threshold))
        return {"distance": 0.05, "simulator_dose_calls": 7}

    x, distances, ancestors, stats = rm.independent_move(
        [[0.2]], [0.1], [11], FixedProposal([[0.8]], old, new),
        0.2, FixedUniform(u), evaluate)
    assert x.tolist() == [[0.8 if moves else 0.2]]
    assert distances.tolist() == [0.05 if moves else 0.1]
    assert ancestors.tolist() == [11]
    assert len(calls) == evaluations
    assert stats == {"attempts": 1, "evaluated": evaluations,
                     "mh_rejected_without_evaluation": 1 - evaluations,
                     "early_rejected": 0, "accepted_moves": moves,
                     "simulator_dose_calls": 7 * evaluations}


@pytest.mark.parametrize("distance,expected_moves", [(0.2, 1), (np.nextafter(0.2, 1), 0)])
def test_mh_pass_is_still_subject_to_exact_conditioning_threshold(distance, expected_moves):
    x, scores, _, stats = rm.independent_move(
        [[0.2]], [0.1], [0], FixedProposal([[0.8]], 1, 1), 0.2,
        FixedUniform(0.5), lambda p, threshold: {"distance": distance, "simulator_dose_calls": 2})
    assert stats["accepted_moves"] == expected_moves
    assert x[0, 0] == (0.8 if expected_moves else 0.2)
    assert scores[0] == (distance if expected_moves else 0.1)


def test_strict_early_rejection_keeps_old_state_and_charges_partial_cost():
    before = np.array([[0.2]])
    x, scores, _, stats = rm.independent_move(
        before, [0.1], [0], FixedProposal([[0.8]], 1, 1), 0.2,
        FixedUniform(0.5), lambda p, threshold: {
            "distance": None, "distance_lower_bound": np.nextafter(threshold, 1),
            "simulator_dose_calls": 6})
    assert np.array_equal(x, before) and before.tolist() == [[0.2]]
    assert scores.tolist() == [0.1]
    assert stats["early_rejected"] == 1 and stats["accepted_moves"] == 0
    assert stats["simulator_dose_calls"] == 6


def test_fixed_proposal_moves_preserve_known_separated_conditional_uniform_target():
    # Exact target: Uniform([0,.2] union [.65,.95]), with masses .4 and .6.
    rng = np.random.default_rng(5294)
    u = rng.uniform(0, 0.5, 18000)
    points = np.where(u < 0.2, u, u + 0.45)[:, None]
    scores = np.zeros(len(points))
    ancestors = np.arange(len(points))
    proposal = BoundedGaussianMixture(1, means=[[0.08]], scales=[[0.06]])

    def evaluate(p, threshold):
        inside = 0 <= p[0] <= 0.2 or 0.65 <= p[0] <= 0.95
        return {"distance": 0.0 if inside else 1.0, "simulator_dose_calls": 1}

    for _ in range(3):
        points, scores, ancestors, _ = rm.independent_move(
            points, scores, ancestors, proposal, 0.5, rng, evaluate)
    assert np.all((points[:, 0] <= 0.2) | ((points[:, 0] >= 0.65) & (points[:, 0] <= 0.95)))
    assert np.all(scores == 0)
    assert np.mean(points[:, 0] < 0.2) == pytest.approx(0.4, abs=0.015)
    assert points[:, 0].mean() == pytest.approx(0.52, abs=0.01)
    # Independent proposal draws alone heavily favor the left mode: the target
    # invariance above depends on MH, not an accidentally well-matched proposal.
    raw, _ = proposal.sample(np.random.default_rng(192), 10000)
    raw_in_target = raw[(raw[:, 0] <= 0.2) | ((raw[:, 0] >= 0.65) & (raw[:, 0] <= 0.95))]
    assert raw_in_target[:, 0].mean() < 0.2


def test_local_widths_ignore_duplicate_neighbor_votes_but_keep_mixture_multiplicity():
    settings = plan(kernel_neighbors=2)
    points = np.array([[0.0], [0.2], [0.8], [1.0]])
    duplicated = np.concatenate([points, np.repeat([[0.2]], 100, axis=0)])
    expected = np.full((4, 1), 0.15)  # Each nearest pair is .2 apart: 1.5*.1.
    assert rm.local_scales(points, settings) == pytest.approx(expected)
    fitted = rm.fit_proposal(duplicated, settings)
    assert len(fitted.means) == 104
    assert fitted.scales == pytest.approx(np.full((104, 1), 0.15))


def test_local_width_floor_and_ceiling_are_explicit():
    identical = np.array([[0.2, 0.7]] * 12)
    assert rm.local_scales(identical, plan()) == pytest.approx(np.full((12, 2), 0.005))
    settings = plan(kernel_neighbors=4, kernel_scale=10, kernel_ceiling=0.2)
    assert rm.local_scales([[0], [0.2], [0.8], [1]], settings) == pytest.approx(np.full((4, 1), 0.2))


def _score(point, threshold):
    distance = float(point[0])
    if threshold is not None and distance > threshold:
        return {"distance": None, "distance_lower_bound": distance, "simulator_dose_calls": 2}
    return {"distance": distance, "simulator_dose_calls": 5}


def test_pilot_thresholds_are_monotone_and_never_cross_final_epsilon():
    settings = plan(levels=8, particles_per_island=32)
    result = rm.train_pilot(_score, 2, 0.1, np.random.SeedSequence(126), settings)
    for island in result["islands"]:
        thresholds = np.array([row["threshold"] for row in island["history"]])
        assert np.all(thresholds >= 0.1) and np.all(np.diff(thresholds) <= 0)
        assert np.all(np.asarray(island["final_distances"]) <= thresholds[-1])
        if island["first_epsilon_level"] is not None:
            assert np.all(thresholds[island["first_epsilon_level"] - 1:] == 0.1)


def test_ties_at_the_threshold_are_all_retained():
    # Every initial score is a tie above epsilon. Choosing an arbitrary half
    # would unnecessarily destroy ancestry and misrepresent the conditioning set.
    result = rm.train_pilot(lambda p, threshold: {"distance": 0.3, "simulator_dose_calls": 1},
                            1, 0.1, np.random.SeedSequence(3), plan())
    for island in result["islands"]:
        assert island["first_epsilon_level"] is None
        for h in island["history"]:
            assert h["threshold"] == 0.3
            assert h["retained_before_resampling"] == 16
            assert not h["resampled"] and h["unique_ancestors"] == 16
    assert not result["all_islands_reached_epsilon"]


def test_reaching_epsilon_does_not_trigger_unnecessary_resampling_or_early_stopping():
    settings = plan(levels=5)
    result = rm.train_pilot(lambda p, threshold: {"distance": 0.05, "simulator_dose_calls": 1},
                            1, 0.1, np.random.SeedSequence(97), settings)
    for island in result["islands"]:
        assert island["first_epsilon_level"] == 1
        assert len(island["history"]) == 5
        assert island["final_ancestors"] == list(range(16))
        assert all(h["threshold"] == 0.1 and not h["resampled"] for h in island["history"])
    assert result["all_islands_reached_epsilon"]
    assert result["attempts"] == 2 * 16 * (1 + 5 * 2)


def test_costs_count_evaluations_and_independent_islands_replay_from_their_streams():
    settings = plan()
    calls = []

    def evaluate(point, threshold):
        calls.append((point.copy(), threshold))
        return {"distance": float(point[0]), "simulator_dose_calls": 5}

    result = rm.train_pilot(evaluate, 2, 0.2, np.random.SeedSequence(293), settings)
    streams = np.random.SeedSequence(293).spawn(2)
    n_evaluated = 2 * 16
    for island, stream in zip(result["islands"], streams):
        expected = np.random.default_rng(stream).uniform(size=(16, 2))
        assert np.array_equal(expected, island["initial_points"])
        assert island["spawn_key"] == list(stream.spawn_key)
        assert island["initial_simulator_dose_calls"] == 16 * 5
        for h in island["history"]:
            assert h["attempts"] == 16 * 2
            assert h["evaluated"] + h["mh_rejected_without_evaluation"] == h["attempts"]
            assert h["accepted_moves"] <= h["evaluated"]
            assert h["simulator_dose_calls"] == h["evaluated"] * 5
            n_evaluated += h["evaluated"]
    assert result["islands"][0]["initial_points"] != result["islands"][1]["initial_points"]
    assert result["attempts"] == 2 * 16 * (1 + 3 * 2)
    assert len(calls) == n_evaluated and result["simulator_dose_calls"] == n_evaluated * 5
    assert len(result["final_points"]) == 32
    replay = rm.train_pilot(evaluate, 2, 0.2, np.random.SeedSequence(293), settings)
    assert replay == result


def test_each_level_uses_one_frozen_proposal_across_all_moves(monkeypatch):
    real_fit, real_move = rm.fit_proposal, rm.independent_move
    fitted, used = [], []

    def fit(points, plan):
        result = real_fit(points, plan)
        fitted.append(result)
        return result

    def move(points, distances, ancestors, proposal, *args):
        used.append(proposal)
        return real_move(points, distances, ancestors, proposal, *args)

    monkeypatch.setattr(rm, "fit_proposal", fit)
    monkeypatch.setattr(rm, "independent_move", move)
    settings = plan(moves_per_level=3)
    rm.train_pilot(_score, 1, 0.2, np.random.SeedSequence(421), settings)
    assert len(fitted) == 2 * 3
    assert len(used) == len(fitted) * 3
    for i, proposal in enumerate(fitted):
        assert all(p is proposal for p in used[i * 3:(i + 1) * 3])


@pytest.mark.parametrize("distance,bound,threshold", [
    (None, 0.2, 0.2), (None, 0.1, 0.2), (None, np.inf, 0.2),
    (None, np.nan, 0.2), (None, None, 0.2), (None, 0.3, None),
    (-1, None, 0.2), (np.nan, None, 0.2), (np.inf, None, 0.2),
])
def test_invalid_scores_or_unproved_early_rejections_are_rejected(distance, bound, threshold):
    with pytest.raises(ValueError):
        rm._evaluation({"distance": distance, "distance_lower_bound": bound,
                        "simulator_dose_calls": 1}, threshold)


@pytest.mark.parametrize("cost", [-1, 0.5, True, np.inf])
def test_evaluation_cost_must_be_an_actual_nonnegative_count(cost):
    with pytest.raises(ValueError):
        rm._evaluation({"distance": 0.1, "simulator_dose_calls": cost}, 0.2)


@pytest.mark.parametrize("change", [
    {"islands": 0}, {"particles_per_island": 1}, {"levels": 1.1},
    {"moves_per_level": True}, {"kernel_neighbors": 0}, {"retained_fraction": 0},
    {"retained_fraction": 1}, {"uniform_weight": np.nan}, {"kernel_scale": -1},
    {"kernel_floor": 0.6, "kernel_ceiling": 0.5},
])
def test_invalid_plan_is_rejected_before_evaluation(change):
    def forbidden(*args):
        pytest.fail("invalid plan reached the simulator")

    with pytest.raises(ValueError):
        rm.train_pilot(forbidden, 1, 0.1, np.random.SeedSequence(7), plan(**change))


@pytest.mark.parametrize("points", [[[np.nan]], [[-0.1]], [[1.1]], [], [0.2]])
def test_invalid_local_kernel_centers_are_rejected(points):
    with pytest.raises(ValueError):
        rm.local_scales(points, plan())
    with pytest.raises(ValueError):
        rm.fit_proposal(points, plan())


@pytest.mark.parametrize("points,scores,ancestors,threshold", [
    ([[0.2]], [0.3], [0], 0.2), ([[0.2]], [-0.1], [0], 0.2),
    ([[0.2]], [np.nan], [0], 0.2), ([[0.2]], [np.inf], [0], 0.2),
    ([[0.2]], [], [0], 0.2), ([[0.2]], [0.1], [], 0.2),
    ([[0.2]], [0.1], [0], np.nan), ([[0.2]], [0], [0], 0),
    ([0.2], [0.1], [0], 0.2),
])
def test_invalid_starting_states_cannot_enter_an_mh_sweep(points, scores, ancestors, threshold):
    def forbidden(*args):
        pytest.fail("invalid starting state reached a target evaluation")

    with pytest.raises(ValueError):
        rm.independent_move(points, scores, ancestors, FixedProposal([[0.8]]),
                            threshold, FixedUniform(0.5), forbidden)
