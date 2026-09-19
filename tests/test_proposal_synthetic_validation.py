"""Known-truth checks without running a learned pilot or the biological model."""

import copy
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import proposal_synthetic_validation as validation  # noqa: E402
from bounded_proposal import BoundedGaussianMixture  # noqa: E402


def quadrature_points(fixture):
    """Degree-three product quadrature, expanded by exact volume proportions.

    These are deterministic integration nodes, never claimed to be independent
    production draws. They make the moment/error-screen unit checks exact.
    """
    nodes = np.array(list(itertools.product([-1 / np.sqrt(3), 1 / np.sqrt(3)], repeat=7)))
    if fixture == "separated_boxes":
        boxes = [(np.zeros(7), np.full(7, .25), 1),
                 (np.array([.25] + [.625] * 6), np.array([1.] + [.875] * 6), 3)]
    else:
        boxes = [(np.zeros(7), np.array([.125] + [1.] * 6), 4)]
    return np.concatenate([np.tile((lo + hi) / 2 + nodes * (hi - lo) / 2, (count, 1))
                           for lo, hi, count in boxes])


@pytest.mark.parametrize("fixture", validation.FIXTURES)
def test_scores_match_exact_regions_at_boundaries_and_random_points(fixture):
    if fixture == "separated_boxes":
        boundaries = np.array(list(itertools.product([0, .25], repeat=7)) +
                              list(itertools.product([.25, 1], *([[.625, .875]] * 6))))
    else:
        boundaries = np.array(list(itertools.product([0, .125], *([[0, 1]] * 6))))
    points = np.concatenate([boundaries, np.random.default_rng(871).uniform(size=(500, 7))])
    regions = validation.membership(points, fixture)
    accepted = np.logical_or.reduce(list(regions.values()))
    scored = [validation.evaluate(point, fixture) for point in points]
    assert np.array([r["distance"] <= 1 for r in scored]).tolist() == accepted.tolist()
    assert accepted[:len(boundaries)].all()
    assert all(r["simulator_dose_calls"] == 0 for r in scored)
    # A supplied intermediate threshold cannot redefine the final target.
    for point in points[:5]:
        assert validation.evaluate(point, fixture, threshold=.5) == validation.evaluate(point, fixture)


@pytest.mark.parametrize("fixture", validation.FIXTURES)
def test_analytic_truth_matches_independent_box_integrals(fixture):
    if fixture == "separated_boxes":
        boxes = [(np.zeros(7), np.full(7, .25)),
                 (np.array([.25] + [.625] * 6), np.array([1.] + [.875] * 6))]
        names = ["A", "B"]
    else:
        boxes = [(np.zeros(7), np.array([.125] + [1.] * 6))]
        names = ["slab"]
    volumes = np.array([np.prod(hi - lo) for lo, hi in boxes])
    proportions = volumes / volumes.sum()
    means = np.array([(lo + hi) / 2 for lo, hi in boxes])
    seconds = np.array([(lo**2 + lo * hi + hi**2) / 3 for lo, hi in boxes])
    expected = {**{f"mean_x{j+1}": proportions @ means[:, j] for j in range(7)},
                "second_x1": proportions @ seconds[:, 0],
                "second_x2": proportions @ seconds[:, 1],
                "cross_x1_x2": proportions @ (means[:, 0] * means[:, 1]),
                "cross_x2_x3": proportions @ (means[:, 1] * means[:, 2])}
    truth = validation.truth(fixture)
    assert truth["mass"] == pytest.approx(volumes.sum(), abs=1e-16)
    assert truth["mode_masses"] == pytest.approx(dict(zip(names, proportions)))
    assert truth["moments"] == pytest.approx(expected, abs=1e-15)


@pytest.mark.parametrize("fixture", validation.FIXTURES)
def test_exact_integration_nodes_pass_moment_and_mass_screens(fixture):
    points = quadrature_points(fixture)
    result = validation.summarize(points, np.full(len(points), -np.log(validation.truth(fixture)["mass"])), fixture)
    assert result["passed"]
    assert result["importance"]["n_attempts"] == len(points)
    assert max(result["absolute_moment_errors"].values()) < 1e-14
    assert max(result["absolute_mode_mass_errors"].values()) < 1e-14
    assert abs(result["relative_mass_error"]) < 1e-14


def test_rejected_production_attempts_remain_in_mass_denominator():
    accepted = quadrature_points("separated_boxes")
    points = np.concatenate([accepted, np.full_like(accepted, .5)])
    result = validation.summarize(points, np.full(len(points), np.log(4096)), "separated_boxes")
    assert result["importance"]["n_attempts"] == len(points)
    assert result["importance"]["n_accepted"] == len(accepted)
    assert result["importance"]["normalizer_estimate"] == pytest.approx(1 / 8192)
    assert not result["truth_checks"]["normalizing_mass"]
    assert result["truth_checks"]["moments"]


def test_inverse_density_weighting_not_raw_component_counts_defines_mode_mass():
    points = np.array([[.125] * 7, [.625] + [.75] * 6])
    result = validation.summarize(points, np.log([3., 1.]), "separated_boxes")
    assert result["mode_counts"] == {"A": 1, "B": 1}
    assert result["mode_masses"] == pytest.approx({"A": .25, "B": .75})
    assert result["moments"]["mean_x2"] == pytest.approx(19 / 32)


def test_zero_hits_are_explicit_failure_with_no_invented_moments():
    result = validation.summarize(np.full((300, 7), .5), np.zeros(300), "separated_boxes")
    assert not result["passed"] and not any(result["truth_checks"].values())
    assert result["moments"] is None and result["mode_masses"] is None
    assert result["importance"]["normalizer_estimate"] == 0
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("seed,usual_pass", [(2026092001, True), (2026092002, True), (2026092003, False)])
def test_required_omitted_region_control_falsifies_weight_screen_reassurance(seed, usual_pass):
    # These negative seeds were inspected during design; none is dropped. The
    # learned positive seeds are distinct and are not exercised in this test.
    result = validation.run_negative(seed)
    assessment = result["assessment"]
    assert assessment["mode_counts"]["A"] == 0
    assert assessment["importance"]["ess"] > 200
    assert all(assessment["usual_checks"].values()) is usual_pass
    assert not any(assessment["truth_checks"].values())
    assert not assessment["passed"]
    assert result["proposal"]["uniform_weight"] == .2


def test_oracle_precision_calibration_is_for_exact_target_proposal_only():
    bound = validation.oracle_precision_bound()
    mass = validation.truth("separated_boxes")["mass"]
    assert bound["accepted_density"] == pytest.approx(.2 + .8 / mass)
    assert bound["acceptance_probability"] == pytest.approx(bound["accepted_density"] * mass)
    assert bound["expected_accepted"] == pytest.approx(6554)
    assert bound["truth_error_probability_bound_per_run"] == pytest.approx(
        sum(bound["bound_terms"].values()))
    assert 0 < bound["truth_error_probability_bound_per_run"] < 3e-12
    assert "not a guarantee for learned proposals" in bound["scope"]


def test_production_stream_is_independent_and_excludes_pilot_points(monkeypatch):
    proposal = BoundedGaussianMixture(7, means=[[.125] * 7], scales=[[.1] * 7])
    pilot_points = [[.125] * 7] * 37
    consumed = []

    def fake_pilot(evaluate, dimension, epsilon, rng, plan):
        consumed.append(np.random.default_rng(rng).uniform(size=5000).sum())
        assert evaluate(pilot_points[0])["distance"] == 0
        return {"final_points": pilot_points, "plan": {"stub": True},
                "all_islands_reached_epsilon": True}

    def fake_fit(points, plan):
        assert points == pilot_points and plan == {"stub": True}
        return proposal

    monkeypatch.setattr(validation, "train_pilot", fake_pilot)
    monkeypatch.setattr(validation, "fit_proposal", fake_fit)
    seed, n = 1871, 128
    result = validation.run_positive("separated_boxes", seed, production_attempts=n)
    production_stream = np.random.SeedSequence(seed).spawn(2)[1]
    points, ids = proposal.sample(np.random.default_rng(production_stream), n)
    expected = validation.summarize(points, proposal.log_density(points), "separated_boxes")
    assert result["assessment"]["importance"] == expected["importance"]
    assert result["assessment"]["moments"] == expected["moments"]
    assert result["assessment"]["importance"]["n_attempts"] == n
    assert result["production_component_counts"] == np.bincount(ids, minlength=2).tolist()
    assert result["rng"] == {"pilot_spawn_key": [0], "production_spawn_key": [1]}
    assert len(consumed) == 1


def test_study_retains_failed_runs_and_renderer_names_failed_checks(monkeypatch):
    def fake_positive(fixture, seed):
        points = quadrature_points(fixture)
        a = validation.summarize(points, np.full(len(points), -np.log(validation.truth(fixture)["mass"])), fixture)
        a["pilot_reached_epsilon"] = False
        a["passed"] = False
        return {"fixture": fixture, "seed": seed, "assessment": a}

    monkeypatch.setattr(validation, "run_positive", fake_positive)
    result = validation.run_study()
    assert len(result["positive_runs"]) == 6 and len(result["negative_controls"]) == 3
    assert not result["passed"] and not result["checks"]["all_learned_runs_passed"]
    assert result["checks"]["all_negative_controls_rejected_by_truth"]
    assert result["checks"]["false_assurance_negative_control_observed"]
    rendered = validation.render(result)
    assert "**Prespecified synthetic checks failed. Failures are retained.**" in rendered
    assert rendered.count("pilot_reached_epsilon.") == 6
    assert "Negative-control seeds were inspected during fixture design" in rendered
    assert "cannot establish coverage of unknown biological acceptance regions" in rendered
    reversed_keys = json.loads(json.dumps(result), object_pairs_hook=lambda pairs: dict(reversed(pairs)))
    assert validation.render(reversed_keys) == rendered


def test_source_provenance_tracks_all_numerical_dependencies():
    hashes = validation.source_hashes()
    assert set(hashes) == {"scripts/proposal_synthetic_validation.py", "scripts/resample_move.py",
                           "scripts/bounded_proposal.py", "scripts/importance_sampling.py"}
    for name, digest in hashes.items():
        assert hashlib.sha256((validation.ROOT / name).read_bytes()).hexdigest() == digest


def test_source_drift_during_study_is_rejected(monkeypatch):
    calls = iter([{"source": "before"}, {"source": "after"}])
    monkeypatch.setattr(validation, "source_hashes", lambda: next(calls))
    monkeypatch.setattr(validation, "run_positive", lambda fixture, seed: {
        "assessment": {"importance": {"ess": 0}, "relative_mass_error": -1, "passed": False}})
    monkeypatch.setattr(validation, "run_negative", lambda seed: {})
    with pytest.raises(RuntimeError, match="sources changed"):
        validation.run_study()


@pytest.fixture
def stored_study(monkeypatch):
    """An artificial pilot trace for replay tests, with no adaptive experiment.

    Small clouds and production count only keep these structural tests cheap;
    the committed validation experiment always uses the untouched fixed plan.
    """
    plan = dict(validation.PILOT_PLAN, particles_per_island=4, levels=2, moves_per_level=1)
    monkeypatch.setattr(validation, "PILOT_PLAN", plan)
    monkeypatch.setattr(validation, "PRODUCTION_ATTEMPTS", 32)
    positives = []
    for fixture in validation.FIXTURES:
        center = [.125] * 7 if fixture == "separated_boxes" else [.0625] + [.5] * 6
        for seed in validation.POSITIVE_SEEDS:
            islands = []
            streams = np.random.SeedSequence(seed).spawn(2)[0].spawn(2)
            for island_id, stream in enumerate(streams):
                initial = np.random.default_rng(stream).uniform(size=(4, 7))
                final = [center] * 4
                islands.append({"island": island_id, "spawn_key": list(stream.spawn_key),
                                "initial_points": initial.tolist(),
                                "initial_distances": [validation.evaluate(p, fixture)["distance"] for p in initial],
                                "initial_simulator_dose_calls": 0,
                                "final_points": final,
                                "final_distances": [validation.evaluate(p, fixture)["distance"] for p in final],
                                "first_epsilon_level": 1,
                                "history": [{"level": level, "threshold": 1., "attempts": 4,
                                             "simulator_dose_calls": 0} for level in (1, 2)]})
            pilot = {"plan": plan, "islands": islands,
                     "final_points": [p for i in islands for p in i["final_points"]],
                     "final_distances": [d for i in islands for d in i["final_distances"]],
                     "all_islands_reached_epsilon": True, "attempts": 24, "simulator_dose_calls": 0}
            proposal = validation.fit_proposal(pilot["final_points"], plan)
            positives.append({"fixture": fixture, "seed": seed, "kind": "learned_proposal",
                              "pilot": pilot, "proposal": proposal.to_dict(), "assessment": {"stale": True},
                              "production_attempts": 32, "production_component_counts": [999],
                              "rng": {"pilot_spawn_key": [0], "production_spawn_key": [1]},
                              "wall_seconds": 17.5})
    return copy.deepcopy({"schema_version": 1, "dimension": 7, "epsilon": 1.,
            "source_hashes": validation.source_hashes(), "runtime": {"recorded": "versions"},
            "pilot_plan": plan, "production_attempts": 32,
            "positive_seeds": list(validation.POSITIVE_SEEDS), "negative_seeds": list(validation.NEGATIVE_SEEDS),
            "gates": dict(validation.GATES), "truth": {"stale": True}, "oracle_precision_bound": {},
            "positive_runs": positives,
            "negative_controls": [validation.run_negative(seed, 32) for seed in validation.NEGATIVE_SEEDS],
            "checks": {"stale": True}, "passed": True})


def test_assembly_replays_frozen_production_and_repairs_all_derived_fields(stored_study):
    original = copy.deepcopy(stored_study)
    stored_study["negative_controls"][0]["assessment"] = {"stale": True}
    result = validation.assemble(stored_study)
    run = result["positive_runs"][0]
    proposal = BoundedGaussianMixture.from_dict(run["proposal"])
    stream = np.random.SeedSequence(run["seed"]).spawn(2)[1]
    points, ids = proposal.sample(np.random.default_rng(stream), 32)
    expected = validation.summarize(points, proposal.log_density(points), run["fixture"])
    assert run["assessment"]["importance"] == expected["importance"]
    assert run["assessment"]["moments"] == expected["moments"]
    assert run["production_component_counts"] == np.bincount(ids, minlength=len(proposal.means) + 1).tolist()
    assert result["truth"] == {name: validation.truth(name) for name in validation.FIXTURES}
    assert result["oracle_precision_bound"] == validation.oracle_precision_bound()
    assert result["checks"] == validation.study_checks(result["positive_runs"], result["negative_controls"])
    assert not result["passed"]
    assert result["runtime"] == original["runtime"] and run["wall_seconds"] == 17.5
    assert stored_study["positive_runs"][0]["assessment"] == {"stale": True}
    assert validation.assemble(result) == result


def test_render_only_reassembles_before_writing_report(stored_study, monkeypatch, tmp_path):
    stored = tmp_path / "result.json"
    report = tmp_path / "report.md"
    stored.write_text(json.dumps(stored_study))
    monkeypatch.setattr(validation, "OUT_JSON", stored)
    monkeypatch.setattr(validation, "OUT_MD", report)
    monkeypatch.setattr(sys, "argv", ["proposal_synthetic_validation.py", "--render-only"])
    assert validation.main() == 0
    reassembled = validation.assemble(stored_study)
    assert report.read_text() == validation.render(reassembled)
    assert json.loads(stored.read_text()) == reassembled


@pytest.mark.parametrize("mutation,match", [
    (lambda r: r["positive_runs"][0]["proposal"]["means"][0].__setitem__(0, .3), "proposal means"),
    (lambda r: r["positive_runs"][0]["proposal"]["scales"][0].__setitem__(0, .2), "proposal scales"),
    (lambda r: r["positive_runs"][0].__setitem__("seed", 51), "run identities"),
    (lambda r: r["positive_runs"].append(copy.deepcopy(r["positive_runs"][0])), "run identities"),
    (lambda r: r["gates"].__setitem__("minimum_ess", 1), "gates"),
    (lambda r: r["pilot_plan"].__setitem__("levels", 1), "pilot_plan"),
    (lambda r: r.__setitem__("epsilon", 2), "epsilon"),
    (lambda r: r["source_hashes"].__setitem__("scripts/resample_move.py", "bad"), "source_hashes"),
    (lambda r: r["positive_runs"][0]["pilot"]["final_distances"].__setitem__(0, .7), "combined final_distances"),
    (lambda r: r["positive_runs"][0]["pilot"]["islands"][0]["final_distances"].__setitem__(0, .7), "final pilot scores"),
    (lambda r: r["positive_runs"][0]["pilot"]["islands"][0]["initial_points"][0].__setitem__(0, .5), "initial pilot points"),
    (lambda r: r["negative_controls"][0]["proposal"].__setitem__("uniform_weight", .9), "proposal uniform_weight"),
])
def test_assembly_rejects_changed_raw_inputs_and_fixed_criteria(stored_study, mutation, match):
    mutation(stored_study)
    with pytest.raises(ValueError, match=match):
        validation.assemble(stored_study)
