"""Revision orchestration tests with explicit artificial pilots, never v2 runs."""

import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bounded_proposal import BoundedGaussianMixture  # noqa: E402
import proposal_synthetic_validation as frozen  # noqa: E402
import proposal_synthetic_validation_v2 as wrapper  # noqa: E402
import synthetic_proposal_study as study  # noqa: E402


def fake_strategy():
    """A trace double, not a competing pilot or synthetic-design experiment."""
    plan = dict(frozen.PILOT_PLAN, particles_per_island=4, levels=2, moves_per_level=1)
    strategy = SimpleNamespace(__name__="explicit_test_strategy", PILOT_PLAN=plan,
                               train_calls=[], validation_calls=[])

    def train(evaluate, dimension, epsilon, rng, plan):
        strategy.train_calls.append(tuple(rng.spawn_key))
        islands = []
        for index, stream in enumerate(rng.spawn(plan["islands"])):
            island_rng = np.random.default_rng(stream)
            initial = island_rng.uniform(size=(plan["particles_per_island"], dimension))
            # Additional pilot consumption must not affect the production stream.
            island_rng.uniform(size=1000)
            scores = [evaluate(point)["distance"] for point in initial]
            first = max(epsilon, float(np.quantile(scores, .5, method="higher")))
            center = [.0625] + [.125] * (dimension - 1)
            final = [center] * plan["particles_per_island"]
            final_scores = [evaluate(point)["distance"] for point in final]
            ancestor = int(np.argmin(scores))
            history = []
            for level, threshold in enumerate([first, epsilon], 1):
                retained = sum(score <= first for score in scores) if level == 1 else len(final)
                history.append({"level": level, "threshold": threshold,
                                "attempts": 4, "evaluated": 4, "mh_rejected_without_evaluation": 0,
                                "early_rejected": 0, "accepted_moves": 4, "simulator_dose_calls": 0,
                                "retained_before_resampling": retained, "resampled": retained < len(final),
                                "unique_ancestors": 1, "unique_particles": 1,
                                "minimum_distance": min(final_scores), "median_distance": float(np.median(final_scores)),
                                "local_attempts": 2, "global_attempts": 2})
            islands.append({"island": index, "spawn_key": list(stream.spawn_key),
                            "initial_points": initial.tolist(), "initial_distances": scores,
                            "initial_simulator_dose_calls": 0, "final_points": final,
                            "final_distances": final_scores, "final_ancestors": [ancestor] * len(final),
                            "history": history, "first_epsilon_level": 1 if first == epsilon else 2})
        return {"plan": copy.deepcopy(plan), "islands": islands,
                "final_points": [point for island in islands for point in island["final_points"]],
                "final_distances": [score for island in islands for score in island["final_distances"]],
                "all_islands_reached_epsilon": True, "attempts": 24, "simulator_dose_calls": 0}

    def fit(points, plan):
        points = np.asarray(points)
        return BoundedGaussianMixture(7, uniform_weight=.2, means=points,
                                      scales=np.full_like(points, .025))

    def validate(pilot, fixture, seed):
        strategy.validation_calls.append((fixture, seed))
        for island in pilot["islands"]:
            for step in island["history"]:
                if step["local_attempts"] + step["global_attempts"] != step["attempts"]:
                    raise ValueError("local/global attempts disagree")

    strategy.train_pilot, strategy.fit_proposal, strategy.validate_pilot_archive = train, fit, validate
    return strategy


@pytest.fixture
def spec():
    return replace(wrapper.SPEC, revision_id="unit-test-revision", strategy_module="explicit_test_strategy",
                   positive_seeds=(10101, 10102, 10103),
                   source_paths=tuple(path for path in wrapper.SPEC.source_paths if path != "scripts/resample_move_local.py") +
                   ("tests/test_synthetic_proposal_study.py",))


@pytest.fixture
def saved(spec):
    strategy = fake_strategy()
    result = study.run_study(spec, strategy)
    return result, strategy


def test_v2_has_distinct_explicit_identity_sources_paths_and_reserved_seeds():
    assert wrapper.SPEC.positive_seeds == (2026092401, 2026092402, 2026092403)
    assert wrapper.SPEC.strategy_module == "resample_move_local"
    assert wrapper.OUT_JSON != frozen.OUT_JSON and wrapper.OUT_MD != frozen.OUT_MD
    assert set(wrapper.SPEC.source_paths) == {"scripts/proposal_synthetic_validation_v2.py",
        "scripts/synthetic_proposal_study.py", "scripts/resample_move_local.py",
        "scripts/proposal_synthetic_validation.py", "scripts/resample_move.py",
        "scripts/bounded_proposal.py", "scripts/importance_sampling.py"}
    assert hashlib.sha256((study.ROOT / wrapper.SPEC.prior_study_path).read_bytes()).hexdigest() == wrapper.SPEC.prior_study_sha256


def test_production_is_independent_replay_without_pilot_points_or_held_back_v2_seeds(saved, spec):
    result, strategy = saved
    assert len(strategy.train_calls) == 6
    assert len(result["positive_runs"]) == 6 and len(result["negative_controls"]) == 3
    for run in result["positive_runs"]:
        proposal = BoundedGaussianMixture.from_dict(run["proposal"])
        stream = np.random.SeedSequence(run["seed"]).spawn(2)[1]
        points, ids = proposal.sample(np.random.default_rng(stream), 8192)
        independent = frozen.summarize(points, proposal.log_density(points), run["fixture"])
        assert run["assessment"]["importance"] == independent["importance"]
        assert run["assessment"]["moments"] == independent["moments"]
        assert run["assessment"]["importance"]["n_attempts"] == 8192
        assert run["production_component_counts"] == np.bincount(ids, minlength=9).tolist()
        assert run["seed"] in spec.positive_seeds and run["seed"] not in wrapper.SPEC.positive_seeds
    assert not result["passed"]  # The deliberately narrow fake proposal is inadequate.
    assert result["checks"]["all_negative_controls_rejected_by_truth"]


def test_reassembly_repairs_cached_results_without_training_or_mutating_history(saved, spec):
    result, strategy = saved
    stale = copy.deepcopy(result)
    stale["positive_runs"][0]["assessment"] = {"stale": True}
    stale["positive_runs"][0]["production_component_counts"] = [999]
    stale["positive_runs"][0]["pilot"]["all_islands_reached_epsilon"] = False
    stale["negative_controls"][0]["assessment"] = {"stale": True}
    stale["truth"] = stale["checks"] = stale["oracle_precision_bound"] = {}
    stale["passed"] = True
    rebuilt = study.assemble(stale, spec, strategy)
    assert rebuilt == result
    assert len(strategy.train_calls) == 6
    assert stale["positive_runs"][0]["assessment"] == {"stale": True}
    assert rebuilt["runtime"] == stale["runtime"]
    assert study.assemble(rebuilt, spec, strategy) == rebuilt


def test_failed_results_and_development_scope_survive_dictionary_reordering(saved):
    result, _ = saved
    rendered = study.render(result)
    assert "**Fixed synthetic checks failed. All outcomes are retained.**" in rendered
    assert "Failed learned-run checks:" in rendered
    assert "Fresh random seeds do not make these" in rendered
    assert "known failure controls, not fresh validation runs" in rendered
    assert result["prior_study"]["sha256"] in rendered
    reordered = json.loads(json.dumps(result), object_pairs_hook=lambda pairs: dict(reversed(pairs)))
    assert study.render(reordered) == rendered


def test_wrapper_lazy_strategy_uses_its_fixed_specification(monkeypatch):
    calls = []
    sentinel = object()

    def importer(name):
        calls.append(name)
        return sentinel

    monkeypatch.setattr(wrapper.importlib, "import_module", importer)
    assert wrapper.strategy() is sentinel
    assert calls == [wrapper.SPEC.strategy_module]


def test_render_only_repairs_both_artifacts_through_explicit_strategy(saved, spec, monkeypatch, tmp_path):
    result, strategy = saved
    stale = copy.deepcopy(result)
    stale["positive_runs"][0]["assessment"] = {"stale": True}
    path, md = tmp_path / "v2.json", tmp_path / "v2.md"
    path.write_text(json.dumps(stale))
    monkeypatch.setattr(wrapper, "SPEC", spec)
    monkeypatch.setattr(wrapper, "strategy", lambda: strategy)
    monkeypatch.setattr(wrapper, "OUT_JSON", path)
    monkeypatch.setattr(wrapper, "OUT_MD", md)
    monkeypatch.setattr(sys, "argv", ["proposal_synthetic_validation_v2.py", "--render-only"])
    assert wrapper.main() == 0
    assert json.loads(path.read_text()) == result
    assert md.read_text() == study.render(result)
    assert len(strategy.train_calls) == 6


@pytest.mark.parametrize("mutation,match", [
    (lambda r: r.__setitem__("revision_id", "other"), "revision_id"),
    (lambda r: r.__setitem__("strategy_module", "arbitrary_module"), "strategy_module"),
    (lambda r: r["positive_seeds"].__setitem__(0, 9), "positive_seeds"),
    (lambda r: r["positive_runs"][0].__setitem__("seed", 9), "run identities"),
    (lambda r: r["positive_runs"].append(copy.deepcopy(r["positive_runs"][0])), "run identities"),
    (lambda r: r["gates"].__setitem__("minimum_ess", 1), "gates"),
    (lambda r: r["pilot_plan"].__setitem__("levels", 100), "pilot_plan"),
    (lambda r: r["prior_study"].__setitem__("sha256", "bad"), "prior_study"),
    (lambda r: r["source_hashes"].__setitem__("scripts/bounded_proposal.py", "bad"), "source_hashes"),
    (lambda r: r["positive_runs"][0]["proposal"]["means"][0].__setitem__(0, .9), "proposal means"),
    (lambda r: r["negative_controls"][0]["proposal"].__setitem__("uniform_weight", .9), "negative proposal"),
])
def test_changed_revision_raw_proposal_or_provenance_is_rejected(saved, spec, mutation, match):
    result, strategy = saved
    mutation(result)
    with pytest.raises(ValueError, match=match):
        study.assemble(result, spec, strategy)


@pytest.mark.parametrize("mutation,match", [
    (lambda p: p["initial_points"][0].__setitem__(0, .9), "initial pilot points"),
    (lambda p: p["final_distances"].__setitem__(0, .9), "final pilot scores"),
    (lambda p: p["final_ancestors"].__setitem__(0, -1), "ancestry"),
    (lambda p: p["history"][0].__setitem__("threshold", .1), "thresholds"),
    (lambda p: p["history"][0].__setitem__("evaluated", 0), "move accounting"),
    (lambda p: p["history"][0].__setitem__("accepted_moves", 100), "accepted moves"),
    (lambda p: p["history"][0].__setitem__("simulator_dose_calls", 1), "analytic move cost"),
    (lambda p: p["history"][0].__setitem__("local_attempts", 99), "local/global attempts"),
    (lambda p: p["history"][0].__setitem__("minimum_distance", np.nan), "distance summaries"),
    (lambda p: p["history"][0].__setitem__("median_distance", np.inf), "distance summaries"),
    (lambda p: p["history"][0].__setitem__("minimum_distance", -.1), "distance summaries"),
    (lambda p: p["history"][0].__setitem__("minimum_distance", p["history"][0]["median_distance"] + .1), "distance summaries"),
    (lambda p: p["history"][0].__setitem__("median_distance", p["history"][0]["threshold"] + .1), "distance summaries"),
    (lambda p: p["history"][0].__setitem__("resampled", 1), "resampling flag must be a boolean"),
    (lambda p: p["history"][-1].__setitem__("unique_particles", 3), "final particle diversity"),
])
def test_common_pilot_checks_and_strategy_hook_reject_inconsistent_trace(saved, spec, mutation, match):
    result, strategy = saved
    mutation(result["positive_runs"][0]["pilot"]["islands"][0])
    with pytest.raises(ValueError, match=match):
        study.assemble(result, spec, strategy)


def test_prior_artifact_or_duplicated_or_reused_seeds_fail_before_training(spec):
    strategy = fake_strategy()
    for bad, match in [(replace(spec, prior_study_sha256="bad"), "prior study SHA256"),
                       (replace(spec, positive_seeds=(1, 1, 2)), "distinct"),
                       (replace(spec, positive_seeds=frozen.POSITIVE_SEEDS), "prior study's seeds")]:
        with pytest.raises(ValueError, match=match):
            study.run_study(bad, strategy)
    assert strategy.train_calls == []


def test_numerical_source_drift_during_a_run_fails_closed(spec, monkeypatch):
    actual = study.source_hashes(spec)
    calls = iter([actual, {**actual, "scripts/bounded_proposal.py": "changed"}])
    monkeypatch.setattr(study, "source_hashes", lambda spec: next(calls))
    with pytest.raises(ValueError, match="specification after study"):
        study.run_study(spec, fake_strategy())


def test_small_real_pilot_adapter_matches_common_archive_contract_without_production():
    import resample_move_local as local

    plan = dict(local.PILOT_PLAN, particles_per_island=8, levels=3, moves_per_level=1)
    adapter = SimpleNamespace(PILOT_PLAN=plan, validate_pilot_archive=local.validate_pilot_archive)
    seed, fixture = 17054, "boundary_slab"
    stream = np.random.SeedSequence(seed).spawn(2)[0]
    pilot = local.train_pilot(lambda point, threshold=None: frozen.evaluate(point, fixture, threshold),
                              7, 1., stream, plan=plan)
    reached = study.check_pilot_archive(pilot, fixture, seed, adapter)
    assert reached == pilot["all_islands_reached_epsilon"]
    assert pilot["attempts"] == 2 * (8 + 3 * 8)
    assert pilot["simulator_dose_calls"] == 0
    proposal = local.fit_proposal(pilot["final_points"], plan)
    restored = frozen._checked_proposal(proposal.to_dict(), local.fit_proposal(pilot["final_points"], plan))
    assert restored.log_density(pilot["final_points"]) == pytest.approx(proposal.log_density(pilot["final_points"]))
