"""Runner checks with a toy target and artificial archives, never pilot training.

The prospective geometries and their reserved study seeds are not evaluated.
"""

import copy
import gzip
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import proposal_coverage_challenges as study  # noqa: E402


@pytest.fixture(autouse=True)
def forbid_pilot_training(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("unit tests must not train any prospective proposal")

    monkeypatch.setattr(study.strategy, "train_pilot", forbidden)


@pytest.fixture
def toy(monkeypatch):
    """Uniform target [0, 3/4] with an analytically known two-region split."""
    monkeypatch.setattr(study, "ATTEMPTS", 400)
    monkeypatch.setattr(study.geometry, "membership", lambda p, f: {
        "left": p[:, 0] < 0.5,
        "right": (p[:, 0] >= 0.5) & (p[:, 0] <= 0.75),
    })
    monkeypatch.setattr(study.geometry, "scores", lambda p, f: np.where(p[:, 0] <= 0.75, 0.5, 2.0))
    monkeypatch.setattr(study.geometry, "features", lambda p, f: {"mean": p[:, 0]})
    monkeypatch.setattr(study.geometry, "truth", lambda f: {
        "mass": 0.75, "region_masses": {"left": 2 / 3, "right": 1 / 3},
        "moments": {"mean": 3 / 8},
    })
    points = np.full((study.ATTEMPTS, 7), 0.5)
    points[:, 0] = (np.arange(study.ATTEMPTS) + 0.5) / study.ATTEMPTS
    return points


def test_summarize_counts_rejected_attempts_in_mass_and_mcse(toy):
    result = study.summarize(toy, np.zeros(len(toy)), "toy")
    assert result["passed"]
    assert result["region_counts"] == {"left": 200, "right": 100}
    assert result["importance"]["n_accepted"] == 300
    assert result["importance"]["ess"] == pytest.approx(300)
    assert result["importance"]["normalizer_estimate"] == pytest.approx(0.75)
    assert result["importance"]["max_normalized_weight"] == pytest.approx(1 / 300)
    assert result["importance"]["normalizer_relative_mcse"] > 0
    assert result["region_masses"] == pytest.approx({"left": 2 / 3, "right": 1 / 3})
    assert result["moments"] == pytest.approx({"mean": 3 / 8})


def test_support_hole_can_pass_all_weight_checks_while_failing_truth(toy):
    toy[:, 0] *= 0.5
    result = study.summarize(toy, np.full(len(toy), np.log(2)), "toy")
    assert all(result["usual_checks"].values())
    assert not any(result["truth_checks"].values())
    assert not result["passed"]
    assert result["region_counts"] == {"left": 400, "right": 0}
    assert result["importance"]["ess"] == pytest.approx(400)
    assert result["relative_mass_error"] == pytest.approx(-1 / 3)
    assert result["moments"] == pytest.approx({"mean": 0.25})


def test_nonuniform_weights_determine_region_mass_and_moments(toy):
    # 200 accepted points each contribute weight 1, 100 contribute weight 2.
    log_q = np.where(toy[:, 0] < 0.5, 0.0, -np.log(2))
    result = study.summarize(toy, log_q, "toy")
    assert result["importance"]["normalizer_estimate"] == pytest.approx(1.0)
    assert result["importance"]["ess"] == pytest.approx(400**2 / 600)
    assert result["region_masses"] == pytest.approx({"left": 0.5, "right": 0.5})
    assert result["moments"] == pytest.approx({"mean": (0.25 + 0.625) / 2})
    assert not result["passed"]


def test_zero_acceptances_are_an_explicit_failed_result(toy):
    toy[:, 0] = 0.9
    result = study.summarize(toy, np.zeros(len(toy)), "toy")
    assert not result["passed"]
    assert result["importance"]["n_accepted"] == 0
    assert result["relative_mass_error"] == -1
    for key in ("region_masses", "moments", "absolute_moment_errors", "absolute_region_mass_errors"):
        assert result[key] is None
    assert not any(result["usual_checks"].values())
    assert not any(result["truth_checks"].values())


@pytest.mark.parametrize("fault", ["missing_point", "extra_point", "dimension", "nan_point", "outside_cube",
                                  "missing_density", "nan_density", "infinite_density"])
def test_invalid_production_is_rejected(toy, fault):
    log_q = np.zeros(len(toy))
    if fault == "missing_point":
        toy = toy[:-1]
    elif fault == "extra_point":
        toy = np.vstack([toy, toy[0]])
    elif fault == "dimension":
        toy = toy[:, :-1]
    elif fault == "nan_point":
        toy[0, 0] = np.nan
    elif fault == "outside_cube":
        toy[0, 0] = -0.1
    elif fault == "missing_density":
        log_q = log_q[:-1]
    else:
        log_q[0] = np.nan if fault == "nan_density" else np.inf
    with pytest.raises(ValueError, match="production"):
        study.summarize(toy, log_q, "toy")


def test_region_overlap_and_score_disagreement_fail(toy, monkeypatch):
    monkeypatch.setattr(study.geometry, "membership", lambda p, f: {
        "left": np.ones(len(p), dtype=bool), "right": np.ones(len(p), dtype=bool)})
    with pytest.raises(ValueError, match="overlap"):
        study.summarize(toy, np.zeros(len(toy)), "toy")
    monkeypatch.setattr(study.geometry, "membership", lambda p, f: {
        "left": np.zeros(len(p), dtype=bool), "right": np.zeros(len(p), dtype=bool)})
    with pytest.raises(ValueError, match="disagree"):
        study.summarize(toy, np.zeros(len(toy)), "toy")


@pytest.mark.parametrize("kind", ["wrong_name", "wrong_shape", "nonfinite", "outside_bounds"])
def test_invalid_feature_contract_is_rejected(toy, monkeypatch, kind):
    def invalid(p, fixture):
        if kind == "wrong_name":
            return {"other": p[:, 0]}
        if kind == "wrong_shape":
            return {"mean": p[:, :1]}
        return {"mean": np.full(len(p), np.nan if kind == "nonfinite" else 1.1)}

    monkeypatch.setattr(study.geometry, "features", invalid)
    with pytest.raises(ValueError, match="features"):
        study.summarize(toy, np.zeros(len(toy)), "toy")


@pytest.fixture
def archives(monkeypatch, tmp_path, toy):
    """Real production replay around a fixed two-kernel proposal; no fitted pilot."""
    monkeypatch.setattr(study, "ARCHIVES", tmp_path)
    monkeypatch.setattr(study, "LEARNED_SEEDS", {"toy": (41101, 41102)})
    monkeypatch.setattr(study, "ORACLE_SEEDS", {"toy": (41201, 41202)})
    monkeypatch.setattr(study, "NEGATIVE_SEEDS", {"toy": (41301, 41302)})
    spec = {"revision": "unit-test", "source_hashes": {"mock_source": "1234"}}
    monkeypatch.setattr(study, "specification", lambda: copy.deepcopy(spec))
    means = [[0.25] * 7, [0.65] * 7]
    proposal = study.BoundedGaussianMixture(7, 0.2, means, [[0.2] * 7] * 2)
    monkeypatch.setattr(study.strategy, "fit_proposal", lambda *args: proposal)
    monkeypatch.setattr(study, "validate_pilot", lambda *args: True)

    def make(seed=41101):
        rng = np.random.default_rng(np.random.SeedSequence(seed).spawn(2)[1])
        points, ids = proposal.sample(rng, study.ATTEMPTS)
        scores = study.geometry.scores(points, "toy")
        return {"schema_version": 1, "fixture": "toy", "seed": seed, "specification": copy.deepcopy(spec),
                "runtime": {"python": "3.14.6", "numpy": np.__version__, "scipy": study.scipy.__version__,
                            "bit_generator": type(rng.bit_generator).__name__},
                "rng": {"pilot_spawn_key": [0], "production_spawn_key": [1]},
                "pilot": {"final_points": means, "attempts": 10,
                          "islands": [{"island": 0, "first_epsilon_level": 1, "history": []}]},
                "proposal": proposal.to_dict(), "wall_seconds": 1.25,
                "production": {"unit_points": points.tolist(), "component_ids": ids.tolist(),
                               "log_q": proposal.log_density(points).tolist(), "scores": scores.tolist(),
                               "accepted": (scores <= 1).tolist()}}

    def write(raw):
        path = study.archive_path(tmp_path, raw["fixture"], raw["seed"])
        path.write_bytes(gzip.compress(json.dumps(raw).encode(), mtime=0))
        return path

    return make, write, spec


def test_archive_assessment_replays_production_without_mutation(archives):
    make, _, spec = archives
    raw = make()
    before = copy.deepcopy(raw)
    result = study.assess_archive(raw, spec)
    assert result["importance"]["n_accepted"] == sum(raw["production"]["accepted"])
    assert result["pilot_reached_epsilon"]
    assert raw == before


def test_failed_pilot_prevents_pass_even_if_production_passes(archives, monkeypatch):
    make, _, spec = archives
    monkeypatch.setattr(study, "validate_pilot", lambda *args: False)
    monkeypatch.setattr(study, "summarize", lambda *args: {"passed": True})
    result = study.assess_archive(make(), spec)
    assert result == {"passed": False, "pilot_reached_epsilon": False}


def test_learned_runner_separates_pilot_and_production_streams_with_a_fake_pilot(archives, monkeypatch):
    make, _, spec = archives
    template = make()
    calls = []

    def fake_pilot(evaluate, dimension, epsilon, stream, *, plan):
        calls.append((dimension, epsilon, stream.entropy, stream.spawn_key, plan))
        return copy.deepcopy(template["pilot"])

    monkeypatch.setattr(study.strategy, "train_pilot", fake_pilot)
    result = study.run_learned("toy", 41101)
    assert calls == [(7, 1.0, 41101, (0,), study.strategy.PILOT_PLAN)]
    assert result["production"] == template["production"]
    assert result["rng"] == {"pilot_spawn_key": [0], "production_spawn_key": [1]}
    assert result["specification"] == spec


def test_learned_runner_rejects_sources_changed_during_fake_training(archives, monkeypatch):
    make, _, spec = archives
    calls = iter([spec, {**spec, "revision": "changed"}])
    monkeypatch.setattr(study, "specification", lambda: next(calls))
    monkeypatch.setattr(study.strategy, "train_pilot", lambda *a, **kw: make()["pilot"])
    with pytest.raises(ValueError, match="sources after learned run"):
        study.run_learned("toy", 41101)


@pytest.mark.parametrize("fault", ["schema_bool", "fixture", "unplanned_seed", "seed_float", "seed_bool",
    "source_hash", "streams", "elapsed_negative", "elapsed_bool", "elapsed_nan", "bit_generator",
    "proposal_extra", "proposal_dimension", "proposal_uniform", "proposal_mean", "proposal_scale",
    "points", "component_type", "component_value", "log_q", "scores", "accepted"])
def test_archive_corruption_is_rejected(archives, fault):
    make, _, spec = archives
    raw = make()
    if fault == "schema_bool": raw["schema_version"] = True
    elif fault == "fixture": raw["fixture"] = "unexpected"
    elif fault == "unplanned_seed": raw["seed"] = 999
    elif fault == "seed_float": raw["seed"] = 41101.0
    elif fault == "seed_bool": raw["seed"] = True
    elif fault == "source_hash": raw["specification"]["source_hashes"]["mock_source"] = "changed"
    elif fault == "streams": raw["rng"]["production_spawn_key"] = [0]
    elif fault == "elapsed_negative": raw["wall_seconds"] = -1
    elif fault == "elapsed_bool": raw["wall_seconds"] = True
    elif fault == "elapsed_nan": raw["wall_seconds"] = float("nan")
    elif fault == "bit_generator": raw["runtime"]["bit_generator"] = "MT19937"
    elif fault == "proposal_extra": raw["proposal"]["extra"] = 0
    elif fault == "proposal_dimension": raw["proposal"]["dimension"] = 7.0
    elif fault == "proposal_uniform": raw["proposal"]["uniform_weight"] = 0.3
    elif fault == "proposal_mean": raw["proposal"]["means"][0][0] += 0.01
    elif fault == "proposal_scale": raw["proposal"]["scales"][0][0] += 0.01
    elif fault == "points": raw["production"]["unit_points"][0][0] += 0.01
    elif fault == "component_type": raw["production"]["component_ids"][0] = float(raw["production"]["component_ids"][0])
    elif fault == "component_value": raw["production"]["component_ids"][0] += 1
    elif fault == "log_q": raw["production"]["log_q"][0] += 0.01
    elif fault == "scores": raw["production"]["scores"][0] += 0.01
    elif fault == "accepted": raw["production"]["accepted"][0] = not raw["production"]["accepted"][0]
    with pytest.raises(ValueError):
        study.assess_archive(raw, spec)


@pytest.mark.parametrize("fault", ["missing", "extra", "blank", "whitespace", "boolean", "integer", "null", "not_mapping"])
def test_runtime_provenance_requires_exact_nonblank_string_fields(archives, fault):
    make, _, spec = archives
    raw = make()
    if fault == "missing":
        del raw["runtime"]["python"]
    elif fault == "extra":
        raw["runtime"]["undeclared"] = "value"
    elif fault == "not_mapping":
        raw["runtime"] = []
    else:
        raw["runtime"]["python"] = {
            "blank": "", "whitespace": " \t", "boolean": True, "integer": 314, "null": None,
        }[fault]
    with pytest.raises(ValueError, match="runtime"):
        study.assess_archive(raw, spec)


def test_offline_runtime_replay_does_not_require_installed_version_match(archives):
    make, _, spec = archives
    raw = make()
    raw["runtime"].update(python="3.12.0", numpy="1.26.0", scipy="1.14.0")
    assert study.assess_archive(raw, spec)["pilot_reached_epsilon"]


@pytest.mark.parametrize("changed", [None, "prior", *study.SOURCE_PATHS[3:7]])
def test_specification_pins_prior_study_and_unchanged_sampler_sources(monkeypatch, tmp_path, changed):
    monkeypatch.setattr(study, "ROOT", tmp_path)
    prior = tmp_path / study.PRIOR_PATH
    prior.parent.mkdir(parents=True)
    old = {path: "frozen-" + str(i) for i, path in enumerate(study.SOURCE_PATHS[3:7])}
    prior.write_text(json.dumps({"source_hashes": old}))

    def source_hash(path):
        relative = path.relative_to(tmp_path).as_posix()
        if relative == study.PRIOR_PATH:
            return "changed" if changed == "prior" else study.PRIOR_SHA256
        if relative == changed:
            return "changed"
        return old.get(relative, "new-source-" + relative)

    monkeypatch.setattr(study, "_hash", source_hash)
    if changed is not None:
        with pytest.raises(ValueError, match="frozen.*changed"):
            study.specification()
    else:
        result = study.specification()
        assert set(result["source_hashes"]) == set(study.SOURCE_PATHS)
        assert result["prior_study"]["sha256"] == study.PRIOR_SHA256
        assert result["source_hashes"][study.PROTOCOL] == "new-source-" + study.PROTOCOL
        assert result["production_attempts"] == study.ATTEMPTS
        assert result["gates"] == study.GATES
        # Callers must not mutate a later specification through its plan.
        result["pilot_plan"]["levels"] = -1
        assert study.specification()["pilot_plan"]["levels"] > 0


@pytest.mark.parametrize("restricted", [False, True])
def test_control_keeps_its_own_density_and_support_label(toy, monkeypatch, restricted):
    calls = []

    def oracle(rng, n, fixture, *, restricted):
        calls.append((type(rng), n, fixture, restricted))
        p = toy.copy()
        p[:, 0] *= 0.5 if restricted else 0.75
        return p, np.full(n, -np.log(0.5 if restricted else 0.75))

    monkeypatch.setattr(study.geometry, "oracle_sample", oracle)
    monkeypatch.setattr(study.geometry, "oracle_retained_fraction", lambda f, r: 2 / 3 if r else 1.0)
    result = study.run_control("toy", 41401, restricted)
    assert calls == [(np.random.Generator, 400, "toy", restricted)]
    assert result["kind"] == ("support_hole_oracle" if restricted else "full_target_oracle")
    assert result["assessment"]["passed"] is (not restricted)
    assert result["retained_mass_fraction"] == (2 / 3 if restricted else 1.0)


@pytest.mark.parametrize("fault", ["missing", "extra", "misnamed"])
def test_report_requires_exact_archive_set(archives, tmp_path, fault):
    make, write, _ = archives
    first = write(make())
    if fault != "missing":
        write(make(41102))
    if fault == "extra":
        (tmp_path / "unexpected.json.gz").write_bytes(b"extra")
    if fault == "misnamed":
        first.rename(tmp_path / "wrong.json.gz")
    with pytest.raises(ValueError, match="every planned run"):
        study.build_report(tmp_path)


@pytest.mark.parametrize("fault", ["identity", "runtime"])
def test_report_rejects_wrong_file_identity_or_mixed_runtimes(archives, tmp_path, fault):
    make, write, _ = archives
    write(make())
    second = make(41102)
    path = write(second)
    if fault == "identity":
        second["seed"] = 41101
    else:
        second["runtime"]["numpy"] = "different"
    path.write_bytes(gzip.compress(json.dumps(second).encode(), mtime=0))
    with pytest.raises(ValueError, match="identity|runtimes"):
        study.build_report(tmp_path)


def _fake_control(fixture, seed, restricted):
    return {"fixture": fixture, "seed": seed,
            "assessment": {"passed": not restricted, "usual_checks": {"ess": True},
                           "truth_checks": {"support": not restricted}}}


def test_reassembly_discards_cached_results_and_rebuilds_every_assessment(archives, tmp_path, monkeypatch):
    make, write, spec = archives
    for seed in (41101, 41102):
        write(make(seed))
    calls = []

    def assess(raw, fixed):
        assert fixed == spec
        calls.append(raw["seed"])
        return {"passed": False, "recomputed": True}

    monkeypatch.setattr(study, "assess_archive", assess)
    monkeypatch.setattr(study, "run_control", _fake_control)
    stored = {"specification": spec, "archive_sha256": study._archive_hashes(tmp_path),
              "positive_runs": [{"assessment": {"passed": True}}], "passed": True}
    rebuilt = study.assemble(stored)
    assert calls == [41101, 41102]
    assert not rebuilt["passed"]
    assert all(run["assessment"]["recomputed"] for run in rebuilt["positive_runs"])
    assert len(rebuilt["oracle_controls"]) == len(rebuilt["negative_controls"]) == 2
    assert rebuilt["checks"]["all_nine_support_holes_fail_truth"]


@pytest.mark.parametrize("failure", ["learned", "oracle", "negative_usual", "negative_truth"])
def test_each_independent_study_gate_can_fail_the_report(archives, tmp_path, monkeypatch, failure):
    make, write, _ = archives
    for seed in (41101, 41102):
        write(make(seed))
    monkeypatch.setattr(study, "assess_archive", lambda *args: {"passed": failure != "learned"})

    def control(fixture, seed, restricted):
        result = _fake_control(fixture, seed, restricted)
        if failure == "oracle" and not restricted:
            result["assessment"]["passed"] = False
        if failure == "negative_usual" and restricted:
            result["assessment"]["usual_checks"]["ess"] = False
        if failure == "negative_truth" and restricted:
            result["assessment"]["truth_checks"]["support"] = True
        return result

    monkeypatch.setattr(study, "run_control", control)
    result = study.build_report(tmp_path)
    assert not result["passed"]
    assert sum(not passed for passed in result["checks"].values()) == 1
    assert len(result["positive_runs"]) == len(result["oracle_controls"]) == len(result["negative_controls"]) == 2


@pytest.mark.parametrize("field", ["specification", "archive_sha256"])
def test_reassembly_rejects_stale_provenance_before_rebuilding(archives, tmp_path, monkeypatch, field):
    make, write, spec = archives
    for seed in (41101, 41102): write(make(seed))
    stored = {"specification": spec, "archive_sha256": study._archive_hashes(tmp_path)}
    stored[field] = {}
    monkeypatch.setattr(study, "build_report", lambda *args: pytest.fail("stale report must fail first"))
    with pytest.raises(ValueError, match="report"):
        study.assemble(stored)


def test_report_detects_sources_changing_during_reconstruction(archives, tmp_path, monkeypatch):
    make, write, spec = archives
    for seed in (41101, 41102): write(make(seed))
    calls = iter([spec, {**spec, "revision": "changed"}])
    monkeypatch.setattr(study, "specification", lambda: next(calls))
    monkeypatch.setattr(study, "run_control", _fake_control)
    with pytest.raises(ValueError, match="sources after"):
        study.build_report(tmp_path)


def test_report_detects_archive_bytes_changing_during_reconstruction(archives, tmp_path, monkeypatch):
    make, write, _ = archives
    paths = [write(make(seed)) for seed in (41101, 41102)]

    def mutating_control(*args):
        paths[0].write_bytes(paths[0].read_bytes() + b"changed")
        return _fake_control(*args)

    monkeypatch.setattr(study, "run_control", mutating_control)
    with pytest.raises(ValueError, match="archives after"):
        study.build_report(tmp_path)


@pytest.fixture
def cli(archives, monkeypatch, tmp_path):
    monkeypatch.setattr(study, "OUT_JSON", tmp_path / "report.json")
    monkeypatch.setattr(study, "OUT_MD", tmp_path / "report.md")
    study.OUT_JSON.write_text("original JSON")
    study.OUT_MD.write_text("original Markdown")
    monkeypatch.setattr(study, "render", lambda result: "rebuilt Markdown")
    monkeypatch.setattr(study, "build_report", lambda: {"passed": False, "fresh": True})
    return archives


def test_cli_refuses_existing_archive_without_resume(cli, monkeypatch):
    make, write, _ = cli
    path = write(make())
    original = path.read_bytes()
    monkeypatch.setattr(sys, "argv", ["runner"])
    monkeypatch.setattr(study, "run_learned", lambda *args: pytest.fail("must not rerun existing study"))
    with pytest.raises(FileExistsError, match="resume"):
        study.main()
    assert path.read_bytes() == original
    assert study.OUT_JSON.read_text() == "original JSON"
    assert study.OUT_MD.read_text() == "original Markdown"


def test_cli_resume_retains_verified_bytes_and_runs_only_missing_seed(cli, monkeypatch):
    make, write, _ = cli
    path = write(make())
    original = path.read_bytes()
    called = []
    monkeypatch.setattr(sys, "argv", ["runner", "--resume"])

    def fake_run(fixture, seed):
        called.append((fixture, seed))
        return make(seed)

    monkeypatch.setattr(study, "run_learned", fake_run)
    assert study.main() == 0
    assert called == [("toy", 41102)]
    assert path.read_bytes() == original
    second = study.archive_path(study.ARCHIVES, "toy", 41102)
    assert json.loads(gzip.decompress(second.read_bytes()))["seed"] == 41102
    assert json.loads(study.OUT_JSON.read_text()) == {"passed": False, "fresh": True}


def test_cli_resume_rejects_corrupt_archive_without_replacing_outputs(cli, monkeypatch):
    make, write, _ = cli
    raw = make()
    raw["rng"]["production_spawn_key"] = [0]
    path = write(raw)
    before = path.read_bytes()
    monkeypatch.setattr(sys, "argv", ["runner", "--resume"])
    monkeypatch.setattr(study, "run_learned", lambda *args: pytest.fail("do not resume after corrupt archive"))
    with pytest.raises(ValueError, match="streams"):
        study.main()
    assert path.read_bytes() == before
    assert study.OUT_JSON.read_text() == "original JSON"
    assert study.OUT_MD.read_text() == "original Markdown"


def test_cli_exclusive_create_does_not_clobber_concurrently_created_archive(cli, monkeypatch):
    make, _, _ = cli
    path = study.archive_path(study.ARCHIVES, "toy", 41101)
    monkeypatch.setattr(sys, "argv", ["runner"])

    def race(fixture, seed):
        path.write_bytes(b"concurrent writer")
        return make(seed)

    monkeypatch.setattr(study, "run_learned", race)
    with pytest.raises(FileExistsError):
        study.main()
    assert path.read_bytes() == b"concurrent writer"
    assert study.OUT_JSON.read_text() == "original JSON"


def test_cli_render_only_uses_archive_reassembly_not_training(cli, monkeypatch):
    study.OUT_JSON.write_text(json.dumps({"stored": "metadata"}))
    monkeypatch.setattr(sys, "argv", ["runner", "--render-only"])
    seen = []

    def reassemble(stored):
        seen.append(stored)
        return {"passed": False, "reconstructed": True}

    monkeypatch.setattr(study, "assemble", reassemble)
    monkeypatch.setattr(study, "run_learned", lambda *args: pytest.fail("render-only cannot train"))
    assert study.main() == 0
    assert seen == [{"stored": "metadata"}]
    assert json.loads(study.OUT_JSON.read_text())["reconstructed"]


def test_cli_render_failure_preserves_both_existing_reports(cli, monkeypatch):
    make, write, _ = cli
    for seed in (41101, 41102):
        write(make(seed))
    monkeypatch.setattr(sys, "argv", ["runner", "--resume"])

    def broken_render(result):
        raise RuntimeError("render failed")

    monkeypatch.setattr(study, "render", broken_render)
    with pytest.raises(RuntimeError, match="render failed"):
        study.main()
    assert study.OUT_JSON.read_text() == "original JSON"
    assert study.OUT_MD.read_text() == "original Markdown"


def test_render_retains_failed_runs_and_labels_unavailable_summaries(toy):
    toy[:, 0] = 0.9
    assessment = study.summarize(toy, np.zeros(len(toy)), "toy")
    assessment["pilot_reached_epsilon"] = False
    result = {"passed": False, "positive_runs": [{"fixture": "toy", "seed": 41101, "assessment": assessment}],
              "oracle_controls": [], "negative_controls": [], "checks": {"all_nine_learned_runs_pass": False}}
    text = study.render(result)
    assert "Prespecified checks failed" in text
    assert "| toy | 41101 | 0 |" in text
    assert "unavailable | unavailable | fail" in text
    assert "pilot_reached_epsilon" in text
    assert "every_region_observed" in text
    assert "None." not in text
