"""Production revision contracts using toy targets and artificial pilots only.

No prospective geometry, reserved study seed, or actual pilot is evaluated.
"""

import copy
import gzip
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import proposal_correlated_study as study  # noqa: E402


@pytest.fixture(autouse=True)
def forbid_pilot_training(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("unit tests must never train a prospective proposal")

    monkeypatch.setattr(study.strategy, "train_pilot", forbidden)


@pytest.fixture
def toy(monkeypatch):
    """Exact midpoint grid from a uniform proposal on [-1/2, 3/2] × [0,1]^6."""
    monkeypatch.setattr(study, "ATTEMPTS", 800)

    def membership(points, fixture):
        inside = ((points >= 0) & (points <= 1)).all(axis=1)
        return {"left": inside & (points[:, 0] < 0.5),
                "right": inside & (points[:, 0] >= 0.5) & (points[:, 0] <= 0.75)}

    monkeypatch.setattr(study.geometry, "membership", membership)
    # Some outside-cube points have a passing score, testing the explicit cube mask.
    monkeypatch.setattr(study.geometry, "scores", lambda p, f: np.where(p[:, 0] <= 0.75, 0.5, 2.0))
    monkeypatch.setattr(study.geometry, "features", lambda p, f: {"mean": p[:, 0]})
    monkeypatch.setattr(study.geometry, "truth", lambda f: {
        "mass": 0.75, "region_masses": {"left": 2 / 3, "right": 1 / 3},
        "moments": {"mean": 3 / 8},
    })
    points = np.full((study.ATTEMPTS, 7), 0.5)
    points[:, 0] = -0.5 + 2 * (np.arange(study.ATTEMPTS) + 0.5) / study.ATTEMPTS
    return points


def test_outside_cube_attempts_keep_zero_weights_and_remain_in_denominator(toy):
    result = study.summarize(toy, np.full(len(toy), -np.log(2)), "toy")
    diagnostics = result["importance"]
    assert result["passed"]
    assert diagnostics["n_attempts"] == 800
    assert diagnostics["n_accepted"] == 300
    assert diagnostics["normalizer_estimate"] == pytest.approx(0.75)
    assert diagnostics["ess"] == pytest.approx(300)
    assert diagnostics["normalizer_relative_mcse"] == pytest.approx(np.sqrt((800 / 300 - 1) / 799))
    assert result["region_counts"] == {"left": 200, "right": 100}
    assert result["region_masses"] == pytest.approx({"left": 2 / 3, "right": 1 / 3})
    assert result["moments"] == pytest.approx({"mean": 3 / 8})


def test_nonuniform_production_weights_determine_mass_and_moments(toy):
    log_q = np.where(toy[:, 0] < 0.5, -np.log(2), -np.log(4))
    result = study.summarize(toy, log_q, "toy")
    assert result["importance"]["normalizer_estimate"] == pytest.approx(1)
    assert result["importance"]["ess"] == pytest.approx(800**2 / 2400)
    assert result["region_masses"] == pytest.approx({"left": 0.5, "right": 0.5})
    assert result["moments"] == pytest.approx({"mean": (0.25 + 0.625) / 2})
    assert not result["passed"]


def test_zero_hits_remain_a_failed_assessment(toy):
    toy[:, 0] = -0.1
    result = study.summarize(toy, np.zeros(len(toy)), "toy")
    assert not result["passed"]
    assert result["importance"]["n_attempts"] == len(toy)
    assert result["importance"]["n_accepted"] == 0
    assert result["relative_mass_error"] == -1
    assert result["moments"] is None
    assert result["region_masses"] is None


def test_three_proposals_share_centers_and_ablation_retains_covariance_diagonal():
    # A correlated artificial endpoint cloud with duplicate ancestry copies.
    first = np.linspace(0.2, 0.8, 12)
    points = np.full((12, 7), 0.5)
    points[:, 0] = first
    points[:, 1] = 0.15 + 0.7 * first
    points[:, 2] = 0.75 - 0.5 * first
    points = np.vstack([points, points[3], points[3]])
    original = points.copy()
    fitted = study.proposals(points)
    assert set(fitted) == {"bounded_diagonal", "unbounded_diagonal", "correlated"}
    assert np.array_equal(points, original)
    for proposal in fitted.values():
        assert np.array_equal(proposal.means, points)
        assert proposal.uniform_weight == 0.2
    correlated = np.asarray(fitted["correlated"].covariances)
    diagonal = np.asarray(fitted["unbounded_diagonal"].covariances)
    indices = np.arange(7)
    assert np.allclose(diagonal[:, indices, indices], correlated[:, indices, indices], rtol=0, atol=0)
    expected_diagonal = np.zeros_like(correlated)
    expected_diagonal[:, indices, indices] = correlated[:, indices, indices]
    assert np.array_equal(diagonal, expected_diagonal)
    assert np.max(np.abs(correlated - diagonal)) > 1e-4
    assert np.array_equal(correlated[-1], correlated[-2])


@pytest.fixture
def archives(monkeypatch, tmp_path, toy):
    """Small fixed proposal kernels; no pilot trajectory or geometry is fitted."""
    monkeypatch.setattr(study, "ARCHIVES", tmp_path)
    monkeypatch.setattr(study, "LEARNED_SEEDS", {"toy": (41101, 41102)})
    monkeypatch.setattr(study, "ORACLE_SEEDS", {"toy": (41201, 41202)})
    monkeypatch.setattr(study, "NEGATIVE_SEEDS", {"toy": (41301, 41302)})
    spec = {"revision": "toy-unit-test", "source_hashes": {"mock-source": "1234"}}
    monkeypatch.setattr(study, "specification", lambda: copy.deepcopy(spec))
    monkeypatch.setattr(study, "committed_sources", lambda: "a" * 40)
    means = [[0.25] * 7, [0.65] * 7]
    covariance = np.eye(7) * 0.04
    covariance[0, 1] = covariance[1, 0] = 0.015
    fitted = {
        "bounded_diagonal": study.BoundedGaussianMixture(7, 0.2, means, [[0.2] * 7] * 2),
        "unbounded_diagonal": study.GaussianMixture(7, 0.2, means, [np.eye(7) * 0.04] * 2),
        "correlated": study.GaussianMixture(7, 0.2, means, [covariance] * 2),
    }
    monkeypatch.setattr(study, "proposals", lambda points: fitted)
    monkeypatch.setattr(study, "validate_pilot", lambda *args: True)

    def make(seed=41101):
        streams = np.random.SeedSequence(seed).spawn(4)
        arms = {}
        for index, name in enumerate(study.ARMS, 1):
            proposal = fitted[name]
            points, ids = proposal.sample(np.random.default_rng(streams[index]), study.ATTEMPTS)
            scores = study.geometry.scores(points, "toy")
            inside = ((points >= 0) & (points <= 1)).all(axis=1)
            arms[name] = {
                "spawn_key": list(streams[index].spawn_key), "proposal": proposal.to_dict(),
                "production": {"points": points.tolist(), "component_ids": ids.tolist(),
                               "log_q": proposal.log_density(points).tolist(), "scores": scores.tolist(),
                               "accepted": (inside & (scores <= study.EPSILON)).tolist()},
            }
        return {
            "schema_version": 1, "fixture": "toy", "seed": seed,
            "specification": copy.deepcopy(spec), "implementation_commit": "a" * 40,
            "runtime": {"python": "3.14.6", "numpy": np.__version__,
                        "scipy": study.scipy.__version__, "bit_generator": "PCG64"},
            "pilot": {"final_points": means, "attempts": 10},
            "arms": arms, "wall_seconds": 1.25,
        }

    def write(raw):
        path = study.archive_path(tmp_path, raw["fixture"], raw["seed"])
        path.write_bytes(gzip.compress(json.dumps(raw).encode(), mtime=0))
        return path

    return make, write, spec, fitted


def test_archive_replay_rebuilds_all_three_arms_without_mutating_raw_data(archives):
    make, _, spec, _ = archives
    raw = make()
    before = copy.deepcopy(raw)
    results = study.assess_archive(raw, spec)
    assert tuple(results) == study.ARMS
    assert raw == before
    for name, result in results.items():
        assert result["importance"]["n_accepted"] == sum(raw["arms"][name]["production"]["accepted"])
        assert result["importance"]["n_attempts"] == study.ATTEMPTS
        assert result["pilot_reached_epsilon"]
    assert results["bounded_diagonal"]["outside_cube"] == 0
    assert results["unbounded_diagonal"]["outside_cube"] > 0
    assert results["correlated"]["outside_cube"] > 0


def test_full_mixture_density_is_required_instead_of_sampled_component_density(archives):
    from scipy.stats import multivariate_normal

    make, _, spec, fitted = archives
    raw = make()
    production = raw["arms"]["correlated"]["production"]
    points, ids = np.asarray(production["points"]), np.asarray(production["component_ids"])
    selected_log_q = np.full(len(points), np.log(0.2))
    proposal = fitted["correlated"]
    for index, (mean, covariance) in enumerate(zip(proposal.means, proposal.covariances), 1):
        selected = ids == index
        selected_log_q[selected] = np.log(0.4) + multivariate_normal.logpdf(
            points[selected], mean=mean, cov=covariance)
    full_log_q = np.asarray(production["log_q"])
    assert np.all(full_log_q >= selected_log_q - 1e-12)
    assert np.max(full_log_q - selected_log_q) > 0.1
    production["log_q"] = selected_log_q.tolist()
    with pytest.raises(ValueError, match="full mixture density"):
        study.assess_archive(raw, spec)


def test_pilot_and_each_production_arm_use_four_distinct_frozen_streams(archives, monkeypatch):
    make, _, spec, _ = archives
    template = make()
    calls = []

    def fake_pilot(evaluate, dimension, epsilon, stream, *, plan):
        calls.append((dimension, epsilon, stream.entropy, stream.spawn_key, plan))
        return copy.deepcopy(template["pilot"])

    monkeypatch.setattr(study.strategy, "train_pilot", fake_pilot)
    raw = study.run_learned("toy", 41101)
    assert calls == [(7, 1.0, 41101, (0,), study.strategy.PILOT_PLAN)]
    assert raw["arms"] == template["arms"]
    assert raw["specification"] == spec
    assert [raw["arms"][name]["spawn_key"] for name in study.ARMS] == [[1], [2], [3]]
    assert not np.array_equal(raw["arms"]["unbounded_diagonal"]["production"]["points"],
                              raw["arms"]["correlated"]["production"]["points"])


def test_uncommitted_sources_prevent_any_pilot_or_production(archives, monkeypatch):
    def uncommitted():
        raise ValueError("sources are not committed")

    monkeypatch.setattr(study, "committed_sources", uncommitted)
    with pytest.raises(ValueError, match="not committed"):
        study.run_learned("toy", 41101)


def test_source_change_during_fake_pilot_is_rejected(archives, monkeypatch):
    make, _, spec, _ = archives
    calls = iter([spec, {**spec, "revision": "changed"}])
    monkeypatch.setattr(study, "specification", lambda: next(calls))
    monkeypatch.setattr(study.strategy, "train_pilot", lambda *a, **kw: make()["pilot"])
    with pytest.raises(ValueError, match="sources after run"):
        study.run_learned("toy", 41101)


def test_failed_pilot_rejects_all_arms_even_when_production_passes(archives, monkeypatch):
    make, _, spec, _ = archives
    monkeypatch.setattr(study, "validate_pilot", lambda *args: False)
    monkeypatch.setattr(study, "summarize", lambda *args: {"passed": True})
    results = study.assess_archive(make(), spec)
    assert all(result == {"passed": False, "pilot_reached_epsilon": False}
               for result in results.values())


@pytest.mark.parametrize("fault", ["seed", "schema_boolean", "missing_arm", "extra_arm", "stream",
                                  "points", "component_type", "density", "acceptance", "covariance"])
def test_corrupt_archive_is_rejected(archives, fault):
    make, _, spec, _ = archives
    raw = make()
    arm = raw["arms"]["correlated"]
    if fault == "seed": raw["seed"] = 41101.0
    elif fault == "schema_boolean": raw["schema_version"] = True
    elif fault == "missing_arm": del raw["arms"]["bounded_diagonal"]
    elif fault == "extra_arm": raw["arms"]["extra"] = copy.deepcopy(arm)
    elif fault == "stream": arm["spawn_key"] = [1]
    elif fault == "points": arm["production"]["points"][0][0] += 0.01
    elif fault == "component_type": arm["production"]["component_ids"][0] = float(arm["production"]["component_ids"][0])
    elif fault == "density": arm["production"]["log_q"][0] += 0.01
    elif fault == "acceptance": arm["production"]["accepted"][0] = not arm["production"]["accepted"][0]
    elif fault == "covariance": arm["proposal"]["covariances"][0][0][0] *= 1.01
    with pytest.raises(ValueError):
        study.assess_archive(raw, spec)


def _fake_control(fixture, seed, restricted):
    return {"fixture": fixture, "seed": seed,
            "kind": "support_hole_oracle" if restricted else "full_target_oracle",
            "assessment": {"passed": not restricted, "usual_checks": {"ess": True},
                           "truth_checks": {"support": not restricted}}}


@pytest.mark.parametrize("fault", ["missing", "extra", "misnamed"])
def test_report_requires_exact_archive_set(archives, tmp_path, fault):
    make, write, _, _ = archives
    first = write(make())
    if fault != "missing":
        write(make(41102))
    if fault == "extra":
        (tmp_path / "unexpected.json.gz").write_bytes(b"extra")
    elif fault == "misnamed":
        first.rename(tmp_path / "wrong.json.gz")
    with pytest.raises(ValueError, match="every planned run"):
        study.build_report(tmp_path)


def test_report_and_render_order_are_stable_and_cached_assessments_are_rebuilt(archives, tmp_path, monkeypatch):
    make, write, _, _ = archives
    for seed in (41102, 41101):
        write(make(seed))
    monkeypatch.setattr(study, "run_control", _fake_control)
    report = study.build_report(tmp_path)
    assert [(run["seed"], run["arm"]) for run in report["positive_runs"]] == [
        (seed, arm) for seed in (41101, 41102) for arm in study.ARMS]
    shuffled = copy.deepcopy(report)
    for key in ("positive_runs", "oracle_controls", "negative_controls"):
        shuffled[key].reverse()
    assert study.render(shuffled) == study.render(report)
    stale = copy.deepcopy(report)
    stale["passed"] = not report["passed"]
    stale["positive_runs"][0]["assessment"]["importance"]["normalizer_estimate"] = 9876
    assert study.assemble(stale) == report


@pytest.mark.parametrize("fault", ["identity", "runtime", "implementation_commit"])
def test_report_rejects_wrong_identity_or_mixed_provenance(archives, tmp_path, fault):
    make, write, _, _ = archives
    write(make())
    second = make(41102)
    path = write(second)
    if fault == "identity": second["seed"] = 41101
    elif fault == "runtime": second["runtime"]["numpy"] = "other"
    else: second["implementation_commit"] = "b" * 40
    path.write_bytes(gzip.compress(json.dumps(second).encode(), mtime=0))
    with pytest.raises(ValueError, match="identity|runtimes|implementation commit"):
        study.build_report(tmp_path)


@pytest.mark.parametrize("fault", ["density", "stream"])
def test_cli_resume_preserves_reports_and_archives_when_validation_fails(archives, monkeypatch, tmp_path, fault):
    make, write, _, _ = archives
    raw = make()
    if fault == "density": raw["arms"]["correlated"]["production"]["log_q"][0] += 1
    else: raw["arms"]["correlated"]["spawn_key"] = [0]
    path = write(raw)
    original = path.read_bytes()
    monkeypatch.setattr(study, "OUT_JSON", tmp_path / "report.json")
    monkeypatch.setattr(study, "OUT_MD", tmp_path / "report.md")
    study.OUT_JSON.write_text("original JSON")
    study.OUT_MD.write_text("original Markdown")
    monkeypatch.setattr(sys, "argv", ["runner", "--resume"])
    monkeypatch.setattr(study, "run_learned", lambda *args: pytest.fail("invalid archives must fail before another run"))
    with pytest.raises(ValueError, match="density|stream"):
        study.main()
    assert path.read_bytes() == original
    assert study.OUT_JSON.read_text() == "original JSON"
    assert study.OUT_MD.read_text() == "original Markdown"


def test_cli_render_only_preserves_reports_when_archive_bytes_change(archives, monkeypatch, tmp_path):
    make, write, _, _ = archives
    paths = [write(make(seed)) for seed in (41101, 41102)]
    monkeypatch.setattr(study, "run_control", _fake_control)
    report = study.build_report(tmp_path)
    monkeypatch.setattr(study, "OUT_JSON", tmp_path / "report.json")
    monkeypatch.setattr(study, "OUT_MD", tmp_path / "report.md")
    original = json.dumps(report)
    study.OUT_JSON.write_text(original)
    study.OUT_MD.write_text("original Markdown")
    paths[0].write_bytes(paths[0].read_bytes() + b"modified")
    monkeypatch.setattr(sys, "argv", ["runner", "--render-only"])
    with pytest.raises(ValueError, match="report archive hashes"):
        study.main()
    assert study.OUT_JSON.read_text() == original
    assert study.OUT_MD.read_text() == "original Markdown"
