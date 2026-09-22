"""Audit committed bounded-sampler artifacts without the biological extension.

Failures of biological sampling screens are valid scientific results. These
tests require reproducible decisions and provenance, never successful screens.
"""

import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CAL = ROOT / "analysis" / "calibration"
sys.path.insert(0, str(ROOT / "scripts"))
import abc_joint_resample as driver  # noqa: E402
import archived_numerical_sources as provenance  # noqa: E402


@pytest.fixture(autouse=True)
def no_simulator(monkeypatch):
    def forbidden():
        pytest.fail("artifact reconstruction must not load the biological extension")

    monkeypatch.setattr(driver.abc.ck, "_fc", forbidden)


@pytest.fixture(scope="module")
def archives():
    result = {}
    for seed in driver.SEEDS:
        path = CAL / f"joint-resample-run-{seed}.json.gz"
        assert path.is_file(), f"Missing required committed archive: {path}"
        result[seed] = json.loads(gzip.decompress(path.read_bytes()))
    return result


def assert_same(actual, expected, path="report"):
    """Require exact structure/decisions and tolerate last-bit numerical drift."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and actual.keys() == expected.keys(), path
        for key in expected:
            assert_same(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            assert_same(a, b, f"{path}[{i}]")
    elif isinstance(expected, float):
        assert isinstance(actual, (int, float)) and not isinstance(actual, bool), path
        assert math.isfinite(actual) and math.isfinite(expected), path
        assert actual == pytest.approx(expected, rel=1e-12, abs=1e-12), path
    else:
        assert type(actual) is type(expected) and actual == expected, path


def test_report_json_and_markdown_rebuild_without_simulator():
    rebuilt = driver.build_report(CAL)
    stored = json.loads((CAL / "joint-resample-sampling.json").read_text())
    assert_same(rebuilt, stored)
    # Render archived numbers exactly; tiny cross-platform numerical drift in
    # the independent JSON reconstruction must not alter a rounding boundary.
    markdown = (CAL / "joint-resample-sampling.md").read_text()
    assert driver.render(stored) == markdown
    assert "No pooled posterior is published." in markdown
    assert stored["stability"]["passed"] == all(stored["stability"]["checks"].values())
    # Reconstruction retains the original numerical-program identity; it does
    # not claim today's archive-verification code generated the old attempts.
    assert rebuilt["source_hashes"] == stored["source_hashes"]
    assert rebuilt["source_hashes"]["scripts/abc_joint_resample.py"] != hashlib.sha256(
        (ROOT / "scripts/abc_joint_resample.py").read_bytes()).hexdigest()


@pytest.mark.parametrize("seed", driver.SEEDS)
def test_archive_source_hashes_cover_all_frozen_numerical_inputs(archives, seed):
    archive = archives[seed]
    required = {
        "scripts/abc_joint_resample.py", "scripts/resample_move.py", "scripts/bounded_proposal.py",
        "scripts/proposal_synthetic_validation.py", "scripts/resample_move_local.py",
        "scripts/proposal_synthetic_validation_v2.py", "scripts/synthetic_proposal_study.py",
        "scripts/abc_joint_importance.py",
        "scripts/importance_sampling.py", "scripts/abc_joint_posterior.py",
        "scripts/calibrate_kill_switch.py", "scripts/calibrate_erastin.py", "scripts/ctrp_dose_support.py",
        "simulations/ferroptosis-python/src/lib.rs", "simulations/Cargo.lock",
        "simulations/rust-toolchain.toml",
    }
    required.update(str(p.relative_to(ROOT)) for p in
                    (ROOT / "simulations" / "ferroptosis-core" / "src").rglob("*.rs"))
    assert set(archive["source_hashes"]) == required
    for relative, digest in archive["source_hashes"].items():
        original = provenance.resolve_source_bytes(relative, digest, root=ROOT)
        assert hashlib.sha256(original).hexdigest() == digest, relative


def source_registry(tmp_path):
    """Small historical-byte fixture independent of production snapshots."""
    relative = "scripts/original.py"
    original = b"original numerical input\n"
    digest = hashlib.sha256(original).hexdigest()
    current = tmp_path / relative
    current.parent.mkdir(parents=True)
    current.write_bytes(b"new reader or model input\n")
    folder = tmp_path / provenance.SNAPSHOTS
    folder.mkdir(parents=True)
    (folder / f"{digest}.source").write_bytes(original)
    (folder / "manifest.json").write_text(json.dumps({
        "schema_version": 1, "sources": {relative: {digest: {"recovered_from_commit": "a" * 40}}}}))
    return relative, digest, original, folder


def test_matching_current_source_bytes_need_no_historical_registry(tmp_path):
    relative = "scripts/current.py"
    path = tmp_path / relative
    path.parent.mkdir()
    blob = b"still the original source\n"
    path.write_bytes(blob)
    digest = hashlib.sha256(blob).hexdigest()
    assert provenance.resolve_source_bytes(relative, digest, root=tmp_path) == blob
    provenance.verify_source_hashes({relative: digest}, root=tmp_path)
    assert not (tmp_path / provenance.SNAPSHOTS).exists()


def test_registered_historical_bytes_resolve_after_current_source_changes(tmp_path):
    relative, digest, original, _ = source_registry(tmp_path)
    assert provenance.resolve_source_bytes(relative, digest, root=tmp_path) == original
    provenance.verify_source_hashes({relative: digest}, root=tmp_path)


@pytest.mark.parametrize("mutation", ["missing_snapshot", "tampered_snapshot", "unknown_path",
                                     "unknown_hash", "unregistered_pair", "missing_registry", "symlink_snapshot"])
def test_historical_sources_require_registered_matching_bytes(tmp_path, mutation):
    relative, digest, _, folder = source_registry(tmp_path)
    snapshot = folder / f"{digest}.source"
    if mutation == "missing_snapshot":
        snapshot.unlink()
    elif mutation == "tampered_snapshot":
        snapshot.write_bytes(b"changed historical source\n")
    elif mutation == "unknown_path":
        relative = "scripts/another.py"
    elif mutation == "unknown_hash":
        digest = "0" * 64
    elif mutation == "unregistered_pair":
        manifest = json.loads((folder / "manifest.json").read_text())
        manifest["sources"] = {"scripts/another.py": manifest["sources"][relative]}
        (folder / "manifest.json").write_text(json.dumps(manifest))
    elif mutation == "missing_registry":
        (folder / "manifest.json").unlink()
    else:
        original = snapshot.read_bytes()
        snapshot.unlink()
        (tmp_path / "borrowed.source").write_bytes(original)
        snapshot.symlink_to(tmp_path / "borrowed.source")
    with pytest.raises(ValueError):
        provenance.verify_source_hashes({relative: digest}, root=tmp_path)


@pytest.mark.parametrize("relative", ["../original.py", "/tmp/original.py", "scripts/../original.py"])
def test_archived_source_paths_cannot_escape_the_repository(tmp_path, relative):
    with pytest.raises(ValueError, match="invalid archived source identity"):
        provenance.resolve_source_bytes(relative, "0" * 64, root=tmp_path)


@pytest.mark.parametrize("mutation", ["manifest_list", "boolean_schema", "path_list", "path_string",
                                     "record_list", "invalid_commit"])
def test_historical_registry_requires_its_declared_structure(tmp_path, mutation):
    relative, digest, _, folder = source_registry(tmp_path)
    path = folder / "manifest.json"
    manifest = json.loads(path.read_text())
    if mutation == "manifest_list":
        manifest = []
    elif mutation == "boolean_schema":
        manifest["schema_version"] = True
    elif mutation == "path_list":
        manifest["sources"][relative] = [digest]
    elif mutation == "path_string":
        manifest["sources"][relative] = digest
    elif mutation == "record_list":
        manifest["sources"][relative][digest] = []
    else:
        manifest["sources"][relative][digest]["recovered_from_commit"] = "not-a-commit"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        provenance.resolve_source_bytes(relative, digest, root=tmp_path)


@pytest.mark.parametrize("seed", driver.SEEDS)
def test_archived_target_matches_current_supported_curves_and_original_criterion(archives, seed):
    target = archives[seed]["target"]
    csv_path = CAL / "ctrpv2_ferroptosis_curves.csv"
    curves, source = driver.abc.ck.load_target_data(csv_path)
    assert target["source"] == source
    assert source["sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    assert target["priors"] == [list(p) for p in driver.abc.PRIORS]
    assert target["simulation_n"] == driver.abc.SIM_N
    assert target["simulation_seed"] == driver.abc.SIM_SEED
    assert target["reference_vector"] == driver.abc.REFERENCE_VECTOR
    assert target["tolerance_factor"] == driver.abc.TOLERANCE_FACTOR
    assert target["epsilon"] == target["reference_distance"] * target["tolerance_factor"]
    assert target["rsl3_doses_um"] == list(driver.abc.ck.DOSE_GRID_UM)
    assert target["erastin_doses_um"] == list(driver.abc.ce.DOSE_GRID_UM)
    for compound, doses, key in (
        ("ML162", target["rsl3_doses_um"], "empirical_rsl3"),
        ("ERASTIN", target["erastin_doses_um"], "empirical_erastin"),
        ("ML210", target["rsl3_doses_um"], "empirical_heldout"),
    ):
        empirical, support = driver.abc.ck.empirical_target(curves[compound], doses)
        np.testing.assert_allclose(target[key], empirical, rtol=1e-12, atol=1e-12)
        assert target["support"][compound] == support
    previous = json.loads((CAL / "joint-importance-sampling.json").read_text())
    assert target == previous["target"]
    original = json.loads((CAL / "joint-posterior.json").read_text())
    assert original["target_source"] == source and original["target_support"] == target["support"]
    assert original["reference_distance"] == pytest.approx(target["reference_distance"], rel=0, abs=5e-5 + 1e-12)
    assert original["epsilon_joint_distance"] == pytest.approx(target["epsilon"], rel=0, abs=5e-5 + 1e-12)


@pytest.mark.parametrize("seed", driver.SEEDS)
def test_pilot_streams_budgets_costs_and_final_kernel_reconstruct(archives, seed):
    archive = archives[seed]
    assert archive["seed"] == seed and archive["schema_version"] == 2
    assert archive["plan"] == driver.PLAN and archive["gates"] == driver.GATES
    driver.validate_pilot(archive)
    plan, pilot = archive["plan"]["pilot"], archive["pilot"]
    n = plan["particles_per_island"]
    d = len(driver.NAMES)
    island_streams = np.random.SeedSequence(seed).spawn(2)[0].spawn(plan["islands"])
    full_cost = len(archive["target"]["rsl3_doses_um"]) + len(archive["target"]["erastin_doses_um"])
    early_cost = len(archive["target"]["erastin_doses_um"])
    calculated_calls = 0
    final_points = []
    for island, stream in zip(pilot["islands"], island_streams):
        np.testing.assert_array_equal(island["initial_points"], np.random.default_rng(stream).uniform(size=(n, d)))
        assert island["spawn_key"] == list(stream.spawn_key)
        assert island["initial_simulator_dose_calls"] == n * full_cost
        calculated_calls += n * full_cost
        assert len(island["history"]) == plan["levels"]
        thresholds = np.array([h["threshold"] for h in island["history"]])
        assert np.all(np.diff(thresholds) <= 0) and np.all(thresholds >= archive["target"]["epsilon"])
        for row in island["history"]:
            assert row["attempts"] == n * plan["moves_per_level"]
            assert row["attempts"] == row["evaluated"] + row["mh_rejected_without_evaluation"]
            for key in ("local_attempts", "global_attempts", "local_accepted_moves"):
                assert type(row[key]) is int and row[key] >= 0
            assert row["local_attempts"] + row["global_attempts"] == row["attempts"]
            assert row["mh_rejected_without_evaluation"] <= row["global_attempts"]
            global_evaluated = row["evaluated"] - row["local_attempts"]
            assert 0 <= global_evaluated <= row["global_attempts"]
            assert row["local_accepted_moves"] <= min(row["local_attempts"], row["accepted_moves"])
            global_accepted = row["accepted_moves"] - row["local_accepted_moves"]
            assert 0 <= global_accepted <= global_evaluated
            full_evaluations = row["evaluated"] - row["early_rejected"]
            assert 0 <= row["accepted_moves"] <= full_evaluations
            assert row["simulator_dose_calls"] == early_cost * row["early_rejected"] + full_cost * full_evaluations
            calculated_calls += row["simulator_dose_calls"]
        final_points.extend(island["final_points"])
    assert pilot["attempts"] == plan["islands"] * n * (1 + plan["levels"] * plan["moves_per_level"])
    assert pilot["simulator_dose_calls"] == calculated_calls
    assert pilot["final_points"] == final_points
    assert_same(driver.fit_proposal(final_points, plan).to_dict(), archive["proposal"], "proposal")


@pytest.mark.parametrize("seed", driver.SEEDS)
def test_production_replay_and_every_decision_reconstruct_from_retained_curves(archives, seed):
    archive = archives[seed]
    driver.replay_production(archive)
    production, target = archive["production"], archive["target"]
    _, stream = np.random.SeedSequence(seed).spawn(2)
    rng = np.random.default_rng(stream)
    assert archive["rng"]["bit_generator"] == type(rng.bit_generator).__name__
    proposal = driver.BoundedGaussianMixture.from_dict(archive["proposal"])
    points, ids = proposal.sample(rng, archive["plan"]["production_attempts"])
    np.testing.assert_array_equal(production["component_ids"], ids)
    np.testing.assert_allclose(production["unit_points"], points, rtol=1e-12, atol=1e-12)
    log_q = proposal.log_density(np.asarray(production["unit_points"]))
    np.testing.assert_allclose(production["log_q"], log_q, rtol=1e-12, atol=1e-12)
    assert np.all((points >= 0) & (points <= 1))
    assert len(production["records"]) == len(points)
    accepted = []
    actual_calls = 0
    for record in production["records"]:
        assert type(record["accepted"]) is bool
        driver._curve(record["erastin"], len(target["erastin_doses_um"]), "stored erastin")
        erastin_error = driver.abc.ck.rmse(record["erastin"], target["empirical_erastin"])
        if record["status"] == "early_rejected":
            assert record["rsl3"] is None and record["distance"] is None
            assert erastin_error > target["epsilon"] and not record["accepted"]
            assert record["distance_lower_bound"] == pytest.approx(erastin_error, rel=1e-12, abs=1e-12)
            assert record["simulator_dose_calls"] == len(target["erastin_doses_um"])
        else:
            assert record["status"] == "complete"
            driver._curve(record["rsl3"], len(target["rsl3_doses_um"]), "stored RSL3")
            distance = erastin_error + driver.abc.ck.rmse(record["rsl3"], target["empirical_rsl3"])
            assert record["distance"] == pytest.approx(distance, rel=1e-12, abs=1e-12)
            assert record["distance_lower_bound"] == record["distance"]
            assert record["accepted"] == (distance <= target["epsilon"])
            assert record["simulator_dose_calls"] == len(target["rsl3_doses_um"]) + len(target["erastin_doses_um"])
        accepted.append(record["accepted"])
        actual_calls += record["simulator_dose_calls"]
    summary = driver.summarize_run(archive)
    raw_weights = np.where(accepted, np.exp(-np.asarray(production["log_q"])), 0.0)
    assert summary["importance"]["n_accepted"] == sum(accepted)
    assert summary["importance"]["normalizer_estimate"] == pytest.approx(raw_weights.mean(), rel=1e-12, abs=1e-12)
    expected_ess = raw_weights.sum()**2 / np.dot(raw_weights, raw_weights) if any(accepted) else 0
    assert summary["importance"]["ess"] == pytest.approx(expected_ess, rel=1e-12, abs=1e-12)
    assert summary["simulator_dose_calls_including_pilot_reference"] == (
        actual_calls + archive["pilot"]["simulator_dose_calls"] + target["reference_simulator_dose_calls"])


def test_report_links_all_archives_prerequisite_and_historical_baselines(archives):
    report = json.loads((CAL / "joint-resample-sampling.json").read_text())
    expected_hashes = {}
    study_path = CAL / "proposal-synthetic-validation-v2.json"
    study = json.loads(study_path.read_text())
    prerequisite = report["synthetic_validation"]
    assert prerequisite["artifact"] == study_path.name
    assert study["revision_id"] == "local-move-v2"
    assert study["pilot_plan"] == driver.PLAN["pilot"]
    assert prerequisite["sha256"] == hashlib.sha256(study_path.read_bytes()).hexdigest()
    # Passing known-target validation is an execution prerequisite; biological
    # adequacy is deliberately not required anywhere in this module.
    assert prerequisite["passed"] is True and study["passed"] is True
    assert prerequisite["source_hashes"] == study["source_hashes"]
    for relative, digest in prerequisite["source_hashes"].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest
    runtimes = []
    for seed, archive in archives.items():
        path = CAL / f"joint-resample-run-{seed}.json.gz"
        expected_hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        assert archive["target"] == report["target"]
        assert archive["source_hashes"] == report["source_hashes"]
        assert archive["synthetic_validation"] == prerequisite
        runtimes.append(archive["runtime"])
    assert report["archive_sha256"] == expected_hashes
    assert report["runtime"] == runtimes[0]
    assert runtimes[0]["extension_binaries"]
    assert all(re.fullmatch(r"[0-9a-f]{64}", digest)
               for digest in runtimes[0]["extension_binaries"].values())
    for runtime in runtimes[1:]:
        for key in ("extension_binaries", "python", "numpy", "scipy"):
            assert runtime[key] == runtimes[0][key]
    assert report["baselines"]["uniform"]["sha256"] == hashlib.sha256((CAL / "joint-posterior.json").read_bytes()).hexdigest()
    assert report["baselines"]["first_importance"]["sha256"] == hashlib.sha256(
        (CAL / "joint-importance-sampling.json").read_bytes()).hexdigest()
