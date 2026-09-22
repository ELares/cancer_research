"""Rebuild committed importance-sampling reports and verify their provenance.

No simulator or local extension binary is required. Archived binary hashes are
checked for consistency, while original source and CSV hashes are recomputed.
Historical source identities resolve to exact registered snapshots after edits.
These checks verify reproducibility and integrity, not that adequacy screens pass
or that the simulation describes experimental outcomes.
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
import abc_joint_importance as driver  # noqa: E402
import archived_numerical_sources as provenance  # noqa: E402


@pytest.fixture(autouse=True)
def no_simulator(monkeypatch):
    def forbidden():
        pytest.fail("committed artifact checks must not load the compiled extension")

    monkeypatch.setattr(driver.abc.ck, "_fc", forbidden)


@pytest.fixture(scope="module")
def archive_by_seed():
    cache = {}

    def load(seed):
        if seed not in cache:
            path = CAL / f"joint-importance-run-{seed}.json.gz"
            assert path.is_file(), f"Missing committed importance archive: {path}"
            # Read the archive independently of the driver's archive loader.
            cache[seed] = json.loads(gzip.decompress(path.read_bytes()))
        return cache[seed]

    return load


def assert_same_report(actual, expected, path="report"):
    """Preserve exact structure and decisions; allow last-bit math-library drift."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict), path
        assert actual.keys() == expected.keys(), path
        for key in expected:
            assert_same_report(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for i, (left, right) in enumerate(zip(actual, expected)):
            assert_same_report(left, right, f"{path}[{i}]")
    elif isinstance(expected, float):
        assert isinstance(actual, (int, float)) and not isinstance(actual, bool), path
        assert math.isfinite(actual) and math.isfinite(expected), path
        assert actual == pytest.approx(expected, rel=1e-12, abs=1e-12), path
    else:
        assert type(actual) is type(expected) and actual == expected, path


def test_committed_json_and_markdown_rebuild_from_all_raw_archives():
    rebuilt = driver.build_report(CAL)
    committed = json.loads((CAL / "joint-importance-sampling.json").read_text())
    assert_same_report(rebuilt, committed)
    markdown = (CAL / "joint-importance-sampling.md").read_text()
    assert driver.render(rebuilt) == markdown
    assert driver.render(committed) == markdown


@pytest.mark.parametrize("seed", driver.SEEDS)
def test_archived_numerical_source_hashes_match_available_original_bytes(archive_by_seed, seed):
    archive = archive_by_seed(seed)
    required = {
        "scripts/abc_joint_importance.py", "scripts/importance_sampling.py",
        "scripts/abc_joint_posterior.py", "scripts/calibrate_kill_switch.py",
        "scripts/calibrate_erastin.py", "scripts/ctrp_dose_support.py",
        "simulations/ferroptosis-python/src/lib.rs", "simulations/Cargo.lock",
        "simulations/rust-toolchain.toml",
    }
    required.update(str(path.relative_to(ROOT)) for path in
                    (ROOT / "simulations" / "ferroptosis-core" / "src").rglob("*.rs"))
    assert set(archive["source_hashes"]) == required
    for relative, digest in archive["source_hashes"].items():
        original = provenance.resolve_source_bytes(relative, digest, root=ROOT)
        assert hashlib.sha256(original).hexdigest() == digest, relative


@pytest.mark.parametrize("seed", driver.SEEDS)
def test_archived_targets_reproduce_supported_csv_cohorts_and_baseline(archive_by_seed, seed):
    target = archive_by_seed(seed)["target"]
    csv_path = CAL / "ctrpv2_ferroptosis_curves.csv"
    curves, source = driver.abc.ck.load_target_data(csv_path)
    assert target["source"] == source
    assert target["source"]["sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
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
        np.testing.assert_allclose(target[key], empirical, rtol=1e-13, atol=1e-15)
        assert target["support"][compound] == support
    baseline = json.loads((CAL / "joint-posterior.json").read_text())
    assert baseline["target_source"] == target["source"]
    assert baseline["target_support"] == target["support"]
    assert baseline["tolerance_factor"] == target["tolerance_factor"]
    for old_key, new_key in (("reference_distance", "reference_distance"),
                             ("epsilon_joint_distance", "epsilon")):
        assert baseline[old_key] == pytest.approx(target[new_key], rel=0, abs=5e-5 + 1e-12)
    # The original rejection artifact stores targets rounded to four decimals.
    for old, new in (
        (baseline["curves"]["empirical_rsl3_ml162"], target["empirical_rsl3"]),
        (baseline["curves"]["empirical_erastin"], target["empirical_erastin"]),
        (baseline["heldout_posterior_predictive"]["empirical"], target["empirical_heldout"]),
    ):
        np.testing.assert_allclose(old, new, rtol=0, atol=5e-5 + 1e-12)


@pytest.mark.parametrize("seed", driver.SEEDS)
def test_frozen_production_points_and_component_ids_reproduce_the_archived_seed(archive_by_seed, seed):
    archive = archive_by_seed(seed)
    assert archive["seed"] == seed
    assert archive["schema_version"] == 1
    assert archive["plan"] == driver.PLAN and archive["gates"] == driver.GATES
    pilot_stream, production_stream = np.random.SeedSequence(seed).spawn(2)
    assert archive["rng"]["pilot_spawn_key"] == list(pilot_stream.spawn_key)
    assert archive["rng"]["production_spawn_key"] == list(production_stream.spawn_key)
    rng = np.random.default_rng(production_stream)
    assert archive["rng"]["bit_generator"] == type(rng.bit_generator).__name__
    proposal = driver.GaussianMixture.from_dict(archive["proposal"])
    points, ids = proposal.sample(rng, archive["plan"]["production_attempts"])
    production = archive["production"]
    np.testing.assert_array_equal(production["component_ids"], ids)
    np.testing.assert_allclose(production["unit_points"], points, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(production["log_q"], proposal.log_density(points), rtol=1e-12, atol=1e-12)
    assert len(production["records"]) == archive["plan"]["production_attempts"]
    assert len(archive["pilot"]) == archive["plan"]["pilot_rounds"]
    for step, history in enumerate(archive["pilot"], start=1):
        assert history["round"] == step
        assert history["attempts"] == archive["plan"]["pilot_attempts"]


def test_report_hashes_link_every_archive_and_the_unchanged_baseline(archive_by_seed):
    report = json.loads((CAL / "joint-importance-sampling.json").read_text())
    expected = {}
    runtimes = []
    for seed in driver.SEEDS:
        path = CAL / f"joint-importance-run-{seed}.json.gz"
        expected[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        archive = archive_by_seed(seed)
        assert archive["target"] == report["target"]
        assert archive["source_hashes"] == report["source_hashes"]
        runtimes.append(archive["runtime"])
    assert report["archive_sha256"] == expected
    baseline_path = CAL / "joint-posterior.json"
    assert report["baseline"]["sha256"] == hashlib.sha256(baseline_path.read_bytes()).hexdigest()
    assert report["runtime"] == runtimes[0]
    assert runtimes[0]["extension_binaries"], "The run must identify the simulator binary."
    for digest in runtimes[0]["extension_binaries"].values():
        assert re.fullmatch(r"[0-9a-f]{64}", digest)
    # Recorded runtime/build identity is portable metadata; no local .so is read.
    for runtime in runtimes[1:]:
        for key in ("extension_binaries", "python", "numpy", "scipy"):
            assert runtime[key] == runtimes[0][key]
