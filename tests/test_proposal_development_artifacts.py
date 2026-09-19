"""Preserve every development attempt without treating selection as validation."""

import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "analysis/calibration/proposal-development-20260919"
sys.path.insert(0, str(ROOT / "scripts"))
from bounded_proposal import BoundedGaussianMixture
import proposal_synthetic_validation as historical


SEEDS = (2026092301, 2026092302, 2026092303)
CONFIGURATIONS = (
    "independent_k32", "independent_k8", "mixed_k32", "mixed_k8",
    "mixed_gap32", "staged_gap32",
)
ARCHIVES = ("development_results.json.gz", "followup_results.json.gz")
SCRIPTS = (
    "run_development.py.txt", "run_followup.py.txt",
    "check_development.py.txt", "verify_candidate.py.txt",
)
ORIGINAL_SOURCES = (
    "scripts/proposal_synthetic_validation.py", "scripts/resample_move.py",
    "scripts/bounded_proposal.py", "scripts/importance_sampling.py",
)


@pytest.fixture(scope="module")
def metadata():
    return json.loads((FOLDER / "metadata.json").read_text())


@pytest.fixture(scope="module")
def studies():
    return {name: json.loads(gzip.decompress((FOLDER / name).read_bytes()))
            for name in ARCHIVES}


@pytest.fixture(scope="module")
def records(studies):
    return [record for name in ARCHIVES for record in studies[name]["results"]]


def assert_nested_close(actual, expected):
    """Decisions/structure are exact; only floating arithmetic gets tolerance."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict)
        assert actual.keys() == expected.keys()
        for key in expected:
            assert_nested_close(actual[key], expected[key])
    elif isinstance(expected, list):
        assert isinstance(actual, list)
        assert len(actual) == len(expected)
        for observed, stored in zip(actual, expected):
            assert_nested_close(observed, stored)
    elif isinstance(expected, float):
        assert not isinstance(actual, bool)
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    else:
        assert type(actual) is type(expected)
        assert actual == expected


def test_all_eighteen_declared_attempts_remain_distinct_development_records(metadata, studies, records):
    assert metadata["status"] == "development_only"
    assert metadata["selected_configuration"] == "mixed_gap32"
    assert metadata["development_seeds"] == list(SEEDS)
    assert [(r["development_config"], r["seed"]) for r in records] == [
        (configuration, seed) for configuration in CONFIGURATIONS for seed in SEEDS
    ]
    assert [len(studies[name]["results"]) for name in ARCHIVES] == [12, 6]
    assert [d["order"] for d in metadata["declarations"]] == [1, 2]
    assert [d["configurations"] for d in metadata["declarations"]] == [
        list(CONFIGURATIONS[:4]), list(CONFIGURATIONS[4:])
    ]
    assert [d["n_records"] for d in metadata["declarations"]] == [12, 6]
    for record in records:
        assert record["development_only"] is True
        assert record["fixture"] == "separated_boxes"
        assert record["kind"] == "learned_proposal"
        assert record["production_attempts"] == 8192
        assert record["pilot"]["attempts"] == 29184
        assert record["pilot"]["simulator_dose_calls"] == 0
        assert record["rng"] == {"pilot_spawn_key": [0], "production_spawn_key": [1]}


@pytest.mark.parametrize("name", ARCHIVES)
def test_archive_compressed_and_original_bytes_match_recorded_hashes(name, metadata, studies):
    raw = (FOLDER / name).read_bytes()
    info = metadata["archives"][name]
    assert hashlib.sha256(raw).hexdigest() == info["sha256"]
    assert hashlib.sha256(gzip.decompress(raw)).hexdigest() == info["uncompressed_sha256"]
    assert raw[:3] == b"\x1f\x8b\x08"
    assert raw[3] & 0x08 == 0  # No original-filename header.
    assert raw[4:8] == b"\x00\x00\x00\x00"
    assert info["gzip_mtime"] == 0
    assert info["gzip_filename"] == ""
    assert len(studies[name]["results"]) == info["n_records"]
    assert Path(info["original_path"]).is_absolute()
    assert studies[name]["source_hashes"] == metadata["original_source_sha256"]


@pytest.mark.parametrize("name", SCRIPTS)
def test_scratch_code_bytes_are_preserved_as_nonexecuted_text(name, metadata):
    info = metadata["scratch_scripts"][name]
    assert hashlib.sha256((FOLDER / name).read_bytes()).hexdigest() == info["sha256"]
    assert Path(info["original_path"]).is_absolute()
    assert Path(info["original_path"]).name + ".txt" == name


@pytest.mark.parametrize("name", ORIGINAL_SOURCES)
def test_original_source_hashes_identify_unchanged_historical_modules(name, metadata):
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == metadata["original_source_sha256"][name]


@pytest.mark.parametrize("configuration,seed", [
    (configuration, seed) for configuration in CONFIGURATIONS for seed in SEEDS
])
def test_every_production_summary_replays_from_its_frozen_proposal(configuration, seed, records):
    record = next(r for r in records if (r["development_config"], r["seed"]) == (configuration, seed))
    proposal = BoundedGaussianMixture.from_dict(record["proposal"])
    stream = np.random.SeedSequence(seed).spawn(2)[1]
    points, components = proposal.sample(np.random.default_rng(stream), record["production_attempts"])
    assert points.shape == (8192, 7)
    assert np.all((points >= 0) & (points <= 1))
    assert np.bincount(components, minlength=len(proposal.means) + 1).tolist() == record["production_component_counts"]
    rebuilt = historical.summarize(points, proposal.log_density(points), record["fixture"])
    rebuilt["pilot_reached_epsilon"] = record["pilot"]["all_islands_reached_epsilon"]
    rebuilt["passed"] = rebuilt["passed"] and rebuilt["pilot_reached_epsilon"]
    assert_nested_close(rebuilt, record["assessment"])


def test_past_candidate_verification_is_historical_evidence_not_a_current_source_pin(metadata):
    info = metadata["candidate_verification"]
    raw = (FOLDER / info["path"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == info["sha256"]
    verification = json.loads(raw)
    assert verification["development_only"] is True
    assert verification["candidate_source_path"] == "scripts/resample_move_local.py"
    assert len(bytes.fromhex(verification["candidate_source_sha256"])) == 32
    assert [check["seed"] for check in verification["checks"]] == list(SEEDS)
    expected_fields = {
        "seed", "pilot_exact", "island_history_original_fields_exact", "proposal_exact",
        "production_points_ids_density_exact", "assessment_exact",
    }
    for check in verification["checks"]:
        assert set(check) == expected_fields
        assert all(check[field] is True for field in expected_fields - {"seed"})


def test_readme_table_preserves_successes_and_failures_from_all_eighteen_records(records):
    readme = (FOLDER / "README.md").read_text()
    rows = [line for line in readme.splitlines() if any(line.startswith(f"| {name} |") for name in CONFIGURATIONS)]
    expected = []
    for record in records:
        result = record["assessment"]
        failed = [name for name, passed in {**result["usual_checks"], **result["truth_checks"]}.items() if not passed]
        if not result["pilot_reached_epsilon"]:
            failed.append("pilot_reached_epsilon")
        expected.append(
            f"| {record['development_config']} | {str(record['seed'])[-4:]} | "
            f"{'yes' if result['passed'] else 'no'} | {result['importance']['ess']:.2f} | "
            f"{result['mode_masses']['A']:.4f} | {result['relative_mass_error']:+.2%} | "
            f"{', '.join(failed) or 'none'} |"
        )
    assert rows == expected
