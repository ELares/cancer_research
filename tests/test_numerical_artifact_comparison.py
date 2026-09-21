"""Platform float tolerance must not turn frozen experiments into fuzzy data."""
from copy import deepcopy
import json
import math
import sys
from types import SimpleNamespace

import pytest

import test_artifact_freshness as freshness


HISTORICAL_NAMES = ("proposal_synthetic_validation", "proposal_synthetic_validation_v2")
COVERAGE_NAME = "proposal_coverage_challenges"
NAMES = (*HISTORICAL_NAMES, COVERAGE_NAME)


def _artifact():
    assessment = {
        "importance": {
            "ess": 1000.0, "n_accepted": 1333, "n_attempts": 8192,
            "normalizer_estimate": 0.0002444950974831571,
        },
        "mode_counts": {"A": 496, "B": 837},
        "moments": {"mean_x1": 0.4965},
        "quantiles": [0.125, 0.625],
        "passed": True,
        "truth_checks": {"normalizing_mass": True},
        "unavailable": None,
    }
    run = {
        "seed": 17,
        "assessment": assessment,
        "wall_seconds": 0.25,
        "pilot": {"initial_particles": [[0.1, 0.2]], "threshold": 0.5},
        "proposal": {"means": [[0.125, 0.75]], "scales": [[0.04, 0.05]],
                     "uniform_weight": 0.2},
    }
    second, negative = deepcopy(run), deepcopy(run)
    second["seed"] = 18
    negative["seed"] = 19
    negative["assessment"]["passed"] = False
    return {
        "schema_version": 2,
        "positive_runs": [run, second],
        "negative_controls": [negative],
        "passed": True,
        "gates": {"minimum_ess": 200.0},
        "source_hashes": {"scripts/example.py": "a" * 64},
        "truth": {"mass": 0.25},
        # Lookalike paths must not acquire assessment tolerance.
        "assessment": {"moment": 0.25},
        "metadata": {"positive_runs": [{"assessment": {"moment": 0.25}}]},
    }


def _dump(value):
    # Independent spelling of the three writers' public JSON format.
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _replace(value, path, replacement):
    for part in path[:-1]:
        value = value[part]
    value[path[-1]] = replacement


def test_only_the_three_declared_synthetic_reports_receive_numeric_tolerance():
    assert freshness.NUMERICAL_REASSEMBLY == frozenset(NAMES)


def _coverage_artifact():
    result = _artifact()
    oracle = deepcopy(result["positive_runs"][0])
    oracle.update(kind="full_target_oracle", retained_mass_fraction=1.0)
    oracle["assessment"]["region_counts"] = {"quadrant_0": 2048, "quadrant_1": 2048}
    result.update(
        oracle_controls=[oracle],
        archive_sha256={"test-run.json.gz": "c" * 64},
        specification={"source_hashes": {"scripts/example.py": "a" * 64},
                       "production_attempts": 8192, "epsilon": 1.0},
        runtime={"numpy": "test-runtime"},
        checks={"all_nine_full_target_oracles_pass": True},
    )
    result["metadata"]["oracle_controls"] = [{"assessment": {"moment": 0.25}}]
    return result


@pytest.mark.parametrize("name", [
    "abc_joint_resample", "proposal_synthetic_validation_v3",
    "proposal_synthetic_validation_v2_copy",
])
def test_other_generators_remain_byte_exact(name):
    original = _artifact()
    changed = deepcopy(original)
    changed["positive_runs"][0]["assessment"]["importance"]["ess"] += 5e-10
    freshness._assert_reassembled_json(name, _dump(original), _dump(original))
    with pytest.raises(AssertionError):
        freshness._assert_reassembled_json(name, _dump(changed), _dump(original))


@pytest.mark.parametrize("branch", ["positive_runs", "negative_controls"])
@pytest.mark.parametrize("old,new", [
    (0.0, 5e-15),
    (0.25, math.nextafter(0.25, math.inf)),
    (1000.0, 1000.0 + 5e-10),
    (-1000.0, -1000.0 - 5e-10),
    (0.0002444950974831571, 0.0002444950974831572),
])
def test_only_small_finite_derived_float_changes_are_allowed(branch, old, new):
    original = _artifact()
    path = (branch, 0, "assessment", "importance", "ess")
    _replace(original, path, old)
    changed = deepcopy(original)
    _replace(changed, path, new)
    for name in NAMES:
        freshness._assert_reassembled_json(name, _dump(changed), _dump(original))
        freshness._assert_reassembled_json(name, _dump(original), _dump(changed))


@pytest.mark.parametrize("old,new", [
    (0.0, 2e-14), (0.25, 0.25 + 5e-13),
    (1000.0, 1000.0 + 2e-9), (-1000.0, -1000.0 - 2e-9),
])
def test_absolute_and_relative_tolerances_have_material_boundaries(old, new):
    original = _artifact()
    path = ("positive_runs", 0, "assessment", "importance", "ess")
    _replace(original, path, old)
    changed = deepcopy(original)
    _replace(changed, path, new)
    with pytest.raises(AssertionError, match="derived float differs"):
        freshness._assert_reassembled_json(NAMES[0], _dump(changed), _dump(original))


@pytest.mark.parametrize("path,replacement", [
    (("positive_runs", 0, "assessment", "importance", "ess"), 1000.01),
    (("positive_runs", 0, "assessment", "importance", "n_accepted"), 1334),
    (("negative_controls", 0, "assessment", "mode_counts", "A"), 497),
    (("positive_runs", 0, "assessment", "importance", "n_accepted"), 1333.0),
    (("positive_runs", 0, "assessment", "importance", "ess"), 1000),
    (("positive_runs", 0, "assessment", "passed"), 1),
    (("negative_controls", 0, "assessment", "passed"), 0),
    (("positive_runs", 0, "assessment", "passed"), False),
    (("positive_runs", 0, "assessment", "truth_checks", "normalizing_mass"), False),
    (("positive_runs", 0, "assessment", "unavailable"), "null"),
    (("source_hashes", "scripts/example.py"), "b" * 64),
    (("positive_runs", 0, "proposal", "scales", 0, 0), math.nextafter(0.04, math.inf)),
    (("positive_runs", 0, "proposal", "means", 0, 0), math.nextafter(0.125, math.inf)),
    (("positive_runs", 0, "proposal", "uniform_weight"), math.nextafter(0.2, math.inf)),
    (("positive_runs", 0, "pilot", "initial_particles", 0, 0), math.nextafter(0.1, math.inf)),
    (("positive_runs", 0, "pilot", "threshold"), math.nextafter(0.5, math.inf)),
    (("positive_runs", 0, "wall_seconds"), math.nextafter(0.25, math.inf)),
    (("truth", "mass"), math.nextafter(0.25, math.inf)),
    (("gates", "minimum_ess"), math.nextafter(200.0, math.inf)),
    (("assessment", "moment"), math.nextafter(0.25, math.inf)),
    (("metadata", "positive_runs", 0, "assessment", "moment"), math.nextafter(0.25, math.inf)),
])
def test_material_results_types_decisions_provenance_and_raw_values_are_strict(path, replacement):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, path, replacement)
    for name in NAMES:
        with pytest.raises(AssertionError):
            freshness._assert_reassembled_json(name, _dump(changed), _dump(original))


@pytest.mark.parametrize("old,new", [
    (0.0, 5e-15),
    (0.25, math.nextafter(0.25, math.inf)),
    (1000.0, 1000.0 + 5e-10),
    (-1000.0, -1000.0 - 5e-10),
])
def test_only_the_new_study_allows_small_oracle_assessment_float_changes(old, new):
    original = _coverage_artifact()
    path = ("oracle_controls", 0, "assessment", "importance", "ess")
    _replace(original, path, old)
    changed = deepcopy(original)
    _replace(changed, path, new)
    freshness._assert_reassembled_json(COVERAGE_NAME, _dump(changed), _dump(original))
    freshness._assert_reassembled_json(COVERAGE_NAME, _dump(original), _dump(changed))
    for name in HISTORICAL_NAMES:
        with pytest.raises(AssertionError, match="exact JSON value differs"):
            freshness._assert_reassembled_json(name, _dump(changed), _dump(original))


@pytest.mark.parametrize("path,replacement", [
    (("oracle_controls", 0, "assessment", "importance", "ess"), 1000.0 + 2e-9),
    (("oracle_controls", 0, "assessment", "importance", "n_accepted"), 1334),
    (("oracle_controls", 0, "assessment", "importance", "n_accepted"), 1333.0),
    (("oracle_controls", 0, "assessment", "importance", "ess"), 1000),
    (("oracle_controls", 0, "assessment", "region_counts", "quadrant_0"), 2049),
    (("oracle_controls", 0, "assessment", "passed"), False),
    (("oracle_controls", 0, "assessment", "passed"), 1),
    (("oracle_controls", 0, "assessment", "truth_checks", "normalizing_mass"), False),
    (("oracle_controls", 0, "seed"), 18),
    (("oracle_controls", 0, "kind"), "support_hole_oracle"),
    (("oracle_controls", 0, "retained_mass_fraction"), math.nextafter(1.0, 0.0)),
    (("oracle_controls", 0, "wall_seconds"), math.nextafter(0.25, math.inf)),
    (("oracle_controls", 0, "pilot", "threshold"), math.nextafter(0.5, math.inf)),
    (("oracle_controls", 0, "proposal", "means", 0, 0), math.nextafter(0.125, math.inf)),
    (("archive_sha256", "test-run.json.gz"), "d" * 64),
    (("specification", "source_hashes", "scripts/example.py"), "b" * 64),
    (("specification", "production_attempts"), 8193),
    (("specification", "epsilon"), math.nextafter(1.0, math.inf)),
    (("runtime", "numpy"), "different-runtime"),
    (("checks", "all_nine_full_target_oracles_pass"), False),
    (("metadata", "oracle_controls", 0, "assessment", "moment"), math.nextafter(0.25, math.inf)),
])
def test_new_oracle_scope_keeps_counts_flags_raw_values_and_provenance_strict(path, replacement):
    original = _coverage_artifact()
    changed = deepcopy(original)
    _replace(changed, path, replacement)
    with pytest.raises(AssertionError):
        freshness._assert_reassembled_json(COVERAGE_NAME, _dump(changed), _dump(original))


def test_large_oracle_count_cannot_receive_relative_float_tolerance():
    original = _coverage_artifact()
    path = ("oracle_controls", 0, "assessment", "region_counts", "quadrant_0")
    _replace(original, path, 10**15)
    changed = deepcopy(original)
    _replace(changed, path, 10**15 + 1)
    with pytest.raises(AssertionError, match="exact JSON value differs"):
        freshness._assert_reassembled_json(COVERAGE_NAME, _dump(changed), _dump(original))


def test_large_integer_counts_cannot_slip_through_relative_float_tolerance():
    original = _artifact()
    path = ("positive_runs", 0, "assessment", "importance", "n_attempts")
    _replace(original, path, 10**15)
    changed = deepcopy(original)
    _replace(changed, path, 10**15 + 1)
    with pytest.raises(AssertionError):
        freshness._assert_reassembled_json(NAMES[0], _dump(changed), _dump(original))


@pytest.mark.parametrize("mutation", [
    "missing_key", "extra_key", "list_length", "list_order", "run_order", "container_type",
])
def test_structure_and_list_order_remain_exact(mutation):
    original = _artifact()
    changed = deepcopy(original)
    assessment = changed["positive_runs"][0]["assessment"]
    if mutation == "missing_key":
        del assessment["passed"]
    elif mutation == "extra_key":
        assessment["extra"] = 0.0
    elif mutation == "list_length":
        assessment["quantiles"].append(0.875)
    elif mutation == "list_order":
        assessment["quantiles"].reverse()
    elif mutation == "run_order":
        changed["positive_runs"].reverse()
    elif mutation == "container_type":
        assessment["quantiles"] = {"0": 0.125, "1": 0.625}
    with pytest.raises(AssertionError):
        freshness._assert_reassembled_json(NAMES[0], _dump(changed), _dump(original))


def test_raw_signed_zero_is_not_equal_as_serialized_json():
    original = _artifact()
    path = ("positive_runs", 0, "proposal", "means", 0, 0)
    _replace(original, path, 0.0)
    changed = deepcopy(original)
    _replace(changed, path, -0.0)
    with pytest.raises(AssertionError):
        freshness._assert_reassembled_json(NAMES[0], _dump(changed), _dump(original))


@pytest.mark.parametrize("side", ["produced", "committed"])
@pytest.mark.parametrize("formatting", ["indent", "key_order", "missing_newline", "extra_newline"])
def test_tolerant_reports_still_require_the_writers_exact_format(side, formatting):
    original = _artifact()
    canonical = _dump(original)
    alternatives = {
        "indent": json.dumps(original, indent=4, sort_keys=True) + "\n",
        "key_order": json.dumps(original, indent=2, sort_keys=False) + "\n",
        "missing_newline": canonical.rstrip("\n"),
        "extra_newline": canonical + "\n",
    }
    args = [alternatives[formatting], canonical] if side == "produced" else [canonical, alternatives[formatting]]
    with pytest.raises(AssertionError, match="JSON formatting differs"):
        freshness._assert_reassembled_json(NAMES[0], *args)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
@pytest.mark.parametrize("side", ["produced", "committed"])
def test_nonfinite_numbers_are_rejected_on_either_side(bad, side):
    original = _artifact()
    changed = deepcopy(original)
    changed["positive_runs"][0]["assessment"]["importance"]["ess"] = bad
    invalid = json.dumps(changed, indent=2, sort_keys=True) + "\n"
    args = [invalid, _dump(original)] if side == "produced" else [_dump(original), invalid]
    with pytest.raises(AssertionError, match="invalid .* JSON numbers"):
        freshness._assert_reassembled_json(NAMES[0], *args)


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("delta", [5e-10, 0.1])
def test_the_freshness_gate_uses_the_scoped_comparator(name, delta, tmp_path, monkeypatch):
    original = _artifact()
    changed = deepcopy(original)
    changed["positive_runs"][0]["assessment"]["importance"]["ess"] += delta
    committed_path = tmp_path / "artifact.json"
    committed_path.write_text(_dump(original))
    monkeypatch.setitem(sys.modules, name, SimpleNamespace(OUT_JSON=committed_path))
    monkeypatch.setattr(freshness, "_reproduce", lambda mod, script: deepcopy(changed))
    if delta < 1e-9:
        freshness.test_the_committed_json_is_what_the_generator_writes(name, tmp_path, monkeypatch)
    else:
        with pytest.raises(AssertionError, match="derived float differs"):
            freshness.test_the_committed_json_is_what_the_generator_writes(name, tmp_path, monkeypatch)
    assert committed_path.read_text() == _dump(original)


@pytest.mark.parametrize("delta", [5e-10, 0.1])
def test_new_oracle_assessments_use_the_scoped_comparator_through_the_freshness_gate(delta, tmp_path, monkeypatch):
    original = _coverage_artifact()
    changed = deepcopy(original)
    changed["oracle_controls"][0]["assessment"]["importance"]["ess"] += delta
    committed_path = tmp_path / "coverage.json"
    committed_path.write_text(_dump(original))
    monkeypatch.setitem(sys.modules, COVERAGE_NAME, SimpleNamespace(OUT_JSON=committed_path))
    monkeypatch.setattr(freshness, "_reproduce", lambda mod, script: deepcopy(changed))
    if delta < 1e-9:
        freshness.test_the_committed_json_is_what_the_generator_writes(COVERAGE_NAME, tmp_path, monkeypatch)
    else:
        with pytest.raises(AssertionError, match="derived float differs"):
            freshness.test_the_committed_json_is_what_the_generator_writes(COVERAGE_NAME, tmp_path, monkeypatch)
    assert committed_path.read_text() == _dump(original)
