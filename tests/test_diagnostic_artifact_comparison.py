"""Diagnostic rounding tolerance cannot alter frozen inputs or decisions."""

from copy import deepcopy
import json
import math
import sys
from types import SimpleNamespace

import pytest

import test_artifact_freshness as freshness


NAME = "proposal_correlated_diagnostics"
RUN = ("runs", 0, "diagnostics")
FEATURE = (*RUN, "features", "angular_sin_2")
REGION = (*RUN, "decomposition", "regions", "sector_0")
POOLED = ("pilots", 0, "diagnostics", "pooled")
ISLAND = ("pilots", 0, "diagnostics", "islands", 0, "summary")
DERIVED_PATHS = (
    *((*RUN, "weights", field) for field in (
        "ess", "max_normalized_weight", "top_1_weight", "top_5_weight", "top_20_weight")),
    *((*FEATURE, field) for field in (
        "estimate", "signed_error", "conditional_mcse", "error_over_mcse")),
    *((*FEATURE, "single_deletion", field) for field in ("signed_change", "absolute_change")),
    *((*RUN, "decomposition", field) for field in (
        "total_signed_error", "between_regions", "within_regions")),
    *((*REGION, field) for field in (
        "estimated_mass", "conditional_estimate", "between_contribution", "within_contribution")),
    (*POOLED, "feature_means", "angular_sin_2"),
    (*POOLED, "region_fractions", "sector_0"),
    (*ISLAND, "feature_means", "angular_sin_2"),
    (*ISLAND, "region_fractions", "sector_0"),
)


def _artifact():
    endpoint = {
        "n_endpoints": 256, "n_unique_endpoints": 250,
        "feature_means": {"angular_sin_2": 0.55},
        "region_counts": {"sector_0": 40}, "region_fractions": {"sector_0": 0.15625},
    }
    diagnostic = {
        "weights": {
            "n_attempts": 8192, "n_positive": 2000, "ess": 400.0,
            "max_normalized_weight": 0.012, "top_1_weight": 0.012,
            "top_5_weight": 0.04, "top_20_weight": 0.12,
        },
        "features": {"angular_sin_2": {
            "estimate": 0.58, "signed_error": 0.08, "conditional_mcse": 0.025,
            "error_over_mcse": 3.2,
            "single_deletion": {"attempt_index": 50, "signed_change": -0.001,
                                "absolute_change": 0.001},
        }},
        "decomposition": {
            "total_signed_error": 0.08, "between_regions": 0.03, "within_regions": 0.05,
            "regions": {"sector_0": {
                "n_accepted": 260, "estimated_mass": 0.13, "conditional_estimate": 0.8,
                "between_contribution": 0.004, "within_contribution": 0.005,
            }},
        },
    }
    return {
        "schema_version": 1,
        "provenance": {
            "source_report_sha256": "a" * 64,
            "archive_sha256": {"annulus-17.json.gz": "b" * 64},
            "source_hashes": {"scripts/example.py": "c" * 64},
            "specification": {"epsilon": 1.0, "production_attempts_per_arm": 8192},
        },
        "truth": {"annular_cylinder": {
            "moments": {"angular_sin_2": 0.5}, "region_masses": {"sector_0": 0.125},
            "conditional_means": {"sector_0": 0.5 + 1 / math.pi},
        }},
        "selected_features": {"annular_cylinder": "angular_sin_2"},
        "completed_benchmark": {"passed": False, "correlated_pass_count": 10},
        "runs": [{
            "fixture": "annular_cylinder", "seed": 17, "arm": "correlated",
            "archive": "annulus-17.json.gz", "original_passed": False,
            "original_failed_checks": ["moments"], "diagnostics": diagnostic,
            "raw": {"points": [[0.25]], "log_q": [2.5], "wall_seconds": 0.25,
                    "proposal": {"uniform_weight": 0.2}},
        }],
        "pilots": [{
            "fixture": "annular_cylinder", "seed": 17, "archive": "annulus-17.json.gz",
            "diagnostics": {"pooled": deepcopy(endpoint),
                            "islands": [{"island": 0, "summary": deepcopy(endpoint)}]},
        }],
        # Same leaf names at other paths must remain exact.
        "metadata": {"runs": [{"diagnostics": deepcopy(diagnostic)}]},
        "positive_runs": [{"assessment": {"estimate": 0.58}}],
    }


def _dump(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _replace(value, path, replacement):
    for part in path[:-1]:
        value = value[part]
    value[path[-1]] = replacement


def _value(value, path):
    for part in path:
        value = value[part]
    return value


def test_diagnostic_scope_does_not_extend_the_four_frozen_study_names():
    assert freshness.DIAGNOSTIC_GENERATOR == NAME
    assert freshness.NUMERICAL_REASSEMBLY == frozenset({
        "proposal_synthetic_validation", "proposal_synthetic_validation_v2",
        "proposal_coverage_challenges", "proposal_correlated_study",
    })


@pytest.mark.parametrize("path", DERIVED_PATHS)
def test_every_declared_derived_metric_allows_last_bit_replay_differences(path):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, path, math.nextafter(_value(original, path), math.inf))
    freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))
    freshness._assert_reassembled_json(NAME, _dump(original), _dump(changed))


@pytest.mark.parametrize("old,new", [
    (0.0, 5e-15), (0.0, -0.0), (1000.0, 1000.0 + 5e-10),
    (-1000.0, -1000.0 - 5e-10),
])
def test_small_finite_derived_differences_use_the_existing_numeric_limits(old, new):
    original = _artifact()
    _replace(original, (*FEATURE, "signed_error"), old)
    changed = deepcopy(original)
    _replace(changed, (*FEATURE, "signed_error"), new)
    freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("path", DERIVED_PATHS)
def test_material_changes_to_each_derived_metric_are_rejected(path):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, path, _value(original, path) + 0.01)
    with pytest.raises(AssertionError, match="derived float differs"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("old,new", [
    (0.0, 2e-14), (0.25, 0.25 + 5e-13), (1000.0, 1000.0 + 2e-9),
])
def test_absolute_and_relative_tolerances_have_material_boundaries(old, new):
    original = _artifact()
    _replace(original, (*FEATURE, "signed_error"), old)
    changed = deepcopy(original)
    _replace(changed, (*FEATURE, "signed_error"), new)
    with pytest.raises(AssertionError, match="derived float differs"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("path", [
    ("truth", "annular_cylinder", "moments", "angular_sin_2"),
    ("truth", "annular_cylinder", "region_masses", "sector_0"),
    ("truth", "annular_cylinder", "conditional_means", "sector_0"),
    ("provenance", "specification", "epsilon"),
    ("runs", 0, "raw", "points", 0, 0),
    ("runs", 0, "raw", "log_q", 0),
    ("runs", 0, "raw", "wall_seconds"),
    ("runs", 0, "raw", "proposal", "uniform_weight"),
    ("metadata", "runs", 0, "diagnostics", "features", "angular_sin_2", "estimate"),
    ("positive_runs", 0, "assessment", "estimate"),
])
def test_one_ulp_changes_to_truth_raw_inputs_provenance_and_lookalikes_are_rejected(path):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, path, math.nextafter(_value(original, path), math.inf))
    with pytest.raises(AssertionError, match="exact JSON value differs"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("path,replacement", [
    (("provenance", "source_report_sha256"), "d" * 64),
    (("provenance", "archive_sha256", "annulus-17.json.gz"), "d" * 64),
    (("provenance", "source_hashes", "scripts/example.py"), "d" * 64),
    (("selected_features", "annular_cylinder"), "angular_cos_2"),
    (("completed_benchmark", "passed"), True),
    (("completed_benchmark", "correlated_pass_count"), 11),
    (("runs", 0, "fixture"), "unequal_balls"),
    (("runs", 0, "seed"), 18),
    (("runs", 0, "arm"), "unbounded_diagonal"),
    (("runs", 0, "archive"), "annulus-18.json.gz"),
    (("runs", 0, "original_passed"), True),
    (("runs", 0, "original_failed_checks", 0), "region_masses"),
    ((*RUN, "weights", "n_attempts"), 8193),
    ((*RUN, "weights", "n_positive"), 2001),
    ((*FEATURE, "single_deletion", "attempt_index"), 51),
    ((*REGION, "n_accepted"), 261),
    ((*POOLED, "n_endpoints"), 257),
    ((*POOLED, "n_unique_endpoints"), 251),
    ((*POOLED, "region_counts", "sector_0"), 41),
    (("pilots", 0, "diagnostics", "islands", 0, "island"), 1),
    ((*ISLAND, "n_endpoints"), 257),
    ((*ISLAND, "region_counts", "sector_0"), 41),
])
def test_hashes_counts_flags_selected_features_and_indices_remain_exact(path, replacement):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, path, replacement)
    with pytest.raises(AssertionError, match="exact JSON value differs"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("path,replacement", [
    ((*FEATURE, "estimate"), 1),
    ((*FEATURE, "conditional_mcse"), None),
    ((*FEATURE, "single_deletion", "attempt_index"), 50.0),
    ((*RUN, "weights", "n_attempts"), 8192.0),
    (("runs", 0, "original_passed"), 0),
    (("schema_version",), True),
])
def test_json_types_remain_exact_at_derived_and_structural_paths(path, replacement):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, path, replacement)
    with pytest.raises(AssertionError, match="JSON value types differ"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


def test_large_integer_count_does_not_receive_relative_float_tolerance():
    original = _artifact()
    _replace(original, (*RUN, "weights", "n_attempts"), 10**15)
    changed = deepcopy(original)
    _replace(changed, (*RUN, "weights", "n_attempts"), 10**15 + 1)
    with pytest.raises(AssertionError, match="exact JSON value differs"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("path", [
    (*RUN, "weights", "n_attempts"), (*REGION, "n_accepted"),
    (*FEATURE, "single_deletion", "attempt_index"),
    (*POOLED, "n_endpoints"), (*ISLAND, "region_counts", "sector_0"),
])
def test_count_and_index_paths_never_acquire_float_tolerance(path):
    # The path rule must remain exact even if both artifacts have the same
    # malformed scalar type; type equality alone does not establish a count.
    original = _artifact()
    _replace(original, path, 1000.0)
    changed = deepcopy(original)
    _replace(changed, path, math.nextafter(1000.0, math.inf))
    with pytest.raises(AssertionError, match="exact JSON value differs"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("path", [
    (*FEATURE, "truth"), (*FEATURE, "point", 0),
    (*RUN, "weights", "new_metric"), (*REGION, "true_mass"),
    (*ISLAND, "raw_coordinates", 0),
])
def test_unlisted_floats_inside_diagnostics_are_not_implicitly_tolerated(path):
    original = _artifact()
    container = original
    for part in path[:-1]:
        if isinstance(container, dict) and part not in container:
            container[part] = [0.25]
        container = container[part]
    container[path[-1]] = 0.25
    changed = deepcopy(original)
    _replace(changed, path, math.nextafter(0.25, math.inf))
    with pytest.raises(AssertionError, match="exact JSON value differs"):
        freshness._assert_reassembled_json(NAME, _dump(changed), _dump(original))


@pytest.mark.parametrize("name", ["proposal_correlated_diagnostics_copy", "other_diagnostics"])
def test_other_generator_names_do_not_inherit_diagnostic_tolerance(name):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, (*FEATURE, "estimate"), math.nextafter(0.58, math.inf))
    with pytest.raises(AssertionError, match="generated JSON bytes differ"):
        freshness._assert_reassembled_json(name, _dump(changed), _dump(original))


@pytest.mark.parametrize("name", sorted(freshness.NUMERICAL_REASSEMBLY))
def test_frozen_study_comparison_scope_does_not_gain_diagnostic_paths(name):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, (*FEATURE, "estimate"), math.nextafter(0.58, math.inf))
    with pytest.raises(AssertionError, match="exact JSON value differs"):
        freshness._assert_reassembled_json(name, _dump(changed), _dump(original))


@pytest.mark.parametrize("side", ["produced", "committed"])
@pytest.mark.parametrize("formatting", ["compact", "no_newline", "extra_newline", "unsorted"])
def test_formatting_must_match_the_diagnostic_writer(side, formatting):
    value = _artifact()
    canonical = _dump(value)
    malformed = {
        "compact": json.dumps(value, sort_keys=True) + "\n",
        "no_newline": canonical[:-1],
        "extra_newline": canonical + "\n",
        "unsorted": json.dumps(value, indent=2) + "\n",
    }[formatting]
    produced, committed = ((malformed, canonical) if side == "produced" else (canonical, malformed))
    with pytest.raises(AssertionError, match=f"{side} JSON formatting differs"):
        freshness._assert_reassembled_json(NAME, produced, committed)


@pytest.mark.parametrize("nonfinite", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_derived_values_are_rejected(nonfinite):
    original = _artifact()
    changed = deepcopy(original)
    _replace(changed, (*FEATURE, "estimate"), nonfinite)
    invalid = json.dumps(changed, indent=2, sort_keys=True) + "\n"
    with pytest.raises(AssertionError, match="invalid produced JSON numbers"):
        freshness._assert_reassembled_json(NAME, invalid, _dump(original))


@pytest.mark.parametrize("delta", [5e-15, 0.01])
def test_diagnostic_comparator_is_used_by_the_actual_freshness_gate(delta, tmp_path, monkeypatch):
    original = _artifact()
    changed = deepcopy(original)
    changed["runs"][0]["diagnostics"]["features"]["angular_sin_2"]["estimate"] += delta
    path = tmp_path / "diagnostic.json"
    path.write_text(_dump(original))
    monkeypatch.setitem(sys.modules, NAME, SimpleNamespace(OUT_JSON=path))
    monkeypatch.setattr(freshness, "_reproduce", lambda mod, name: deepcopy(changed))
    if delta < 1e-12:
        freshness.test_the_committed_json_is_what_the_generator_writes(NAME, tmp_path, monkeypatch)
    else:
        with pytest.raises(AssertionError, match="derived float differs"):
            freshness.test_the_committed_json_is_what_the_generator_writes(NAME, tmp_path, monkeypatch)
    assert path.read_text() == _dump(original)
