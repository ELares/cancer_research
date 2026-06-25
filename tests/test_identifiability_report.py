"""Tests for the practical-identifiability synthesis (#503).

Pure-Python stdlib (no compiled extension), so these run in CI. They guard the
load-bearing facts (degrees of freedom, the non-identifiable count, that no
headline is point-estimable, that zero headlines are data-conditioned in the
production regime) and that the committed report is not stale.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import identifiability_report as ir  # noqa: E402


def test_degrees_of_freedom_matches_prcc():
    """The DOF count must equal the actual PRCC swept-parameter set (11), so the
    accounting cannot drift away from the real parameter set."""
    r = ir.build()
    assert r["degrees_of_freedom"] == 11
    assert len(r["swept_parameters"]) == 11


def test_non_identifiable_params_are_real_and_count_holds():
    """6 of 11 parameters are non-identifiable from the kill rate, and each named
    one is a genuine PRCC parameter (build() asserts this; re-check here)."""
    r = ir.build()
    sc = r["single_cell_sobol"]
    assert sc["non_identifiable_count"] == 6
    for p in sc["non_identifiable"]:
        assert p in r["swept_parameters"], f"{p} not a swept parameter"
    # The dominant three are also real parameters and dominate.
    for p in sc["dominant_ST"]:
        assert p in r["swept_parameters"]
    assert sc["dominant_ST"]["lp_propagation"] > sc["dominant_ST"]["gpx4_rate"]


def test_zero_headlines_data_conditioned_in_production():
    """The load-bearing finding: no headline is conditioned on data in the
    regime that produces it (the only fit is the disjoint in-vitro switch)."""
    r = ir.build()
    assert r["data_constrained_in_production"] == 0


def test_no_headline_is_point_estimable():
    """Every headline verdict must be directional-only or direction-robust-
    magnitude-not; none may claim a point-estimable magnitude."""
    r = ir.build()
    allowed = {"directional_only", "direction_robust_magnitude_not"}
    for h in r["headlines"]:
        assert h["verdict"] in allowed, f"{h['key']} claims {h['verdict']}"
    # At least the single-cell kill rate and the immune ratio are directional-only.
    by_key = {h["key"]: h for h in r["headlines"]}
    assert by_key["single_cell_kill_rate"]["verdict"] == "directional_only"
    assert by_key["immune_amplification_ratio"]["verdict"] == "directional_only"


def test_report_is_deterministic_and_committed_json_is_fresh():
    """build() is deterministic, and the committed JSON must match a fresh build
    (so the committed report cannot silently go stale)."""
    import json

    r1 = ir.build()
    r2 = ir.build()
    assert r1 == r2
    committed = json.loads((REPO_ROOT / "analysis" / "identifiability-report.json").read_text())
    assert committed == r1, "analysis/identifiability-report.json is stale; re-run scripts/identifiability_report.py"
