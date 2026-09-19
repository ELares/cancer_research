"""Tests for the practical-identifiability synthesis (#503).

Pure-Python stdlib (no compiled extension), so these run in CI. They guard the
load-bearing facts (degrees of freedom, the non-identifiable count, that no
headline is point-estimable, that zero headlines are data-conditioned in the
production regime) and that the committed report is not stale.
"""

import sys
import pytest
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
    """Uncertain direction is permitted; no headline is point-estimable."""
    r = ir.build()
    allowed = {"directional_only", "direction_robust_magnitude_not",
               "direction_and_magnitude_uncertain"}
    for h in r["headlines"]:
        assert h["verdict"] in allowed, f"{h['key']} claims {h['verdict']}"
    # At least the single-cell kill rate and the immune ratio are directional-only.
    by_key = {h["key"]: h for h in r["headlines"]}
    assert by_key["single_cell_kill_rate"]["verdict"] == "directional_only"
    assert by_key["immune_amplification_ratio"]["verdict"] == "directional_only"


def test_bliss_verdict_tracks_the_sampled_range_not_a_rounded_null(tmp_path, monkeypatch):
    """Exercise both sides of the inference, including an exactly additive draw."""
    source = tmp_path / "uncertainty.md"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    for minimum, lower, expected in (
        (0.953, 1.000, "direction_and_magnitude_uncertain"),
        (1.000, 1.100, "direction_and_magnitude_uncertain"),
        (1.050, 1.100, "direction_robust_magnitude_not"),
    ):
        source.write_text(
            "| default point estimate (no perturbation) | 1.992 |\n"
            "| prior-predictive median | 1.348 |\n"
            f"| 95% prior-predictive interval | [{lower:.3f}, 5.242] |\n"
            f"| full range (min, max) | [{minimum:.3f}, 7.800] |\n"
        )
        report = ir.build()
        bliss = next(h for h in report["headlines"] if h["key"] == "bliss_synergy")
        assert bliss["verdict"] == expected
        assert f"[{minimum:.3f}, 7.800]" in bliss["prior_predictive"]
        assert ("not established" in report["overall"]) == (minimum <= 1.0)
        if minimum == 1.0:
            assert "precision does not establish" in bliss["rationale"]
            assert "is not uniform" not in bliss["rationale"]
        assert "not equivalent to exceeding Bliss independence" in bliss["rationale"]


def test_missing_bliss_evidence_cannot_reuse_a_hardcoded_verdict(tmp_path, monkeypatch):
    source = tmp_path / "uncertainty.md"
    source.write_text("| full range (min, max) | [0.953, 7.800] |\n")
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    with pytest.raises(ValueError, match="Missing or malformed Bliss uncertainty row"):
        ir.build()


def test_report_is_deterministic_and_committed_json_is_fresh():
    """build() is deterministic, and the committed JSON must match a fresh build
    (so the committed report cannot silently go stale)."""
    import json

    r1 = ir.build()
    r2 = ir.build()
    assert r1 == r2
    committed = json.loads((REPO_ROOT / "analysis" / "identifiability-report.json").read_text())
    assert committed == r1, "analysis/identifiability-report.json is stale; re-run scripts/identifiability_report.py"
