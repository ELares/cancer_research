"""
Integration tests for pipeline functions.

These test actual function logic with known inputs and expected
outputs (using tolerance bands for floats). They catch formula
errors, keyword matching regressions, and scoring bugs.

Run: pytest tests/test_integration.py -v
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


# ============================================================
# Tagging integration
# ============================================================

class TestTaggingIntegration:
    """Tagging functions should produce correct, deterministic results."""

    def test_match_keywords_deterministic(self):
        """Same input must always produce the same output."""
        from tag_articles import match_keywords
        from config import MECHANISM_KEYWORDS

        text = "immunotherapy checkpoint inhibitor combined with sonodynamic therapy"
        r1 = match_keywords(text, MECHANISM_KEYWORDS)
        r2 = match_keywords(text, MECHANISM_KEYWORDS)
        assert r1 == r2

    def test_match_keywords_crispr(self):
        """CRISPR-related text should match the crispr mechanism."""
        from tag_articles import match_keywords
        from config import MECHANISM_KEYWORDS

        text = "this study uses CRISPR-Cas9 gene editing to target tumor cells"
        result = match_keywords(text, MECHANISM_KEYWORDS)
        assert "crispr" in result, f"Expected 'crispr' in {result}"

    def test_match_keywords_no_match(self):
        """Text with no mechanism-related content should return empty."""
        from tag_articles import match_keywords
        from config import MECHANISM_KEYWORDS

        text = "this paper is about cooking recipes and gardening tips"
        result = match_keywords(text, MECHANISM_KEYWORDS)
        assert len(result) == 0, f"Expected no matches, got {result}"

    def test_match_keywords_multi_mechanism(self):
        """Text mentioning multiple mechanisms should return all of them."""
        from tag_articles import match_keywords
        from config import MECHANISM_KEYWORDS

        text = "nanoparticle-based delivery of checkpoint immunotherapy with CAR-T cell combination"
        result = match_keywords(text, MECHANISM_KEYWORDS)
        assert len(result) >= 2, f"Expected 2+ mechanisms, got {result}"
        assert "nanoparticle" in result
        assert "immunotherapy" in result


# ============================================================
# Evidence weight integration
# ============================================================

class TestEvidenceWeightIntegration:
    """Evidence weight function should return consistent values for known inputs.

    All comparisons use tolerance bands, not exact floating-point matching.
    """

    def test_phase3_high_citation_recent(self):
        """Phase III with high citation and recent year should score ~18.8."""
        from analyze_corpus import evidence_weight

        entry = {"evidence_level": "phase3-clinical", "icite_percentile": 90, "year": 2025}
        w = evidence_weight(entry)
        # phase3(12) × citation(1.45) × recency(~1.082) ≈ 18.82
        assert 18.0 < w < 20.0, f"Expected ~18.8, got {w:.2f}"

    def test_preclinical_low_citation_old(self):
        """Preclinical in-vitro with low citation and old year should score ~0.96."""
        from analyze_corpus import evidence_weight

        entry = {"evidence_level": "preclinical-invitro", "icite_percentile": 10, "year": 2016}
        w = evidence_weight(entry)
        # invitro(1.0) × citation(1.05) × recency(~0.918) ≈ 0.96
        assert 0.9 < w < 1.1, f"Expected ~0.96, got {w:.2f}"

    def test_missing_percentile_defaults_gracefully(self):
        """Missing citation percentile should still produce a positive weight."""
        from analyze_corpus import evidence_weight

        entry = {"evidence_level": "phase1-clinical", "year": 2022}
        w = evidence_weight(entry)
        assert w > 0, f"Should produce positive weight with missing percentile, got {w}"

    def test_missing_year_defaults_gracefully(self):
        """Missing year should still produce a positive weight."""
        from analyze_corpus import evidence_weight

        entry = {"evidence_level": "phase2-clinical", "icite_percentile": 50}
        w = evidence_weight(entry)
        assert w > 0, f"Should produce positive weight with missing year, got {w}"


# ============================================================
# News pipeline integration
# ============================================================

class TestNewsScoreIntegration:
    """News scoring should produce bounded results for edge-case inputs."""

    def test_compute_score_within_bounds(self):
        """Score for any valid input must be in [0, 100]."""
        from score_news import compute_score

        fm = {
            "tier": 1,
            "claims": [],
            "author_credentialed": False,
            "author": "",
            "date_published": None,
        }
        score = compute_score(fm)
        assert 0 <= score <= 100, f"Score {score} out of [0, 100] range"
