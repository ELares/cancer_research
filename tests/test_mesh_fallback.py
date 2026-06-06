#!/usr/bin/env python3
"""Unit tests for the MeSH-descriptor evidence fallback (#346).

The fallback is an off-by-default (`FERRO_MESH_EXPANSION=1`) controlled-vocabulary
layer in `tag_articles.match_evidence_level`: when the pub_type + keyword passes
find no evidence level, it infers one from an article's MeSH descriptors via
exact set-membership. These tests lock the precision-safety invariants
(conflict order, the editorial veto, the empty-MeSH floor, pub_type precedence,
the review guard) and, crucially, that the layer is byte-reversible (flag OFF
reproduces the pre-#346 baseline).

Run: pytest tests/test_mesh_fallback.py -v
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def _fm(mesh=None, pub_types=None):
    """Minimal article front matter for the matcher."""
    return {"mesh_terms": list(mesh or []), "pub_types": list(pub_types or [])}


class TestMatchEvidenceMesh:
    """The pure `match_evidence_mesh(fm)` helper (independent of the flag)."""

    def test_invivo_from_animal_descriptor(self):
        from tag_articles import match_evidence_mesh
        assert match_evidence_mesh(_fm(["Animals", "Mice"])) == "preclinical-invivo"

    def test_invitro_from_cell_line_descriptor(self):
        from tag_articles import match_evidence_mesh
        assert match_evidence_mesh(_fm(["Cell Line, Tumor"])) == "preclinical-invitro"

    def test_invivo_beats_invitro(self):
        """A paper with both animal and cell-line MeSH is an in-vivo study."""
        from tag_articles import match_evidence_mesh
        assert (
            match_evidence_mesh(_fm(["Cell Line, Tumor", "Animals"]))
            == "preclinical-invivo"
        )

    def test_clinical_other_from_study_design_descriptor(self):
        from tag_articles import match_evidence_mesh
        assert (
            match_evidence_mesh(_fm(["Retrospective Studies"], ["Journal Article"]))
            == "clinical-other"
        )

    def test_clinical_other_vetoed_by_editorial_pubtype(self):
        """Regression guard (PMID 32167722): an editorial citing a cohort is
        not itself evidence, so the clinical-other branch must not fire."""
        from tag_articles import match_evidence_mesh
        assert match_evidence_mesh(_fm(["Prospective Studies"], ["Editorial"])) == ""
        assert match_evidence_mesh(_fm(["Retrospective Studies"], ["Comment"])) == ""

    def test_empty_mesh_returns_empty_floor(self):
        from tag_articles import match_evidence_mesh
        assert match_evidence_mesh(_fm([])) == ""

    def test_non_discriminative_mesh_returns_empty(self):
        """Bare clinical-context MeSH (no study-design / organism descriptor)
        must not be promoted — the precision-safety core."""
        from tag_articles import match_evidence_mesh
        assert (
            match_evidence_mesh(_fm(["Humans", "Stomach Neoplasms", "Middle Aged"]))
            == ""
        )

    def test_deterministic(self):
        from tag_articles import match_evidence_mesh
        fm = _fm(["Animals", "Disease Models, Animal"])
        assert match_evidence_mesh(fm) == match_evidence_mesh(fm) == "preclinical-invivo"


class TestEvidenceLevelFallbackWiring:
    """`match_evidence_level` integration: flag gating + precedence + guards."""

    def test_flag_off_reproduces_baseline(self, monkeypatch):
        """Off by default: a MeSH-only article resolves to "" exactly as the
        pre-#346 tagger did (the byte-reversibility invariant)."""
        import tag_articles
        monkeypatch.setattr(tag_articles, "EVIDENCE_USE_MESH_FALLBACK", False)
        fm = _fm(["Animals", "Mice"], ["Journal Article"])
        assert tag_articles.match_evidence_level(fm, "") == ""

    def test_flag_on_enables_mesh_fallback(self, monkeypatch):
        import tag_articles
        monkeypatch.setattr(tag_articles, "EVIDENCE_USE_MESH_FALLBACK", True)
        fm = _fm(["Animals", "Mice"], ["Journal Article"])
        assert tag_articles.match_evidence_level(fm, "") == "preclinical-invivo"

    def test_pubtype_precedence_over_mesh(self, monkeypatch):
        """An authoritative phase-3 pub_type must win over animal MeSH."""
        import tag_articles
        monkeypatch.setattr(tag_articles, "EVIDENCE_USE_MESH_FALLBACK", True)
        fm = _fm(["Animals"], ["Clinical Trial, Phase III"])
        assert tag_articles.match_evidence_level(fm, "") == "phase3-clinical"

    def test_review_guard_runs_before_mesh(self, monkeypatch):
        """A review-like article is never MeSH-tagged (the upstream guard wins)."""
        import tag_articles
        monkeypatch.setattr(tag_articles, "EVIDENCE_USE_MESH_FALLBACK", True)
        fm = _fm(["Animals", "Mice"], ["Review"])
        assert tag_articles.match_evidence_level(fm, "") == ""

    def test_keyword_pass_still_wins_when_present(self, monkeypatch):
        """If the existing keyword pass resolves a level, MeSH never runs."""
        import tag_articles
        from config import EVIDENCE_LEVEL_KEYWORDS
        monkeypatch.setattr(tag_articles, "EVIDENCE_USE_MESH_FALLBACK", True)
        # Pick a preclinical-invitro keyword and an article whose MeSH would say invivo.
        invitro_kw = next(iter(EVIDENCE_LEVEL_KEYWORDS["preclinical-invitro"]))
        fm = _fm(["Animals"], ["Journal Article"])
        text = f"this study used {invitro_kw.lower()} extensively"
        # Keyword pass resolves invitro; the MeSH invivo fallback must not override.
        assert tag_articles.match_evidence_level(fm, text) == "preclinical-invitro"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
