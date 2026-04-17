"""
Smoke tests for the Python corpus pipeline.

These verify that core scripts import cleanly and produce expected outputs
without running the full pipeline (which takes minutes). Each test is
designed to catch obvious breakage: missing imports, broken function
signatures, and empty outputs.

Run: pytest tests/test_pipeline_smoke.py -v
"""

import importlib
import json
import sys
from pathlib import Path

import pytest

# Add scripts/ to path so we can import pipeline modules
REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

CORPUS_INDEX = REPO_ROOT / "corpus" / "INDEX.jsonl"
ANALYSIS_DIR = REPO_ROOT / "analysis"


# ============================================================
# Import tests: verify modules load without error
# ============================================================

class TestImports:
    def test_config_imports(self):
        config = importlib.import_module("config")
        assert hasattr(config, "MECHANISM_KEYWORDS")
        assert hasattr(config, "CANCER_TYPE_KEYWORDS")
        assert hasattr(config, "DIAGNOSTIC_THERAPY_KEYWORDS")
        assert hasattr(config, "TISSUE_CATEGORY_ORDER")
        assert len(config.MECHANISM_KEYWORDS) > 0

    def test_tag_articles_imports(self):
        mod = importlib.import_module("tag_articles")
        assert hasattr(mod, "match_keywords")
        assert hasattr(mod, "match_diagnostic_therapy_links")
        assert hasattr(mod, "match_evidence_level")

    def test_evidence_utils_imports(self):
        mod = importlib.import_module("evidence_utils")
        assert hasattr(mod, "normalize_text")
        assert hasattr(mod, "is_review_like")

    def test_analyze_corpus_imports(self):
        mod = importlib.import_module("analyze_corpus")
        assert hasattr(mod, "build_mechanism_matrix")
        assert hasattr(mod, "build_diagnostic_therapy_audit")

    def test_generate_figures_imports(self):
        # This imports matplotlib which sets Agg backend
        mod = importlib.import_module("generate_figures")
        assert hasattr(mod, "load_corpus")
        assert hasattr(mod, "load_index")

    def test_generate_latex_imports(self):
        # generate_latex is a script, not a module with functions,
        # but we can verify it parses without SyntaxError
        path = SCRIPTS_DIR / "generate_latex.py"
        assert path.exists()
        compile(path.read_text(), str(path), "exec")


# ============================================================
# Keyword/config consistency tests
# ============================================================

class TestConfigConsistency:
    def test_diagnostic_therapy_order_matches_keywords(self):
        from config import DIAGNOSTIC_THERAPY_ORDER, DIAGNOSTIC_THERAPY_KEYWORDS
        assert set(DIAGNOSTIC_THERAPY_ORDER) == set(DIAGNOSTIC_THERAPY_KEYWORDS.keys())

    def test_tissue_order_matches_mapping(self):
        from config import TISSUE_CATEGORY_ORDER, CANCER_TYPE_TO_TISSUE
        tissue_values = set(CANCER_TYPE_TO_TISSUE.values())
        assert tissue_values.issubset(set(TISSUE_CATEGORY_ORDER))

    def test_cancer_subtype_order_matches_keywords(self):
        from config import CANCER_SUBTYPE_ORDER, CANCER_SUBTYPE_KEYWORDS
        assert set(CANCER_SUBTYPE_ORDER) == set(CANCER_SUBTYPE_KEYWORDS.keys())

    def test_all_mechanisms_have_keywords(self):
        from config import MECHANISM_KEYWORDS
        for mech, keywords in MECHANISM_KEYWORDS.items():
            assert len(keywords) > 0, f"Mechanism {mech} has no keywords"

    def test_diagnostic_therapy_chains_have_all_fields(self):
        from config import DIAGNOSTIC_THERAPY_KEYWORDS
        for chain_id, chain in DIAGNOSTIC_THERAPY_KEYWORDS.items():
            assert "diagnostic" in chain, f"{chain_id} missing 'diagnostic'"
            assert "feature" in chain, f"{chain_id} missing 'feature'"
            assert "intervention" in chain, f"{chain_id} missing 'intervention'"


# ============================================================
# Corpus index tests
# ============================================================

class TestCorpusIndex:
    @pytest.fixture
    def index_entries(self):
        if not CORPUS_INDEX.exists():
            pytest.skip("INDEX.jsonl not found")
        entries = []
        with open(CORPUS_INDEX) as f:
            for line in f:
                entries.append(json.loads(line))
        return entries

    def test_index_is_nonempty(self, index_entries):
        assert len(index_entries) > 4000, f"Expected 4000+ entries, got {len(index_entries)}"

    def test_index_has_required_fields(self, index_entries):
        required = {"pmid", "title", "mechanisms", "cancer_types", "evidence_level"}
        sample = index_entries[0]
        for field in required:
            assert field in sample, f"Missing field: {field}"

    def test_index_has_diagnostic_therapy_field(self, index_entries):
        # At least some entries should have the new field
        with_links = [e for e in index_entries if e.get("diagnostic_therapy_links")]
        assert len(with_links) > 50, f"Expected 50+ entries with diagnostic_therapy_links, got {len(with_links)}"

    def test_index_has_tissue_categories(self, index_entries):
        with_tissue = [e for e in index_entries if e.get("tissue_categories")]
        assert len(with_tissue) > 2000, f"Expected 2000+ entries with tissue_categories, got {len(with_tissue)}"


# ============================================================
# Analysis output tests
# ============================================================

class TestAnalysisOutputs:
    EXPECTED_OUTPUTS = [
        "mechanism-matrix.md",
        "convergence-map.md",
        "evidence-tiers.md",
        "gap-analysis.md",
        "evidence-coverage-audit.md",
        "tissue-mechanism-summary.md",
        "tissue-evidence-summary.md",
        "weighted-evidence-summary.md",
        "diagnostic-therapy-audit.md",
        "sarcoma-subtype-audit.md",
        "pathway-target-audit.md",
        "radioligand-audit.md",
        "designed-combinations.md",
    ]

    @pytest.mark.parametrize("filename", EXPECTED_OUTPUTS)
    def test_analysis_output_exists_and_nonempty(self, filename):
        path = ANALYSIS_DIR / filename
        if not path.exists():
            pytest.skip(f"{filename} not found")
        content = path.read_text()
        assert len(content) > 100, f"{filename} is too small ({len(content)} chars)"

    def test_evidence_gold_eval_exists(self):
        path = ANALYSIS_DIR / "evidence-gold-eval.md"
        if not path.exists():
            pytest.skip("evidence-gold-eval.md not found")
        content = path.read_text()
        assert "46/100" in content or "46.0%" in content, "Gold-set eval should report 46% accuracy"


# ============================================================
# Matching function tests
# ============================================================

class TestMatchingFunctions:
    def test_match_keywords_basic(self):
        from tag_articles import match_keywords
        from config import MECHANISM_KEYWORDS
        text = "this paper describes immunotherapy with checkpoint inhibitors for breast cancer treatment"
        matches = match_keywords(text, MECHANISM_KEYWORDS)
        assert "immunotherapy" in matches

    def test_match_diagnostic_therapy_requires_intervention(self):
        from tag_articles import match_diagnostic_therapy_links
        # Only diagnostic language, no intervention — should NOT match
        text = "pd-l1 immunohistochemistry showed high pd-l1 expression in the tumor"
        matches = match_diagnostic_therapy_links(text)
        assert len(matches) == 0, "Should not match without intervention keyword"

    def test_match_diagnostic_therapy_with_intervention(self):
        from tag_articles import match_diagnostic_therapy_links
        # Diagnostic + intervention — should match
        text = "pd-l1 ihc showed high expression and the patient received pembrolizumab"
        matches = match_diagnostic_therapy_links(text)
        assert "pdl1-ihc-to-checkpoint" in matches

    def test_normalize_text(self):
        from evidence_utils import normalize_text
        result = normalize_text("  Hello   World  ")
        assert result == "hello world"
