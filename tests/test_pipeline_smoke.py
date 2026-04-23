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
        assert path.exists(), f"Analysis output missing: {filename}"
        content = path.read_text()
        assert len(content) > 100, f"{filename} is too small ({len(content)} chars)"

    def test_evidence_gold_eval_has_expected_structure(self):
        path = ANALYSIS_DIR / "evidence-gold-eval.md"
        assert path.exists(), "evidence-gold-eval.md missing"
        content = path.read_text()
        assert "## Overall Metrics" in content, "Gold-set eval should have Overall Metrics section"
        assert "Exact-label accuracy:" in content, "Gold-set eval should report exact-label accuracy"
        assert "## Per-Label Metrics" in content, "Gold-set eval should have Per-Label Metrics section"


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


class TestNewsPipeline:
    """Smoke tests for the news integration pipeline (issue #99)."""

    def test_news_imports(self):
        from fetch_news import classify_source, slugify
        from extract_claims import split_sentences, detect_factual_markers
        from score_news import compute_score
        from build_news_index import build_index

    def test_classify_source_tier1(self):
        from fetch_news import classify_source
        tier, name, domain = classify_source("https://cancer.gov/news-events/some-article")
        assert tier == 1
        assert domain == "cancer.gov"

    def test_classify_source_tier2(self):
        from fetch_news import classify_source
        tier, name, domain = classify_source("https://statnews.com/2026/02/article")
        assert tier == 2

    def test_classify_source_excluded(self):
        from fetch_news import classify_source
        tier, name, domain = classify_source("https://twitter.com/post/123")
        assert tier == 0
        assert name == "Excluded"

    def test_classify_source_longest_prefix(self):
        from fetch_news import classify_source
        # nature.com/news is Tier 1; nature.com/articles is NOT in config
        tier1, _, _ = classify_source("https://nature.com/news/some-article")
        tier0, _, _ = classify_source("https://nature.com/articles/s41586-paper")
        assert tier1 == 1
        assert tier0 == 0

    def test_split_sentences(self):
        from extract_claims import split_sentences
        result = split_sentences("Dr. Smith found a 43.2% response rate. This is important.")
        assert len(result) == 2
        assert "43.2%" in result[0]

    def test_detect_factual_markers(self):
        from extract_claims import detect_factual_markers
        markers = detect_factual_markers("The Phase 3 trial showed a 43% response rate.")
        assert len(markers) >= 2  # percentage + phase

    def test_detect_no_factual_markers(self):
        from extract_claims import detect_factual_markers
        markers = detect_factual_markers("This is an interesting development.")
        assert len(markers) == 0

    def test_classify_claim_type(self):
        from extract_claims import classify_claim_type
        assert classify_claim_type("The FDA approved the drug.") == "event"
        assert classify_claim_type("This could lead to new treatments.") == "speculation"

    def test_score_formula_tier1(self):
        from score_news import compute_score
        fm = {
            "tier": 1,
            "claims": [{"category": "FACTUAL", "verification_status": "verified"}] * 5,
            "author_credentialed": True,
            "author": "Test Author",
            "date_published": "2026-04-01",
        }
        score = compute_score(fm)
        # tier_weight=1.0, verified_ratio=1.0, author=1.0, recency=1.0, cross=0.0
        # 1.0 * (40 + 30 + 20 + 0) = 90
        assert 85 <= score <= 95

    def test_score_zero_factual_claims(self):
        from score_news import compute_score
        fm = {
            "tier": 3,
            "claims": [{"category": "INTERPRETIVE", "verification_status": None}],
            "author_credentialed": True,
            "author": "Expert",
            "date_published": "2026-04-01",
        }
        score = compute_score(fm)
        assert score > 0  # Should not be zero — verified_ratio defaults to 1.0

    def test_slugify(self):
        from fetch_news import slugify
        result = slugify("Key Study of Grail's Cancer Detection Test Fails!")
        assert result.islower() or result.replace("-", "").isalnum()
        assert len(result) <= 60

    def test_verify_imports(self):
        from verify_news_claims import search_corpus, extract_search_terms

    def test_extract_search_terms(self):
        from verify_news_claims import extract_search_terms
        terms = extract_search_terms("The Phase 3 trial of pembrolizumab showed 43% response rate.")
        assert len(terms) > 0

    def test_extract_opinion_claim(self):
        from extract_claims import split_sentences, detect_factual_markers
        from config import CLAIM_OPINION_TRIGGERS
        sentence = "According to Dr. Smith, the results are encouraging."
        # No factual markers
        assert len(detect_factual_markers(sentence)) == 0
        # But should match opinion trigger
        assert any(t in sentence.lower() for t in CLAIM_OPINION_TRIGGERS)

    def test_extract_speculation_claim(self):
        from extract_claims import split_sentences, detect_factual_markers
        from config import CLAIM_SPECULATION_TRIGGERS
        sentence = "This could lead to new treatments within the next decade."
        assert len(detect_factual_markers(sentence)) == 0
        assert any(t in sentence.lower() for t in CLAIM_SPECULATION_TRIGGERS)

    def test_config_has_news_additions(self):
        from config import NEWS_RATE, NEWS_DIR, CLAIM_FACTUAL_MARKERS, CLAIM_TYPE_MARKERS
        assert NEWS_RATE is not None
        assert NEWS_DIR is not None
        assert len(CLAIM_FACTUAL_MARKERS) >= 10
        assert "event" in CLAIM_TYPE_MARKERS
        assert "speculation" in CLAIM_TYPE_MARKERS
