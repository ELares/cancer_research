"""Unit tests for the 2026-06 data-quality scripts (#566).

Offline / stdlib-only (no network, no corpus, no compiled extension), so they run
in CI. They lock the load-bearing helpers of the scripts added during the hardening
review: the duplicate-title normaliser (#535), the non-OA full-text strip regex
(#526, incl. the #572 non-greedy fix), and the collaborator-ranking frontmatter
parser (#541).
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))


class TestDetectCorpusDuplicates:
    def test_norm_title_collapses_punctuation_and_case(self):
        import detect_corpus_duplicates as d
        a = d.norm_title("Epigenetic Targeting of PGBD5-Dependent DNA Damage!")
        b = d.norm_title("epigenetic targeting of pgbd5 dependent dna damage")
        assert a == b
        assert " " not in a and "-" not in a

    def test_norm_title_handles_none(self):
        import detect_corpus_duplicates as d
        assert d.norm_title(None) == ""

    def test_preprint_regex_matches_servers_not_journals(self):
        import detect_corpus_duplicates as d
        assert d.PREPRINT_RE.search("bioRxiv : the preprint server for biology")
        assert d.PREPRINT_RE.search("medRxiv")
        assert not d.PREPRINT_RE.search("The Journal of clinical investigation")


class TestStripNonOaFulltext:
    def test_strips_full_text_keeps_abstract(self):
        import strip_non_oa_fulltext as s
        text = "---\npmid: '1'\n---\n\n# Title\n\n## Abstract\n\nthe abstract.\n\n## Full Text\n\nbody.\n"
        out = s.FULLTEXT_HEADING.sub(s.NOTICE.rstrip("\n"), text)
        assert "## Abstract" in out and "the abstract." in out
        assert "body." not in out
        assert "Full text removed (#526" in out

    def test_non_greedy_preserves_trailing_section(self):
        # The #572 fix: a section AFTER Full Text (e.g. a Source footer) must survive.
        import strip_non_oa_fulltext as s
        text = "# T\n\n## Abstract\n\nabs.\n\n## Full Text\n\nbody.\n\n## Source\n\nhttps://doi/x\n"
        out = s.FULLTEXT_HEADING.sub(s.NOTICE.rstrip("\n"), text)
        assert "## Source" in out and "https://doi/x" in out
        assert "body." not in out

    def test_idempotent_marker_present(self):
        import strip_non_oa_fulltext as s
        assert "Full text removed (#526" in s.NOTICE


class TestRankCollaboratorCandidates:
    def test_parse_frontmatter_reads_author_block_and_scalars(self):
        import rank_collaborator_candidates as r
        text = (
            "---\n"
            "pmid: '12345'\n"
            "title: A study\n"
            "authors:\n"
            "- Smith John\n"
            "- Doe Jane\n"
            "journal: Nature\n"
            "year: 2024\n"
            "cited_by_count: 42\n"
            "---\n\n"
            "## Abstract\n\nbody\n"
        )
        fm = r.parse_frontmatter(text)
        assert fm["authors"] == ["Smith John", "Doe Jane"]
        assert fm["journal"] == "Nature"
        assert fm["year"] == "2024"
        assert fm["cited_by_count"] == "42"

    def test_parse_frontmatter_no_frontmatter_returns_empty(self):
        import rank_collaborator_candidates as r
        assert r.parse_frontmatter("no frontmatter here") == {}

    def test_relevance_matches_ferroptosis_terms(self):
        import rank_collaborator_candidates as r
        assert r.RELEVANCE.search("a study of GPX4 and FSP1 in ferroptosis")
        assert r.RELEVANCE.search("erastin induces system xc- inhibition")  # SLC7A11/xCT path
        assert not r.RELEVANCE.search("a study of unrelated kinase signalling")

    def test_score_rewards_recent_and_frequent(self):
        import rank_collaborator_candidates as r
        many_recent = {"papers": 5, "recent": 5, "latest_year": 2025,
                       "citations": 100, "journals": {"a", "b"}, "pmids": []}
        one_old = {"papers": 1, "recent": 0, "latest_year": 2010,
                   "citations": 5, "journals": {"a"}, "pmids": []}
        assert r.score(many_recent) > r.score(one_old)
