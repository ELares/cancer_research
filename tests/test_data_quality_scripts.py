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


class TestScoreNewsCrossCitation:
    """#571: the 10% cross-citation term now measures the count of DISTINCT corpus
    PMIDs cited across all of an article's claims, saturating at 3 — replacing the
    fraction-of-anchored-claims term that was redundant with verified_ratio.

    Each fixture pins the other three terms deterministically so the score isolates
    the cross-citation contribution: tier 1 (weight 1.0), author_credentialed
    (author 1.0), no FACTUAL claims (verified_ratio 1.0), date_published None
    (recency 0.5, date-independent). So score == 40 + 30 + 10 + 10*cross == 80 +
    10*cross, and cross == min(1, distinct_corpus / 3).
    """

    @staticmethod
    def _score(monkeypatch, linked_lists, corpus):
        import score_news
        monkeypatch.setattr(score_news, "_corpus_pmids", lambda: {str(p) for p in corpus})
        fm = {
            "tier": 1,
            "author_credentialed": True,
            "date_published": None,
            # category != FACTUAL => verified_ratio 1.0, claims still feed the
            # distinct-corpus comprehension (which scans ALL claims).
            "claims": [{"category": "BACKGROUND", "linked_pmids": ll} for ll in linked_lists],
        }
        return score_news.compute_score(fm)

    def test_saturates_at_three_distinct_corpus_pmids(self, monkeypatch):
        corpus = {"1", "2", "3", "4"}
        assert self._score(monkeypatch, [["1"]], corpus) == 83.3              # 1 distinct -> 1/3
        assert self._score(monkeypatch, [["1"], ["2"]], corpus) == 86.7        # 2 distinct -> 2/3
        assert self._score(monkeypatch, [["1"], ["2"], ["3"]], corpus) == 90.0  # 3 distinct -> 1.0
        assert self._score(monkeypatch, [["1", "2", "3", "4"]], corpus) == 90.0  # 4 -> capped at 1.0

    def test_dedups_ignores_noncorpus_and_handles_none(self, monkeypatch):
        corpus = {"1", "2", "3"}
        # same corpus PMID across three claims counts once -> 1 distinct -> 1/3
        assert self._score(monkeypatch, [["1"], ["1"], ["1"]], corpus) == 83.3
        # non-corpus PMIDs ignored; a None linked_pmids is safe -> 2 distinct -> 2/3
        assert self._score(monkeypatch, [["1", "999"], None, ["2", "888"]], corpus) == 86.7
        # no corpus citations at all -> cross 0 -> 80.0
        assert self._score(monkeypatch, [["999"], None], corpus) == 80.0

    def test_int_pmids_reconcile_with_string_corpus(self, monkeypatch):
        corpus = {"1", "2", "3"}
        # linked_pmids given as ints must match the string corpus via str(): 2 distinct -> 2/3
        assert self._score(monkeypatch, [[1, 2]], corpus) == 86.7
