#!/usr/bin/env python3
"""Regression guards for the v2 evidence tagger (#TAGGER-V2).

Two things must hold at once, and they pull in opposite directions:

1. CONTAINMENT. With `FERRO_EVIDENCE_V2` unset the tagger must be byte-identical
   to the frozen corpus, because `corpus/INDEX.jsonl` and every manuscript
   number depend on it. Same contract as the #346 MeSH fallback.

2. THE GAIN MUST SURVIVE. With the flag on, the measured error reduction on the
   human/full-text CONSENSUS subset must not regress. That subset is the
   conservative measurement: both annotators had to agree on the label.

Thresholds are floors, not point pins. The consensus subset has 77 records, so
one record is ~1.3 points; pinning to the measured 75.3% would false-fail on a
benign single-record change. The floors below sit roughly 5 points under the
measured values, which resolves a real regression while tolerating noise.

Run: pytest tests/test_evidence_v2.py -v
"""

import csv
import json
import random
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import evidence_sections  # noqa: E402
from evidence_sections import CITED, DROP, SELF  # noqa: E402

PMID_DIR = REPO_ROOT / "corpus" / "by-pmid"
INDEX = REPO_ROOT / "corpus" / "INDEX.jsonl"
HUMAN = REPO_ROOT / "analysis" / "evidence-gold-labels-v1.csv"
FULLTEXT = REPO_ROOT / "analysis" / "evidence-gold-set-v3-fulltext.csv"

# Measured on the consensus subset: v1 51.9% exact, v2 75.3% (1.95x error
# reduction). Floors sit ~5 points below to tolerate single-record noise.
MIN_V2_EXACT = 0.70
MIN_ERROR_REDUCTION = 1.70
MIN_V2_PRECISION = 0.90


# --------------------------------------------------------------------------
# Section splitter
# --------------------------------------------------------------------------

def test_heading_classification():
    assert evidence_sections.classify_heading("Materials and methods") == SELF
    assert evidence_sections.classify_heading("2.1 Methods") == SELF
    assert evidence_sections.classify_heading("Results") == SELF
    assert evidence_sections.classify_heading("Trial design") == SELF
    assert evidence_sections.classify_heading("Introduction") == CITED
    assert evidence_sections.classify_heading("Discussion") == CITED
    assert evidence_sections.classify_heading("Conclusions") == CITED
    assert evidence_sections.classify_heading("References") == DROP
    assert evidence_sections.classify_heading("Competing interests") == DROP
    # Unrecognised headings must return None so callers can inherit the
    # enclosing class -- Methods subsections vary far too much to enumerate.
    assert evidence_sections.classify_heading("Calcein assay") is None


def test_unknown_subsection_inherits_enclosing_class():
    body = (
        "## Full Text\n\nTitle line\n\n"
        "Introduction\n\ncited prose here\n\n"
        "Materials and methods\n\nmethods prose\n\n"
        "Calcein assay\n\nmore methods prose\n\n"
        "Discussion\n\ndiscussion prose\n"
    )
    got = {h: k for h, k, _ in evidence_sections.split_sections(body)}
    assert got["Introduction"] == CITED
    assert got["Materials and methods"] == SELF
    assert got["Discussion"] == CITED
    # 'Calcein assay' is unrecognised, so it is absorbed into Methods rather
    # than starting a section of its own.
    assert "Calcein assay" not in got
    self_txt = evidence_sections.self_text(body)
    assert "more methods prose" in self_txt
    assert "discussion prose" not in self_txt
    assert "cited prose here" not in self_txt


@pytest.mark.skipif(not (PMID_DIR / "40700574.md").exists(), reason="corpus record absent")
def test_discussion_citation_is_excluded_from_self_text():
    """The worked case that motivated section scoping.

    PMID 40700574 is a mouse study whose Discussion cites someone else's phase 3
    trial. Reading full text whole promotes it to phase3-clinical; SELF-scoped
    text must not contain that sentence.
    """
    body = (PMID_DIR / "40700574.md").read_text(errors="ignore")
    assert "phase 3" in body.lower()
    assert "phase 3" not in evidence_sections.self_text(body).lower()


# --------------------------------------------------------------------------
# Containment: flag OFF must not move the frozen corpus
# --------------------------------------------------------------------------

@pytest.mark.skipif(not INDEX.exists(), reason="corpus index absent")
def test_flag_off_is_byte_identical_to_frozen_corpus():
    import tag_articles
    from article_io import load_article

    assert not tag_articles.EVIDENCE_USE_V2, "v2 must default to OFF"
    assert not tag_articles.EVIDENCE_USE_MESH_FALLBACK, "MeSH fallback must default to OFF"

    rows = [json.loads(line) for line in INDEX.open()]
    random.Random(11).shuffle(rows)
    checked = diffs = 0
    for row in rows[:250]:
        path = PMID_DIR / f"{row['pmid']}.md"
        if not path.exists():
            continue
        fm, body = load_article(path)
        got = tag_articles.match_evidence_level(fm, tag_articles.get_searchable_text(fm, body))
        checked += 1
        if (got or "") != (row.get("evidence_level", "") or ""):
            diffs += 1
    assert checked > 100, "too few records checked to be meaningful"
    assert diffs == 0, f"{diffs}/{checked} frozen evidence_level values changed with the flag OFF"


# --------------------------------------------------------------------------
# The gain
# --------------------------------------------------------------------------

def _consensus_labels():
    human = {r["pmid"]: r["gold_evidence_level"] for r in csv.DictReader(HUMAN.open())}
    with FULLTEXT.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    ft = {r["pmid"]: r["tier_fulltext"] for r in rows}
    return {p: human[p] for p in human if ft.get(p) == human[p]}


def _score(gold, use_v2):
    """Import tag_articles fresh under the requested flag state."""
    import tag_articles
    from article_io import load_article

    prev_v2 = tag_articles.EVIDENCE_USE_V2
    prev_mesh = tag_articles.EVIDENCE_USE_MESH_FALLBACK
    tag_articles.EVIDENCE_USE_V2 = use_v2
    tag_articles.EVIDENCE_USE_MESH_FALLBACK = use_v2
    try:
        exact = tp = fp = 0
        n = 0
        for pmid, want in gold.items():
            path = PMID_DIR / f"{pmid}.md"
            if not path.exists():
                continue
            fm, body = load_article(path)
            got = tag_articles.match_evidence_level(fm, tag_articles.get_searchable_text(fm, body))
            n += 1
            if (got or "none-applicable") == want:
                exact += 1
            if got and want != "none-applicable":
                tp += 1
            elif got and want == "none-applicable":
                fp += 1
        return n, exact / n, (tp / (tp + fp) if tp + fp else 0.0)
    finally:
        tag_articles.EVIDENCE_USE_V2 = prev_v2
        tag_articles.EVIDENCE_USE_MESH_FALLBACK = prev_mesh


@pytest.mark.skipif(not FULLTEXT.exists() or not PMID_DIR.exists(),
                    reason="gold set or corpus absent")
def test_v2_error_reduction_on_consensus_subset():
    gold = _consensus_labels()
    assert len(gold) >= 60, f"consensus subset unexpectedly small: {len(gold)}"

    n, v1_exact, _ = _score(gold, use_v2=False)
    n2, v2_exact, v2_prec = _score(gold, use_v2=True)
    assert n == n2

    assert v2_exact >= MIN_V2_EXACT, (
        f"v2 exact-label accuracy {v2_exact:.1%} fell below the {MIN_V2_EXACT:.0%} floor")
    assert v2_prec >= MIN_V2_PRECISION, (
        f"v2 precision {v2_prec:.1%} fell below the {MIN_V2_PRECISION:.0%} floor -- "
        "the v2 vocabularies are over-firing on none-applicable records")

    reduction = (1 - v1_exact) / (1 - v2_exact)
    assert reduction >= MIN_ERROR_REDUCTION, (
        f"error reduction {reduction:.2f}x fell below the {MIN_ERROR_REDUCTION}x floor "
        f"(v1 {v1_exact:.1%} -> v2 {v2_exact:.1%})")


@pytest.mark.skipif(not PMID_DIR.exists(), reason="corpus absent")
def test_opinion_pubtypes_are_never_primary_evidence():
    """An editorial discussing a phase III trial must not inherit that tier."""
    import tag_articles

    prev = tag_articles.EVIDENCE_USE_V2
    tag_articles.EVIDENCE_USE_V2 = True
    try:
        fm = {"title": "Some commentary", "pub_types": ["Editorial", "Comment"]}
        text = "this phase 3 randomized controlled trial enrolled 500 patients"
        assert tag_articles.match_evidence_level(fm, text) == ""
    finally:
        tag_articles.EVIDENCE_USE_V2 = prev


@pytest.mark.skipif(not PMID_DIR.exists(), reason="corpus absent")
def test_invivo_reagent_context_is_discounted():
    """'murine ... cell line' is an in-vitro reagent, not an animal experiment."""
    import tag_articles

    reagent_only = ("we used the murine lewis lung carcinoma cell line ll/2 and "
                    "cells were cultured in dmem")
    genuine = "tumor-bearing mice were treated and tumor volume was measured weekly"
    assert tag_articles._invivo_hits(reagent_only) == 0
    assert tag_articles._invivo_hits(genuine) > 0
