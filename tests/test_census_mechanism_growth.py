"""Guards for the census growth analysis and the manuscript sentences on it.

THE DEFECT THIS FILE EXISTS FOR was found by reviewing my own draft rather than
by any test: the manuscript called electrochemical therapy "the one mechanism
whose literature is measurably shrinking", an extremum asserted over a set
nobody had enumerated. Epigenetic work is smaller too. The set is derived in
the generator now, and these guards keep the prose and the derivation together.

The second concern is the denominator. This analysis exists because
`manuscript_vs_census` compared a retrieved corpus against the WHOLE cancer
literature, which a corpus of emerging-therapy papers outgrows whether or not
anything unusual happened. Reverting to that comparison would restore a
conclusion the arithmetic does not support, so the matched row is pinned.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-mechanism-growth.json"
MD = REPO / "analysis/census-mechanism-growth.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_shrinking_set_is_derived_not_described(d):
    """Recompute it from the rows rather than trusting the stored field."""
    expected = sorted(r["mechanism"] for r in d["rows"]
                      if r["growth"] is not None and r["growth"] < 1.0)
    assert sorted(d["shrinking"]) == expected
    expected_recent = sorted(r["mechanism"] for r in d["rows"]
                             if r["recent_pct"] is not None and r["recent_pct"] < 0)
    assert sorted(d["recent_shrinking"]) == expected_recent


def test_the_report_states_how_many_shrink(d):
    md = MD.read_text()
    n = len(d["shrinking"])
    assert f"{n} mechanism(s) are SMALLER" in md
    for m in d["shrinking"]:
        assert f"`{m}`" in md, f"{m} shrinks and the report does not name it"


def test_the_manuscript_does_not_claim_a_lone_shrinking_mechanism(d):
    """The exact wording that was wrong, banned; and the correct count required.

    Banning the phrase alone would be a guard satisfied by a reword, so the
    manuscript must also state a count that MATCHES the derived set.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    assert "the one mechanism in this book whose literature is measurably shrinking" not in txt
    n = len(d["shrinking"])
    words = {1: "one", 2: "two", 3: "three", 4: "four"}
    assert f"only {words[n]} mechanisms" in txt or f"only {words[n]} mechanism" in txt, (
        f"{n} mechanisms shrink and the manuscript does not say so; it must "
        "not describe the set without counting it")
    for m in d["shrinking"]:
        stem = m.split("-")[0]
        assert stem in txt.lower(), f"{m} shrinks and the manuscript does not name it"


def test_the_recent_window_claim_names_its_base_threshold(d):
    """A "fastest-growing" claim off a base of 24 needs its qualifier.

    Without it, an even smaller mechanism takes the title on a handful of
    papers -- the extremum-over-a-self-selected-set defect one step along.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    if "fastest-growing" in txt:
        assert "at least 20 articles" in txt, (
            "the manuscript claims a fastest-growing mechanism without stating "
            "the base threshold that decides which mechanisms were eligible")


def test_the_matched_denominator_is_present_and_beats_the_field(d):
    """The row this analysis exists to supply, and the comparison it replaces."""
    assert d["union_growth"] and d["field_growth"]
    assert d["mechanisms_over_field"] == pytest.approx(
        d["union_growth"] / d["field_growth"], abs=0.01)
    md = MD.read_text()
    assert "denominator a growth claim needs" in md
    assert f"x{d['union_growth']}" in md and f"x{d['field_growth']}" in md


def test_a_growth_ratio_is_withheld_below_the_base_floor(d):
    """A ratio off a handful of articles measures the handful."""
    for r in d["rows"]:
        if r["start"] < d["min_base"]:
            assert r["growth"] is None, (
                f"{r['mechanism']} starts at {r['start']} and still carries a "
                f"ratio; below {d['min_base']} the ratio measures the base")


def test_render_only_reassembles_rather_than_trusting_stored_fields():
    """Adding a derived column must not strand --render-only on an old artifact.

    It also means the stored derived fields are checked against a fresh
    derivation on every render instead of being carried forward untested.
    """
    src = (REPO / "scripts/census_mechanism_growth.py").read_text()
    block = src[src.index("if a.render_only:"):src.index("    else:", src.index("if a.render_only:"))]
    assert "assemble(" in block, (
        "--render-only reads the stored JSON without re-deriving it, so a new "
        "column will crash it against an artifact written before that column")
