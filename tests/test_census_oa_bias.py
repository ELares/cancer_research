"""Guards for the open-access ranking-sensitivity test (Section 3.3.1).

The analysis this replaced compared 4,830 full-text records against 5,586
abstract-only ones with the full-text side 98.7% open access -- which is not a
contrast at all, since both arms were drawn from a retrieval that had already
selected for availability. The census version splits 936,347 against 3,467,647
on identical expert descriptors, so only availability differs.

THE CONFOUND IS THE POINT OF THE GUARDS, not the ranking. PMC deposition rose
steeply over the same period the newer mechanisms grew, so a mechanism with a
recent median year has a high OA rate for reasons that have nothing to do with
its subject. If the median-year column is ever dropped, the OA column starts
reading as a fact about access, and this page becomes an attribution it has no
design to support.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-oa-bias.json"
MD = REPO / "analysis/census-oa-bias.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"
# PINNED HERE AS A LITERAL, not read from the artifact. A first version took
# the threshold from the JSON, which made every guard that used it compare the
# generator to itself: raising it from 3 to 6 emptied the "moved" set, changed
# the headline, and passed the whole file. A judgement that can be loosened
# until the result disappears is exactly the kind that needs two files to
# change.
BIG_SHIFT = 3


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def _rows(d):
    return d["mechanisms"]


def test_both_arms_partition_the_same_articles(d):
    """The comparison is only meaningful if the two arms are complements."""
    for r in _rows(d):
        assert r["with_fulltext"] + r["without"] == r["total"], (
            f"{r['mechanism']}: the two arms do not sum to its census count, so "
            "they are not a partition and the ranking contrast is meaningless")
        # Stored as a FRACTION; the render multiplies by 100. Comparing against
        # a percentage here would have passed only for a mechanism at 1%.
        assert r["oa_rate"] == pytest.approx(
            r["with_fulltext"] / r["total"], abs=0.0005)


def test_the_shift_column_is_derived_from_the_two_rankings(d):
    rows = _rows(d)
    oa = {r["mechanism"]: i for i, r in enumerate(
        sorted(rows, key=lambda r: -r["with_fulltext"]), 1)}
    non = {r["mechanism"]: i for i, r in enumerate(
        sorted(rows, key=lambda r: -r["without"]), 1)}
    for r in rows:
        assert r["rank_oa"] == oa[r["mechanism"]], r["mechanism"]
        assert r["rank_non_oa"] == non[r["mechanism"]], r["mechanism"]
        assert r["shift"] == non[r["mechanism"]] - oa[r["mechanism"]], r["mechanism"]


def test_the_headline_count_matches_the_shifts_it_counts(d):
    """An extremum or a count stated over a set nobody enumerated is how this
    repo has been wrong before.

    The moved SET is now stored in the artifact, so the count in the prose, the
    set in the JSON and the bolded rows in the table are three renderings of one
    derivation rather than three chances to disagree.
    """
    assert d["big_shift_threshold"] == BIG_SHIFT, (
        f"the generator now calls a shift 'large' at {d['big_shift_threshold']} "
        f"places where this guard pins {BIG_SHIFT}. That threshold decides how "
        "many mechanisms the headline reports as moving, so changing it must be "
        "a deliberate edit in both places rather than a quiet one here")
    rows = _rows(d)
    big = [r for r in rows if abs(r["shift"]) >= BIG_SHIFT]
    assert sorted(r["mechanism"] for r in big) == sorted(d["moved"])
    md = MD.read_text()
    assert f"**{len(big)} of {len(rows)} mechanisms shift" in md, (
        f"{len(big)} mechanisms shift {BIG_SHIFT}+ places and the report does "
        "not say so")
    for r in big:
        assert f"**{r['shift']:+d}**" in md, (
            f"{r['mechanism']} shifts {r['shift']} and is not marked in the table")


def test_the_era_confound_is_printed_beside_the_oa_column(d):
    """Read alone, the OA column reads as a fact about access. It is not one.

    Every row must carry a median year, and the report must say the confound is
    uncontrolled -- otherwise a reader takes an association for a cause.
    """
    for r in _rows(d):
        assert isinstance(r.get("median_year"), (int, float)), (
            f"{r['mechanism']} has no median year, so its OA rate has nothing "
            "to be read against")
    md = MD.read_text()
    assert "not controlled here" in md
    assert "median year" in md.lower()
    # The worked example: the lowest-OA mechanism should also be among the
    # oldest, which is what makes the confound concrete rather than abstract.
    rows = _rows(d)
    lowest = min(rows, key=lambda r: r["oa_rate"])
    median_years = sorted(r["median_year"] for r in rows)
    assert lowest["median_year"] <= sorted(median_years)[len(median_years) // 2], (
        "the lowest-OA mechanism is no longer among the older half, so the "
        "report's worked example of the confound no longer holds -- check the "
        "prose before trusting it")


def test_it_does_not_claim_the_ranking_is_unaffected(d):
    md = MD.read_text()
    for overclaim in ("ranking is stable", "no effect on the ranking",
                      "availability does not matter"):
        assert overclaim not in md.lower()


def test_the_manuscript_quotes_this_analysis_not_the_superseded_one():
    txt = " ".join(MANUSCRIPT.read_text().split())
    d = json.loads(JSON.read_text())
    rows = _rows(d)
    big = [r for r in rows if abs(r["shift"]) >= BIG_SHIFT]
    assert "census-oa-bias.md" in txt, (
        "Section 3.3.1 does not cite the census version of this test")
    assert f"{len(big)} of {len(rows)} mechanisms shift" in txt, (
        "the manuscript's rank-shift count is not the one the artifact derives")
    # The superseded contrast may still be NAMED in a retraction; what it must
    # not do is stand as the section's measurement. So the guard requires the
    # census figures to be the ones quoted, rather than banning a string that a
    # legitimate retraction also contains.
    assert f"{d['with_pmcid']:,}" in txt, (
        "Section 3.3.1 does not quote the census availability split, so it may "
        "still be resting on the superseded 4,830-vs-5,586 contrast whose two "
        "arms were both drawn from a retrieval that selected for availability")
