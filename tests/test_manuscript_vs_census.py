"""Guards for testing the manuscript's own claims against the census.

WHAT IT DOES
------------
Two published claims the frozen 4,830-article corpus could not test are now
testable against the whole indexed cancer literature with expert MeSH labels.
Both survive, and one is understated by the manuscript.

WHY A CONFIRMATION NEEDS GUARDING AS MUCH AS A REFUTATION
----------------------------------------------------------
The pressure on this analysis runs the opposite way to the rest of the atlas.
Everywhere else the risk is claiming too much; here the risk is a document that
exists to check the manuscript quietly becoming a document that endorses it. So
these guards pin the things that would make a confirmation hollow:

  * that the comparison is between RATIOS, not absolute counts, because the
    census is cancer-restricted and the manuscript's figures were not, so the
    absolute counts differ several-fold by construction
  * that rows too small to support a ratio are reported as unmeasurable rather
    than folded into the verdict
  * that a row whose concept has no MeSH descriptor is reported as untestable
    rather than as a zero that agrees with the manuscript
  * that the verdict FLIPS if the measurement does

The last is the one that matters. A test which can only ever print "the claim
survives" is not a test.
"""

import importlib.util
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JSON = REPO_ROOT / "analysis" / "manuscript-vs-census.json"
DOC = REPO_ROOT / "analysis" / "manuscript-vs-census.md"
SCRIPT = REPO_ROOT / "scripts" / "manuscript_vs_census.py"
MANUSCRIPT = REPO_ROOT / "article" / "drafts" / "v1.md"


def d() -> dict:
    return json.loads(JSON.read_text())


def flat() -> str:
    return " ".join(DOC.read_text().split())


def mod():
    spec = importlib.util.spec_from_file_location("manuscript_vs_census", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --- the claims must still be the manuscript's ------------------------------

def test_the_quoted_manuscript_figures_are_the_manuscript_s():
    """A claim test is worthless if it tests a claim the manuscript stopped making.

    Every figure attributed to the manuscript is checked against v1.md itself,
    so an edit there fails this rather than silently changing what is being
    tested.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    m = mod()
    for name, desc, ferro, _icd in m.MODALITIES:
        if desc is None:
            continue
        assert f"| **{ferro}** |" in txt or f"| {ferro} |" in txt, (
            f"the manuscript no longer states {ferro} for {name}; the modality "
            "table has been edited and this analysis is testing a claim that "
            "is no longer made")
    assert f"from {m.CORPUS_GROWTH_START} full-text articles" in txt, (
        "the manuscript's growth start is not the one this analysis quotes")
    assert f"to {m.CORPUS_GROWTH_END:,} in {m.GROWTH_END}" in txt, (
        "the manuscript's growth end is not the one this analysis quotes")


def test_the_verdict_would_flip_if_the_measurement_did():
    """A test that can only print "survives" is not a test.

    The first version recomputed each verdict from the same JSON fields the
    generator wrote. Because both claims genuinely hold, hardcoding the verdict
    to True was invisible: the recomputation agreed with the constant. It is
    the precise failure this test was written to prevent, committed inside it.

    The verdicts are functions now, and they are handed inputs where the answer
    MUST be False.
    """
    m = mod()
    # understatement: the census ratio must exceed the manuscript's, and the
    # comparison must be refused outright when the ratio is not measurable.
    assert m.understates(4.75, 2.93, True) is True
    assert m.understates(2.93, 4.75, True) is False, (
        "a SMALLER census ratio is still reported as understatement")
    assert m.understates(4.75, 4.75, True) is False, (
        "an equal ratio is reported as understatement")
    assert m.understates(4.75, 2.93, False) is False, (
        "an unmeasurable ratio still yields a verdict")
    assert m.understates(None, 2.93, True) is False
    # growth
    assert m.outgrew(30.71, 1.1) is True
    assert m.outgrew(1.1, 30.71) is False, (
        "a corpus growing SLOWER than the field is reported as outgrowing it")
    assert m.outgrew(1.1, 1.1) is False
    assert m.outgrew(None, 1.1) is False
    # and the floor
    assert m.ratio_is_measurable(152, 32) is True
    assert m.ratio_is_measurable(152, 1) is False, (
        "a ratio with one side below the floor is called measurable")
    assert m.ratio_is_measurable(0, 0) is False

    # the shipped verdicts must BE those functions applied to the shipped inputs
    r = d()
    mt, g = r["modality_table"], r["growth"]
    assert mt["census_exceeds_manuscript"] == m.understates(
        mt["census_pdt_sdt_ratio"], mt["manuscript_pdt_sdt_ratio"],
        mt["ratio_is_measurable"])
    assert g["corpus_exceeds_field"] == m.outgrew(
        g["corpus_growth"], g["census_growth"])


# --- what would make a confirmation hollow ----------------------------------

def test_the_absolute_counts_are_not_presented_as_agreement():
    """The census is cancer-restricted; the manuscript's queries were not.

    So absolute counts differ several-fold BY CONSTRUCTION. If the document
    ever reads them as agreement or disagreement, its central comparison is
    invalid.
    """
    r, txt = d(), flat()
    pdt = next(x for x in r["modality_table"]["rows"] if x["modality"] == "PDT")
    assert pdt["census_ferroptosis"] != pdt["manuscript_ferroptosis"], (
        "the census and manuscript counts now agree exactly, which would mean "
        "the scope difference this document is built around has gone away")
    assert "by construction" in txt.lower(), (
        "the report no longer states that the absolute counts differ by "
        "construction, so a reader may read the gap as a discrepancy")
    assert "Only ratios within a table, and directions, are comparable" in txt, (
        "the report no longer says which quantities are comparable")


def test_rows_too_small_for_a_ratio_are_reported_as_unmeasurable():
    r = d()
    mt = r["modality_table"]
    m = mod()
    assert mt["ratio_is_measurable"] == m.ratio_is_measurable(
        next(x["census_ferroptosis"] for x in mt["rows"] if x["modality"] == "PDT"),
        next(x["census_ferroptosis"] for x in mt["rows"] if x["modality"] == "SDT")), (
        "the measurability flag is not the function applied to the shipped counts")
    assert mt["rows_below_ratio_floor"], (
        "no row falls below the ratio floor any more; if every modality now "
        "carries census signal the report's unmeasurable section is wrong")
    for name in mt["rows_below_ratio_floor"]:
        row = next(x for x in mt["rows"] if x["modality"] == name)
        assert row["census_ferroptosis"] < mt["min_for_a_ratio"]
    assert not mt["icd_column_is_measurable"], (
        "the ICD column is now measurable at census scale; the report says it "
        "is not and computes no ratio from it, so that has to be revisited")
    txt = flat()
    for name in mt["rows_below_ratio_floor"]:
        assert name in txt, f"{name} is not named as unmeasurable in the report"


def test_a_row_with_no_descriptor_is_untestable_not_a_zero():
    """TTFields has no MeSH descriptor. Reporting its census count as 0 would
    read as agreeing with the manuscript's 0, which is not what happened."""
    r = d()
    un = r["modality_table"]["untestable_rows"]
    assert un, "no row is reported as untestable any more"
    for name in un:
        row = next(x for x in r["modality_table"]["rows"] if x["modality"] == name)
        assert row["census_ferroptosis"] is None, (
            f"{name} has no descriptor but ships a census count of "
            f"{row['census_ferroptosis']}, which a reader would compare with "
            "the manuscript's figure")
    assert "_no MeSH descriptor_" in flat(), (
        "the untestable row is not marked as such in the shipped table")


def test_the_growth_claim_rules_out_both_confounds_it_names():
    """The corpus was retrieved once with fixed queries, so open-access growth
    is the only year-dependent retrieval effect. Both must be reported."""
    g = d()["growth"]
    for k in ("census_growth", "open_access_growth", "corpus_growth"):
        assert g[k] and g[k] > 0, f"{k} is missing or zero"
    assert g["corpus_growth"] > g["open_access_growth"] > g["census_growth"], (
        f"the ordering has changed (corpus {g['corpus_growth']}, open access "
        f"{g['open_access_growth']}, field {g['census_growth']}); the argument "
        "is that neither confound accounts for the corpus, so a different "
        "ordering needs a different argument")
    assert g["unexplained_by_availability"] > 1.0
    assert len(g["per_year"]) == g["end_year"] - g["start_year"] + 1, (
        "the per-year table does not cover the span the growth figure spans")


def test_the_report_states_what_a_surviving_claim_does_not_establish():
    txt = flat()
    assert "Surviving a census test is not validation" in txt, (
        "the report no longer distinguishes removing a confound from showing "
        "the underlying biology is true, which is the misreading a "
        "confirmation invites")
    assert "MeSH indexing lag" in txt or "MeSH indexing lags" in txt, (
        "the report no longer states that recent years are undercounted, so "
        "every growth figure reads as exact rather than as a lower bound")
    assert "broader than" in txt, (
        "the report no longer states that two descriptors are broader than the "
        "modalities they stand for")


def test_the_headline_agrees_with_the_verdicts_it_summarises():
    """It was a STATIC STRING, and this guard used to pin it in place.

    A review forced `understates()` to return False and got a document whose
    section body read "the census ratio is smaller" while its opening line
    still read "Both claims tested here survive". The guard asserted the fixed
    words "survive" and "understated", so it held the endorsement steady
    against the measurement -- the exact failure it was written to prevent,
    landing on the one sentence a reader is most likely to read.

    It now checks AGREEMENT with the computed verdicts, and the headline is
    generated from them.
    """
    m, r = mod(), d()
    txt = DOC.read_text()
    head = txt[:txt.index("## The scope difference")]
    mt, g = r["modality_table"], r["growth"]

    headline = m._headline(r)
    assert headline in " ".join(head.split()), (
        "the shipped summary is not the sentence the verdicts generate; it has "
        "gone back to being a literal")

    # and the sentence must SAY what the verdicts hold
    survived = (mt["census_exceeds_manuscript"] or mt["direction_holds"]) \
        and g["corpus_exceeds_field"]
    assert ("survive" in headline.lower()) == bool(survived), (
        f"the headline says {headline!r} while the verdicts are "
        f"8.2={mt['direction_holds']} 3.7={g['corpus_exceeds_field']}")
    assert ("understate" in headline.lower()) == bool(
        mt["census_exceeds_manuscript"]), (
        "the headline's understatement clause does not track the verdict")

    # the generator must flip it when the verdicts flip
    flipped = m._headline({**r, "modality_table": {**mt,
                                                  "census_exceeds_manuscript": False,
                                                  "direction_holds": False},
                           "growth": {**g, "corpus_exceeds_field": False}})
    assert "survive" not in flipped.lower() or "Nothing" in flipped, (
        f"with both verdicts false the headline still reads {flipped!r}")
    assert "understate" not in flipped.lower()


def test_the_descriptor_choice_is_reported_as_a_choice():
    """A single pair ships a point estimate whose margin depends on it."""
    mt = d()["modality_table"]
    vs = mt["descriptor_variants"]
    assert len(vs) >= 3, "the sensitivity analysis has been reduced to one pair"
    assert mt["ratio_range"] and mt["ratio_range"][0] < mt["ratio_range"][1], (
        "the variants no longer span a range, so the sensitivity is not shown")
    assert mt["understatement_holds_under_every_variant"], (
        f"the understatement no longer holds across the whole range "
        f"{mt['ratio_range']}; the conclusion then depends on the descriptor "
        "pair chosen and must not be stated unconditionally")
    txt = flat()
    for v in vs:
        assert f"| {v['variant']} |" in txt, (
            f"the {v['variant']!r} variant is computed but not shipped")


def test_the_named_invalidator_is_measured_not_just_named():
    """"the ratio is only as good as their relative over-estimation" sat beside
    an unconditional verdict, unmeasured. A named invalidator that is never
    measured reads as diligence while the verdict reads as settled."""
    om = d()["modality_table"]["on_modality_and_tumour"]
    for k in ("pdt_pct", "sdt_pct", "gap_points", "filtered_ratio"):
        assert om.get(k) is not None, f"{k} is not measured"
    assert om["symmetric_within_5_points"], (
        f"the over-estimation is NOT symmetric ({om['pdt_pct']}% against "
        f"{om['sdt_pct']}%, a {om['gap_points']}-point gap), so the raw ratio "
        "is inflated and the verdict is unearned as stated")
    assert 0 < om["pdt_pct"] <= 100 and 0 < om["sdt_pct"] <= 100
    assert "symmetric" in flat(), (
        "the symmetry result is measured but not reported")
