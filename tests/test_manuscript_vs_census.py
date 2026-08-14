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


def test_the_confirmation_is_reported_as_prominently_as_a_failure():
    """The specific risk of a self-check: burying the good news, or the bad."""
    txt = DOC.read_text()
    head = txt[:txt.index("## The scope difference")]
    assert "survive" in head, (
        "the summary does not say the claims survived; a reader has to reach "
        "the tables to learn the outcome")
    assert "understated" in head, (
        "the summary does not say the manuscript understated one of its own "
        "claims, which is the more surprising half of the result")
