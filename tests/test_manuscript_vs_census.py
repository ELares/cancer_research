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


def test_the_shipped_document_is_the_generator_s_output():
    """Every other guard here is a substring check, which defends nothing.

    A review rewrote the verdict to "survives decisively", DELETED the caveat
    the newest section is built on, moved the central narrowing paragraph into
    an appendix below the limitations, and weakened "the ratio survives it" to
    "the ratio is unaffected by it" -- and the whole suite passed. The
    substring guards survived because later sections happened to supply
    duplicate matches for the phrases that were their only defence.

    Re-rendering from the shipped JSON and comparing byte-for-byte means prose
    can only change by changing the generator, and every other guard here then
    only has to defend the generator's text.
    """
    m = mod()
    assert m.render(d()) == DOC.read_text(), (
        "the shipped markdown is not what the generator produces from the "
        "shipped JSON, so the prose has been edited by hand or the two are "
        "out of step; regenerate rather than editing the document")


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
    # AGREEMENT, for these two as well. The commit that converted the variant
    # and symmetry guards to agreement left these asserting the favourable
    # outcome -- and ADDED the render() branches for the unfavourable case,
    # making them unreachable. A census in which the ICD column became
    # measurable would produce a correct, self-consistent document and a RED
    # suite.
    txt = flat()
    icd_yes = "The ICD column IS now measurable" in txt
    icd_no = "The ICD column is not measurable" in txt
    assert icd_yes != icd_no, "the report states both or neither ICD branch"
    assert icd_yes == mt["icd_column_is_measurable"], (
        f"the report says the ICD column "
        f"{'IS' if icd_yes else 'is not'} measurable while the measurement "
        f"says {mt['icd_column_is_measurable']}")
    # the floor list must be exactly the rows that are under it
    under = [x["modality"] for x in mt["rows"]
             if x["testable"] and x["census_ferroptosis"] < mt["min_for_a_ratio"]]
    assert mt["rows_below_ratio_floor"] == under, (
        f"the reported floor list {mt['rows_below_ratio_floor']} is not the "
        f"rows actually under the floor {under}")
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


def test_the_base_fragility_is_reported_and_the_finding_survives_it():
    """A 31-fold ratio off a base of 38 is fragile as a FIGURE.

    Stated so the number is not read as precise, and measured so the reader can
    see the fragility does not reach the finding. If it ever does, the
    attribution in section 3.7 needs revisiting rather than this guard
    relaxing.
    """
    g = d()["growth"]
    bs = g["base_sensitivity"]
    assert len(bs) >= 4, "the sensitivity band has been reduced to a point"
    assert any(b["base"] != g["corpus_start"] for b in bs), (
        "every row uses the shipped base, so nothing is being varied")
    assert all(b["still_outgrows_field"] for b in bs), (
        "the corpus stops outgrowing the field somewhere in the base band, so "
        "the fragility reaches the finding and section 3.7's attribution has "
        "to be revisited")
    # the band must actually span, or "fragile" is being asserted not shown
    ratios = [b["ratio"] for b in bs]
    assert max(ratios) - min(ratios) > 1.0, (
        "the sensitivity band is flat, so the report calls the figure fragile "
        "while showing nothing that moves")
    txt = flat()
    assert "should not be read as precise" in txt
    # The break-even, which a chosen band cannot be accused of flattering.
    assert g["break_even_base"] and g["break_even_base"] > g["corpus_start"], (
        "the break-even base is not above the actual base, which would mean "
        "the conclusion already fails")
    assert g["break_even_multiple_of_actual"] > 2.0, (
        f"the conclusion breaks at only "
        f"{g['break_even_multiple_of_actual']}x the actual base, so it is "
        "sensitive enough that the band should not be presented as reassuring")
    assert f"{g['break_even_base']:,}" in txt, (
        "the break-even base is computed but not reported")
    # and it must BE the break-even: at that base the corpus must no longer win
    m = mod()
    assert not m.outgrew(g["corpus_end"] / g["break_even_base"],
                         g["census_growth"]), (
        "at the reported break-even base the corpus still outgrows the field, "
        "so the figure is not the break-even it claims to be")
    assert m.outgrew(g["corpus_end"] / (g["break_even_base"] - 1),
                     g["census_growth"]), (
        "one article below the reported break-even the corpus already fails, "
        "so the break-even is off by more than rounding")


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
    # PER CLAIM. Equating "the word survive appears" with "both held" raised a
    # false alarm on a state the generator handles correctly (8.2 failed, 3.7
    # held), and the likely repair under CI pressure is to loosen the guard.
    for name, ok in (("section 8.2", mt["census_exceeds_manuscript"]
                      or mt["direction_holds"]),
                     ("section 3.7", g["corpus_exceeds_field"])):
        says_failed = f"{name} does not survive" in headline.lower()
        assert says_failed == (not ok), (
            f"the headline {'says' if says_failed else 'does not say'} "
            f"{name} failed while its verdict is {ok}: {headline!r}")
    assert ("understated by the manuscript" in headline.lower()) == bool(
        mt["census_exceeds_manuscript"]), (
        "the headline's understatement clause does not track the verdict")

    # EVERY state, and through render() so the headline is checked against the
    # BODY it must agree with. Driving only the shipped state and an all-false
    # fixture left the commit's central fix untested: swapping the branch order
    # back -- the original defect verbatim -- kept the suite green, because
    # neither of the two states that reach the changed branch was ever driven.
    for meas in (True, False):
        for exceeds in (True, False):
            for direction in (True, False):
                for grew in (True, False):
                    if exceeds and not direction:
                        continue          # cannot exceed without holding
                    v = {**r,
                         "modality_table": {**mt, "ratio_is_measurable": meas,
                                            "census_exceeds_manuscript": meas and exceeds,
                                            "direction_holds": meas and direction},
                         "growth": {**g, "corpus_exceeds_field": grew}}
                    doc = m.render(v).lower()
                    head = doc[:doc.index("## the scope difference")]
                    label = f"meas={meas} exceeds={exceeds} dir={direction} grew={grew}"
                    if not meas:
                        assert "cannot be decided" in head, (
                            f"[{label}] an undecidable ratio is not reported as "
                            f"undecidable in the headline: {head.strip()[:160]!r}")
                        assert "section 8.2 does not survive" not in head, (
                            f"[{label}] undecidable is reported as refuted, "
                            "which the document says three times it does not do")
                    else:
                        assert (("section 8.2 does not survive" in head)
                                == (not direction)), (
                            f"[{label}] the headline's 8.2 clause disagrees "
                            f"with its verdict: {head.strip()[:160]!r}")
                        assert (("understated by the manuscript" in head)
                                == bool(exceeds)), (
                            f"[{label}] the understatement clause does not "
                            "track the verdict")
                    assert (("section 3.7 does not survive" in head)
                            == (not grew)), (
                        f"[{label}] the headline's 3.7 clause disagrees with "
                        "its verdict")


def test_the_descriptor_choice_is_reported_as_a_choice():
    """A single pair ships a point estimate whose margin depends on it."""
    mt = d()["modality_table"]
    vs = mt["descriptor_variants"]
    assert len(vs) >= 3, "the sensitivity analysis has been reduced to one pair"
    assert mt["ratio_range"] and mt["ratio_range"][0] < mt["ratio_range"][1], (
        "the variants no longer span a range, so the sensitivity is not shown")
    # AGREEMENT, not the favourable outcome. Asserting the flag is true is the
    # pinning shape the headline guard was rewritten to stop doing, and render()
    # already carries a correct branch for the unfavourable case -- which an
    # assertion on the flag makes unreachable, so a census that refuted the
    # claim would surface as a broken test rather than as a document reporting
    # the refutation.
    txt_ = flat()
    favourable = "Every variant exceeds the manuscript's" in txt_
    unfavourable = "Not every variant exceeds the manuscript's" in txt_
    assert favourable != unfavourable, (
        "the document states both or neither branch of the variant verdict")
    assert favourable == mt["understatement_holds_under_every_variant"], (
        f"the document says the understatement "
        f"{'holds' if favourable else 'does not hold'} across the range while "
        f"the measurement says "
        f"{mt['understatement_holds_under_every_variant']}")
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
    txt_ = flat()
    sym = "the over-estimation is symmetric" in txt_
    asym = "one descriptor is materially broader" in txt_
    assert sym != asym, "the document states both or neither symmetry branch"
    assert sym == om["symmetric_within_5_points"], (
        f"the document says the over-estimation is "
        f"{'symmetric' if sym else 'asymmetric'} while the measurement "
        f"({om['pdt_pct']}% against {om['sdt_pct']}%, a "
        f"{om['gap_points']}-point gap) says "
        f"{om['symmetric_within_5_points']}")
    assert 0 < om["pdt_pct"] <= 100 and 0 < om["sdt_pct"] <= 100
    assert "symmetric" in flat(), (
        "the symmetry result is measured but not reported")
