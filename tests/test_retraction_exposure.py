"""Guards for the retraction-exposure measurement.

WHAT IT MEASURES
----------------
Retracted papers stay indexed and machine-extracted relation graphs keep the
assertions they contributed, so anyone building on a PubTator-derived graph
inherits them silently. This bounds that for the cancer census using PubMed's own
`Retracted Publication` type, which this repository has stored per record all
along and never read.

THE TWO TRAPS THIS GUARDS, both of which the project has fallen into before
------------------------------------------------------------------------
1. CONTAINMENT REPORTED AS DAMAGE. `atlas_ambiguity_impact.py` found about half of
   rows touching a contested identifier against a far smaller share actually at risk -- a
   39-fold gap. The same shape applies here, so the report must carry both.
2. THE WRONG DENOMINATOR. 70.2% of the graph's pairs rest on a single paper,
   which is the population that a retraction can actually strand. Quoting the
   rate against all pairs understates it; quoting it against single-paper pairs
   alone hides that the 70% has nothing to do with retraction.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JSON = REPO_ROOT / "analysis" / "atlas-retraction-exposure.json"
DOC = REPO_ROOT / "analysis" / "atlas-retraction-exposure.md"
SCRIPT = REPO_ROOT / "scripts" / "atlas_retraction_exposure.py"


def d() -> dict:
    return json.loads(JSON.read_text())


def flat() -> str:
    """Document text with whitespace collapsed.

    Content assertions must not depend on where a line happens to wrap. The
    predicate-skew guard first looked for "hypothesis rather than a result" and
    failed on prose that reads exactly that way but wraps mid-phrase -- a true
    claim reported as missing, which is the failure that trains people to
    weaken guards. Same fix as tests/test_preregistration_ordering.py.
    """
    return re.sub(r"\s+", " ", DOC.read_text())


def test_the_census_denominator_matches_the_repos_own_figure():
    """4,403,994 is the census size every other analysis uses."""
    assert d()["census_records"] == 4_403_994, (
        f"census total is {d()['census_records']:,}; every other analysis in the "
        "repo uses 4,403,994, so either this scan or those are wrong")


def test_containment_and_damage_are_both_reported():
    r = d()
    assert r["pairs_touched"] > r["pairs_uncorroborated"], (
        "touched should exceed uncorroborated; if they are equal the "
        "corroboration check is not doing anything")
    txt = DOC.read_text()
    assert f"{r['pairs_touched']:,}" in txt and f"{r['pairs_uncorroborated']:,}" in txt, (
        "the document does not carry both the containment and the damage figure")
    assert "touched" in txt and "only" in txt.lower()


def test_the_single_paper_denominator_is_reported():
    """The population a retraction can actually strand."""
    r = d()
    assert r["single_paper_pairs"] > 0
    frac = r["single_paper_pairs"] / r["pairs"]
    assert 0.5 < frac < 0.95, (
        f"single-paper share is {frac:.1%}; the document's framing assumes it is "
        "the large majority of the graph")
    txt = DOC.read_text()
    assert f"{r['single_paper_pairs']:,}" in txt, (
        "the at-risk denominator is not in the document, so the uncorroborated "
        "rate is quoted against all pairs only")
    both = f"{100.0*r['pairs_uncorroborated']/r['single_paper_pairs']:.3f}%"
    assert both in txt, "the rate against the single-paper denominator is missing"


def test_only_the_retracted_article_is_treated_as_tainted():
    """A retraction notice or erratum is a legitimate record."""
    src = SCRIPT.read_text()
    assert 'TAINTED = "Retracted Publication"' in src
    fam = d()["retraction_family_types"]
    assert "Retracted Publication" in fam
    # the notice must be SEEN but not counted as tainted
    assert any("Notice" in k or "Erratum" in k for k in fam), (
        "no retraction notice or erratum type was found at all, which is what the "
        "first version of this script reported when it guessed the type string")


def test_the_family_types_are_discovered_not_guessed():
    """The zero that was a wrong constant.

    The first version looked for `Retraction of Publication` -- a plausible
    string PubMed does not use -- and reported ZERO notices across 4.4M records.
    A zero of that shape reads as a finding. Types are now matched by pattern.
    """
    src = SCRIPT.read_text()
    # The ASSIGNMENT and the USE, not the bare name. Renaming the constant left
    # "FAMILY" in its own explanatory comment and at the use site, so a bare
    # membership check passed a mutation that removed the discovery entirely --
    # the same shape as a guard satisfied by its own description.
    assert "FAMILY = (" in src, "the retraction-family pattern tuple is gone"
    assert "in FAMILY" in src, (
        "FAMILY is defined but never consulted, so types are no longer discovered")
    assert "retract" in src
    assert 'NOTICE = "Retraction of Publication"' not in src, (
        "the guessed type string is back")
    assert len(d()["retraction_family_types"]) >= 3, (
        "fewer retraction-family types than expected; the discovery pattern may "
        "have narrowed")


def test_the_report_states_it_is_a_lower_bound():
    txt = DOC.read_text().lower()
    # In the LIMITATIONS section specifically. The phrase appears twice, so a
    # whole-document check survived deleting the one beside the headline figure.
    lim = txt[txt.index("what this cannot say"):] if "what this cannot say" in txt else ""
    assert lim, "the report has no limitations section"
    assert "lower bound" in lim, (
        "the limitations section does not say retraction indexing lags, so every "
        "count is a lower bound")
    assert "not necessarily false" in txt, (
        "the report no longer distinguishes withdrawn SUPPORT from falsity")
    assert "extractor" in txt or "pubtator" in txt, (
        "the report does not bound itself by the extractor's own error rate")


def test_the_predicate_skew_is_offered_as_a_hypothesis():
    """Small counts, uncontrolled. It must not read as a result."""
    txt = flat()
    if "Which predicates carry it" in txt:
        assert "hypothesis rather than a result" in txt, (
            "the predicate skew is stated without marking it as uncontrolled")
