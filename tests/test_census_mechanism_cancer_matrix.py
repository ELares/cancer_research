"""Guards for the mechanism-by-site matrix.

The analysis exists to replace a measure, so the failure that matters is
quietly keeping the old one: reporting that zero cells fell from 22.5% to 2.1%
reads as good news about the field, when what it establishes is that COUNTING
ZEROS NO LONGER WORKS. A gap measure returning 2% on the full literature is not
detecting gaps.

Two further failures are specific to what replaces it. A ratio against expected
invites a p-value, and the null here -- that a mechanism's literature ignores
which site it is -- is known false before any data is read, so a test against it
would return significance nearly everywhere and answer a question nobody asked.
And a depleted cell invites being called an opportunity, which the counts
cannot support: an ablative modality is depleted in disseminated disease
because there is nothing to ablate.

OFFLINE: these read only the committed artifact.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-mechanism-cancer-matrix.json"
MD = REPO / "analysis/census-mechanism-cancer-matrix.md"
SCRIPT = REPO / "scripts/census_mechanism_cancer_matrix.py"
# Pinned literals, not read from the artifact: a threshold taken from the file
# the generator wrote compares the generator to itself.
MIN_EXPECTED = 20.0
MIN_CELLS_FOR_SPREAD = 8


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_universe_is_articles_carrying_both_axes(d):
    """An expectation computed against a population most of the matrix could
    never enter is not an expectation.

    If the universe were the whole census, every expected value would be a
    fraction of a population that has no mechanism tag and no site, and every
    observed cell would look enormously enriched.
    """
    assert d["universe"] < d["census"] / 10, (
        "the universe is close to the whole census, so it is probably not "
        "restricted to articles carrying BOTH a mechanism and a site")
    src = SCRIPT.read_text()
    assert "THE UNIVERSE IS ARTICLES CARRYING BOTH" in src
    # Marginals must be counted over the same universe, or the expectation
    # mixes two populations.
    for m, n in d["mechanism_totals"].items():
        assert n <= d["universe"], f"{m} exceeds the universe it is counted in"
    for s, n in d["site_totals"].items():
        assert n <= d["universe"], f"{s} exceeds the universe it is counted in"


def test_expected_and_ratio_recompute_from_the_marginals(d):
    U = d["universe"]
    for r in d["rows"]:
        exp = d["mechanism_totals"][r["mechanism"]] * d["site_totals"][r["site"]] / U
        assert r["expected"] == pytest.approx(exp, abs=0.15), r
        if r["ratio"] is not None:
            assert r["ratio"] == pytest.approx(r["observed"] / exp, abs=0.02)
        assert r["interpretable"] == (exp >= MIN_EXPECTED)


def test_the_zero_count_is_reported_as_a_retirement_not_a_result(d):
    """The trap this whole analysis is built around.

    "Empty cells fell from 22.5% to 2.1%" reads as a finding about the field.
    It is a finding about the MEASURE: at census scale nearly every pair has
    been written about once, so zeros no longer separate a neglected
    combination from a well-studied one.
    """
    zeros = [r for r in d["rows"] if r["observed"] == 0]
    assert d["n_zero"] == len(zeros)
    assert sorted(d["zero_cells"]) == sorted(
        f"{r['mechanism']}/{r['site']}" for r in zeros)
    assert d["census_zero_share"] == pytest.approx(
        100 * len(zeros) / len(d["rows"]), abs=0.1)
    md = MD.read_text()
    assert "RETIRES the measure" in md or "RETIRES ZERO AS A MEASURE" in md
    assert "is not detecting gaps" in md, (
        "the report gives the falling zero count without saying the measure "
        "has stopped working, which is the finding")


def test_no_significance_test_is_computed_or_implied(d):
    md = MD.read_text()
    assert "## No p-value, deliberately" in md
    assert "known false" in md.lower()
    for banned in ("p =", "p-value of", "p <", "chi-square", "statistically significant"):
        assert banned not in md.lower().replace("no p-value", ""), (
            f"the report reports {banned!r} against a null that is false a "
            "priori, which tests a hypothesis nobody holds")


def test_no_depleted_cell_is_named_an_opportunity(d):
    """Separating "nobody has tried this" from "this cannot work here" needs
    knowledge of the mechanism, and the counts do not carry it."""
    md = MD.read_text()
    assert "A depleted cell is not a gap" in md
    assert "working as designed" in md
    for overclaim in ("untapped", "an opportunity", "should be studied",
                      "promising gap", "neglected combination that"):
        assert overclaim not in md.lower(), (
            f"the report names a cell an opportunity ({overclaim!r}); the "
            "counts cannot distinguish neglect from impossibility")


def test_the_size_control_is_derived_and_its_verdict_follows_the_number(d):
    """Without it, the tails could be small mechanisms with noisy ratios and
    the ranking would be a property of the denominator."""
    import math

    r = d["size_vs_spread_r"]
    if r is None:
        pytest.skip("too few mechanisms for the control")
    xs = [math.log10(x["articles"]) for x in d["spread"]]
    ys = [x["spread"] for x in d["spread"]]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    assert r == pytest.approx(num / den, abs=0.02)
    md = MD.read_text()
    assert f"**{r:+.2f}**" in md
    # The verdict must track the number rather than standing beside it.
    if abs(r) < 0.4:
        assert "not a small-number artifact" in md
    else:
        assert "partly a measurement of mechanism size" in md


def test_a_mechanism_in_both_tails_is_flagged_as_one_fact(d):
    """A row's ratios are constrained: concentrating a fixed literature into a
    few sites FORCES depletion in the rest. Reading the two tables as
    independent evidence double-counts a single pattern."""
    dep = {r["mechanism"] for r in d["most_depleted"]}
    enr = {r["mechanism"] for r in d["most_enriched"]}
    assert sorted(d["in_both_tails"]) == sorted(dep & enr)
    if d["in_both_tails"]:
        md = MD.read_text()
        assert "one fact rather than two" in md
        assert "arithmetic shadow" in md
        for m in d["in_both_tails"]:
            assert f"`{m}`" in md


def test_spread_needs_enough_cells_to_have_a_shape(d):
    """A standard deviation over two numbers is not a spread. Microbiome
    clears the expectation floor in 2 of 18 sites."""
    counted = {r["mechanism"] for r in d["spread"]}
    for r in d["spread"]:
        assert r["cells"] >= MIN_CELLS_FOR_SPREAD
    thin = {m for m in d["mechanism_totals"]
            if len([x for x in d["rows"]
                    if x["mechanism"] == m and x["interpretable"]])
            < MIN_CELLS_FOR_SPREAD}
    assert not (counted & thin), (
        f"{sorted(counted & thin)} have too few interpretable cells for a "
        "spread and are ranked anyway")


def test_mechanisms_mesh_cannot_express_are_absent_not_zero(d):
    """The failure a falling zero count would otherwise hide.

    TTFields and bioelectric modulation have no descriptor, so they cannot
    appear in this matrix at all -- and an absent row contributes no zeros,
    which means the headline "only 6 empty cells" is silent about them.
    """
    for absent in ("ttfields", "bioelectric", "cold-atmospheric-plasma"):
        assert absent not in d["mechanism_totals"], (
            f"{absent} has no MeSH descriptor and must not appear as a row")
    md = MD.read_text()
    assert "ABSENT rather than empty" in md, (
        "the report reports the zero count falling without saying that the "
        "mechanisms MeSH cannot express contribute no zeros to it")
