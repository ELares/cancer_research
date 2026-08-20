"""Guards for the external check of the census build.

This is the only artifact in the repo that compares the census against a source
outside it, so its failure mode is specific: a check that agrees for the WRONG
REASON certifies a build nobody has actually verified. Three query properties
do all the work, and getting any one wrong produces a comfortable-looking
number:

* Explosion. `X[mh]` includes the whole MeSH subtree while the census matches
  DescriptorName exactly. Left exploded, `Ultrasonic Therapy` returns 4,371
  against the census's 2,513 and silently includes the HIFU records this
  project counts separately -- a 74% "disagreement" that is entirely the query.
* The baseline date. PubMed keeps growing past the snapshot, so an uncapped
  query makes the census read low on every recent mechanism.
* The cancer restriction. The census admits nine adjacent descriptors that
  `neoplasms[mh]` cannot return, so the comparison must use the C04 core.

OFFLINE: these read only the committed artifact and never touch the network.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-external-check.json"
MD = REPO / "analysis/census-external-check.md"
SCRIPT = REPO / "scripts/census_external_check.py"
# Pinned as literals, NOT read from the artifact. A threshold taken from the
# file the generator wrote compares the generator to itself -- the defect a
# mutation sweep found in the sibling OA analysis, where raising a cut-point
# emptied the reported set and passed every test.
FLAG_AT = 0.10
BASELINE_DATE = "2026/06/17"


@pytest.fixture(scope="module")
def d():
    if not JSON.exists():
        pytest.skip("external check has not been run (needs network)")
    return json.loads(JSON.read_text())


def test_every_query_is_unexploded_dated_and_cancer_restricted(d):
    """The three properties that make the comparison mean anything."""
    for r in d["rows"]:
        q = r["query"]
        if q is None:
            # A mechanism with no descriptor has no query to check. It is NOT
            # queried, deliberately: building one from an empty list yields
            # `()`, which PubMed answers, and that answer would enter the table
            # as a disagreement.
            assert r["mechanism"] in d["not_comparable"]
            continue
        assert "[mh:noexp]" in q, (
            f"{r['mechanism']}'s query explodes the MeSH tree, so it counts a "
            "subtree the census never matched")
        assert "[mh]" not in q.replace("[mh:noexp]", "").replace(
            "neoplasms[mh]", ""), (
            f"{r['mechanism']}'s query has an unqualified [mh] clause")
        assert "neoplasms[mh]" in q, (
            f"{r['mechanism']}'s query is not cancer-restricted")
        assert BASELINE_DATE in q, (
            f"{r['mechanism']}'s query is not capped at the baseline, so "
            "PubMed's growth since the snapshot reads as a build defect")


def test_the_comparison_uses_the_c04_core_not_the_whole_census(d):
    """The adjacent extension admits records with NO C04 descriptor, which
    `neoplasms[mh]` cannot return. Including them would guarantee the census
    reads high and the gap would be a property of the comparison."""
    assert d["c04_core"] < d["census_records"], (
        "the C04 core is not smaller than the census, so the adjacent "
        "extension is being folded into a comparison that cannot see it")
    src = SCRIPT.read_text()
    assert 'r.get("cancer_basis") != "C04"' in src


def test_gaps_and_flags_recompute_from_the_two_counts(d):
    for r in d["rows"]:
        if r["pubmed"] is None:
            assert r["ratio"] is None and r["rel_gap"] is None
            assert not r["flagged"]
            continue
        assert r["ratio"] == pytest.approx(r["census"] / r["pubmed"], abs=0.002)
        assert r["rel_gap"] == pytest.approx(
            abs(r["census"] - r["pubmed"]) / r["pubmed"], abs=0.002)
        assert r["flagged"] == (r["rel_gap"] > FLAG_AT)
    assert sorted(d["flagged"]) == sorted(
        r["mechanism"] for r in d["rows"] if r["flagged"])


def test_the_flag_threshold_is_the_pinned_one(d):
    assert d["flag_threshold"] == FLAG_AT, (
        f"the generator flags at {d['flag_threshold']} where this guard pins "
        f"{FLAG_AT}. That threshold decides how many mechanisms are reported as "
        "disagreeing, so loosening it must take a deliberate edit in two files")


def test_the_direction_verdict_follows_the_counts(d):
    """A scatter around zero is noise; a one-sided gap is systematic and needs
    an explanation rather than a tolerance. The report must not call a
    one-sided gap 'balanced'."""
    higher = sum(1 for r in d["rows"] if r["ratio"] and r["ratio"] > 1)
    lower = sum(1 for r in d["rows"] if r["ratio"] and r["ratio"] <= 1)
    assert d["census_higher"] == higher
    assert d["census_lower"] == lower
    md = MD.read_text()
    balanced = abs(higher - lower) <= 2
    if balanced:
        assert "close to balanced" in md
        assert "ONE-SIDED" not in md
    else:
        assert "ONE-SIDED" in md, (
            f"the census reads higher on {higher} and lower on {lower}, which "
            "noise does not produce, and the report calls it balanced")


def test_a_failed_request_is_not_pooled_with_a_missing_descriptor(d):
    """Two different reasons a row has no count, and pooling them would let a
    network outage read as a taxonomy gap -- or the reverse, which is worse: a
    mechanism MeSH cannot express reported as a request that happened to fail,
    and therefore as something a retry would fix.
    """
    no_query = {r["mechanism"] for r in d["rows"] if r["query"] is None}
    failed = {r["mechanism"] for r in d["rows"]
              if r["pubmed"] is None and r["query"] is not None}
    assert set(d["not_comparable"]) == no_query
    assert set(d["unresolved"]) == failed
    assert not (no_query & failed)
    assert d["compared"] == len(d["rows"]) - len(no_query) - len(failed), (
        "the compared count does not account for both kinds of missing row, so "
        "a failure could shrink the denominator into a better-looking median")
    md = MD.read_text()
    for m in no_query:
        assert f"| {m} |" in md and "*not comparable*" in md
    for m in failed:
        assert f"| {m} |" in md and "*unresolved*" in md


def test_the_one_sided_gap_is_explained_and_the_explanation_is_tested(d):
    """A one-sided gap is the finding this check exists to surface.

    Leaving it flagged but unexplained would be the honest half of the job; the
    other half is a prediction that could have failed. The candidate cause --
    MeSH indexing applied after the snapshot to records that entered PubMed
    before it -- predicts the shortfall grows with a mechanism's recency, and
    the correlation must be computed rather than asserted.
    """
    balanced = abs(d["census_higher"] - d["census_lower"]) <= 2
    if balanced:
        pytest.skip("the gap is balanced, so no systematic cause is claimed")
    rt = d["recency_test"]
    assert rt["prediction"].startswith("positive"), (
        "the prediction must be stated in the artifact so it can be read "
        "against the result rather than fitted to it")
    rho = rt["spearman_year_vs_gap"]
    assert rho is not None and rt["n"] >= 4
    assert rt["supported"] == (rho >= 0.5)
    md = MD.read_text()
    assert f"**{rho:+.2f}**" in md
    if rt["supported"]:
        assert "the prediction holds" in md
        # The consequence a reader needs: a bound on growth claims, not a
        # correction to apply.
        assert "LOWER bound" in md
    else:
        assert "the prediction FAILS" in md, (
            "the recency prediction failed, so the one-sided gap is "
            "unexplained and the report must say so rather than moving on")


def test_the_recency_correlation_recomputes(d):
    """Derived, not stored-and-trusted."""
    rt = d["recency_test"]
    if rt["spearman_year_vs_gap"] is None:
        pytest.skip("no correlation computed")
    pairs = [(r["median_year"], r["rel_gap"]) for r in d["rows"]
             if r.get("median_year") and r.get("rel_gap") is not None]
    assert len(pairs) == rt["n"]
    # Recompute with an INDEPENDENT implementation rather than importing the
    # generator's, which would compare it to itself.
    def rank(vals):
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        out = [0.0] * len(vals)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            for k in range(i, j + 1):
                out[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return out
    xs, ys = rank([p[0] for p in pairs]), rank([p[1] for p in pairs])
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    assert rt["spearman_year_vs_gap"] == pytest.approx(num / den, abs=0.02)


def test_the_manuscript_carries_the_measured_bound():
    """The growth section says recent years are under-counted. It must say by
    how much, or the caveat is unactionable."""
    md = (REPO / "article/drafts/v1.md").read_text()
    if "lower bound" not in md.lower():
        pytest.skip("the manuscript no longer states the growth bound")
    d = json.loads(JSON.read_text())
    assert f"{100 * d['median_rel_gap']:.1f}%" in md, (
        "the manuscript states the direction of the under-count but not the "
        "measured size")
    assert "census-external-check.md" in md


def test_it_does_not_claim_the_descriptors_mean_what_we_take_them_to_mean():
    """Agreement on counts is agreement on ADMISSION. The breadth problem --
    whether `Plasma Gases` means cold atmospheric plasma -- is untouched by it,
    and a reader who takes this page for a validation of the taxonomy would be
    taking it for something it cannot be."""
    md = MD.read_text()
    assert "agreement on ADMISSION, not on content" in md
    for overclaim in ("validates the taxonomy", "the census is correct",
                      "confirms the descriptors"):
        assert overclaim not in md.lower()


def test_the_worked_explosion_example_is_present():
    """The concrete number is what stops someone 'simplifying' the query.

    Without it, `[mh:noexp]` reads as a stylistic choice rather than as the
    difference between 4,371 and 2,513.
    """
    md = MD.read_text()
    assert "4,371" in md and "[mh:noexp]" in SCRIPT.read_text()
    assert re.search(r"EXPLODES", md)
