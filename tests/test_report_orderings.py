"""The specific orderings this campaign repaired, pinned individually.

`tests/test_artifact_freshness.py` gates the CLASS: it shuffles an artifact's
dicts and requires the rendered text not to move, which catches an order
INHERITED from serialisation. What it cannot do is tell a correct rank from a
wrong-but-explicit one -- a renderer sorting the wrong way is deterministic,
reproducible, and green.

That distinction is not theoretical. Re-introducing this campaign's own
regression (an ascending sort where the sweep steps OUTWARD from a bound)
leaves the class gate passing with the bound on the last row, directly under
prose describing it as the first. And the regression that started all of this
-- a verdict flipping from "3 of 3" to "2 of 3" -- would likewise pass, because
the wrong order was perfectly deterministic.

So each repaired ordering is pinned here against the thing that makes it
correct: the prose beside it, or the header above it, or the constant that
defines the sequence. These assertions are about MEANING, which is why they
are hand-written and why a class gate cannot replace them.

OFFLINE: reads only committed artifacts.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
A = REPO / "analysis"


def _rows(md: str, pattern: str) -> list:
    return [l for l in md.splitlines() if re.match(pattern, l)]


def _nums(line: str) -> list:
    return [int(x.replace(",", "")) for x in re.findall(r"\b\d[\d,]*\b", line)]


def test_the_penetration_verdict_is_computed_over_the_declared_tissue_order():
    """The regression that started this: a MONOTONE test read dict order.

    It was correct only while the scan inserted well -> poorly -> cns; sorting
    alphabetically reversed the sequence and reported the gradient BROKEN for
    `default`, the one admissible set, whose kills are 12.10 > 2.60 > 1.80.
    """
    md = " ".join((A / "headline-at-fitted-cascade.md").read_text().split())
    assert "ordering is preserved in 3 of 3" in md, (
        "the penetration-ordering verdict is not 3 of 3; if the data changed "
        "say so, but check first that the tissue sequence is still in "
        "penetration order rather than alphabetical")
    body = (A / "headline-at-fitted-cascade.md").read_text()
    order = [l.split("|")[1].strip() for l in _rows(body, r"\| (well|poorly|cns)")]
    assert order == ["well_vascularized", "poorly_vascularized", "cns_bbb"], order
    # And the order must come from the constant that builds the dict, not a
    # copy typed beside it.
    src = (REPO / "scripts/headline_at_fitted.py").read_text()
    assert "PENETRATION_ORDER = tuple(k for k, _ in hs.PENETRATION_TISSUES)" in src, (
        "PENETRATION_ORDER is no longer derived from the constant that "
        "determines the keys, so the two can drift apart")


def test_the_bound_leads_each_abc_sweep_as_its_prose_says():
    """Each sweep steps OUTWARD from its bound (3->2->1, 6->8->10).

    An ascending sort restores `hill` and silently reverses `k_erastin`,
    putting "(prior low bound)" on the last row under a sentence reading
    "pushing `k_erastin` BELOW its bound".
    """
    md = (A / "calibration/abc-acceptance-diagnostic.md").read_text()
    ke = [l for l in md.splitlines() if l.startswith("| `k_erastin`")]
    hl = [l for l in md.splitlines() if l.startswith("| `hill`")]
    assert ke and hl
    assert "prior low bound" in ke[0], f"k_erastin must lead with its bound: {ke[0]}"
    assert "prior high bound" in hl[0], f"hill must lead with its bound: {hl[0]}"
    vals = [float(l.split("|")[2].strip()) for l in ke]
    assert vals == sorted(vals, reverse=True), (
        f"k_erastin must step DOWN away from its bound, got {vals}")


def test_the_retraction_family_table_is_count_ranked():
    """`Expression of Concern` (172) above `Retracted Publication` (7,900) is
    what rendering a round-tripped Counter produced."""
    md = (A / "atlas-retraction-exposure.md").read_text()
    rows = _rows(md, r"\| `[A-Z]")
    counts = [_nums(r)[-1] for r in rows if _nums(r)]
    assert counts == sorted(counts, reverse=True), (
        f"retraction-family table is not count-descending: {counts}")


def test_the_variant_spelling_table_leads_with_its_top_row():
    """Headed "top 6 shown", and the sentence below names the 712-row leader."""
    md = (A / "atlas-variant-drug-map.md").read_text()
    rows = _rows(md, r"\| `[cp]\.")
    counts = [_nums(r)[0] for r in rows if _nums(r)]
    assert counts, "no spelling rows found"
    assert counts == sorted(counts, reverse=True), (
        f"the top-N spelling table is not a rank: {counts}")
    assert counts[0] == max(counts)


def test_the_variant_predicate_columns_are_count_ranked():
    md = (A / "atlas-variant-drug-map.md").read_text()
    bad = []
    for line in md.splitlines():
        m = re.findall(r"`([a-z_]+)` (\d+)", line)
        if len(m) >= 2 and "|" in line:
            ns = [int(n) for _, n in m]
            if ns != sorted(ns, reverse=True):
                bad.append(line[:90])
    assert not bad, f"predicate columns not count-ranked: {bad[:3]}"


def test_the_untagged_partner_page_follows_its_own_opening_sentence():
    """It names radiotherapy, chemotherapy and surgery in that order."""
    md = (A / "atlas-untagged-partner.md").read_text()
    rows = [l.split("|")[1].strip() for l in _rows(md, r"\| (radiotherapy|chemotherapy|surgery) ")]
    assert rows[:3] == ["radiotherapy", "chemotherapy", "surgery"], (
        f"the tables contradict the page's own first paragraph: {rows[:3]}")
    src = (REPO / "scripts/atlas_untagged_partner.py").read_text()
    assert "_in_candidate_order" in src, (
        "the declared-order helper is gone, so the page will alphabetise again")


@pytest.mark.parametrize("artifact,pattern", [
    ("census-normal-tissue.md", r"\| (acute kidney|chemical and drug|cardiotox)"),
])
def test_ranked_tables_stay_ranked(artifact, pattern):
    md = (A / artifact).read_text()
    counts = [_nums(r)[-1] for r in _rows(md, pattern) if _nums(r)]
    assert counts == sorted(counts, reverse=True), f"{artifact}: {counts}"
