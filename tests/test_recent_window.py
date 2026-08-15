"""Guards for the recent-window analysis.

WHAT THIS ANALYSIS CLAIMS, AND WHY IT IS EASY TO GET WRONG
----------------------------------------------------------
It says that about a third of the MeSH descriptors which rise significantly
against the whole census are no longer rising against the most recent complete
year alone -- i.e. that the obvious comparison measures the era rather than the
topic. That is a claim about a METHOD, so the ways it can be wrong are all
methodological:

1. THE COMPARATOR. Comparing against the newest year inflates the share,
   because the newest year is only partly indexed. The report must use the most
   recent COMPLETE year, the same exclusion the thesis-leg section applies.

2. AN EMPTY RESULT MUST NOT RENDER. The first version of the generator read a
   field named `mesh_terms`, which does not exist in these records -- the field
   is `mesh` -- and produced a report reading "Of **0 MeSH descriptors** ...
   **0 (0.0%)** are flat or lower". A missing field rendered as a measurement.
   This repository has hit that exact shape before, where a search of a field
   the frozen index does not carry returned 5 mentions instead of 38.

3. THE VERB. "Flat or falling" is a point estimate over a pool selected for
   having risen, so regression to the mean is expected. Only a fraction are
   demonstrably falling at the same confidence the pool was selected at, and
   the report must not claim decline.

4. THE SHARE IS NOT A CONSTANT. It ranges from single digits against a
   four-year-old comparator to over 40% against the newest, so the report has
   to show that dependence rather than quote one number.
"""

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_recent_window.py"
MD = REPO_ROOT / "analysis" / "atlas-recent-window.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-recent-window.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def test_the_headline_share_is_the_one_in_the_artifact():
    """Prose and JSON must agree, and the prose must not carry its own number."""
    d = _doc()
    comp = d["composition"]
    years = sorted(comp["by_comparator"], key=int, reverse=True)
    complete = years[1] if len(years) > 1 else years[0]
    row = comp["by_comparator"][complete]
    md = MD.read_text()
    assert f"{row['flat_or_down']:,} ({row['flat_share']}%)" in md, (
        f"the report does not state the {complete} share "
        f"{row['flat_or_down']:,} ({row['flat_share']}%) that its own JSON holds")
    assert f"comparator is {complete} alone" in md, (
        f"the report's headline comparator is not {complete}")


def test_the_comparator_is_a_complete_year_not_the_newest():
    """The newest year is partly indexed, which inflates the share.

    This is the same incomplete-year effect the thesis-leg section strips out,
    and using the newest year here would contradict it three screens apart.
    """
    d = _doc()
    comp = d["composition"]
    years = sorted(comp["by_comparator"], key=int, reverse=True)
    assert len(years) > 1, "only one comparator year; the dependence is untestable"
    newest, complete = years[0], years[1]
    md = MD.read_text()
    assert f"comparator is {newest} alone" not in md, (
        f"the headline compares against {newest}, which is the incomplete "
        "trailing year the report itself excludes elsewhere")
    # and the trailing year must actually be the inflating one, or the
    # exclusion is cargo-culted rather than doing work
    assert comp["by_comparator"][newest]["flat_share"] > \
        comp["by_comparator"][complete]["flat_share"], (
        "the trailing year no longer inflates the share, so the reason given "
        "for excluding it is not the reason it is excluded")
    assert d["legs"]["complete_through"] == int(complete), (
        f"the legs section treats {d['legs']['complete_through']} as the last "
        f"complete year while the composition section uses {complete}")


def test_an_empty_pool_refuses_to_render():
    """The missing-field failure, pinned. A wrong field name must raise."""
    src = SCRIPT.read_text()
    assert 'r.get("mesh")' in src, (
        "the generator no longer reads the `mesh` field; `mesh_terms` does not "
        "exist in these records and silently yields an empty pool")
    assert 'r.get("mesh_terms")' not in src, "the nonexistent field is back"
    assert "raise SystemExit" in src and "is not a finding" in src, (
        "an empty descriptor pool no longer refuses to render, so a wrong "
        "field name would print '0 descriptors (0.0%)' as a measurement")


def test_the_report_does_not_claim_decline():
    """Selection is on having risen, so regression to the mean is expected."""
    d = _doc()
    comp = d["composition"]
    for year, row in comp["by_comparator"].items():
        assert row["demonstrably_falling"] <= row["flat_or_down"], (
            f"{year}: more descriptors demonstrably falling "
            f"({row['demonstrably_falling']}) than flat-or-down "
            f"({row['flat_or_down']}), which is arithmetically impossible")
        assert row["flat_or_down"] <= comp["pool_size"], (
            f"{year}: flat-or-down exceeds the pool it was drawn from")
    md = MD.read_text()
    assert "no longer demonstrably *rising*" in md, (
        "the report no longer distinguishes 'stopped rising' from 'declining'")


def test_the_share_is_shown_to_depend_on_the_comparator():
    """One number here would be a claim the data does not support."""
    d = _doc()
    shares = {y: r["flat_share"] for y, r in d["composition"]["by_comparator"].items()}
    assert len(shares) >= 3, f"only {len(shares)} comparators; the dependence is not shown"
    lo, hi = min(shares.values()), max(shares.values())
    assert hi > lo * 1.5, (
        f"the share barely moves across comparators ({lo}% to {hi}%), so the "
        "claim that it depends on the comparator is not supported by this table")
    md = MD.read_text()
    for y in shares:
        assert f"| {y} |" in md, f"comparator {y} is in the JSON but not the table"


def test_the_unindexed_pool_collapse_is_real():
    """The claim is that resolution collapses with age, not that it is low."""
    idx = _doc()["indexing"]
    rows = {r["year"]: r["rate_pct"] for r in idx["rows"]}
    assert len(rows) >= 4, "too few cohorts to show a collapse"
    newest = max(rows)
    old = [r for y, r in rows.items() if y <= newest - 3]
    assert old, "no cohort three or more years old; the collapse is untestable"
    assert rows[newest] > 10 * max(old), (
        f"the newest cohort resolves at {rows[newest]}% against a maximum of "
        f"{max(old)}% among cohorts 3+ years old; that is not a collapse")
    assert idx["settled_years"] >= 2, (
        "fewer than two cohorts resolve below 1%, so 'permanently un-indexed' "
        "rests on a single year")


def test_the_thesis_legs_are_reported_on_complete_years():
    """The raw gain is a partial year refilling, not the field moving."""
    lg = _doc()["legs"]
    assert lg["gain_trailing_share_pct"] > 50, (
        "the trailing year no longer dominates the gain, so the complete-year "
        "restriction is not doing the work the report says it does")
    assert lg["ferroptosis_complete_after"] >= lg["ferroptosis_complete_before"]
    md = MD.read_text()
    assert f"through {lg['complete_through']}" in md, (
        "the report does not say which years the leg comparison is restricted to")


def test_the_generator_can_rebuild_the_report_from_the_artifact():
    """--render-only must not need the 2GB census to redraw the prose."""
    res = subprocess.run(
        [sys.executable, str(SCRIPT), "--render-only"],
        cwd=REPO_ROOT, capture_output=True, text=True)
    assert res.returncode == 0, (
        f"--render-only failed, so the report cannot be checked without the "
        f"bulk census:\n{res.stdout}\n{res.stderr}")
