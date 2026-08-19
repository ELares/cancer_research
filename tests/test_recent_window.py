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
import re
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


def test_the_cohort_sentence_counts_the_cohorts_it_says_it_counts():
    """"29 of the cohorts shown resolve at under 1%" -- 12 are shown.

    The count came from the FULL row list while the table printed `rows[:12]`,
    so the sentence was arithmetically impossible against the table beneath it.
    """
    idx, md = _doc()["indexing"], MD.read_text()
    allr = idx.get("rows_all")
    assert allr, (
        "only the printed cohorts are committed, so no count taken over every "
        "cohort can be checked -- which is how the impossible sentence shipped")
    assert idx["n_cohorts"] == len(allr)
    assert idx["rows"] == allr[:len(idx["rows"])], (
        "the printed rows are not the head of the full list")
    assert idx["n_cohorts_shown"] == len(idx["rows"])
    assert idx["n_shown_under_1pct"] == sum(
        1 for r in idx["rows"] if r["rate_pct"] < 1.0), (
        "the shown-cohort count is not computed over the shown cohorts")
    assert idx["n_older_cohorts"] == sum(
        1 for r in allr if r["year"] <= idx["older_cutoff_year"])
    assert idx["settled_years"] == sum(
        1 for r in allr
        if r["year"] <= idx["older_cutoff_year"] and r["rate_pct"] < 1.0)
    assert (f"Of the {idx['n_cohorts_shown']} shown, "
            f"{idx['n_shown_under_1pct']} resolve at under 1%") in md
    # any count quoted over ALL cohorts must name its own denominator
    assert (f"across all {idx['n_cohorts']} cohorts") in md, (
        "a count over every cohort is quoted without saying so, which is how "
        "29 came to sit beside a 12-row table")


def test_no_maximum_is_taken_over_a_set_selected_by_that_maximum():
    """`settled_max_rate_pct` was max over rows already filtered to <1%.

    It could not exceed 1% whatever the data did, and it was quoted as the
    evidence that older cohorts settle.
    """
    d, idx, md = _doc(), _doc()["indexing"], MD.read_text()
    assert "settled_max_rate_pct" not in json.dumps(idx), (
        "a maximum over a set selected for being below the threshold it is "
        "quoted against is back")
    src = SCRIPT.read_text()
    for m in re.finditer("settled_max_rate_pct", src):
        w = src[max(0, m.start() - 400):m.end() + 200]
        assert "USED TO BE" in w, (
            "the self-selecting maximum is computed again rather than only "
            "described in the note explaining why it was removed")
    allr = idx["rows_all"]
    older = [r["rate_pct"] for r in allr if r["year"] <= idx["older_cutoff_year"]]
    assert older, "no cohort old enough for the claim to be about"
    # RECOMPUTED over every cohort, including the ones the table truncates
    assert abs(idx["older_max_rate_pct"] - max(older)) < 1e-9, (
        f"the reported maximum is {idx['older_max_rate_pct']}% and the "
        f"cohorts give {max(older)}%")
    assert idx["older_max_rate_year"] == max(
        (r for r in allr if r["year"] <= idx["older_cutoff_year"]),
        key=lambda r: r["rate_pct"])["year"]
    assert f"**{idx['older_max_rate_pct']}%**" in md


def test_the_resolution_claim_is_not_monotone_in_age_if_the_data_is_not():
    """"a cohort not indexed by then generally never will be" is contradicted
    by the table printed directly beneath it, and by a cohort the 12-row
    truncation hid entirely.
    """
    idx, md = _doc()["indexing"], MD.read_text()
    faster = idx.get("older_than_and_faster_than_reference")
    assert faster is not None, "the monotonicity check was removed"
    ref = idx.get("reference_rate_pct")
    # RECOMPUTED from every cohort, so narrowing the window cannot hide one
    allr = idx["rows_all"]
    want = sorted(({"year": r["year"], "rate_pct": r["rate_pct"]}
                   for r in allr
                   if r["year"] < idx["reference_year"] and r["rate_pct"] > ref),
                  key=lambda r: -r["rate_pct"])
    assert faster == want, (
        f"the report lists {[r['year'] for r in faster]} as older-and-faster "
        f"and the cohorts give {[r['year'] for r in want]}")
    for r in faster:
        assert r["rate_pct"] > ref, "a listed cohort is not actually faster"
        assert r["year"] < idx["reference_year"], "a listed cohort is not older"
    if faster:
        assert "not monotone in age" in md, (
            f"{len(faster)} cohorts older than {idx['reference_year']} resolve "
            f"faster than it does, and the report does not say so")
        assert "is WITHDRAWN" in md
        for r in faster[:4]:
            assert f"{r['year']} at {r['rate_pct']}%" in md
    for gone in ("generally never will be, so the recent blind spot",
                 "mostly permanent, not a backlog"):
        assert gone not in md, f"the withdrawn claim {gone!r} is back"
    assert "bounds a per-window rate and not a lifetime" in md, (
        "the measurement covers one update window and the page does not say "
        "what that bounds")


def test_the_exit_rule_asymmetry_is_priced():
    """Enter on a 99% interval above 1.3, leave on a point at 1.0.

    That is an interval test on the way in and a point test on the way out,
    and the headline is what the mismatch produces.
    """
    d, md = _doc(), MD.read_text()
    comp = d["composition"]
    best = comp.get("by_comparator")
    assert best, "the comparator block is gone"
    for y, r in best.items():
        ev = r.get("exit_rule_variants") or {}
        assert len(ev) >= 3, f"{y}: the exit rule is reported without variants"
        shipped = [v for k, v in ev.items() if "shipped" in k]
        assert shipped and shipped[0] == r["flat_or_down"], (
            f"{y}: the shipped variant does not equal the headline count")
        mirror = [v for k, v in ev.items() if "mirrors admission" in k]
        assert mirror, f"{y}: the mirror of the admission rule is not computed"
        # A STRICTER TEST CANNOT ADMIT MORE. Ordered by strictness, the four
        # variants must be monotone -- an interval excluding 1.3 is the
        # hardest bar to clear, a point estimate at 1.3 the easiest.
        order = ["mirrors admission", "interval_excludes_1.0",
                 "point_le_1.0 (shipped)", "point_le_"]
        got = []
        for tag in order:
            hit = [v for k, v in ev.items()
                   if (tag in k) and not (tag == "point_le_"
                                          and "1.0" in k)]
            assert hit, f"{y}: no variant matching {tag!r}"
            got.append(hit[0])
        assert got == sorted(got), (
            f"{y}: the exit variants are not monotone in strictness {got}, so "
            "at least one is not the test its label claims")
        assert r.get("n_within_10pct_of_cut") is not None
        assert len(r.get("point_cut_sweep") or {}) >= 3
    # and the page must render them for the comparator it leads with
    # the comparator the page LEADS with is the second-newest, matching the
    # renderer's own choice -- quoting variants for a different one would
    # compare the table to a paragraph about another year
    years = sorted((int(y) for y in best), reverse=True)
    lead = str(years[1] if len(years) > 1 else years[0])
    r = best[lead]
    for k, v in (r.get("exit_rule_variants") or {}).items():
        assert f"| {k} | {v:,} |" in md, f"exit variant {k!r} is not rendered"
    assert "HOLDS ONE SIDE TO AN INTERVAL AND THE OTHER TO A POINT" in md
    assert "at the same z but a LOWER floor" in md, (
        "the report still says 'at the same confidence', which is the claim "
        "the asymmetry falsifies")


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
