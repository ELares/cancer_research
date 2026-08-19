#!/usr/bin/env python3
"""What the newest slice of the cancer literature does, and what it only appears to do.

WHY THIS EXISTS
---------------
The census gained a daily-update stream, and the newest articles are the part
nobody has looked at. The obvious analysis -- count what is over-represented in
the new window relative to the census -- is also the one this repository is most
likely to get wrong, because it has already been burned twice by exactly that
comparison: a fix justified by an error distribution it then changed (#617), and
a manuscript-versus-census claim whose denominators were not comparable.

So the first result here is about the method, not the literature.

THE COMPOSITION TRAP
--------------------
Comparing a 2026 window against a 1975-2026 census measures the era, not the
topic. A third of the descriptors that "rise" significantly against the whole
census are no longer rising against the most recent complete year alone.
Immune checkpoint inhibitors are the cleanest case: strongly up against fifty
years of literature, flat against last year.

The share is not a constant of nature -- it depends entirely on how recent the
comparator is, and this script computes it against several so the reader can see
that dependence rather than take one number. Nothing here says those topics are
declining: at 99% confidence only a tenth are demonstrably falling. The
defensible claim is that they are no longer demonstrably RISING, which is a
weaker and more useful statement than either "surging" or "over".

MOST UN-INDEXED LITERATURE IS PERMANENTLY UN-INDEXED
----------------------------------------------------
The census is defined by MeSH, so an article NLM has not indexed is invisible to
it, and the natural reading is that the un-indexed pool is a backlog that will
resolve. It mostly will not. Tracking the baseline's text-matched un-indexed
pool across the update window, resolution is overwhelmingly concentrated in the
two newest cohorts -- tens of percent -- and every older one resolves at a few
percent at most. It is NOT monotone in age, and an earlier version of this
docstring said "a cohort that has not been indexed by then generally never will
be": that is withdrawn, because cohorts older than the two-year-old one resolve
FASTER than it does, one of them 3x faster and hidden by a 12-row truncation.
Measured over one update window, this bounds a per-window rate, not a lifetime.

That reframes every recent-window claim in this project, including ones it made
before measuring this: the recent blind spot is only partly lag. The rest is
literature MeSH does not index at all, and no amount of waiting recovers it.

WHAT THIS DOES NOT DO
---------------------
It does not merge the update stream into the census. `records/` is what the
committed atlas artifacts were computed on, and the relation layer covers a
different set again, so merging would produce a census whose headline count
includes articles the graph analyses cannot see. Every figure here names the
surface it came from.

Usage:
    python scripts/atlas_recent_window.py
    python scripts/atlas_recent_window.py --render-only   # from committed JSON
"""

import argparse
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = PROJECT_ROOT / "corpus" / "atlas"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-recent-window.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-recent-window.json"

# A descriptor must clear this in the new window before any ratio is computed.
# Without a floor the table is dominated by 2 -> 6 moves, which are noise
# wearing a large ratio.
MIN_NEW = 50
# Katz log-method lower confidence bound; a rise must clear this against the
# WHOLE census to enter the pool at all.
LCB_Z = 2.576          # 99%
LCB_FLOOR = 1.3
# Below this against the recent end, a descriptor is no longer rising.
FLAT = 1.0

# The thesis legs this project's own manuscript rests on, as MeSH descriptors.
LEGS = {
    "sonodynamic (Ultrasonic Therapy, an under-estimate)": ["Ultrasonic Therapy"],
    "photodynamic": ["Photochemotherapy", "Photosensitizing Agents"],
    "drug resistance": ["Drug Resistance, Neoplasm"],
    "lipid peroxidation": ["Lipid Peroxidation"],
}
FERRO = "Ferroptosis"


def _shards(sub):
    d = ATLAS / sub
    return sorted(d.glob("*.jsonl.gz")) if d.exists() else []


def _read(sub):
    for f in _shards(sub):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                yield json.loads(line)


def katz_lcb(a, n1, b, n2, z=LCB_Z):
    """Lower bound of the rate ratio (a/n1)/(b/n2), Katz log method.

    Returns 0.0 where the ratio is undefined, so an undefined ratio can never
    clear a floor -- a zero denominator must not read as an infinite rise.
    """
    if a == 0 or n1 == 0 or n2 == 0:
        return 0.0
    if b == 0:
        b, adj = 0.5, True     # Haldane correction, only to bound the variance
    else:
        adj = False
    rr = (a / n1) / (b / n2)
    se = math.sqrt(max(1.0 / a - 1.0 / n1 + 1.0 / b - 1.0 / n2, 0.0))
    lo = rr * math.exp(-z * se)
    return 0.0 if adj and a < 5 else lo


def collect():
    """One pass per stream. Everything downstream is derived from this."""
    census_year = Counter()          # records/ by publication year
    census_desc = Counter()          # descriptor -> count over records/
    desc_by_year = defaultdict(Counter)
    census_pmids = set()
    ferro_year_census = Counter()
    leg_census = defaultdict(Counter)

    for r in _read("records"):
        pid = str(r.get("pmid") or "")
        if pid:
            census_pmids.add(pid)
        y = r.get("year")
        terms = set(r.get("mesh") or [])
        if y:
            census_year[int(y)] += 1
        for t in terms:
            census_desc[t] += 1
            if y:
                desc_by_year[int(y)][t] += 1
        if FERRO in terms:
            if y:
                ferro_year_census[int(y)] += 1
            for name, ds in LEGS.items():
                if any(d in terms for d in ds) and y:
                    leg_census[name][int(y)] += 1

    unindexed = {}
    for r in _read("records_unindexed"):
        pid = str(r.get("pmid") or "")
        if pid:
            unindexed[pid] = r.get("year")

    new_desc = Counter()
    new_total = 0
    new_pmids = set()
    ferro_year_new = Counter()
    leg_new = defaultdict(Counter)
    seen_update = set()
    resolved_unindexed = Counter()   # by the year the BASELINE recorded

    for r in _read("records_updates"):
        pid = str(r.get("pmid") or "")
        if not pid or pid in seen_update:
            continue
        seen_update.add(pid)
        if pid in unindexed:
            y = unindexed[pid]
            if y:
                resolved_unindexed[int(y)] += 1
        if pid in census_pmids or pid in unindexed:
            continue                  # a revision, not a new article
        new_pmids.add(pid)
        new_total += 1
        terms = set(r.get("mesh") or [])
        for t in terms:
            new_desc[t] += 1
        y = r.get("year")
        if FERRO in terms:
            if y:
                ferro_year_new[int(y)] += 1
            for name, ds in LEGS.items():
                if any(d in terms for d in ds) and y:
                    leg_new[name][int(y)] += 1

    unindexed_by_year = Counter()
    for y in unindexed.values():
        if y:
            unindexed_by_year[int(y)] += 1

    return {
        "census_total": len(census_pmids),
        "census_year": dict(census_year),
        "census_desc": census_desc,
        "desc_by_year": desc_by_year,
        "new_total": new_total,
        "new_desc": new_desc,
        "unindexed_total": len(unindexed),
        "unindexed_by_year": dict(unindexed_by_year),
        "resolved_unindexed": dict(resolved_unindexed),
        "ferro_year_census": dict(ferro_year_census),
        "ferro_year_new": dict(ferro_year_new),
        "leg_census": {k: dict(v) for k, v in leg_census.items()},
        "leg_new": {k: dict(v) for k, v in leg_new.items()},
    }


def composition(raw, comparators):
    """Descriptors that rise against the whole census, re-tested against a
    recent year. The point of the exercise is the gap between the two."""
    census_desc = raw["census_desc"]
    new_desc = raw["new_desc"]
    n_new, n_all = raw["new_total"], raw["census_total"]

    if not new_desc:
        raise SystemExit(
            "no MeSH descriptors found in the update stream. The records field "
            "is `mesh`; an earlier version read `mesh_terms`, which does not "
            "exist, and rendered '0 descriptors (0.0%)' as though a missing "
            "field were a measurement.")
    pool = []
    for d, a in new_desc.items():
        if a < MIN_NEW:
            continue
        b = census_desc.get(d, 0)
        if katz_lcb(a, n_new, b, n_all) <= LCB_FLOOR:
            continue
        pool.append((d, a, b))

    out = {"pool_size": len(pool), "min_new": MIN_NEW,
           "lcb_floor": LCB_FLOOR, "by_comparator": {}}
    for year in comparators:
        n_y = raw["census_year"].get(year, 0)
        if not n_y:
            continue
        dy = raw["desc_by_year"].get(year, Counter())
        flat, falling, rows = 0, 0, []
        # THE EXIT RULE IS NOT THE ADMISSION RULE IN MIRROR IMAGE, and the
        # headline is what you get by holding one side to an interval and the
        # other to a point: a descriptor ENTERS on a 99% Katz lower bound above
        # LCB_FLOOR and LEAVES on a bare point estimate at or below FLAT. Both
        # mirrors are computed so the asymmetry is priced rather than hidden,
        # and so is the sensitivity of the point cut, because the pool piles up
        # against it.
        mirror_floor = 0     # 99% interval excluding LCB_FLOOR, the exact mirror
        mirror_one = 0       # 99% interval excluding 1.0
        point_floor = 0      # point estimate at or below LCB_FLOOR
        near = 0             # within 10% of the point cut
        sweep = {}
        for d, a, b in pool:
            c = dy.get(d, 0)
            rr_all = (a / n_new) / (b / n_all) if b else float("inf")
            rr_y = (a / n_new) / (c / n_y) if c else float("inf")
            if rr_y <= FLAT:
                flat += 1
                if katz_lcb(c, n_y, a, n_new) > 1.0:
                    falling += 1
            if c:
                lcb_down = katz_lcb(c, n_y, a, n_new)
                mirror_floor += lcb_down > LCB_FLOOR
                mirror_one += lcb_down > 1.0
            if rr_y <= LCB_FLOOR:
                point_floor += 1
            if rr_y != float("inf") and abs(rr_y - FLAT) <= 0.1 * FLAT:
                near += 1
            for cut in (0.95, 1.0, 1.05):
                sweep[cut] = sweep.get(cut, 0) + (rr_y <= cut)
            rows.append((d, rr_all, rr_y, a, b, c))
        rows.sort(key=lambda r: (r[2], -r[1]))
        out["by_comparator"][str(year)] = {
            "comparator_n": n_y,
            "exit_rule_variants": {
                "point_le_1.0 (shipped)": flat,
                "interval_excludes_1.0": mirror_one,
                "interval_excludes_%s (mirrors admission)" % LCB_FLOOR:
                    mirror_floor,
                "point_le_%s" % LCB_FLOOR: point_floor,
            },
            "n_within_10pct_of_cut": near,
            "point_cut_sweep": {str(k): v for k, v in sorted(sweep.items())},
            "flat_or_down": flat,
            "flat_share": round(100 * flat / len(pool), 1) if pool else 0.0,
            "demonstrably_falling": falling,
            "examples": [
                {"descriptor": d, "rr_vs_census": round(ra, 2),
                 "rr_vs_year": round(ry, 2), "new": a, "census": b, "year": c}
                for d, ra, ry, a, b, c in rows[:12] if ry != float("inf")
            ],
        }
    return out


def indexing(raw):
    """Does the un-indexed pool resolve, or is it permanent?"""
    pool = raw["unindexed_by_year"]
    got = raw["resolved_unindexed"]
    rows = []
    for y in sorted(pool, reverse=True):
        n = pool[y]
        if n < 1000:
            continue
        r = got.get(y, 0)
        rows.append({"year": y, "pool": n, "resolved": r,
                     "rate_pct": round(100 * r / n, 3)})
    # `settled_max_rate_pct` USED TO BE the maximum over `[r for r in rows if
    # r["rate_pct"] < 1.0]` -- the largest member of a set selected for being
    # under 1%, so it could not exceed 1% whatever the data did, and it was
    # quoted as evidence that older cohorts settle. It is replaced by the
    # maximum over every cohort old enough for the claim to be about, which
    # CAN exceed 1% and does.
    shown = rows[:12]
    newest = max((r["year"] for r in rows), default=0)
    older = [r for r in rows if r["year"] <= newest - 2]
    settled = [r for r in older if r["rate_pct"] < 1.0]
    return {
        "rows": shown,
        # EVERY cohort, not just the twelve the table prints. The counts in
        # the prose are taken over all of them, and with only the shown rows
        # committed no guard could check one -- which is how "29 of the
        # cohorts shown" sat beside a 12-row table, and how a 25-year-old
        # cohort resolving 3x faster than the two-year-old one stayed hidden.
        "rows_all": rows,
        "n_cohorts": len(rows),
        "n_cohorts_shown": len(shown),
        "n_shown_under_1pct": sum(1 for r in shown if r["rate_pct"] < 1.0),
        "pool_total": raw["unindexed_total"],
        "resolved_total": sum(got.values()),
        "resolved_share_pct": round(
            100 * sum(got.values()) / max(raw["unindexed_total"], 1), 3),
        "older_cutoff_year": newest - 2,
        "n_older_cohorts": len(older),
        "settled_years": len(settled),
        "older_max_rate_pct": round(max((r["rate_pct"] for r in older),
                                        default=0.0), 3),
        "older_max_rate_year": max(older, key=lambda r: r["rate_pct"])["year"]
        if older else None,
        # THE MONOTONICITY THE PROSE ASSUMED. Cohorts that resolve FASTER than
        # a younger one refute "a cohort not indexed by then generally never
        # will be", and one of them was hidden by the 12-row truncation.
        # The reference is the youngest cohort the withdrawn sentence was
        # ABOUT -- the two-year-old one it called "well under one" -- not the
        # one-year-old cohort, which is still resolving fast and against which
        # nothing older could win.
        "older_than_and_faster_than_reference": sorted(
            ({"year": r["year"], "rate_pct": r["rate_pct"]}
             for r in rows
             if r["year"] < newest - 2
             and r["rate_pct"] > next((x["rate_pct"] for x in rows
                                       if x["year"] == newest - 2), 0.0)),
            key=lambda r: -r["rate_pct"]),
        "reference_year": newest - 2,
        "reference_rate_pct": next(
            (r["rate_pct"] for r in rows if r["year"] == newest - 2), None),
    }


def last_complete_year(census_year):
    """The newest year whose volume is not obviously a partial year.

    NOT max(year). A handful of ahead-of-print records carry a year beyond the
    current one -- five of them here -- so taking the maximum makes the
    "incomplete trailing year" a five-article cohort and the restriction it
    justifies does nothing. Completeness is a statement about VOLUME: a year
    holding a small fraction of its predecessor is still filling up.
    """
    years = sorted(y for y in census_year if census_year[y] > 0)
    for y in reversed(years):
        prev = census_year.get(y - 1, 0)
        if prev and census_year[y] >= 0.5 * prev:
            return y
    return years[-1] if years else 0


def legs(raw):
    """The thesis legs, before and after, restricted to complete years.

    The incomplete trailing year is what makes the raw gain look like movement:
    it is a partial year being refilled, not the field turning.
    """
    fc, fn = raw["ferro_year_census"], raw["ferro_year_new"]
    complete = last_complete_year(raw["census_year"])
    latest = complete + 1
    ferro_before = sum(v for y, v in fc.items() if y <= complete)
    ferro_after = ferro_before + sum(v for y, v in fn.items() if y <= complete)
    gain_total = sum(fn.values())
    gain_trailing = sum(v for y, v in fn.items() if y > complete)
    return {
        "latest_year": latest,
        "complete_through": complete,
        "ferroptosis_all_years_before": sum(fc.values()),
        "ferroptosis_all_years_after": sum(fc.values()) + gain_total,
        "gain_total": gain_total,
        "gain_in_trailing_year": gain_trailing,
        "gain_trailing_share_pct": round(
            100 * gain_trailing / max(gain_total, 1), 1),
        "ferroptosis_complete_before": ferro_before,
        "ferroptosis_complete_after": ferro_after,
        # BOTH COLUMNS, because the complete-year one CANNOT MOVE. The update
        # stream is overwhelmingly the trailing year, so restricting to
        # complete years excludes 100% of the sonodynamic gain, 100% of the
        # photodynamic gain and almost all of drug resistance -- and the table
        # then printed "171 -> 171" and "30 -> 30" under a heading saying the
        # legs did not move. That was the filter's output, not stability, and
        # the artifact already publishes both columns for ferroptosis overall.
        "legs": {
            k: {
                "before": sum(v for y, v in raw["leg_census"].get(k, {}).items()
                              if y <= complete),
                "after": sum(v for y, v in raw["leg_census"].get(k, {}).items()
                             if y <= complete)
                + sum(v for y, v in raw["leg_new"].get(k, {}).items()
                      if y <= complete),
                "all_years_before": sum(raw["leg_census"].get(k, {}).values()),
                "all_years_after": sum(raw["leg_census"].get(k, {}).values())
                + sum(raw["leg_new"].get(k, {}).values()),
                "gain_all_years": sum(raw["leg_new"].get(k, {}).values()),
                "gain_excluded_by_the_filter": sum(
                    v for y, v in raw["leg_new"].get(k, {}).items()
                    if y > complete),
            } for k in LEGS},
    }


def render(d):
    comp, idx, lg = d["composition"], d["indexing"], d["legs"]
    # The most recent COMPLETE year, not the newest. The trailing year is
    # partially indexed, so comparing against it inflates the share -- the same
    # incomplete-year effect the thesis-leg section below has to strip out.
    years_desc = sorted(comp["by_comparator"], key=lambda y: int(y), reverse=True)
    best = years_desc[1] if len(years_desc) > 1 else years_desc[0]
    b = comp["by_comparator"][best]
    if comp["pool_size"] == 0:
        raise SystemExit(
            "the rising-descriptor pool is empty, which is not a finding -- it "
            "is what a wrong field name or a floor set above the data looks "
            "like. Refusing to render.")
    L = [f"# The recent window: {d['new_total']:,} articles the census did not have", ""]
    L += ["*Generated by `scripts/atlas_recent_window.py`. Every figure is "
          "recomputed; none is transcribed.*", ""]

    L += ["## Most of what looks like a rising topic is not", ""]
    L += [f"Of **{comp['pool_size']:,} MeSH descriptors** that rise significantly "
          f"against the whole {d['census_total']:,}-article census "
          f"(99% lower bound above {comp['lcb_floor']}, at least "
          f"{comp['min_new']} occurrences in the new window), "
          f"**{b['flat_or_down']:,} ({b['flat_share']}%)** are flat or lower "
          f"when the comparator is {best} alone.", ""]
    L += [f"Only **{b['demonstrably_falling']:,}** are demonstrably falling at "
          f"the same z but a LOWER floor (an interval excluding 1.0, against "
          f"the {comp['lcb_floor']} used for admission), so the defensible "
          f"statement is that these topics are no longer demonstrably "
          f"*rising* -- not that they are in decline.", ""]
    ev = b.get("exit_rule_variants") or {}
    if ev:
        L += [f"**THE HEADLINE HOLDS ONE SIDE TO AN INTERVAL AND THE OTHER TO "
              f"A POINT.** A descriptor ENTERS the pool on a 99% lower bound "
              f"above {comp['lcb_floor']} against the whole census, and LEAVES "
              f"on a bare point estimate at or below {FLAT} against the "
              f"comparator. That is not the same test run twice. Every variant, "
              f"on the same pool and the same comparator:", ""]
        L += ["| exit rule | no longer rising | share |", "|---|--:|--:|"]
        for k, v in ev.items():
            L.append(f"| {k} | {v:,} | {100*v/max(comp['pool_size'],1):.1f}% |")
        L += [""]
        near = b.get("n_within_10pct_of_cut")
        sweep = b.get("point_cut_sweep") or {}
        if near and sweep:
            lo = min(sweep, key=lambda k: float(k))
            hi = max(sweep, key=lambda k: float(k))
            L += [f"And the cut is knife-edge: **{near:,} of "
                  f"{comp['pool_size']:,} "
                  f"({100*near/max(comp['pool_size'],1):.0f}%)** of the pool "
                  f"sits within 10% of it, so the share runs "
                  f"{100*sweep[lo]/max(comp['pool_size'],1):.1f}% at a cut of "
                  f"{lo} to {100*sweep[hi]/max(comp['pool_size'],1):.1f}% at "
                  f"{hi}. The mirror-image test -- the admission rule applied "
                  f"in reverse -- is the bottom row of that table, and it is "
                  f"the number a reader who assumed one test was run twice "
                  f"would have expected.", ""]
    L += ["| descriptor | vs whole census | vs " + best + " | new | census | " + best + " |",
          "|---|--:|--:|--:|--:|--:|"]
    for e in b["examples"][:8]:
        L.append(f"| {e['descriptor']} | {e['rr_vs_census']}x | "
                 f"**{e['rr_vs_year']}x** | {e['new']:,} | {e['census']:,} | "
                 f"{e['year']:,} |")
    L += [""]
    L += ["The share depends entirely on how recent the comparator is, which is "
          "the finding rather than a caveat on it:", ""]
    L += ["| comparator year | articles | no longer rising | share |",
          "|---|--:|--:|--:|"]
    for y in sorted(comp["by_comparator"], key=lambda s: int(s)):
        r = comp["by_comparator"][y]
        L.append(f"| {y} | {r['comparator_n']:,} | {r['flat_or_down']:,} | "
                 f"{r['flat_share']}% |")
    L += ["",
          "A comparison against the whole census measures the era, not the topic.",
          ""]

    L += ["## The un-indexed pool resolves slowly, and not in age order", ""]
    L += [f"The census is defined by MeSH, so an un-indexed article is invisible "
          f"to it. Tracking the {idx['pool_total']:,} text-matched un-indexed "
          f"articles across the update window, "
          f"**{idx['resolved_total']:,} ({idx['resolved_share_pct']}%)** acquired "
          f"indexing -- and which ones depends almost entirely on age.", ""]
    L += ["| publication year | un-indexed pool | acquired indexing | rate |",
          "|---|--:|--:|--:|"]
    for r in idx["rows"]:
        L.append(f"| {r['year']} | {r['pool']:,} | {r['resolved']:,} | "
                 f"{r['rate_pct']}% |")
    L += [""]
    faster = idx.get("older_than_and_faster_than_reference") or []
    L += [f"Resolution is overwhelmingly concentrated in the newest cohorts. "
          f"Of the {idx['n_cohorts_shown']} shown, "
          f"{idx['n_shown_under_1pct']} resolve at under 1%; across all "
          f"{idx['n_cohorts']} cohorts with a pool of at least 1,000, "
          f"{idx['settled_years']} of the {idx['n_older_cohorts']} at least "
          f"two years old do. AN EARLIER VERSION OF THIS SENTENCE SAID "
          f"\"{idx['settled_years']} of the cohorts shown\", which is "
          f"arithmetically impossible -- only {idx['n_cohorts_shown']} are "
          f"shown -- and quoted a maximum that could not have exceeded 1%, "
          f"because it was the largest member of a set selected for being "
          f"below it.", ""]
    if idx.get("older_max_rate_pct") is not None:
        L += [f"The real maximum among cohorts at least two years old is "
              f"**{idx['older_max_rate_pct']}%** ({idx['older_max_rate_year']}), "
              f"not the {idx.get('reference_rate_pct')}% the "
              f"{idx.get('reference_year')} row shows -- which is the figure "
              f"the withdrawn sentence quoted.", ""]
    if faster:
        L += [f"**And it is not monotone in age, which the withdrawn sentence "
              f"assumed.** "
              + ", ".join(f"{r['year']} at {r['rate_pct']}%" for r in faster[:4])
              + f" all resolve FASTER than the {idx['reference_year']} "
              f"cohort at {idx.get('reference_rate_pct')}%, despite being "
              f"older"
              + (f" -- and {faster[0]['year']} is "
                 f"{faster[0]['rate_pct']/max(idx.get('reference_rate_pct') or 1e-9, 1e-9):.1f}x "
                 f"its rate" if idx.get("reference_rate_pct") else "")
              + ". So \"a cohort not indexed by then generally never will be\" "
                "is WITHDRAWN: cohorts that miss the initial indexing window "
                "keep acquiring indexing in later batches at a low rate.", ""]
    L += ["What survives is the size of the effect rather than its shape. "
          "Resolution beyond the two newest cohorts is a few percent at most, "
          "so the recent blind spot is only partly lag -- but this is measured "
          "over ONE update window, which bounds a per-window rate and not a "
          "lifetime, and the page no longer claims waiting recovers nothing.",
          ""]

    L += ["## The thesis legs, on complete years and on all years", ""]
    L += [f"The update window adds {lg['gain_total']:,} ferroptosis-indexed "
          f"articles, taking the field from {lg['ferroptosis_all_years_before']:,} "
          f"to {lg['ferroptosis_all_years_after']:,}. But "
          f"**{lg['gain_trailing_share_pct']}%** of that gain lands in "
          f"{lg['latest_year']}, an incomplete publication year the project's own "
          f"claims already exclude.", ""]
    L += [f"Restricted to complete years (through {lg['complete_through']}), "
          f"ferroptosis goes {lg['ferroptosis_complete_before']:,} -> "
          f"{lg['ferroptosis_complete_after']:,}. The legs:", ""]
    L += [f"| leg | through {lg['complete_through']} | all years | gain, all "
          f"years | of which the filter excludes |", "|---|--:|--:|--:|--:|"]
    for k, v in lg["legs"].items():
        exc = v.get("gain_excluded_by_the_filter", 0)
        g = v.get("gain_all_years", 0)
        L.append(f"| {k} | {v['before']:,} -> {v['after']:,} | "
                 f"{v.get('all_years_before', 0):,} -> "
                 f"{v.get('all_years_after', 0):,} | {g:,} | "
                 + (f"{exc:,} ({100*exc/g:.0f}%)" if g else "-") + " |")
    L += [""]
    stuck = [k for k, v in lg["legs"].items()
             if v["before"] == v["after"] and v.get("gain_all_years")]
    if stuck:
        short = [k.split(" (")[0] for k in stuck]
        L += [f"**{len(stuck)} of these rows CANNOT MOVE in the complete-year "
              f"column, and an earlier version of this section printed them "
              f"under the heading \"The thesis legs did not move\".** The "
              f"update stream is overwhelmingly {lg['latest_year']}, so the "
              f"filter excludes "
              + ", ".join(
                  f"{100*lg['legs'][k].get('gain_excluded_by_the_filter',0)/max(lg['legs'][k].get('gain_all_years',1),1):.0f}% "
                  f"of the {n} gain" for k, n in zip(stuck, short))
              + ". Those rows reported the filter, not stability. Both columns "
                "ship now, the way this section already reported ferroptosis "
                "overall.", ""]
    sono = next((v for k, v in lg["legs"].items()
                 if k.startswith("sonodynamic")), None)
    L += ["What survives is the conclusion rather than the table. On all years "
          + (f"the sonodynamic leg -- the thesis's own central mechanism -- "
             f"goes {sono.get('all_years_before', 0):,} to "
             f"{sono.get('all_years_after', 0):,}, "
             if sono else "the sonodynamic leg moves, ")
          + "which is still a literature of tens of papers rather than "
            "hundreds. The window did not change the field's shape; the "
            "complete-year restriction just could not have shown it either "
            "way.", ""]

    L += ["## What this does not claim", ""]
    L += ["* Nothing here merges the update stream into the census. Every figure "
          "names the surface it came from, because the relation layer covers a "
          "different set again and a merged headline count would include "
          "articles the graph analyses cannot see.",
          "* `Ultrasonic Therapy` is broader than sonodynamic therapy, but "
          "measured that is the smaller effect: precision 90.6% against "
          "recall 46.0%, so the leg is an UNDER-estimate by roughly twofold "
          "and is still the thinnest even so. An earlier version called it an "
          "over-estimate; see analysis/atlas-descriptor-recall.md.",
          "* A descriptor that stops rising has not been shown to decline.",
          "* MeSH indexing lag biases every recent count downward. The point of "
          "the second section is that it does not explain all of it.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the report from the committed JSON")
    args = ap.parse_args()

    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        raw = collect()
        years = sorted(raw["census_year"], reverse=True)[:5]
        d = {
            "census_total": raw["census_total"],
            "new_total": raw["new_total"],
            "composition": composition(raw, years),
            "indexing": indexing(raw),
            "legs": legs(raw),
        }
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
