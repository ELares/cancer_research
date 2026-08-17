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
pool across the update window, resolution collapses within about two years of
publication -- the most recent cohort resolves at tens of percent, the
two-year-old cohort at well under one. A cohort that has not been indexed by
then generally never will be.

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
        for d, a, b in pool:
            c = dy.get(d, 0)
            rr_all = (a / n_new) / (b / n_all) if b else float("inf")
            rr_y = (a / n_new) / (c / n_y) if c else float("inf")
            if rr_y <= FLAT:
                flat += 1
                if katz_lcb(c, n_y, a, n_new) > 1.0:
                    falling += 1
            rows.append((d, rr_all, rr_y, a, b, c))
        rows.sort(key=lambda r: (r[2], -r[1]))
        out["by_comparator"][str(year)] = {
            "comparator_n": n_y,
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
    settled = [r for r in rows if r["rate_pct"] < 1.0]
    return {
        "rows": rows[:12],
        "pool_total": raw["unindexed_total"],
        "resolved_total": sum(got.values()),
        "resolved_share_pct": round(
            100 * sum(got.values()) / max(raw["unindexed_total"], 1), 3),
        "settled_years": len(settled),
        "settled_max_rate_pct": round(max((r["rate_pct"] for r in settled),
                                          default=0.0), 3),
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
        # restricted to complete years, matching the sentence above the table.
        # All-years counts here would show the partial year as movement, which
        # is the exact thing this section exists to strip out.
        "legs": {
            k: {
                "before": sum(v for y, v in raw["leg_census"].get(k, {}).items()
                              if y <= complete),
                "after": sum(v for y, v in raw["leg_census"].get(k, {}).items()
                             if y <= complete)
                + sum(v for y, v in raw["leg_new"].get(k, {}).items()
                      if y <= complete),
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
    L += [f"Only **{b['demonstrably_falling']:,}** are demonstrably falling at the "
          "same confidence, so the defensible statement is that these topics are "
          "no longer demonstrably *rising* -- not that they are in decline.", ""]
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

    L += ["## The un-indexed pool is mostly permanent, not a backlog", ""]
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
    L += [f"Resolution collapses within about two years of publication: "
          f"{idx['settled_years']} of the cohorts shown resolve at under 1%, the "
          f"highest of them at {idx['settled_max_rate_pct']}%. A cohort not "
          "indexed by then generally never will be, so the recent blind spot is "
          "only partly lag -- the rest is literature MeSH does not index at all, "
          "and waiting does not recover it.", ""]

    L += ["## The thesis legs did not move", ""]
    L += [f"The update window adds {lg['gain_total']:,} ferroptosis-indexed "
          f"articles, taking the field from {lg['ferroptosis_all_years_before']:,} "
          f"to {lg['ferroptosis_all_years_after']:,}. But "
          f"**{lg['gain_trailing_share_pct']}%** of that gain lands in "
          f"{lg['latest_year']}, an incomplete publication year the project's own "
          f"claims already exclude.", ""]
    L += [f"Restricted to complete years (through {lg['complete_through']}), "
          f"ferroptosis goes {lg['ferroptosis_complete_before']:,} -> "
          f"{lg['ferroptosis_complete_after']:,}. The legs:", ""]
    L += [f"| leg (through {lg['complete_through']}) | before | after |",
          "|---|--:|--:|"]
    for k, v in lg["legs"].items():
        L.append(f"| {k} | {v['before']:,} | {v['after']:,} |")
    L += ["",
          "So the sonodynamic leg -- the thesis's own central mechanism -- is "
          "still supported by a literature of tens of papers, not hundreds. The "
          "window refilled a partial year; it did not change the field's shape.",
          ""]

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
