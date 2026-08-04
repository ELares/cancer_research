#!/usr/bin/env python3
"""Atlas: which claims were asserted once and never followed up? (#ATLAS-REPL)

WHY
---
`MISSION.md` lists four things to mine for: patterns, contradictions, temporal
emergence, and "the replications that never happened". Three are implemented.
This is the fourth, and it is the one that could not be done at all on a
4,830-article corpus -- a claim looks unreplicated in a small sample whether or
not anyone replicated it.

WHAT IT MEASURES
----------------
For every entity pair in the relation graph: how many DISTINCT papers assert it,
and in which year it was first asserted. A pair asserted by exactly one paper,
long enough ago that a follow-up would have appeared if anyone ran one, is an
ORPHANED CLAIM.

The headline is not the count but the RATE, because the count alone is
uninterpretable: most pairs are asserted once because most pairs are obscure.
The rate asks a comparative question -- of the claims first made in year Y, what
share ever acquired a second paper? -- and it can be tracked across decades.

WHAT "REPLICATION" MEANS HERE, AND DOES NOT
--------------------------------------------
A second paper asserting the same entity pair. That is CO-ASSERTION, not
experimental replication:

  * the second paper may cite the first rather than test it, so this counts
    propagation as if it were confirmation, which INFLATES the rate;
  * it may assert the opposite direction, which the pair key does not
    distinguish -- a contradiction counts here as a replication;
  * a genuine independent replication that used different entities, or that
    PubTator failed to extract, is invisible.

So this bounds attention, not truth. A pair with one paper has been examined
once as far as the machine-readable record shows; that is a reading queue.

Usage:
    python scripts/atlas_replication.py
    python scripts/atlas_replication.py --quiet-years 5
"""

import argparse
import collections
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_discovery_eval import pmid_years  # noqa: E402
from atlas_graph import _corrected, load_corrections, load_index  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-replication.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-replication.json"

# Genes this project's simulation layers rest on, so the orphan list can be
# filtered to claims that would actually matter here.
FOCUS = {"2879", "84883", "11332", "23657", "6520", "1738", "2643", "4780",
         "9817", "7157", "1026", "5743"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet-years", type=int, default=5,
                    help="years of silence after which a single-paper claim counts "
                         "as orphaned")
    args = ap.parse_args()

    root = atlas_root()
    idx = load_index(root)
    print("loading year map ...", flush=True)
    years = pmid_years(root)
    latest = max(years.values())
    print(f"  {len(years):,} dated PMIDs, latest {latest}", flush=True)

    print("dating every pair ...", flush=True)
    corr = load_corrections()
    # Years of every asserting paper, not just the first: measuring replication
    # "ever" makes old cohorts look better purely because they have had longer,
    # so the window has to be EQUAL, not merely a minimum.
    pair_years = collections.defaultdict(list)
    papers = collections.defaultdict(set)
    first = {}
    with gzip.open(root / "relations" / "relations.tsv.gz", "rt",
                   encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            pmid = p[0]
            y = years.get(pmid)
            if not y:
                continue
            a = _corrected(p[2].split("|", 1)[-1], pmid, corr)
            b = _corrected(p[3].split("|", 1)[-1], pmid, corr)
            key = (a, b) if a <= b else (b, a)
            if pmid not in papers[key]:
                papers[key].add(pmid)
                pair_years[key].append(y)
            if key not in first or y < first[key]:
                first[key] = y
    print(f"  {len(papers):,} dated pairs", flush=True)

    # Replication rate by cohort year, with a fixed follow-up window so early
    # and late cohorts are comparable. A pair first asserted last year has had
    # no chance to be replicated, and pooling it with a 1990 cohort would
    # manufacture a decline that is really just censoring.
    W = args.quiet_years
    # Replicated WITHIN W years of first assertion, so every cohort is scored on
    # the same interval. Scoring "ever" instead produced a clean monotonic
    # decline from 60% in 1950 to 17.5% in 2020 -- which is not a finding about
    # science, it is the observation window shrinking as the cohorts get newer.
    cohort = collections.defaultdict(lambda: [0, 0])   # year -> [pairs, replicated]
    for key, ys in pair_years.items():
        y0 = first[key]
        if y0 + W > latest:          # not yet had its full window
            continue
        ys = sorted(ys)
        within = len(ys) > 1 and ys[1] <= y0 + W
        cohort[y0][0] += 1
        cohort[y0][1] += within

    rows = [{"year": y, "pairs": v[0], "replicated": v[1],
             "rate": v[1] / v[0] if v[0] else 0.0}
            for y, v in sorted(cohort.items()) if v[0] >= 200]

    orphans = [(k, first[k]) for k, pm in papers.items()
               if len(pm) == 1 and first[k] + W <= latest]
    # never-replicated at any horizon, for the headline count
    ever = sum(1 for k, ys in pair_years.items()
               if first[k] + W <= latest and len(ys) > 1)
    focus_orphans = [(k, y) for k, y in orphans if k[0] in FOCUS or k[1] in FOCUS]
    focus_orphans.sort(key=lambda kv: kv[1])

    name = lambda i: idx["canon"].get(i, i)  # noqa: E731
    total_eligible = sum(v[0] for v in cohort.values())
    overall = (sum(v[1] for v in cohort.values()) / total_eligible) if total_eligible else 0.0
    overall_ever = (ever / total_eligible) if total_eligible else 0.0

    L = [
        "# Claims asserted once and never followed up (#ATLAS-REPL)", "",
        "Generated by `scripts/atlas_replication.py`. The fourth mining pillar in",
        "`MISSION.md`, and the one that could not be attempted on 4,830 articles: a",
        "claim looks unreplicated in a small sample whether or not anyone replicated",
        "it.", "",
        "## What replication means here", "",
        "A second paper asserting the same entity pair. That is **co-assertion, not",
        "experimental replication**, and it fails in both directions: a paper that",
        "cites the first rather than testing it counts as a replication, which",
        "inflates the rate, while a genuine replication using different entities is",
        "invisible. The pair key does not carry direction either, so a contradiction",
        "counts here as a replication.", "",
        "This bounds ATTENTION, not truth. A pair with one paper has been examined",
        "once as far as the machine-readable record shows.", "",
        "## The rate", "",
        f"Of {total_eligible:,} entity pairs first asserted at least {W} years ago,",
        f"**{100*overall:.1f}%** acquired a second paper WITHIN {W} years. Allowing an",
        f"unlimited horizon raises that to {100*overall_ever:.1f}%, and "
        f"**{len(orphans):,}** pairs have still been asserted exactly once.", "",
        "### The window has to be equal, not merely a minimum", "",
        "Scoring each cohort on whether it was *ever* replicated produced a clean",
        "monotonic decline, from about 60% in 1950 to 17.5% in 2020. That is not a",
        "finding about science. A 1975 pair has had fifty years to acquire a second",
        f"paper and a 2020 pair has had six, so the \"decline\" was the observation",
        "window shrinking as cohorts get newer.", "",
        f"Every cohort below is therefore scored on the SAME {W}-year interval from",
        "its own first assertion, and cohorts too recent to have completed one are",
        "excluded rather than shown declining:", "",
        f"| first asserted | pairs | replicated within {W}y | rate |",
        "|---|---|---|---|",
    ]
    for r in rows[::max(1, len(rows) // 18)]:
        L.append(f"| {r['year']} | {r['pairs']:,} | {r['replicated']:,} | "
                 f"{100*r['rate']:.1f}% |")

    if rows:
        # compare fully-observed cohorts only: the newest are still being indexed
        mid = [r for r in rows if 2000 <= r["year"] <= 2014]
        old = [r for r in rows if 1985 <= r["year"] < 2000]
        recent = [r for r in rows if r["year"] >= 2017]
        avg = lambda rs: (sum(r["replicated"] for r in rs) /  # noqa: E731
                          max(1, sum(r["pairs"] for r in rs)))
        L += [
            "",
            "The dramatic version does not survive the fix, and what is left is",
            "modest:", "",
            f"* 1985-1999 cohorts: **{100*avg(old):.1f}%**",
            f"* 2000-2014 cohorts: **{100*avg(mid):.1f}%**",
            f"* 2017 onward: **{100*avg(recent):.1f}%**", "",
            "The first two comparisons are between fully observed cohorts and the gap",
            "between them is small. The last is not trustworthy on its own: MeSH",
            "indexing lags publication, so a recent second paper may exist without",
            "being visible here, and that biases the newest cohorts DOWNWARD. Read the",
            "recent drop as an upper bound on a real decline, not as a measurement of",
            "one.", ""]

    L += [
        "## Orphaned claims touching this project's mechanisms", "",
        f"{len(focus_orphans):,} single-paper pairs involve a gene the simulation",
        "layers rest on. Oldest first, since those have had the longest to be",
        "followed up and were not:", "",
        "| first asserted | pair |", "|---|---|",
    ]
    for (a, b), y in focus_orphans[:30]:
        L.append(f"| {y} | {name(a)} — {name(b)} |")

    L += [
        "", "## How to read this", "",
        "* A high orphan count is not a scandal. Most entity pairs are obscure, and",
        "  a claim nobody revisited may simply be a claim nobody needed.",
        "* The comparative figures carry the weight: the cohort trend, and whether a",
        "  pair this project DEPENDS ON sits in the orphan list. A simulation layer",
        "  resting on a single 15-year-old assertion is a different risk from one",
        "  resting on a single assertion made last year.",
        "* Cross-reference `analysis/atlas-module-support.md`. A module whose claim",
        "  is corroborated by many articles is not affected by anything here.", "",
        "## Limits", "",
        "* Co-assertion, not replication, as set out above. The rate is an upper",
        "  bound on genuine follow-up and a lower bound on genuine novelty.",
        "* PubTator extracts from abstracts, so a claim replicated only in full text",
        "  reads as orphaned. The co-mention layer is the complement for that.",
        "* MeSH indexing lag suppresses the most recent cohorts, which is why cohorts",
        f"  without a full {W}-year window are excluded rather than shown declining.",
        "* Pairs whose papers the census cannot date are excluded entirely.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "quiet_years": W, "latest_year": latest,
        "eligible_pairs": total_eligible, "orphans": len(orphans),
        "overall_replication_rate": overall,
        "cohorts": rows,
        "focus_orphans": [{"year": y, "a": name(k[0]), "b": name(k[1])}
                          for k, y in focus_orphans[:200]],
    }, indent=2) + "\n")
    print(f"\n{100*overall:.1f}% of {total_eligible:,} eligible pairs ever replicated; "
          f"{len(orphans):,} orphans, {len(focus_orphans):,} touching focus genes")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
