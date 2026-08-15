#!/usr/bin/env python3
"""Atlas: find what the cancer literature has only recently started saying (#ATLAS).

WHY
---
A census is only worth building if it answers questions a sample cannot. This is
one of them: which entity relations are NEW?

Every paper's year is known, and every relation carries the PMIDs asserting it,
so for any pair the graph gives a publication-year distribution. A relation whose
support is concentrated in the last few years is an emerging claim; one spread
evenly across decades is settled background. A sample of 4,830 articles cannot
distinguish these, because it lacks both the volume and the time depth.

WHAT IT REPORTS
---------------
For each pair with enough support, the median year of its asserting papers and
the share published in the recent window. Ranked by recent share, then volume,
so a pair with 40 papers all since 2022 outranks one with 8.

LIMITS
------
Emergence is not importance. A relation can be new because the entity was only
recently named (a 2023 drug cannot have 1990s support), because a technique made
it newly measurable, or because a field is faddish. The report gives the year
distribution so a reader can tell those apart, not a verdict.

The PMID list per pair is capped when the index is built, so `recent_share` is
computed on a bounded sample of each pair's support, not all of it. Pairs with
very large support are therefore estimated, and the cap is reported.

MEASURED ACCURACY
-----------------
The recent share is computed from the index's 60-PMID sample, not from every
asserting paper, so for well-supported pairs it is an ESTIMATE and the >=80%
rule turns that estimate into a yes or no. `scripts/atlas_emergence_error.py`
compares it against the exact share computed from every dated asserting paper:

  * 89.4% of examined pairs carry no more papers than the sample holds, so their
    share is exact rather than estimated;
  * on the 10.6% genuinely estimated, median absolute error in the share is
    0.017 and the threshold decision is 86.4% precise at 93.2% recall;
  * across every examined pair, precision 99.0% and recall 99.6% -- the figure a
    reader of this report is actually getting.

Error grows with support, as sampling error should (median 0.013 at 61-120
asserting papers, 0.035 above 2,000). So do not read a single pair's share as
exact when that pair has hundreds of papers behind it.

Usage:
    python scripts/atlas_emergence.py
    python scripts/atlas_emergence.py --since 2021 --min-papers 8 --top 40
"""

import argparse
import collections
import glob
import gzip
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_graph import load_index  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-emergence.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-emergence.json"


# Every stream that carries a `year`. NOT a cosmetic list: the two names this
# tuple used to hold were `records` and `records_c04only`, and c04only is a
# strict SUBSET of records -- measured, it adds exactly zero PMIDs -- so a loop
# that read like a merge merged nothing, while the two streams holding the most
# recent literature were never named. The layer whose subject is recency was
# blind to 814,015 dated articles, overwhelmingly from 2021 onward.
YEAR_STREAMS = ("records", "records_c04only", "records_unindexed", "records_updates")


def pmid_years(root: Path) -> tuple:
    """PMID -> publication year, merged across every census directory present.

    MUST merge rather than pick one. PubMed baseline files are ordered
    CHRONOLOGICALLY, so a partially-rebuilt census contains only the OLDEST
    literature, and dating relations against it makes every pair look ancient
    and reports zero emerging claims. Merging all available sources, and
    reporting the resulting year span, makes that failure visible instead of
    silent.
    """
    years = {}
    used = []
    for d in YEAR_STREAMS:
        files = sorted(glob.glob(str(root / d / "*.jsonl.gz")))
        if not files:
            continue
        used.append(f"{d} ({len(files)} shards)")
        for f in files:
            with gzip.open(f, "rt", encoding="utf-8") as fh:
                for line in fh:
                    r = json.loads(line)
                    y = r.get("year")
                    if y:
                        years[r["pmid"]] = y
    return years, ", ".join(used) if used else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--since", type=int, default=2021, help="start of the 'recent' window")
    ap.add_argument("--min-papers", type=int, default=6,
                    help="minimum sampled asserting papers with a known year")
    ap.add_argument("--top", type=int, default=40)
    args = ap.parse_args()

    root = atlas_root()
    idx = load_index(root)
    years, source = pmid_years(root)
    if not years:
        raise SystemExit("no census records found; run scripts/atlas_baseline.py first")
    span = (min(years.values()), max(years.values()))
    print(f"year map: {len(years):,} PMIDs spanning {span[0]}-{span[1]} (from {source})",
          flush=True)
    if span[1] < args.since:
        raise SystemExit(
            f"the year map ends at {span[1]}, before the recent window opens at "
            f"{args.since}. The census is probably still rebuilding: baseline files are "
            "ordered chronologically, so a partial census holds only the oldest "
            "literature and every pair would look settled. Wait for the ingest to "
            "finish, or pass --since below that year.")

    rows = []
    for key, pmids in idx["pmids"].items():
        ys = [years[p] for p in pmids if p in years]
        if len(ys) < args.min_papers:
            continue
        recent = sum(1 for y in ys if y >= args.since)
        preds = idx["edges"].get(key, {})
        rows.append(dict(
            a=idx["canon"].get(key[0], key[0]), b=idx["canon"].get(key[1], key[1]),
            n_sampled=len(ys), median_year=int(statistics.median(ys)),
            min_year=min(ys), max_year=max(ys),
            recent=recent, recent_share=recent / len(ys),
            total_relations=sum(preds.values()),
            predicates=dict(sorted(preds.items(), key=lambda kv: -kv[1])[:3]),
            pmids=[p for p in pmids if p in years][:5]))

    rows.sort(key=lambda r: (-r["recent_share"], -r["n_sampled"]))
    emerging = [r for r in rows if r["recent_share"] >= 0.8]
    RAW.write_text(json.dumps(rows[:1000], indent=1), encoding="utf-8")

    L = [
        "# What the cancer literature has only recently started saying (#ATLAS)", "",
        "Generated by `scripts/atlas_emergence.py`.", "",
        f"Every relation carries the PMIDs asserting it and every paper carries a year, so",
        f"each pair has a publication-year distribution. Pairs whose support is",
        f"concentrated since **{args.since}** are emerging claims; pairs spread across",
        "decades are settled background. A 4,830-article sample cannot tell these apart --",
        "it has neither the volume nor the time depth.", "",
        f"Scanned {len(rows):,} pairs with at least {args.min_papers} dated asserting papers. "
        f"**{len(emerging):,}** have {'≥'}80% of their support in the recent window.", "",
        "## Limits", "",
        "Emergence is not importance. A relation can be new because the entity was only",
        "recently named (a 2023 drug cannot have 1990s support), because a technique made",
        "it newly measurable, or because a field is faddish. The year range is shown so a",
        "reader can tell those apart.",
        "",
        "The PMID list per pair is capped when the index is built, so `recent share` is",
        "computed on a bounded sample of each pair's support rather than all of it. Pairs",
        "with very large support are estimated.", "",
        f"## Most recently-emerged relations", "",
        "| pair | sampled papers | years | recent share | predicates |",
        "|---|---|---|---|---|",
    ]
    for r in rows[:args.top]:
        preds = ", ".join(f"{k} {v}" for k, v in r["predicates"].items())
        L.append(f"| {r['a']} — {r['b']} | {r['n_sampled']} | {r['min_year']}-{r['max_year']} | "
                 f"{r['recent_share']:.0%} | {preds} |")

    oldest = sorted(rows, key=lambda r: r["median_year"])[:15]
    L += ["", "## For contrast: the most settled relations", "",
          "| pair | sampled papers | years | median year | predicates |",
          "|---|---|---|---|---|"]
    for r in oldest:
        preds = ", ".join(f"{k} {v}" for k, v in r["predicates"].items())
        L.append(f"| {r['a']} — {r['b']} | {r['n_sampled']} | {r['min_year']}-{r['max_year']} | "
                 f"{r['median_year']} | {preds} |")

    L += ["", "## How to use this", "",
          "An emerging relation with real volume is a candidate for the simulation suite to",
          "test, and a candidate the manuscript's literature chapters would have missed:",
          "the frozen corpus's newest-mechanism queries were date-capped and 500-record",
          "capped, so it systematically under-samples exactly this.", ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"{len(rows):,} pairs scanned, {len(emerging):,} emerging (>=80% since {args.since})")
    for r in rows[:12]:
        print(f"  {r['a'][:26]:<28}{r['b'][:24]:<26}n={r['n_sampled']:<4}"
              f"{r['min_year']}-{r['max_year']}  {r['recent_share']:.0%}")


if __name__ == "__main__":
    main()
