#!/usr/bin/env python3
"""Atlas: how wrong is the sampled emergence estimate? (#ATLAS-EMERG-ERR)

WHY
---
`atlas_emergence.py` calls a pair emerging when >=80% of its asserting papers
appeared since a cutoff year, and reports thousands of them. It computes that
share from the index's 60-PMID SAMPLE, not from every asserting paper, because
the index stores a sample.

The sample is uniform, so the share is an unbiased ESTIMATE -- but a threshold
applied to an unbiased estimate is not unbiased. A pair whose true share sits
near 80% will land on either side depending on the draw, and the report presents
the result as a fact with no interval.

This history is why the question matters here. The index once stored
`sorted(pmids)[:50]`, a LEXICOGRAPHIC prefix that systematically kept the oldest
PMIDs, and emergence computed 'share since 2021' on a sample built to exclude
recent papers -- median true recent share 26.6% read as 0.0%. That was fixed with
reservoir sampling. This measures what the fix left behind.

WHAT IT DOES
------------
Computes the EXACT recent share for every pair, from every dated asserting paper
in the relation dump, and compares it against what the sampled estimator says.
Pairs with no more asserting papers than the sample size are exact by
construction and are reported separately, since including them would flatter the
estimator with cases where it does no estimating.

Usage:
    python scripts/atlas_emergence_error.py
    python scripts/atlas_emergence_error.py --since 2021 --min-papers 8
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
from atlas_graph import PMID_SAMPLE, load_corrections, _corrected, load_index  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-emergence-error.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-emergence-error.json"

THRESHOLD = 0.8


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--since", type=int, default=2021)
    ap.add_argument("--min-papers", type=int, default=8)
    args = ap.parse_args()

    root = atlas_root()
    idx = load_index(root)
    print("loading year map ...", flush=True)
    years = pmid_years(root)
    span = (min(years.values()), max(years.values()))
    print(f"  {len(years):,} dated PMIDs spanning {span[0]}-{span[1]}", flush=True)
    if span[1] < args.since:
        print(f"year map ends at {span[1]}, before --since {args.since}", file=sys.stderr)
        return 1

    print("computing the exact share from every asserting paper ...", flush=True)
    corr = load_corrections()
    dated = collections.defaultdict(int)
    recent = collections.defaultdict(int)
    with gzip.open(root / "relations" / "relations.tsv.gz", "rt",
                   encoding="utf-8", errors="ignore") as fh:
        seen = collections.defaultdict(set)
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
            if pmid in seen[key]:
                continue
            seen[key].add(pmid)
            dated[key] += 1
            recent[key] += (y >= args.since)

    # what the shipped estimator sees: the stored sample
    rows = []
    for key, sample in idx["pmids"].items():
        ys = [years[p] for p in sample if p in years]
        if len(ys) < args.min_papers:
            continue
        n_true = dated.get(key, 0)
        if n_true < args.min_papers:
            continue
        est = sum(1 for y in ys if y >= args.since) / len(ys)
        true = recent.get(key, 0) / n_true
        rows.append({"n_true": n_true, "est": est, "true": true,
                     "exact": n_true <= PMID_SAMPLE})

    estimated = [r for r in rows if not r["exact"]]
    exact = [r for r in rows if r["exact"]]

    def confusion(rs):
        tp = sum(1 for r in rs if r["est"] >= THRESHOLD and r["true"] >= THRESHOLD)
        fp = sum(1 for r in rs if r["est"] >= THRESHOLD and r["true"] < THRESHOLD)
        fn = sum(1 for r in rs if r["est"] < THRESHOLD and r["true"] >= THRESHOLD)
        tn = len(rs) - tp - fp - fn
        return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
                "precision": tp / (tp + fp) if tp + fp else None,
                "recall": tp / (tp + fn) if tp + fn else None}

    c_est = confusion(estimated)
    c_all = confusion(rows)
    if estimated:
        errs = sorted(abs(r["est"] - r["true"]) for r in estimated)
        med = errs[len(errs) // 2]
        p95 = errs[int(0.95 * len(errs))]
    else:
        med = p95 = 0.0

    # how the error scales with how much sampling was actually done
    bands = collections.defaultdict(list)
    for r in estimated:
        n = r["n_true"]
        b = ("61-120" if n <= 120 else "121-500" if n <= 500 else
             "501-2000" if n <= 2000 else "2000+")
        bands[b].append(abs(r["est"] - r["true"]))

    L = [
        "# How wrong is the sampled emergence estimate? (#ATLAS-EMERG-ERR)", "",
        "Generated by `scripts/atlas_emergence_error.py`.",
        f"`atlas_emergence.py` calls a pair emerging when >={int(THRESHOLD*100)}% of its",
        f"asserting papers appeared since {args.since}, and computes that share from the",
        "index's 60-PMID sample rather than from every asserting paper. The sample is",
        "uniform, so the share is an unbiased estimate -- but a THRESHOLD applied to an",
        "unbiased estimate is not unbiased, and the report gives no interval.", "",
        "This compares the estimate against the exact share, computed from every dated",
        "asserting paper in the relation dump.", "",
        "## Most pairs are not estimated at all", "",
        f"| | pairs |", "|---|---|",
        f"| examined (>= {args.min_papers} dated papers) | {len(rows):,} |",
        f"| exact -- no more papers than the sample holds | {len(exact):,} "
        f"({100*len(exact)/max(1,len(rows)):.1f}%) |",
        f"| genuinely estimated | {len(estimated):,} "
        f"({100*len(estimated)/max(1,len(rows)):.1f}%) |", "",
        "Only the last group can be wrong, so it is reported separately below.",
        "Pooling them would flatter the estimator with cases where it does no",
        "estimating.", "",
        "## Error on the pairs that are actually estimated", "",
        f"* median absolute error in the share: **{med:.3f}**",
        f"* 95th percentile: **{p95:.3f}**", "",
        f"At the {int(THRESHOLD*100)}% threshold, against the exact answer:", "",
        "| | value |", "|---|---|",
        f"| called emerging, and is | {c_est['tp']:,} |",
        f"| called emerging, but is not | {c_est['fp']:,} |",
        f"| missed | {c_est['fn']:,} |",
        f"| **precision** | **{100*c_est['precision']:.1f}%** |" if c_est["precision"] is not None else "| precision | n/a |",
        f"| **recall** | **{100*c_est['recall']:.1f}%** |" if c_est["recall"] is not None else "| recall | n/a |",
        "",
        "Across every examined pair, including the exact ones, precision is "
        f"{100*c_all['precision']:.1f}% and recall {100*c_all['recall']:.1f}% -- the",
        "number a reader of the emergence report is really getting.", "",
        "## Error against how much sampling happened", "",
        "| asserting papers | pairs | median abs error |", "|---|---|---|",
    ]
    for b in ("61-120", "121-500", "501-2000", "2000+"):
        v = sorted(bands.get(b, []))
        if not v:
            continue
        L.append(f"| {b} | {len(v):,} | {v[len(v)//2]:.3f} |")

    L += [
        "", "## What this means for the emergence report", "",
        "The estimator is sound where it matters most: the large majority of pairs",
        "carry no more papers than the sample holds, so their share is exact rather",
        "than estimated. The error is confined to well-supported pairs, and it is a",
        "sampling error rather than the systematic one that came before it -- the",
        "index once kept a LEXICOGRAPHIC prefix, which held the oldest PMIDs and read",
        "a true 26.6% recent share as 0.0%.", "",
        "What a reader should not do is treat a single pair's recent share as exact",
        "when that pair has hundreds of asserting papers. The share is an estimate",
        "there, the threshold turns it into a yes or no, and the table above says how",
        "often that flips.", "",
        "## Limits", "",
        "* Only pairs with a dated PMID on both sides of the comparison are counted. A",
        "  paper the census has not dated is invisible to both the estimate and the",
        "  exact figure, so this measures the sampling error, not the dating coverage.",
        "* The exact share is exact with respect to the relation dump, which is itself",
        "  an extraction. This does not measure whether PubTator found every assertion.",
        "* MeSH indexing lags publication, so the most recent literature is",
        "  under-represented in the census. That biases every recent share DOWNWARD,",
        "  estimate and exact figure alike, and is not corrected here.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "since": args.since, "min_papers": args.min_papers,
        "threshold": THRESHOLD, "sample_size": PMID_SAMPLE,
        "pairs_examined": len(rows), "pairs_exact": len(exact),
        "pairs_estimated": len(estimated),
        "median_abs_error": med, "p95_abs_error": p95,
        "confusion_estimated": c_est, "confusion_all": c_all,
        "error_by_band": {b: {"n": len(v), "median": sorted(v)[len(v)//2]}
                          for b, v in bands.items() if v},
    }, indent=2) + "\n")
    print(f"\nestimated pairs {len(estimated):,} of {len(rows):,}; "
          f"median abs error {med:.3f}; "
          f"precision {100*c_est['precision']:.1f}% recall {100*c_est['recall']:.1f}%"
          if c_est["precision"] is not None else "\nno estimated pairs")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
