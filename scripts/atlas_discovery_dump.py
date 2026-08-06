#!/usr/bin/env python3
"""Dump every discovery candidate with its features and its outcome.

WHY THIS EXISTS
---------------
`atlas_discovery_eval.py` builds a candidate pool per seed, ranks it seven ways,
and scores only the top 20 of each. Everything below rank 20 is computed and
thrown away, and the outcome label -- did the literature go on to assert this
pair -- is only ever evaluated for those 140 rows per seed.

That makes every follow-up question require a fresh multi-minute run over a
2.2 GB graph: any re-ranking, any stratification, any "what if we scored it this
way instead". Dumping the pool ONCE with its features and its hit label makes
all of those offline and instant, and it fits the repo's contract that CI and
downstream analyses read committed artifacts rather than rebuilding the census.

WHAT IS DUMPED
--------------
One row per (seed, candidate) over the SAME split, seeds and pools the precision
table uses: the four ranking features (`n` bridges, `deg_c`, hypergeometric `p`,
Adamic-Adar, resource allocation, Jaccard) and `hit`, whether that pair is first
asserted in the split year or later. Per-seed constants that the features cannot
be recomputed without -- `deg_a_usable` and the split's `n_nodes` -- go in a
sidecar, because `degrees[seed]` is the UNFILTERED degree and using it to
re-derive `p` silently gives a different number.

WHY IT RECOMPUTES THE FEATURES AND THEN CHECKS ITSELF
------------------------------------------------------
`rank_all` returns orderings, not the feature rows behind them, so this rebuilds
the pool with the same construction. Rebuilding invites divergence, so it does
not ask to be trusted: for every seed it asserts that ranking its own rows
reproduces `rank_all`'s ordering EXACTLY, for all seven methods, ties included.
If the reconstruction drifts, the dump refuses to write rather than shipping
plausible numbers.

Usage:
    python scripts/atlas_discovery_dump.py                # 2018 split, 200 seeds
    python scripts/atlas_discovery_dump.py --seeds 25     # quick check
"""

import argparse
import collections
import csv
import gzip
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_discovery import (  # noqa: E402
    HUB_PERCENTILE, MIN_BRIDGES, MIN_CANDIDATE_DEGREE,
)
from atlas_discovery_eval import (  # noqa: E402
    SEED_DEGREE_MAX, SEED_DEGREE_MIN, SEED_RNG, load_corrections, load_index,
    pair_first_year, pmid_years, rank_all,
)
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-discovery-candidates.csv.gz"
SIDE = PROJECT_ROOT / "analysis" / "atlas-discovery-candidates-seeds.json"

FIELDS = ["seed", "c", "n", "deg_c", "p", "aa", "ra", "jac", "hit"]


def build_pool(adj, degrees, cutoff, n_nodes, sid):
    """The candidate pool for one seed, with every ranking feature.

    Mirrors `rank_all`'s construction exactly, including the sort by identifier
    that makes ties reproducible under salted string hashing.
    """
    from scipy.stats import hypergeom

    a_nb = adj.get(sid, set())
    usable = {b for b in a_nb if degrees.get(b, 0) <= cutoff}
    if not usable:
        return None, 0
    deg_a = len(usable)

    bridges = collections.defaultdict(set)
    for b in usable:
        for c in adj.get(b, ()):
            if c == sid or c in a_nb:
                continue
            bridges[c].add(b)

    cands = []
    for c, bs in bridges.items():
        n = len(bs)
        if n < MIN_BRIDGES:
            continue
        deg_c = degrees.get(c, 0)
        if deg_c < MIN_CANDIDATE_DEGREE:
            continue
        aa = sum(1.0 / math.log(max(degrees.get(b, 1), 2)) for b in bs)
        ra = sum(1.0 / max(degrees.get(b, 1), 1) for b in bs)
        union = deg_a + deg_c - n
        cands.append({"c": c, "n": n, "deg_c": deg_c,
                      "p": float(hypergeom.sf(n - 1, n_nodes, deg_a, deg_c)),
                      "aa": aa, "ra": ra,
                      "jac": (n / union) if union > 0 else 0.0})
    cands.sort(key=lambda r: r["c"])
    return (cands or None), deg_a


def orderings(cands, rng):
    """The same seven orderings `rank_all` produces, from local rows."""
    order = lambda key: [r["c"] for r in sorted(cands, key=key)]  # noqa: E731
    rnd = [r["c"] for r in cands]
    rng.shuffle(rnd)
    return {
        "abc": order(lambda r: (r["p"], -r["n"])),
        "popularity": order(lambda r: -r["deg_c"]),
        "adamic_adar": order(lambda r: -r["aa"]),
        "resource_alloc": order(lambda r: -r["ra"]),
        "jaccard": order(lambda r: -r["jac"]),
        "bridges": order(lambda r: -r["n"]),
        "random": rnd,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", type=int, default=2018)
    ap.add_argument("--seeds", type=int, default=200)
    args = ap.parse_args()

    root = atlas_root()
    idx = load_index(root)
    first = pair_first_year(root, pmid_years(root), load_corrections())

    before = {k for k, y in first.items() if y < args.split}
    after = {k for k, y in first.items() if y >= args.split}
    adj = collections.defaultdict(set)
    for a, b in before:
        adj[a].add(b)
        adj[b].add(a)
    degrees = {k: len(v) for k, v in adj.items()}
    cutoff = sorted(degrees.values())[int(len(degrees) * HUB_PERCENTILE)]
    n_nodes = len(degrees)
    print(f"  split {args.split}: {len(before):,} pairs before, {len(after):,} after",
          flush=True)

    rng = random.Random(SEED_RNG)
    pool = sorted(k for k, d in degrees.items()
                  if SEED_DEGREE_MIN <= d <= SEED_DEGREE_MAX)
    seeds = rng.sample(pool, min(args.seeds, len(pool)))

    # Two RNGs advanced in lockstep: `rank_all` consumes draws for its random
    # baseline, so a shared generator would desynchronise the check below.
    rng_ref = random.Random(SEED_RNG)
    rng_ref.sample(pool, min(args.seeds, len(pool)))
    rng_mine = random.Random(SEED_RNG)
    rng_mine.sample(pool, min(args.seeds, len(pool)))

    rows, sidecar, checked = [], [], 0
    for sid in seeds:
        ref = rank_all(adj, degrees, cutoff, n_nodes, sid, rng_ref)
        cands, deg_a = build_pool(adj, degrees, cutoff, n_nodes, sid)
        if not ref or not cands:
            continue
        mine = orderings(cands, rng_mine)
        for m in ref:
            assert mine[m] == ref[m], (
                f"reconstruction diverges from rank_all for seed {sid}, "
                f"method {m}; the dump would not describe the evaluated pools")
        checked += 1
        h = 0
        for r in cands:
            key = (sid, r["c"]) if sid <= r["c"] else (r["c"], sid)
            hit = int(key in after)
            h += hit
            rows.append([sid, r["c"], r["n"], r["deg_c"], f"{r['p']:.6g}",
                         f"{r['aa']:.6g}", f"{r['ra']:.6g}", f"{r['jac']:.6g}", hit])
        sidecar.append({"seed": sid, "seed_name": idx["canon"].get(sid, sid),
                        "degree": degrees.get(sid, 0), "deg_a_usable": deg_a,
                        "candidates": len(cands), "pool_hits": h})

    with gzip.open(OUT, "wt", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(FIELDS)
        w.writerows(rows)
    SIDE.write_text(json.dumps({
        "split": args.split, "n_nodes": n_nodes, "hub_cutoff": cutoff,
        # The BUILD fingerprint, not just the split. Split year and k matched
        # between two graphs that differed by 11% in edge count, and a guard
        # comparing only those passed while an analysis quoted selectivity from
        # one build beside precision from the other.
        "pairs_before": len(before), "pairs_after": len(after),
        "seeds_evaluated": len(sidecar), "rows": len(rows), "seeds": sidecar,
    }, indent=2) + "\n", encoding="utf-8")

    hits = sum(r[-1] for r in rows)
    print(f"  {checked} seeds reproduced rank_all exactly (all 7 orderings)")
    print(f"  {len(rows):,} candidate rows, {hits:,} hits "
          f"({100*hits/max(1,len(rows)):.2f}% pool-wide base rate)")
    print(f"wrote {OUT} ({OUT.stat().st_size/1e6:.1f} MB)\nwrote {SIDE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
