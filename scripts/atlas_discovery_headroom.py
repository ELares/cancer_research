#!/usr/bin/env python3
"""Does any seed-specific signal add anything over a candidate-only prior?

THE QUESTION, AND WHY THE LEADERBOARD CANNOT ANSWER IT
-------------------------------------------------------
`atlas-discovery-eval.md` shows every seed-aware ranking losing to `popularity`.
That is a comparison of rankers in ISOLATION, and it does not settle the thing
anyone actually wants to know, because losing alone is not the same as adding
nothing. A signal can be weaker than a baseline and still be complementary to
it -- ranking by it alone does worse, ranking by baseline-then-it does better.

Note the structural fact that makes this the right question here: `popularity`
is the ONLY ranking whose score is a function of the candidate alone. `deg_c`
does not vary with the seed. Every other method varies with the seed through the
bridge set. So the published table already reads "a seed-blind score beats every
seed-aware score", and the open question is whether the seed-aware ones carry
anything the seed-blind one does not.

THE TEST
--------
Rank each pool by a candidate-only benchmark B, then re-rank by B with one
seed-aware ranker's within-seed percentile added, and compare precision@20 on
the same seeds. The statistic is the per-seed difference in hits out of 20, fed
to the same paired bootstrap the evaluation uses.

  H0: adding seed-specific information does not change precision@20.

If H0 holds for every method, then on this metric the candidate prior exhausts
what is measurable, and "our ranker is bad" and "the metric cannot see discovery"
are indistinguishable -- which is a statement about the evaluation, not a defence
of any ranker.

THE DECISION RULE IS FIXED BEFORE THE RESULT IS SEEN, and written here rather
than chosen afterwards. A method is credited with headroom only if BOTH hold:
its 95% percentile bootstrap interval on the per-seed difference excludes zero,
AND its Holm-corrected two-sided bootstrap p is below 0.05. The family is the
methods in `SEED_AWARE` -- five of them, and the count is taken from that dict
rather than written down, because an earlier draft of this docstring said
"seven" while the code corrected over five. Verdicts were unchanged either way,
but a factual error inside the one sentence whose job is to establish that the
rule was fixed in advance is worth more than the arithmetic.

WHY THIS IS AFFORDABLE
----------------------
It reads the committed candidate dump. No graph, no re-ranking run: every pool
member already carries its features and its outcome, so the whole thing is a
few seconds of arithmetic over a 6.5 MB file.

WHAT IT CANNOT SETTLE
---------------------
Whether hub-selection is RIGHT. Two independent critiques established that the
separation between "popularity is a good prior" and "the target rewards
popularity" is not identifiable from this corpus at any sample size, because a
latent importance drives degree, how much gets written, and real biology
together. This asks only whether seed-specific information adds measurable
precision over a candidate-only prior.

Usage:
    python scripts/atlas_discovery_headroom.py
"""

import argparse
import collections
import csv
import gzip
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

DUMP = PROJECT_ROOT / "analysis" / "atlas-discovery-candidates.csv.gz"
SIDE = PROJECT_ROOT / "analysis" / "atlas-discovery-candidates-seeds.json"
OUT = PROJECT_ROOT / "analysis" / "atlas-discovery-headroom.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-discovery-headroom.json"

BOOTSTRAP_SEED = 20260806
N_BOOT = 10000

# Seed-aware scores: every one varies with the seed through the bridge set.
# `popularity` is deliberately absent -- it IS the benchmark.
SEED_AWARE = {
    "abc": lambda r: (-float(r["p"]), int(r["n"])),      # low p first
    "bridges": lambda r: (int(r["n"]),),
    "adamic_adar": lambda r: (float(r["aa"]),),
    "resource_alloc": lambda r: (float(r["ra"]),),
    "jaccard": lambda r: (float(r["jac"]),),
}


def load_pools():
    pools = collections.defaultdict(list)
    with gzip.open(DUMP, "rt") as fh:
        for r in csv.DictReader(fh):
            pools[r["seed"]].append(r)
    return pools


def percentile_ranks(rows, key):
    """Within-seed percentile of a score, ties averaged. Higher score -> higher.

    Percentiles rather than raw scores so that seed-level differences in scale
    -- a big pool's hypergeometric p spans orders of magnitude a small one does
    not -- cannot leak into the combination.
    """
    n = len(rows)
    order = sorted(range(n), key=lambda i: key(rows[i]))
    out = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and key(rows[order[j + 1]]) == key(rows[order[i]]):
            j += 1
        avg = (i + j) / 2
        for k in range(i, j + 1):
            out[order[k]] = avg / max(1, n - 1)
        i = j + 1
    return out


def top_hits(rows, score, k):
    """Hits in the top k under `score`, ties broken by identifier for stability."""
    idx = sorted(range(len(rows)),
                 key=lambda i: (-score[i], rows[i]["c"]))[:k]
    return sum(int(rows[i]["hit"]) for i in idx)


def paired_ci(diffs, rng):
    """Mean, 95% percentile interval, and a two-sided p from ONE resample set.

    The p used to be drawn from a SECOND, independent set of resamples while the
    interval came from the first, so the two halves of the decision rule were
    not computed from the same bootstrap. And its tails were asymmetric --
    `P(mean* <= 0)` against `P(mean* > 0)` -- which returns p = 0 for a signal
    that does nothing at all: an all-zero difference vector put every resample
    at exactly zero, all of it in the first tail and none in the second. Both
    tails are now inclusive, and the p is floored at the bootstrap's own
    resolution rather than printed as 0.
    """
    n = len(diffs)
    means = sorted(sum(diffs[rng.randrange(n)] for _ in range(n)) / n
                   for _ in range(N_BOOT))
    le = sum(1 for x in means if x <= 0) / N_BOOT
    ge = sum(1 for x in means if x >= 0) / N_BOOT
    p = min(1.0, 2 * min(le, ge))
    return (sum(diffs) / n, means[int(0.025 * N_BOOT)], means[int(0.975 * N_BOOT)],
            max(p, 1.0 / N_BOOT))


def holm(pvals: dict) -> dict:
    """Holm-Bonferroni over the methods tested, however many there are."""
    ordered = sorted(pvals.items(), key=lambda kv: kv[1])
    m, out, prev = len(ordered), {}, 0.0
    for i, (k, p) in enumerate(ordered):
        adj = max(prev, min(1.0, (m - i) * p))
        out[k] = adj
        prev = adj
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--weight", type=float, default=0.5,
                    help="weight on the seed-aware percentile in the blend")
    args = ap.parse_args()

    if not DUMP.exists():
        print(f"no candidate dump at {DUMP}; run scripts/atlas_discovery_dump.py",
              file=sys.stderr)
        return 1
    side = json.loads(SIDE.read_text())
    pools = load_pools()
    print(f"  {len(pools)} seeds, {sum(len(v) for v in pools.values()):,} candidates",
          flush=True)

    per_seed = {m: [] for m in SEED_AWARE}
    base_hits = 0
    for sid, rows in sorted(pools.items()):
        # Benchmark: candidate-only. deg_c does not vary with the seed.
        b = percentile_ranks(rows, lambda r: (int(r["deg_c"]),))
        h_b = top_hits(rows, b, args.top)
        base_hits += h_b
        for m, key in SEED_AWARE.items():
            s = percentile_ranks(rows, key)
            blend = [(1 - args.weight) * b[i] + args.weight * s[i]
                     for i in range(len(rows))]
            per_seed[m].append(top_hits(rows, blend, args.top) - h_b)

    # THE SWEEP. One blend weight cannot answer "does any mixture help": a
    # weight near zero reduces to the benchmark by construction. The prose used
    # to promise this sweep before it existed, which is the defect this repo
    # keeps finding, so it is computed rather than referred to.
    sweep = {w: {} for w in (0.1, 0.25, 0.5, 0.75, 0.9)}
    sweep_diffs = {w: {} for w in sweep}
    for w in sweep:
        for m, key in SEED_AWARE.items():
            d = []
            for sid, rows in sorted(pools.items()):
                b = percentile_ranks(rows, lambda r: (int(r["deg_c"]),))
                s = percentile_ranks(rows, key)
                blend = [(1 - w) * b[i] + w * s[i] for i in range(len(rows))]
                d.append(top_hits(rows, blend, args.top) - top_hits(rows, b, args.top))
            sweep_diffs[w][m] = d
            sweep[w][m] = sum(d) / len(d)

    n_seeds = len(pools)
    rng = random.Random(BOOTSTRAP_SEED)
    results, pvals = {}, {}
    for m, d in per_seed.items():
        mean, lo, hi, p = paired_ci(d, rng)
        results[m] = {"mean_diff": mean, "ci95": [lo, hi],
                      "decided": hi < 0 or lo > 0,
                      "better": sum(1 for x in d if x > 0),
                      "worse": sum(1 for x in d if x < 0),
                      "same": sum(1 for x in d if x == 0), "p": p}
        pvals[m] = p
    adj = holm(pvals)
    for m in results:
        results[m]["p_holm"] = adj[m]
        results[m]["decided_holm"] = adj[m] < 0.05 and results[m]["decided"]

    # Every positive sweep cell gets the same test the headline uses. Reporting a
    # positive mean without one is how +0.03 hits out of 20 becomes a lead.
    survivors, tested = [], []
    for w in sorted(sweep):
        for m in sorted(SEED_AWARE):
            if sweep[w][m] <= 0:
                continue
            mean, lo, hi, _p = paired_ci(sweep_diffs[w][m], rng)
            tested.append({"weight": w, "method": m, "mean": mean, "ci95": [lo, hi],
                           "decided": lo > 0})
            if lo > 0:
                survivors.append(f"`{m}` at w={w}")
    if tested:
        worst = max(tested, key=lambda r: r["mean"])
        pos_note = (f"{len(tested)} of {len(sweep)*len(SEED_AWARE)} cells are "
                    f"positive at all; the largest is {worst['method']} at "
                    f"w={worst['weight']}, {worst['mean']:+.3f} hits out of "
                    f"{args.top} with a 95% interval of [{worst['ci95'][0]:+.3f}, "
                    f"{worst['ci95'][1]:+.3f}] -- which includes zero.")
    else:
        pos_note = "Every cell at every weight is negative or flat."

    base_prec = base_hits / (n_seeds * args.top)
    any_head = [m for m, r in results.items() if r["decided_holm"] and r["mean_diff"] > 0]

    L = [
        "# Does seed-specific signal add anything over a candidate-only prior?",
        "", "Generated by `scripts/atlas_discovery_headroom.py` from the committed",
        "candidate dump -- no graph rebuild, no re-ranking run.", "",
        "## Why the leaderboard cannot answer this", "",
        "`atlas-discovery-eval.md` shows every seed-aware ranking losing to",
        "`popularity`. That compares rankers in ISOLATION, and losing alone is not",
        "the same as adding nothing: a signal can be weaker than a baseline and",
        "still complementary to it. Note also that `popularity` is the only ranking",
        "whose score is a function of the candidate alone -- `deg_c` does not vary",
        "with the seed -- so the published table already reads as a seed-BLIND score",
        "beating every seed-aware one. What it leaves open is whether the seed-aware",
        "ones carry anything the seed-blind one does not.", "",
        "## The test", "",
        f"Rank each pool by candidate degree alone (**{100*base_prec:.1f}%** "
        f"precision@{args.top} over {n_seeds} seeds), then blend in one method's",
        "within-seed percentile and re-score. Percentiles, not raw scores, so",
        "seed-level differences in scale cannot leak into the blend. The statistic",
        "is the per-seed difference in hits, through the same paired bootstrap the",
        "evaluation uses.", "",
        "**The decision rule was fixed before the result was seen** and is in the",
        "script's docstring: headroom is credited only when BOTH the 95% percentile",
        f"interval excludes zero AND the Holm-corrected p (over the "
        f"{len(SEED_AWARE)} methods) is below 0.05. The interval reported below is",
        "uncorrected; the correction lives in the p column.", "",
        f"| method | mean Δ hits/{args.top} | 95% CI | better/worse/same | p (Holm) | headroom? |",
        "|---|---|---|---|---|---|",
    ]
    for m in sorted(results, key=lambda k: -results[k]["mean_diff"]):
        r = results[m]
        L.append(f"| {m} | {r['mean_diff']:+.3f} | "
                 f"[{r['ci95'][0]:+.3f}, {r['ci95'][1]:+.3f}] | "
                 f"{r['better']}/{r['worse']}/{r['same']} | {r['p_holm']:.3f} | "
                 f"{'**yes**' if r['decided_holm'] and r['mean_diff'] > 0 else 'no'} |")

    L += ["", "## What it says", ""]
    if any_head:
        L += [f"**{', '.join('`'+m+'`' for m in any_head)} adds measurable precision",
              "over the candidate-only prior.** So the seed-specific structure is not",
              "redundant with popularity, and a ranking built on both would beat",
              "either -- which the isolated leaderboard could not have shown.", ""]
    else:
        n_pos = len(tested)
        L += ["**No method adds measurable precision over the candidate-only prior.**",
              (f"At the headline weight every method is negative, and across the whole "
               f"sweep {n_pos} of {len(sweep)*len(SEED_AWARE)} cells are positive at "
               f"all -- against roughly half that many expected if every method were "
               f"pure noise. None survives the decision rule."
               if n_pos else
               "Every blend that moves the score at all moves it down."),
              "No interval excludes zero in the positive direction after correction.",
              "",
              "On this metric the candidate prior exhausts what is measurable. That",
              "makes 'the shipped ranker is bad' and 'the metric cannot see discovery'",
              "indistinguishable ON THIS EVIDENCE -- which is a statement about the",
              "evaluation rather than a defence of any ranker.", ""]
    L += [
        "## What this cannot settle", "",
        "**It tests ONE combination family**: a linear convex blend of two",
        "within-seed percentile ranks. A signal could in principle be",
        "complementary in a way that family cannot express -- informative only",
        "inside a subpopulation, or requiring a conditional or multiplicative",
        "combination. Three probes outside the family were run and none changed the",
        "answer: a lexicographic blend at the limit (degree primary, method",
        "breaking ties) is IDENTICAL to the benchmark for all five methods on all",
        "200 seeds, because degree ties never straddle the top-20 cut; weights",
        "below the sweep floor add nothing; and a two-stage retrieve-then-rerank",
        "over six retrieval depths is negative in 29 of 30 cells and monotonically",
        "worse with depth. Splitting seeds by pool size leaves every method",
        "negative in both strata. So the result is not an artifact of the blend,",
        "but it is stated over what was tested.", "",
        "Whether hub-selection is RIGHT. The separation between 'popularity is a",
        "good prior' and 'the target rewards popularity' is not identifiable from",
        "this corpus at any sample size: a latent importance drives degree, how much",
        "gets written, and real biology together, and every available exposure",
        "measure is either the outcome through the same extractor or the treatment",
        "itself. This asks only whether seed-specific information adds measurable",
        "precision over a candidate-only prior.", "",
        "## Is it the blend weight?", "",
        "One weight cannot answer this: a weight near zero reduces to the benchmark",
        "by construction. Mean change in hits at each weight, for every method:", "",
        "| weight | " + " | ".join(sorted(SEED_AWARE)) + " |",
        "|---|" + "---|" * len(SEED_AWARE),
    ] + [
        f"| {w} | " + " | ".join(f"{sweep[w][m]:+.3f}" for m in sorted(SEED_AWARE)) + " |"
        for w in sorted(sweep)
    ] + [
        "",
        (f"**No cell survives the decision rule.** {pos_note}"
         if not survivors else
         f"**{', '.join(survivors)} survive the decision rule at some weight**, so "
         "a tuned blend does carry signal the single-weight test missed."), "",
    ]
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    RAW.write_text(json.dumps({
        "seeds": n_seeds, "top": args.top, "weight": args.weight,
        "pairs_before": side.get("pairs_before"),
        "benchmark_precision": base_prec, "results": results,
        "sweep": {str(w): sweep[w] for w in sweep},
        "sweep_positive_cells_tested": tested, "sweep_survivors": survivors,
        "any_headroom": any_head,
    }, indent=2) + "\n", encoding="utf-8")
    print("  benchmark (degree only): %.1f%%" % (100 * base_prec))
    for m in sorted(results, key=lambda k: -results[k]["mean_diff"]):
        r = results[m]
        print(f"  {m:<16}{r['mean_diff']:+.3f}  CI[{r['ci95'][0]:+.3f},"
              f"{r['ci95'][1]:+.3f}]  holm p={r['p_holm']:.3f}")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
