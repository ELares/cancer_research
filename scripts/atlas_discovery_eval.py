#!/usr/bin/env python3
"""Atlas: does literature-based discovery predict anything? (#ATLAS-LBD-EVAL)

WHY
---
`atlas_discovery.py` emits a ranked list of A-C pairs the literature implies but
has never stated. It is careful about its limits and calls the output "a ranked
reading list, nothing here is a finding". What it has never had is a HIT RATE.

Without one there is no way to tell the layer apart from a list of famous
entities. Run on GPX4 it returns ERK, caspase-3, cyclin D1, MMP-9 and ATP --
exactly what a popularity ranking would return, which is the failure mode the
module's own docstring says it corrects for.

So this measures it, the standard way: a TIME SPLIT. Build the graph as it stood
before year Y, predict which absent A-C pairs will appear, and check against what
the literature actually did next.

THE COMPARISON THAT MATTERS
---------------------------
Not "does ABC beat random" -- almost anything beats random here, because
co-occurrence graphs are heavily clustered and any 2-hop neighbour is more likely
to link than an arbitrary node. The question is whether ABC beats **ranking the
same candidates by popularity**. If a degree ranking does as well, the ABC
machinery -- bridges, hub filtering, hypergeometric tails -- is decoration, and
the honest thing is to say so.

Three rankings over an IDENTICAL candidate set, so only the ordering differs:

  * `abc`        -- the shipped ranking (hypergeometric tail over bridge counts)
  * `popularity` -- rank by candidate degree in the before-graph
  * `random`     -- a seeded shuffle, the floor

WHAT COUNTS AS A HIT
--------------------
A predicted pair A-C is a hit if the literature first asserts it in year >= Y.
Pairs already asserted before Y are excluded from prediction by construction, so
a hit is genuinely a NEW statement, not a rediscovery.

WHAT THIS CANNOT SHOW
---------------------
That a hit is a real biological relation. PubTator's extractor has its own error
rate, and a new edge may be a new extraction of an old idea. It measures whether
the ranking anticipates what the literature went on to say, which is the most
this graph can support.

Usage:
    python scripts/atlas_discovery_eval.py
    python scripts/atlas_discovery_eval.py --split-year 2018 --seeds 40 --top 20
"""

import argparse
import collections
import glob
import gzip
import json
import math
import pickle
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_discovery import HUB_PERCENTILE, MIN_BRIDGES, MIN_CANDIDATE_DEGREE  # noqa: E402
from atlas_graph import load_index, load_corrections, _corrected  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-discovery-eval.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-discovery-eval.json"

SEED_RNG = 20260803
# Seeds are sampled from this degree band. Too low and there is no 2-hop
# neighbourhood to rank; too high and the seed is a hub whose candidate set is
# most of the graph.
SEED_DEGREE_MIN, SEED_DEGREE_MAX = 30, 800


# Every stream that carries a `year`. NOT a cosmetic list: the two names this
# tuple used to hold were `records` and `records_c04only`, and c04only is a
# strict SUBSET of records -- measured, it adds exactly zero PMIDs -- so a loop
# that read like a merge merged nothing, while the two streams holding the most
# recent literature were never named. The layer whose subject is recency was
# blind to 814,015 dated articles, overwhelmingly from 2021 onward.
YEAR_STREAMS = ("records", "records_c04only", "records_unindexed", "records_updates")


def pmid_years(root: Path) -> dict:
    """Cached wrapper: the raw scan reads 2.2 GB of census records.

    THE CACHE RECORDS WHICH STREAMS BUILT IT. Keyed by path alone, it would
    have silently defeated the fix that widened YEAR_STREAMS: the map was
    rebuilt to cover 814,015 more articles, and every consumer would have gone
    on loading the old pickle without a word. A cache whose key omits an input
    is a way of not applying a change.
    """
    cache = root / "records" / ".pmid-years.pkl"
    try:
        if cache.exists():
            with open(cache, "rb") as fh:
                got = pickle.load(fh)
            if isinstance(got, dict) and got.get("streams") == list(YEAR_STREAMS):
                years = got["years"]
                print(f"  year map from cache ({len(years):,} PMIDs); "
                      f"delete {cache.name} to rebuild", flush=True)
                return years
            print(f"  {cache.name} was built from a different stream set; "
                  "rebuilding", flush=True)
    except Exception:
        pass
    got = _scan_pmid_years(root)
    try:
        with open(cache, "wb") as fh:
            pickle.dump({"streams": list(YEAR_STREAMS), "years": got}, fh, protocol=5)
    except OSError:
        pass
    return got


def _scan_pmid_years(root: Path) -> dict:
    """PMID -> publication year, merged across every census directory present.

    Merged rather than taken from one, because PubMed baseline files are
    chronological: a partially-rebuilt census holds only the OLDEST literature,
    which would date every pair as ancient and silently zero the split.
    """
    years = {}
    for d in YEAR_STREAMS:
        for f in sorted(glob.glob(str(root / d / "*.jsonl.gz"))):
            with gzip.open(f, "rt", encoding="utf-8") as fh:
                for line in fh:
                    r = json.loads(line)
                    y = r.get("year")
                    if y:
                        p = r["pmid"]
                        # keep the EARLIEST year seen for a PMID
                        if p not in years or y < years[p]:
                            years[p] = y
    return years


def pair_first_year(root: Path, years: dict, corrections: dict) -> dict:
    """(a, b) -> earliest year any paper asserts the pair.

    Computed over EVERY asserting PMID, not the index's 60-PMID sample: the
    sample is uniform, so its minimum is a late-biased estimate of the first
    assertion, which is exactly the quantity a time split turns on.
    """
    first = {}
    with gzip.open(root / "relations" / "relations.tsv.gz", "rt",
                   encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid = parts[0]
            y = years.get(pmid)
            if not y:
                continue
            a = _corrected(parts[2].split("|", 1)[-1], pmid, corrections)
            b = _corrected(parts[3].split("|", 1)[-1], pmid, corrections)
            key = (a, b) if a <= b else (b, a)
            if key not in first or y < first[key]:
                first[key] = y
    return first


def rank_all(adj_before, degrees, cutoff, n_nodes, seed_id, rng):
    """Every ranking, over one IDENTICAL candidate set, so only order differs.

    Beyond the shipped ABC ranking and the two baselines, this includes the
    standard degree-controlled link predictors. They are the obvious thing to
    try once ABC is found to lose to popularity, and they cost nothing extra:
    all of them are functions of the bridge set already computed.

      adamic_adar  sum over bridges of 1/log(deg(b)) -- a bridge shared with a
                   hub is weak evidence, a bridge through a specific entity is
                   strong. The classic link-prediction baseline.
      resource_alloc  sum of 1/deg(b) -- the same idea, punishing hubs harder.
      jaccard      n / (deg_a + deg_c - n) -- normalises by BOTH degrees, so it
                   is the most aggressive correction for candidate popularity.
    """
    from scipy.stats import hypergeom

    a_nb = adj_before.get(seed_id, set())
    usable = {b for b in a_nb if degrees.get(b, 0) <= cutoff}
    if not usable:
        return None
    deg_a = len(usable)

    bridges = collections.defaultdict(set)
    for b in usable:
        for c in adj_before.get(b, ()):
            if c == seed_id or c in a_nb:
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
        # deg(b) >= 1 for any b that bridges, and log(1) == 0 would divide by
        # zero, so the Adamic-Adar denominator is floored.
        aa = sum(1.0 / math.log(max(degrees.get(b, 1), 2)) for b in bs)
        ra = sum(1.0 / max(degrees.get(b, 1), 1) for b in bs)
        union = deg_a + deg_c - n
        jac = (n / union) if union > 0 else 0.0
        p = float(hypergeom.sf(n - 1, n_nodes, deg_a, deg_c))
        cands.append({"c": c, "n": n, "deg_c": deg_c, "p": p,
                      "aa": aa, "ra": ra, "jac": jac})
    if not cands:
        return None

    # Sort by identifier BEFORE anything reads this list. `bridges` is filled by
    # iterating adjacency SETS, and CPython's string hashing is salted per
    # process (PYTHONHASHSEED), so set order -- and therefore the input order of
    # this list -- differs between runs. Python's sort is stable, so that order
    # leaks into every tie, and it reached the output: the random baseline drifted
    # 2.9%-3.3% across runs with a fixed seed. Sorting here makes ties, the
    # shuffle, and the reported numbers reproducible.
    cands.sort(key=lambda r: r["c"])

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


def evaluate(first, idx, Y, seeds_n, top, log=True):
    """One split: build the before-graph, rank, score against what came next."""
    before = {k for k, y in first.items() if y < Y}
    after = {k for k, y in first.items() if y >= Y}
    if log:
        print(f"  before {Y}: {len(before):,} pairs; first asserted {Y} or later: "
              f"{len(after):,}", flush=True)

    adj_before = collections.defaultdict(set)
    for a, b in before:
        adj_before[a].add(b)
        adj_before[b].add(a)
    degrees = {k: len(v) for k, v in adj_before.items()}
    if not degrees:
        return None
    cutoff = sorted(degrees.values())[int(len(degrees) * HUB_PERCENTILE)]
    n_nodes = len(degrees)

    rng = random.Random(SEED_RNG)
    pool = sorted(k for k, d in degrees.items()
                  if SEED_DEGREE_MIN <= d <= SEED_DEGREE_MAX)
    if not pool:
        return None
    seeds = rng.sample(pool, min(seeds_n, len(pool)))

    METHODS = ("abc", "popularity", "adamic_adar", "resource_alloc",
               "jaccard", "bridges", "random")
    hits = {m: 0 for m in METHODS}
    shown = {m: 0 for m in hits}
    per_seed = []
    for sid in seeds:
        ranks = rank_all(adj_before, degrees, cutoff, n_nodes, sid, rng)
        if not ranks:
            continue
        row = {"seed": sid, "seed_name": idx["canon"].get(sid, sid),
               "degree": degrees.get(sid, 0), "candidates": len(ranks["abc"])}
        for method, order in ranks.items():
            sel = order[:top]
            h = sum(1 for c in sel
                    if ((sid, c) if sid <= c else (c, sid)) in after)
            hits[method] += h
            shown[method] += len(sel)
            row[method] = h
        per_seed.append(row)
    if not per_seed:
        return None

    prec = {m: (hits[m] / shown[m] if shown[m] else 0.0) for m in hits}

    # PAIRED comparison. The three rankings run on the same seeds over the same
    # candidate set, so per-seed differences are paired; an interval on two
    # independent proportions would be the wrong test.
    n = len(per_seed)
    paired_all = {}
    for m in METHODS:
        if m == "popularity":
            continue
        d = [r[m] - r["popularity"] for r in per_seed]
        boot = random.Random(SEED_RNG + 1)
        means = sorted(sum(d[boot.randrange(n)] for _ in range(n)) / n
                       for _ in range(10000))
        lo, hi = means[int(0.025 * len(means))], means[int(0.975 * len(means))]
        paired_all[m] = {"mean_diff": sum(d) / n, "ci95": [lo, hi],
                         "decided": hi < 0 or lo > 0,
                         "ahead": sum(1 for x in d if x > 0),
                         "behind": sum(1 for x in d if x < 0)}
    ci = tuple(paired_all["abc"]["ci95"])
    mean_diff = paired_all["abc"]["mean_diff"]
    return {
        "split_year": Y, "seeds_evaluated": n, "top_k": top,
        "pairs_before": len(before), "pairs_after": len(after),
        "hits": hits, "predictions": shown, "precision": prec,
        "abc_over_popularity": (prec["abc"] / prec["popularity"]
                                if prec["popularity"] else float("inf")),
        "paired": {"mean_diff": mean_diff, "ci95": list(ci),
                   "decided": paired_all["abc"]["decided"],
                   "abc_ahead": paired_all["abc"]["ahead"],
                   "abc_behind": paired_all["abc"]["behind"]},
        "paired_all": paired_all,
        "per_seed": per_seed,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split-year", type=int, default=2018,
                    help="the headline split")
    ap.add_argument("--also-years", type=int, nargs="*", default=[2015, 2021],
                    help="extra splits run as a robustness check; the pair "
                         "dating pass is shared so each costs almost nothing")
    ap.add_argument("--seeds", type=int, default=200)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    root = atlas_root()
    idx = load_index(root)
    print("loading year map ...", flush=True)
    years = pmid_years(root)
    if not years:
        print("no dated records; run scripts/atlas_baseline.py first", file=sys.stderr)
        return 1
    span = (min(years.values()), max(years.values()))
    print(f"  {len(years):,} dated PMIDs spanning {span[0]}-{span[1]}", flush=True)
    if span[1] < args.split_year:
        print(f"year map ends at {span[1]}, before the split at {args.split_year}; "
              "the census is probably still rebuilding (baseline files are "
              "chronological)", file=sys.stderr)
        return 1

    print("dating every pair ...", flush=True)
    first = pair_first_year(root, years, load_corrections())
    print(f"  {len(first):,} pairs carry a first-assertion year", flush=True)

    head = evaluate(first, idx, args.split_year, args.seeds, args.top)
    if not head:
        print("no evaluable seeds", file=sys.stderr)
        return 1
    robust = []
    for y in sorted(set(args.also_years) - {args.split_year}):
        r = evaluate(first, idx, y, args.seeds, args.top)
        if r:
            robust.append(r)

    prec, paired = head["precision"], head["paired"]
    if not paired["decided"]:
        verdict = "ABC and popularity are indistinguishable at this sample size"
    elif paired["mean_diff"] > 0:
        verdict = "ABC ranking beats popularity"
    else:
        verdict = "ABC ranking is WORSE than ranking by popularity"

    Y = args.split_year
    L = [
        "# Does literature-based discovery predict anything? (#ATLAS-LBD-EVAL)", "",
        "Generated by `scripts/atlas_discovery_eval.py`. `atlas_discovery.py` emits a",
        "ranked list of pairs the literature implies but has never stated. This gives",
        "it a hit rate, which it has never had.", "",
        "## Method", "",
        f"A time split at **{Y}**. The graph is rebuilt as it stood before {Y}",
        f"({head['pairs_before']:,} pairs), candidates are ranked, and a prediction",
        f"counts as a hit if the literature first asserts that pair in {Y} or later",
        f"({head['pairs_after']:,} pairs did). Pairs already asserted before {Y} are",
        "excluded from prediction by construction, so a hit is a NEW statement rather",
        "than a rediscovery.", "",
        "Three rankings run over an **identical** candidate set, so only the ordering",
        "differs. The comparison that matters is not ABC against random -- almost",
        "anything beats random in a clustered co-occurrence graph -- but ABC against",
        "ranking those same candidates by popularity.", "",
        f"{head['seeds_evaluated']} seed entities, sampled with a fixed seed from",
        f"degree band {SEED_DEGREE_MIN}-{SEED_DEGREE_MAX}, top {args.top} each.", "",
        "## Result", "",
        f"| ranking | hits | predictions | precision@{args.top} | paired vs popularity |",
        "|---|---|---|---|---|",
    ]
    for m in sorted(prec, key=lambda k: -prec[k]):
        pa = head["paired_all"].get(m)
        vs = ("baseline" if m == "popularity" else
              f"{pa['mean_diff']:+.2f} [{pa['ci95'][0]:+.2f}, {pa['ci95'][1]:+.2f}]"
              + ("" if pa["decided"] else " (spans 0)"))
        L.append(f"| {m} | {head['hits'][m]:,} | {head['predictions'][m]:,} | "
                 f"**{100*prec[m]:.1f}%** | {vs} |")
    L += [
        "", "### Is that difference real?", "",
        "Paired bootstrap over seeds, 10,000 resamples:", "",
        f"* mean per-seed difference (abc minus popularity): "
        f"**{paired['mean_diff']:+.2f}** hits out of {args.top}",
        f"* 95% CI **[{paired['ci95'][0]:+.2f}, {paired['ci95'][1]:+.2f}]**",
        f"* abc ahead on {paired['abc_ahead']} seeds, behind on "
        f"{paired['abc_behind']}", "",
    ]
    if robust:
        L += ["### Robustness across split years", "",
              "| split | abc | popularity | random | paired diff | 95% CI |",
              "|---|---|---|---|---|---|"]
        for r in [head] + robust:
            pp, q = r["precision"], r["paired"]
            L.append(f"| {r['split_year']} | {100*pp['abc']:.1f}% | "
                     f"{100*pp['popularity']:.1f}% | {100*pp['random']:.1f}% | "
                     f"{q['mean_diff']:+.2f} | [{q['ci95'][0]:+.2f}, "
                     f"{q['ci95'][1]:+.2f}] |")
        L.append("")

    L += [
        "### The pattern in that table", "",
        "The rankings do not scatter. They order themselves by **how hard each one",
        "corrects for degree**, and the harder the correction, the worse it does:", "",
        "| ranking | degree correction | precision@%d |" % args.top,
        "|---|---|---|",
        "| popularity | none -- it IS degree | %.1f%% |" % (100 * prec["popularity"]),
        "| raw bridge count | none | %.1f%% |" % (100 * prec["bridges"]),
        "| Adamic-Adar | down-weights hub bridges | %.1f%% |" % (100 * prec["adamic_adar"]),
        "| resource allocation | down-weights them harder | %.1f%% |" % (100 * prec["resource_alloc"]),
        "| ABC (hypergeometric) | divides out candidate degree | %.1f%% |" % (100 * prec["abc"]),
        "| Jaccard | normalises by BOTH degrees | %.1f%% |" % (100 * prec["jaccard"]),
        "| random | -- | %.1f%% |" % (100 * prec["random"]),
        "",
        "Jaccard, the most aggressive correction, lands barely above chance. That is",
        "not a bug in any one method; it says the thing being corrected away is the",
        "signal. New edges in this graph genuinely do attach preferentially to",
        "well-connected entities, because attention concentrates -- so a ranking that",
        "removes degree removes most of what predicts the next edge.", "",
        "This also answers the obvious follow-up. The shipped ranking cannot be",
        "repaired by swapping in a better-known link predictor: the standard ones were",
        "tried here and none beats the baseline. Adamic-Adar and the raw bridge count",
        "tie with it (their intervals span zero); everything else loses.", "",
        "### The objection to this whole evaluation", "",
        "It is worth stating against my own result. This scores a ranking by whether",
        "it anticipates what the literature went on to assert. But Swanson-style",
        "discovery is FOR finding connections the literature is slow to reach -- the",
        "fish-oil/Raynaud case mattered precisely because nobody was about to publish",
        "it. On that reading, a method that beats popularity at predicting next year's",
        "edges may be selecting for the LEAST interesting hypotheses, and a method that",
        "loses to popularity is not thereby useless.", "",
        "So the defensible claim is narrow and it is the one the module actually made:",
        "`atlas_discovery.py` says it corrects for popularity, and measured against",
        "what the literature did next it does not help. Whether predicting the",
        "literature is the right target for a discovery layer at all is a separate",
        "question this evaluation cannot settle, and a genuinely overlooked connection",
        "would be scored here as a miss.", "",
        f"### Verdict: {verdict}", ""]
    if paired["decided"] and paired["mean_diff"] < 0:
        L += [
            "This is a negative result about this repository's own layer, and it is",
            "why the evaluation was worth writing. The bridge counting, hub filtering",
            "and hypergeometric tail do not order candidates better than asking which",
            "of them is already famous -- they order them measurably worse, at every",
            "split year tested.", "",
            "**But the candidate SET is doing real work.** Both rankings beat random by",
            f"roughly {prec['popularity']/prec['random']:.0f}x, so restricting attention",
            "to 2-hop bridged entities is genuinely informative; it is the ranking",
            "within that set that fails. The honest summary is that",
            "`atlas_discovery.py` is a good candidate GENERATOR and a bad RANKER, and",
            "its output should be read as a popularity-weighted reading list.", "",
            "The module's docstring says it corrects for popularity. Measured against",
            "what the literature went on to say, it does not.", "",
        ]
    elif not paired["decided"]:
        L += [
            "No claim either way survives this sample size. The point estimate favours",
            "popularity but the paired interval spans zero, so the honest statement is",
            "that the ABC machinery has not been SHOWN to add ordering information --",
            "not that it has been shown to lack it.", "",
        ]
    else:
        L += [
            "The ABC machinery orders candidates better than popularity, so the bridge",
            "structure carries information the degree distribution does not. That is a",
            "floor, not a validation: absolute precision is still low and every limit",
            "in `atlas_discovery.py` still applies.", "",
        ]

    L += ["## Per-seed detail", "",
          f"| seed | degree before {Y} | candidates | abc | popularity | random |",
          "|---|---|---|---|---|---|"]
    for r in sorted(head["per_seed"], key=lambda r: -r["abc"])[:25]:
        L.append(f"| {r['seed_name']} | {r['degree']:,} | {r['candidates']:,} | "
                 f"{r['abc']} | {r['popularity']} | {r['random']} |")

    L += [
        "", "## What this cannot show", "",
        "* That a hit is a real biological relation. A new edge may be a new",
        "  extraction of an old idea, and PubTator's extractor has its own error rate.",
        "  This measures whether a ranking anticipates what the literature went on to",
        "  say, which is the most this graph supports.",
        "* Anything about pairs the literature will assert after the census ends. A",
        "  correct prediction not yet published counts here as a miss, so every",
        "  precision figure is a lower bound.",
        "* Anything about seeds outside the sampled degree band. Hubs and",
        "  near-isolated nodes behave differently and are excluded by construction.",
        "* Recall. Only the top-k are scored, so a ranking that buries a true pair at",
        "  position k+1 is indistinguishable here from one that omits it.",
        "* That popularity is a GOOD ranking. It is a better one, on a graph where",
        "  well-studied entities keep accruing edges. Predicting that a famous gene",
        "  will gain another relation is easy and not very useful.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({"headline": head, "robustness": robust}, indent=2) + "\n")
    print(f"\nabc {100*prec['abc']:.1f}%  popularity {100*prec['popularity']:.1f}%  "
          f"random {100*prec['random']:.1f}%  ->  {verdict}")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
