#!/usr/bin/env python3
"""Atlas: does literature-based discovery predict anything? (#ATLAS-LBD-EVAL)

WHY
---
`atlas_discovery.py` emits a ranked list of A-C pairs absent from its relation
graph. It calls the output "a ranked reading list, nothing here is a finding".
This evaluation measures a hit rate against later dated, observed assertions.

Without one there is no way to tell the layer apart from a list of famous
entities. Run on GPX4 it returns ERK, caspase-3, cyclin D1, MMP-9 and ATP --
exactly what a popularity ranking would return, which is the failure mode the
module's own docstring says it corrects for.

Use a TIME SPLIT: reconstruct the graph from dated assertions before year Y,
predict which absent A-C pairs will appear, and check the remaining observations.

THE COMPARISON THAT MATTERS
---------------------------
The question is whether ABC beats **ranking the same candidates by popularity**.
The random baseline also samples that same candidate pool. Beating it measures
ordering within the pool; it cannot establish the value of candidate generation.

Seven rankings use an IDENTICAL candidate set, including:

  * `abc`        -- the shipped ranking (hypergeometric tail over bridge counts)
  * `popularity` -- rank by candidate degree in the before-graph
  * `random`     -- a seeded shuffle within the pool

WHAT COUNTS AS A HIT
--------------------
A predicted pair A-C is a hit if its earliest dated, observed assertion is in
year >= Y. Undated or missing assertions can conceal prior knowledge, so this
does not establish novelty.

WHAT THIS CANNOT SHOW
---------------------
That a hit is a real biological relation. PubTator's extractor has its own error
rate, and a new edge may be a new extraction of an old idea. It measures whether
the ranking anticipates what the literature went on to say, which is the most
this graph can support.

Usage:
    python scripts/atlas_discovery_eval.py
    python scripts/atlas_discovery_eval.py --split-year 2018 --seeds 40 --top 20
    python scripts/atlas_discovery_eval.py --render-only
"""

import argparse
import collections
import functools
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

OUT_MD = PROJECT_ROOT / "analysis" / "atlas-discovery-eval.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-discovery-eval.json"
OUT, RAW = OUT_MD, OUT_JSON  # Historical import names.

SEED_RNG = 20260803
METHODS = ("abc", "popularity", "adamic_adar", "resource_alloc",
           "jaccard", "bridges", "random")
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



def _integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


@functools.lru_cache(maxsize=128)
def _paired_summary(differences):
    """The original bootstrap, cached only by its immutable ordered input.

    Re-rendering must not change the RNG draws, seed order, or interval method.
    Returning an immutable tuple keeps callers from mutating the cached result.
    """
    n = len(differences)
    boot = random.Random(SEED_RNG + 1)
    means = sorted(sum(differences[boot.randrange(n)] for _ in range(n)) / n
                   for _ in range(10000))
    lo, hi = means[int(0.025 * len(means))], means[int(0.975 * len(means))]
    return (sum(differences) / n, lo, hi, hi < 0 or lo > 0,
            sum(1 for x in differences if x > 0),
            sum(1 for x in differences if x < 0))


def _assemble_split(raw):
    if not isinstance(raw, dict):
        raise ValueError("each split must be an object")
    metadata = {}
    for key, minimum in (("split_year", 1), ("top_k", 1),
                         ("pairs_before", 0), ("pairs_after", 0)):
        metadata[key] = _integer(raw.get(key), key, minimum)
    rows = raw.get("per_seed")
    if not isinstance(rows, list) or not rows:
        raise ValueError("per_seed must contain at least one evaluable seed")
    per_seed, seen = [], set()
    top = metadata["top_k"]
    for position, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"per_seed[{position}] must be an object")
        sid = row.get("seed")
        if not isinstance(sid, str) or not sid.strip():
            raise ValueError(f"per_seed[{position}].seed must be a nonempty string")
        if sid in seen:
            raise ValueError(f"duplicate seed {sid!r} in split {metadata['split_year']}")
        seen.add(sid)
        name = row.get("seed_name")
        if not isinstance(name, str):
            raise ValueError(f"seed_name for {sid!r} must be a string")
        degree = _integer(row.get("degree"), f"degree for {sid!r}")
        candidates = _integer(row.get("candidates"), f"candidates for {sid!r}", 1)
        copied = {"seed": sid, "seed_name": name, "degree": degree,
                  "candidates": candidates}
        selected = min(top, candidates)
        for method in METHODS:
            count = _integer(row.get(method), f"{method} hits for {sid!r}")
            if count > selected:
                raise ValueError(f"{method} hits for {sid!r} exceed {selected} predictions")
            copied[method] = count
        per_seed.append(copied)

    hits = {m: sum(r[m] for r in per_seed) for m in METHODS}
    denominator = sum(min(top, r["candidates"]) for r in per_seed)
    shown = {m: denominator for m in METHODS}
    precision = {m: hits[m] / denominator for m in METHODS}
    paired_all = {}
    for method in METHODS:
        if method == "popularity":
            continue
        differences = tuple(r[method] - r["popularity"] for r in per_seed)
        mean, lo, hi, decided, ahead, behind = _paired_summary(differences)
        paired_all[method] = {"mean_diff": mean, "ci95": [lo, hi],
                              "decided": decided, "ahead": ahead, "behind": behind}
    abc = paired_all["abc"]
    return {
        "split_year": metadata["split_year"], "seeds_evaluated": len(per_seed),
        "top_k": top, "pairs_before": metadata["pairs_before"],
        "pairs_after": metadata["pairs_after"],
        "hits": hits, "predictions": shown, "precision": precision,
        "abc_over_popularity": (precision["abc"] / precision["popularity"]
                                if precision["popularity"] else None),
        "paired": {"mean_diff": abc["mean_diff"], "ci95": list(abc["ci95"]),
                   "decided": abc["decided"], "abc_ahead": abc["ahead"],
                   "abc_behind": abc["behind"]},
        "paired_all": paired_all, "per_seed": per_seed,
    }


def assemble(raw):
    """Validate observations and recompute summaries without mutating the input.

    Stored hits, precision, ratios, verdicts, and intervals are derived fields;
    none is trusted as evidence. The historical per-seed order is retained.
    """
    if not isinstance(raw, dict):
        raise ValueError("evaluation must be an object")
    robustness = raw.get("robustness", [])
    if not isinstance(robustness, list):
        raise ValueError("robustness must be a list")
    splits = [_assemble_split(raw.get("headline"))]
    splits.extend(_assemble_split(split) for split in robustness)
    seen = set()
    for split in splits:
        year = split["split_year"]
        if year in seen:
            raise ValueError(f"duplicate split year {year}")
        seen.add(year)
    return {"headline": splits[0], "robustness": splits[1:]}


def evaluate(first, idx, Y, seeds_n, top, log=True):
    """One split: build the before-graph, rank, score against later observations."""
    before = {k for k, y in first.items() if y < Y}
    after = {k for k, y in first.items() if y >= Y}
    if log:
        print(f"  before {Y}: {len(before):,} pairs; earliest dated assertion "
              f"{Y} or later: {len(after):,}", flush=True)

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

    per_seed = []
    for sid in seeds:
        ranks = rank_all(adj_before, degrees, cutoff, n_nodes, sid, rng)
        if not ranks:
            continue
        row = {"seed": sid, "seed_name": idx["canon"].get(sid, sid),
               "degree": degrees.get(sid, 0), "candidates": len(ranks["abc"])}
        for method, order in ranks.items():
            sel = order[:top]
            row[method] = sum(1 for c in sel
                              if ((sid, c) if sid <= c else (c, sid)) in after)
        per_seed.append(row)
    if not per_seed:
        return None
    return _assemble_split({
        "split_year": Y, "top_k": top, "pairs_before": len(before),
        "pairs_after": len(after), "per_seed": per_seed,
    })


def _direction(paired):
    lo, hi = paired["ci95"]
    if lo > 0:
        return "higher than popularity"
    if hi < 0:
        return "lower than popularity"
    return "interval includes zero"


def _verdict(paired):
    direction = _direction(paired)
    if direction == "higher than popularity":
        return "ABC ranking beats popularity on observed future assertions"
    if direction == "lower than popularity":
        return "ABC ranking performs worse than popularity on observed future assertions"
    return "The ABC versus popularity comparison is inconclusive"


def _cell(value):
    """Keep a canonical entity name within its Markdown table cell."""
    return " ".join(str(value).split()).replace("|", "\\|")


def render(raw):
    """Render only recomputed results; this function performs no file or graph I/O."""
    data = assemble(raw)
    head, robust = data["headline"], data["robustness"]
    year, top = head["split_year"], head["top_k"]
    prec, paired = head["precision"], head["paired"]
    short = sum(r["candidates"] < top for r in head["per_seed"])
    lines = [
        "# Does literature-based discovery predict anything? (#ATLAS-LBD-EVAL)", "",
        "Generated by `scripts/atlas_discovery_eval.py`. This evaluates how rankings",
        "anticipate later dated, observed assertions within a shared candidate pool.", "",
        "## Method", "",
        f"A time split at **{year}** uses {head['pairs_before']:,} pairs whose earliest",
        f"dated, observed assertion is before {year}. A prediction is a hit when the",
        f"pair's earliest dated, observed assertion is in {year} or later",
        f"({head['pairs_after']:,} pairs across the graph). Pairs observed in dated",
        "pre-split records are excluded from prediction. Missing or undated earlier",
        "assertions can conceal prior knowledge, so a hit does not establish novelty.", "",
        "All seven rankings use an **identical candidate set** for each seed.",
        "Candidates are eligible two-hop neighbors in the pre-split graph; the",
        "existing bridge, hub, and candidate-degree filters apply. Popularity ranks",
        "by pre-split candidate degree. Random is one seeded shuffle **within the",
        "same candidate pool**, not a sample from unrestricted graph pairs.", "",
        f"The report includes {head['seeds_evaluated']} evaluable seed entities from",
        f"the fixed-seed sample in degree band {SEED_DEGREE_MIN}-{SEED_DEGREE_MAX}.",
        f"Each method selects up to {top} candidates per seed: min(k, candidate count).",
        f"{short} seeds have fewer than {top} candidates. Precision is total hits",
        "divided by the actual number of predictions, pooled over those seeds.", "",
        "## Result", "",
        f"| ranking | hits | predictions | precision@{top} | paired vs popularity |",
        "|---|---|---|---|---|",
    ]
    for method in sorted(prec, key=lambda m: -prec[m]):
        comparison = "baseline"
        if method != "popularity":
            pa = head["paired_all"][method]
            lo, hi = pa["ci95"]
            comparison = (f"{pa['mean_diff']:+.2f} [{lo:+.2f}, {hi:+.2f}]; "
                          f"{_direction(pa)}")
        lines.append(f"| {method} | {head['hits'][method]:,} | "
                     f"{head['predictions'][method]:,} | **{100*prec[method]:.1f}%** | "
                     f"{comparison} |")

    lines.extend([
        "", "### Paired uncertainty", "",
        "Paired bootstrap over seeds, 10,000 resamples. Differences are hits per",
        f"seed, with at most {top} predictions per seed; they are not percentage points.",
        f"The mean difference (ABC minus popularity) is **{paired['mean_diff']:+.2f}**",
        f"with 95% percentile interval **[{paired['ci95'][0]:+.2f}, "
        f"{paired['ci95'][1]:+.2f}]**. ABC is ahead on {paired['abc_ahead']} seeds",
        f"and behind on {paired['abc_behind']}.", "",
        "Intervals are unadjusted comparisons against popularity, not simultaneous",
        "95% guarantees for all methods. An interval including zero leaves the",
        "difference unresolved; it does not establish that two rankings are equal.", "",
        f"### Verdict: {_verdict(paired)}", "",
    ])
    if not paired["decided"]:
        mean = paired["mean_diff"]
        if mean > 0:
            lines.append("The point estimate favors ABC, but its interval includes zero.")
        elif mean < 0:
            lines.append("The point estimate favors popularity, but its interval includes zero.")
        else:
            lines.append("The point estimates are equal, and the interval includes zero.")
        lines.append("")
    lines.extend([
        "These comparisons describe ordering within the evaluated candidate pools.",
        "Outperforming the within-pool random baseline does **not** establish that",
        "the candidate generator improves on unrestricted candidate selection.",
        "No such control is evaluated here. Nor does this comparison establish",
        "biological validity or the usefulness of an overlooked hypothesis.", "",
    ])
    if robust:
        lines.extend([
            "### Comparisons across split years", "",
            "| split | k | abc | popularity | random | paired diff | 95% CI | conclusion |",
            "|---|---|---|---|---|---|---|---|",
        ])
        for split in [head] + robust:
            pp, pa = split["precision"], split["paired"]
            lo, hi = pa["ci95"]
            lines.append(f"| {split['split_year']} | {split['top_k']} | "
                         f"{100*pp['abc']:.1f}% | {100*pp['popularity']:.1f}% | "
                         f"{100*pp['random']:.1f}% | {pa['mean_diff']:+.2f} | "
                         f"[{lo:+.2f}, {hi:+.2f}] | {_direction(pa)} |")
        directions = {_direction(split["paired"]) for split in [head] + robust}
        if len(directions) == 1 and paired["decided"]:
            lines.extend(["", f"ABC is {_direction(paired)} in every evaluated split."])
        else:
            lines.extend(["", "The evaluated splits do not all resolve the ABC comparison in the same direction."])
        lines.extend([
            "Each split follows assertions through the available corpus snapshot;",
            "follow-up lengths and sampled seeds can differ. This is not a comparison",
            "at a common prediction horizon. Overlapping years are not independent replications.", "",
        ])

    lines.extend([
        "## Per-seed detail", "",
        "Up to 25 seeds with the most ABC hits are shown; aggregate metrics use every evaluable seed.", "",
        f"| seed | degree before {year} | candidates | abc | popularity | random |",
        "|---|---|---|---|---|---|",
    ])
    for row in sorted(head["per_seed"], key=lambda r: -r["abc"])[:25]:
        lines.append(f"| {_cell(row['seed_name'])} | {row['degree']:,} | "
                     f"{row['candidates']:,} | {row['abc']} | "
                     f"{row['popularity']} | {row['random']} |")
    lines.extend([
        "", "## What this cannot show", "",
        "* Biological validity or genuinely new knowledge. The extractor can miss or",
        "  misidentify assertions, and absent dated evidence is not absence of prior knowledge.",
        "* Eventual precision beyond the available corpus. Unobserved future assertions",
        "  count as misses here. The stored results do not record an observation end year.",
        "* Uncertainty for a new corpus or extraction process. The bootstrap resamples",
        "  seeds and assumes exchangeability; overlapping entities, papers, and pairs",
        "  can violate independence. Its intervals are conditional on this corpus and",
        "  the observed per-seed outcomes, not calibrated biological confidence bounds.",
        "* Performance outside the sampled degree band or for seeds without eligible",
        "  candidates. Those seeds do not contribute to the reported denominator.",
        "* Recall or candidate-generator benefit. Only the selected top-k hits are",
        "  evaluated, and all ranking baselines share the same candidate pool.",
        "* Whether popularity reflects biology, publication attention, or both. This",
        "  design does not separate those explanations or validate discovery utility.", "",
        "An objection to the target itself: discovery can aim at overlooked",
        "connections that the literature is slow to reach. Those hypotheses count",
        "as misses here until a dated assertion is observed, so predicting later",
        "assertions and identifying useful discoveries are different objectives.", "",
    ])
    return "\n".join(lines)


def _positive_int(value):
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError("must be a positive integer") from None
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split-year", type=_positive_int, default=2018,
                    help="the headline split")
    ap.add_argument("--also-years", type=_positive_int, nargs="*", default=[2015, 2021],
                    help="additional requested splits; all must be evaluable")
    ap.add_argument("--seeds", type=_positive_int, default=200)
    ap.add_argument("--top", type=_positive_int, default=20)
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild Markdown from stored per-seed counts without graph I/O")
    args = ap.parse_args(argv)

    try:
        if args.render_only:
            raw = json.loads(OUT_JSON.read_text(encoding="utf-8"))
            data = assemble(raw)
        else:
            root = atlas_root()
            idx = load_index(root)
            print("loading year map ...", flush=True)
            years = pmid_years(root)
            if not years:
                raise ValueError("no dated records; run scripts/atlas_baseline.py first "
                                 "or use --render-only for stored results")
            for value in years.values():
                _integer(value, "publication year", 1)
            span = (min(years.values()), max(years.values()))
            print(f"  {len(years):,} dated PMIDs spanning {span[0]}-{span[1]}", flush=True)
            requested = [args.split_year] + sorted(set(args.also_years) - {args.split_year})
            unavailable = [year for year in requested if not span[0] < year <= span[1]]
            if unavailable:
                raise ValueError(f"requested split years {unavailable} unavailable in dated "
                                 f"record span {span[0]}-{span[1]}; each needs observations "
                                 "before and at or after the split")
            print("dating every pair ...", flush=True)
            first = pair_first_year(root, years, load_corrections())
            print(f"  {len(first):,} pairs carry an earliest dated assertion year", flush=True)
            results = []
            for year in requested:
                result = evaluate(first, idx, year, args.seeds, args.top)
                if not result:
                    raise ValueError(f"no evaluable seeds for requested split {year}; "
                                     "reports were not updated")
                results.append(result)
            raw = {"headline": results[0], "robustness": results[1:]}
            data = assemble(raw)
    except (OSError, ValueError, TypeError) as exc:
        print(f"discovery evaluation unavailable: {exc}", file=sys.stderr)
        return 1

    # Both payloads must be valid before either artifact is replaced. Replay
    # deliberately preserves the original JSON bytes and historical observations.
    json_text = json.dumps(data, indent=2, allow_nan=False) + "\n"
    markdown = render(data)
    if not args.render_only:
        OUT_JSON.write_text(json_text, encoding="utf-8")
    OUT_MD.write_text(markdown, encoding="utf-8")
    head = data["headline"]
    print(f"abc {100*head['precision']['abc']:.1f}%  "
          f"popularity {100*head['precision']['popularity']:.1f}%  "
          f"random {100*head['precision']['random']:.1f}%  -> {_verdict(head['paired'])}")
    print(f"wrote {OUT_MD}" + ("" if args.render_only else f"\nwrote {OUT_JSON}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
