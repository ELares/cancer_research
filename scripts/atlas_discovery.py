#!/usr/bin/env python3
"""Atlas: literature-based discovery over the relation graph (#ATLAS).

WHY
---
This is the analysis the census was built for. Swanson's ABC model: if A relates
to B, and B relates to C, but A and C have never been discussed together, the
A-C link is a candidate hypothesis that exists implicitly in the literature and
that nobody has stated. Swanson's own worked cases -- fish oil and Raynaud's
syndrome (1986), magnesium and migraine (1988) -- were later supported
experimentally.

A 4,830-article corpus cannot do this. Discovery needs both literatures present
at once, and the whole point is that the two literatures do not cite each other.

WHAT IT DOES
------------
For a seed entity A:
  1. collect B = everything A is directly related to;
  2. collect C = everything those Bs relate to, EXCLUDING anything already
     related to A;
  3. rank each C by how many distinct Bs bridge it to A.

THE HARD PART: NOT REDISCOVERING THAT TP53 IS FAMOUS
----------------------------------------------------
Raw bridge counts rank by popularity. TP53 connects to everything, so it bridges
to everything, and a naive ABC search returns a list of the most-studied entities
in oncology every time.

Two corrections are applied:

  * a HYPERGEOMETRIC TAIL TEST, not a ratio. The first attempt scored
    observed/expected under a degree-preserving null, and it was degenerate: for
    an entity of degree 2 the expectation is near zero, so every obscure node
    saturated at the same maximal "lift" and the ranking filled with
    degree-2 variant identifiers. That is the opposite failure from ranking by
    popularity and no more useful. The tail probability of sharing at least the
    observed number of neighbours, under Hypergeometric(N nodes, K=deg(A),
    n=deg(C)), handles both ends correctly: sharing 2 of 2 is unremarkable,
    sharing 40 of 200 is not.
  * a HUB FILTER on the bridges themselves. A B connected to tens of thousands
    of entities carries no information, so bridges above a degree percentile are
    discarded before counting.
  * a MINIMUM DEGREE for candidates, below which nothing is testable.

LIMITS -- READ BEFORE BELIEVING ANY ROW
---------------------------------------
Every output is a HYPOTHESIS, and the base rate for these is poor. Specifically:

  * ABSENCE OF AN A-C EDGE IS NOT ABSENCE OF KNOWLEDGE. It may be missing
    because the extractor failed, because the pair is discussed only in
    full text, or because the relation is stated in a non-cancer paper outside
    this census.
  * The graph does not record direction of effect, so an A-B-C chain may compose
    two relations whose signs cancel. "A relates to B relates to C" does not
    imply A affects C, or in which direction.
  * No context is carried, so a chain can be true in one tissue and meaningless
    in another.
  * Many high-scoring pairs will be trivially related (synonyms, members of the
    same complex, a drug and its own target class) rather than novel.

Treat the output as a ranked reading list. Nothing here is a finding.

Usage:
    python scripts/atlas_discovery.py --seed GPX4
    python scripts/atlas_discovery.py --seed 2879 --top 30
"""

import argparse
import collections
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_graph import load_index, resolve  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT_DIR = PROJECT_ROOT / "analysis"

# Bridges above this degree percentile carry no information and are dropped.
HUB_PERCENTILE = 0.995
MIN_BRIDGES = 3
# Below this, a candidate has too few partners for the tail test to say anything:
# sharing 2 of 2 neighbours is not evidence.
MIN_CANDIDATE_DEGREE = 15


def _comention(a: str, b: str) -> int:
    """PubMed abstracts mentioning BOTH terms. -1 if the query failed.

    This is the check that decides whether a candidate is a discovery at all.
    The ABC premise is that A and C have never been discussed together; if they
    co-occur in PubMed, the missing A-C edge is an extraction failure, not an
    undiscovered link.
    """
    if len(b) > 40 or b.startswith(("RS#", "HGVS", "CorrespondingGene")):
        return -1
    term = f'{a}[tiab] AND "{b}"[tiab]'
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           "?db=pubmed&retmode=json&rettype=count&term=" + urllib.parse.quote(term))
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            n = int(json.load(r)["esearchresult"]["count"])
    except Exception:
        return -1
    time.sleep(0.35)
    return n


def adjacency(idx: dict):
    adj = collections.defaultdict(set)
    for a, b in idx["edges"]:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def discover(idx: dict, adj: dict, seed_id: str, top: int = 30):
    degrees = {k: len(v) for k, v in adj.items()}
    total_edges = len(idx["edges"])
    cutoff = sorted(degrees.values())[int(len(degrees) * HUB_PERCENTILE)]

    a_nb = adj.get(seed_id, set())
    if not a_nb:
        return None, []

    # The null must use the SAME neighbour set the observed overlap is counted
    # over. Counting bridges only over non-hub neighbours while passing the
    # seed's FULL degree as K inflates the expectation and understates every
    # enrichment by the hub fraction.
    usable_nb = {b for b in a_nb if degrees.get(b, 0) <= cutoff}
    deg_a = len(usable_nb)

    bridges = collections.defaultdict(set)
    for b in usable_nb:
        for c in adj.get(b, ()):
            if c == seed_id or c in a_nb:   # already known, not a discovery
                continue
            bridges[c].add(b)

    n_nodes = len(degrees)
    rows = []
    for c, bs in bridges.items():
        n = len(bs)
        if n < MIN_BRIDGES:
            continue
        deg_c = degrees.get(c, 0)
        if deg_c < MIN_CANDIDATE_DEGREE:
            continue
        # P(share >= n) under Hypergeometric(N=n_nodes, K=deg_a, n=deg_c)
        p = float(hypergeom.sf(n - 1, n_nodes, deg_a, deg_c))
        expected = deg_a * deg_c / n_nodes
        rows.append(dict(
            c=c, c_name=idx["canon"].get(c, c), bridges=n, deg_c=deg_c,
            expected=expected, p=p,
            enrichment=(n / expected) if expected > 0 else float("inf"),
            via=[idx["canon"].get(b, b) for b in sorted(bs, key=lambda x: degrees.get(x, 0))][:6]))

    # Benjamini-Hochberg, because this is tens of thousands of simultaneous tests
    rows.sort(key=lambda r: r["p"])
    m = len(rows)
    for i, r in enumerate(rows, 1):
        r["q"] = min(1.0, r["p"] * m / i)
    for i in range(m - 2, -1, -1):
        rows[i]["q"] = min(rows[i]["q"], rows[i + 1]["q"])
    meta = dict(seed=seed_id, seed_name=idx["canon"].get(seed_id, seed_id),
                degree=deg_a, hub_cutoff=cutoff, candidates=len(rows))
    return meta, rows[:top]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", required=True, help="entity symbol or identifier")
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--no-validate", action="store_true",
                    help="skip the PubMed co-mention check (NOT recommended)")
    ap.add_argument("--validate-n", type=int, default=15)
    args = ap.parse_args()

    idx = load_index(atlas_root())
    seed_id = resolve(idx, args.seed)
    if not seed_id:
        raise SystemExit(f"could not resolve {args.seed!r}")

    print("building adjacency ...", flush=True)
    adj = adjacency(idx)
    meta, rows = discover(idx, adj, seed_id, args.top)
    if meta is None:
        raise SystemExit(f"{args.seed} has no relations in the graph")

    if not args.no_validate:
        print(f"validating the top {args.validate_n} against PubMed co-mention ...", flush=True)
        checked = documented = 0
        for r in rows[:args.validate_n]:
            n = _comention(meta["seed_name"], r["c_name"])
            r["comentions"] = n
            if n >= 0:
                checked += 1
                documented += (n > 0)
        meta["validated"] = checked
        meta["documented"] = documented
        meta["documented_share"] = documented / checked if checked else None
        print(f"  {documented}/{checked} candidates are ALREADY co-mentioned in PubMed")

    slug = meta["seed_name"].lower().replace(" ", "-").replace("/", "-")
    out = OUT_DIR / f"atlas-discovery-{slug}.md"
    raw = OUT_DIR / f"atlas-discovery-{slug}.json"
    raw.write_text(json.dumps(dict(meta=meta, rows=rows), indent=1), encoding="utf-8")

    L = [
        f"# Literature-based discovery from {meta['seed_name']} (#ATLAS)", "",
        "Generated by `scripts/atlas_discovery.py`.", "",
        "Swanson's ABC model: A relates to B, B relates to C, but A and C have never",
        "been discussed together, so the A-C link exists implicitly in the literature",
        "and nobody has stated it. Swanson's own cases -- fish oil and Raynaud's (1986),",
        "magnesium and migraine (1988) -- were later supported experimentally.", "",
        f"Seed **{meta['seed_name']}** (`{meta['seed']}`) has {meta['degree']:,} direct partners.",
        f"{meta['candidates']:,} indirectly-linked entities were found with at least",
        f"{MIN_BRIDGES} bridges, after discarding bridge nodes above degree",
        f"{meta['hub_cutoff']:,} (the {HUB_PERCENTILE:.1%} percentile) as uninformative.", "",
        "## Ranking", "",
        "Raw bridge counts rank by popularity -- TP53 bridges to everything. A ratio",
        "score fails at the other end: for a degree-2 entity the expectation is near",
        "zero, so every obscurity saturates at the same maximal value. The ranking is",
        "therefore the **hypergeometric tail probability** of sharing at least the",
        "observed number of neighbours given both degrees, Benjamini-Hochberg corrected",
        "across all candidates. Sharing 2 of 2 is unremarkable; sharing 40 of 200 is not.", "",
        f"Candidates below degree {MIN_CANDIDATE_DEGREE} are dropped as untestable.", "",
        "## Candidates", "",
        "| candidate | bridges | expected | enrichment | q | PubMed co-mentions | via |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        cm = r.get("comentions")
        cm_s = "-" if cm is None or cm < 0 else (f"**{cm:,}**" if cm > 0 else "0")
        L.append(f"| {r['c_name']} | {r['bridges']} | {r['expected']:.1f} | "
                 f"{r['enrichment']:.1f}x | {r['q']:.1e} | {cm_s} | "
                 f"{', '.join(r['via'][:4])} |")

    if meta.get("validated"):
        share = meta["documented_share"]
        L += ["", "## Validation: does the ABC premise hold?", "",
              f"Of the top {meta['validated']} candidates checked against PubMed, "
              f"**{meta['documented']} ({share:.0%}) are ALREADY co-mentioned** in "
              "published abstracts.", ""]
        if share and share > 0.5:
            L += ["> **This seed produces no usable discoveries.** The ABC premise is that A",
                  "> and C have never been discussed together. Here that premise is false for",
                  "> most candidates, so the missing A-C edge reflects the relation",
                  "> extractor's RECALL, not an undiscovered link. GPX4 and caspase-3 share",
                  "> 236 PubMed abstracts and no graph edge.",
                  ">",
                  "> Literature-based discovery needs high edge recall to work at all, because",
                  "> it reasons from ABSENCE. This graph is built from abstract-level",
                  "> extraction, and its absences are not informative. Making LBD viable here",
                  "> would need full-text relation extraction or a curated knowledge base, and",
                  "> until then this tool is better read as a RECALL DIAGNOSTIC: the top rows",
                  "> are the edges the graph should already have.", ""]

    L += ["", "## Every row here is a hypothesis, not a finding", "",
          "* **Absence of an A-C edge is not absence of knowledge.** It may be missing",
          "  because the extractor failed, because the pair is discussed only in full",
          "  text, or because the relation is stated outside this cancer census.",
          "* **Direction of effect is not carried.** An A-B-C chain may compose two",
          "  relations whose signs cancel. It does not imply A affects C, or how.",
          "* **No context is carried**, so a chain can hold in one tissue and be",
          "  meaningless in another.",
          "* Many high-lift pairs will be trivially related -- synonyms, members of one",
          "  complex, a drug and its own target class -- rather than novel.", "",
          "The honest use is as a ranked reading list, and the honest next step for any",
          "row is to read the bridging papers and decide whether the chain composes.", ""]

    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print(f"seed {meta['seed_name']} deg={meta['degree']:,}, {meta['candidates']:,} candidates")
    for r in rows[:12]:
        print(f"  {r['c_name'][:32]:<34}bridges={r['bridges']:<4}exp={r['expected']:<7.1f}"
              f"{r['enrichment']:>6.1f}x  q={r['q']:.1e}  deg={r['deg_c']:<6}"
              f" via {', '.join(r['via'][:2])}")


if __name__ == "__main__":
    main()
