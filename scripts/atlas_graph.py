#!/usr/bin/env python3
"""Atlas: query the cancer literature's typed relation graph (#ATLAS).

WHAT THIS IS FOR
----------------
The simulation suite has ~30 off-by-default ferroptosis "realism layers", and
each was added on the strength of one or two hand-picked papers -- the module
docs cite them by PMID. That is the weakest form of literature grounding: it
cannot distinguish a mechanism the field has replicated a hundred times from one
asserted once in a single paper the author happened to read.

The atlas can now answer that with a denominator. It holds 7,951,325 typed,
normalized relations over 1,603,105 cancer-article PMIDs, drawn from NCBI's
PubTator3 bulk release, with a FIXED predicate vocabulary (associate, treat,
cause, inhibit, stimulate, positive_correlate, negative_correlate, cotreat,
interact, compare, prevent, drug_interact).

So for any claim of the form "X relates to Y", this reports how many distinct
cancer articles assert it, under which predicates, and gives the PMIDs.

WHAT IT IS NOT
--------------
Co-occurrence in an abstract is not proof of mechanism, and PubTator's relation
extraction has its own error rate (BioREx scores ~79.6 F1 on BioRED). A high
count means the field talks about the pair, not that the pair is true. Roughly
half of all relations here are the weakest predicate, `associate`, which is
closer to co-mention than to knowledge -- the per-predicate breakdown is
reported precisely so that weight can be discounted.

Usage:
    python scripts/atlas_graph.py --build              # one-time index (~2 min)
    python scripts/atlas_graph.py --pair GPX4 ferroptosis
    python scripts/atlas_graph.py --entity GPX4 --top 25
"""

import argparse
import collections
import gzip
import json
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402

MIN_MENTIONS = 3  # a surface form must appear this often to enter the alias map


def paths(root: Path):
    return {
        "relations": root / "relations" / "relations.tsv.gz",
        "index": root / "relations" / "graph-index.pkl",
        "entities": root / "entities",
    }


def build_index(root: Path) -> dict:
    """Alias -> id map, id -> canonical name, and the relation edge list."""
    p = paths(root)
    alias: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    canon: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)

    for kind in ("gene", "chemical", "disease"):
        f = p["entities"] / f"{kind}.tsv.gz"
        if not f.exists():
            continue
        print(f"  reading {kind} ...", flush=True)
        with gzip.open(f, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                _pmid, _typ, ident, mention = parts[0], parts[1], parts[2], parts[3]
                if not ident or ident in ("-", "None"):
                    continue
                # PubTator packs synonyms into the mention with '|'
                for surface in mention.split("|"):
                    s = surface.strip().lower()
                    if s:
                        alias[s][ident] += 1
                        canon[ident][surface.strip()] += 1

    # keep only reasonably-attested surface forms, and resolve each to its
    # most frequent identifier (surface forms are genuinely ambiguous)
    alias_res = {a: c.most_common(1)[0][0] for a, c in alias.items()
                 if sum(c.values()) >= MIN_MENTIONS}
    canon_res = {i: c.most_common(1)[0][0] for i, c in canon.items()}
    print(f"  aliases {len(alias_res):,}, identifiers {len(canon_res):,}", flush=True)

    print("  reading relations ...", flush=True)
    edges: dict[tuple, collections.Counter] = collections.defaultdict(collections.Counter)
    pmids: dict[tuple, set] = collections.defaultdict(set)
    with gzip.open(p["relations"], "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, rel, a, b = parts[0], parts[1], parts[2], parts[3]
            ida = a.split("|", 1)[-1]
            idb = b.split("|", 1)[-1]
            key = tuple(sorted((ida, idb)))
            edges[key][rel] += 1
            if len(pmids[key]) < 200:  # cap: enough to cite, bounded memory
                pmids[key].add(pmid)

    idx = {"alias": alias_res, "canon": canon_res,
           "edges": {k: dict(v) for k, v in edges.items()},
           "pmids": {k: sorted(v)[:50] for k, v in pmids.items()}}
    p["index"].parent.mkdir(parents=True, exist_ok=True)
    with open(p["index"], "wb") as fh:
        pickle.dump(idx, fh, protocol=5)
    print(f"  wrote {p['index']} ({p['index'].stat().st_size/1e6:.0f} MB), "
          f"{len(idx['edges']):,} distinct entity pairs", flush=True)
    return idx


def load_index(root: Path) -> dict:
    p = paths(root)["index"]
    if not p.exists():
        raise SystemExit(f"no index at {p}; run: python scripts/atlas_graph.py --build")
    with open(p, "rb") as fh:
        return pickle.load(fh)


def resolve(idx: dict, name: str):
    """Surface form -> identifier, or the identifier itself if given directly."""
    n = name.strip()
    if n in idx["canon"]:
        return n
    return idx["alias"].get(n.lower())


def support(idx: dict, a: str, b: str):
    ida, idb = resolve(idx, a), resolve(idx, b)
    if not ida or not idb:
        return None
    key = tuple(sorted((ida, idb)))
    return {"a": ida, "a_name": idx["canon"].get(ida, ida),
            "b": idb, "b_name": idx["canon"].get(idb, idb),
            "predicates": idx["edges"].get(key, {}),
            "total": sum(idx["edges"].get(key, {}).values()),
            "pmids": idx["pmids"].get(key, [])}


def neighbours(idx: dict, name: str, top: int = 20):
    ident = resolve(idx, name)
    if not ident:
        return None, []
    out = []
    for key, preds in idx["edges"].items():
        if ident in key:
            other = key[0] if key[1] == ident else key[1]
            out.append((sum(preds.values()), other, idx["canon"].get(other, other), preds))
    out.sort(reverse=True)
    return ident, out[:top]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--pair", nargs=2, metavar=("A", "B"))
    ap.add_argument("--entity", metavar="NAME")
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    root = atlas_root()
    if args.build:
        build_index(root)
        return

    idx = load_index(root)
    if args.pair:
        r = support(idx, *args.pair)
        if not r:
            raise SystemExit(f"could not resolve one of: {args.pair}")
        print(f"{r['a_name']} ({r['a']})  <->  {r['b_name']} ({r['b']})")
        print(f"  asserting relations: {r['total']:,}")
        for k, v in sorted(r["predicates"].items(), key=lambda kv: -kv[1]):
            print(f"    {k:22} {v:>8,}")
        if r["pmids"]:
            print(f"  example PMIDs: {', '.join(r['pmids'][:10])}")
    elif args.entity:
        ident, ns = neighbours(idx, args.entity, args.top)
        if not ident:
            raise SystemExit(f"could not resolve {args.entity!r}")
        print(f"{idx['canon'].get(ident, ident)} ({ident}) — top {len(ns)} partners\n")
        print(f"  {'partner':<34}{'n':>8}  predicates")
        for n, _oid, oname, preds in ns:
            top = ", ".join(f"{k}={v}" for k, v in
                            sorted(preds.items(), key=lambda kv: -kv[1])[:3])
            print(f"  {oname[:32]:<34}{n:>8,}  {top}")
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
