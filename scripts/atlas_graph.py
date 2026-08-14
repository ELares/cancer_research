#!/usr/bin/env python3
"""Atlas: query the cancer literature's typed relation graph (#ATLAS).

WHAT THIS IS FOR
----------------
The simulation suite has ~30 off-by-default ferroptosis "realism layers", and
each was added on the strength of one or two hand-picked papers -- the module
docs cite them by PMID. That is the weakest form of literature grounding: it
cannot distinguish a mechanism the field has replicated a hundred times from one
asserted once in a single paper the author happened to read.

The atlas can now answer that with a denominator. It holds the census's typed,
normalized relation layer over the cancer-article PMIDs it was extracted from
(the row count is in `corpus/atlas/relations/manifest.json`, and is not
repeated here because it grew by about a third when the layer was re-ingested
and this sentence did not), drawn from NCBI's
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
import hashlib
import json
import pickle
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402

MIN_MENTIONS = 3  # a surface form must appear this often to enter the alias map

# Sense collisions measured by scripts/atlas_ambiguity.py. resolve() refuses to
# guess on these rather than returning a plausible-looking wrong gene.
AMBIGUITY_JSON = Path(__file__).resolve().parents[1] / "analysis" / "atlas-ambiguity.json"
_AMBIG = None

# Per-pair PMIDs are sampled, not truncated. The first version kept
# `sorted(v)[:50]` over PMID STRINGS, which is a lexicographic (first-digit)
# order, so the retained PMIDs were systematically the oldest -- and
# atlas_emergence.py then computed 'share of support since 2021' on a sample
# built to exclude recent papers. Measured on the affected pairs, median true
# recent share 26.6% collapsed to 0.0% as stored. A uniform reservoir sample is
# unbiased, and `n_pmids` records the TRUE support size so any share computed
# from the sample has a real denominator.
PMID_SAMPLE = 60
_RESERVOIR_SEED = 20260803


# Per-paper sense corrections from scripts/atlas_disambiguate.py, applied when
# the relation edges are read so every consumer of the index sees the right
# gene. Without this the manuscript's own GPX4+FSP1 claim has ZERO typed
# relations: PubTator extracted them, then filed all of them under ATL1.
DISAMBIGUATION_JSON = (Path(__file__).resolve().parents[1]
                       / "analysis" / "atlas-disambiguation.json")
# sense -> the identifier a corrected mention should carry
_SENSE_ID = {"AIFM2": "84883", "S100A4": "6275", "ATL1": "51062"}
# Only the ATL1 identifiers may be moved, and the restriction is deliberate.
# relations.tsv.gz records identifiers, not the surface form that produced
# them, so a paper discussing BOTH cancer-associated fibroblasts and
# ferroptosis would have its genuine S100A4 edge rewritten to AIFM2 if every
# colliding id were correctable. ATL1 is safe to move because it is a
# hereditary spastic-paraplegia gene that only 1.9% of these papers mention at
# all, whereas a paper that really means S100A4 or AIFM2 gets that identifier
# from its own unambiguous surface form. This trades recall for precision: the
# 33 papers corrected S100A4 -> AIFM2 keep their original edges.
_CORRECTABLE = {"51062", "73991"}


def load_corrections() -> dict:
    """pmid -> corrected identifier, for the papers the disambiguator decided."""
    try:
        with open(DISAMBIGUATION_JSON) as fh:
            raw = json.load(fh).get("corrections", {})
    except (OSError, ValueError):
        return {}
    out = {}
    for pmid, v in raw.items():
        ident = _SENSE_ID.get(v.get("corrected"))
        if ident:
            out[pmid] = ident
    return out


def _corrected(ident: str, pmid: str, corrections: dict) -> str:
    """Remap a colliding identifier to the sense this paper actually uses."""
    if ident in _CORRECTABLE:
        return corrections.get(pmid, ident)
    return ident


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
    # How often each surface form appears at all. Computed here already and
    # previously discarded, which cost the co-mention layer dearly: `MIN_MENTIONS`
    # of 3 admits a generic English phrase mis-annotated in a handful of
    # abstracts, and the majority vote then hands it an identifier. `tumor cells`
    # resolves to Glucagonoma on 312 mentions, `overall survival` to Prosthesis
    # Failure on 37, `et al` to Multiple Myeloma on 11. Support separates those
    # cleanly from real matches -- measured on a labelled sample, correct
    # matches have a median support of 1,067 and wrong ones 16 -- so consumers
    # that read running text need it. See analysis/atlas-comention-audit.md.
    alias_support = {a: sum(c.values()) for a, c in alias.items()
                     if sum(c.values()) >= MIN_MENTIONS}
    # Total mentions per identifier, so a consumer can ask what SHARE of an
    # entity's mentions a given surface form accounts for. That share is the
    # signal that separates a real name from a mis-annotation without needing
    # any authority-name lookup: `tumor cells` is 2.8% of Glucagonoma's
    # mentions and `et al` is 0.006% of Multiple Myeloma's, while `gpx4` is
    # 63% of GPX4's and `erastin` 98.5% of erastin's.
    ident_mentions = {i: sum(c.values()) for i, c in canon.items()}
    # The count of a form FOR THE IDENTIFIER IT RESOLVES TO, which is the only
    # numerator that makes the ratio above an actual share.
    #
    # `alias_support` sums a form across every sense it carries, while
    # `ident_mentions` sums an identifier across every form. Dividing one by the
    # other is not a share and can exceed 1 -- measured at 274% for `as`, 133%
    # for `gp`, 104% for `tss`. Worse, it is biased in exactly the wrong
    # direction: an ambiguous generic word collects a LARGER numerator from its
    # other senses, so a minimum-share filter admits it more readily than a
    # specific name. That is how `as`, `treatment` and `effects` passed a filter
    # written to exclude them (analysis/atlas-comention-audit.md).
    alias_ident_support = {a: c.most_common(1)[0][1] for a, c in alias.items()
                           if sum(c.values()) >= MIN_MENTIONS}
    canon_res = {i: c.most_common(1)[0][0] for i, c in canon.items()}
    print(f"  aliases {len(alias_res):,}, identifiers {len(canon_res):,}", flush=True)

    corrections = load_corrections()
    if corrections:
        print(f"  applying {len(corrections):,} per-paper sense corrections", flush=True)

    print("  reading relations ...", flush=True)
    edges: dict[tuple, collections.Counter] = collections.defaultdict(collections.Counter)
    sample: dict[tuple, list] = collections.defaultdict(list)
    seen_pmids: dict[tuple, set] = collections.defaultdict(set)
    # One RNG PER PAIR, seeded from the pair itself. A single shared RNG is
    # consumed in file order, so changing any pair's membership shifts every
    # later draw and a rebuild reshuffles the sample for unrelated pairs too --
    # applying 441 sense corrections rewrote the example PMIDs of MDM2-p53,
    # which shares no paper with any of them. That makes rebuild diffs
    # unreadable and hides which pairs actually moved. Seeding per pair means a
    # pair's sample changes only when that pair's own evidence changes.
    rngs: dict[tuple, random.Random] = {}

    def pair_rng(key: tuple) -> random.Random:
        r = rngs.get(key)
        if r is None:
            # hashlib, not hash(): PYTHONHASHSEED salts str hashing per process,
            # so hash() would make the sample irreproducible across runs.
            digest = hashlib.blake2b(
                f"{_RESERVOIR_SEED}|{key[0]}|{key[1]}".encode(),
                digest_size=8).digest()
            r = random.Random(int.from_bytes(digest, "big"))
            rngs[key] = r
        return r

    with gzip.open(p["relations"], "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, rel, a, b = parts[0], parts[1], parts[2], parts[3]
            ida = _corrected(a.split("|", 1)[-1], pmid, corrections)
            idb = _corrected(b.split("|", 1)[-1], pmid, corrections)
            key = tuple(sorted((ida, idb)))
            edges[key][rel] += 1
            # one pair can assert several predicates from the same paper; count
            # each PMID once so `n_pmids` is a distinct-article support size
            seen = seen_pmids[key]
            if pmid in seen:
                continue
            seen.add(pmid)
            # Algorithm R: uniform reservoir over the DISTINCT PMIDs
            res = sample[key]
            n = len(seen)
            if len(res) < PMID_SAMPLE:
                res.append(pmid)
            else:
                j = pair_rng(key).randrange(n)
                if j < PMID_SAMPLE:
                    res[j] = pmid

    idx = {"alias": alias_res, "alias_support": alias_support,
           "alias_ident_support": alias_ident_support,
           "ident_mentions": ident_mentions, "canon": canon_res,
           "edges": {k: dict(v) for k, v in edges.items()},
           "pmids": {k: sorted(v, key=int) for k, v in sample.items()},
           "n_pmids": {k: len(v) for k, v in seen_pmids.items()}}
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


def _ambiguity():
    """The committed sense-collision blocklist, loaded once. Empty if absent."""
    global _AMBIG
    if _AMBIG is None:
        try:
            with open(AMBIGUITY_JSON) as fh:
                d = json.load(fh)
            _AMBIG = (set(d.get("blocklist", [])), d.get("domain_sense", {}))
        except (OSError, ValueError):
            _AMBIG = (set(), {})
    return _AMBIG


def resolve(idx: dict, name: str, allow_domain_sense: bool = False):
    """Surface form -> identifier, or the identifier itself if given directly.

    Returns None for a symbol that `atlas_ambiguity.py` measured as a genuine
    SENSE collision (different genes, not merely different species). The alias
    map resolves those by majority vote, which returns one identifier with no
    signal that anything was discarded -- so `ER` silently becomes epiregulin
    rather than the estrogen receptor, and `FSP1` becomes a spastic-paraplegia
    gene. Failing loudly is the whole point; call `resolve_reason` for the
    explanation, or pass allow_domain_sense=True to accept the curated
    cancer-domain sense where one is defensible.
    """
    n = name.strip()
    key = n.lower()
    blocked, domain = _ambiguity()
    # Checked BEFORE the canon shortcut: PubTator's canonical name for gene
    # 51062 is itself the string "FSP1", so a blocked symbol can also be a
    # canon key and would otherwise resolve straight through the block.
    if key in blocked:
        if allow_domain_sense and key in domain:
            return domain[key]["id"]
        return None
    if n in idx["canon"]:
        return n
    return idx["alias"].get(key)


def resolve_majority(idx: dict, name: str):
    """The raw majority-vote identifier, blocklist ignored.

    For callers whose job is to REPORT what the vote does -- the entity audit --
    rather than to rely on it. Analysis code should use resolve().
    """
    n = name.strip()
    if n in idx["canon"]:
        return n
    return idx["alias"].get(n.lower())


def resolve_reason(idx: dict, name: str) -> str:
    """Why resolve() returned None, in words a caller can print."""
    key = name.strip().lower()
    blocked, domain = _ambiguity()
    if key in blocked:
        majority = idx["alias"].get(key)
        # Name the identifier, not just PubTator's label for it: the label for
        # gene 51062 is the string "FSP1", so "returns FSP1" reads as agreement
        # when it is in fact the collision.
        got = (f"{idx['canon'].get(majority, majority)} (NCBI Gene {majority})"
               if majority else "nothing")
        msg = (f"`{name}` is a measured sense collision: the majority vote "
               f"returns {got}, which may not be the sense you mean "
               f"(analysis/atlas-ambiguity.md).")
        if key in domain:
            d = domain[key]
            msg += (f" The cancer-domain sense is {d['symbol']} ({d['id']}): "
                    f"{d['why']} Pass allow_domain_sense=True to accept it.")
        else:
            msg += " No domain default is defensible; disambiguate per paper."
        return msg
    if not idx["alias"].get(key):
        return f"`{name}` is not in the alias map (needs >= {MIN_MENTIONS} mentions)."
    return ""


def support(idx: dict, a: str, b: str):
    ida, idb = resolve(idx, a), resolve(idx, b)
    if not ida or not idb:
        return None
    key = tuple(sorted((ida, idb)))
    return {"a": ida, "a_name": idx["canon"].get(ida, ida),
            "b": idb, "b_name": idx["canon"].get(idb, idb),
            "predicates": idx["edges"].get(key, {}),
            "total": sum(idx["edges"].get(key, {}).values()),
            # distinct articles asserting ANY relation for this pair
            "n_articles": idx.get("n_pmids", {}).get(key, 0),
            # a UNIFORM SAMPLE of those articles, not a prefix
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
