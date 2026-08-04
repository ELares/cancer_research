#!/usr/bin/env python3
"""An OFFLINE authority-name table for the co-mention identifiers (#628).

WHY
---
#628's discriminator asks whether a surface form is a NAME of the entity
according to an authority, rather than merely a form somebody annotated. The
index cannot answer that -- its `canon` field is the most frequent surface form,
derived from the same annotations, so testing against it is circular.

Two things blocked shipping the check. It needed a live NLM call per identifier,
which breaks this repository's offline contract, and NCBI Gene identifiers have
no MeSH descriptor at all, so the rule cut every one of them unconditionally.
That is not a corner case: of the 26,382 identifiers the alias map resolves to,
**12,906 are genes** -- in a repository whose subject is GPX4 and ACSL4.

This builds the table once, from the authorities, and commits the derived
result. CI and every downstream check then read a file.

SOURCES
-------
MeSH   NLM's SPARQL endpoint, the same one `atlas_baseline` already uses for the
       C04 cancer definition. One paginated query for descriptors and one for
       supplementary concept records, rather than 13,396 individual lookups.
NCBI   E-utilities esummary, batched. Each gene yields its official symbol, its
Gene   full description, and its listed aliases -- all three count as names, so
       `xCT` resolves as a name of SLC7A11 and `PHGPx` as one of GPX4.

WHAT IS AND IS NOT COMMITTED
----------------------------
Only the derived table, and only for identifiers this corpus actually uses. The
raw downloads are not committed and CI never fetches, matching the contract the
CTRPv2 and cBioPortal legs already follow.

A NOTE ON WHAT THIS DOES NOT SETTLE
-----------------------------------
Having gene names does not make the discriminator correct for genes. It makes it
TESTABLE for genes, which is the precondition. Whether a gene symbol matching its
official symbol predicts a true match is a separate measurement, and the judged
samples are the place to make it.

Usage:
    python scripts/build_label_source.py            # build, ~2 minutes
    python scripts/build_label_source.py --stats    # report the committed table
"""

import argparse
import gzip
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_comention import build_alias_map  # noqa: E402
from atlas_graph import load_index  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "comention" / "authority-labels.tsv.gz"
MESH_SPARQL = "https://id.nlm.nih.gov/mesh/sparql"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
UA = "cancer_research-atlas/1.0 (https://github.com/ELares/cancer_research)"


def _get(url: str, timeout: int = 180) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf-8", "replace")


def mesh_labels(uis: list) -> dict:
    """UI -> label, asked for the identifiers we need rather than all of MeSH.

    Dumping every concept of a class and paginating by offset returns HTTP 500:
    an unfiltered ORDER BY over roughly thirty thousand rows is more than the
    endpoint will do. A VALUES block naming the wanted identifiers is both
    reliable and far smaller, and it covers descriptors and supplementary
    records in one query without having to know which is which.
    """
    found = {}
    batch = 150
    for i in range(0, len(uis), batch):
        chunk = uis[i:i + batch]
        values = " ".join(f"mesh:{u}" for u in chunk)
        query = (
            "PREFIX mesh: <http://id.nlm.nih.gov/mesh/> "
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> "
            f"SELECT ?d ?label WHERE {{ VALUES ?d {{ {values} }} ?d rdfs:label ?label }}"
        )
        try:
            params = urllib.parse.urlencode({"query": query, "format": "JSON"})
            rows = json.loads(_get(f"{MESH_SPARQL}?{params}"))["results"]["bindings"]
        except Exception as exc:
            print(f"\n    ! MeSH batch at {i} failed: {exc}", file=sys.stderr)
            time.sleep(2)
            continue
        for b in rows:
            found[b["d"]["value"].rsplit("/", 1)[-1]] = b["label"]["value"]
        print(f"    mesh: {len(found):,}/{len(uis):,}", end="\r", flush=True)
        time.sleep(0.3)
    print(f"    mesh: {len(found):,}/{len(uis):,}          ")
    return found


def gene_labels(gene_ids: list) -> dict:
    """Gene id -> its official symbol, description and aliases.

    All three are names by NCBI's own account, so a form matching any of them is
    a name of that gene. Restricting to the official symbol would reject `xCT`
    for SLC7A11 and `PHGPx` for GPX4, which are how the literature writes them.
    """
    found = {}
    for i in range(0, len(gene_ids), 250):
        batch = gene_ids[i:i + 250]
        params = urllib.parse.urlencode(
            {"db": "gene", "id": ",".join(batch), "retmode": "json"})
        try:
            res = json.loads(_get(f"{ESUMMARY}?{params}"))["result"]
        except Exception as exc:
            print(f"\n    ! batch at {i} failed: {exc}", file=sys.stderr)
            time.sleep(2)
            continue
        for uid in res.get("uids", []):
            r = res[uid]
            names = [r.get("name"), r.get("description")]
            names += [a.strip() for a in (r.get("otheraliases") or "").split(",")]
            found[uid] = [n for n in names if n]
        print(f"    genes: {len(found):,}/{len(gene_ids):,}", end="\r", flush=True)
        time.sleep(0.34)
    print(f"    genes: {len(found):,}/{len(gene_ids):,}          ")
    return found


def top_unreachable(root, reachable: set, limit: int) -> set:
    """The most-annotated identifiers no alias form resolves to.

    Included so the table can answer questions about entities the layer misses,
    not only about the ones it finds.
    """
    import collections
    import gzip as _gzip

    mentions = collections.Counter()
    for kind in ("gene", "chemical", "disease"):
        f = root / "entities" / f"{kind}.tsv.gz"
        if not f.exists():
            continue
        with _gzip.open(f, "rt", errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3 and p[2] and p[2] != "-" and p[2] not in reachable:
                    mentions[p[2]] += 1
    return {i for i, _ in mentions.most_common(limit)}


def load_table() -> dict:
    """identifier -> [names], from the committed table."""
    if not OUT.exists():
        return {}
    out = {}
    with gzip.open(OUT, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ident, names = line.rstrip("\n").split("\t", 1)
            out[ident] = names.split("|")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stats", action="store_true")
    args = ap.parse_args()

    if args.stats:
        t = load_table()
        if not t:
            print(f"no table at {OUT}", file=sys.stderr)
            return 1
        genes = sum(1 for k in t if not k.startswith(("MESH:", "OMIM:")))
        print(f"{OUT}\n  identifiers : {len(t):,}\n  genes       : {genes:,}"
              f"\n  mesh/omim   : {len(t)-genes:,}"
              f"\n  names total : {sum(len(v) for v in t.values()):,}")
        return 0

    idx = load_index(atlas_root())
    alias, _ = build_alias_map(idx)
    wanted = set(alias.values())
    # Also cover the most-annotated identifiers the alias map CANNOT reach.
    # Without them the table can only describe entities the layer already
    # finds, which makes "is this entity unreachable by its own name?"
    # unanswerable by construction -- the question that motivated the table.
    wanted |= top_unreachable(atlas_root(), wanted, limit=3000)
    mesh_ids = {i.split(":", 1)[1] for i in wanted if i.startswith("MESH:")}
    gene_ids = sorted(i for i in wanted if not i.startswith(("MESH:", "OMIM:")))
    print(f"identifiers in use: {len(wanted):,} "
          f"({len(mesh_ids):,} MeSH, {len(gene_ids):,} gene)")

    print("  fetching MeSH labels ...")
    mesh = mesh_labels(sorted(mesh_ids))

    print("  fetching gene names ...")
    genes = gene_labels(gene_ids)

    rows = {}
    for ident in sorted(wanted):
        if ident.startswith("MESH:"):
            lab = mesh.get(ident.split(":", 1)[1])
            if lab:
                rows[ident] = [lab]
        elif not ident.startswith("OMIM:"):
            names = genes.get(ident)
            if names:
                rows[ident] = names

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT, "wt", encoding="utf-8") as fh:
        fh.write("# Authority names for the identifiers the co-mention alias map "
                 "resolves to (#628).\n")
        fh.write("# MeSH: NLM SPARQL endpoint. Genes: NCBI E-utilities esummary "
                 "(symbol, description, aliases).\n")
        fh.write("# Regenerate: python scripts/build_label_source.py\n")
        fh.write(f"# identifiers: {len(rows)} of {len(wanted)} in use\n")
        for ident, names in rows.items():
            fh.write(ident + "\t" + "|".join(names) + "\n")

    miss_m = len(mesh_ids) - sum(1 for k in rows if k.startswith("MESH:"))
    miss_g = len(gene_ids) - sum(1 for k in rows if not k.startswith(("MESH:", "OMIM:")))
    print(f"\nwrote {OUT}: {len(rows):,} identifiers")
    print(f"  unresolved MeSH : {miss_m:,}")
    print(f"  unresolved gene : {miss_g:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
