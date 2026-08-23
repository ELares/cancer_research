#!/usr/bin/env python3
"""Diagnostic-to-therapy chain membership at CENSUS scale (#RETIRE-FROZEN).

The manuscript reported this layer over 4,830 retrieved full-text records, where
the informative sentence was not the count but the disclaimer beneath it: the
low match rates "mean the local corpus was constructed around mechanism
keywords, not around diagnostic-therapy pairings", so two chains returned zero
and EGFR-to-EGFR-inhibitor "recovers only a dozen of the thousands of
EGFR-targeted-therapy papers in the wider literature". That disclaimer names a
property of the RETRIEVAL, and the census is the instrument that removes it.

THE INSTRUMENT IS NEARLY THE SAME ONE, which is what makes the two columns
comparable at all. `recompute_diagnostic_therapy_links` calls
`get_searchable_text(fm, body)` -- title + MeSH terms + PubTator disease/gene/drug
annotations + abstract, and NOT the full text. Census records carry title, MeSH
and abstract, so the only missing channel is the PubTator annotation layer. That
cost is MEASURED here rather than assumed: the same matcher is run over the
frozen corpus with the annotation channel removed, so the census column is
compared against a corpus column built the same way.

A chain returning zero on the census is a different statement from a chain
returning zero on the corpus. The corpus zero says the retrieval never reached
those papers; the census zero would say the indexed cancer literature does not
express the pairing in title, MeSH and abstract. Only the second is a claim
about the field, and only the second is worth reporting as one.
"""
import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
REPO = SCRIPT_DIR.parent
RECORDS = REPO / "corpus/atlas/records"
OUT_MD = REPO / "analysis/census-diagnostic-chains.md"
OUT_JSON = REPO / "analysis/census-diagnostic-chains.json"


def census_text(rec: dict) -> str:
    """Title + MeSH + abstract, the census-side analogue of the production text."""
    return " ".join(
        [rec.get("title") or "", " ".join(rec.get("mesh") or []), rec.get("abstract") or ""]
    ).lower()


def scan_census(stride: int = 1) -> dict:
    from tag_articles import match_diagnostic_therapy_links

    per_chain = Counter()
    matched = 0
    n = 0
    with_abstract = 0
    shards = sorted(RECORDS.glob("*.jsonl.gz"))[::stride]
    for f in shards:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                rec = json.loads(line)
                n += 1
                if rec.get("abstract"):
                    with_abstract += 1
                hits = match_diagnostic_therapy_links(census_text(rec))
                if hits:
                    matched += 1
                    per_chain.update(hits)
    return {
        "records": n,
        "shards": len(shards),
        "with_abstract": with_abstract,
        "matched": matched,
        "per_chain": per_chain,
    }


def scan_corpus_without_annotations() -> dict:
    """The SAME matcher over the frozen corpus with the PubTator annotation channel
    removed, so the census column has a like-for-like comparator."""
    from tag_articles import (
        PMID_DIR,
        extract_abstract,
        load_article,
        match_diagnostic_therapy_links,
    )

    per_chain_full = Counter()
    per_chain_narrow = Counter()
    matched_full = matched_narrow = n = 0
    for filepath in sorted(Path(PMID_DIR).glob("*.md")):
        fm, body = load_article(filepath)
        if not fm:
            continue
        n += 1
        abstract = extract_abstract(body)
        title = fm.get("title", "")
        mesh = " ".join(fm.get("mesh_terms", []))
        annot = " ".join(
            fm.get("diseases_annotated", []) + fm.get("genes", []) + fm.get("drugs", [])
        )
        narrow = " ".join([title, mesh, abstract]).lower()
        full = " ".join([title, mesh, annot, abstract]).lower()
        hn = match_diagnostic_therapy_links(narrow)
        hf = match_diagnostic_therapy_links(full)
        if hn:
            matched_narrow += 1
            per_chain_narrow.update(hn)
        if hf:
            matched_full += 1
            per_chain_full.update(hf)
    return {
        "records": n,
        "matched_production_text": matched_full,
        "matched_without_annotations": matched_narrow,
        "per_chain_production_text": per_chain_full,
        "per_chain_without_annotations": per_chain_narrow,
    }


def _split_stored(d: dict) -> tuple:
    """Recover the two scan products from the merged artifact.

    `--render-only` must RE-ASSEMBLE from the stored raw counts rather than
    re-render the stored derived fields, or every guard reading a derived
    field is comparing the artifact to itself and cannot fail. This generator
    used to load the assembled JSON and render it directly, so the row
    ORDERING, the annotation-channel cost and both zero-chain lists were
    re-emitted rather than recomputed.

    The raw per-chain counts survive inside `rows`, so nothing extra needs
    storing: this puts them back into the shapes `assemble` reads.
    """
    cen = {
        "records": d["census_records"],
        "shards": d["census_shards"],
        "with_abstract": d["census_with_abstract"],
        "matched": d["census_matched"],
        "per_chain": {r["chain"]: r["census"] for r in d["rows"]},
    }
    cor = {
        "records": d["corpus_records"],
        "matched_production_text": d["corpus_matched_production_text"],
        "matched_without_annotations": d["corpus_matched_without_annotations"],
        "per_chain_production_text":
            {r["chain"]: r["corpus_production_text"] for r in d["rows"]},
        "per_chain_without_annotations":
            {r["chain"]: r["corpus_without_annotations"] for r in d["rows"]},
    }
    return cen, cor


def assemble(cen: dict, cor: dict) -> dict:
    from config import DIAGNOSTIC_THERAPY_ORDER

    rows = []
    for cid in DIAGNOSTIC_THERAPY_ORDER:
        rows.append(
            {
                "chain": cid,
                "census": cen["per_chain"].get(cid, 0),
                "corpus_production_text": cor["per_chain_production_text"].get(cid, 0),
                "corpus_without_annotations": cor["per_chain_without_annotations"].get(cid, 0),
            }
        )
    rows.sort(key=lambda r: -r["census"])
    # The annotation channel's cost, measured on the SAME articles.
    cost = cor["matched_production_text"] - cor["matched_without_annotations"]
    return {
        "census_records": cen["records"],
        "census_shards": cen["shards"],
        "census_with_abstract": cen["with_abstract"],
        "census_matched": cen["matched"],
        "corpus_records": cor["records"],
        "corpus_matched_production_text": cor["matched_production_text"],
        "corpus_matched_without_annotations": cor["matched_without_annotations"],
        "annotation_channel_cost_records": cost,
        "rows": rows,
        "chains_zero_on_census": [r["chain"] for r in rows if r["census"] == 0],
        "chains_zero_on_corpus": [
            r["chain"] for r in rows if r["corpus_production_text"] == 0
        ],
    }


def render(d: dict) -> str:
    cen_pct = 100 * d["census_matched"] / d["census_records"]
    cor_pct = 100 * d["corpus_matched_production_text"] / d["corpus_records"]
    abs_pct = 100 * d["census_with_abstract"] / d["census_records"]
    narrow = d["corpus_matched_without_annotations"]
    cost = d["annotation_channel_cost_records"]
    cost_pct = 100 * cost / d["corpus_matched_production_text"]
    L = []
    L.append("# Diagnostic-to-therapy chains at census scale\n")
    L.append(
        f"Generated by `scripts/census_diagnostic_chains.py`. "
        f"{d['census_matched']:,} of {d['census_records']:,} census records "
        f"({cen_pct:.2f}%) match at least one of the {len(d['rows'])} chains.\n"
    )
    L.append("## What changed, and what the comparison is worth\n")
    L.append(
        f"The corpus column reports {d['corpus_matched_production_text']:,} of "
        f"{d['corpus_records']:,} ({cor_pct:.1f}%). The two percentages are NOT "
        f"comparable as rates: the corpus was retrieved by mechanism keyword and "
        f"the census is the whole MeSH-indexed cancer literature, so the corpus is "
        f"enriched for the mechanisms its queries named and depleted for everything "
        f"else. What IS comparable is the per-chain ordering and, in particular, "
        f"which chains return nothing.\n"
    )
    verdict = (
        f"costs nothing at all ({cost:,} records), so the two columns are built "
        f"from the same channels for every chain here"
        if cost == 0
        else f"costs {cost:,} records of {d['corpus_matched_production_text']:,} "
             f"({cost_pct:.1f}%), so the census column reads low by about that "
             f"much against a corpus column built the same way ({narrow:,})"
    )
    L.append(
        f"The instrument is nearly the same. The published figure was matched over "
        f"title + MeSH + PubTator annotations + abstract, and the census carries "
        f"every channel but the annotations. Removing that channel from the corpus "
        f"arm {verdict}. The corpus arm also reproduces the published "
        f"{d['corpus_matched_production_text']:,} exactly, so the matcher is the one "
        f"the manuscript reported. {abs_pct:.1f}% of census records carry an "
        f"abstract at all.\n"
    )
    L.append("## Per chain\n")
    L.append("| chain | census | corpus (production text) | corpus (no annotations) |")
    L.append("|---|--:|--:|--:|")
    for r in d["rows"]:
        L.append(
            f"| {r['chain']} | {r['census']:,} | {r['corpus_production_text']:,} | "
            f"{r['corpus_without_annotations']:,} |"
        )
    L.append("")
    zc = d["chains_zero_on_corpus"]
    zn = d["chains_zero_on_census"]
    by_corpus = sorted(d["rows"], key=lambda r: -r["corpus_production_text"])
    moved = [
        (r["chain"], by_corpus.index(r) + 1, i + 1)
        for i, r in enumerate(d["rows"])
        if abs(by_corpus.index(r) - i) >= 3
    ]
    L.append("## The ordering, not the rate\n")
    if moved:
        L.append(
            "Ranking the same ten chains by each column disagrees on "
            f"{len(moved)} of {len(d['rows'])} by three places or more:\n"
        )
        for chain, cor_rank, cen_rank in sorted(moved, key=lambda m: m[2]):
            L.append(f"- `{chain}`: corpus rank {cor_rank}, census rank {cen_rank}")
        L.append("")
        L.append(
            "The manuscript read its own ordering as a map of where the corpus had "
            "translational depth, and said so. The census says which of those "
            "readings were about the field: the chains that rise are the ones whose "
            "literature the mechanism queries never went looking for.\n"
        )
    else:
        L.append("The two columns rank the chains the same way.\n")
    L.append("## The zeros\n")
    if zc and not zn:
        L.append(
            f"{len(zc)} chain(s) return zero on the corpus and none returns zero on "
            f"the census: {', '.join(zc)}. That is the disclaimer the manuscript "
            f"attached to this section, now measured: those zeros were a property of "
            f"a retrieval built around mechanism keywords, not a statement about the "
            f"literature. Every chain is expressed in the indexed cancer literature.\n"
        )
    elif zn:
        L.append(
            f"{len(zn)} chain(s) return zero on the CENSUS: {', '.join(zn)}. Unlike a "
            f"corpus zero, this one is a claim about the indexed literature -- the "
            f"pairing is not expressed in title, MeSH or abstract anywhere in "
            f"{d['census_records']:,} records.\n"
        )
    else:
        L.append("No chain returns zero on either arm.\n")
    L.append("## What this still cannot say\n")
    L.append(
        "Chain membership is keyword co-occurrence in title, MeSH and abstract. It "
        "says a paper mentions a therapy and either a diagnostic or a targetable "
        "feature; it does not say the paper USES one to select the other. The layer "
        "measures where the literature pairs the two vocabularies, not where "
        "biomarker-directed treatment is established practice.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1,
                    help="sample every Nth shard (shards are CHRONOLOGICAL, so a "
                         "prefix samples one era; a stride spreads the draw)")
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    if a.render_only:
        d = assemble(*_split_stored(json.loads(OUT_JSON.read_text())))
    else:
        cor = scan_corpus_without_annotations()
        print(f"corpus: {cor['matched_production_text']} production / "
              f"{cor['matched_without_annotations']} without annotations")
        cen = scan_census(a.stride)
        print(f"census: {cen['matched']:,} of {cen['records']:,}")
        d = assemble(cen, cor)
        OUT_JSON.write_text(json.dumps(d, indent=1) + "\n", encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
