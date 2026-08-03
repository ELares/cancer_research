#!/usr/bin/env python3
"""Atlas: normalized entities and typed relations for the cancer census (#ATLAS).

WHY
---
The frozen corpus stores entities as raw surface strings, because
`scripts/enrich_metadata.py` calls PubTator3 and then discards
`infons['identifier']`. The result is 3,260 distinct gene spellings for roughly
1,200 genes, `PD-1` and `PD1` as different entities, and `MTT` -- an assay
reagent -- filed as a drug. Nothing joins to HGNC, MeSH or anything else, and no
relation between entities is extracted anywhere in the project.

NCBI already publishes, for the whole of PubMed, exactly what is needed: entity
annotations resolved to stable identifiers, and typed relations between them,
public domain. This ingests those bulk files and keeps the cancer slice.

    relation2pubtator3.gz   PMID \\t type \\t Type|ID \\t Type|ID
    gene2pubtator3.gz       PMID \\t NCBI Gene id \\t mention \\t resolver
    disease2pubtator3.gz    PMID \\t MeSH id      \\t mention \\t resolver
    chemical2pubtator3.gz   PMID \\t MeSH id      \\t mention \\t resolver

Relation types are a fixed vocabulary (associate, treat, cause, inhibit,
stimulate, positive_correlate, negative_correlate, cotreat, interact, compare,
prevent, drug_interact), which is what makes the result queryable rather than a
pile of free-text predicates.

WHAT IT WRITES
--------------
    <root>/relations/relations.tsv.gz   cancer-PMID relations, source columns kept
    <root>/entities/<kind>.tsv.gz       cancer-PMID entity annotations
    <root>/relations/manifest.json      per-source state so runs resume

Requires the census from `scripts/atlas_baseline.py` -- the cancer PMID set is
what the bulk files are filtered against.

Stdlib only. Downloads resume via HTTP Range.

Usage:
    python scripts/atlas_relations.py --only relation --sample 200000   # smoke test
    python scripts/atlas_relations.py                                   # full ingest
"""

import argparse
import gzip
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402

BASE = "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/"
USER_AGENT = "cancer_research-atlas/1.0 (https://github.com/ELares/cancer_research)"

SOURCES = {
    "relation": "relation2pubtator3.gz",
    "gene": "gene2pubtator3.gz",
    "disease": "disease2pubtator3.gz",
    "chemical": "chemical2pubtator3.gz",
}


def load_cancer_pmids(root: Path) -> set:
    """The census PMID set, as ints (~150 MB at full scale, vs ~400 MB as str)."""
    rec_dir = root / "records"
    files = sorted(rec_dir.glob("*.jsonl.gz"))
    if not files:
        raise SystemExit(f"no census records under {rec_dir}; run scripts/atlas_baseline.py first")
    pmids = set()
    for f in files:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                p = json.loads(line).get("pmid", "")
                if p.isdigit():
                    pmids.add(int(p))
    return pmids


def download_resumable(name: str, dest: Path, quiet: bool = False) -> Path:
    """Download with HTTP Range resume. These files are multi-GB."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    pos = dest.stat().st_size if dest.exists() else 0
    headers = {"User-Agent": USER_AGENT}
    if pos:
        headers["Range"] = f"bytes={pos}-"
    req = urllib.request.Request(BASE + name, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=600) as resp:
            if resp.status == 416:  # already complete
                return dest
            mode = "ab" if pos and resp.status == 206 else "wb"
            if mode == "wb":
                pos = 0
            with open(dest, mode) as fh:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
                    pos += len(chunk)
                    if not quiet and pos % (200 << 20) < (1 << 20):
                        print(f"    {pos / (1 << 30):.1f} GB", flush=True)
    except urllib.error.HTTPError as e:
        if e.code != 416:
            raise
    return dest


def filter_source(kind: str, src: Path, out: Path, pmids: set, sample: int = 0) -> dict:
    """Stream a bulk file, keep rows whose PMID is in the cancer census."""
    kept = seen = 0
    types = {}
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rt", encoding="utf-8", errors="ignore") as fh, \
            gzip.open(out, "wt", encoding="utf-8") as w:
        for line in fh:
            seen += 1
            if sample and seen > sample:
                break
            tab = line.find("\t")
            if tab <= 0:
                continue
            head = line[:tab]
            if not head.isdigit() or int(head) not in pmids:
                continue
            w.write(line)
            kept += 1
            if kind == "relation":
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    types[parts[1]] = types.get(parts[1], 0) + 1
    return {"seen": seen, "kept": kept, "types": types}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--only", choices=sorted(SOURCES), help="ingest a single source")
    ap.add_argument("--sample", type=int, default=0,
                    help="stop after N input lines (smoke test)")
    ap.add_argument("--keep-raw", action="store_true",
                    help="keep the downloaded bulk file (default: delete after filtering)")
    args = ap.parse_args()

    root = atlas_root()
    print("loading census PMIDs ...", flush=True)
    pmids = load_cancer_pmids(root)
    print(f"cancer PMIDs in census: {len(pmids):,}")

    raw = root / "raw"
    man_path = root / "relations" / "manifest.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"sources": {}}

    todo = [args.only] if args.only else list(SOURCES)
    for kind in todo:
        name = SOURCES[kind]
        print(f"\n[{kind}] {name}")
        t0 = time.time()
        src = raw / name
        download_resumable(name, src)
        out = (root / ("relations" if kind == "relation" else "entities")
               / (("relations" if kind == "relation" else kind) + ".tsv.gz"))
        stats = filter_source(kind, src, out, pmids, sample=args.sample)
        if not args.keep_raw and not args.sample:
            src.unlink(missing_ok=True)
        stats["census_pmids"] = len(pmids)
        stats["seconds"] = round(time.time() - t0, 1)
        man["sources"][kind] = stats
        man_path.parent.mkdir(parents=True, exist_ok=True)
        man_path.write_text(json.dumps(man, indent=1, sort_keys=True), encoding="utf-8")
        print(f"  scanned {stats['seen']:,} rows -> kept {stats['kept']:,} "
              f"({stats['kept']/max(stats['seen'],1):.1%}) in {stats['seconds']}s")
        if stats.get("types"):
            top = sorted(stats["types"].items(), key=lambda kv: -kv[1])[:12]
            print("  relation types: " + ", ".join(f"{k}={v:,}" for k, v in top))


if __name__ == "__main__":
    main()
