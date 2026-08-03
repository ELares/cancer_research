#!/usr/bin/env python3
"""Atlas: ingest the whole cancer literature from the PubMed baseline (#ATLAS).

WHY THIS EXISTS
---------------
The frozen corpus is 4,830 full-text articles. PubMed indexes 4,276,707 articles
under the MeSH tree C04 (Neoplasms). The corpus is therefore **0.11% of the
cancer literature**, and it was assembled from 19 hand-written keyword queries
with a 500-record-per-query cap and inconsistent date windows. Every statement
the project makes about "where research is concentrated" or "where the gaps are"
is really a statement about that retrieval design, not about the field.

This script replaces the keyword-query corpus with a census.

WHY THE BASELINE AND NOT THE API
--------------------------------
E-utilities cannot deliver the set. `esearch` over `neoplasms[mh]` reports all
4,276,707 hits, but paging is capped at retstart 9,999 even with WebEnv history
(verified: retstart=10000 returns an error document, not records). Date
partitioning to stay under the cap needs thousands of requests and still drifts
as records are reindexed. The annual baseline has no such limit: 1,334 gzipped
XML files, ~19 MB each, ~30,000 records each, ~25 GB total, and it is the same
data NLM ships to every other consumer.

WHY MeSH AND NOT KEYWORDS
-------------------------
"Cancer" is defined here as *any MeSH descriptor in tree C04*, resolved from
NLM's own MeSH SPARQL endpoint and cached in `mesh/c04-descriptors.tsv`. That is
a controlled-vocabulary definition maintained by NLM indexers, reproducible from
a committed file, and it does not inherit the coverage holes of a hand-written
keyword list -- the current corpus has no photodynamic-therapy query and no
ferroptosis query at all, despite both being central to the project's own thesis.

WHAT IT WRITES
--------------
    <root>/mesh/c04-descriptors.tsv     the cancer definition (committed, small)
    <root>/raw/pubmed26nNNNN.xml.gz     downloaded baseline files (NOT committed)
    <root>/records/part-NNNN.jsonl.gz   parsed cancer records, one JSON per line
    <root>/manifest.json                per-file state, so runs resume

Nothing under `corpus/by-pmid/`, `corpus/INDEX.jsonl` or `tags/` is touched. This
is a new surface alongside the frozen corpus and the living review.

The data root defaults to `corpus/atlas/` and is overridable with
`FERRO_ATLAS_ROOT`, so the bulk can live on external storage while the small
committed artifacts stay in the repo.

Stdlib only: urllib, gzip, xml.etree, json. No new dependencies.

Usage:
    python scripts/atlas_baseline.py --limit 2        # smoke test on 2 files
    python scripts/atlas_baseline.py                  # full ingest (resumable)
    python scripts/atlas_baseline.py --status         # progress without fetching
"""

import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

BASELINE_URL = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"
MESH_SPARQL = "https://id.nlm.nih.gov/mesh/sparql"
USER_AGENT = "cancer_research-atlas/1.0 (https://github.com/ELares/cancer_research)"

DEFAULT_ROOT = PROJECT_ROOT / "corpus" / "atlas"


def atlas_root() -> Path:
    return Path(os.getenv("FERRO_ATLAS_ROOT", str(DEFAULT_ROOT)))


def _get(url: str, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


# --------------------------------------------------------------------------
# The cancer definition: MeSH tree C04
# --------------------------------------------------------------------------

def fetch_c04_descriptors(dest: Path, force: bool = False) -> dict:
    """Descriptor UI -> label for every topical descriptor under MeSH tree C04.

    Cached to `dest`. NLM's SPARQL endpoint paginates, so this walks offsets
    until a page comes back short.
    """
    if dest.exists() and not force:
        out = {}
        for line in dest.read_text(encoding="utf-8").splitlines():
            if line and not line.startswith("#"):
                ui, label = line.split("\t", 1)
                out[ui] = label
        return out

    query = (
        "PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#> "
        "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> "
        "SELECT DISTINCT ?d ?label WHERE { "
        "?d a meshv:TopicalDescriptor . ?d meshv:treeNumber ?t . ?d rdfs:label ?label . "
        'FILTER(STRSTARTS(STR(?t), "http://id.nlm.nih.gov/mesh/C04")) } ORDER BY ?d'
    )
    found: dict[str, str] = {}
    offset, page = 0, 500
    while True:
        params = urllib.parse.urlencode({
            "query": query, "format": "JSON", "inference": "true",
            "limit": page, "offset": offset,
        })
        data = json.loads(_get(f"{MESH_SPARQL}?{params}", timeout=180))
        rows = data["results"]["bindings"]
        for b in rows:
            found[b["d"]["value"].rsplit("/", 1)[-1]] = b["label"]["value"]
        if len(rows) < page:
            break
        offset += page
        time.sleep(0.3)

    dest.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# MeSH tree C04 (Neoplasms) topical descriptors -- the ATLAS cancer definition.\n"
        "# Source: NLM MeSH SPARQL endpoint, https://id.nlm.nih.gov/mesh/sparql (inference on).\n"
        "# Regenerate: python scripts/atlas_baseline.py --refresh-mesh\n"
        f"# descriptors: {len(found)}\n"
    )
    dest.write_text(header + "\n".join(f"{ui}\t{lab}" for ui, lab in sorted(found.items())) + "\n",
                    encoding="utf-8")
    return found


# --------------------------------------------------------------------------
# Baseline file inventory + download
# --------------------------------------------------------------------------

_FILE_RE = re.compile(rb"pubmed\d+n\d+\.xml\.gz")


def list_baseline_files() -> list[str]:
    html = _get(BASELINE_URL, timeout=180)
    return sorted({m.decode() for m in _FILE_RE.findall(html)})


def download(name: str, raw_dir: Path) -> Path:
    dest = raw_dir / name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    raw_dir.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(_get(BASELINE_URL + name, timeout=600))
    tmp.rename(dest)
    return dest


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def _text(node) -> str:
    """Flattened text of an element including tail-bearing children."""
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip() if node is not None else ""


def parse_articles(path: Path, c04: dict):
    """Yield one record per CANCER article in a baseline file.

    An article is cancer if any of its MeSH descriptor UIs is in the C04 set.
    Articles with no MeSH headings at all (very recent, not yet indexed) cannot
    be classified this way and are skipped -- `--status` reports how many, so
    the blind spot is visible rather than silent.
    """
    with gzip.open(path, "rb") as fh:
        for _event, elem in ET.iterparse(fh, events=("end",)):
            if not elem.tag.endswith("PubmedArticle"):
                continue
            try:
                cit = elem.find("MedlineCitation")
                if cit is None:
                    continue
                mesh_uis, mesh_labels, major = [], [], []
                for mh in cit.findall("./MeshHeadingList/MeshHeading/DescriptorName"):
                    ui = mh.get("UI", "")
                    mesh_uis.append(ui)
                    mesh_labels.append(mh.text or "")
                    if mh.get("MajorTopicYN") == "Y":
                        major.append(ui)
                hits = [u for u in mesh_uis if u in c04]
                if not hits:
                    continue

                art = cit.find("Article")
                pmid_el = cit.find("PMID")
                journal = art.find("./Journal/Title") if art is not None else None
                year = art.find("./Journal/JournalIssue/PubDate/Year") if art is not None else None
                medline_date = (art.find("./Journal/JournalIssue/PubDate/MedlineDate")
                                if art is not None else None)
                yr = None
                if year is not None and (year.text or "").isdigit():
                    yr = int(year.text)
                elif medline_date is not None and medline_date.text:
                    m = re.search(r"\d{4}", medline_date.text)
                    yr = int(m.group()) if m else None

                ids = {a.get("IdType"): (a.text or "")
                       for a in elem.findall("./PubmedData/ArticleIdList/ArticleId")}

                yield {
                    "pmid": pmid_el.text if pmid_el is not None else "",
                    "title": _text(art.find("ArticleTitle")) if art is not None else "",
                    "abstract": _text(art.find("Abstract")) if art is not None else "",
                    "journal": (journal.text or "") if journal is not None else "",
                    "year": yr,
                    "doi": ids.get("doi", ""),
                    "pmcid": ids.get("pmc", ""),
                    "mesh_ui": mesh_uis,
                    "mesh": mesh_labels,
                    "mesh_major": major,
                    "cancer_ui": hits,
                    "pub_types": [
                        (p.text or "") for p in
                        (art.findall("./PublicationTypeList/PublicationType") if art is not None else [])
                    ],
                }
            finally:
                elem.clear()


def count_articles(path: Path) -> tuple[int, int]:
    """(total articles, articles with no MeSH headings) in a baseline file."""
    total = nomesh = 0
    with gzip.open(path, "rb") as fh:
        for _event, elem in ET.iterparse(fh, events=("end",)):
            if elem.tag.endswith("PubmedArticle"):
                total += 1
                if elem.find("./MedlineCitation/MeshHeadingList") is None:
                    nomesh += 1
                elem.clear()
    return total, nomesh


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def load_manifest(root: Path) -> dict:
    p = root / "manifest.json"
    if p.exists():
        return json.loads(p.read_text(encoding="utf-8"))
    return {"files": {}, "source": BASELINE_URL}


def save_manifest(root: Path, man: dict) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "manifest.json").write_text(json.dumps(man, indent=1, sort_keys=True), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--limit", type=int, default=0, help="process at most N baseline files")
    ap.add_argument("--status", action="store_true", help="report progress and exit")
    ap.add_argument("--refresh-mesh", action="store_true", help="re-fetch the C04 descriptor set")
    ap.add_argument("--keep-raw", action="store_true",
                    help="keep downloaded XML (default: delete after parsing to save disk)")
    args = ap.parse_args()

    root = atlas_root()
    man = load_manifest(root)

    if args.status:
        done = [f for f, s in man["files"].items() if s.get("parsed")]
        recs = sum(s.get("cancer", 0) for s in man["files"].values())
        seen = sum(s.get("total", 0) for s in man["files"].values())
        nomesh = sum(s.get("no_mesh", 0) for s in man["files"].values())
        print(f"atlas root      : {root}")
        print(f"files parsed    : {len(done)}")
        print(f"articles seen   : {seen:,}")
        print(f"cancer articles : {recs:,}" + (f"  ({recs/seen:.1%} of seen)" if seen else ""))
        print(f"no MeSH headings: {nomesh:,}" + (f"  ({nomesh/seen:.1%}, unclassifiable)" if seen else ""))
        return

    c04 = fetch_c04_descriptors(root / "mesh" / "c04-descriptors.tsv", force=args.refresh_mesh)
    print(f"cancer definition: {len(c04)} MeSH C04 descriptors")

    files = list_baseline_files()
    print(f"baseline files available: {len(files)}")
    todo = [f for f in files if not man["files"].get(f, {}).get("parsed")]
    if args.limit:
        todo = todo[:args.limit]
    print(f"to process this run: {len(todo)}")

    raw_dir, rec_dir = root / "raw", root / "records"
    rec_dir.mkdir(parents=True, exist_ok=True)

    for i, name in enumerate(todo, 1):
        t0 = time.time()
        path = download(name, raw_dir)
        total, nomesh = count_articles(path)
        out = rec_dir / f"{name.replace('.xml.gz', '')}.jsonl.gz"
        n = 0
        with gzip.open(out, "wt", encoding="utf-8") as fh:
            for rec in parse_articles(path, c04):
                fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
                n += 1
        if not args.keep_raw:
            path.unlink(missing_ok=True)
        man["files"][name] = {"parsed": True, "total": total, "cancer": n,
                              "no_mesh": nomesh, "records": out.name}
        save_manifest(root, man)
        print(f"  [{i}/{len(todo)}] {name}: {total:,} articles -> {n:,} cancer "
              f"({n/total:.1%}), {time.time()-t0:.0f}s")

    done = sum(1 for s in man["files"].values() if s.get("parsed"))
    recs = sum(s.get("cancer", 0) for s in man["files"].values())
    print(f"\nparsed {done}/{len(files)} baseline files; {recs:,} cancer articles so far")


if __name__ == "__main__":
    main()
