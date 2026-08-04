#!/usr/bin/env python3
"""Atlas: recover the cancer literature MeSH has not indexed yet (#ATLAS).

THE BLIND SPOT
--------------
`scripts/atlas_baseline.py` defines cancer as membership in MeSH tree C04. That
is the right definition -- a controlled vocabulary maintained by NLM indexers,
not a keyword list -- but it can only classify articles that HAVE been indexed,
and indexing lags publication.

Measured across the baseline, the unindexed share trends sharply upward with
recency (not monotonically -- individual files fluctuate):

    file 0001 (1975)   0.0% unindexed
    file 0409          5.2%
    file 0613         12.6%
    file 0817         18.4%
    file 0959         37.6%

The newest files are the worst. A pure MeSH census therefore loses the most
recent literature -- precisely the part that matters most for spotting what is
emerging. This recovers it.

WHAT IT DOES
------------
Re-scans the baseline for articles with NO MeSH headings at all, and keeps those
whose title or abstract matches a high-precision cancer term set. Output is
written to a SEPARATE stream so a text-matched record can never be mistaken for
a MeSH-indexed one:

    <root>/records_unindexed/<file>.jsonl.gz    (source: "text-match")

VALIDATION
----------
The text matcher's own accuracy is measurable, because the MeSH-indexed portion
of the corpus is ground truth for it. `--validate` runs the matcher over indexed
articles and reports precision and recall against the C04 label. That number is
written into the manifest and must be quoted alongside any count derived from
this layer.

Usage:
    python scripts/atlas_unindexed.py --validate        # measure the matcher first
    python scripts/atlas_unindexed.py --limit 3         # smoke test
    python scripts/atlas_unindexed.py                   # full recovery pass
"""

import argparse
import gzip
import json
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import xml.etree.ElementTree as ET  # noqa: E402

from atlas_baseline import (  # noqa: E402
    atlas_root, download, fetch_c04_descriptors, list_baseline_files, _text,
)

# High-precision cancer vocabulary. Deliberately narrow: these terms are almost
# never used outside oncology in a title or abstract. Broad words that would
# raise recall at a real precision cost ("growth", "mass", "lesion", "benign")
# are excluded on purpose -- this layer's job is to recover what MeSH has not
# reached yet, not to redefine what cancer means.
CANCER_TERMS = re.compile(
    r"\b("
    r"cancers?|carcinomas?|neoplas(m|ms|tic|ia)|tumou?rs?|oncolog(y|ic|ical)"
    r"|malignan(t|cy|cies)|metasta(sis|ses|tic|size)"
    r"|leukemi(a|as)|leukaemi(a|as)|lymphomas?|melanomas?|sarcomas?|myelomas?"
    r"|gliomas?|glioblastomas?|blastomas?|adenocarcinomas?|mesotheliomas?"
    r"|chemotherap(y|eutic)|radiotherap(y|eutic)|immunotherap(y|eutic)"
    r"|antineoplastic|antitumou?r|anti-tumou?r|cytotoxic therapy"
    r"|tumou?rigenes(is|es)|carcinogenes(is|es)|oncogenes?|tumou?r suppressor"
    r")\b",
    re.I,
)


def is_cancer_text(title: str, abstract: str) -> bool:
    return bool(CANCER_TERMS.search(f"{title} {abstract}"))


def scan(path: Path, c04: dict, want_unindexed: bool):
    """Yield (record, has_mesh, is_c04) for every article in a baseline file."""
    with gzip.open(path, "rb") as fh:
        for _ev, elem in ET.iterparse(fh, events=("end",)):
            if not elem.tag.endswith("PubmedArticle"):
                continue
            try:
                cit = elem.find("MedlineCitation")
                if cit is None:
                    continue
                descs = cit.findall("./MeshHeadingList/MeshHeading/DescriptorName")
                has_mesh = bool(descs)
                if want_unindexed and has_mesh:
                    continue
                uis = [d.get("UI", "") for d in descs]
                is_c04 = any(u in c04 for u in uis)

                art = cit.find("Article")
                pmid_el = cit.find("PMID")
                title = _text(art.find("ArticleTitle")) if art is not None else ""
                abstract = _text(art.find("Abstract")) if art is not None else ""
                year_el = art.find("./Journal/JournalIssue/PubDate/Year") if art is not None else None
                yr = int(year_el.text) if year_el is not None and (year_el.text or "").isdigit() else None
                journal = art.find("./Journal/Title") if art is not None else None
                ids = {a.get("IdType"): (a.text or "")
                       for a in elem.findall("./PubmedData/ArticleIdList/ArticleId")}
                yield ({
                    "pmid": pmid_el.text if pmid_el is not None else "",
                    "title": title,
                    "abstract": abstract,
                    "journal": (journal.text or "") if journal is not None else "",
                    "year": yr,
                    "doi": ids.get("doi", ""),
                    "pmcid": ids.get("pmc", ""),
                    "source": "text-match",
                }, has_mesh, is_c04)
            finally:
                elem.clear()


def validate(files, c04, root: Path, n_files: int = 6) -> dict:
    """Measure the text matcher against MeSH truth on INDEXED articles."""
    tp = fp = fn = tn = 0
    raw = root / "raw"
    for name in files[:n_files]:
        p = download(name, raw)
        for rec, has_mesh, is_c04 in scan(p, c04, want_unindexed=False):
            if not has_mesh:
                continue
            pred = is_cancer_text(rec["title"], rec["abstract"])
            if is_c04 and pred:
                tp += 1
            elif is_c04 and not pred:
                fn += 1
            elif not is_c04 and pred:
                fp += 1
            else:
                tn += 1
        p.unlink(missing_ok=True)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec_ = tp / (tp + fn) if tp + fn else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "precision": prec, "recall": rec_,
            "f1": (2 * prec * rec_ / (prec + rec_)) if prec + rec_ else 0.0,
            "files_used": n_files}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--validate", action="store_true",
                    help="measure the text matcher against MeSH truth and exit")
    ap.add_argument("--validate-files", type=int, default=6)
    args = ap.parse_args()

    root = atlas_root()
    c04 = fetch_c04_descriptors(root / "mesh" / "c04-descriptors.tsv")
    files = list_baseline_files()

    man_path = root / "unindexed-manifest.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"files": {}}

    if args.validate:
        # Validate on the RECENT files, where the unindexed problem actually is.
        v = validate(files[-args.validate_files:], c04, root, args.validate_files)
        man["validation"] = v
        man_path.parent.mkdir(parents=True, exist_ok=True)
        man_path.write_text(json.dumps(man, indent=1, sort_keys=True), encoding="utf-8")
        print("text-matcher accuracy vs MeSH C04 on indexed articles "
              f"({v['files_used']} recent baseline files):")
        print(f"  precision {v['precision']:.1%}   recall {v['recall']:.1%}   F1 {v['f1']:.3f}")
        print(f"  tp={v['tp']:,} fp={v['fp']:,} fn={v['fn']:,} tn={v['tn']:,}")
        return

    out_dir = root / "records_unindexed"
    out_dir.mkdir(parents=True, exist_ok=True)
    todo = [f for f in files if not man["files"].get(f, {}).get("done")]
    if args.limit:
        todo = todo[:args.limit]
    print(f"recovery pass over {len(todo)} baseline files")

    raw = root / "raw"
    for i, name in enumerate(todo, 1):
        t0 = time.time()
        p = download(name, raw)
        out = out_dir / f"{name.replace('.xml.gz', '')}.jsonl.gz"
        n_un = n_keep = 0
        with gzip.open(out, "wt", encoding="utf-8") as fh:
            for rec, _hm, _c in scan(p, c04, want_unindexed=True):
                n_un += 1
                if is_cancer_text(rec["title"], rec["abstract"]):
                    fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
                    n_keep += 1
        p.unlink(missing_ok=True)
        man["files"][name] = {"done": True, "unindexed": n_un, "cancer_text": n_keep}
        man_path.write_text(json.dumps(man, indent=1, sort_keys=True), encoding="utf-8")
        print(f"  [{i}/{len(todo)}] {name}: {n_un:,} unindexed -> {n_keep:,} cancer "
              f"by text, {time.time()-t0:.0f}s", flush=True)

    tot = sum(v.get("cancer_text", 0) for v in man["files"].values())
    print(f"\nrecovered {tot:,} text-matched cancer articles MeSH has not indexed")


if __name__ == "__main__":
    main()
