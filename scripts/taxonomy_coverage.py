#!/usr/bin/env python3
"""Expanded-taxonomy coverage report (#347).

The manuscript's quantitative results use a FROZEN 19-mechanism corpus. The
working taxonomy in `scripts/config.py::MECHANISM_KEYWORDS` has since grown a
next-pass scaffold of under-covered mechanisms (radioligand therapy, targeted
protein degradation, oncolytic virus, mRNA vaccine, TTFields, bispecifics, cold
atmospheric plasma, and now cuproptosis + disulfidptosis). This script versions
that expansion concretely: for every mechanism it reports the PubMed-wide
literature size (E-utilities esearch, cancer-scoped) and the local-corpus
coverage, flagging the scaffold mechanisms added beyond the frozen 19.

NON-RETROACTIVITY: this is a SCAFFOLD report. Adding mechanisms to the working
taxonomy does NOT re-tag the frozen corpus or change any quantitative result in
the manuscript (the corpus `INDEX.jsonl` is not regenerated here). It documents
which expanded mechanisms exist and how much literature they have, so the
manuscript's "next-pass scaffold" paragraph (Section 3.2) links to real numbers.

Usage:
  python3 scripts/taxonomy_coverage.py            # offline corpus coverage only
  python3 scripts/taxonomy_coverage.py --pubmed   # also query PubMed counts (network)
"""

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from config import MECHANISM_KEYWORDS  # noqa: E402

INDEX = REPO / "corpus" / "INDEX.jsonl"
ABSTRACT_DIR = REPO / "corpus" / "abstracts" / "by-pmid"
REPORT = REPO / "analysis" / "taxonomy-coverage-report.md"

# The 19 mechanisms the manuscript's quantitative results are frozen on. Anything
# in MECHANISM_KEYWORDS but NOT here is a next-pass scaffold mechanism.
FROZEN_19 = {
    "ttfields", "immunotherapy", "car-t", "crispr", "nanoparticle",
    "metabolic-targeting", "oncolytic-virus", "mRNA-vaccine", "synthetic-lethality",
    "bioelectric", "electrolysis", "sonodynamic", "cold-atmospheric-plasma", "hifu",
    "electrochemical-therapy", "epigenetic", "microbiome", "frequency-therapy",
    "antibody-drug-conjugate",
}
# Mechanisms first added by THIS issue (#347).
NEW_THIS_PR = {"cuproptosis", "disulfidptosis"}


def local_tagged_count(mech):
    """Full-text records already TAGGED with this mechanism in INDEX.jsonl."""
    n = 0
    if not INDEX.exists():
        return 0
    for line in INDEX.read_text().splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if mech in (rec.get("mechanisms") or []):
            n += 1
    return n


def local_keyword_hits(keywords):
    """Records (full-text titles via INDEX + abstract files) whose text contains a
    keyword. A read-only preview of how many corpus records the mechanism WOULD
    cover if the corpus were re-tagged (it is NOT)."""
    kws = [k.lower() for k in keywords]
    hits = 0
    # full-text titles (cheap proxy from the index)
    if INDEX.exists():
        for line in INDEX.read_text().splitlines():
            if not line.strip():
                continue
            try:
                t = (json.loads(line).get("title") or "").lower()
            except json.JSONDecodeError:
                continue
            if any(k in t for k in kws):
                hits += 1
    # abstract records (title + body)
    for f in ABSTRACT_DIR.glob("*.md"):
        txt = f.read_text(errors="ignore").lower()
        if any(k in txt for k in kws):
            hits += 1
    return hits


def pubmed_count(keywords):
    """Cancer-scoped PubMed count for the OR'd seed terms (esearch)."""
    terms = " OR ".join(f'"{k}"' for k in keywords)
    q = f"({terms}) AND cancer"
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed"
           f"&term={urllib.parse.quote(q)}&retmode=json&retmax=0")
    try:
        d = json.loads(urllib.request.urlopen(url, timeout=30).read())
        return int(d["esearchresult"]["count"])
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pubmed", action="store_true", help="also query PubMed counts (network)")
    args = ap.parse_args()

    rows = []
    for mech, kws in sorted(MECHANISM_KEYWORDS.items()):
        scaffold = mech not in FROZEN_19
        tagged = local_tagged_count(mech)
        hits = local_keyword_hits(kws) if scaffold else tagged
        pm = pubmed_count(kws) if args.pubmed else None
        if args.pubmed:
            time.sleep(0.4)
        rows.append((mech, scaffold, mech in NEW_THIS_PR, tagged, hits, pm))
        print(f"  {mech:26s} {'scaffold' if scaffold else 'frozen-19':10s} "
              f"tagged={tagged:4d} text-hits={hits:4d} pubmed={pm}")

    n_total = len(MECHANISM_KEYWORDS)
    n_scaffold = sum(1 for _, s, *_ in rows if s)
    lines = [
        "# Expanded-taxonomy coverage report (#347)", "",
        "Generated by `scripts/taxonomy_coverage.py`.",
        "",
        f"The working taxonomy now has **{n_total} mechanisms** ({n_total - n_scaffold} "
        f"frozen + {n_scaffold} next-pass scaffold). **This is a scaffold report: it does "
        "NOT re-tag the frozen corpus or change any manuscript quantitative result** "
        "(the corpus `INDEX.jsonl` is not regenerated here). `tagged` = records already "
        "tagged with the mechanism in the frozen index; `text-hits` (scaffold rows) = "
        "records whose title/abstract text matches the seeds, a preview of coverage IF the "
        "corpus were re-tagged.",
        "",
        "| Mechanism | Tier | Tagged (frozen) | Text hits (preview) | PubMed (cancer) |",
        "|---|---|---|---|---|",
    ]
    for mech, scaffold, new, tagged, hits, pm in rows:
        tier = "**new (#347)**" if new else ("scaffold" if scaffold else "frozen-19")
        pm_s = "n/a" if pm is None else str(pm)
        lines.append(f"| {mech} | {tier} | {tagged} | {hits if scaffold else '-'} | {pm_s} |")
    lines += [
        "",
        "**cuproptosis** and **disulfidptosis** (the emerging regulated-cell-death "
        "modalities named in manuscript Section 3.2) are added here. The corpus already "
        "contains matching records (see Text hits); they are NOT mechanism-tagged in the "
        "frozen index, so adding the mechanisms makes the taxonomy current without altering "
        "the frozen 19-mechanism results.",
    ]
    # Only write the committed report when the PubMed column is real. A default
    # offline run (no --pubmed) would fill every PubMed cell with n/a and silently
    # clobber the checked-in report's verified counts, so in that mode we print the
    # coverage to stdout and leave the file untouched.
    if args.pubmed:
        REPORT.write_text("\n".join(lines) + "\n")
        print(f"wrote {REPORT}")
    else:
        print(f"offline mode: {REPORT} left unchanged (re-run with --pubmed to refresh it)")


if __name__ == "__main__":
    main()
