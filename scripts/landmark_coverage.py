#!/usr/bin/env python3
"""Landmark-coverage guardrail + recovery for the local corpus (#345).

Field-defining papers missing from the local corpus distort mechanism-maturity
claims (e.g. radioligand therapy looks immature partly because the VISION trial
is absent). This script:

  1. holds a CURATED registry of landmark / trial-defining papers per mechanism
     (every PMID verified against PubMed by hand before adding here);
  2. checks each landmark's membership in the corpus
     (`corpus/by-pmid/<pmid>.md` = full text, `corpus/abstracts/by-pmid/<pmid>.md`
     = abstract-only, or absent);
  3. writes a per-mechanism coverage report
     (`analysis/landmark-coverage-report.md`);
  4. with `--recover-missing`, fetches any absent landmark's metadata + abstract
     from NCBI E-utilities (+ iCite for the citation metrics) and writes it as a
     tagged abstract record in the corpus format, so the field-definer is no
     longer silently missing.

Honesty: recovered records are ABSTRACT-ONLY. They do NOT change the frozen
full-text quantitative results in the manuscript (which analyze the 4,830
full-text records); they fix the *coverage guardrail* so absence claims are not
corpus artifacts. See `analysis/landmark-corpus-gaps.md`.

Usage:
  python3 scripts/landmark_coverage.py                 # report only
  python3 scripts/landmark_coverage.py --recover-missing  # also fetch+write missing
"""

import argparse
import json
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FULLTEXT_DIR = REPO / "corpus" / "by-pmid"
ABSTRACT_DIR = REPO / "corpus" / "abstracts" / "by-pmid"
REPORT = REPO / "analysis" / "landmark-coverage-report.md"

# Curated registry: every PMID hand-verified against PubMed. Each entry carries
# the tags the corpus analysis reads (mechanism / cancer type / evidence level)
# so a recovered record is consistent with the pipeline's tagging.
LANDMARKS = [
    dict(pmid="34161051", mechanism="radioligand-therapy", cancer="prostate",
         evidence="phase3-clinical", note="VISION: 177Lu-PSMA-617 for mCRPC (NEJM 2021)"),
    dict(pmid="28076709", mechanism="radioligand-therapy", cancer="neuroendocrine",
         evidence="phase3-clinical", note="NETTER-1: 177Lu-DOTATATE for midgut NET (NEJM 2017) "
         "— the SSTR/neuroendocrine radioligand pivotal, complementing the PSMA/prostate "
         "VISION landmark so the registry covers both deployed radioligand targets (#536)"),
    dict(pmid="40448572", mechanism="ttfields", cancer="pancreatic",
         evidence="phase3-clinical", note="PANOVA-3: TTFields + gem/nab-paclitaxel (JCO 2025)"),
    dict(pmid="33016924", mechanism="mRNA-vaccine", cancer="gastrointestinal",
         evidence="phase2-clinical", note="mRNA neoantigen vaccine, GI cancers (JCI 2020)"),
    dict(pmid="36027916", mechanism="mRNA-vaccine", cancer="lung",
         evidence="phase2-clinical", note="NEO-PV-01 + chemo + anti-PD-1, NSCLC (Cancer Cell 2022)"),
    dict(pmid="35970920", mechanism="mRNA-vaccine", cancer="solid-tumor",
         evidence="phase2-clinical", note="ChAd/samRNA individualized neoantigen vaccine (Nat Med 2022)"),
]


def coverage(pmid):
    if (FULLTEXT_DIR / f"{pmid}.md").exists():
        return "full-text"
    if (ABSTRACT_DIR / f"{pmid}.md").exists():
        return "abstract-only"
    return "MISSING"


def _get(url):
    return urllib.request.urlopen(url, timeout=30).read()


def fetch_pubmed(pmid):
    """efetch metadata + abstract; esummary-free (efetch carries it all)."""
    xml = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml")
    art = ET.fromstring(xml).find(".//PubmedArticle")
    if art is None:
        raise RuntimeError(f"PMID {pmid} not found at efetch")
    # itertext() (not findtext) so titles with inline markup (e.g. the
    # <sup>177</sup>Lu in NETTER-1) are captured in full, not truncated at the
    # first child element.
    _title_el = art.find(".//ArticleTitle")
    title = "".join(_title_el.itertext()).strip() if _title_el is not None else ""
    # itertext() (not .text) so structured abstracts with inline formula elements
    # (e.g. <sup>177</sup>Lu) are captured in full, not truncated at the first child.
    abstract = " ".join("".join(t.itertext()).strip()
                        for t in art.findall(".//Abstract/AbstractText")).strip()
    journal = art.findtext(".//Journal/Title") or ""
    year = art.findtext(".//JournalIssue/PubDate/Year") or ""
    month = art.findtext(".//JournalIssue/PubDate/Month") or ""
    vol = art.findtext(".//JournalIssue/Volume") or ""
    issue = art.findtext(".//JournalIssue/Issue") or ""
    pages = art.findtext(".//Pagination/MedlinePgn") or ""
    mesh = [m.findtext("DescriptorName") for m in art.findall(".//MeshHeading") if m.findtext("DescriptorName")]
    ptypes = [p.text for p in art.findall(".//PublicationType") if p.text]
    authors = [f"{a.findtext('LastName')} {a.findtext('Initials')}".strip()
               for a in art.findall(".//Author") if a.findtext("LastName")]
    doi = ""
    for idn in art.findall(".//ELocationID"):
        if idn.get("EIdType") == "doi":
            doi = idn.text
    return dict(title=title, abstract=abstract, journal=journal, year=year, month=month,
                volume=vol, issue=issue, pages=pages, mesh=mesh, ptypes=ptypes,
                authors=authors, doi=doi)


def fetch_icite(pmid):
    try:
        d = json.loads(_get(f"https://icite.od.nih.gov/api/pubs?pmids={pmid}"))["data"][0]
        return dict(rcr=d.get("relative_citation_ratio"), citations=d.get("citation_count"),
                    is_clinical=d.get("is_clinical"))
    except Exception:
        return dict(rcr=None, citations=None, is_clinical=None)


def _yaml_list(items):
    return "".join(f"- {i}\n" for i in items) if items else "[]\n"


def write_abstract_record(entry):
    """Fetch a landmark and write it as a tagged abstract record (corpus format)."""
    pmid = entry["pmid"]
    pm = fetch_pubmed(pmid)
    ic = fetch_icite(pmid)
    fm = []
    fm.append(f"pmid: '{pmid}'")
    fm.append(f"doi: {pm['doi']}")
    fm.append("pmcid: ''")
    fm.append(f"title: {json.dumps(pm['title'])}")
    fm.append("authors:")
    fm += [f"- {a}" for a in pm["authors"]]
    fm.append(f"journal: {json.dumps(pm['journal'])}")
    fm.append(f"year: {pm['year'] or 'null'}")
    fm.append("mesh_terms:")
    fm += [f"- {m}" for m in pm["mesh"]]
    fm.append("pub_types:")
    fm += [f"- {p}" for p in pm["ptypes"]]
    fm.append("mechanisms:")
    fm.append(f"- {entry['mechanism']}")
    fm.append("cancer_types:")
    fm.append(f"- {entry['cancer']}")
    fm.append(f"evidence_level: {entry['evidence']}")
    fm.append(f"icite_rcr: {ic['rcr'] if ic['rcr'] is not None else 'null'}")
    # iCite citation count goes in icite_citation_count (matching existing records);
    # cited_by_count is the OpenAlex field, which this offline recovery does not query.
    fm.append(f"icite_citation_count: {ic['citations'] if ic['citations'] is not None else 'null'}")
    fm.append("cited_by_count: null")
    fm.append(f"icite_is_clinical: {str(bool(ic['is_clinical'])).lower() if ic['is_clinical'] is not None else 'null'}")
    fm.append("source: landmark-recovery (#345)")
    body = pm["abstract"] or "(no abstract available from PubMed)"
    text = "---\n" + "\n".join(fm) + "\n---\n\n" + body + "\n"
    ABSTRACT_DIR.mkdir(parents=True, exist_ok=True)
    (ABSTRACT_DIR / f"{pmid}.md").write_text(text)
    return pm["title"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recover-missing", action="store_true",
                    help="fetch + write abstract records for any MISSING landmark")
    args = ap.parse_args()

    rows = []
    for e in LANDMARKS:
        rows.append((e, coverage(e["pmid"])))

    missing = [e for e, c in rows if c == "MISSING"]
    if args.recover_missing and missing:
        for e in missing:
            title = write_abstract_record(e)
            print(f"recovered {e['pmid']} ({e['mechanism']}): {title[:70]}")
            time.sleep(1)
        rows = [(e, coverage(e["pmid"])) for e in LANDMARKS]  # re-check

    # Per-mechanism report
    by_mech = {}
    for e, c in rows:
        by_mech.setdefault(e["mechanism"], []).append((e, c))
    lines = ["# Landmark coverage report (#345)", "",
             "Generated by `scripts/landmark_coverage.py`. Reusable guardrail: flags",
             "field-defining papers absent from the local corpus so mechanism-level",
             "absence claims are not corpus artifacts. Recovered records are",
             "ABSTRACT-ONLY and do NOT change the frozen full-text quantitative results.",
             "", "| Mechanism | PMID | Coverage | Landmark |", "|---|---|---|---|"]
    for mech in sorted(by_mech):
        for e, c in by_mech[mech]:
            lines.append(f"| {mech} | {e['pmid']} | {c} | {e['note']} |")
    n_missing = sum(1 for _, c in rows if c == "MISSING")
    lines += ["", f"**{len(rows)} landmarks tracked; {n_missing} still MISSING.**",
              "Run `python3 scripts/landmark_coverage.py --recover-missing` to fetch absent ones."]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT}")
    for e, c in rows:
        print(f"  {c:14s} {e['pmid']}  {e['mechanism']}")
    return 1 if n_missing and not args.recover_missing else 0


if __name__ == "__main__":
    sys.exit(main())
