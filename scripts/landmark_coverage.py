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

Accounting: recovered records are OUT-OF-BAND BY PROVENANCE from the PRISMA
systematic search. By the scoping protocol (`analysis/prisma-scr-protocol.md`),
the systematic-search totals (10,415 screened -> 4,830 full-text + 5,586
abstract-only) count only records returned by the 19 mechanism queries; landmark
field-definers are curated by hand, NOT returned by those queries. They are
written to `corpus/abstracts/by-pmid/`, so `oa_bias_analysis.py` (which globs that
directory) DOES count them in the live archive total -- but they never touch the
4,830 full-text quantitative corpus, and the §3.3.1 open-access CONCLUSIONS are
unaffected: regenerating `analysis/oa-bias-report.md` after recovery moves only
the raw abstract counts (the archive total and a few mechanism tallies) by the
handful of added records, while every rounded mechanism SHARE and RANK --
including the manuscript-cited 34.4->28.7% / 14.7->22.4% shifts and the
bioelectric 14->3 reordering -- is identical. The manuscript's frozen PRISMA
figure (5,586 abstract-only) is the archive snapshot at manuscript-freeze time
(which folded in the first two recoveries, VISION + NETTER-1).

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
    # CAR-T pivotals (#568): the two field-defining registrational single-arm
    # phase-2 trials behind the first FDA CAR-T approvals. car-t is well
    # represented in the corpus (474 records), so unlike radioligand these do
    # not fix a maturity-distortion — they make the registry's field-definer
    # coverage explicit for the headline cell-therapy mechanism.
    dict(pmid="29226797", mechanism="car-t", cancer="lymphoma",
         evidence="phase2-clinical", note="ZUMA-1: axicabtagene ciloleucel (axi-cel) for refractory "
         "large B-cell lymphoma (NEJM 2017) — the pivotal single-arm phase-2 behind the first DLBCL "
         "CAR-T approval (#568)"),
    dict(pmid="29385370", mechanism="car-t", cancer="leukemia",
         evidence="phase2-clinical", note="ELIANA: tisagenlecleucel for relapsed/refractory pediatric "
         "and young-adult B-cell ALL (NEJM 2018) — the pivotal single-arm phase-2 behind the first "
         "CAR-T approval of any kind (#568)"),
    # ADC pivotals (#568): the trastuzumab-deruxtecan HER2 randomised phase-3s
    # that redefined the antibody-drug-conjugate class (incl. the HER2-low
    # expansion that created a new treatable population).
    dict(pmid="35320644", mechanism="antibody-drug-conjugate", cancer="breast",
         evidence="phase3-clinical", note="DESTINY-Breast03: trastuzumab deruxtecan (T-DXd) vs T-DM1 in "
         "HER2-positive metastatic breast cancer (NEJM 2022) — head-to-head ADC superiority (#568)"),
    dict(pmid="35665782", mechanism="antibody-drug-conjugate", cancer="breast",
         evidence="phase3-clinical", note="DESTINY-Breast04: trastuzumab deruxtecan in previously "
         "treated HER2-LOW metastatic breast cancer (NEJM 2022) — established HER2-low as a new ADC "
         "target population (#568)"),
    # Checkpoint-immunotherapy landmarks (#568): three of the trials that defined
    # the modern checkpoint era (melanoma anti-PD-1 and dual blockade, plus the
    # chemo-IO combination that became first-line NSCLC). immunotherapy is the
    # largest corpus mechanism (2,297 records); these anchor its field-definers.
    dict(pmid="25891173", mechanism="immunotherapy", cancer="melanoma",
         evidence="phase3-clinical", note="KEYNOTE-006: pembrolizumab vs ipilimumab in advanced "
         "melanoma (NEJM 2015) — anti-PD-1 superiority over the prior anti-CTLA-4 standard (#568)"),
    dict(pmid="26027431", mechanism="immunotherapy", cancer="melanoma",
         evidence="phase3-clinical", note="CheckMate-067: nivolumab + ipilimumab or monotherapy in "
         "untreated melanoma (NEJM 2015) — the dual-checkpoint-blockade landmark (#568)"),
    dict(pmid="29658856", mechanism="immunotherapy", cancer="lung",
         evidence="phase3-clinical", note="KEYNOTE-189: pembrolizumab + chemotherapy in metastatic "
         "nonsquamous NSCLC (NEJM 2018) — established first-line chemo-immunotherapy (#568)"),
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
             "They are out-of-band BY PROVENANCE (curated by hand, not returned by the",
             "19 mechanism queries). They are stored in `corpus/abstracts/by-pmid/`, so",
             "the §3.3.1 OA-bias tally counts them in the live archive total, but they",
             "move only raw counts -- every rounded mechanism share and rank is identical.",
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
