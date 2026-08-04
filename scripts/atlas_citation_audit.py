#!/usr/bin/env python3
"""Atlas: do the module citations point at the papers they claim? (#ATLAS-CITE)

WHY
---
`atlas_module_support.py` exists because each simulation layer rests on "one or
two papers the author read, cited by PMID in the module docs". It measures how
many OTHER articles assert the same relation. It never checks the citation
itself.

That check turned out to matter. Three module citations point at papers with
nothing to do with their claim -- the `dhodh` layer cites a Nature news item
about US fetal-tissue policy, `prom2` cites a Theriogenology paper on embryo
vitrification, and the copper/FIN synergy note cites a study of electric-scooter
fractures. A citation that does not exist in the census reads in the support
table as "cited-absent", which looks like an extraction gap and is not.

WHAT IT CHECKS
--------------
For each claim, whether the cited paper's own title, abstract and MeSH mention
either entity in the claim. A paper can describe a mechanism without naming a
gene, so a miss is a REVIEW FLAG rather than proof of error -- but a paper that
mentions neither entity, in a claim built entirely on it, is worth a human look.

Requires network (NCBI E-utilities). The verdicts are committed so the result is
readable offline.

Usage:
    python scripts/atlas_citation_audit.py
"""

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_module_support import CLAIMS  # noqa: E402

from config import PROJECT_ROOT  # noqa: E402

# A paper can establish a mechanism without naming the gene symbol this project
# picked as its proxy -- Dixon 2012 defines System Xc- inhibition without writing
# "SLC7A11". So a missing entity name alone cannot distinguish a WRONG PAPER from
# a naming mismatch. This vocabulary does: a citation that mentions neither its
# entities nor any term from the field is in the wrong subject area entirely.
DOMAIN_TERMS = re.compile(
    r"ferroptos|lipid peroxidation|glutathione|gpx4|erastin|rsl3|lipoxygenase"
    r"|system xc|cystine|labile iron|ferritin|reactive oxygen|peroxidation"
    r"|antioxidant|redox|tumor|tumour|cancer|carcinoma|neoplas")

OUT = PROJECT_ROOT / "analysis" / "atlas-citation-audit.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-citation-audit.json"


def fetch(pmids: list) -> dict:
    out = {}
    for i in range(0, len(pmids), 100):
        batch = pmids[i:i + 100]
        data = urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}).encode()
        req = urllib.request.Request(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", data=data)
        try:
            xml = urllib.request.urlopen(req, timeout=180).read().decode("utf-8", "replace")
        except Exception as exc:
            print(f"  ! efetch failed: {exc}", file=sys.stderr)
            continue
        for art in xml.split("<PubmedArticle>")[1:]:
            m = re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
            if not m:
                continue
            t = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", art, re.S)
            j = re.search(r"<Title>(.*?)</Title>", art, re.S)
            out[m.group(1)] = {
                "title": re.sub(r"<[^>]+>", "", t.group(1)).strip() if t else "",
                "journal": re.sub(r"<[^>]+>", "", j.group(1)).strip() if j else "",
                "text": re.sub(r"<[^>]+>", " ", art).lower(),
            }
        time.sleep(0.4)
    return out


def prose_citations() -> dict:
    """PMID -> the documents citing it, for prose that carries no entity pair.

    Module claims name two entities, so a citation can be checked against them.
    The manuscript and the repo guide cite PMIDs in running text with no such
    handle, so those are checked only for existing and for being in the right
    subject area at all -- which is enough to catch a citation pointing at
    speech-language pathology or electric-scooter injuries.
    """
    out = {}
    for rel in ("article/drafts/v1.md", "CLAUDE.md"):
        f = PROJECT_ROOT / rel
        if not f.exists():
            continue
        for pm in set(re.findall(r"PMID:?\s*(\d{7,8})", f.read_text())):
            out.setdefault(pm, []).append(rel)
    return out


def main() -> int:
    pmids = sorted({c[3] for c in CLAIMS})
    print(f"checking {len(pmids)} cited PMIDs behind {len(CLAIMS)} claims ...", flush=True)
    recs = fetch(pmids)

    rows = []
    for module, a, b, pmid, claim in CLAIMS:
        rec = recs.get(pmid)
        if rec is None:
            rows.append({"module": module, "pmid": pmid, "claim": claim,
                         "status": "unresolvable", "title": "", "journal": "",
                         "mentions": []})
            continue
        text = rec["text"]
        hits = [e for e in (a, b) if e.lower() in text]
        if hits:
            status = "ok"
        elif DOMAIN_TERMS.search(text):
            status = "entity-not-named"
        else:
            status = "wrong-subject"
        rows.append({"module": module, "pmid": pmid, "claim": claim,
                     "status": status, "title": rec["title"],
                     "journal": rec["journal"], "mentions": hits})

    wrong = [r for r in rows if r["status"] in ("wrong-subject", "unresolvable")]
    soft = [r for r in rows if r["status"] == "entity-not-named"]
    bad = wrong + soft

    # The same existence-and-subject check over prose citations.
    prose = prose_citations()
    prose_only = sorted(set(prose) - set(pmids))
    print(f"checking {len(prose_only)} further PMIDs cited in prose ...", flush=True)
    precs = fetch(prose_only)
    prose_rows = []
    for pm in prose_only:
        rec = precs.get(pm)
        if rec is None:
            prose_rows.append({"pmid": pm, "cited_in": prose[pm], "status": "unresolvable",
                               "title": "", "journal": ""})
            continue
        st = "ok" if DOMAIN_TERMS.search(rec["text"]) else "off-topic"
        prose_rows.append({"pmid": pm, "cited_in": prose[pm], "status": st,
                           "title": rec["title"], "journal": rec["journal"]})
    prose_bad = [r for r in prose_rows if r["status"] != "ok"]
    L = [
        "# Do the module citations point at the papers they claim? (#ATLAS-CITE)", "",
        "Generated by `scripts/atlas_citation_audit.py`.",
        "`atlas_module_support.py` measures how many OTHER articles assert each",
        "module's relation. It never checked the citation itself, and that turned out",
        "to matter.", "",
        "## Method", "",
        "For each claim, whether the cited paper's own title, abstract and MeSH mention",
        "either entity in the claim. A paper can describe a mechanism without naming a",
        "gene, so a miss is a **review flag, not proof of error** -- but a claim built",
        "entirely on a paper that mentions neither of its entities deserves a look.", "",
        f"## Result: {len(wrong)} broken citations, {len(soft)} naming mismatches", "",
        "The two are different problems and must not be pooled. A paper that mentions",
        "neither entity but is plainly in the field is very likely the right paper",
        "under a different name -- Dixon 2012 defines System Xc- inhibition without",
        "ever writing `SLC7A11`. A paper that mentions neither the entities nor any",
        "term from the field is a different paper altogether.", "",
        "| module | PMID | journal | title | mentions |", "|---|---|---|---|---|",
    ]
    for r in rows:
        mark = "" if r["status"] == "ok" else " **FLAG**"
        L.append(f"| `{r['module']}`{mark} | {r['pmid']} | {r['journal'][:28]} | "
                 f"{r['title'][:64]} | {', '.join(r['mentions']) or '-'} |")

    if wrong:
        L += ["", "## Broken citations -- these point at unrelated papers", ""]
        for r in wrong:
            L += [f"### `{r['module']}` -> PMID {r['pmid']}", "",
                  f"* **claim**: {r['claim']}",
                  f"* **cited paper**: {r['title'] or '(unresolvable)'}"
                  + (f" _{r['journal']}_" if r["journal"] else ""),
                  "* mentions neither of its entities NOR any term from this field"
                  if r["status"] == "wrong-subject"
                  else "* PMID does not resolve at PubMed", ""]
    if soft:
        L += ["## Naming mismatches -- probably the right paper, wrong proxy name", "",
              "These are in the field but do not use the gene symbol this project",
              "chose as the claim's proxy. Reported for completeness, not as errors.", ""]
        for r in soft:
            L += [f"* `{r['module']}` -> {r['pmid']}: {r['title'][:80]}"]
        L += [""]

    if wrong:
        L += [
            "## What this means", "",
            "These are not extraction failures. A citation pointing at an unrelated",
            "paper reads in `atlas-module-support.md` as `cited-absent`, which looks",
            "like the graph failed to find a real paper -- when in fact the paper being",
            "looked for is not the one the module meant.", "",
            "Correcting them needs a human: the right PMID has to be identified from the",
            "mechanism, and guessing one would replace a visible error with an invisible",
            "one. They are listed here rather than patched.", "",
        ]

    L += [f"## Prose citations: {len(prose_bad)} of {len(prose_rows)} flagged", "",
          "The manuscript and repo guide cite PMIDs in running text with no entity",
          "pair to check against, so these are tested only for existing and for being",
          "in the right subject area. That is a low bar, and the point is that a",
          "citation can fail it.", ""]
    if prose_bad:
        L += ["| PMID | cited in | journal | title |", "|---|---|---|---|"]
        for r in sorted(prose_bad, key=lambda r: r["status"]):
            L.append(f"| {r['pmid']} | {', '.join(r['cited_in'])} | "
                     f"{r['journal'][:24]} | {r['title'][:60] or '(unresolvable)'} |")
        L += ["",
              "A subject-area miss is not automatically an error -- a methods or",
              "statistics citation legitimately mentions no cancer term. Each needs a",
              "human look, and they are listed rather than patched.", ""]

    L += [
        "## Limits", "",
        "* Entity mention is a weak proxy for relevance. A paper can establish a",
        "  mechanism without naming the gene this project chose as its proxy, so an",
        "  `ok` here is not confirmation that the citation supports the claim -- only",
        "  that it is in the right subject area.",
        "* Only the 20 module claims are checked, not every PMID in the manuscript.",
        "* Title, abstract and MeSH only. A gene named solely in the body of the paper",
        "  counts as a miss.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({"claims": len(rows), "flagged": len(bad),
                               "rows": rows,
                               "prose_checked": len(prose_rows),
                               "prose_flagged": len(prose_bad),
                               "prose_rows": prose_rows}, indent=2) + "\n")
    print(f"\nprose: {len(prose_bad)} of {len(prose_rows)} flagged")
    for r in prose_bad:
        print(f"  {r['status']:<14}{r['pmid']}  {r['title'][:56] or '(unresolvable)'}")
    print(f"\n{len(bad)} of {len(rows)} claims flagged")
    for r in bad:
        print(f"  {r['module']:<16} {r['pmid']}  {r['title'][:56] or '(unresolvable)'}")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
