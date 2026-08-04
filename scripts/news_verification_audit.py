#!/usr/bin/env python3
"""Do the news pipeline's "verified" PMIDs support the claims? (#NEWS-VERIFY)

WHY
---
`news/NEWS_INDEX.jsonl` marked 44 claims `verification_status: "verified"` and
attached `linked_pmids` as the evidence. Two of those links reached the
manuscript as footnotes reading *"Verified: PMID:X"*.

They did not survive being looked at. The claim about electric fields and brain
cancer was "verified" against five papers on freshwater fish biodiversity,
school vision programs, printable colouration, post-vaccine inflammation, and a
survey of speech-language pathologists. Not one concerned electric fields, brain
tumours, or immunotherapy.

The pattern said what happened: the linked identifiers clustered in a narrow
numeric band (66 of 175 shared the prefix `4202`), which is what you get from
PubMed identifiers assigned in the same indexing window. The linking was picking
up records adjacent in time, not in content. That measurement -- 55.7% of 212
resolved pairs sharing no content word -- is kept below as BASELINE so the
delta stays readable.

`verify_news_claims.py` has since been fixed (sentence-initial capitals no
longer read as proper nouns; a candidate must share >=2 content words with the
claim; a claim with no searchable term is marked unverified instead of keeping
a stale verdict; withdrawn evidence is cleared rather than left attached) and
re-run over every article. This script is what says whether that worked.

WHAT IT MEASURES
----------------
Content-word overlap between the news article and each linked paper's title. It
is a deliberately crude test -- a real supporting citation shares vocabulary with
the claim it supports, and a citation sharing NOTHING is not evidence of
anything.

Crude cuts both ways, so the result is reported as a band rather than a verdict:
overlap can miss a genuine supporting paper that uses different words, and can
credit an unrelated one that happens to share a common term. It measures the
FLOOR -- links that cannot be support -- and never establishes that a surviving
link IS support.

Requires network (NCBI E-utilities). Verdicts are committed.

Usage:
    python scripts/news_verification_audit.py
"""

import collections
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

INDEX = PROJECT_ROOT / "news" / "NEWS_INDEX.jsonl"
OUT = PROJECT_ROOT / "analysis" / "news-verification-audit.md"
RAW = PROJECT_ROOT / "analysis" / "news-verification-audit.json"

# Words too common to count as shared subject matter.
STOP = set(
    "the a an of and or in on for with to by at from as is are was were this that "
    "have has had been be can could may might will would new study studies research "
    "using via into over under more most less than then when what which who how why "
    "its it their they them there here about after before during between within "
    "results result found find shows show showed报".split())

# Oncology boilerplate: real words, but every second cancer paper carries them,
# so an overlap made ONLY of these clears the >=2 bar without establishing that
# the paper is about the same thing as the claim. Reported separately so the
# headline share is not read as "62% of the links are support".
GENERIC = set(
    "cancer cancers cancerous tumor tumors tumour tumours patient patients "
    "clinical trial trials phase therapy therapies therapeutic treatment "
    "treatments cell cells survival prognosis prognostic disease diseases "
    "center centre centers multicenter multicentre randomized randomised "
    "efficacy safety analysis outcomes united states europe medicine medical "
    "oncology oncologic".split())

# The pre-fix measurement, kept so the delta is readable from the report itself
# rather than from a commit message. Measured 2026-06 on the same script.
BASELINE = {
    "claims_with_links": 44,
    "distinct_pmids": 175,
    "pairs_resolved": 212,
    "zero_overlap": 118,
    "one_word": 41,
    "two_plus": 53,
    "dominant_prefix": "4202",
    "dominant_prefix_count": 66,
}


def titles(pmids: list) -> dict:
    out = {}
    for i in range(0, len(pmids), 150):
        batch = pmids[i:i + 150]
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
               "?db=pubmed&retmode=json&id=" + ",".join(batch))
        try:
            res = json.load(urllib.request.urlopen(url, timeout=90))["result"]
        except Exception as exc:
            print(f"  ! esummary failed: {exc}", file=sys.stderr)
            continue
        for p in batch:
            if p in res:
                out[p] = res[p].get("title", "")
        time.sleep(0.4)
    return out


def content_words(s: str) -> set:
    return {w for w in re.findall(r"[a-z]{4,}", (s or "").lower()) if w not in STOP}


def main() -> int:
    if not INDEX.exists():
        print(f"missing {INDEX}", file=sys.stderr)
        return 1
    rows = [json.loads(l) for l in INDEX.read_text().splitlines() if l.strip()]
    linked = [r for r in rows if r.get("linked_pmids")]
    pmids = sorted({p for r in linked for p in r["linked_pmids"]})
    print(f"{len(rows)} claims, {len(linked)} with linked PMIDs, "
          f"{len(pmids)} distinct identifiers", flush=True)
    tmap = titles(pmids)

    pairs = []
    for r in linked:
        news = content_words(r.get("article_title", "")) | content_words(r.get("claim_text", ""))
        for p in r["linked_pmids"]:
            t = tmap.get(p)
            if t is None:
                continue
            shared = news & content_words(t)
            pairs.append({"pmid": p, "overlap": len(shared),
                          "shared": sorted(shared)[:8],
                          "generic_only": bool(shared) and shared <= GENERIC,
                          "article": r.get("article_title", "")[:70],
                          "paper": t[:70],
                          "claim_id": r.get("claim_id", "")})
    zero = [x for x in pairs if x["overlap"] == 0]
    weak = [x for x in pairs if x["overlap"] == 1]
    ok = [x for x in pairs if x["overlap"] >= 2]
    generic = [x for x in ok if x["generic_only"]]

    prefix = collections.Counter(p[:4] for p in pmids)
    top_prefix, top_n = prefix.most_common(1)[0]

    def pct(n, d):
        return 100 * n / max(1, d)

    b = BASELINE
    withdrawn = b["claims_with_links"] - len(linked)
    L = [
        "# Do the news pipeline's \"verified\" PMIDs support the claims? (#NEWS-VERIFY)", "",
        "Generated by `scripts/news_verification_audit.py`.", "",
        f"`news/NEWS_INDEX.jsonl` now marks **{len(linked)}** claims `verified` with",
        f"`linked_pmids` attached as the evidence, against **{b['claims_with_links']}** before",
        "`verify_news_claims.py` was fixed and re-run.", "",
        "## The measurement", "",
        "Content-word overlap between each news article and each linked paper's title.",
        "A real supporting citation shares vocabulary with the claim it supports; one",
        "sharing nothing is not evidence of anything. This is a FLOOR test: it says",
        "which links cannot be support, never that a surviving link is support.", "",
        "| overlap with the news article | before | share | after | share |",
        "|---|---|---|---|---|",
        f"| **none** | {b['zero_overlap']:,} | **{pct(b['zero_overlap'], b['pairs_resolved']):.1f}%** "
        f"| {len(zero):,} | **{pct(len(zero), len(pairs)):.1f}%** |",
        f"| one word | {b['one_word']:,} | {pct(b['one_word'], b['pairs_resolved']):.1f}% "
        f"| {len(weak):,} | {pct(len(weak), len(pairs)):.1f}% |",
        f"| two or more | {b['two_plus']:,} | {pct(b['two_plus'], b['pairs_resolved']):.1f}% "
        f"| {len(ok):,} | {pct(len(ok), len(pairs)):.1f}% |",
        f"| total resolved pairs | {b['pairs_resolved']:,} | | {len(pairs):,} | |", "",
        "## What the first measurement found", "",
        f"The linked identifiers clustered in a narrow numeric band -- {b['dominant_prefix_count']} of",
        f"{b['distinct_pmids']} shared the prefix `{b['dominant_prefix']}`. PubMed identifiers are",
        "assigned sequentially as records are indexed, so that is the signature of",
        "matching on **when a record was indexed rather than what it says**.", "",
        "The clearest case: a claim about electric fields treating brain cancer, linked",
        "to papers on freshwater fish biodiversity, school-based vision programs,",
        "printable colouration, post-vaccine inflammation tracking, and a survey of",
        "speech-language pathologists.", "",
        "Two root causes, both now fixed in `verify_news_claims.py`: a sentence-initial",
        "capital was read as a proper noun (\"Seven of these 26 patients\" -> the query",
        "`Seven`, 835,973 records, five newest returned), and ANY non-empty search",
        "result was accepted as verification.", "",
        "## What the re-run changed, and what it did not", "",
        f"The no-overlap share falls from {pct(b['zero_overlap'], b['pairs_resolved']):.1f}% to "
        f"{pct(len(zero), len(pairs)):.1f}%, and the index-adjacency signature is gone",
        f"({top_n} of {len(pmids)} identifiers now share the leading prefix `{top_prefix}`,",
        f"against {b['dominant_prefix_count']}/{b['distinct_pmids']} before).", "",
        "**Read that as withdrawal, not repair.** The denominator fell with it: "
        f"{b['pairs_resolved']:,} pairs to {len(pairs):,}, because {withdrawn} of the "
        f"{b['claims_with_links']} verifications were",
        "dropped outright and their evidence cleared. The bad links were removed; they",
        "were not replaced with good ones.", "",
        "## The residual: two shared words is a low bar", "",
        f"Of the {len(ok):,} pairs clearing the two-word bar, **{len(generic):,}** clear it on",
        "oncology boilerplate alone -- words like *cancer*, *trial*, *phase*, *patients*,",
        "*united states* -- which every second cancer paper carries. Those pairs are not",
        "distinguishable from coincidence by this measurement.", "",
        "| claim | linked paper | shared words |", "|---|---|---|",
    ]
    for x in generic[:12]:
        L.append(f"| {x['article']} | {x['paper']} | {', '.join(x['shared'])} |")

    L += [
        "", "The other unfixed asymmetry is that `supports_claim()` gates only the PubMed",
        "fallback. The local-corpus path accepts a hit whenever two extracted terms",
        "appear as substrings of a corpus title, and it runs FIRST, so a weak corpus",
        "match short-circuits the search before PubMed is ever asked.", "",
    ]

    if zero:
        L += ["## Remaining pairs with no shared content word", "",
              "| news article | linked paper |", "|---|---|"]
        for x in zero[:15]:
            L.append(f"| {x['article']} | {x['paper']} |")
        L.append("")

    L += [
        "## What this does not say", "",
        "Word overlap is a crude test and it cuts both ways. It can miss a genuine",
        "supporting paper that describes the same finding in different vocabulary, and",
        "it can credit an unrelated paper that happens to share a common term. So the",
        "figure above is a band, not a verdict on any individual claim.", "",
        "In particular it does NOT establish that the surviving `verified` claims are",
        "verified. It establishes that the specific failure it was built to detect --",
        "links chosen by indexing date -- is no longer the dominant pattern.", "",
        "## What should happen", "",
        "* A `verification_status: \"verified\"` label is now worth reading, but a claim",
        "  quoted as verified anywhere load-bearing -- a manuscript footnote above all --",
        "  should still be checked against the linked paper by hand first.",
        "* The corpus path needs the same content bar as the PubMed path.",
        "* Query construction is still the weak step: multi-word capitalised phrases",
        "  produce queries like `This Phase United States`, which return recently",
        "  indexed records for the same reason the original bug did.",
        "* Re-running the linker remains a separate, reviewable job from this audit,",
        "  which only measures.", "",
        "## Limits", "",
        "* Titles only. A paper whose abstract supports a claim its title does not",
        "  mention counts here as unrelated, so the no-overlap share is an over-estimate",
        "  of the true error rate.",
        "* Only claims that carry `linked_pmids` are examined; claims with none are",
        "  outside this measurement entirely. After the re-run that exclusion is much",
        f"  larger than it was: {withdrawn} claims left the measured set by losing their links.",
        "* The stoplist and the boilerplate list are both hand-written, so the counts",
        "  shift slightly with them.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "claims_total": len(rows), "claims_with_links": len(linked),
        "distinct_pmids": len(pmids), "pairs_resolved": len(pairs),
        "zero_overlap": len(zero), "one_word": len(weak), "two_plus": len(ok),
        "two_plus_generic_only": len(generic),
        "dominant_prefix": top_prefix, "dominant_prefix_count": top_n,
        "baseline": b,
        "examples_zero_overlap": zero[:40],
        "examples_generic_only": generic[:40],
    }, indent=2) + "\n")
    print(f"\nno shared content word: {len(zero)}/{len(pairs)} "
          f"({100*len(zero)/max(1,len(pairs)):.1f}%)")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
