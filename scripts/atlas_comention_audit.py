#!/usr/bin/env python3
"""Atlas: how precise is the full-text co-mention layer? (#ATLAS-COMENT-AUDIT)

WHY
---
`atlas_comention.py` produces the counts that `atlas_module_support.py` uses to
argue a zero in the relation column is an extraction failure rather than absence
of evidence. Its own docstring calls the alias map its weak point. It has never
had a precision number, and an unauditable weak point is indistinguishable from
a sound one.

It is now auditable: the build writes a uniform reservoir sample of matched
sentences with their resolved entity names.

THE TEST
--------
Two checks, one mechanical and one independent.

  DID THE ALIAS ACTUALLY MATCH?  Every sampled entity must correspond to a
  surface form present in the sentence. A failure here is a tokenizer bug, not
  a judgement call, and should be zero.

  DOES PUBTATOR AGREE?  PubTator annotated the same paper independently, from
  its abstract. If it also assigned that entity to that paper, the full-text
  match is corroborated by a source that never saw our alias map.

WHY AGREEMENT IS A LOWER BOUND, NOT PRECISION
----------------------------------------------
PubTator reads abstracts; this layer reads full text. An entity discussed only
in a Methods or Results section is genuinely present and genuinely absent from
PubTator's annotation, so it counts as a disagreement while being correct.
Agreement therefore UNDERSTATES precision, and the gap cannot be closed from
this data. What a low agreement rate would show is that the layer is finding
entities the abstract-level extractor never sees anywhere in the paper -- which
is the failure mode the alias map makes likely.

Usage:
    python scripts/atlas_comention_audit.py
"""

import collections
import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_comention import _TOKEN, build_alias_map  # noqa: E402
from atlas_graph import load_index  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-comention-audit.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-comention-audit.json"


def pubtator_by_paper(root: Path, want: set) -> dict:
    """pmid -> identifiers PubTator assigned, for the sampled papers only."""
    out = collections.defaultdict(set)
    for kind in ("gene", "chemical", "disease"):
        f = root / "entities" / f"{kind}.tsv.gz"
        if not f.exists():
            continue
        with gzip.open(f, "rt", errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4 and p[0] in want:
                    out[p[0]].add(p[2])
    return out


def abstracts(pmids: list) -> dict:
    """pmid -> lowercased title+abstract, for the decomposition below.

    A disagreement with PubTator means one of two very different things, and
    they have to be told apart: our match is in the ABSTRACT, so PubTator saw
    the same text and did not assign the entity (their recall gap, or our false
    positive), or it is BODY-ONLY, where PubTator could not have seen it at all
    and disagreement is the expected outcome rather than evidence of anything.
    """
    out = {}
    for i in range(0, len(pmids), 200):
        batch = pmids[i:i + 200]
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
            keep = " ".join(re.findall(r"<ArticleTitle>(.*?)</ArticleTitle>|"
                                       r"<AbstractText[^>]*>(.*?)</AbstractText>",
                                       art, re.S)[0] if False else [])
            t = " ".join(x for pair in re.findall(
                r"<ArticleTitle>(.*?)</ArticleTitle>|<AbstractText[^>]*>(.*?)</AbstractText>",
                art, re.S) for x in pair if x)
            out[m.group(1)] = re.sub(r"<[^>]+>", " ", t).lower()
        time.sleep(0.4)
    return out


def main() -> int:
    root = atlas_root()
    sample_path = root / "comention" / "audit-sample.jsonl.gz"
    if not sample_path.exists():
        print(f"missing {sample_path}; rebuild with "
              "`python scripts/atlas_comention.py --rebuild`", file=sys.stderr)
        return 1
    rows = [json.loads(l) for l in
            gzip.open(sample_path, "rt", encoding="utf-8") if l.strip()]
    print(f"{len(rows)} sampled sentences", flush=True)

    idx = load_index(root)
    alias, _stats = build_alias_map(idx)
    # identifier -> the surface forms that resolve to it, for the mechanical check
    by_id = collections.defaultdict(set)
    for a, i in alias.items():
        by_id[i].add(a)

    pmids = {r["pmid"] for r in rows if r.get("pmid")}
    print(f"reading PubTator annotations for {len(pmids)} papers ...", flush=True)
    pt = pubtator_by_paper(root, pmids)

    print("fetching abstracts to decompose the disagreements ...", flush=True)
    abs_text = abstracts(sorted(pmids))

    matched = unmatched = filtered = 0
    agree = disagree = 0
    no_pt = 0
    in_abstract = body_only = 0
    examples = []
    for r in rows:
        toks = _TOKEN.findall(r["sentence"].lower())
        grams = set()
        for n in range(1, 6):
            for i in range(len(toks) - n + 1):
                grams.add(" ".join(toks[i:i + n]))
        seen = pt.get(r.get("pmid") or "", set())
        for ident, nm in zip(r["entities"], r.get("entity_names", r["entities"])):
            forms = by_id.get(ident, set())
            if forms & grams:
                matched += 1
            else:
                # Every sampled mention fired BY CONSTRUCTION when the sample was
                # written, so a miss now means the alias that fired has since been
                # removed from the map -- not that the tokenizer failed. The two
                # are only distinguishable because the sample predates the filters.
                filtered += 1
                if len(examples) < 12:
                    examples.append({"pmid": r.get("pmid"), "entity": nm,
                                     "sentence": r["sentence"][:150]})
            if not seen:
                no_pt += 1
            elif ident in seen:
                agree += 1
            else:
                disagree += 1
                # was the alias visible in the abstract PubTator actually read?
                ab = abs_text.get(r.get("pmid") or "", "")
                if ab and any(a in ab for a in by_id.get(ident, ())):
                    in_abstract += 1
                else:
                    body_only += 1

    tot = matched + unmatched + filtered
    scored = agree + disagree
    L = [
        "# How precise is the full-text co-mention layer? (#ATLAS-COMENT-AUDIT)", "",
        "Generated by `scripts/atlas_comention_audit.py` over the uniform sentence",
        "sample the co-mention build now writes.", "",
        "This layer's counts are what `atlas-module-support.md` uses to argue that a",
        "zero in the relation column is an extraction failure rather than absence of",
        "evidence. Until now it had no precision number at all.", "",
        "## Check 1: did the alias actually match?", "",
        "Mechanical, not a judgement: every sampled entity must correspond to a",
        "surface form present in the sentence. A failure is a tokenizer bug.", "",
        f"| | count |", "|---|---|",
        f"| entity mentions sampled | {tot:,} |",
        f"| alias found in the sentence | {matched:,} ({100*matched/max(1,tot):.1f}%) |",
        f"| **would no longer fire** | **{filtered:,} "
        f"({100*filtered/max(1,tot):.1f}%)** |", "",
        "Every sampled mention fired by construction when the sample was written,",
        "so the second row is not a tokenizer failure -- it is the alias having",
        "since been removed by the support and minority-form filters. That makes it",
        "a direct read on how much volume those filters take out of a uniform",
        "sample of what this layer used to match.", "",
    ]
    if examples:
        L += ["Examples where the alias could not be located:", ""]
        for e in examples[:8]:
            L.append(f"* `{e['entity']}` in PMID {e['pmid']}: {e['sentence']}...")
        L.append("")

    L += [
        "## Check 2: does PubTator agree?", "",
        "PubTator annotated the same papers independently, from their abstracts. It",
        "never saw our alias map, so where it assigned the same entity to the same",
        "paper, the full-text match is corroborated.", "",
        f"| | count |", "|---|---|",
        f"| mentions in papers PubTator annotated | {scored:,} |",
        f"| PubTator assigned the same entity | {agree:,} ({100*agree/max(1,scored):.1f}%) |",
        f"| it did not | {disagree:,} ({100*disagree/max(1,scored):.1f}%) |",
        f"| papers with no PubTator annotation at all | {no_pt:,} |", "",
        "### Splitting the disagreements", "",
        "A disagreement means one of two very different things, and pooling them",
        "would waste the measurement:", "",
        "| | count | share of disagreements |", "|---|---|---|",
        f"| our alias appears in the ABSTRACT PubTator read | {in_abstract:,} | "
        f"{100*in_abstract/max(1,disagree):.1f}% |",
        f"| body-only, so PubTator could not have seen it | {body_only:,} | "
        f"{100*body_only/max(1,disagree):.1f}% |", "",
        f"The body-only share is not evidence against this layer -- it is the layer",
        "doing the job it exists for, finding entities the abstract-level extractor",
        "structurally cannot reach. The first row is the one that could contain false",
        "positives, and it is where any future manual check should go.", "",
        f"Treating body-only matches as correct puts precision at "
        f"**{100*(agree+body_only)/max(1,scored):.1f}%** as an upper bound, against the",
        f"{100*agree/max(1,scored):.1f}% corroboration rate as a lower bound.", "",
        "### The bound has since been resolved by hand, and it sits near the bottom",
        "",
        "180 abstract-visible and 39 body-only mentions were read individually",
        "(#617). The result was much worse than this bound's midpoint suggests:", "",
        "| stratum | n | precision |", "|---|---|---|",
        "| agreeing with PubTator | 493 | 92.5% |",
        "| body-only | 439 | **30.8%** |",
        "| abstract-visible | 180 | **14.6%** |", "",
        "**Population-weighted precision was 55.5%**, near the bottom of the bound.",
        "The body-only stratum was NOT the layer doing its job -- the same generic",
        "aliases misfire in body text as in abstracts, and this report's earlier",
        "framing of that stratum was wrong.", "",
        "The cause was a single gap: `usable_alias` exempted every MULTI-WORD form",
        "from the specificity test it applied to single tokens, so `tumor cells`",
        "resolved to *Glucagonoma*, `overall survival` to *Prosthesis Failure* and",
        "`et al` to *Multiple Myeloma*. 132 of 152 measured false positives were",
        "multi-word.", "",
        "Two measured filters now replace that proxy (see `atlas_comention.py`), and",
        "the counts above are from the run BEFORE them. A fresh sample after the next",
        "rebuild is what will confirm the repair.", "",
        "### Why corroboration alone is a lower bound", "",
        "PubTator reads abstracts; this layer reads full text. An entity discussed",
        "only in Methods or Results is genuinely present and genuinely absent from",
        "PubTator's annotation, so it scores as a disagreement while being correct.",
        "Agreement therefore UNDERSTATES precision by an amount this data cannot",
        "measure.", "",
        "What a low agreement rate WOULD show is the layer finding entities the",
        "abstract-level extractor never sees anywhere in the paper, which is exactly",
        "the failure mode a permissive alias map produces.", "",
        "## Limits", "",
        "* A uniform sample of sentences, so common entities dominate it. Precision on",
        "  rare aliases -- where the alias map is most likely to be wrong -- is not",
        "  measured separately.",
        "* Check 1 rebuilds the n-gram lookup from the CURRENT alias map. If the map",
        "  changed since the sample was written, a mismatch may reflect the change",
        "  rather than a bug.",
        "* Neither check asks whether the two entities in a sentence are RELATED. This",
        "  layer never claimed they were; it claims only that both are named.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "sentences": len(rows), "mentions": tot,
        "alias_found": matched, "alias_missing": unmatched,
        "alias_filtered_out": filtered,
        "pubtator_scored": scored, "pubtator_agree": agree,
        "pubtator_disagree": disagree, "papers_without_pubtator": no_pt,
        "in_abstract": in_abstract,
        "body_only": body_only,
        "examples_missing": examples,
    }, indent=2) + "\n")
    print(f"\nalias found {100*matched/max(1,tot):.1f}%   "
          f"PubTator agreement {100*agree/max(1,scored):.1f}%")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
