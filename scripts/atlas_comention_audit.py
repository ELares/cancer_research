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

import argparse
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


def _measured_strata_lines() -> list:
    """Report the hand-judged strata, ordered by the data.

    An earlier version hardcoded "20.0% ... the worst of the three". Body-only
    is the SECOND worst -- abstract-visible is lower -- and the two are not
    distinguishable at these sample sizes. A superlative written by hand beside
    numbers written by hand is exactly the drift this report keeps finding
    elsewhere, so the ordering is now read from the measurement.
    """
    f = PROJECT_ROOT / "analysis" / "comention-regression.json"
    try:
        d = json.loads(f.read_text())
        m = d["measured_strata"]
    except (OSError, ValueError, KeyError):
        return ["Hand-judged strata precisions are not available; run "
                "`scripts/comention_regression.py`.", ""]
    rows = [
        ("body-only", m["body_only"]["precision"], m["body_only"]["n"],
         m["body_only"]["ci"]),
        ("abstract-visible", d["after"]["abstract_precision"],
         d["after"]["judged_n"], d["after"]["abstract_precision_ci"]),
        ("corroborated", m["agree"]["precision"], m["agree"]["n"], m["agree"]["ci"]),
    ]
    rows.sort(key=lambda r: r[1])
    body = next(r for r in rows if r[0] == "body-only")
    lowest = rows[0]
    overlap = body[3][0] <= lowest[3][1] and lowest[3][0] <= body[3][1]
    L = ["Hand judging puts the three strata at "
         "(`analysis/comention/*-judgements.csv`):", "",
         "| stratum | precision | 95% CI | n |", "|---|---|---|---|"]
    for name, prec, n_, ci in rows:
        L.append(f"| {name} | **{100*prec:.1f}%** | "
                 f"[{100*ci[0]:.1f}%, {100*ci[1]:.1f}%] | {n_} |")
    tail = (", and it is not distinguishable from the lowest at these sample sizes "
            "-- the intervals overlap across most of their range."
            if overlap and body is not lowest else ".")
    L += ["",
          f"Body-only is the {'lowest' if body is lowest else 'second lowest'}{tail}",
          "It carries the most volume of the three, which is what makes it matter.",
          ""]
    return L


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
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dump-strata", metavar="DIR",
                    help="write body-only and corroborated mentions to DIR as "
                         "JSONL, so those strata can be hand-judged instead of "
                         "carrying assumed precisions")
    ap.add_argument("--dump-abstract-visible", metavar="PATH",
                    help="write the abstract-visible disagreements to PATH as "
                         "JSONL, so the stratum can be hand-judged")
    args = ap.parse_args()

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
    abstract_visible = []  # the only stratum PubTator could have contradicted
    body_only_rows, corroborated_rows = [], []
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
            rec = {"pmid": r.get("pmid"), "identifier": ident, "entity": nm,
                   "sentence": r["sentence"][:600]}
            if not seen:
                no_pt += 1
            elif ident in seen:
                agree += 1
                if len(corroborated_rows) < 5000:
                    corroborated_rows.append(rec)
            else:
                disagree += 1
                # was the alias visible in the abstract PubTator actually read?
                ab = abs_text.get(r.get("pmid") or "", "")
                if ab and any(a in ab for a in by_id.get(ident, ())):
                    in_abstract += 1
                    if len(abstract_visible) < 5000:
                        abstract_visible.append({
                            "pmid": r.get("pmid"), "identifier": ident,
                            "entity": nm, "sentence": r["sentence"][:400]})
                else:
                    body_only += 1
                    if len(body_only_rows) < 5000:
                        body_only_rows.append(rec)

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
        "A body-only match is one PubTator could not have contradicted, because it",
        "never read that text. That is what this layer exists to find -- and it is",
        "NOT a reason to score the stratum as correct. This report asserted BOTH",
        "readings at once -- that the body-only share was \"the layer doing the job",
        "it exists for\" here, and that it was \"NOT the layer doing its job\" in the",
        "#617 section below. The second is the right one.", "",
    ] + _measured_strata_lines() + [
        "The same generic forms misfire in body text as in abstracts, and there is",
        "no abstract-level extractor to disagree with them.", "",
        f"So read the {100*(agree+body_only)/max(1,scored):.1f}% that treating",
        "body-only as correct would give as a bound that is now known to be far",
        f"above the truth, and the {100*agree/max(1,scored):.1f}% corroboration rate",
        "as a lower bound that is close to it. All three strata are measured in",
        "`analysis/comention-regression.md`.", "",
        "### The PRE-FILTER measurement (#617), kept for comparison",
        "",
        "180 abstract-visible and 39 body-only mentions were read individually",
        "(#617), on the sample from BEFORE the filter change. Retained because it is",
        "the only measurement of that run, and because its body-only n is an order of",
        "magnitude larger than the post-filter one. **The two are not like-for-like**:",
        "this pass judged whether the sentence contained the matched string, while the",
        "post-filter pass asked whether the sentence discusses the entity the",
        "IDENTIFIER denotes, which is stricter. Read the 30.8% below as an upper bound",
        "on what a strict re-judging of that run would have given.", "",
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
    if args.dump_strata:
        d = Path(args.dump_strata)
        d.mkdir(parents=True, exist_ok=True)
        for name, rowset in (("body-only", body_only_rows),
                             ("corroborated", corroborated_rows)):
            f = d / f"{name}.jsonl"
            f.write_text("\n".join(json.dumps(e) for e in rowset) + "\n")
            print(f"wrote {f}: {len(rowset)} {name} mentions for hand judging")

    if args.dump_abstract_visible:
        dest = Path(args.dump_abstract_visible)
        dest.write_text("\n".join(json.dumps(e) for e in abstract_visible) + "\n")
        print(f"wrote {dest}: {len(abstract_visible)} abstract-visible disagreements "
              f"for hand judging")

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
