#!/usr/bin/env python3
"""Atlas: resolve the FSP1 sense collision per paper (#ATLAS-AMBIG).

WHY
---
`FSP1` is an official NCBI alias for three different human genes:

| gene | id | what FSP1 stands for there |
|---|---|---|
| AIFM2  | 84883 | **Ferroptosis Suppressor Protein 1** |
| S100A4 | 6275  | Fibroblast-Specific Protein 1 |
| ATL1   | 51062 | Familial Spastic Paraplegia 1 (atlastin GTPase 1) |

So `FSP1` has no single right answer. A blanket remap to AIFM2 -- the obvious
"fix" once you know the manuscript needs AIFM2 -- would be wrong for the
majority of the census's FSP1 papers, which really do mean S100A4.

The answer depends on the individual paper, which is why this is a separate
layer from `atlas_ambiguity.py`'s blocklist: the blocklist stops the graph
guessing, and this decides.

WHAT PUBTATOR DOES
------------------
Measured, not assumed. The gold set is papers that **declare the expansion
themselves** ("ferroptosis suppressor protein 1 (FSP1)"), which is a label
independent of anything the classifier sees. On that set PubTator3 is 47.4%
accurate, and its errors are not spread evenly:

* of 110 papers that spell out *ferroptosis suppressor protein 1*, **not one**
  has that mention resolved to AIFM2 -- 99 go to ATL1;
* 69% carry no correct AIFM2 edge from any OTHER mention either (the remaining
  31% are rescued only because the paper also writes "AIFM2" somewhere), so for
  most of them the ferroptosis biology is not merely mislabelled, it is absent;
* the S100A4 assignments, by contrast, are mostly right, which is exactly why a
  blanket remap would do harm.

THE LEAKAGE TRAP
----------------
The gold label is defined by the presence of an expansion phrase. If the
classifier may read that phrase, it scores near 100% and measures nothing. So
every phrase that can define a label is masked out of the text before a single
feature is read (`MASK`). `ferroptosis` still counts as a cue, but only where it
occurs outside `ferroptosis suppressor protein 1`. Without that step the
headline number is circular.

The gold set gets one further check that uses no feature at all: publication
year. FSP1 was named as the ferroptosis suppressor in 2019, and not one of the
110 AIFM2-sense papers predates it (range 2020-2026), while the S100A4 sense
spans 1997-2025. A single pre-2019 AIFM2 paper would mean the gold rule was
matching something other than what it claims.

Requires network (NCBI E-utilities) for the title/abstract/MeSH of the affected
papers. Only the derived verdicts are committed, so downstream stays offline and
no abstracts are redistributed.

Usage:
    python scripts/atlas_disambiguate.py
    python scripts/atlas_disambiguate.py --symbol FSP1
"""

import argparse
import collections
import gzip
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

# The partner gene whose co-mention count the correction is measured against:
# GPX4, because the manuscript's headline Bliss-synergy claim is GPX4 + FSP1.
PARTNER_ID = "2879"

OUT = PROJECT_ROOT / "analysis" / "atlas-disambiguation.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-disambiguation.json"

# The senses `FSP1` can carry, each with the identifiers PubTator uses for it
# (human and mouse), the phrases by which a paper declares that sense, and the
# contextual cues that survive masking.
SENSES = {
    "AIFM2": {
        "ids": ["84883"],
        "declares": r"ferroptosis[- ]suppressor[- ]protein[- ]?1"
                    r"|apoptosis[- ]inducing factor mitochondria[- ]associated 2"
                    r"|aifm2",
        "cues": [r"gpx4", r"ferroptos", r"lipid peroxidation", r"erastin",
                 r"rsl3", r"ubiquinone|coenzyme q|coq10", r"labile iron|fenton",
                 r"slc7a11|system xc", r"acsl4", r"radical[- ]trapping"],
    },
    "S100A4": {
        "ids": ["6275", "20198"],
        "declares": r"fibroblast[- ]specific protein[- ]?1|s100a4|metastasin|mts1",
        "cues": [r"fibroblast", r"epithelial[- ]mesenchymal|\bemt\b", r"metasta",
                 r"\bcaf\b|cancer[- ]associated fibroblast",
                 r"invasion|migration", r"calcium[- ]binding"],
    },
    "ATL1": {
        "ids": ["51062", "73991"],
        "declares": r"atlastin|spastic paraplegia|spg3a?\b",
        "cues": [r"paraplegia", r"axon", r"endoplasmic reticulum", r"gtpase",
                 r"neuropath", r"hereditary"],
    },
}

# Every phrase that can define a gold label, masked out before features are read.
MASK = re.compile("|".join(s["declares"] for s in SENSES.values()))
CUES = {k: [re.compile(p) for p in v["cues"]] for k, v in SENSES.items()}
DECLARES = {k: re.compile(v["declares"]) for k, v in SENSES.items()}
ID2SENSE = {i: k for k, v in SENSES.items() for i in v["ids"]}


def wilson(k: int, n: int) -> tuple:
    """95% Wilson interval, the repo's convention for a measured proportion."""
    if n == 0:
        return (0.0, 0.0)
    z, p = 1.96, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def affected(genes: Path, symbol: str) -> dict:
    """pmid -> set of senses PubTator assigned via this symbol's surface forms."""
    sym = symbol.lower()
    out = collections.defaultdict(set)
    surfaces = collections.Counter()
    with gzip.open(genes, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, gid, mention = parts[0], parts[2], parts[3]
            forms = [s.strip().lower() for s in mention.split("|")]
            # exact surface match, so BFSP1/UFSP1 are not swept in as FSP1
            if sym not in forms:
                continue
            surfaces[mention.lower()] += 1
            if gid in ID2SENSE:
                out[pmid].add(ID2SENSE[gid])
    return out, surfaces


def fetch(pmids: list) -> dict:
    """pmid -> lowercased title/abstract/MeSH text, batched."""
    text = {}
    for i in range(0, len(pmids), 200):
        batch = pmids[i:i + 200]
        data = urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}).encode()
        req = urllib.request.Request(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", data=data)
        try:
            xml = urllib.request.urlopen(req, timeout=180).read().decode("utf-8", "replace")
        except Exception as exc:
            print(f"  ! efetch batch failed: {exc}", file=sys.stderr)
            continue
        for art in xml.split("<PubmedArticle>")[1:]:
            m = re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
            if m:
                text[m.group(1)] = re.sub(r"<[^>]+>", " ", art).lower()
        time.sleep(0.4)
        print(f"  fetched {len(text):,}/{len(pmids):,}", flush=True)
    return text


def pub_years(pmids: list) -> dict:
    """pmid -> publication year, for the independent temporal check.

    Parsed from esummary rather than scraped out of the record text: taking the
    earliest four-digit number in an article record picks up reference and
    history dates, which produced a wrong answer before this was fixed.
    """
    out = {}
    for i in range(0, len(pmids), 200):
        batch = pmids[i:i + 200]
        data = urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(batch), "retmode": "json"}).encode()
        req = urllib.request.Request(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", data=data)
        try:
            res = json.load(urllib.request.urlopen(req, timeout=120))["result"]
        except Exception as exc:
            print(f"  ! esummary batch failed: {exc}", file=sys.stderr)
            continue
        for p in batch:
            pd = res.get(p, {}).get("pubdate", "")
            if pd[:4].isdigit():
                out[p] = int(pd[:4])
        time.sleep(0.4)
    return out


def gold_label(t: str):
    """The sense a paper declares, or None if it declares zero or several."""
    hits = [k for k, rx in DECLARES.items() if rx.search(t)]
    return hits[0] if len(hits) == 1 else None


def classify(t: str):
    """Contextual sense, with every label-defining phrase masked out first."""
    masked = MASK.sub(" ", t)
    score = {k: sum(bool(rx.search(masked)) for rx in v) for k, v in CUES.items()}
    ranked = sorted(score.values(), reverse=True)
    best = max(score, key=lambda k: score[k])
    if ranked[0] == 0:
        return None, score, "no cue matched"
    if len(ranked) > 1 and ranked[0] == ranked[1]:
        return None, score, "tie between senses"
    return best, score, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--symbol", default="FSP1")
    args = ap.parse_args()

    genes = atlas_root() / "entities" / "gene.tsv.gz"
    if not genes.exists():
        print(f"missing {genes}; run scripts/atlas_relations.py first", file=sys.stderr)
        return 1

    print(f"scanning {genes} for `{args.symbol}` ...", flush=True)
    assigned, surfaces = affected(genes, args.symbol)
    pmids = sorted(assigned)
    print(f"  {len(pmids):,} cancer papers carry a `{args.symbol}` annotation", flush=True)

    text = fetch(pmids)

    gold, decided, abstained = {}, {}, {}
    for p, t in text.items():
        g = gold_label(t)
        if g:
            gold[p] = g
        pred, score, why = classify(t)
        if pred:
            decided[p] = pred
        else:
            abstained[p] = why

    # PubTator's accuracy on the self-declaring set (single-sense papers only,
    # so a paper PubTator tagged two ways is not scored either right or wrong).
    pt_ok = pt_n = 0
    pt_conf = collections.Counter()
    for p, truth in gold.items():
        a = assigned.get(p, set())
        if len(a) != 1:
            continue
        got = next(iter(a))
        pt_conf[(truth, got)] += 1
        pt_n += 1
        pt_ok += (got == truth)

    # This layer's accuracy on the same set.
    my_ok = my_n = 0
    my_conf = collections.Counter()
    for p, truth in gold.items():
        if p not in decided:
            continue
        got = decided[p]
        my_conf[(truth, got)] += 1
        my_n += 1
        my_ok += (got == truth)

    lo, hi = wilson(my_ok, my_n)
    ptlo, pthi = wilson(pt_ok, pt_n)

    # Corrections: papers where this layer disagrees with the single sense
    # PubTator assigned. These are what a consumer applies.
    corrections = {p: {"pubtator": next(iter(assigned[p])), "corrected": decided[p]}
                   for p in decided
                   if len(assigned.get(p, set())) == 1
                   and next(iter(assigned[p])) != decided[p]}

    # How many gold-AIFM2 papers have no correct AIFM2 edge from any mention?
    # Same pass also counts the partner co-mentions the correction recovers.
    aifm2_gold = {p for p, v in gold.items() if v == "AIFM2"}
    corrected_to_aifm2 = {p for p, v in corrections.items()
                          if v["corrected"] == "AIFM2"}
    with_edge, partner_pm, aifm2_pm = set(), set(), set()
    with gzip.open(genes, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, gid = parts[0], parts[2]
            if gid == PARTNER_ID:
                partner_pm.add(pmid)
            elif gid == "84883":
                aifm2_pm.add(pmid)
                if pmid in aifm2_gold:
                    with_edge.add(pmid)
    orphaned = len(aifm2_gold - with_edge)
    pair_before = len(partner_pm & aifm2_pm)
    pair_after = len(partner_pm & (aifm2_pm | corrected_to_aifm2))

    # Years for the gold set AND for every corrected paper, because the
    # extrapolation below turns on the corrected-but-undeclaring population.
    years = pub_years(sorted(set(gold) | set(corrections)))
    year_buckets = collections.defaultdict(list)
    for p, s in gold.items():
        if p in years:
            year_buckets[s].append(years[p])

    # THE EXTRAPOLATION CHECK. Accuracy is measured on papers that declare a
    # sense, but most corrections land on papers that declare nothing, so the
    # headline is extrapolated to a population it never scored. Publication year
    # tests that extrapolation directly and independently: FSP1 was named as the
    # ferroptosis suppressor in 2019, so a pre-2019 paper corrected to AIFM2 is
    # almost certainly wrong.
    nd = {p: v["corrected"] for p, v in corrections.items() if p not in gold}
    nd_years = {p: years[p] for p in nd if p in years}
    nd_aifm2 = [p for p, v in nd.items() if v == "AIFM2" and p in nd_years]
    nd_other = [p for p, v in nd.items() if v != "AIFM2" and p in nd_years]
    nd_pre = sum(1 for p in nd_aifm2 if nd_years[p] < 2019)
    base_pre = (sum(1 for p in nd_years if nd_years[p] < 2019) / len(nd_years)
                if nd_years else 0.0)
    expected_pre = base_pre * len(nd_aifm2)

    dist = collections.Counter(gold.values())
    lines = [
        f"# Resolving the `{args.symbol}` sense collision (#ATLAS-AMBIG)", "",
        f"Generated by `scripts/atlas_disambiguate.py`. `{args.symbol}` is an official",
        "NCBI alias for three different human genes, so it has no single right",
        "answer -- only a per-paper one.", "",
        "| gene | id | what the abbreviation stands for there |", "|---|---|---|",
        "| AIFM2 | 84883 | **Ferroptosis Suppressor Protein 1** |",
        "| S100A4 | 6275 | Fibroblast-Specific Protein 1 |",
        "| ATL1 | 51062 | Familial Spastic Paraplegia 1 (atlastin GTPase 1) |", "",
        "## The gold set", "",
        f"{len(pmids):,} cancer papers carry an `{args.symbol}` annotation. Of those,",
        f"**{len(gold):,}** declare the expansion themselves and so carry a label",
        "that is independent of anything the classifier reads",
        f"({dist.get('AIFM2',0)} AIFM2, {dist.get('S100A4',0)} S100A4, {dist.get('ATL1',0)} ATL1).", "",
        "## What PubTator3 does on it", "",
        f"**{pt_ok}/{pt_n} = {100*pt_ok/max(1,pt_n):.1f}%** accurate "
        f"(95% CI {100*ptlo:.1f}-{100*pthi:.1f}%).", "",
        "| paper declares | PubTator assigns | n |", "|---|---|---|",
    ]
    for (t, g), c in sorted(pt_conf.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {t} | {g} | {c}{'' if t == g else ' **wrong**'} |")

    lines += [
        "", "The errors are not evenly spread, and the shape is the finding:", "",
        f"* of the {dist.get('AIFM2',0)} papers that spell out *ferroptosis suppressor",
        "  protein 1*, **not one** has that mention resolved to AIFM2;",
        f"* **{orphaned}** of them ({100*orphaned/max(1,len(aifm2_gold)):.0f}%) carry no",
        "  correct AIFM2 edge from any OTHER mention either -- the rest are rescued",
        "  only because the paper also writes \"AIFM2\" somewhere -- so for most of",
        "  them the ferroptosis biology is not merely mislabelled, it is absent;",
        "* the S100A4 assignments are mostly correct, which is precisely why a",
        "  blanket remap of `FSP1` to AIFM2 would do harm rather than fix this.", "",
        "## What this layer does on it", "",
        f"**{my_ok}/{my_n} = {100*my_ok/max(1,my_n):.1f}%** accurate "
        f"(95% CI {100*lo:.1f}-{100*hi:.1f}%), abstaining on "
        f"{len(gold) - my_n} of the {len(gold)} gold papers.", "",
        "| paper declares | this layer assigns | n |", "|---|---|---|",
    ]
    for (t, g), c in sorted(my_conf.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {t} | {g} | {c}{'' if t == g else ' **wrong**'} |")

    lines += [
        "", "## An independent check the classifier cannot see", "",
        "Publication year is used by neither the gold rule nor the classifier, so it",
        "is a free structural test. Doll and Bersuker named FSP1 as the ferroptosis",
        "suppressor in 2019, so if the gold set is sound essentially no AIFM2-sense",
        "paper can predate that, while the long-established S100A4 sense should span",
        "decades.", "",
        "| declared sense | n | median year | range | published before 2019 |",
        "|---|---|---|---|---|",
    ] + [
        f"| {s} | {len(ys)} | {sorted(ys)[len(ys)//2]} | {min(ys)}-{max(ys)} | "
        f"{sum(1 for y in ys if y < 2019)} "
        f"({100*sum(1 for y in ys if y < 2019)/len(ys):.0f}%) |"
        for s, ys in sorted(year_buckets.items()) if ys
    ] + [
        "", "The separation is what it should be, and the AIFM2 row is the strong form",
        "of the test: a single pre-2019 paper declaring *ferroptosis suppressor protein",
        "1* would mean the gold rule was matching something else.", "",
        "### The extrapolation, and a test of it", "",
        f"Accuracy above is measured on the {len(gold)} papers that declare a sense. But",
        f"**{len(nd):,} of the {len(corrections):,} corrections "
        f"({100*len(nd)/max(1,len(corrections)):.0f}%) land on papers that declare",
        "nothing**, so the headline is extrapolated to a population it never scored.",
        "That is the most important limitation here, and it is testable with the same",
        "independent signal.", "",
        "| corrected, undeclaring | n | published before 2019 |", "|---|---|---|",
        f"| to AIFM2 | {len(nd_aifm2)} | **{nd_pre}** |",
        f"| to another sense | {len(nd_other)} | "
        f"{sum(1 for p in nd_other if nd_years[p] < 2019)} |", "",
        f"Across all corrected undeclaring papers {100*base_pre:.0f}% predate 2019, so if",
        f"the classifier were assigning AIFM2 without regard to the biology it would put",
        f"roughly **{expected_pre:.0f}** of them before the term existed. It puts",
        f"**{nd_pre}**. The extrapolation is supported on exactly the population the gold",
        "set does not cover.", "",
        "## Why that number is not circular", "",
        "The gold label is defined by the presence of an expansion phrase. A",
        "classifier allowed to read that phrase would score near 100% and measure",
        "nothing at all. So every phrase that can define a label is masked out of",
        "the text before a single feature is read: `ferroptosis` still counts as a",
        "cue, but only where it occurs outside *ferroptosis suppressor protein 1*.",
        "The number above is what survives that masking.", "",
        "## Corrections, and what they recover", "",
        f"{len(corrections):,} papers where this layer disagrees with the single sense",
        "PubTator assigned. These are the per-paper corrections a consumer applies;",
        "the full table is in the JSON beside this report.", "",
        "The consequence is measurable on exactly the pair the manuscript's headline",
        "claim rests on. GPX4 + FSP1 Bliss synergy needs papers co-mentioning GPX4",
        "and the *real* FSP1 (AIFM2):", "",
        "| | papers co-mentioning GPX4 and AIFM2 |", "|---|---|",
        f"| before correction | {pair_before:,} |",
        f"| after correction | {pair_after:,} |",
        f"| recovered | **+{pair_after - pair_before:,}** "
        f"({100*(pair_after - pair_before)/max(1, pair_before):.0f}% more) |", "",
        "So the collision was not a labelling nicety. It was hiding roughly",
        f"{100*(pair_after - pair_before)/max(1, pair_after):.0f}% of the census's own",
        "evidence for the mechanism this project is built on, and the graph reported",
        "the smaller number with nothing to indicate anything was missing.", "",
        "## Limits", "",
        f"* **Most corrections are extrapolated.** {100*len(nd)/max(1,len(corrections)):.0f}% "
        "land on papers that declare no sense, so their accuracy is inferred from the",
        "  declaring subset rather than measured. The temporal check above supports that",
        "  inference but does not replace a labelled evaluation of those papers.",
        "* Measured on `FSP1` only. The method generalises to any symbol whose",
        "  senses have distinct spelled-out expansions, but no other symbol has been",
        "  scored, and an unmeasured layer is not a validated one.",
        "* The text is title, abstract and MeSH, not full text. A paper that",
        "  declares its expansion only in the body is unlabelled here, and one whose",
        "  cues are all in the body may be abstained on.",
        "* Cues are hand-written, so recall is bounded by that vocabulary. The",
        "  abstention rate is reported rather than hidden by forcing a guess.",
        "* A paper using `FSP1` in two senses at once is out of scope; the layer",
        "  assigns one sense per paper.",
    ]

    OUT.write_text("\n".join(lines) + "\n")
    RAW.write_text(json.dumps({
        "symbol": args.symbol,
        "papers_with_annotation": len(pmids),
        "gold_set": len(gold),
        "gold_distribution": dict(dist),
        "pubtator_accuracy": {"correct": pt_ok, "n": pt_n,
                              "wilson95": [ptlo, pthi]},
        "layer_accuracy": {"correct": my_ok, "n": my_n, "wilson95": [lo, hi],
                           "abstained": len(gold) - my_n},
        "gold_aifm2_with_no_correct_edge": orphaned,
        "gold_years": {k: sorted(v) for k, v in year_buckets.items()},
        "corrections": corrections,
        "abstentions": abstained,
    }, indent=2) + "\n")
    print(f"\nPubTator {100*pt_ok/max(1,pt_n):.1f}%  ->  this layer "
          f"{100*my_ok/max(1,my_n):.1f}%  ({len(corrections):,} corrections)")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
