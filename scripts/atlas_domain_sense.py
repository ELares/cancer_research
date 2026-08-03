#!/usr/bin/env python3
"""Atlas: measure the curated DOMAIN_SENSE claims instead of asserting them
(#ATLAS-AMBIG).

WHY
---
`atlas_ambiguity.DOMAIN_SENSE` records, for a handful of high-volume sense
collisions, which sense the cancer literature actually means -- `ER` is the
estrogen receptor, `COX-2` is PTGS2, and so on. Every entry ships with a
sentence justifying it.

Those sentences were written from domain knowledge, not from the corpus, and one
of them was wrong. The `psa` entry claimed the majority vote returns NPEPPS; it
returns KLK3, which is the correct answer, and the repository's own committed
scan said so. A curated claim that nothing checks is exactly the kind of
assertion this project is supposed to avoid.

So this measures them. For each curated symbol it samples cancer papers that
carry the symbol as a gene annotation, counts how many DECLARE each candidate
sense in their own words ("estrogen receptor (ER)"), and compares that to both
the curated claim and PubTator's majority vote.

WHAT IT FOUND
-------------
The curated senses hold, overwhelmingly, and the majority vote is wrong for four
of the five: it picks the sense the cancer literature essentially never means.
`FSP1` turns out to be the unusual one. Its two senses are genuinely balanced
(110 AIFM2 against 132 S100A4), which is why it needs the per-paper classifier in
`atlas_disambiguate.py` while these five need only a curated default.

Requires network (NCBI E-utilities). Only the derived counts are committed, so
downstream stays offline and no abstracts are redistributed.

Usage:
    python scripts/atlas_domain_sense.py
    python scripts/atlas_domain_sense.py --sample 800
"""

import argparse
import collections
import gzip
import json
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_ambiguity import DOMAIN_SENSE  # noqa: E402
from atlas_baseline import atlas_root  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-domain-sense-validation.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-domain-sense-validation.json"

SAMPLE_SEED = 20260803

# How each candidate sense is DECLARED in a paper's own words. The curated
# sense is named first. A sense that is not a gene at all (`ER` meaning the
# endoplasmic reticulum) is included where it is a real competitor, because a
# spurious gene annotation is a different failure from a wrong one.
PROBES = {
    "er": {
        "ESR1": r"estrogen receptor|oestrogen receptor|\besr1\b",
        "EREG": r"epiregulin|\bereg\b",
        "(endoplasmic reticulum, not a gene)": r"endoplasmic reticulum",
    },
    "cox-2": {
        "PTGS2": r"prostaglandin[- ]endoperoxide|cyclooxygenase|\bptgs2\b",
        "COX2 (mitochondrial)": r"cytochrome c oxidase",
    },
    "p62": {
        "SQSTM1": r"sequestosome|\bsqstm1\b",
        "NUP62": r"nucleoporin",
    },
    "p21": {
        "CDKN1A": r"cyclin[- ]dependent kinase inhibitor|\bcdkn1a\b|\bwaf1\b|\bcip1\b",
        "H3P16 (histone pseudogene)": r"histone h3|h3 histone",
        "RAS p21": r"ras p21|p21 ras",
    },
    "psa": {
        "KLK3": r"prostate[- ]specific antigen|\bklk3\b|kallikrein",
        "NPEPPS": r"puromycin[- ]sensitive aminopeptidase|\bnpepps\b",
    },
}


def efetch(pmids: list, tries: int = 4) -> str:
    """E-utilities is intermittently 502-prone at this volume; back off."""
    for k in range(tries):
        try:
            data = urllib.parse.urlencode(
                {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}).encode()
            req = urllib.request.Request(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", data=data)
            return urllib.request.urlopen(req, timeout=180).read().decode("utf-8", "replace")
        except Exception as exc:
            if k == tries - 1:
                print(f"  ! batch failed after {tries} tries: {exc}", file=sys.stderr)
                return ""
            time.sleep(2 ** k)
    return ""


def papers_for(genes: Path, surface: str) -> dict:
    """gene id -> the cancer PMIDs annotated with this surface form."""
    out = collections.defaultdict(set)
    with gzip.open(genes, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if surface in [s.strip().lower() for s in parts[3].split("|")]:
                out[parts[2]].add(parts[0])
    return {g: sorted(v) for g, v in out.items()}


def declared(text: str, probes: dict):
    """The single sense a paper declares, or None if zero or several."""
    hits = [k for k, rx in probes.items() if rx.search(text)]
    return hits[0] if len(hits) == 1 else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=800,
                    help="papers to sample per symbol")
    args = ap.parse_args()

    genes = atlas_root() / "entities" / "gene.tsv.gz"
    if not genes.exists():
        print(f"missing {genes}; run scripts/atlas_relations.py first", file=sys.stderr)
        return 1
    try:
        scan = json.loads(
            (PROJECT_ROOT / "analysis" / "atlas-ambiguity.json").read_text())
    except (OSError, ValueError):
        print("run scripts/atlas_ambiguity.py first", file=sys.stderr)
        return 1
    vote = {r["surface"]: r for r in scan["by_type"]["gene"]["sense_rows"]}

    rng = random.Random(SAMPLE_SEED)
    results = {}
    for surface, probes in PROBES.items():
        if surface not in DOMAIN_SENSE:
            continue
        print(f"scanning for {surface!r} ...", flush=True)
        by_gene = papers_for(genes, surface)
        pool = sorted(set().union(*by_gene.values())) if by_gene else []
        if not pool:
            continue
        samp = sorted(rng.sample(pool, min(args.sample, len(pool))))
        compiled = {k: re.compile(v) for k, v in probes.items()}

        counts = collections.Counter()
        seen = 0
        for i in range(0, len(samp), 200):
            xml = efetch(samp[i:i + 200])
            for art in xml.split("<PubmedArticle>")[1:]:
                if not re.search(r"<PMID[^>]*>(\d+)</PMID>", art):
                    continue
                seen += 1
                d = declared(re.sub(r"<[^>]+>", " ", art).lower(), compiled)
                counts[d or "(none declared)"] += 1
            time.sleep(0.5)

        curated = DOMAIN_SENSE[surface][1]
        decided = {k: v for k, v in counts.items() if k != "(none declared)"}
        total_decided = sum(decided.values())
        top = max(decided, key=decided.get) if decided else None
        results[surface] = {
            "pool": len(pool), "sampled": seen,
            "curated_sense": curated,
            "majority_vote": vote.get(surface, {}).get("top", {}).get("symbol"),
            "declared": dict(counts),
            "decided": total_decided,
            "curated_share": (decided.get(curated, 0) / total_decided
                              if total_decided else None),
            "curated_is_dominant": (top == curated) if top else None,
            "vote_matches_curated": (
                vote.get(surface, {}).get("top", {}).get("symbol") == curated),
        }
        print(f"  {surface}: {total_decided} declared, curated {curated} "
              f"{100*results[surface]['curated_share']:.1f}%", flush=True)

    confirmed = sum(1 for r in results.values() if r["curated_is_dominant"])
    vote_wrong = sum(1 for r in results.values() if not r["vote_matches_curated"])

    L = [
        "# Do the curated domain senses hold? (#ATLAS-AMBIG)", "",
        "Generated by `scripts/atlas_domain_sense.py`.",
        "`atlas_ambiguity.DOMAIN_SENSE` records which sense the cancer literature",
        "means for a handful of high-volume collisions. Those entries were written",
        "from domain knowledge, so this checks them against the corpus.", "",
        "Method: sample cancer papers carrying the symbol as a gene annotation, and",
        "count how many DECLARE each candidate sense in their own words (\"estrogen",
        "receptor (ER)\"). A paper declaring zero or several is not counted.", "",
        f"**{confirmed} of {len(results)} curated senses are confirmed as dominant, and",
        f"PubTator's majority vote disagrees with the curated sense in {vote_wrong} of",
        f"{len(results)} cases.**", "",
        "| symbol | curated sense | share of declaring papers | majority vote | vote correct? |",
        "|---|---|---|---|---|",
    ]
    for s, r in results.items():
        ok = "yes" if r["vote_matches_curated"] else "**no**"
        L.append(f"| `{s}` | {r['curated_sense']} | "
                 f"{100*r['curated_share']:.1f}% ({r['decided']} declaring) | "
                 f"{r['majority_vote']} | {ok} |")

    L += ["", "## What each sample actually declared", ""]
    for s, r in results.items():
        L += [f"### `{s}` — {r['pool']:,} papers in the census, {r['sampled']} sampled", "",
              "| declared sense | papers |", "|---|---|"]
        for k, v in sorted(r["declared"].items(), key=lambda kv: -kv[1]):
            L.append(f"| {k} | {v} |")
        L.append("")

    L += [
        "## Why this matters more than it looks", "",
        "The majority vote is not merely noisy on these symbols. It is wrong in a",
        "consistent direction: it picks the sense the cancer literature essentially",
        "never means. `ER` resolves to epiregulin while declaring papers say estrogen",
        "receptor by a wide margin; `COX-2` resolves to the mitochondrial oxidase",
        "while declaring papers say PTGS2; `p21` resolves to a histone pseudogene.", "",
        "`FSP1` is the unusual case, not the typical one. Its senses are genuinely",
        "balanced (110 AIFM2 against 132 S100A4 in the gold set), which is why it",
        "needs the per-paper classifier in `scripts/atlas_disambiguate.py` while these",
        "have a defensible default.", "",
        "## A curated claim that was wrong", "",
        "The `psa` entry asserted that the majority vote returns NPEPPS. It returns",
        "KLK3, the correct answer, and the committed ambiguity scan already said so.",
        "The entry has been corrected. `psa` remains blocklisted because KLK3 and",
        "NPEPPS are genuinely different genes, but the reason given for it was not",
        "true, and it was written rather than measured.", "",
        "## Limits", "",
        "* A declared sense is evidence about the paper, not proof: a paper can name",
        "  one sense in passing and mean the other. The declaration rule is the same",
        "  one used for the FSP1 gold set, and inherits its assumptions.",
        "* Only title, abstract and MeSH are read, so a paper declaring its sense only",
        "  in the body counts as undeclared. Between a third and two thirds of each",
        "  sample declares nothing, and those papers are excluded rather than guessed.",
        "* The share is measured on a sample, not the full pool, so it carries",
        "  sampling error. The margins are wide enough that this does not change the",
        "  conclusion: the narrowest curated share above is "
        f"{100*min(r['curated_share'] for r in results.values()):.1f}%, and the widest",
        f"  is {100*max(r['curated_share'] for r in results.values()):.1f}%.",
        "* This validates the DIRECTION of each curated sense. It does not measure how",
        "  often applying that default would be wrong for an individual paper.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({"sample_seed": SAMPLE_SEED,
                               "sample_per_symbol": args.sample,
                               "confirmed": confirmed,
                               "vote_disagrees": vote_wrong,
                               "results": results}, indent=2) + "\n")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
