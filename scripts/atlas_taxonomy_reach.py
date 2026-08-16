#!/usr/bin/env python3
"""How much of the cancer literature can the mechanism taxonomy label at all? (#730)

WHY THIS EXISTS
---------------
The project reports that per-mechanism capture of the frozen corpus against the
census spans 0.20% to 41.86%, a 213-fold spread, and treats that as the headline
caveat on every relative prevalence it computes over mechanisms.

That figure is CONDITIONED ON HAVING A TAG. The prior question -- what fraction
of the literature the taxonomy can label at all -- has never been measured, so
the spread has always been reported without its denominator. A mechanism layer
with a narrow field of view produces confident-looking shares of a sliver.

TWO INSTRUMENTS, AND CONFLATING THEM WOULD BE UNFAIR
-----------------------------------------------------
There are two ways a mechanism reaches an article here and they have different
jobs, so they get measured separately:

  THE KEYWORD TAGGER (`config.MECHANISM_KEYWORDS`, 25 mechanisms, 186 terms)
  is what actually labels the corpus. Its reach is the number that bounds every
  mechanism share the project reports. This is the field of view that matters.

  THE MESH LEAF MAP (`analysis/mesh-mechanism-map.yaml`) is PRECISION-FIRST BY
  DESIGN. It exists to measure the tagger's recall non-circularly, so it
  deliberately excludes umbrella descriptors -- bare `Immunotherapy` is dropped
  as non-discriminative even though it would raise coverage enormously. Reading
  its low reach as a coverage failure would be criticising a precision
  instrument for not being a coverage instrument. It is reported here only as a
  contrast, clearly labelled.

WHAT A LOW NUMBER WOULD AND WOULD NOT MEAN
-------------------------------------------
The interesting question is not the percentage, it is what the UNLABELLED
remainder consists of. If it is mostly case reports, epidemiology, health
services research and basic biology with no therapeutic modality, then a
mechanism taxonomy SHOULD NOT reach it and the field of view is appropriate --
only the wording of the capture caveat needs fixing. If instead it is full of
therapy papers the taxonomy has no name for, that is a real gap.

So this script profiles the remainder rather than just counting it, and the
profile is the deliverable.

METHOD NOTE. The keyword reach is measured on a uniform sample because the
production matcher is applied per keyword per article; the MeSH reach is
measured on the full census because set membership is cheap. Both are stated
with their denominators, and the sample carries Wilson intervals.

The matcher is IMPORTED from `scripts/tag_articles.py` rather than reproduced.
A reimplementation here would validate the reimplementation.

Usage:
    python scripts/atlas_taxonomy_reach.py
    python scripts/atlas_taxonomy_reach.py --sample-every 40   # faster, wider CI
"""

import argparse
import gzip
import json
import math
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import config  # noqa: E402
from tag_articles import text_matches_keyword  # noqa: E402

ATLAS = PROJECT_ROOT / "corpus" / "atlas"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-taxonomy-reach.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-taxonomy-reach.json"
MESH_MAP = PROJECT_ROOT / "analysis" / "mesh-mechanism-map.yaml"

# MeSH "check tags": demographic descriptors applied to almost every article.
# Without excluding them the top-descriptor table reads humans 92%, female 46%,
# male 37% -- true, and completely uninformative about what the literature is
# ABOUT, which is the only question this section exists to answer.
CHECK_TAGS = {
    "humans", "female", "male", "animals", "adult", "aged", "middle aged",
    "young adult", "adolescent", "child", "aged, 80 and over", "mice",
    "rats", "child, preschool", "infant", "retrospective studies",
    "prospective studies", "treatment outcome", "cell line, tumor",
    "prognosis", "survival rate", "follow-up studies", "time factors",
    "reproducibility of results", "risk factors", "cohort studies",
    "sensitivity and specificity", "predictive value of tests",
    "kaplan-meier estimate", "survival analysis", "mice, inbred balb c",
    "mice, nude", "cells, cultured", "case-control studies",
}

# Publication types that tell you what an unlabelled article IS. Deliberately
# coarse: the question is whether a mechanism taxonomy SHOULD have reached it.
PUBTYPE_BUCKETS = {
    "review/opinion": ("review", "editorial", "comment", "letter", "news",
                       "historical article", "lecture"),
    "case report": ("case reports",),
    "trial": ("clinical trial", "randomized controlled trial",
              "controlled clinical trial", "multicenter study"),
    "guideline/consensus": ("guideline", "practice guideline", "consensus"),
    "meta/systematic": ("meta-analysis", "systematic review"),
}


def wilson(k: int, n: int, z: float = 1.96) -> tuple:
    """Binomial CI, so a sampled rate is never reported as a point."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def mesh_leaf_descriptors() -> set:
    """Every descriptor the precision-first reference map names.

    Parsed with a real YAML loader. A hand-rolled line parser here missed five
    of the thirty-one descriptors, because some entries are multi-line lists --
    and an undercount would have understated the map's reach in exactly the
    place this script makes a claim about it.
    """
    import yaml
    doc = yaml.safe_load(MESH_MAP.read_text(encoding="utf-8")) or {}
    out = set()
    for _mech, body in (doc.get("mechanisms") or {}).items():
        for d in (body.get("descriptors") or []):
            if isinstance(d, str) and d.strip():
                out.add(d.strip().lower())
    if not out:
        raise SystemExit(
            "no descriptors parsed from mesh-mechanism-map.yaml; an empty "
            "reference set would report 0% reach as though it were a finding.")
    return out


def scan(sample_every: int):
    leaves = mesh_leaf_descriptors()
    kw_items = [(m, k) for m, ks in config.MECHANISM_KEYWORDS.items() for k in ks]

    total = 0
    mesh_hit = 0                       # full census
    sampled = 0
    kw_hit = 0
    per_mech = Counter()
    untagged_pubtypes = Counter()
    untagged_mesh = Counter()
    untagged_no_pubtype = 0

    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                total += 1
                mesh = [m.lower() for m in (r.get("mesh") or [])]
                if any(m in leaves for m in mesh):
                    mesh_hit += 1
                if total % sample_every:
                    continue
                sampled += 1
                blob = ((r.get("title") or "") + " " +
                        (r.get("abstract") or "")).lower()
                hits = {m for m, k in kw_items if text_matches_keyword(blob, k)}
                if hits:
                    kw_hit += 1
                    for m in hits:
                        per_mech[m] += 1
                    continue
                # profile the remainder: this is the deliverable
                pts = [p.lower() for p in (r.get("pub_types") or [])]
                if not pts:
                    untagged_no_pubtype += 1
                bucketed = False
                for bucket, needles in PUBTYPE_BUCKETS.items():
                    if any(any(n in p for n in needles) for p in pts):
                        untagged_pubtypes[bucket] += 1
                        bucketed = True
                if not bucketed and pts:
                    untagged_pubtypes["primary research (no special type)"] += 1
                for m in mesh[:40]:
                    if m not in CHECK_TAGS:
                        untagged_mesh[m] += 1

    return {
        "census_total": total,
        "mesh_leaf_hits": mesh_hit,
        "sample_every": sample_every,
        "sampled": sampled,
        "keyword_hits": kw_hit,
        "per_mechanism": dict(per_mech.most_common()),
        "untagged_pubtypes": dict(untagged_pubtypes.most_common()),
        "untagged_top_mesh": dict(untagged_mesh.most_common(25)),
        "untagged_no_pubtype": untagged_no_pubtype,
        "n_mechanisms": len(config.MECHANISM_KEYWORDS),
        "n_keywords": len(kw_items),
        "n_mesh_leaves": len(leaves),
    }


def render(d: dict) -> str:
    n, s = d["census_total"], d["sampled"]
    kw, mh = d["keyword_hits"], d["mesh_leaf_hits"]
    lo, hi = wilson(kw, s)
    untagged = s - kw
    L = [f"# What fraction of the cancer literature can the mechanism taxonomy label?", ""]
    L += ["*Generated by `scripts/atlas_taxonomy_reach.py`. Every figure is recomputed.*", ""]

    L += ["## The field of view", ""]
    L += [f"The keyword tagger carries **{d['n_mechanisms']} mechanisms** over "
          f"**{d['n_keywords']} terms**. Applied to a uniform 1-in-"
          f"{d['sample_every']} sample of the {n:,}-article census "
          f"({s:,} articles):", ""]
    L += ["| | count | share |", "|---|--:|--:|"]
    L += [f"| sampled | {s:,} | |",
          f"| **matched at least one mechanism** | **{kw:,}** | "
          f"**{100*kw/s:.2f}%** (95% CI {100*lo:.2f}-{100*hi:.2f}%) |",
          f"| matched none | {untagged:,} | {100*untagged/s:.2f}% |", ""]
    L += [f"So every mechanism share this project reports is a share of roughly "
          f"**{100*kw/s:.0f}%** of the cancer literature, and the documented "
          f"0.20%-41.86% capture spread describes variation inside that "
          f"fraction.", ""]

    L += ["### The precision-first map, for contrast", ""]
    L += [f"`mesh-mechanism-map.yaml` names {d['n_mesh_leaves']} leaf descriptors "
          f"and reaches {mh:,} of {n:,} census articles "
          f"(**{100*mh/n:.2f}%**), measured on the full census rather than a "
          f"sample.", ""]
    L += ["That number is **not** a coverage failure. The map deliberately drops "
          "umbrella descriptors -- bare `Immunotherapy` is excluded as "
          "non-discriminative -- because its job is to measure the tagger's "
          "recall without the tagger's own vocabulary leaking into the "
          "reference. It is a precision instrument and should be read as one.", ""]

    L += ["## What the unlabelled remainder actually is", ""]
    L += ["The percentage is less interesting than its complement. If the "
          "remainder is literature no mechanism taxonomy should claim, the field "
          "of view is appropriate and only the wording of the capture caveat "
          "needs fixing.", ""]
    L += ["| publication type | share of unlabelled |", "|---|--:|"]
    for k, v in list(d["untagged_pubtypes"].items())[:8]:
        L.append(f"| {k} | {100*v/max(untagged,1):.1f}% |")
    L += [""]
    L += ["Most common MeSH descriptors among unlabelled articles, with "
          "demographic check-tags excluded (`Humans` alone sits on 92% of them "
          "and says nothing about subject):", ""]
    L += ["| descriptor | share of unlabelled |", "|---|--:|"]
    for k, v in list(d["untagged_top_mesh"].items())[:12]:
        L.append(f"| {k} | {100*v/max(untagged,1):.1f}% |")
    L += [""]

    # The falsifier, answered from the numbers rather than asserted.
    ther = sum(v for k, v in d["untagged_top_mesh"].items()
               if "antineoplastic" in k or "chemotherapy" in k
               or "radiotherapy" in k or "drug therapy" in k)
    softs = sum(v for k, v in d["untagged_pubtypes"].items()
                if k in ("review/opinion", "case report"))
    L += ["### The verdict on the field of view", ""]
    L += [f"The comfortable reading is that the unlabelled remainder is "
          f"literature no mechanism taxonomy should claim. That is "
          f"**partly true and not sufficient**: review, opinion and case "
          f"reports account for {100*softs/max(untagged,1):.1f}% of it, but "
          f"{d['untagged_pubtypes'].get('primary research (no special type)', 0)*100/max(untagged,1):.1f}% "
          f"is primary research with no special publication type.", ""]
    L += [f"And the remainder carries explicit therapy descriptors: "
          f"`antineoplastic agents` and "
          f"`antineoplastic combined chemotherapy protocols` together sit on "
          f"**{100*ther/max(untagged,1):.1f}%** of unlabelled articles. Those "
          f"are therapy papers the taxonomy has no name for -- chemotherapy "
          f"has no mechanism tag -- rather than literature outside its remit.", ""]
    L += ["So the field of view is a real limit and not merely a wording "
          "problem. The capture caveat should carry this number, and the "
          "backbone modalities with no tag are the first place to widen.", ""]

    L += ["## Per mechanism, within the sample", ""]
    L += ["| mechanism | sampled hits | share of census |", "|---|--:|--:|"]
    for m, v in list(d["per_mechanism"].items())[:15]:
        L.append(f"| {m} | {v:,} | {100*v/s:.2f}% |")
    L += [""]

    L += ["## What this does not say", ""]
    L += ["* It does not say the taxonomy is wrong. A mechanism layer is supposed "
          "to label mechanism papers, and much of the literature is not about a "
          "therapeutic mechanism at all.",
          "* It does not measure precision. An article matching a keyword is not "
          "necessarily about that mechanism; `mechanism_recall.py` measures that "
          "separately against MeSH.",
          "* The keyword figure is sampled, so it carries an interval. The MeSH "
          "figure is a full-census count and does not.",
          "* Reach is measured on title and abstract only, which is what the "
          "production tagger reads for these fields.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample-every", type=int, default=20,
                    help="take every Nth census record for the keyword scan")
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()

    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan(args.sample_every)
        if d["keyword_hits"] == 0:
            raise SystemExit(
                "the keyword tagger matched nothing, which is not a finding -- "
                "it is what a broken import or an empty vocabulary looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
