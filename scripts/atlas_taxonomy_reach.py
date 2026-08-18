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
from tag_articles import (  # noqa: E402
    text_matches_keyword, match_mechanisms, normalize_text)

ATLAS = PROJECT_ROOT / "corpus" / "atlas"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-taxonomy-reach.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-taxonomy-reach.json"
MESH_MAP = PROJECT_ROOT / "analysis" / "mesh-mechanism-map.yaml"

# MeSH "check tags": demographic descriptors applied to almost every article.
# Without excluding them the top-descriptor table reads humans 92%, female 46%,
# male 37% -- true, and completely uninformative about what the literature is
# ABOUT, which is the only question this section exists to answer.
# SPLIT, because the report said "demographic check-tags excluded" while the
# set also removed STUDY-DESIGN descriptors -- and those are exactly what
# answers the section's question about what the unlabelled remainder IS.
# `cell line, tumor`, `retrospective studies`, `prognosis`, `treatment
# outcome` and `risk factors` were all silently dropped from a table headed
# "most common descriptors among unlabelled articles".
DEMOGRAPHIC_TAGS = {
    "humans", "female", "male", "animals", "adult", "aged", "middle aged",
    "young adult", "adolescent", "child", "aged, 80 and over", "mice",
    "rats", "child, preschool", "infant", "mice, inbred balb c",
    "mice, nude", "cells, cultured",
}
STUDY_DESIGN_TAGS = {
    "retrospective studies", "prospective studies", "treatment outcome",
    "cell line, tumor", "prognosis", "survival rate", "follow-up studies",
    "time factors", "reproducibility of results", "risk factors",
    "cohort studies", "sensitivity and specificity",
    "predictive value of tests", "kaplan-meier estimate",
    "survival analysis", "case-control studies",
}
CHECK_TAGS = DEMOGRAPHIC_TAGS | STUDY_DESIGN_TAGS

# The descriptors the 'therapy papers the taxonomy cannot name'
# sentence is about. Named here so the measured union and the
# rendered sentence cannot drift apart.
THERAPY_DESCRIPTORS = {
    "antineoplastic agents",
    "antineoplastic combined chemotherapy protocols",
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
    prod_hit = 0
    prod_ta_hit = 0
    per_mech = Counter()
    untagged_pubtypes = Counter()
    untagged_mesh = Counter()
    untagged_design = Counter()
    # UNIONS, because two overlapping buckets summed is not a share of
    # anything. `antineoplastic agents` and the chemotherapy-protocols
    # descriptor co-occur, and so do review/opinion and case report.
    untagged_therapy_union = 0
    untagged_soft_union = 0
    untagged_no_pubtype = 0
    multi_bucket = 0

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
                # THE PRODUCTION MATCHER, measured beside the raw loop. The
                # published 5.98% came from looping the low-level
                # `text_matches_keyword`; production calls `match_mechanisms`,
                # which opens with a cancer-context gate and reads MeSH as
                # well as title and abstract. The two are not the same field
                # of view and the report claimed they were.
                title_t = normalize_text(r.get("title") or "")
                if match_mechanisms(blob, title_t):
                    prod_ta_hit += 1
                prod_text = normalize_text(
                    " ".join([r.get("title") or "", " ".join(mesh),
                              r.get("abstract") or ""]))
                if match_mechanisms(prod_text, title_t):
                    prod_hit += 1
                if hits:
                    kw_hit += 1
                    for m in hits:
                        per_mech[m] += 1
                    continue
                # profile the remainder: this is the deliverable
                pts = [p.lower() for p in (r.get("pub_types") or [])]
                if not pts:
                    untagged_no_pubtype += 1
                # FIRST MATCH WINS, so the buckets partition. Without the
                # break a record landed in several and the table -- presented
                # as a breakdown of the remainder -- summed past 100%.
                bucketed = False
                for bucket, needles in PUBTYPE_BUCKETS.items():
                    if any(any(n in p for n in needles) for p in pts):
                        untagged_pubtypes[bucket] += 1
                        bucketed = True
                        break
                if bucketed:
                    # how often the old multi-count would have fired, kept so
                    # the correction is visible rather than silent
                    extra = sum(
                        1 for b, ns in PUBTYPE_BUCKETS.items()
                        if b != bucket and any(any(n in p for n in ns) for p in pts))
                    if extra:
                        multi_bucket += 1
                if not bucketed and pts:
                    untagged_pubtypes["primary research (no special type)"] += 1
                # EXACTLY the two descriptors the report's sentence names.
                # A wider stem set here would make the number and the prose
                # describe different things -- the sum-versus-union defect in
                # another form.
                if THERAPY_DESCRIPTORS.intersection(mesh):
                    untagged_therapy_union += 1
                if any(any(n in pp for n in
                           PUBTYPE_BUCKETS.get("review/opinion", []) +
                           PUBTYPE_BUCKETS.get("case report", []))
                       for pp in pts):
                    untagged_soft_union += 1
                for m in mesh[:40]:
                    if m in DEMOGRAPHIC_TAGS:
                        continue
                    if m in STUDY_DESIGN_TAGS:
                        untagged_design[m] += 1
                        continue
                    untagged_mesh[m] += 1

    return {
        "census_total": total,
        "mesh_leaf_hits": mesh_hit,
        "sample_every": sample_every,
        "sampled": sampled,
        "keyword_hits": kw_hit,
        "production_hits": prod_hit,
        "production_title_abstract_hits": prod_ta_hit,
        # ORDERED LISTS OF PAIRS. A dict here is reordered by
        # `json.dumps(sort_keys=True)`, so the committed artifact was
        # ALPHABETICAL and `--render-only` produced a different
        # document from the documented command. The same defect is
        # already recorded for atlas_untagged_partner.py.
        "per_mechanism": [[k, v] for k, v in per_mech.most_common()],
        "untagged_pubtypes": [[k, v] for k, v in untagged_pubtypes.most_common()],
        "untagged_top_mesh": [[k, v] for k, v in untagged_mesh.most_common(25)],
        "untagged_top_study_design": [[k, v] for k, v in untagged_design.most_common(12)],
        "untagged_no_pubtype": untagged_no_pubtype,
        "untagged_multi_bucket": multi_bucket,
        "untagged_therapy_union": untagged_therapy_union,
        "untagged_soft_union": untagged_soft_union,
        "n_mechanisms": len(config.MECHANISM_KEYWORDS),
        "n_keywords": len(kw_items),
        "n_mesh_leaves": len(leaves),
    }


def _pairs(v):
    """Ordered (key, value) pairs from either shape.

    The artifact stores ORDERED LISTS now, because `json.dumps(sort_keys=True)`
    reordered the dicts and made `--render-only` produce an alphabetical
    document while the documented command produced a ranked one. A dict read
    from an older artifact is re-sorted by value here rather than trusted, so
    the two paths cannot diverge again.
    """
    if isinstance(v, dict):
        return sorted(v.items(), key=lambda kv: -kv[1])
    return [(k, n) for k, n in v]


def _get(v, key, default=0):
    for k, n in _pairs(v):
        if k == key:
            return n
    return default


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
          f"| matched none | {untagged:,} | {100*untagged/s:.2f}% |"]
    if not d.get("production_hits"):
        L += [""]
    prod = d.get("production_hits")
    prod_ta = d.get("production_title_abstract_hits")
    if prod:
        plo, phi = wilson(prod, s)
        L += [f"| **the PRODUCTION matcher** | **{prod:,}** | "
              f"**{100*prod/s:.2f}%** (95% CI {100*plo:.2f}-{100*phi:.2f}%) |",
              ""]
        L += [f"The `matched at least one mechanism` row is a raw loop over "
              f"`text_matches_keyword` on title and abstract. **Production is "
              f"the {100*prod/s:.2f}% row**: `match_mechanisms` opens with a "
              f"cancer-context gate and reads the MeSH descriptors as well as "
              f"the title and abstract. On the same articles the gate alone "
              f"costs {100*(prod_ta - kw)/s:+.2f} points and the MeSH channel "
              f"adds {100*(prod - prod_ta)/s:+.2f}. An earlier version "
              f"published the raw figure as the production field of view.", ""]
    L += [f"So every mechanism share this project reports is a share of "
          f"**{100*(prod or kw)/s:.1f}%** of the cancer literature.", ""]
    # NOT the capture spread. That is computed entirely by the MeSH map, whose
    # reach this same page reports separately -- attaching it here conflates
    # the two instruments the docstring forbids conflating.
    L += ["The documented 0.20%-41.86% per-mechanism capture spread is NOT "
          "variation inside this fraction: it is computed by the MeSH leaf "
          "map on both sides, a different instrument whose reach is reported "
          "below. An earlier version of this sentence attached the spread to "
          "the keyword reach, which is the conflation this script's own "
          "docstring forbids.", ""]

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
    for k, v in _pairs(d["untagged_pubtypes"])[:8]:
        L.append(f"| {k} | {100*v/max(untagged,1):.1f}% |")
    L += [""]
    L += ["Most common MeSH descriptors among unlabelled articles, with "
          "demographic check-tags excluded (`Humans` alone sits on 92% of them "
          "and says nothing about subject):", ""]
    L += ["| descriptor | share of unlabelled |", "|---|--:|"]
    for k, v in _pairs(d["untagged_top_mesh"])[:12]:
        L.append(f"| {k} | {100*v/max(untagged,1):.1f}% |")
    L += [""]

    # The falsifier, answered from the numbers rather than asserted.
    # UNION, measured in the scan. Summing the two descriptor rows
    # double-counts every article carrying both.
    ther = d.get("untagged_therapy_union")
    if ther is None:  # older artifact
        ther = sum(v for k, v in _pairs(d["untagged_top_mesh"])
                   if "antineoplastic" in k or "chemotherapy" in k
                   or "radiotherapy" in k or "drug therapy" in k)
    softs = d.get("untagged_soft_union")
    if softs is None:
        softs = sum(v for k, v in _pairs(d["untagged_pubtypes"])
                    if k in ("review/opinion", "case report"))
    L += ["### The verdict on the field of view", ""]
    L += [f"The comfortable reading is that the unlabelled remainder is "
          f"literature no mechanism taxonomy should claim. That is "
          f"**partly true and not sufficient**: review, opinion and case "
          f"reports account for {100*softs/max(untagged,1):.1f}% of it, but "
          f"{_get(d['untagged_pubtypes'], 'primary research (no special type)')*100/max(untagged,1):.1f}% "
          f"is primary research with no special publication type.", ""]
    L += [f"And the remainder carries explicit therapy descriptors: "
          f"`antineoplastic agents` and "
          f"`antineoplastic combined chemotherapy protocols` sit on "
          f"**{100*ther/max(untagged,1):.1f}%** of unlabelled articles. Those "
          f"are therapy papers the taxonomy has no name for -- chemotherapy "
          f"has no mechanism tag -- rather than literature outside its remit. "
          f"That is the UNION of the therapy descriptors, not the sum of two "
          f"overlapping rows, which an earlier version published.", ""]
    L += ["So the field of view is a real limit and not merely a wording "
          "problem. The capture caveat should carry this number, and the "
          "backbone modalities with no tag are the first place to widen.", ""]

    L += ["## Per mechanism, within the sample", ""]
    L += [f"The {min(15, len(_pairs(d['per_mechanism'])))} largest of "
          f"{len(_pairs(d['per_mechanism']))} mechanisms with any hit, "
          f"BY COUNT. An earlier version sliced an alphabetically "
          f"reordered dict and omitted `nanoparticle`, the "
          f"second-largest.", ""]
    L += ["| mechanism | sampled hits | share of census |", "|---|--:|--:|"]
    # top 15 BY COUNT, and the truncation is stated: an alphabetical
    # slice omitted `nanoparticle`, the second-largest mechanism.
    _pm = _pairs(d["per_mechanism"])
    for m, v in _pm[:15]:
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
          "* Two reaches are reported: a raw keyword loop over title and "
          "abstract, and the production matcher, which also reads MeSH and "
          "applies a cancer-context gate. An earlier version reported only "
          "the first and described it as what production reads.",
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
