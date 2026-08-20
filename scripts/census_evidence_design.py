#!/usr/bin/env python3
"""What study design each census article reports, on labels this project did not assign (#RETIRE-FROZEN).

WHY THIS REPLACES THE EVIDENCE LADDER
-------------------------------------
The manuscript's seven-tier evidence ladder was produced by a keyword tagger
reading the full text of 4,830 articles, and scored against a 100-article gold
set. Neither input describes the census: there is no full text for three
quarters of it, and a gold set drawn from a retrieval says nothing about the
literature that retrieval did not select.

What the census carries instead is BETTER on the axis that matters most and
WORSE on another, and this page states both rather than trading one for the
other quietly:

  BETTER -- NLM's own `PublicationType` is assigned by professional indexers,
  not by this project. A trial labelled `Clinical Trial, Phase III` here is
  labelled that by the National Library of Medicine. There is no tagger to
  score, so no accuracy figure is owed and none is claimed.

  WORSE -- publication types are not evidence tiers. NLM has no concept of
  preclinical in-vivo versus in-vitro, so that distinction is approximated from
  MeSH check tags, which answer "what organism or material appears in this
  study" rather than "what kind of study is this". The two are close but not
  the same question, and the classes are named for what they measure.

COVERAGE IS THE HEADLINE, NOT A CAVEAT
--------------------------------------
A large share of the census carries `Journal Article` and nothing more
discriminating. That share is reported first, because a design distribution
computed over the classifiable remainder and presented as a distribution over
the literature would be the same error this project has made before: quoting a
number over one population against another.

Usage:
    python scripts/census_evidence_design.py
    python scripts/census_evidence_design.py --render-only
"""

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RECORDS = PROJECT_ROOT / "corpus" / "atlas" / "records"
OUT_MD = PROJECT_ROOT / "analysis" / "census-evidence-design.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "census-evidence-design.json"

# NLM publication types. Trial-specific only: `Multicenter Study` and
# `Comparative Study` are excluded because they apply to observational work
# too -- the same exclusion `atlas_landscape.CLINICAL_TYPES` already makes.
TRIAL = {
    "Clinical Trial", "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial",
}
PHASED = {"Clinical Trial, Phase I", "Clinical Trial, Phase II",
          "Clinical Trial, Phase III", "Clinical Trial, Phase IV"}
# Human study that is not a trial.
CLINICAL_OTHER = {"Observational Study", "Case Reports", "Clinical Study",
                  "Twin Study", "Validation Study"}
# Literature ABOUT studies rather than reporting one. `atlas_evidence_check`
# measured 73.8% of the untagged pile as this, and leaving it untagged was the
# tagger working rather than failing.
NON_PRIMARY = {"Review", "Systematic Review", "Meta-Analysis", "Editorial",
               "Comment", "Letter", "News", "Historical Article",
               "Practice Guideline", "Guideline", "Consensus Development Conference",
               "Scoping Review", "Narrative"}
# `Journal Article` and the funding tags carry no design information at all.
UNINFORMATIVE = {"Journal Article", "English Abstract", "Multicenter Study",
                 "Comparative Study", "Evaluation Study"}

# MeSH check tags. These answer "what organism or material is in this study",
# NOT "what kind of study is this", and the class names say so.
IN_VIVO = {"xenograft model antitumor assays", "disease models, animal",
           "mice, nude", "mice, inbred balb c", "mice, inbred c57bl",
           "mice, transgenic", "mice, scid", "rats, nude", "neoplasm transplantation"}
ANIMAL_ANY = {"animals", "mice", "rats", "zebrafish", "dogs", "rabbits", "swine"}
IN_VITRO = {"cell line, tumor", "in vitro techniques", "cells, cultured",
            "tumor cells, cultured", "spheroids, cellular", "cell culture techniques"}
HUMAN = "humans"


def classify(pub_types, mesh) -> str:
    """One design class per record, most specific first.

    Precedence is deliberate and is the judgement this file makes: a record
    carrying both a trial type and a mouse check tag is a trial that also used
    a model, not a preclinical study.
    """
    pt = set(pub_types or [])
    if pt & TRIAL:
        return "trial"
    if pt & NON_PRIMARY:
        return "non-primary"
    if pt & CLINICAL_OTHER:
        return "clinical-other"
    ms = {m.lower() for m in (mesh or [])}
    if ms & IN_VIVO:
        return "animal-model"
    if ms & IN_VITRO:
        # a study with both an animal model and cultured cells was caught above
        return "cell-culture"
    if ms & ANIMAL_ANY and HUMAN not in ms:
        return "animal-other"
    return "undetermined"


ORDER = ["trial", "clinical-other", "animal-model", "cell-culture",
         "animal-other", "non-primary", "undetermined"]

WHAT_IT_MEASURES = {
    "trial": "NLM assigned a trial publication type",
    "clinical-other": "NLM assigned a non-trial human-study type",
    "animal-model": "carries a tumour-model or animal-disease-model descriptor",
    "cell-culture": "carries a cultured-cell or in-vitro descriptor and no animal model",
    "animal-other": "carries an animal check tag and no `Humans` tag",
    "non-primary": "review, editorial, comment or guideline -- reports no study of its own",
    "undetermined": "no publication type or descriptor that discriminates design",
}


def scan() -> dict:
    cls = Counter()
    phase = Counter()
    by_year_trial = Counter()
    by_year_total = Counter()
    bare = 0
    n = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                pt = r.get("pub_types") or []
                c = classify(pt, r.get("mesh"))
                cls[c] += 1
                if set(pt) <= UNINFORMATIVE or not pt:
                    bare += 1
                for p in set(pt) & PHASED:
                    phase[p] += 1
                y = r.get("year")
                if isinstance(y, int) and 1975 <= y <= 2026:
                    by_year_total[y] += 1
                    if c == "trial":
                        by_year_trial[y] += 1
    classifiable = n - cls["undetermined"]
    return {
        "census": n,
        "classes": {k: cls.get(k, 0) for k in ORDER},
        "classifiable": classifiable,
        "bare_or_uninformative_pub_types": bare,
        "phased_trials": dict(phase.most_common()),
        "trial_share_by_year": {str(y): [by_year_trial.get(y, 0), by_year_total[y]]
                                for y in sorted(by_year_total)},
    }


def render(d: dict) -> str:
    n = d["census"]
    c = d["classes"]
    L = ["# What study design the census reports", ""]
    L += [f"*Generated by `scripts/census_evidence_design.py` over {n:,} "
          f"MeSH-indexed census articles. Every label here is NLM's, not this "
          f"project's: publication types are assigned by professional indexers "
          f"and MeSH check tags by the same. There is no tagger of ours to "
          f"score, so no accuracy figure is owed and none is claimed.*", ""]

    und = c["undetermined"]
    L += ["## Coverage first", ""]
    L += [f"**{und:,} of {n:,} ({100*und/n:.1f}%) carry nothing that "
          f"discriminates study design** -- no trial or review publication "
          f"type, no model or culture descriptor. A distribution computed over "
          f"the remaining {d['classifiable']:,} and presented as a "
          f"distribution over the literature would be a number quoted against "
          f"a population it was not measured on, which is a mistake this "
          f"project has made before. Both denominators are given in the table.",
          ""]
    L += [f"{d['bare_or_uninformative_pub_types']:,} "
          f"({100*d['bare_or_uninformative_pub_types']/n:.1f}%) carry only "
          f"`Journal Article` or funding and language tags, which say nothing "
          f"about design at all.", ""]

    L += ["| class | what the label actually means | records | of census | "
          "of classifiable |", "|---|---|--:|--:|--:|"]
    for k in ORDER:
        v = c[k]
        cl = f"{100*v/d['classifiable']:.1f}%" if k != "undetermined" else "-"
        L.append(f"| `{k}` | {WHAT_IT_MEASURES[k]} | {v:,} | "
                 f"{100*v/n:.1f}% | {cl} |")
    L += [""]

    ph = d.get("phased_trials") or {}
    if ph:
        L += ["## The trial column is the strong one", ""]
        L += [f"Of the {c['trial']:,} trials, these carry an explicit phase, "
              f"assigned by NLM:", ""]
        L += ["| phase | records |", "|---|--:|"]
        for k, v in sorted(ph.items()):
            L.append(f"| {k} | {v:,} |")
        L += ["", "This is the one place the census is strictly better than the "
              "retrieval it replaces: a phase label here is the National "
              "Library of Medicine's judgement, where the superseded ladder "
              "inferred phase from keywords in article text and was measured "
              "at 46% exact-label accuracy against its own gold set.", ""]

    L += ["## What this does NOT measure", ""]
    L += ["* **It is not an evidence hierarchy.** The classes are study "
          "designs as NLM labels them. Ordering them by strength is a "
          "judgement this page does not make.",
          "* **`animal-model` and `cell-culture` are check-tag inferences.** "
          "MeSH check tags answer what organism or material appears in a "
          "study, not what kind of study it is. A clinical paper that "
          "mentions a cell line in passing can carry the descriptor. These "
          "two classes are the weakest rows here and are named for the "
          "descriptor rather than for a tier.",
          "* **Nothing distinguishes theoretical or computational work.** The "
          "superseded ladder had a tier for it; no publication type or "
          "descriptor identifies it reliably, so it is not approximated.",
          "* **Recent literature is undercounted in every class.** MeSH "
          "indexing lags publication, and `analysis/atlas-recent-window.md` "
          "measures how much.",
          "* **The 783,271 text-recovered census records are excluded** -- "
          "they carry no MeSH and no publication types, so no label of this "
          "kind exists for them. The denominator here is the "
          f"{n:,} indexed stream, not the {n + 783271:,} full census.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan()
        if d["classes"]["trial"] == 0:
            raise SystemExit(
                "no trial publication types matched, which is not a finding -- "
                "it is what a field-name or case mismatch looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
        d = json.loads(OUT_JSON.read_text())
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}")
    for k in ORDER:
        print(f"  {k:16s} {d['classes'][k]:>9,}  {100*d['classes'][k]/d['census']:>5.1f}%")


if __name__ == "__main__":
    main()
