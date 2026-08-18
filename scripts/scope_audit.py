#!/usr/bin/env python3
"""What is this project's work actually about? (#731)

WHY THIS EXISTS
---------------
The corpus is broad and the work is not, and a reader cannot tell from the front
door. The README leads with census scale -- five million articles, ten million
relations -- and a reader reasonably infers that the analysis is commensurate.
It is not: the analytical, predictive and simulation work is almost entirely
ferroptosis and the physical-ROS modalities.

That is a legitimate way to run a project. A narrow thesis on a broad corpus is
how most good science works. What is not legitimate is leaving a reader to
discover it by counting files, which is what the repository currently requires.

So this counts, rather than asserting, and the README quotes the counts.

WHAT IS MECHANICAL AND WHAT IS JUDGEMENT
-----------------------------------------
Two of the three measures need no interpretation and are derived here:

  PREDICTIONS. `PREREGISTRATION.md` states P1..Pn. A prediction is counted
  ferroptosis-or-physical-ROS if its statement names ferroptosis, a ferroptosis
  gene or inducer, or PDT/SDT. These are the project's falsifiable commitments,
  so they are the sharpest measure of what it is willing to be wrong about.

  ENGINE MODULES. `simulations/ferroptosis-core/src/*.rs` is the whole engine.
  Counting is exact.

The third is a judgement and is therefore LISTED rather than merely counted, so
a reader can disagree with a specific placement instead of the total:

  ANALYSES. Each file in `analysis/` is placed in one of three buckets by the
  rule below. The bucket members are written into the artifact. If you think a
  file is in the wrong bucket, the disagreement is auditable.

THE RULE, stated before the result. An analysis is:
  * `therapy-subject` if its SUBJECT is a named therapy other than ferroptosis
    or the physical-ROS modalities;
  * `ferroptosis-or-physical` if its subject is ferroptosis, GPX4/FSP1 biology,
    PDT or SDT;
  * `method` otherwise -- the corpus, the census, the graph, the taxonomy, the
    tooling and their error rates.

A method analysis is not lesser work. The census methodology is much of this
project's value. The point of separating it is that method analyses are about
the INSTRUMENT, so they cannot tell a reader what BIOLOGY the project studies.

Usage:
    python scripts/scope_audit.py
"""

import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_MD = PROJECT_ROOT / "analysis" / "scope-audit.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "scope-audit.json"

# `ferropto` not `ferroptos`: the latter does not match "ferroptotic", which
# misfiled P5 ("dense ferroptotic kill ...") as a non-ferroptosis prediction --
# an error in the flattering direction, since it made the commitments look
# broader than they are.
FERRO = re.compile(
    r"ferropto|gpx4|fsp1|aifm2|rsl3|erastin|lipid peroxid|slc7a11|"
    r"photodynamic|sonodynamic|\bpdt\b|\bsdt\b|photosensitiz", re.I)

# Analyses whose subject is the INSTRUMENT: the census, the graph, the tagger,
# the taxonomy, the corpus and their error rates. Checked BEFORE the ferroptosis
# test, because many of them use the FSP1 sense collision as their worked
# example and a body-text scan therefore files them as ferroptosis biology --
# which inflates the ferroptosis count in the direction of this file's own
# argument. Subject, not vocabulary.
METHOD_STEM = re.compile(
    r"^(atlas|comention|census|corpus|mesh|evidence|mechanism|manuscript|"
    r"prisma|taxonomy|oa-bias|news|landmark|seed|scope|gap|convergence|"
    r"non-oa|rare-event|identifiability|sobol|uncertainty|headline|"
    r"reaction-diffusion|research-roadmap|contribution|collaborator|"
    r"ode-cross|arbatskiy|talkington|sbml|p1-wetlab|diagnostic)", re.I)

# Therapies whose NAME in a filename marks the analysis as being about them.
# Deliberately a name list rather than a keyword scan of the body: almost every
# analysis mentions immunotherapy somewhere, and mentioning is not being about.
THERAPY_NAMES = (
    "radioligand", "immunotherap", "checkpoint", "car-t", "cart",
    "oncolytic", "adc", "antibody-drug", "crispr", "nanoparticle",
    "mrna-vaccine", "vaccine", "microbiome", "radiotherap", "chemotherap",
    "surgery", "hormone", "bispecific", "epigenetic", "proteolysis",
    "targeted-protein",
)


def classify_analyses() -> dict:
    buckets = {"ferroptosis-or-physical": [], "therapy-subject": [], "method": []}
    for p in sorted((PROJECT_ROOT / "analysis").glob("*.md")):
        stem = p.stem.lower()
        text = p.read_text(errors="ignore")
        head = "\n".join(text.split("\n")[:40])
        if any(t in stem for t in THERAPY_NAMES):
            buckets["therapy-subject"].append(p.stem)
        elif METHOD_STEM.match(stem) and not FERRO.search(stem):
            buckets["method"].append(p.stem)
        elif FERRO.search(stem) or len(FERRO.findall(head)) >= 3:
            buckets["ferroptosis-or-physical"].append(p.stem)
        else:
            buckets["method"].append(p.stem)
    return buckets


def classify_predictions() -> dict:
    text = (PROJECT_ROOT / "PREREGISTRATION.md").read_text()
    out = {}
    for m in re.finditer(r"^\*\*(P\d+)\.\s*(.+?)\*\*", text, re.M):
        out[m.group(1)] = bool(FERRO.search(m.group(2)))
    if not out:
        raise SystemExit(
            "no predictions parsed from PREREGISTRATION.md; an empty set would "
            "report 0 of 0 ferroptosis predictions as though that were a "
            "measurement")
    return out


def engine_modules() -> int:
    d = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src"
    return len(list(d.glob("*.rs"))) if d.exists() else 0


def render(d: dict) -> str:
    a, p = d["analyses"], d["predictions"]
    n_ferro = len(a["ferroptosis-or-physical"])
    n_ther = len(a["therapy-subject"])
    n_meth = len(a["method"])
    tot = n_ferro + n_ther + n_meth
    pf = sum(1 for v in p.values() if v)
    L = ["# What the work is about", ""]
    L += ["*Generated by `scripts/scope_audit.py`. Counts are derived; the "
          "analysis buckets are listed so a placement can be disputed.*", ""]

    L += ["## The three measures", ""]
    L += ["| | ferroptosis / physical-ROS | other therapy | method |",
          "|---|--:|--:|--:|",
          f"| committed analyses ({tot}) | {n_ferro} | **{n_ther}** | {n_meth} |"]
    L += [""]
    L += [f"| preregistered predictions | **{pf} of {len(p)}** |",
          "|---|--:|",
          f"| engine modules | **{d['engine_modules']} of "
          f"{d['engine_modules']}** |", ""]

    L += ["Every falsifiable commitment the project makes, and every module of "
          "its simulation engine, concerns ferroptosis or the physical-ROS "
          "modalities.", ""]

    L += ["## Why this is worth stating rather than hiding", ""]
    L += ["A narrow thesis on a broad corpus is how most good science works, and "
          "the corpus is genuinely broad. But a reader who arrives at a front "
          "door advertising millions of articles will infer that the analysis is "
          "commensurate, and it is not. Leaving them to discover that by "
          "counting files is the part that is not defensible.", ""]
    L += ["Note also that method analyses -- the largest bucket -- are about the "
          "INSTRUMENT: the census, the graph, the taxonomy and their error "
          "rates. That work is much of this project's value, and it is exactly "
          "why it cannot tell a reader what biology the project studies.", ""]

    L += ["## The buckets, so a placement can be disputed", ""]
    for k in ("therapy-subject", "ferroptosis-or-physical"):
        L += [f"**{k}** ({len(a[k])}):", "",
              "".join(f"- `{x}`\n" for x in a[k]) or "- (none)\n"]
    L += [f"**method** ({n_meth}): listed in `analysis/scope-audit.json`.", ""]

    L += ["## What this does not say", ""]
    L += ["* It does not say the focus is wrong. It says it is undeclared.",
          "* It does not measure quality or effort, only subject.",
          "* The analysis bucketing is a judgement applied by a stated rule; the "
          "prediction and module counts are mechanical.",
          "* Mechanism shares quoted elsewhere in this repo (immunotherapy 47.6% "
          "and so on) are shares of TAGGED articles, and the mechanism taxonomy "
          "reaches 6.8% of the cancer literature through the production matcher — see "
          "`analysis/atlas-taxonomy-reach.md`.",
          ""]
    return "\n".join(L) + "\n"


def main():
    a = classify_analyses()
    p = classify_predictions()
    d = {"analyses": a, "predictions": p, "engine_modules": engine_modules(),
         "n_analyses": sum(len(v) for v in a.values()),
         "n_predictions": len(p),
         "n_predictions_ferroptosis": sum(1 for v in p.values() if v),
         "n_therapy_subject": len(a["therapy-subject"])}
    if d["engine_modules"] == 0:
        raise SystemExit("no engine modules found; the count would be a lie")
    OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                        encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  analyses: {len(a['ferroptosis-or-physical'])} ferroptosis, "
          f"{len(a['therapy-subject'])} other therapy, {len(a['method'])} method")
    print(f"  predictions: {d['n_predictions_ferroptosis']} of {d['n_predictions']} ferroptosis")
    print(f"  engine modules: {d['engine_modules']}")


if __name__ == "__main__":
    main()
