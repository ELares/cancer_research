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
# ONE vocabulary for both admission rules, built from THERAPY_NAMES so the
# body route and the filename route cannot differ in what they look for.
# An earlier version used a separate 14-term list against the filename
# route's 21, so the published contrast moved two things at once.
def _therapy_re():
    import re as _re
    return _re.compile("|".join(_re.escape(x) for x in
                                sorted(THERAPY_NAMES, key=len, reverse=True)),
                       _re.I)

THERAPY_NAMES = (
    "radioligand", "immunotherap", "checkpoint", "car-t", "cart",
    "oncolytic", "adc", "antibody-drug", "crispr", "nanoparticle",
    "mrna-vaccine", "vaccine", "microbiome", "radiotherap", "chemotherap",
    "surgery", "hormone", "bispecific", "epigenetic", "proteolysis",
    "targeted-protein",
)


def classify_analyses() -> dict:
    buckets = {"ferroptosis-or-physical": [], "therapy-subject": [], "method": []}
    body_therapy = []
    for p in sorted((PROJECT_ROOT / "analysis").glob("*.md")):
        stem = p.stem.lower()
        text = p.read_text(errors="ignore")
        # BLOCKQUOTE BANNERS ARE NOT CONTENT. The window is 40 lines of a
        # document's own prose; a provenance banner prepended to the top of a
        # page pushes real lines out of it and moves this audit's numbers
        # without any subject changing. That happened a third time when five
        # analysis pages gained retired-corpus banners: therapy-by-body fell
        # 46 -> 44 and the ferroptosis overlap 7 -> 5, purely from displacement.
        # The page's own prose already warns the figure has moved on this
        # window rather than on content, so the window now skips them.
        lines = [ln for ln in text.split("\n") if not ln.lstrip().startswith(">")]
        head = "\n".join(lines[:40])
        # THE TWO BUCKETS USE DIFFERENT ADMISSION RULES, and that is the
        # finding rather than a bug to paper over. Therapy is FILENAME ONLY;
        # ferroptosis is filename OR body. So the therapy count moves on a
        # rename with contents unchanged, and an empty file with a therapy
        # word in its name is filed as therapy-subject.
        #
        # The naive symmetric fix -- giving therapy the same >=3-head-hits
        # body route -- is NOT the answer: most files it would admit are
        # genuinely instrument analyses that merely mention a therapy as a
        # worked example, so it over-claims in the opposite direction. Both
        # numbers are reported, with what each rule is.
        if len(_therapy_re().findall(head)) >= 3 and not FERRO.search(stem):
            body_therapy.append(p.stem)
        if any(t in stem for t in THERAPY_NAMES):
            buckets["therapy-subject"].append(p.stem)
        elif METHOD_STEM.match(stem) and not FERRO.search(stem):
            buckets["method"].append(p.stem)
        elif FERRO.search(stem) or len(FERRO.findall(head)) >= 3:
            buckets["ferroptosis-or-physical"].append(p.stem)
        else:
            buckets["method"].append(p.stem)
    return buckets, sorted(body_therapy)


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


def mechanism_denominators() -> dict:
    """The three denominators a mechanism share can be quoted against.

    They are not interchangeable, and both this page and the README said the
    shares were "of TAGGED articles" when the published figure is the share of
    the CORPUS. Derived here so the sentence cannot drift from the numbers.
    """
    idx = PROJECT_ROOT / "corpus" / "INDEX.jsonl"
    if not idx.exists():
        return {}
    n_rec = tagged = tags = top_n = 0
    top = None
    counts = {}
    with idx.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except ValueError:
                continue
            n_rec += 1
            ms = r.get("mechanisms") or []
            if ms:
                tagged += 1
            tags += len(ms)
            for m in ms:
                counts[m] = counts.get(m, 0) + 1
    if not n_rec or not counts:
        return {}
    top, top_n = max(counts.items(), key=lambda kv: kv[1])
    return {"corpus_records": n_rec, "tagged": tagged, "tags": tags,
            "top_mechanism": top, "top_n": top_n,
            "share_of_corpus": top_n / n_rec,
            "share_of_tagged": top_n / tagged if tagged else None,
            "share_of_tags": top_n / tags if tags else None}


def engine_modules() -> dict:
    """How many engine modules actually CONCERN ferroptosis or physical ROS.

    An earlier version returned a file COUNT and the report printed
    `f"{n} of {n}"` -- the same value on both sides of "of", arithmetic that
    cannot come out any other way. It supported the sentence "every module of
    its simulation engine concerns ferroptosis or the physical-ROS
    modalities", which is FALSE: three modules contain zero matches of this
    file's own FERRO regex anywhere in their text.

    `lib.rs` is the crate root rather than a module and is excluded, so the
    denominator is the module count, not the file count.
    """
    d = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src"
    if not d.exists():
        return {"modules": 0, "mention": 0, "in_code": 0,
                "in_production_code": 0, "silent": []}
    mods = [f for f in sorted(d.glob("*.rs")) if f.stem != "lib"]
    mention, in_code, in_prod, silent = 0, 0, 0, []
    for f in mods:
        text = f.read_text(errors="ignore")
        if FERRO.search(text):
            mention += 1
        else:
            silent.append(f.name)
        # non-comment lines only: a module can cite ferroptosis in a doc
        # comment while its code is about geometry
        code = "\n".join(ln for ln in text.split("\n")
                          if not ln.strip().startswith(("//", "*", "/*")))
        if FERRO.search(code):
            in_code += 1
        # STRICTER STILL: drop the #[cfg(test)] block. Four modules pass the
        # in-code check only via test code -- a field name in a byte-identity
        # assert, a string literal in a CSV writer -- while the sentence
        # beside the row invites the reader to take it as production code.
        i = text.find("#[cfg(test)]")
        prod = text if i < 0 else text[:i]
        prod_code = "\n".join(ln for ln in prod.split("\n")
                               if not ln.strip().startswith(("//", "*", "/*")))
        if FERRO.search(prod_code):
            in_prod += 1
    return {"modules": len(mods), "mention": mention, "in_code": in_code,
            "in_production_code": in_prod, "silent": silent}


def _denominator_note(d: dict) -> str:
    """The mechanism-share denominator, derived.

    Both this page and README.md said the shares were "of TAGGED articles".
    They are shares of the CORPUS -- a different denominator, and the corpus
    carries three that are not interchangeable.
    """
    m = d.get("mechanism_denominators") or {}
    if not m:
        return ("* Mechanism shares quoted elsewhere in this repo are shares "
                "of the frozen corpus, not of the cancer literature; the "
                "mechanism taxonomy matches a small fraction of the census "
                "(see `analysis/atlas-taxonomy-reach.md`).")
    return (
        f"* Mechanism shares quoted elsewhere in this repo are shares of the "
        f"**{m['corpus_records']:,}-record frozen corpus** -- {m['top_mechanism']} is "
        f"{m['top_n']:,} of {m['corpus_records']:,} = "
        f"**{100*m['share_of_corpus']:.1f}%** -- not shares of the cancer "
        f"literature. An earlier version of this bullet, and of README.md, "
        f"called them shares of TAGGED articles; they are not, and the corpus "
        f"carries three denominators that are not interchangeable: "
        f"{100*m['share_of_corpus']:.1f}% of corpus articles, "
        f"{100*m['share_of_tagged']:.1f}% of the {m['tagged']:,} tagged "
        f"records, {100*m['share_of_tags']:.1f}% of all {m['tags']:,} "
        f"mechanism tags. The taxonomy itself matches about 6.8% of the "
        f"census through the production matcher -- see "
        f"`analysis/atlas-taxonomy-reach.md`.")


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
    em = d["engine_modules"]
    if isinstance(em, int):          # legacy artifact
        em = {"modules": em, "mention": em, "in_code": em, "silent": []}
    L += [f"| preregistered predictions | **{pf} of {len(p)}** |",
          "|---|--:|",
          f"| engine modules mentioning it anywhere | "
          f"**{em['mention']} of {em['modules']}** |",
          f"| engine modules mentioning it in code | "
          f"**{em['in_code']} of {em['modules']}** |",
          f"| engine modules mentioning it in PRODUCTION code | "
          f"**{em.get('in_production_code', em['in_code'])} of "
          f"{em['modules']}** |", ""]
    L += ["The module rows used to read `N of N` -- the same number on both "
          "sides of \"of\", arithmetic that cannot come out any other way, "
          "produced by counting `.rs` files without opening one. They are "
          "content measurements now, and the LAST is the one to read: a "
          "module can cite ferroptosis in a doc comment while its code is "
          "about geometry, and four modules pass the in-code check only via "
          "their `#[cfg(test)]` block -- a field name in a byte-identity "
          "assert, a string literal in a CSV writer.", ""]
    if em["silent"]:
        L += [f"**{len(em['silent'])} modules mention neither ferroptosis nor a "
              f"physical-ROS modality anywhere in their text**: "
              f"{', '.join(f'`{x}`' for x in em['silent'])}. An earlier "
              f"version of this page said every module of the engine concerns "
              f"one or the other. It does not.", ""]

    # DERIVED, because the sentence here was unconditional: a doctored input
    # with half the predictions non-ferroptosis still printed "Every
    # falsifiable commitment ...". This repo already fixed the same shape in
    # manuscript_vs_census.py, where a static headline could not flip.
    preds_all = (pf == len(p))
    # PRODUCTION CODE, not a mention. A module citing a PMID in a doc
    # comment does not "concern" ferroptosis in the sense the sentence
    # claims, and branching on `mention` meant three comment edits could
    # restore the universal claim this page retracts.
    mods_all = (em.get("in_production_code", em["mention"]) == em["modules"])
    if preds_all and mods_all:
        L += ["Every falsifiable commitment the project makes, and every "
              "module of its simulation engine, concerns ferroptosis or the "
              "physical-ROS modalities.", ""]
    elif preds_all:
        L += [f"Every falsifiable commitment the project makes concerns "
              f"ferroptosis or the physical-ROS modalities. The engine is "
              f"mostly but not entirely about them: {em['mention']} of "
              f"{em['modules']} modules mention one, and only "
              f"{em['in_code']} do so in code.", ""]
    else:
        L += [f"{pf} of {len(p)} preregistered predictions and "
              f"{em['mention']} of {em['modules']} engine modules concern "
              f"ferroptosis or the physical-ROS modalities.", ""]

    body_t = d.get("therapy_by_body") or []
    if body_t:
        overlap = {k: sorted(set(body_t) & set(a[k]))
                   for k in ("ferroptosis-or-physical", "therapy-subject",
                             "method")}
        L += [f"### The '{n_ther}' is a filename marker, not a subject "
              f"measurement", ""]
        L += [f"The two buckets do not use the same admission rule. Therapy is "
              f"matched on the FILENAME only; ferroptosis is matched on the "
              f"filename OR the first 40 lines of body text. So the therapy "
              f"count moves on a rename with contents unchanged, and an empty "
              f"file with a therapy word in its name is filed as a therapy "
              f"analysis.", ""]
        L += [f"Applying the ferroptosis bucket's body rule to the SAME "
              f"vocabulary admits **{len(body_t)}** analyses. That is not an "
              f"upper bound on 'other therapy', and an earlier version of "
              f"this page published it as one. Where those "
              f"{len(body_t)} already sit in the table above:", ""]
        L += ["| already classified as | body-route matches |", "|---|--:|"]
        for k, v in overlap.items():
            L.append(f"| {k} | {len(v)} |")
        L += [""]
        n_f = len(overlap["ferroptosis-or-physical"])
        n_m = len(overlap["method"])
        L += [f"**{n_f} of them are in this page's own FERROPTOSIS column** "
              f"({', '.join(f'`{x}`' for x in overlap['ferroptosis-or-physical'])}), "
              f"so under a mutually exclusive bucketing "
              f"{len(body_t)} cannot bound the therapy count. A further "
              f"{n_m} are method analyses -- instrument work that cites a "
              f"therapy as a worked example rather than studying one.", ""]
        L += [f"So the body route does not measure 'other therapy' either. "
              f"What it measures is how many analyses MENTION a therapy, "
              f"which is a different question, and the honest reading of the "
              f"table above is that the therapy figure is small and its exact "
              f"value is not established by either rule. An earlier version "
              # A QUOTATION FROM A PAST COMMIT, so it is a literal. Deriving
              # it from the CURRENT count rewrote history. An earlier note
              # here also got the history wrong twice over: the first
              # published range really did say 37, and the count moved
              # because THIS project's own census-findings page was
              # rewritten, not an unrelated one. Measured across commits the
              # body-route figure has read 37, then 38, then 37.
              f"published 'the true figure lies between 1 and 37', a range "
              f"whose upper end is not a bound. That figure is quoted as a "
              f"fact about a past commit and is deliberately NOT derived -- "
              f"it has read 37, 38 and 37 across commits, moving TWICE on "
              f"this script's 40-line head window and for opposite reasons -- "
              f"once when a sibling page's head grew INTO the window past the "
              f"three-hit threshold, once when this page's own therapy "
              f"mention was pushed OUT of it -- rather than when any "
              f"subject changed. The live "
              f"figures are {n_ther} and {len(body_t)}.", ""]

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
    L += ["* It does not say the focus is wrong. It says it is UNDECLARED, "
          "which is a different criticism and the only one this page makes.",
          "* It does not measure quality or effort, only subject. A method "
          "analysis is not lesser work; it is work about the instrument.",
          _denominator_note(d), ""]
    return "\n".join(L) + "\n"


def main():
    a, body_therapy = classify_analyses()
    p = classify_predictions()
    d = {"analyses": a, "therapy_by_body": body_therapy,
         "predictions": p, "engine_modules": engine_modules(),
         "mechanism_denominators": mechanism_denominators(),
         "n_analyses": sum(len(v) for v in a.values()),
         "n_predictions": len(p),
         "n_predictions_ferroptosis": sum(1 for v in p.values() if v),
         "n_therapy_subject": len(a["therapy-subject"])}
    if not d["engine_modules"]["modules"]:
        raise SystemExit("no engine modules found; the count would be a lie")
    OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                        encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  analyses: {len(a['ferroptosis-or-physical'])} ferroptosis, "
          f"{len(a['therapy-subject'])} other therapy, {len(a['method'])} method")
    print(f"  predictions: {d['n_predictions_ferroptosis']} of {d['n_predictions']} ferroptosis")
    print(f"  engine modules: {d['engine_modules']['mention']} of "
          f"{d['engine_modules']['modules']} mention; "
          f"{d['engine_modules']['in_code']} in code")


if __name__ == "__main__":
    main()
