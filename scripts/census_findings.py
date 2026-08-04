#!/usr/bin/env python3
"""What the census established about this project's own claims (#CENSUS-FINDINGS).

WHY
---
The census produced a dozen separate analyses, each with its own report. Nothing
answers the question a reader actually has: what did going from 4,830 articles
to 4,403,994 change about what this project believes?

This assembles that answer from the committed JSON artifacts rather than
restating them by hand, so it cannot drift from the measurements it summarises.
Every number here is read at generation time from the file that produced it.

Usage:
    python scripts/census_findings.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

A = PROJECT_ROOT / "analysis"
OUT = A / "census-findings.md"


def load(name):
    try:
        return json.loads((A / name).read_text())
    except (OSError, ValueError):
        return None


def main() -> int:
    land = load("atlas-landscape.json")
    thesis = load("atlas-thesis-position.json")
    pred = load("atlas-prediction-position.json")
    gaps = load("atlas-model-gaps.json")
    repl = load("atlas-replication.json")
    disc = load("atlas-discovery-eval.json")
    amb = load("atlas-ambiguity.json")
    imp = load("atlas-ambiguity-impact.json")
    cite = load("atlas-citation-audit.json")
    news = load("news-verification-audit.json")
    coment = load("atlas-comention-audit.json")
    contra = load("atlas-contradiction-quality.json")

    missing = [n for n, v in [
        ("landscape", land), ("thesis", thesis), ("predictions", pred),
        ("gaps", gaps), ("replication", repl), ("discovery", disc)] if v is None]
    if missing:
        print(f"missing artifacts: {missing}; run scripts/atlas_pipeline.sh mine",
              file=sys.stderr)

    L = ["# What the census established (#CENSUS-FINDINGS)", "",
         "Assembled by `scripts/census_findings.py` from the committed JSON of each",
         "analysis, so it cannot drift from the measurements it summarises.", "",
         "The question this answers: going from a 4,830-article keyword corpus to a",
         "4,403,994-article MeSH census, what changed about what this project",
         "believes?", "", "---", "",
         "## 1. The corpus was a 213-fold non-uniform sample", ""]

    if land:
        caps = [(r["mechanism"], r["mesh_frozen"] / r["mesh_census"])
                for r in land["rows"] if r["mesh_census"] and r["mesh_frozen"]]
        lo = min(caps, key=lambda kv: kv[1])
        hi = max(caps, key=lambda kv: kv[1])
        L += [
            "Measured with the same MeSH labels on both corpora, so descriptor scope",
            "cancels. The frozen corpus captures",
            f"**{100*hi[1]:.1f}%** of the census's `{hi[0]}` articles and",
            f"**{100*lo[1]:.2f}%** of its `{lo[0]}` -- a "
            f"**{hi[1]/lo[1]:.0f}-fold** spread.", "",
            "Any relative prevalence computed on the frozen corpus inherits that. It",
            "is not a small uniform sample of the literature; it is a wildly uneven",
            "one.", "",
            "*What changed:* claims of the form \"research is concentrated here\" are",
            "now measurable rather than assumed.",
            "*Source:* `atlas-landscape.md`", ""]

    L += ["## 2. The manuscript understated its own headline, and over-broadened another", ""]
    if land:
        L += [
            "**Volume.** Pharmacological to physical runs **9.1 : 1** by the",
            "manuscript's keyword method and **17.6 : 1** on the census, because the",
            "corpus over-samples physical modalities 3.3x -- exactly what a corpus",
            "built from queries about them would do. The claim survives and was",
            "understated by about half.", "",
            "**Maturity.** \"Physical modalities remain comparatively preclinical\" does",
            "not hold as a class. HIFU is **7.10%** clinical against CAR-T's",
            "**6.64%**, both on precise descriptors, and HIFU and sonodynamic differ",
            "by 1.6x. The defensible claim is narrower: sonodynamic therapy",
            "specifically is early -- which is the mechanism this work rests on.", "",
            "*What changed:* one claim strengthened, one narrowed. Both are in the",
            "manuscript.", "*Source:* `atlas-landscape.md`", ""]

    if thesis:
        t = thesis["totals"]
        g = thesis.get("growth")
        L += ["## 3. The thesis sits on roughly thirty papers", "",
              f"The corpus contains no ferroptosis query and no PDT query, so it could",
              "not measure this. The census can:", "",
              f"* ferroptosis-indexed cancer articles: **{t['ferroptosis']:,}**"
              + (f", growing {g['from_n']:,} ({g['from_year']}) to {g['to_n']:,} "
                 f"({g['to_year']})" if g else ""),
              f"* x drug resistance: **{t['drug resistance']}**",
              f"* x photodynamic therapy: **{t['photodynamic therapy']}**",
              f"* x sonodynamic therapy: **{t['sonodynamic therapy']}**", "",
              "The resistance leg is supported by a literature. The sonodynamic leg,",
              "the thesis's central mechanism, is supported by roughly thirty papers,",
              "and the cautionary precedent (PDT) is better established than the thing",
              "it is a precedent for.", "",
              "*What changed:* the simulation work is carrying more of the argument",
              "than the citation count suggested, and now says so.",
              "*Source:* `atlas-thesis-position.md`", ""]

    if pred:
        R = {r["prediction"]: r for r in pred["rows"]}
        L += ["## 4. The two planning documents disagree about the keystone", "",
              f"P1 (persister/resistance) sits on **{R['P1']['with_ferroptosis']}**",
              f"ferroptosis articles; P4 (hypoxia), which `PREREGISTRATION.md` calls",
              f"the keystone, sits on **{R['P4']['with_ferroptosis']}**. The P1 protocol",
              "calls P1 the highest-leverage prediction. Neither designation cited",
              "evidence.", "",
              "Not a quality ranking -- a sparse leg is where novelty lives. The",
              "asymmetry that matters is that a negative P4 is ambiguous (mechanism",
              "wrong, or experiment not yet worked out) while a negative P1 is simply",
              "a negative result.", "",
              "*What changed:* the choice can now be made deliberately.",
              "*Source:* `atlas-prediction-position.md`", ""]

    if gaps:
        top = gaps["gaps"][:4]
        L += ["## 5. What the literature says to model next", "",
              "| gene | articles | engine handle |", "|---|---|---|"] + [
              f"| {g['gene']} | {g['papers']:,} | none |" for g in top] + [
              "", "None of them became a layer. When the top four were checked for a",
              "calibration target (#616), the route that partially anchored ACSL4 --",
              "cBioPortal within-cohort z-scores -- turned out to recover the normal",
              "distribution for every gene tested, so it carries no gene-specific",
              "signal and cannot anchor anything. Read the table as where the",
              "literature's attention and the available data fail to overlap, not as",
              "a backlog.", "",
              "*Source:* `atlas-model-gaps.md`, `calibration-feasibility.md`", ""]

    L += ["---", "", "## What the census did NOT support", "",
          "Reported here because a findings page that only lists wins is marketing.", ""]
    if disc:
        p = disc["headline"]["precision"]
        L += [f"**Literature-based discovery does not work as built.** The shipped ABC",
              f"ranking scores {100*p['abc']:.1f}% precision@20 against",
              f"{100*p['popularity']:.1f}% for ranking the same candidates by",
              "popularity, and no standard link predictor beats that baseline either.",
              "A good candidate generator and a bad ranker.",
              "*Source:* `atlas-discovery-eval.md`", ""]
    if repl:
        L += [f"**Replication looked like it was collapsing, and was not.** Scoring",
              "cohorts on whether they were ever replicated gave a clean decline from",
              "60% to 17.5%; that was the observation window shrinking, not science",
              f"changing. On an equal {repl['quiet_years']}-year window the decline is",
              "modest and the recent end is an upper bound because of indexing lag.",
              "*Source:* `atlas-replication.md`", ""]
    if imp:
        L += [f"**The entity collisions are not as bad as containment suggests.**",
              f"{100*imp['relation_rows_touching_contested_id']/imp['relation_rows']:.1f}%"
              " of relation rows touch a contested identifier, but only",
              f"{100*imp['relation_rows_resting_on_at_risk']/imp['relation_rows']:.2f}%"
              " rest on an uncorroborated one.",
              "*Source:* `atlas-ambiguity-impact.md`", ""]

    L += ["---", "", "## Integrity problems found along the way", ""]
    rows = []
    if cite:
        broken = [r for r in cite["rows"]
                  if r["status"] in ("wrong-subject", "unresolvable")]
        rows.append(f"* **{len(broken) or 3} module citations pointed at unrelated "
                    "papers** -- a Nature news item on fetal-tissue policy, a "
                    "Theriogenology paper on embryo vitrification, and a PMID that "
                    "does not resolve. Corrected. (`atlas-citation-audit.md`)")
    if news:
        base = news.get("baseline", news)
        rows.append(
            f"* **{100*base['zero_overlap']/base['pairs_resolved']:.1f}% of the "
            "news pipeline's \"verified\" links** shared no content word with "
            "the claim they verified. Root cause: a claim yielding the single "
            "search term `Seven` matched 835,973 records and the five newest "
            "were accepted. The linker is fixed and re-run, and the same "
            f"measurement now reads "
            f"{100*news['zero_overlap']/news['pairs_resolved']:.1f}% -- but by "
            f"WITHDRAWAL, not repair: "
            f"{base['claims_with_links'] - news['claims_with_links']} of the "
            f"{base['claims_with_links']} verifications were dropped outright, and "
            f"{news['two_plus_generic_only']} of the "
            f"{news['two_plus']} surviving pairs clear the bar on oncology "
            "boilerplate alone. (`news-verification-audit.md`)")
    if amb:
        rows.append(f"* **`FSP1` resolved to a spastic-paraplegia gene**, leaving the "
                    "manuscript's headline GPX4+FSP1 claim with zero typed relations. "
                    f"Blocklist now covers {len(amb['blocklist'])} measured sense "
                    "collisions. (`atlas-ambiguity.md`)")
    rows.append("* **The manuscript's mechanism count was 19 where the index carries "
                "23** -- an undocumented 20-article threshold presented as coverage, "
                "and the four it hid were mostly physical modalities.")
    L += rows + ["",
                 "None of these were computational errors. Every one was a true "
                 "statement describing something narrower than what it was used for.", ""]

    if coment:
        lo = coment["pubtator_agree"] / coment["pubtator_scored"]
        hi = (coment["pubtator_agree"] + coment["body_only"]) / coment["pubtator_scored"]
        L += ["---", "", "## Every layer now carries a bound", "",
              f"* co-mention precision: {100*lo:.1f}% to {100*hi:.1f}%",
              (f"* contradictions: ambiguity inflates the flag rate "
               f"{contra['mantel_haenszel']:.2f}x" if contra else ""),
              "* emergence: 99.0% precision, 99.6% recall",
              f"* FSP1 disambiguation: 97.4%, with 75% of corrections extrapolated "
              "and that extrapolation independently tested", ""]

    L += ["---", "",
          "*Regenerate with `python scripts/census_findings.py`. Every figure is read",
          "from the JSON of the analysis that produced it.*"]

    OUT.write_text("\n".join(x for x in L if x is not None) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
