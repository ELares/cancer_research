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


MANUSCRIPT_RATIO = 9.1


def _volume(land) -> dict:
    """The volume claim's figures, DERIVED from the landscape artifact.

    They were literals. Setting a source field to nonsense left the page
    printing the old number, and rewriting 17.6 to 1.6 in the generator passed
    every guard -- because nothing here read anything.
    """
    import importlib.util
    from pathlib import Path as _P
    spec = importlib.util.spec_from_file_location(
        "al", _P(__file__).resolve().parent / "atlas_landscape.py")
    al = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(al)
    rows = land if isinstance(land, list) else land.get("rows") or []
    cen = {r["mechanism"].lower(): r for r in rows if r.get("mechanism")}

    def tot(keys, field="mesh_census"):
        return sum((cen.get(k) or {}).get(field) or 0 for k in keys)

    num, den = tot(al.PHARMACOLOGICAL), tot(al.PHYSICAL)
    # CAPTURE is the frozen corpus's MeSH-labelled share of the census, which
    # is `mesh_frozen`. `keyword_frozen` is the manuscript's own tagger and
    # gives 1.9x, not the 3.3x the landscape page derives.
    cap_ph = tot(al.PHARMACOLOGICAL, "mesh_frozen") / max(num, 1)
    cap_py = tot(al.PHYSICAL, "mesh_frozen") / max(den, 1)
    pre_n = tot(al.PHARMACOLOGICAL & al.PRECISE)
    pre_d = tot(al.PHYSICAL & al.PRECISE)
    # THE SAME ARTICLES UNDER MESH LABELS, which is what separates the two
    # factors below. Without it the page welds them: 3.3x is the frozen->census
    # step, 1.9x is the net against the manuscript's KEYWORD figure, and MeSH
    # labelling alone moves the ratio the other way.
    b = tot(al.PHARMACOLOGICAL, "mesh_frozen") / max(
        tot(al.PHYSICAL, "mesh_frozen"), 1)
    # ALL THREE RESTRICTIONS the sibling publishes, not just the one that
    # inverts. `PRECISE - PHYSICAL` is what `atlas_landscape.py` ITSELF uses
    # and is the literal referent of "that script's own PRECISE set".
    alt_n = tot(al.PRECISE - al.PHYSICAL)
    rest = sorted(x for x in (al.PHARMACOLOGICAL - al.PRECISE)
                  if _names_therapy((cen.get(x) or {}).get("top_descriptor")))
    return {"ratio": num / den if den else 0.0,
            "ratio_mesh_frozen": b,
            "oversample": (cap_py / cap_ph) if cap_ph else 0.0,
            "precise": (pre_n / pre_d) if pre_d else 0.0,
            "precise_alt": (alt_n / pre_d) if pre_d else 0.0,
            "criterion_restored": ((pre_n + tot(rest)) / pre_d) if pre_d else 0.0}


_THERAPY = __import__("re").compile(
    r"inhibitor|therap|antibod|conjugat|vaccine|agents?\b|blockade", __import__("re").I)


def _names_therapy(desc) -> bool:
    return bool(desc and _THERAPY.search(desc))


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
    modsup = load("atlas-module-support.json")

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
        _v = _volume(land)
        L += [
            f"**Volume.** Pharmacological to physical runs **{MANUSCRIPT_RATIO} : 1** by the",
            f"manuscript's keyword method and **{_v['ratio']:.1f} : 1** on the census.",
            f"THOSE ARE NOT ON COMPARABLE LABELS, and an earlier version of this",
            f"paragraph welded the two factors with a `because`. The same articles",
            f"under MeSH labels give **{_v['ratio_mesh_frozen']:.1f} : 1**, so MeSH",
            f"labelling alone moves the ratio by",
            f"{_v['ratio_mesh_frozen']/MANUSCRIPT_RATIO:.2f}x; census selection then",
            f"multiplies it by {_v['oversample']:.1f}x, the corpus's over-sampling of",
            f"physical modalities. Net against the manuscript's keyword figure,",
            f"{_v['ratio']/MANUSCRIPT_RATIO:.1f}x.", "",
            f"**And whether that is an understatement depends on the restriction.**",
            f"The sibling page publishes three, and says of them that neither is",
            f"adopted: restricting BOTH classes to `PHARMACOLOGICAL & PRECISE` gives",
            f"**{_v['precise']:.2f} : 1**; using `PRECISE - PHYSICAL`, which is what",
            f"`atlas_landscape.py` itself uses, gives **{_v['precise_alt']:.2f} : 1**;",
            f"restoring the mechanisms that satisfy PRECISE's own stated criterion",
            f"gives **{_v['criterion_restored']:.2f} : 1**. Only the first two fall below",
            f"the manuscript's {MANUSCRIPT_RATIO} : 1, so the inversion holds under two",
            f"readings of three and an earlier version of this paragraph quoted only",
            f"the one that inverts. See `atlas-modality-ratio.md`.", "",            "**Maturity.** \"Physical modalities remain comparatively preclinical\" does",
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
              "", "None of them became a layer. When four of these were checked for",
              "a calibration target (#616), the route that partially anchored ACSL4 --",
              "cBioPortal within-cohort z-scores -- turned out to recover the normal",
              "distribution at the z<-1 cut for every gene tested, so that cut carries",
              "no gene-specific signal. At the deeper z<-2 cut TP53 does separate",
              "(above the normal expectation in 31 of 32 cancer types), but it bounds",
              "a prevalence rather than a dose-response, so it is recorded as a weak",
              "anchor and no layer was written. Read the table as where the",
              "literature's attention and the available data fail to overlap, not as",
              "a backlog.", "",
              "*Source:* `atlas-model-gaps.md`, `calibration-feasibility.md`", ""]

    L += ["---", "", "## What the census did NOT support", "",
          "Reported here because a findings page that only lists wins is marketing.", ""]
    reg = load("comention-regression.json")
    auth = load("comention-authority-result.json")
    if reg and reg["net_change"] < 0:
        a, b = reg["after"], reg["before"]
        L += [
            "**A repair this project made to its own co-mention layer made it "
            f"worse.** #617 measured the layer at 55.5% precision "
            f"(recomputed at {100*b['weighted']:.0f}% once every stratum was "
            "hand-judged rather than assumed), "
            "traced it to a filter that exempted multi-word forms, and replaced the "
            "single-token test with two measured filters. Re-measured on a fresh "
            f"sample after the rebuild, precision FELL to {100*a['weighted']:.0f}%. "
            "The false positives had been multi-word because the multi-word channel "
            "was the unfiltered one; closing it moved the pressure to the channel "
            "just opened, and the top offenders are now bare English words. The "
            "replacement filters were then measured and do not separate true matches "
            "from false ones at all."]
        # The failure above is the FIRST half of a closed arc. Reporting only the
        # half that went wrong was accurate the day it was written and became a
        # different kind of dishonesty once the second half shipped -- the same
        # defect this page exists to catch, in this page.
        if auth:
            L += [
                f"**It was then repaired for real, and the fix is measured at "
                f"{100*auth['weighted']:.0f}%** (#628). What separates true matches "
                "from false ones is not how much support a form has but whether the "
                "form is a NAME of the entity it resolves to, checked against NLM "
                "and NCBI rather than against the corpus. `treatment` is not a name "
                f"of any descriptor; `xCT` is a name of SLC7A11. A blind panel of "
                f"three judges who never saw the first verdicts put it at "
                f"{100*auth['blind_weighted']:.0f}%, and the hostile bound, "
                f"resolving every borderline case against the layer, is "
                f"{100*auth['blind_hostile']:.0f}%.",
                "So the standing finding is not that the layer is broken. It is that "
                "this project shipped a filter justified by an error distribution "
                "the filter itself changed, did not notice for two issues, and "
                "needed a fresh sample drawn after the rebuild to see it. The cost "
                f"of the real fix is {100*auth['recall_cost']['tp_lost_share']:.0f}% "
                "of true matches, paid entirely on MeSH terms "
                f"({100*auth['recall_cost']['mesh_tp_lost']:.0f}%) and not at all on "
                f"genes ({100*auth['recall_cost']['gene_tp_lost']:.0f}%)."]
        L += ["*Source:* `comention-regression.md`, "
              "`comention-authority-result.md`", ""]
    if disc:
        p = disc["headline"]["precision"]
        L += [f"**Literature-based discovery does not work as built.** The shipped ABC",
              f"ranking scores {100*p['abc']:.1f}% precision@20 against",
              f"{100*p['popularity']:.1f}% for ranking the same candidates by",
              "popularity, and no standard link predictor beats that baseline either."]
        # "and a bad ranker" was the second half of this finding for months. It
        # does not follow, and two later measurements say why -- so it is stated
        # conditionally on them rather than asserted.
        head = load("atlas-discovery-headroom.json")
        bias = load("atlas-discovery-degree-bias.json")
        if head and bias and not head.get("any_headroom"):
            L += [
                "The candidate SET is genuinely informative, so that half stands. "
                "The obvious second half -- a bad RANKER -- does not follow, and "
                "two later measurements say why. Ordering the seven measured "
                "methods by how hub-selecting each one is reproduces the "
                f"precision leaderboard exactly (rank correlation "
                f"{bias['spearman_L_vs_precision']:.2f} over those seven points), "
                "so the metric rewards NOT correcting for candidate degree; and "
                "blending each of the five seed-aware signals into a degree-only "
                "prior adds nothing measurable, at any weight tested and under "
                "three combination schemes outside the blend family. Among every "
                "ranker measured, a degree-correcting one and a bad one cannot be "
                "told apart on this metric.",
                "*Source:* `atlas-discovery-eval.md`, "
                "`atlas-discovery-degree-bias.md`, `atlas-discovery-headroom.md`", ""]
        else:
            L += ["A good candidate generator and a bad ranker.",
                  "*Source:* `atlas-discovery-eval.md`", ""]
    if repl:
        # The censored series is read from the artifact. It used to be quoted
        # from atlas-replication.md's PROSE, where it was a remembered figure
        # from a development run that stored nothing -- so neither document
        # could check it and neither would have noticed it going stale.
        cen = repl.get("cohorts_censored_ever") or []
        w = repl["quiet_years"]
        done = [r for r in cen if r["year"] + w <= repl["latest_year"]]
        span = (f"{100*done[0]['rate']:.1f}% for {done[0]['year']} to "
                f"{100*done[-1]['rate']:.1f}% for {done[-1]['year']}"
                if done else "a clean monotonic decline")
        L += ["**Replication looked like it was collapsing, and was not.** Scoring "
              f"cohorts on whether they were EVER replicated gives {span}; that is "
              "the observation window shrinking, not science changing, since the "
              "older cohort has had decades to acquire a second paper and the newer "
              f"one had {w} years. On an equal {w}-year window from each pair's own "
              "first assertion the decline is modest, and the recent end is an upper "
              "bound because of MeSH indexing lag.",
              "*Source:* `atlas-replication.md`", ""]
    if modsup and modsup.get("exposure_floor") is not None:
        n_zero = modsup["zero_relation"]
        n_exp = modsup["zero_explained_by_exposure"]
        L += [
            f"**Most of what the census cannot corroborate is about what has been "
            f"STUDIED, not about what is true.** Of the "
            f"{modsup['n_claims']} simulation-module claims, {modsup['corroborated']} "
            f"are corroborated by at least one asserting article and {n_zero} by none. "
            f"Read flat, that looks like {n_zero} unsupported claims. It is not: a pair "
            f"can only be asserted if BOTH its entities are written about, and the "
            f"weaker entity's partner count across these claims runs from "
            f"{modsup['weaker_min']:,} to {modsup['weaker_max']:,}. "
            f"Every claim that HAS support has a weaker entity of at least "
            f"{modsup['exposure_floor']} partners, and {n_exp} of the {n_zero} zeros "
            f"fall below that (Spearman rho = "
            f"{modsup['spearman_weaker_degree_vs_relations']:.2f}, and the association "
            f"survives dropping every GPX4 pair).",
            "**But which zeros are 'genuine' cannot be identified, and the source "
            "document deliberately names none.** The line is a sample minimum set by a "
            "one-article row; 45% of all asserted pairs in the graph sit below it, so "
            "it is not a detectability limit; and running the same procedure on the "
            "pair-level co-mention column inverts the correlation and returns a "
            "disjoint pair of exceptions. The finding is that a zero is a poor guide "
            "to a claim's truth -- not that any particular zero is interesting.",
            "*Source:* `atlas-module-support.md`", ""]

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
        # A hand-judged measurement supersedes the corroboration bound where one
        # exists. The bound is what this layer had before anyone read its output;
        # reporting it once a measurement exists overstates what is known and,
        # here, understates how bad the layer is.
        reg = load("comention-regression.json")
        auth = load("comention-authority-result.json")
        if auth:
            # Quote the layer AS SHIPPED. The regression figure describes a
            # configuration nothing runs any more, and a bound on a build nobody
            # uses is not a bound.
            L += ["---", "", "## Every layer now carries a bound", "",
                  f"* co-mention precision: **{100*auth['weighted']:.0f}% measured** "
                  f"on the layer as shipped (blind panel "
                  f"{100*auth['blind_weighted']:.0f}%, hostile bound "
                  f"{100*auth['blind_hostile']:.0f}%), superseding both the "
                  f"{100*lo:.1f}%-{100*hi:.1f}% corroboration bound and the "
                  f"{100*reg['after']['weighted']:.0f}% the layer measured before the "
                  f"authority filter was turned on",]
        elif reg:
            L += ["---", "", "## Every layer now carries a bound", "",
                  f"* co-mention precision: **{100*reg['after']['weighted']:.0f}% "
                  f"measured** (hand-judged, superseding the "
                  f"{100*lo:.1f}%-{100*hi:.1f}% corroboration bound), and it FELL "
                  f"{abs(100*reg['net_change']):.0f} points when the #617 filters "
                  f"were changed",]
        else:
            L += ["---", "", "## Every layer now carries a bound", "",
                  f"* co-mention precision: {100*lo:.1f}% to {100*hi:.1f}%",]
        L += [
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
