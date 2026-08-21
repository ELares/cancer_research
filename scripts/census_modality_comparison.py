#!/usr/bin/env python3
"""Section 9.4's shared-axis comparison, for the axes the census can measure.

The manuscript proposes evaluating the physical-ROS hypothesis against its
alternatives on four shared axes -- regulatory maturity, evidence depth,
delivery burden, and applicability to residual disease -- and names radioligand
therapy as the most informative comparison, on the argument that if it meets
the same clinical need through a more mature path then the physical-ROS
hypothesis is academically interesting and clinically unnecessary.

TWO OF THE FOUR AXES ARE NOW MEASURABLE and two are not, so this does half the
job and says which half. Evidence depth and regulatory maturity both reduce to
what NLM records: whether an article carries a clinical-trial publication type,
and which phase it is assigned. Delivery burden and applicability to residual
disease do not reduce to anything in a bibliographic record, and no arrangement
of counts will produce them.

THE INSTRUMENT DIFFERS BY MODALITY AND THAT DECIDES THE DESIGN. Radioligand
therapy -- the comparison the section turns on -- has no MeSH descriptor, so it
is measurable only through this project's keyword vocabulary, while
immunotherapy and CAR-T have precise descriptors. Comparing a text-arm number
against a descriptor-arm number is an instrument mismatch, and it would fall
exactly where the argument is: the alternative would be measured one way and
the thesis modality another.

So every modality is measured on BOTH arms wherever both exist, each row states
which arms it has, and the comparison the section asks for is drawn on the arm
that covers all of them -- with the descriptor arm shown beside it as a check
rather than as the answer.

WHAT A TRIAL SHARE IS NOT: evidence that a modality WORKS. It counts trials
that happened, not trials that succeeded, and a modality with a low barrier to
a first-in-human study can accumulate them without any of them reading out
well. It measures how far along a development path a literature is, which is
the axis Section 9.4 asks for, and nothing more.
"""
import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-modality-comparison.md"
OUT_JSON = REPO / "analysis/census-modality-comparison.json"

TRIAL_TYPES = {
    "Clinical Trial", "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial",
}
PHASES = ["Clinical Trial, Phase I", "Clinical Trial, Phase II",
          "Clinical Trial, Phase III", "Clinical Trial, Phase IV"]
# The modalities Section 9.4 puts on the same axes: the thesis modality, the
# alternative it names, and the physical and pharmacological comparators it
# invokes. Named here rather than swept, because the section names them.
COMPARED = ["sonodynamic", "hifu", "electrochemical-therapy",
            "radioligand-therapy", "ttfields", "immunotherapy", "car-t",
            "antibody-drug-conjugate"]
# EVERY mechanism carrying both instruments, for the validity sweep. The
# sonodynamic artifact was found by accident on a subset chosen for a different
# reason, which is not a way to learn whether other descriptors are borrowed
# too. Filled at scan time from the intersection of the two vocabularies.
SWEEP_ALL = True
# A share off a handful of trials is a statement about the handful.
MIN_TRIALS = 10
# THE FLOOR THAT MATTERS FOR THE VALIDITY TEST IS ON ARTICLES, NOT TRIALS, and
# a first version got this exactly backwards. Requiring both arms to clear ten
# TRIALS excluded sonodynamic -- the case the whole test was built from --
# because a descriptor that inflates a share does so by piling trials onto the
# descriptor arm while the term's own arm stays thin. So the worst artifacts
# are precisely the ones a trial floor hides.
#
# The share ratio is trial-ratio over article-ratio, so what it needs is enough
# ARTICLES on both arms to estimate two shares, plus enough trials on at least
# one arm for the numerator to mean something. Sonodynamic has 2,513 and 1,379
# articles: amply testable, and its 114-against-4 trials is the finding rather
# than a reason to skip it.
MIN_ARTICLES_FOR_VALIDITY = 200


def _mesh_sets():
    import yaml

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    return {k.lower(): {x.lower() for x in v["descriptors"]}
            for k, v in mp.items() if v["descriptors"]}


def _text_patterns():
    import sys

    sys.path.insert(0, str(REPO / "scripts"))
    import config

    # WORD-BOUNDED, for the reason the translation-lag analysis established:
    # `car t` matches inside "scar tissue". A share is more forgiving of that
    # than a first-appearance year, but the fix costs nothing.
    return {k.lower(): re.compile(
        "|".join(rf"\b{re.escape(t.lower())}\b" for t in terms))
        for k, terms in config.MECHANISM_KEYWORDS.items() if terms}


def scan(stride: int = 1) -> dict:
    mesh = _mesh_sets()
    text = _text_patterns()
    # The sweep needs every mechanism with both arms; the Section 9.4 view is a
    # named subset of the same scan rather than a second pass.
    want = set(COMPARED) | (set(mesh) & set(text) if SWEEP_ALL else set())
    mesh = {k: v for k, v in mesh.items() if k in want}
    text = {k: v for k, v in text.items() if k in want}
    prefilter = re.compile("|".join(p.pattern for p in text.values()))

    arms = {"mesh": defaultdict(Counter), "text": defaultdict(Counter)}
    n = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                pts = set(r.get("pub_types") or [])
                is_trial = bool(pts & TRIAL_TYPES)
                phases = [p for p in PHASES if p in pts]

                ms = {m.lower() for m in (r.get("mesh") or [])}
                hits = {"mesh": [k for k, d in mesh.items() if ms & d] if ms else []}
                blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
                hits["text"] = ([k for k, p in text.items() if p.search(blob)]
                                if blob and prefilter.search(blob) else [])
                for arm, ks in hits.items():
                    for k in ks:
                        arms[arm][k]["articles"] += 1
                        if is_trial:
                            arms[arm][k]["trials"] += 1
                        for p in phases:
                            arms[arm][k][p] += 1
    return {
        "census": n,
        "min_trials": MIN_TRIALS,
        "compared": COMPARED,
        "swept": sorted(set(mesh) & set(text)),
        "mesh_measurable": sorted(mesh),
        "text_measurable": sorted(text),
        "arms": {a: {k: dict(v) for k, v in arms[a].items()} for a in arms},
    }


def assemble(d: dict) -> dict:
    rows = []
    for k in sorted(set(d["compared"]) | set(d.get("swept", []))):
        row = {"modality": k,
               "mesh_measurable": k in d["mesh_measurable"],
               "text_measurable": k in d["text_measurable"]}
        for arm in ("mesh", "text"):
            c = d["arms"][arm].get(k, {})
            art, tri = c.get("articles", 0), c.get("trials", 0)
            row[f"{arm}_articles"] = art
            row[f"{arm}_trials"] = tri
            row[f"{arm}_trial_share"] = round(100 * tri / art, 2) if art else None
            phased = {p.replace("Clinical Trial, Phase ", ""): c.get(p, 0)
                      for p in PHASES}
            row[f"{arm}_phases"] = phased
            n_ph = sum(phased.values())
            row[f"{arm}_phased"] = n_ph
            # LATE-PHASE SHARE is the regulatory-maturity axis: a literature
            # with trials but none past Phase II is at a different point on a
            # development path from one with Phase III work, and the raw trial
            # count cannot tell them apart.
            row[f"{arm}_late_phase"] = phased["III"] + phased["IV"]
            row[f"{arm}_late_share"] = (
                round(100 * row[f"{arm}_late_phase"] / n_ph, 1) if n_ph else None)
            row[f"{arm}_interpretable"] = tri >= d["min_trials"]
        # ARM AGREEMENT AS A DESCRIPTOR-VALIDITY TEST. Where a descriptor names
        # what the words name, the two arms should land close: they read the
        # same literature through different instruments. A large divergence
        # says the descriptor is measuring something else, and it says so
        # QUANTITATIVELY rather than as the qualitative "this descriptor is
        # broad" caveat this project has carried in prose.
        a, b = row["text_trial_share"], row["mesh_trial_share"]
        row["arm_ratio"] = (round(max(a, b) / min(a, b), 1)
                            if a and b else None)
        row["arms_disagree"] = bool(row["arm_ratio"] and row["arm_ratio"] >= 3)
        rows.append(row)
    out = dict(d)
    out["rows"] = rows
    # The arm every compared modality has, which is the one the comparison must
    # be drawn on.
    out["common_arm"] = ("text" if all(r["text_measurable"] for r in rows)
                         else None)
    out["mesh_missing"] = sorted(r["modality"] for r in rows
                                 if not r["mesh_measurable"])
    out["arms_disagree"] = sorted(r["modality"] for r in rows
                                  if r["arms_disagree"])
    # THE SWEEP. Only rows with both arms and enough trials on each to make a
    # share mean anything; a mechanism with three trials on one arm produces a
    # large ratio for reasons that have nothing to do with its descriptor.
    testable = [r for r in rows
                if r["mesh_articles"] >= MIN_ARTICLES_FOR_VALIDITY
                and r["text_articles"] >= MIN_ARTICLES_FOR_VALIDITY
                and max(r["mesh_trials"], r["text_trials"]) >= MIN_TRIALS
                and r["arm_ratio"] is not None]
    out["validity_tested"] = sorted(r["modality"] for r in testable)
    out["validity_failed"] = sorted(r["modality"] for r in testable
                                    if r["arms_disagree"])
    out["section_94_view"] = list(COMPARED)
    out["min_articles_for_validity"] = MIN_ARTICLES_FOR_VALIDITY
    # Rows a trial-based floor WOULD have hidden, kept visible because the
    # motivating case was one of them.
    out["hidden_by_a_trial_floor"] = sorted(
        r["modality"] for r in testable
        if min(r["mesh_trials"], r["text_trials"]) < MIN_TRIALS)
    return out


def render(d: dict) -> str:
    arm = d["common_arm"] or "text"
    L = ["# Modalities on the axes a census can measure\n"]
    L.append(
        f"Generated by `scripts/census_modality_comparison.py` over "
        f"{d['census']:,} census records. Section 9.4 proposes evaluating the "
        f"physical-ROS hypothesis against its alternatives on four shared axes "
        f"-- regulatory maturity, evidence depth, delivery burden, and "
        f"applicability to residual disease -- and names radioligand therapy "
        f"as the comparison that matters most.\n"
    )
    L.append(
        "**Two of those axes are measurable here and two are not.** Evidence "
        "depth and regulatory maturity both reduce to what NLM records: "
        "whether an article carries a clinical-trial publication type, and "
        "which phase it is assigned. Delivery burden and applicability to "
        "residual disease reduce to nothing in a bibliographic record, and no "
        "arrangement of counts produces them. This is half of the comparison "
        "the section asks for.\n"
    )
    if d["mesh_missing"]:
        L.append(
            f"**The instrument differs by modality, which decides the design.** "
            + ", ".join(f"`{m}`" for m in d["mesh_missing"])
            + " have no MeSH descriptor, so they are measurable only through "
              "this project's keyword vocabulary, while immunotherapy and CAR-T "
              "have precise descriptors. Comparing a text-arm number against a "
              "descriptor-arm number would be an instrument mismatch falling "
              "exactly where the argument is -- the alternative measured one "
              f"way and the thesis modality another. The comparison is drawn "
              f"on the **{arm} arm**, which covers every modality, with the "
              f"descriptor arm beside it as a check rather than as the "
              f"answer.\n"
        )
    L.append(f"## Evidence depth and regulatory maturity ({arm} arm)\n")
    L.append("| modality | articles | trials | trial share | phased | "
             "Ph I | Ph II | Ph III | Ph IV | late-phase share |")
    L.append("|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for r in sorted(d["rows"], key=lambda r: -(r[f"{arm}_trial_share"] or 0)):
        ph = r[f"{arm}_phases"]
        thin = "" if r[f"{arm}_interpretable"] else " *"
        ls = (f"{r[f'{arm}_late_share']}%" if r[f"{arm}_late_share"] is not None
              else "-")
        L.append(
            f"| {r['modality']}{thin} | {r[f'{arm}_articles']:,} | "
            f"{r[f'{arm}_trials']:,} | {r[f'{arm}_trial_share']}% | "
            f"{r[f'{arm}_phased']:,} | {ph['I']:,} | {ph['II']:,} | "
            f"{ph['III']:,} | {ph['IV']:,} | {ls} |")
    L.append("")
    thin = [r["modality"] for r in d["rows"] if not r[f"{arm}_interpretable"]]
    if thin:
        L.append(f"\\* fewer than {d['min_trials']} trials, so the share "
                 f"describes a handful: " + ", ".join(f"`{m}`" for m in thin)
                 + ".\n")
    L.append("## Arm agreement is a descriptor-validity test\n")
    L.append(
        "Where a MeSH descriptor names what the words name, the two arms should "
        "land close -- they read the same literature through different "
        "instruments. A large divergence says the descriptor is measuring "
        "something else, and says it as a NUMBER rather than as the "
        "qualitative \"this descriptor is broad\" caveat this project has "
        "carried in prose.\n"
    )
    L.append("| modality | text share | descriptor share | ratio |")
    L.append("|---|--:|--:|--:|")
    for r in sorted(d["rows"], key=lambda r: -(r["arm_ratio"] or 0)):
        if r["mesh_trial_share"] is None:
            continue
        mark = "**" if r["arms_disagree"] else ""
        L.append(f"| {r['modality']} | {r['text_trial_share']}% | "
                 f"{r['mesh_trial_share']}% | {mark}{r['arm_ratio']}x{mark} |")
    L.append("")
    L.append(
        f"Swept across every mechanism carrying both instruments: "
        f"**{len(d['validity_failed'])} of {len(d['validity_tested'])} fail**, "
        f"which is what makes a failure interpretable rather than a property "
        f"of the method. Testability requires "
        f"{d['min_articles_for_validity']} articles on each arm and "
        f"{d['min_trials']} trials on at least one.\n"
    )
    if d["hidden_by_a_trial_floor"]:
        L.append(
            f"**The floor is on ARTICLES, not trials, and a first version had "
            f"it backwards.** Requiring ten trials on BOTH arms excluded "
            + ", ".join(f"`{m}`" for m in d["hidden_by_a_trial_floor"])
            + " -- including the case this test was built from. A descriptor "
              "that inflates a share does it by piling trials onto the "
              "descriptor arm while the term's own arm stays thin, so the "
              "worst artifacts are exactly the ones a trial floor hides. The "
              "test would have been blind to the thing that motivated it.\n"
        )
    DIAGNOSES = {
        "sonodynamic": (
            "**The descriptor is too BROAD.** Of the 114 census articles "
            "carrying `Ultrasonic Therapy` and a trial publication type, "
            "exactly **2** mention sonodynamic therapy or a sonosensitiser "
            "anywhere in title or abstract. The other 112 are ultrasound "
            "HYPERTHERMIA trials, mostly from the 1980s and 1990s. Neither "
            "reading favours the descriptor here: hyperthermia is not "
            "sonodynamic therapy under any definition, so 4.54% is 98% "
            "borrowed and 0.29% is the figure. THE MANUSCRIPT WAS CORRECTED."),
        "nanoparticle": (
            "**The descriptor is too NARROW, and the disagreement CONFIRMS a "
            "caveat the manuscript already makes.** 81% of the text arm's "
            "1,308 trials match only on `liposome`: approved liposomal "
            "medicines are indexed by drug name, not under `Nanoparticles`. "
            "So the two arms answer different questions -- 0.49% is how "
            "clinical the literature that calls itself nanoparticle research "
            "is, and 1.52% is how clinical nanoparticle-based delivery is once "
            "liposomal drugs are counted. The manuscript's claim is about "
            "novel experimental nanoplatforms, so 0.49% is the right figure "
            "FOR THAT CLAIM and no correction follows -- but the caveat is now "
            "a measurement rather than an assertion."),
        "mrna-vaccine": (
            "**The descriptor is too YOUNG.** `mRNA Vaccines` was minted in "
            "2020 and its earliest census record is 2020, so it cannot see a "
            "trial that predates it however relevant. The text arm reaches "
            "earlier neoantigen-vaccine work. This is the descriptor-"
            "introduction effect the translation-lag analysis measures, "
            "showing up in a share rather than in a first-appearance year."),
    }
    if d["arms_disagree"]:
        L.append("### What each failure turned out to be\n")
        L.append(
            "Three failures, three DIFFERENT causes, and the direction "
            "differs too -- the descriptor arm reads high in one and low in "
            "two. A test that flags mismatch does not tell you which arm is "
            "right; that needs looking, and each was traced rather than "
            "labelled.\n"
        )
        for m in d["arms_disagree"]:
            r = next(x for x in d["rows"] if x["modality"] == m)
            L.append(f"**`{m}`** — text {r['text_trial_share']}%, descriptor "
                     f"{r['mesh_trial_share']}%, {r['arm_ratio']}x. "
                     + DIAGNOSES.get(m, "Not yet traced.") + "\n")
        L.append(
            "The other mechanisms pass, which is what makes a failure "
            "interpretable rather than a property of the method: HIFU's two "
            "arms agree to within a tenth of a point, and the descriptor arm "
            "is not systematically higher or lower across the set.\n"
        )
    L.append("## Validity is a property of the descriptor AND the question\n")
    L.append(
        "The framing above -- \"`Ultrasonic Therapy` is too broad\" -- is "
        "not quite right, and the correction matters because this book's "
        "central thesis figure rests on the same descriptor.\n"
    )
    L.append(
        "Section 8.2 counts 32 census articles carrying both `Ferroptosis` and "
        "`Ultrasonic Therapy`, and reads them as the sonodynamic leg of this "
        "project's own hypothesis. If that descriptor were simply borrowed, "
        "that number would be borrowed too. **It is not: 30 of the 32 (94%) "
        "mention sonodynamic therapy or a sonosensitiser in title or abstract, "
        "and NONE mentions hyperthermia.** The same descriptor that fails "
        "badly on trial share passes cleanly here.\n"
    )
    L.append(
        "The reason is measurable rather than hand-waved. The contaminating "
        "literature is OLD -- the 114 trials run 1981 to 2025 with a median of "
        "2005 and 94% predating 2020 -- while `Ferroptosis` was minted in 2020 "
        "and its intersection with `Ultrasonic Therapy` runs 2021 to 2026 with "
        "a median of 2025 and NOTHING before 2020. Intersecting with a young "
        "descriptor imposes an era filter that removes exactly the "
        "contamination.\n"
    )
    L.append(
        "So descriptor validity is not a property a descriptor has or lacks. "
        "It is a property of the descriptor AND the question asked of it, and "
        "a single verdict attached to a descriptor name would be wrong in one "
        "of these two directions whichever way it was written. The test has to "
        "be re-run per question, which is cheap, and this project's "
        "long-standing prose caveat -- \"that descriptor is broad\" -- was "
        "not merely imprecise but unattached to the question it was "
        "qualifying.\n"
    )
    L.append("## The comparison Section 9.4 turns on\n")
    rl = next((r for r in d["rows"] if r["modality"] == "radioligand-therapy"),
              None)
    sdt = next((r for r in d["rows"] if r["modality"] == "sonodynamic"), None)
    if rl and sdt:
        L.append(
            f"The section argues that if radioligand therapy meets the same "
            f"clinical need through a more mature path, the physical-ROS "
            f"hypothesis is academically interesting and clinically "
            f"unnecessary. On these two axes radioligand therapy is further "
            f"along: {rl[f'{arm}_trial_share']}% of its literature carries a "
            f"trial publication type against sonodynamic therapy's "
            f"{sdt[f'{arm}_trial_share']}%, and it has "
            f"{rl[f'{arm}_late_phase']:,} Phase III or IV records against "
            f"{sdt[f'{arm}_late_phase']:,}.\n"
        )
        L.append(
            "That is a fact about development stage, and the section's "
            "argument needs more than it. Being further along a path does not "
            "establish meeting the same need: radioligand therapy requires a "
            "targetable receptor, and the residual-disease setting this "
            "project is concerned with is defined by a metabolic state rather "
            "than by receptor expression. Whether those coincide is a "
            "biological question the counts do not touch, and it is precisely "
            "the applicability axis this analysis cannot measure.\n"
        )
    L.append("## What a trial share is not\n")
    L.append(
        "Evidence that a modality WORKS. It counts trials that happened, not "
        "trials that read out well, and a modality with a low barrier to a "
        "first-in-human study accumulates them without any of them "
        "succeeding. The late-phase share is the closer proxy for progress "
        "along a development path, because reaching Phase III generally "
        "requires earlier phases to have been survived -- but a literature can "
        "also stall at Phase II for a decade and this table will not show it.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only
                 else scan(a.stride))
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    for r in d["rows"]:
        print(f"  {r['modality']:26s} text {r['text_trial_share']}% "
              f"({r['text_trials']:,} trials)  mesh {r['mesh_trial_share']}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
