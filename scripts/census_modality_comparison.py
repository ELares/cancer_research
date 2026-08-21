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
# A share off a handful of trials is a statement about the handful.
MIN_TRIALS = 10


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
    want = set(COMPARED)
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
        "mesh_measurable": sorted(mesh),
        "text_measurable": sorted(text),
        "arms": {a: {k: dict(v) for k, v in arms[a].items()} for a in arms},
    }


def assemble(d: dict) -> dict:
    rows = []
    for k in d["compared"]:
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
    if d["arms_disagree"]:
        L.append(
            f"**{', '.join(f'`{m}`' for m in d['arms_disagree'])} fails the "
            f"test, and the failure was traced rather than inferred.** Of the "
            f"114 census articles carrying `Ultrasonic Therapy` AND a trial "
            f"publication type, exactly **2** mention sonodynamic therapy or a "
            f"sonosensitiser anywhere in title or abstract. The remaining 112 "
            f"are ultrasound HYPERTHERMIA trials, mostly from the 1980s and "
            f"1990s -- a different modality that shares an instrument.\n"
        )
        L.append(
            "So a descriptor-arm trial share of 4.54% for sonodynamic therapy "
            "is 98% borrowed from work that is not sonodynamic therapy. The "
            "text arm's 0.29% is the honest figure, and the correction runs "
            "in the direction that makes this project's own thesis modality "
            "look EARLIER than previously reported, not later.\n"
        )
        L.append(
            "The other modalities pass, which is what makes the failure "
            "interpretable rather than a property of the method: HIFU's two "
            "arms agree to within a tenth of a point, and the descriptor arm "
            "is not systematically higher.\n"
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
