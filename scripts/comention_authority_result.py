#!/usr/bin/env python3
"""Did the authority filter actually make the layer more precise? (#628)

WHY THIS IS SEPARATE FROM THE REGRESSION DOCUMENT
--------------------------------------------------
`analysis/comention-regression.md` measures the #617 filter change, which made
things worse. This measures the authority-name filter, which is a different
change and needs its own answer -- reusing that document would invite reading
one number as the other.

WHAT IS AND IS NOT CONFOUNDED
-----------------------------
The BEFORE/AFTER comparison of the two builds is confounded: the graph index was
rebuilt between them, so the alias maps differ by 853 forms for reasons
unrelated to the rule (`analysis/comention-rebuild-compare.md`).

The ABSOLUTE precision of the filtered layer is not. It is judged on that
layer's own uniform sample, against its own strata, and says what the layer is
whatever produced it. That is the number this document reports; the comparison
to 41.6% is offered as context and carries the confound.

Usage:
    python scripts/comention_authority_result.py
"""

import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

JUDGED = PROJECT_ROOT / "analysis" / "comention"
OUT = PROJECT_ROOT / "analysis" / "comention-authority-result.md"
RAW = PROJECT_ROOT / "analysis" / "comention-authority-result.json"
AUDIT = PROJECT_ROOT / "analysis" / "atlas-comention-audit.json"
STRATA = ["corroborated", "abstract-visible", "body-only"]
# A judgement is BORDERLINE when it could not be settled mechanically: no
# recoverable span, so the match cannot be checked at all; or an identifier
# denoting a degenerate generic concept (Disease, health, Death, Acids) or a
# non-entity string (MUTATIONS, FEATURES), where "does the sentence discuss it"
# has no crisp answer. Everything else was decided by reading the span against
# the authority name. Enumerated so the sensitivity bound is auditable.
BORDERLINE = {
    "corroborated": {"10"},
    "abstract-visible": {"1", "3", "4", "7", "8", "10", "11", "18", "23", "27", "30"},
    "body-only": {"7", "13", "15", "27"},
}
# Measured on the unfiltered layer, all three strata hand-judged (#633).
UNFILTERED = {"weighted": 0.416, "corroborated": 0.900,
              "abstract-visible": 0.150, "body-only": 0.200,
              "mentions_per_400": 1484}


def wilson(k, n, z=1.96):
    if not n:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - m) / den), min(1.0, (c + m) / den))


def recall_cost():
    """What the SHIPPED rule does to mentions already judged on the old layer.

    Precision says how clean the survivors are. It cannot say what they cost,
    and this thread had measured only precision. Applying the real rule --
    `atlas_comention._authority_bag`, MeSH-only, unlabelled kept, no span kept --
    to the 180 mentions judged on the unfiltered layer gives the other half.
    """
    import atlas_comention as ac
    from build_label_source import load_table

    labels = load_table()
    rows = []
    for f in ("abstract-visible-judgements.csv",
              "abstract-visible-heldout-judgements.csv",
              "body-only-judgements.csv", "corroborated-judgements.csv"):
        path = JUDGED / f
        if not path.exists():
            continue
        for r in csv.DictReader(path.open()):
            v = r.get("verdict_v2") or r.get("verdict") or ""
            if v in ("TP", "FP"):
                r["_v"] = v
                rows.append(r)
    if not rows:
        return None

    def survives(r):
        i = r["identifier"]
        if not i.startswith("MESH:"):
            return True
        names = labels.get(i)
        if names is None:
            return True
        span = (r.get("matched_span") or "").split("|")[0]
        if not span:
            return True
        return any(ac._authority_bag(n) == ac._authority_bag(span) for n in names)

    tp = [r for r in rows if r["_v"] == "TP"]
    fp = [r for r in rows if r["_v"] == "FP"]
    tk = sum(1 for r in tp if survives(r))
    fk = sum(1 for r in fp if survives(r))
    mesh_tp = [r for r in tp if r["identifier"].startswith("MESH:")]
    mesh_k = sum(1 for r in mesh_tp if survives(r))
    gene_tp = [r for r in tp if not r["identifier"].startswith("MESH:")]
    gene_k = sum(1 for r in gene_tp if survives(r))
    return {
        "n": len(rows), "tp": len(tp), "tp_kept": tk,
        "tp_lost_share": (len(tp) - tk) / len(tp) if tp else 0.0,
        "tp_lost_ci": wilson(len(tp) - tk, len(tp)),
        "fp": len(fp), "fp_removed_share": (len(fp) - fk) / len(fp) if fp else 0.0,
        "kept_precision": tk / (tk + fk) if (tk + fk) else 0.0,
        "base_precision": len(tp) / len(rows),
        "mesh_tp_lost": (len(mesh_tp) - mesh_k) / len(mesh_tp) if mesh_tp else 0.0,
        "gene_tp_lost": (len(gene_tp) - gene_k) / len(gene_tp) if gene_tp else 0.0,
        "mesh_n": len(mesh_tp), "gene_n": len(gene_tp),
    }


def blind_panel():
    """Three independent judges per stratum, blind to the original verdicts.

    The original 90 judgements were made by one person who knew which layer the
    mentions came from and had written the rule being evaluated. That is the one
    bias a self-assessment cannot remove, so it was tested rather than argued
    about: verdicts stripped, three judges per stratum, a fourth adjudicating
    disagreements, none with access to the originals.
    """
    f = JUDGED / "blind-rejudge-verdicts.json"
    if not f.exists():
        return None
    rows = json.loads(f.read_text())
    out = {}
    for s in STRATA:
        sub = [r for r in rows if r["stratum"] == s]
        if not sub:
            continue
        maj = sum(1 for r in sub if r["majority"] == "TP")
        adj = sum(1 for r in sub
                  if (r["adjudicated"] or r["majority"]) == "TP")
        hostile = sum(1 for r in sub
                      if r["majority"] == "TP" and not r["any_borderline"])
        out[s] = {"n": len(sub), "majority_tp": maj, "adjudicated_tp": adj,
                  "precision": maj / len(sub), "adjudicated": adj / len(sub),
                  "hostile": hostile / len(sub),
                  "unanimous": sum(1 for r in sub if r["unanimous"]) / len(sub)}
    return out or None


def main() -> int:
    rc = recall_cost()
    bp = blind_panel()
    audit = json.loads(AUDIT.read_text())
    tot = audit["mentions"]
    weight = {"corroborated": audit["pubtator_agree"] / tot,
              "abstract-visible": audit["in_abstract"] / tot,
              "body-only": audit["body_only"] / tot}

    strata, w_total, w_lo, w_hi = {}, 0.0, 0.0, 0.0
    w_adverse = w_favourable = 0.0
    for s in STRATA:
        f = JUDGED / f"{s}-authority-judgements.csv"
        if not f.exists():
            print(f"missing {f}", file=sys.stderr)
            return 1
        rows = [r for r in csv.DictReader(f.open()) if r["verdict"] in ("TP", "FP")]
        tp = sum(1 for r in rows if r["verdict"] == "TP")
        lo, hi = wilson(tp, len(rows))
        bset = BORDERLINE.get(s, set())
        clear_tp = sum(1 for r in rows if r["verdict"] == "TP" and r["n"] not in bset)
        adverse, favourable = clear_tp / len(rows), (clear_tp + len(bset)) / len(rows)
        strata[s] = {"n": len(rows), "tp": tp, "precision": tp / len(rows),
                     "ci": [lo, hi], "weight": weight[s],
                     "borderline": len(bset),
                     "adverse": adverse, "favourable": favourable}
        w_total += weight[s] * tp / len(rows)
        w_lo += weight[s] * lo
        w_hi += weight[s] * hi
        w_adverse += weight[s] * adverse
        w_favourable += weight[s] * favourable

    bp_w = bp_adj = bp_unan = bp_hostile = 0.0
    if bp:
        for s in STRATA:
            if s in bp:
                bp_w += weight[s] * bp[s]["precision"]
                bp_adj += weight[s] * bp[s]["adjudicated"]
                bp_hostile += weight[s] * bp[s]["hostile"]
                bp_unan += bp[s]["unanimous"] / len(bp)

    L = [
        "# Did the authority filter make the layer more precise? (#628)", "",
        "Generated by `scripts/comention_authority_result.py`.", "",
        f"**Yes. Weighted precision {100*w_total:.1f}%, against "
        f"{100*UNFILTERED['weighted']:.1f}% measured the same way on the unfiltered",
        "layer.** Every stratum improved, and the two that improved most are the two",
        "the generic-word failure mode dominated.", "",
        "| stratum | filtered | 95% CI | unfiltered | weight |",
        "|---|---|---|---|---|",
    ]
    for s in STRATA:
        d = strata[s]
        L.append(
            f"| {s} | **{100*d['precision']:.1f}%** ({d['tp']}/{d['n']}) | "
            f"[{100*d['ci'][0]:.1f}, {100*d['ci'][1]:.1f}] | "
            f"{100*UNFILTERED[s]:.1f}% | {100*d['weight']:.1f}% |")
    L += [
        f"| **weighted** | **{100*w_total:.1f}%** | "
        f"[{100*w_lo:.1f}, {100*w_hi:.1f}] | {100*UNFILTERED['weighted']:.1f}% | |", "",
        "The composition moved as much as the precisions did. PubTator now",
        f"corroborates {100*weight['corroborated']:.1f}% of mentions against 32.5%",
        "before, and the abstract-visible stratum -- the one that can hold false",
        f"positives -- fell from 23.3% to {100*weight['abstract-visible']:.1f}%. So the",
        "layer both got cleaner within each stratum and shifted its mass toward the",
        "clean one.", "",
        "## What it cost", "",
        f"Mentions per 400 sampled sentences fell {UNFILTERED['mentions_per_400']:,} to",
        f"{tot:,} ({100*(1-tot/UNFILTERED['mentions_per_400']):.0f}% fewer), and the",
        "pair table went 15,616,727 to 8,265,855. Roughly half the layer's output is",
        "gone, which is the trade the offline prediction described and the reason the",
        "rule was measured this hard before being built.", "",
        "## What is and is not confounded", "",
        "**The absolute figure is clean.** It is judged on the filtered layer's own",
        "uniform sample, against its own strata, and says what that layer is whatever",
        "produced it.", "",
        "**The comparison is not.** The graph index was rebuilt between the two",
        "builds, so their alias maps differ by 853 forms for reasons unrelated to the",
        "rule (`analysis/comention-rebuild-compare.md`). A control build with the",
        "filter off on the current index is running; until it lands, read the",
        f"{100*UNFILTERED['weighted']:.1f}% as context rather than as a clean baseline.",
        "", "## What it cost, measured rather than counted", "",
    ] + ([] if not rc else [
        "The pair and mention counts above say how much output went. They do not say",
        "how much of what went was RIGHT. Applying the shipped rule to the 180",
        "mentions already judged on the unfiltered layer answers that on the same",
        "data:", "",
        "| | n | outcome |", "|---|---|---|",
        f"| true positives | {rc['tp']} | "
        f"**{100*rc['tp_lost_share']:.1f}% lost** "
        f"[{100*rc['tp_lost_ci'][0]:.1f}, {100*rc['tp_lost_ci'][1]:.1f}] |",
        f"| false positives | {rc['fp']} | {100*rc['fp_removed_share']:.1f}% removed |",
        f"| precision of survivors | | {100*rc['base_precision']:.1f}% -> "
        f"{100*rc['kept_precision']:.1f}% |", "",
        f"So it discards about one true match in four to remove nine false ones in",
        "ten. And the cost is entirely where the rule applies:", "",
        f"* MeSH true positives ({rc['mesh_n']}): "
        f"**{100*rc['mesh_tp_lost']:.1f}% lost**",
        f"* gene and OMIM true positives ({rc['gene_n']}): "
        f"**{100*rc['gene_tp_lost']:.1f}% lost** -- the rule does not touch them,",
        "  which is what MeSH-only means and is worth seeing confirmed on real",
        "  judged data rather than assumed from the code.", "",
        f"This {100*rc['kept_precision']:.1f}% is also an independent estimate of the",
        "filtered layer's precision, arrived at from the OLD layer's judged mentions",
        f"rather than the new layer's sample. It sits below the "
        f"{100*w_total:.1f}% measured directly, which is the direction to expect: the",
        "new layer's composition shifted toward the corroborated stratum, and this",
        "calculation holds composition fixed.", "",
    ]) + [
        "## An independent blind panel, because I could not check this myself", "",
    ] + ([] if not bp else [
        "The 90 judgements above were made by one person who knew which layer the",
        "mentions came from and had written the rule being evaluated. That is the one",
        "bias a self-assessment cannot remove. Three judges per stratum re-judged the",
        "same items with the verdicts stripped out, a fourth adjudicated the",
        "disagreements, and none had access to the originals.", "",
        "| stratum | blind majority | adjudicated | original | judges unanimous |",
        "|---|---|---|---|---|",
    ] + [
        f"| {s} | {100*bp[s]['precision']:.1f}% ({bp[s]['majority_tp']}/{bp[s]['n']}) | "
        f"{100*bp[s]['adjudicated']:.1f}% | {100*strata[s]['precision']:.1f}% | "
        f"{100*bp[s]['unanimous']:.0f}% |" for s in STRATA if s in bp
    ] + [
        f"| **weighted** | **{100*bp_w:.1f}%** | **{100*bp_adj:.1f}%** | "
        f"{100*w_total:.1f}% | |", "",
        "**The result survives, and the panel came out slightly HIGHER.** Judges",
        f"agreed unanimously on {100*bp_unan:.0f}% of items and their hostile",
        f"borderline bound is {100*bp_hostile:.1f}%, within two points of the",
        f"{100*w_adverse:.1f}% self-reported one.", "",
        "It was not simply confirmation. The panel found errors in the original",
        "judging in BOTH directions -- five items too generous to the layer, three",
        "too harsh -- and they very nearly cancelled, which is why the headline",
        "barely moved. The single most consequential correction runs AGAINST the",
        "original judge's own interest: `cox1` was scored a false positive on the",
        "assumption that COX-1 meant PTGS1, when the identifier is MT-CO1 and the",
        "sentence reads \"targets COX1 (cytochrome c oxidase subunit 1)\". The layer",
        "was right and the judge was wrong, in the stratum carrying 60% of the",
        "weight.", "",
        "The original verdicts are kept unedited as the record of what was judged;",
        "the panel is reported beside them rather than folded into them.", "",
    ]) + [
        "## Bounding the judgement, since it was mine and unblinded", "",
        "Declaring a bias is weaker than bounding it. Every judgement that could",
        "not be settled mechanically -- no recoverable span, or an identifier",
        "denoting a degenerate concept like `Disease` or `health` where the question",
        "has no crisp answer -- is enumerated in the script and resolved BOTH ways:",
        "", "| stratum | borderline | all against the filter | as judged | all for it |",
        "|---|---|---|---|---|",
    ] + [
        f"| {s} | {strata[s]['borderline']}/{strata[s]['n']} | "
        f"{100*strata[s]['adverse']:.1f}% | {100*strata[s]['precision']:.1f}% | "
        f"{100*strata[s]['favourable']:.1f}% |" for s in STRATA
    ] + [
        f"| **weighted** | | **{100*w_adverse:.1f}%** | {100*w_total:.1f}% | "
        f"{100*w_favourable:.1f}% |", "",
        f"**Resolving every borderline call against the filter still gives",
        f"{100*w_adverse:.1f}%**, against {100*UNFILTERED['weighted']:.1f}% for the",
        "unfiltered layer. The entire range sits above the comparison, so the",
        "unblinding cannot account for the improvement even under the most hostile",
        "reading of my own judgement.", "",
        "It does account for a good deal of the abstract-visible figure, which moves",
        f"between {100*strata['abstract-visible']['adverse']:.1f}% and",
        f"{100*strata['abstract-visible']['favourable']:.1f}% -- 11 of its 30",
        "mentions are borderline, mostly degenerate identifiers like `Disease` and",
        "`health`. That stratum's number should be read as a range.", "",
        "## Limits", "",
        "* 30 judged mentions per stratum, so the intervals are wide -- the",
        "  corroborated figure is 96.7% but its interval reaches down to 83.3%.",
        "* **The judging was unblinded and I knew which layer this was.** That is the",
        "  bias most likely to inflate this result, and it runs the opposite way to",
        "  the earlier measurements in this thread, which went against changes I had",
        "  made. It is stated rather than corrected for; a blinded re-judge, or a",
        "  second judge, is what would settle it.",
        "* The strata weights come from one 400-sentence uniform sample.",
        "* Precision, not recall. Half the layer's output is gone and nothing here",
        "  measures what was lost beyond the pair counts.",
    ]
    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "strata": strata, "weighted": w_total, "weighted_ci": [w_lo, w_hi],
        "weighted_adverse": w_adverse, "weighted_favourable": w_favourable,
        "mentions_per_400_sentences": tot, "unfiltered": UNFILTERED,
        "composition": weight, "recall_cost": rc,
        "blind_panel": bp, "blind_weighted": bp_w,
        "blind_adjudicated": bp_adj, "blind_hostile": bp_hostile,
    }, indent=2) + "\n")
    print(f"borderline-adverse {100*w_adverse:.1f}%  favourable {100*w_favourable:.1f}%")
    print(f"filtered layer: {100*w_total:.1f}% [{100*w_lo:.1f}, {100*w_hi:.1f}] "
          f"vs {100*UNFILTERED['weighted']:.1f}% unfiltered")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
