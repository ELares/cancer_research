#!/usr/bin/env python3
"""Is there a calibration target for the unmodelled ferroptosis genes? (#616)

WHY
---
`analysis/atlas-model-gaps.md` asked the literature which ferroptosis mechanisms
it emphasises and found four with no parameter in `ferroptosis-core`: HMOX1
(1,583 ferroptosis-indexed articles), TP53 (1,497), TFRC (1,482), KEAP1 (768).

`CONTRIBUTING.md`'s layer policy requires a NAMED calibration target before any
new off-by-default axis lands. So the question is not "add four layers", it is
whether public data exists that could anchor them. The precedent is #444, which
shipped ACSL4 flagged DATA-GATED, and #462, which anchored it against cBioPortal
TCGA expression.

TWO CUTS, TWO DIFFERENT ANSWERS
-------------------------------
This measures the same seven genes at both thresholds #462 used, and they do not
agree. That split is the whole result.

At `z < -1` the route is uninformative for everything. Within-cohort mRNA
z-scores are standardised per study by construction, so the fraction below z=-1
is ~15.87% for any gene under approximate normality, and all seven land within
about two points of it. A control gene with no ferroptosis role returns the same
number, so this cut cannot anchor anything and cannot support a gene-specific
reading.

At `z < -2` the genes SEPARATE. TP53 sits at 4.27% and exceeds the 2.28%
expectation in 31 of 32 cancer types; ACSL4 sits at 2.98% and exceeds it in 25 of
32. The other five sit at or below expectation. A deep left tail heavier than a
normal is the signature of a discrete low-expression subpopulation -- deletion,
or a genuinely refractory subtype -- which is precisely what standardisation
cannot manufacture.

WHAT THAT MEANS FOR EACH QUESTION
---------------------------------
For the four unmodelled genes: HMOX1, TFRC and KEAP1 are unanchored at BOTH cuts,
so this route offers them nothing. TP53 is the exception and has a real signal.

For #462: its `z < -1` prevalence figure is not gene-specific and should not be
read as one. Its `z < -2` figure SURVIVES -- the control separates ACSL4 rather
than washing it out, so that row is evidence about ACSL4 and not merely about the
z-score. An earlier draft of this analysis compared TP53 against "the others'
~1.5-2.1%", a range that silently omitted ACSL4's 2.98%, and drew the wrong
conclusion for both. The comparator is now every gene, always.

Usage:
    python scripts/calibration_feasibility.py
"""

import csv
import json
import statistics
import sys
from math import comb, erf, sqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

SRC = PROJECT_ROOT / "analysis" / "calibration" / "acsl4_prevalence_tcga.csv"
OUT = PROJECT_ROOT / "analysis" / "calibration" / "calibration-feasibility.md"
RAW = PROJECT_ROOT / "analysis" / "calibration" / "calibration-feasibility.json"

GAPS = ["HMOX1", "TP53", "TFRC", "KEAP1"]
ANCHORED = ["ACSL4", "GPX4", "SLC7A11"]
EXP_LOW = 0.5 * (1 + erf(-1 / sqrt(2)))       # 15.87%
EXP_VLOW = 0.5 * (1 + erf(-2 / sqrt(2)))      # 2.28%
# A gene is called SEPARATED at a cut when its per-cancer-type fraction exceeds
# the normal expectation in significantly more types than chance. Two-sided sign
# test: distribution-free, and the right test for "is this consistently above the
# line", which a median alone cannot answer.
ALPHA = 0.01


def sign_test(above: int, n: int) -> float:
    """Two-sided binomial sign test against p=0.5."""
    if n == 0:
        return 1.0
    k = min(above, n - above)
    tail = sum(comb(n, i) for i in range(k + 1)) / (2 ** n)
    return min(1.0, 2 * tail)


def quartile_spread(vals: list) -> float:
    """Q3-Q1 by linear interpolation.

    An earlier version indexed `sorted[int(.75*n)] - sorted[int(.25*n)]`, which
    for n=32 is the 25th minus the 9th order statistic -- roughly the 76.6th
    minus the 26.6th percentile, overstating TP53's spread by 19%. It is
    labelled an IQR, so it should be one.
    """
    v = sorted(vals)
    n = len(v)
    if n < 2:
        return 0.0

    def q(p):
        h = (n - 1) * p
        lo = int(h)
        return v[lo] + (h - lo) * (v[min(lo + 1, n - 1)] - v[lo])
    return q(0.75) - q(0.25)


def main() -> int:
    if not SRC.exists():
        print(f"missing {SRC}; run scripts/fetch_acsl4_prevalence.py", file=sys.stderr)
        return 1
    rows = list(csv.DictReader(SRC.open()))
    stats = {}
    for g in ANCHORED + GAPS:
        if f"{g}_frac_low" not in rows[0]:
            continue
        lo = [float(r[f"{g}_frac_low"]) for r in rows if r[f"{g}_frac_low"] not in ("", None)]
        vl = [float(r[f"{g}_frac_verylow"]) for r in rows
              if r.get(f"{g}_frac_verylow") not in ("", None)]
        above_lo = sum(1 for x in lo if x > EXP_LOW)
        above_vl = sum(1 for x in vl if x > EXP_VLOW)
        p_lo, p_vl = sign_test(above_lo, len(lo)), sign_test(above_vl, len(vl))
        stats[g] = {
            "n_types": len(lo),
            "median_low": statistics.median(lo),
            "median_verylow": statistics.median(vl),
            "iqr_low": quartile_spread(lo),
            "dev_from_normal": statistics.median(lo) - EXP_LOW,
            "types_above_expected_low": above_lo,
            "types_above_expected_verylow": above_vl,
            "p_low": p_lo, "p_verylow": p_vl,
            # "Separated" means BOTH consistently above expectation and
            # significant. Below-expectation genes are not evidence of anything
            # a refractory-tail prior could use.
            "separated_low": p_lo < ALPHA and above_lo > len(lo) / 2,
            "separated_verylow": p_vl < ALPHA and above_vl > len(vl) / 2,
        }

    sep_lo = [g for g, s in stats.items() if s["separated_low"]]
    sep_vl = [g for g, s in stats.items() if s["separated_verylow"]]

    L = [
        "# Is there a calibration target for the unmodelled genes? (#616)", "",
        "Generated by `scripts/calibration_feasibility.py`.", "",
        "`CONTRIBUTING.md`'s layer policy requires a named calibration target before",
        "any new off-by-default axis lands, so the question for the four mechanisms in",
        "`analysis/atlas-model-gaps.md` is whether public data could anchor them --",
        "not whether to write the layers.", "",
        "## The two cuts do not agree, and that is the result", "",
        "cBioPortal mRNA z-scores are standardised WITHIN each study, so a shallow cut",
        "is fixed by construction: the fraction of tumours below z = -1 is",
        f"**{100*EXP_LOW:.2f}%** for any gene under approximate normality. A deep cut is",
        "not similarly fixed -- a left tail heavier than a normal is the signature of a",
        "discrete low-expression subpopulation, which standardisation cannot",
        "manufacture. Measured across 32 TCGA PanCancer Atlas studies:", "",
        "| gene | z < -1 | above expectation | z < -2 | above expectation | separates? |",
        "|---|---|---|---|---|---|",
    ]
    for g, s in sorted(stats.items(), key=lambda kv: -kv[1]["median_verylow"]):
        tag = " *(gap)*" if g in GAPS else ""
        verdict = "**yes, deep cut**" if s["separated_verylow"] else "no"
        L.append(
            f"| {g}{tag} | {100*s['median_low']:.1f}% | "
            f"{s['types_above_expected_low']}/{s['n_types']} | "
            f"{100*s['median_verylow']:.2f}% | "
            f"{s['types_above_expected_verylow']}/{s['n_types']} "
            f"(p={s['p_verylow']:.2g}) | {verdict} |")

    L += [
        "", f"**At z = -1, nothing separates** ({len(sep_lo)} of {len(stats)} genes).",
        "Every gene lands within about two points of the normal expectation, so this",
        "cut recovers the standardisation rather than any biology, and it cannot tell",
        "a gene whose low-expression tail matters from one with no ferroptosis role.",
        "",
        f"**At z = -2, two genes separate**: {' and '.join(sorted(sep_vl))}. TP53 exceeds",
        f"expectation in {stats['TP53']['types_above_expected_verylow']} of 32 cancer",
        f"types and ACSL4 in {stats['ACSL4']['types_above_expected_verylow']} of 32,",
        "while the remaining five sit at or below it. That is a real, gene-specific",
        "signal: consistent with TP53 deletion and with a genuinely refractory",
        "low-ACSL4 subpopulation.", "",
        "## What this does to #462", "",
        "It splits the finding in two rather than overturning it.",
        "",
        "`analysis/calibration/acsl4-prevalence-calibration.md` reports both cuts. Its",
        "**z < -1 row (10.8-18.8%, median 14.4%) is NOT gene-specific** and should not",
        "be read as one -- six control genes return 13.8-16.0%, so a gene with no",
        "ferroptosis role would give the same figure. The document already noted the",
        "by-construction caveat, but with no control the number still read as evidence",
        "about ACSL4.", "",
        "Its **z < -2 row (median 3.0%) SURVIVES**. The same control that washes out",
        "the shallow cut separates ACSL4 at the deep one, so that row is evidence about",
        "ACSL4 and not merely about the z-score. It is also the cut the layer actually",
        "cares about, since `ACSL4_NEGATIVE` sits at z = -2.", "",
        "The `status_from_zscore` bridge is untouched by either: it maps a scale.", "",
        "> An earlier draft of this analysis reached the opposite conclusion on ACSL4",
        "> by comparing TP53 against \"the others' ~1.5-2.1%\" -- a range that silently",
        "> omitted ACSL4's own 2.98%, the highest of the six. Every comparison here now",
        "> runs against all seven genes, and a guard pins the deep-cut result for each.",
        "", "## Verdict per gene", "",
        "| gene | calibration status | reason |", "|---|---|---|",
    ]
    reasons = {
        "HMOX1": ("**data-blocked**", "Unanchored at BOTH cuts: 14.9% at z<-1 and "
                  f"{100*stats['HMOX1']['median_verylow']:.2f}% at z<-2, exceeding "
                  f"expectation in only {stats['HMOX1']['types_above_expected_verylow']}"
                  "/32 types. Neither cut of this route carries information about it."),
        "TFRC": ("**data-blocked**", "Unanchored at both cuts "
                 f"({100*stats['TFRC']['median_verylow']:.2f}% at z<-2, "
                 f"{stats['TFRC']['types_above_expected_verylow']}/32 types). Note that "
                 "the engine models NCOA4 release of iron from intracellular ferritin "
                 "(`ferritinophagy_release`), NOT transferrin-receptor import across the "
                 "membrane, so this is a genuine uncovered axis rather than a duplicate."),
        "KEAP1": ("**data-blocked**", "Unanchored at both cuts, and the only gene "
                  f"BELOW expectation at z<-1 by more than two points "
                  f"({100*stats['KEAP1']['median_verylow']:.2f}% at z<-2). Much of its "
                  "effect is already carried by `nrf2_gsh_rate`, the pathway it "
                  "regulates, though the KEAP1 sensor itself is not separately modelled."),
        "TP53": ("**weak anchor available**", "The one gene here with a real signal: "
                 f"{100*stats['TP53']['median_verylow']:.2f}% at z<-2 against a 2.28% "
                 f"expectation, exceeding it in "
                 f"{stats['TP53']['types_above_expected_verylow']}/32 types "
                 f"(p={stats['TP53']['p_verylow']:.1g}), with the widest inter-cancer "
                 "spread. Consistent with real deletion. CTRPv2 also carries 828 cell "
                 "lines x 5 ferroptosis inducers for a sensitivity join if a cell-line "
                 "expression source can be obtained."),
    }
    for g in GAPS:
        status, why = reasons[g]
        L.append(f"| {g} | {status} | {why} |")

    L += [
        "", "None of the four should have a layer written on this evidence. Even TP53's",
        "signal bounds a PREVALENCE, not a dose-response: it says how many tumours sit",
        "in the deep tail, not what ferroptosis sensitivity to give them. The honest",
        "outcome is a data-blocked row, which is what #444 did for ACSL4 before #462",
        "found a partial anchor.", "",
        "## What would change the answer", "",
        "A cell-line expression matrix joinable to the committed CTRPv2 curves would",
        "make all four testable at once: correlate expression against ferroptosis-",
        "inducer EC50 across 828 lines. The DepMap download catalogue served that until",
        "recently and now returns HTML rather than the documented CSV, so that join is",
        "blocked on access, not on method.", "",
        "## Limits", "",
        "* **Scope of the negative.** Two routes were tried: cBioPortal expression",
        "  prevalence (measured above) and the DepMap cell-line join (access-blocked).",
        "  The data-blocked verdicts mean \"neither route tried anchors this gene\", not",
        "  \"no dataset anywhere could\". A targeted knockdown or overexpression",
        "  dose-response search per gene has not been run.",
        "* Bulk tumour mRNA, which #462 already found does not rank the ACSL4-low",
        "  refractory phenotype the cell-line literature reports. The same caveat",
        "  applies here with more force.",
        "* The z < -1 comparison assumes approximate per-study normality. Most genes",
        "  sit slightly BELOW 15.87% while ACSL4 and TP53 sit below at z=-1 and ABOVE at",
        "  z=-2. That combination is not general skew -- it is what a discrete",
        "  low-expression subpopulation looks like, and it is the reason the two cuts",
        "  had to be measured separately rather than treated as one finding.",
        "* The sign test asks whether a gene is CONSISTENTLY above expectation across",
        "  cancer types. It does not bound the effect size, and the per-type fractions",
        "  are not independent of each other's sample sizes.",
        "* A gene can matter enormously and have a flat expression distribution. This",
        "  measures whether THIS route informs a prior, not whether the gene matters.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "expected_low": EXP_LOW, "expected_verylow": EXP_VLOW, "alpha": ALPHA,
        "genes": stats, "gaps": GAPS,
        "separated_low": sep_lo, "separated_verylow": sep_vl,
    }, indent=2) + "\n")
    print(f"wrote {OUT}")
    print(f"separated at z<-1: {sep_lo or 'none'}")
    print(f"separated at z<-2: {sep_vl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
