#!/usr/bin/env python3
"""Did the #617 co-mention filters work? Measured answer: no (#ATLAS-COMENT-REG).

WHY THIS EXISTS
---------------
#617 measured the full-text co-mention layer at 55.5% population-weighted
precision and traced it to one gap: `usable_alias` exempted every MULTI-WORD
form from the specificity test it applied to single tokens, so `tumor cells`
resolved to Glucagonoma and `et al` to Multiple Myeloma. 132 of 152 measured
false positives were multi-word.

The fix removed the single-token shape proxy entirely and replaced it with two
measured filters, census support and share of an identifier's mentions. The
stated justification was that the proxy "kills the same 141 of 152 false
positives and retains strictly more true ones".

That was measured on the pre-fix sample. A fresh uniform sample after the
rebuild says the change made precision WORSE, and this script is the accounting.

WHAT IT FOUND
-------------
1. Net population-weighted precision fell by roughly 9 points once every
   stratum was hand-judged rather than assumed: 50.6% -> 41.6%. The exact
   levels depend on which body-only measurement is placed on the before side,
   and both bases are tabulated in the report; the DELTA is robust to that
   choice and the levels are not. The layer produces 1.33x more mentions per
   sentence and PubTator corroborates a smaller fraction (44.3% -> 32.5%).

2. The abstract-visible stratum -- the only one that can contain false
   positives -- nearly tripled as a share of all mentions, 8.2% -> 23.7%. Its
   own precision, judged at the IDENTIFIER level over 60 hand-read mentions, is
   15.0%; it did not measurably improve, and the stratum tripled. All three
   strata are now judged: corroborated 90.0% (36/40), body-only 20.0% (8/40).

3. The replacement filters DO NOT SEPARATE true from false matches on that
   stratum. Measured over the judged sample, support ranges 77-355,572 for
   false positives and 297-100,758 for true ones; share ranges 11.5%-274% for
   false and 20.2%-68.9% for true. The distributions overlap almost entirely,
   so no threshold on either filter recovers the loss. This is a design
   failure, not a tuning problem.

4. A concrete bug explains part of it. `alias_support` sums a form's mentions
   across every sense it carries; `ident_mentions` sums an identifier's across
   every form. Their ratio is not a share and can exceed 1 -- 274% for `as`,
   133% for `gp`, 104% for `tss`. It is also biased the wrong way: an ambiguous
   generic word collects a larger numerator from its other senses, so a
   MINIMUM-share filter admits it more readily than a specific name. The filter
   written to exclude `as` and `treatment` was structurally disposed to keep
   them. `atlas_graph.build_index` now records `alias_ident_support`, the
   numerator the ratio always needed -- but re-scoring the judged forms against
   the corrected share shows it does NOT rescue the filter: `as` still reads
   124.5%, and at the shipped 5% threshold the corrected share cuts 0 of 37
   false positives and 0 of 13 true ones. The filter is inert on this stratum
   either way.

THE GENERAL LESSON, WHICH IS THE POINT
--------------------------------------
The #617 justification was true about the sample it was measured on and wrong
about the population it was applied to. Removing the single-token test looked
free because the false positives of the day were multi-word -- but they were
multi-word BECAUSE the multi-word channel was unfiltered. Closing that channel
moved the pressure to the one that had just been opened.

Usage:
    python scripts/comention_regression.py
"""

import csv
import json
import math
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

AUDIT = PROJECT_ROOT / "analysis" / "atlas-comention-audit.json"
OUT = PROJECT_ROOT / "analysis" / "comention-regression.md"
RAW = PROJECT_ROOT / "analysis" / "comention-regression.json"

# Every stratum is now hand-judged at the identifier level. These were assumed
# values carried over from #617 (0.925 and 0.308) until #628 measured them, and
# the assumption mattered: body-only was optimistic by half.
JUDGED = {
    "abstract": PROJECT_ROOT / "analysis" / "comention" / "abstract-visible-judgements.csv",
    "heldout": PROJECT_ROOT / "analysis" / "comention" / "abstract-visible-heldout-judgements.csv",
    "body": PROJECT_ROOT / "analysis" / "comention" / "body-only-judgements.csv",
    "agree": PROJECT_ROOT / "analysis" / "comention" / "corroborated-judgements.csv",
}
PRIOR_ABSTRACT_PRECISION = 0.146  # #617's figure, lenient criterion


# Reference lists, table rows and PMC section markers. A co-mention harvested
# from a cited paper's title is recorded against the CITING pmid, which is a
# weaker claim than the layer is usually read as making.
_STRUCTURAL = re.compile(r"====|\t|\bdoi\b|10\.\d{4}/|\b(?:19|20)\d\d\b\s+\d+\s+\d+")


def structural_share(path) -> int:
    """How many TRUE positives in this stratum sit in structural text."""
    return sum(1 for r in csv.DictReader(path.open())
               if r.get("verdict") == "TP" and _STRUCTURAL.search(r["sentence"]))


def prose_precision(path, column="verdict"):
    """Precision restricted to the paper's own prose.

    A co-mention harvested from a cited reference's title is recorded against
    the CITING pmid, which is a weaker claim than the layer is read as making.
    Excluding those says whether the confound moves the headline.
    """
    rows = [r for r in csv.DictReader(path.open())
            if (r.get(column) or r.get("verdict_v2")) in ("TP", "FP")
            and not _STRUCTURAL.search(r["sentence"])]
    if not rows:
        return None
    tp = sum(1 for r in rows
             if (r.get(column) or r.get("verdict_v2")) == "TP")
    return {"n": len(rows), "tp": tp, "precision": tp / len(rows),
            "ci": wilson(tp, len(rows))}


def judged_precision(path, column="verdict"):
    """Hand-judged precision for a stratum, with its Wilson interval."""
    rows = [r for r in csv.DictReader(path.open()) if r.get(column)]
    tp = sum(1 for r in rows if r[column] == "TP")
    return tp / len(rows), len(rows), tp, wilson(tp, len(rows))


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - m) / den), min(1.0, (c + m) / den))


# The pre-rebuild audit, pinned to the commit it came from. Reading it from
# HEAD was wrong once the post-rebuild audit was committed: the script then
# compared a run against itself and refused to produce anything. The historical
# run cannot change, so it is a constant with its provenance attached.
PRE_REBUILD_COMMIT = "4389be80a1163ed294afa4ad2da71a8780bd6390"
PRE_REBUILD_AUDIT = {"mentions": 1112, "pubtator_agree": 493,
                     "in_abstract": 91, "body_only": 528}
# The AFTER side, pinned for the same reason the BEFORE side is.
#
# This document measures ONE historical change: the #617 filter swap. Its
# "after" used to be read live from `analysis/atlas-comention-audit.json`, which
# was correct exactly until that artifact was regenerated for a DIFFERENT build.
# It was, when the authority filter shipped (#646), and the audit now describes
# the filtered layer -- so re-running this generator silently recomputed a
# historical regression against a build it does not describe and reported
# precision RISING 10.9 points where the finding is that it FELL 9. Nothing
# failed; the document simply rewrote its own conclusion. Found by accident,
# months of commits after it became possible.
#
# Recovered from the audit as committed at ed6b10ec (#633), the last state that
# describes the post-#617 pre-authority build.
POST_617_AUDIT = {"mentions": 1484, "pubtator_agree": 482,
                  "in_abstract": 346, "body_only": 656}


def prior_audit():
    """The audit as it stood before the rebuild.

    Verified against git when the repository still has that commit, so the
    pinned numbers cannot drift silently; falls back to the constants when it
    does not (a shallow clone, for instance).
    """
    r = subprocess.run(
        ["git", "show", f"{PRE_REBUILD_COMMIT}:analysis/atlas-comention-audit.json"],
        capture_output=True, text=True, cwd=PROJECT_ROOT)
    if r.returncode == 0:
        from_git = json.loads(r.stdout)
        for k, v in PRE_REBUILD_AUDIT.items():
            if from_git.get(k) != v:
                raise SystemExit(
                    f"pinned pre-rebuild {k}={v} does not match commit "
                    f"{PRE_REBUILD_COMMIT[:12]} ({from_git.get(k)}); "
                    "the comparison baseline has moved")
    return dict(PRE_REBUILD_AUDIT)


def strata(a):
    n = a["mentions"]
    return {"agree": a["pubtator_agree"] / n,
            "abstract": a["in_abstract"] / n,
            "body": a["body_only"] / n}


def weighted(s, abstract_precision, agree_p, body_p):
    return (s["agree"] * agree_p
            + s["abstract"] * abstract_precision
            + s["body"] * body_p)


def main() -> int:
    if not AUDIT.exists() or not all(f.exists() for f in JUDGED.values()):
        print("run scripts/atlas_comention_audit.py first", file=sys.stderr)
        return 1
    # Both sides pinned: this is a historical comparison and neither end may
    # drift with whatever the live artifact currently describes.
    after = dict(POST_617_AUDIT)
    before = prior_audit()
    if before["mentions"] == after["mentions"]:
        print("the current audit equals the pre-rebuild one; nothing to compare",
              file=sys.stderr)
        return 1
    if AUDIT.exists():
        live = json.loads(AUDIT.read_text())
        if live["mentions"] != after["mentions"]:
            print(f"  note: the live audit describes a different build "
                  f"({live['mentions']:,} mentions vs the pinned "
                  f"{after['mentions']:,}); this document measures the #617 "
                  f"change and uses the pinned values.", file=sys.stderr)

    judged = list(csv.DictReader(JUDGED["abstract"].open()))
    n = len(judged)
    # verdict_v2 judges at the IDENTIFIER level: does the sentence discuss the
    # entity the identifier DENOTES, per NLM, rather than merely contain the
    # matched string. `verdict` is the first pass and is retained so the
    # correction below is checkable.
    key = "verdict_v2" if "verdict_v2" in judged[0] else "verdict"
    tp = sum(1 for r in judged if r[key] == "TP")
    tp_v1 = sum(1 for r in judged if r["verdict"] == "TP")
    kept = [r for r in judged if r.get("is_authority_name") == "True"]
    kept_tp = sum(1 for r in kept if r[key] == "TP")
    klo, khi = wilson(kept_tp, len(kept)) if kept else (0.0, 0.0)
    v1_prec = tp_v1 / n
    struct_tp, struct_body = structural_share(JUDGED["agree"]), structural_share(JUDGED["body"])
    ho = list(csv.DictReader(JUDGED["heldout"].open()))
    ho_n = len(ho)
    ho_tp = sum(1 for r in ho if r["verdict"] == "TP")
    ho_prec = ho_tp / ho_n
    ho_lo, ho_hi = wilson(ho_tp, ho_n)

    def _bag(s):
        return frozenset(w for w in re.split(r"[^a-z0-9]+", (s or "").lower()) if w)

    ho_keep = [r for r in ho if r["authority_name"]
               and _bag(r["entity"]) == _bag(r["authority_name"])]
    ho_keep_n = len(ho_keep)
    ho_keep_tp = sum(1 for r in ho_keep if r["verdict"] == "TP")
    ho_keep_prec = ho_keep_tp / ho_keep_n if ho_keep_n else 0.0
    ho_fp = ho_n - ho_tp
    ho_cut_fp = ((ho_fp - (ho_keep_n - ho_keep_tp)) / ho_fp) if ho_fp else 0.0
    ho_cut_tp = ((ho_tp - ho_keep_tp) / ho_tp) if ho_tp else 0.0
    ho_span_keep = sum(1 for r in ho if r["authority_name"] and
                       _bag((r["matched_span"] or "").split("|")[0]) == _bag(r["authority_name"]))
    abs_prec = tp / n
    lo, hi = wilson(tp, n)

    sb, sa = strata(before), strata(after)
    agree_p, agree_n, agree_tp, agree_ci = judged_precision(JUDGED["agree"])
    body_p, body_n, body_tp, body_ci = judged_precision(JUDGED["body"])
    # The strata precisions are measured on the AFTER sample. Applying them to
    # the BEFORE run assumes the failure modes within a stratum did not change,
    # only the stratum sizes -- which is the point of stratifying, but it is an
    # assumption and it is stated in Limits.
    wb = weighted(sb, PRIOR_ABSTRACT_PRECISION, agree_p, body_p)
    wa = weighted(sa, abs_prec, agree_p, body_p)
    prose_agree = prose_precision(JUDGED["agree"])
    prose_body = prose_precision(JUDGED["body"])
    prose_abs = prose_precision(JUDGED["abstract"], "verdict_v2")
    w_prose = weighted(
        sa,
        prose_abs["precision"] if prose_abs else abs_prec,
        prose_agree["precision"] if prose_agree else agree_p,
        prose_body["precision"] if prose_body else body_p)

    # Extremes over the three strata intervals. Reported because a point
    # estimate quoted five times in a downstream document reads as precision the
    # samples do not support.
    wa_lo = weighted(sa, lo, agree_ci[0], body_ci[0])
    wa_hi = weighted(sa, hi, agree_ci[1], body_ci[1])
    wb_lo = weighted(sb, PRIOR_ABSTRACT_PRECISION, agree_ci[0], body_ci[0])
    wb_hi = weighted(sb, PRIOR_ABSTRACT_PRECISION, agree_ci[1], body_ci[1])

    L = [
        "# Did the #617 co-mention filters work? (#ATLAS-COMENT-REG)", "",
        "Generated by `scripts/comention_regression.py`.", "",
        "**No. Measured on a fresh uniform sample, population-weighted precision",
        f"fell from {100*wb:.1f}% to {100*wa:.1f}%, about "
        f"{100*(wb-wa):.0f} points.**", "",
        "#617 found the layer at 55.5% precision and traced it to one gap:",
        "`usable_alias` exempted every MULTI-WORD form from the specificity test",
        "it applied to single tokens, and 132 of 152 measured false positives were",
        "multi-word. The fix removed the single-token shape proxy entirely and",
        "replaced it with two measured filters. This is the follow-up measurement",
        "that fix asked for.", "",
        "## What changed, per 400 sampled sentences", "",
        "| | before | after |", "|---|---|---|",
        f"| entity mentions | {before['mentions']:,} | {after['mentions']:,} "
        f"({after['mentions']/before['mentions']:.2f}x) |",
        f"| PubTator corroborates | {100*sb['agree']:.1f}% | {100*sa['agree']:.1f}% |",
        f"| abstract-visible disagreement | {100*sb['abstract']:.1f}% | "
        f"**{100*sa['abstract']:.1f}%** |",
        f"| body-only disagreement | {100*sb['body']:.1f}% | {100*sa['body']:.1f}% |", "",
        "The layer now produces a third more mentions per sentence and PubTator",
        "corroborates a smaller share of them. The abstract-visible stratum -- the",
        "only one that can contain false positives, since PubTator read that text --",
        f"nearly tripled.", "",
        "## A correction to how this was judged", "",
        "The first pass judged the SURFACE FORM: does the sentence contain the",
        "string that matched? That is the wrong question, and it flattered the",
        f"layer. It gave {100*tp_v1/n:.1f}% on this sample.", "",
        "The right question is whether the sentence discusses the entity the",
        "IDENTIFIER denotes, because that identifier is what the co-mention pair is",
        "recorded against. Checked against NLM's own descriptor labels, the matched",
        "forms routinely resolve to something else entirely:", "",
        "| matched form | what the identifier actually denotes |", "|---|---|",
        "| `apoptosis` | Malformations of Cortical Development, Group I |",
        "| `node` | Sick Sinus Syndrome |",
        "| `chemotherapy` | Chemotherapy-Related Cognitive Impairment |",
        "| `resistance` | Disease Resistance (host resistance, not drug) |",
        "| `Migration` | Tooth Migration |",
        "| `Aging` | Aging, Premature |",
        "| `artery` | Renal Artery Obstruction |", "",
        "Every one of those was scored correct on the first pass. None is.", "",
        "## Hand-judged precision on that stratum", "",
        f"{n} mentions were read individually "
        f"(`analysis/comention/abstract-visible-judgements.csv`, seeded sample, "
        f"with each identifier's NLM label attached so the judgement is checkable):",
        f"**{tp}/{n} = {100*abs_prec:.1f}%** correct at the identifier level, 95% CI "
        f"[{100*lo:.1f}%, {100*hi:.1f}%].", "",
        f"That is against {100*PRIOR_ABSTRACT_PRECISION:.1f}% before the filter",
        "change -- but the two are NOT known to be methodologically comparable. The",
        "earlier figure was produced under #617 and may have used the same",
        "surface-form criterion this document has just abandoned. Read the",
        "comparison as indicative; the identifier-level number is the sound one.", "",
        "So the stratum did not measurably get cleaner, and it tripled in size.",
        "The weighted total:", "",
        "| stratum | before | after |", "|---|---|---|",
        f"| corroborated | {100*sb['agree']:.1f}% x {100*agree_p:.1f}% | "
        f"{100*sa['agree']:.1f}% x {100*agree_p:.1f}% |",
        f"| abstract-visible | {100*sb['abstract']:.1f}% x "
        f"{100*PRIOR_ABSTRACT_PRECISION:.1f}% | {100*sa['abstract']:.1f}% x "
        f"{100*abs_prec:.1f}% |",
        f"| body-only | {100*sb['body']:.1f}% x {100*body_p:.1f}% | "
        f"{100*sa['body']:.1f}% x {100*body_p:.1f}% |",
        f"| **weighted** | **{100*wb:.1f}%** | **{100*wa:.1f}%** |", "",
        f"Propagating the three strata intervals to their extremes, the after total",
        f"spans **{100*wa_lo:.1f}% to {100*wa_hi:.1f}%** and the before total "
        f"**{100*wb_lo:.1f}% to {100*wb_hi:.1f}%**. The point estimates are what this",
        "document reports elsewhere for readability, and they are not the finding:",
        "the separation between the two runs is, and it survives because the strata",
        "WEIGHTS are exact counts while only the precisions carry sampling error.", "",
        "**All three rows are now hand-judged**, which they were not when this was",
        "first written. The corroborated and body-only precisions were assumed at",
        "92.5% and 30.8%, carried over from #617, and the assumption mattered:", "",
        f"* corroborated measures **{100*agree_p:.1f}%** ({agree_tp}/{agree_n}, 95% CI "
        f"[{100*agree_ci[0]:.1f}%, {100*agree_ci[1]:.1f}%]). The 92.5% assumption",
        "  holds. Its failures are a generic span landing on an over-specific",
        "  descriptor: `hypoxic` to Hypoxia BRAIN, `metabolic` to Metabolic DISEASES.",
        f"* body-only measures **{100*body_p:.1f}%** ({body_tp}/{body_n}, 95% CI "
        f"[{100*body_ci[0]:.1f}%, {100*body_ci[1]:.1f}%]). The 30.8% assumption was",
        "  optimistic by about a third, and this is the stratum carrying the most",
        "  volume. `evaluation` resolves to Phobia, Social; `world` to Cutaneous",
        "  Leishmaniasis; `mPTP`, the mitochondrial permeability transition pore, to",
        "  the neurotoxin MPTP.", "",
        "### Which body-only precision belongs on the BEFORE side", "",
        "The table above uses the after-run 20.0% on both sides. #617 separately",
        "read 39 body-only mentions from the BEFORE run and got 30.8%, so that is a",
        "measurement of that run, not an assumption, and using it raises the before",
        "level:", "",
        "| basis for the before run | before | after | change |", "|---|---|---|---|",
        f"| after-run precisions on both sides (as reported) | {100*wb:.1f}% | "
        f"{100*wa:.1f}% | {100*(wa-wb):+.1f} |",
        f"| each run's own body-only measurement | "
        f"{100*weighted(sb, PRIOR_ABSTRACT_PRECISION, agree_p, 0.308):.1f}% | "
        f"{100*wa:.1f}% | "
        f"{100*(wa-weighted(sb, PRIOR_ABSTRACT_PRECISION, agree_p, 0.308)):+.1f} |", "",
        "The two body-only figures are not separable at these sample sizes -- 8/40",
        "against 12/39 is a two-sided Fisher p of about 0.31 -- so this is a choice",
        "between defensible options rather than a correction. **The reported row is",
        "the conservative one**: it gives the smaller regression. The delta is robust",
        "to the choice (the two runs' body-only shares, 47.5% and 44.2%, nearly",
        "cancel); the LEVEL is not, and quoting 41.6% without saying which basis",
        "produced its counterpart would hide that.", "",
        "## Why retuning the thresholds will not fix it", "",
        "The replacement filters do not separate true matches from false ones on",
        "this stratum. Over the judged sample:", "",
        "| | census support | share of identifier |", "|---|---|---|",
        "| false positives | 77 - 355,572 | 11.5% - 274% |",
        "| true positives | 297 - 100,758 | 20.2% - 68.9% |", "",
        "The distributions overlap almost entirely. No threshold on either filter",
        "recovers the loss without taking the true matches with it, so this is a",
        "design failure rather than a tuning problem.", "",
        "## A concrete bug, and it is biased the wrong way", "",
        "`alias_support` sums a form's mentions across every sense it carries.",
        "`ident_mentions` sums an identifier's across every form. Their ratio is",
        "therefore not a share, and it can exceed 1 -- measured at **274% for**",
        "**`as`**, 133% for `gp`, 104% for `tss`.", "",
        "The direction of the bias is the damaging part. An ambiguous generic word",
        "collects a larger numerator from its other senses, so a MINIMUM-share",
        "filter admits it more readily than a specific name would be admitted. The",
        "filter written to exclude `as` and `treatment` was structurally disposed",
        "to keep them.", "",
        "`atlas_graph.build_index` now records `alias_ident_support`, the count of",
        "a form for the identifier its own majority vote picks, which is the",
        "numerator the ratio always needed.", "",
        "**It is not enough, and that was measured rather than assumed.** Rebuilding",
        "the index and re-scoring the judged forms against the corrected share:", "",
        "* `as` still reads 124.5% (down from 274%), because the co-mention layer can",
        "  redirect a blocklisted form to a curated sense whose identifier is not the",
        "  one the majority vote counted. Numerator and denominator are still not",
        "  guaranteed commensurable.",
        "* More decisively, at the shipped 5% threshold the corrected share cuts",
        "  **nothing**: 0 of 37 judged false positives and 0 of 13 true positives.",
        "  The filter is inert on this stratum whether it is computed correctly or",
        "  not.", "",
        "So the bug is real and worth fixing on its own terms, but fixing it does not",
        "recover the lost precision. That is the same conclusion the overlap above",
        "reaches, arrived at independently: the share filter cannot do this job.", "",
        "## The general lesson", "",
        "The #617 justification was true about the sample it was measured on and",
        "wrong about the population it was applied to. Removing the single-token",
        "test looked free because that day's false positives were multi-word -- but",
        "they were multi-word BECAUSE the multi-word channel was the unfiltered",
        "one. Closing it moved the pressure to the channel that had just been",
        "opened, and the top offenders are now bare English words: `treatment`",
        "(16 of the sampled disagreements), `effects` (15), `as` (13), `advanced`",
        "(10), `function` (7), `other` (6).", "",
        "A filter justified by the error distribution it was measured against",
        "needs re-measuring once it changes that distribution.", "",
        "## What this does NOT change", "",
        "* The multi-word repair was real. `tumor cells`, `overall survival` and",
        "  `et al` no longer reach the matcher, and none appears in the fresh",
        "  sample's disagreements.",
        "* The layer's PURPOSE is unaffected. Body-only matches are still entities",
        "  an abstract-level extractor structurally cannot reach.",
        "* No manuscript number depends on this layer. It feeds",
        "  `atlas-module-support.md`, which argues a zero in the relation column is",
        "  an extraction failure rather than absence of evidence -- a claim that",
        f"  should now be read against {100*wa:.0f}% precision, not 55%.", "",
        "## The headline compares two criteria, so here are all three", "",
        "The bold figure pairs a LENIENT before-abstract (14.6%, measured under #617)",
        "with a STRICT after-abstract. Like for like:", "",
        "| comparison | before | after | change |", "|---|---|---|---|",
        f"| both lenient | {100*wb:.1f}% | "
        f"{100*weighted(sa, v1_prec, agree_p, body_p):.1f}% | "
        f"{100*(weighted(sa, v1_prec, agree_p, body_p)-wb):+.1f} |",
        f"| both strict (before scaled by the same {abs_prec/max(1e-9,v1_prec):.2f} "
        f"factor) | "
        f"{100*weighted(sb, PRIOR_ABSTRACT_PRECISION*abs_prec/max(1e-9,v1_prec), agree_p, body_p):.1f}% | "
        f"{100*wa:.1f}% | "
        f"{100*(wa-weighted(sb, PRIOR_ABSTRACT_PRECISION*abs_prec/max(1e-9,v1_prec), agree_p, body_p)):+.1f} |",
        f"| **as reported (mixed)** | **{100*wb:.1f}%** | **{100*wa:.1f}%** | "
        f"**{100*(wa-wb):+.1f}** |", "",
        "The reported figure is the most pessimistic of the three. All three are",
        "negative, so the direction does not turn on the choice, but the magnitude",
        "is a 7-to-9 point band rather than a single number.", "",
        "## Does the reference-list confound move the headline? (measured: no)", "",
        "A co-mention harvested from a cited title is recorded against the CITING",
        "paper, which is a weaker claim than the layer is usually read as making, and",
        f"it is asymmetric: {struct_tp} of the {agree_tp} corroborated true positives",
        f"sit in structural text against {struct_body} of {body_tp} body-only. That",
        "was disclosed and uncorrected. Restricting every stratum to the paper's own",
        "prose:", "",
        "| stratum | all | prose only | n (prose) |", "|---|---|---|---|",
    ] + [
        f"| {nm} | {100*allp:.1f}% | {100*pr['precision']:.1f}% | {pr['n']} |"
        for nm, allp, pr in (
            ("corroborated", agree_p, prose_agree),
            ("body-only", body_p, prose_body),
            ("abstract-visible", abs_prec, prose_abs))
        if pr
    ] + [
        "",
        f"Weighted, that is **{100*w_prose:.1f}%** against the reported",
        f"{100*wa:.1f}% -- a difference of {abs(100*(w_prose-wa)):.1f} points. **The",
        "confound is real and does not move the headline.**", "",
        "It also does not run the way it was assumed to. Structural text is MORE",
        "precise in the corroborated stratum, not less: a reference title naming two",
        "entities usually names them correctly, because a title is written to be",
        "precise. The effect is a wash rather than an inflation, and the stratum it",
        "helps is the one that needed help least.", "",
        "This holds the strata WEIGHTS fixed and varies only the precisions, so it",
        "answers whether the confound biases the precision estimate. It does not",
        "model a layer that excluded structural text at build time, which would",
        "shrink volume as well.", "",
        "## Held-out validation of the discriminator (#628 criterion 2)", "",
        "Every judged mention above was used to SELECT the authority-name rule, so",
        "none of them can validate it. A second sample was drawn the same way with a",
        f"different seed, disjoint from the first ({ho_n} mentions,",
        "`scripts/comention_judge_prep.py --exclude`), and judged before its numbers",
        "were computed.", "",
        f"**The stratum itself generalises**: held-out precision is {ho_tp}/{ho_n} = "
        f"{100*ho_prec:.1f}% (CI [{100*ho_lo:.1f}%, {100*ho_hi:.1f}%]) against "
        f"{100*abs_prec:.1f}% on the selection sample. The two are the same number",
        "within noise, which is the first thing that had to be true.", "",
        # Derived, not asserted. This was a fixed label claiming the estimate
        # was optimistic -- printed unconditionally, and a test asserted the
        # string to prove the document "reports that kept precision fell". If
        # the held-out sample had come back HIGHER the label would still have
        # printed and the test would still have passed.
        f"**The rule generalises in direction, and its precision was "
        f"{'optimistic' if ho_keep_prec < 0.667 else 'if anything conservative'}:**",
        "",
        "| | selection | held-out |", "|---|---|---|",
        f"| false positives removed | 96% | {100*ho_cut_fp:.0f}% |",
        f"| true positives removed | 60% | {100*ho_cut_tp:.0f}% |",
        f"| precision of what it keeps | 66.7% (4/6) | "
        f"{100*ho_keep_prec:.1f}% ({ho_keep_tp}/{ho_keep_n}) |", "",
        "The direction survives -- it still removes almost every false positive at a",
        "consistent 60% cost in true ones. The kept-precision estimate falls by a",
        "third, which is exactly the selection optimism a held-out sample exists to",
        "expose, and neither figure is well constrained: 4 and "
        f"{ho_keep_n} surviving mentions give intervals so wide they overlap almost",
        "entirely. What is established is that the rule separates; what is not is by",
        "how much.", "",
        "One further observation the held-out sample makes available. Applying the",
        "rule to the MATCHED SPAN rather than the identifier's canonical name is",
        f"stronger here -- it keeps {ho_span_keep} mentions, all correct, removing",
        "every false positive -- but on a sample that small the difference is not",
        "distinguishable from the identifier-level version, and the span was not",
        "recorded when the first sample was judged, so this is a lead rather than a",
        "result.", "",
        "## Does the regression survive if BOTH sides were judged leniently?", "",
        "This is the weakest point in the comparison, so it is worth working out",
        "rather than hedging. The carried-over precisions were measured under #617",
        "and may share the surface-form criterion this document abandoned. If they",
        "do, the `before` total is overstated too, and the regression could in",
        "principle be an artifact of correcting only one side.", "",
        "It is not. Suppose both sides' impure strata are overstated by the same",
        "factor, so their true precisions are `k` times the reported ones for some",
        "`k` in [0,1]. The corroborated stratum is left alone, since PubTator",
        "assigned the same identifier there and that is an identifier-level check",
        "already. Then:", "",
        "| k | before | after | change |", "|---|---|---|---|",
    ]
    for k in (1.0, 0.8, 0.6, 0.4, 0.2, 0.0):
        B = (sb["agree"] * agree_p
             + k * (sb["abstract"] * PRIOR_ABSTRACT_PRECISION
                    + sb["body"] * body_p))
        A = (sa["agree"] * agree_p
             + k * (sa["abstract"] * abs_prec + sa["body"] * body_p))
        L.append(f"| {k:.1f} | {100*B:.1f}% | {100*A:.1f}% | {100*(A-B):+.1f} |")

    L += [
        "", "The change is negative for every `k`, and it gets MORE negative as `k`",
        "falls. Correcting both sides therefore WIDENS the gap rather than closing",
        "it, because the `after` run carries more of its weight in the impure",
        "strata. The regression does not depend on the carried-over numbers being",
        "right; it depends only on their being wrong by a similar factor on both",
        "sides, which is the assumption that would have to hold for the objection",
        "to bite.", "",
        "## A signal that does separate, and what it actually is (#628)", "",
        "The filters above fail because they measure how OFTEN a form appears, and",
        "the offenders are common words. The signal that separates them is whether",
        "the entity's own name looks like a name at all, checked against NLM's",
        "descriptor labels, which are external to this corpus.", "",
        "**Read the next table carefully, because an earlier draft of this section",
        "described it wrongly and the correction changes what it is.** The audit",
        "sample does not record the span that matched. `atlas_comention` writes",
        "`entity_names` as `canon[identifier]` -- the identifier's most frequent",
        "surface form ACROSS THE WHOLE CORPUS -- so the column is constant for a",
        "given identifier and is absent from the sentence in 20 of these 60 rows",
        "(`left ventricular dysfunction` sits against \"She was diagnosed with left",
        "breast cancer\"; the span that fired was presumably `left`).", "",
        "So this is an IDENTIFIER-LEVEL check, not a per-mention one: does this",
        "identifier's dominant corpus name coincide with its authority label? It",
        "cannot distinguish two mentions of the same identifier, and it would be",
        "applied when building the alias map rather than when matching text. That is",
        "a cheaper filter than the one first described, and a blunter one.", "",
        "Compared as a word bag so `squamous cell carcinoma` matches `Carcinoma,",
        "Squamous Cell`:", "",
        "| | true | false |", "|---|---|---|",
        f"| kept | {kept_tp} | {len(kept)-kept_tp} |",
        f"| cut | {tp-kept_tp} | {(n-tp)-(len(kept)-kept_tp)} |", "",
        f"It removes **{100*((n-tp)-(len(kept)-kept_tp))/max(1,n-tp):.0f}% of false",
        f"positives** and {100*(tp-kept_tp)/max(1,tp):.0f}% of true ones, lifting",
        f"precision on what it keeps from {100*abs_prec:.1f}% to "
        f"{100*kept_tp/max(1,len(kept)):.1f}%. Against filters that cut nothing at",
        "all, that is a discriminator rather than a threshold.", "",
        f"It is aggressive. It discards {n-len(kept)} of {n} mentions on this",
        f"stratum ({100*(n-len(kept))/n:.0f}%), and the kept precision carries a 95% CI of",
        f"[{100*klo:.1f}%, {100*khi:.1f}%] on {len(kept)} survivors -- the point",
        "estimate is not the finding, the separation is.", "",
        "**Three reasons it is NOT a recommendation**, beyond the sample size:", "",
        "* It was selected on this sample, so its numbers are optimistic by",
        "  construction. A held-out judged sample is the next requirement.",
        "* **It removes every gene.** NCBI Gene identifiers have no MeSH descriptor,",
        "  so the rule cuts them unconditionally. This sample is 58 of 60 MeSH, so",
        "  the 96% figure says nothing about the gene half of a layer built over",
        "  gene, chemical and disease annotations -- in a repository whose subject is",
        "  GPX4 and ACSL4. A label source covering NCBI Gene is not an optimisation",
        "  here, it is a precondition.",
        "* It is tested only on the abstract-visible stratum. The corroborated",
        f"  stratum is {100*sa['agree']:.0f}% of volume at {100*agree_p:.0f}%",
        "  precision, and the trade would have to be re-judged against that, where",
        "  the filter has much more to lose.", "",
        "The word-bag comparison is also NOT clearly better than plain equality on",
        "this sample: equality fires on 4 rows, the word bag on 6, and the two extra",
        "are one true match (`squamous cell carcinoma`) and one false one (`left",
        "ventricular dysfunction` against `Ventricular Dysfunction, Left`). It adds a",
        "true positive and a false positive, and LOWERS kept precision from 75% to",
        f"{100*kept_tp/max(1,len(kept)):.1f}%. The favourable example was the one",
        "quoted; the unfavourable one is now quoted beside it.", "",
        "## Limits", "",
        f"* {n} judged mentions gives a wide interval "
        f"([{100*lo:.1f}%, {100*hi:.1f}%]); the direction of the net change is",
        "  robust to it, the magnitude is not.",
        f"* **A confound concentrated in the best stratum.** {struct_tp} of the "
        f"{agree_tp} corroborated true positives come from reference-list entries,",
        "  table rows or PMC structural markup rather than the paper's own prose,",
        f"  against {struct_body} of {body_tp} in body-only. A co-mention harvested",
        "  from a cited title is recorded against the CITING paper, which is a",
        "  different claim from \"this paper discusses X and Y\". It is asymmetric,",
        "  it sits in the stratum with the highest precision and a third of the",
        "  volume, and it is not corrected for here.",
        "* Judgement is mine and unblinded. I knew which run each sentence came",
        "  from, which is exactly the bias that would flatter a fix I wrote, and",
        "  the result runs against my own prior change rather than for it.",
        "* The first pass of that judgement was WRONG in a way that flattered the",
        "  layer, and it was caught only because an unrelated experiment surfaced",
        "  the authority labels. The corrected verdicts are committed beside the",
        "  originals so the correction itself can be audited.",
        "* **The strata precisions are measured on the AFTER sample and applied to",
        "  the BEFORE run.** That assumes the failure modes WITHIN a stratum did not",
        "  change, only the stratum sizes -- which is the point of stratifying, but",
        "  it is an assumption and it moves the levels. #617 separately measured the",
        "  before run's body-only stratum at 30.8% on 39 mentions, so that figure is",
        "  a MEASUREMENT of that run rather than an assumption, and substituting the",
        "  after-run 20.0% is the conservative choice. Both bases are tabulated",
        "  above; the delta is robust to the choice and the level is not.",

    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "before": {"mentions": before["mentions"], "strata": sb,
                   "abstract_precision": PRIOR_ABSTRACT_PRECISION, "weighted": wb},
        "after": {"mentions": after["mentions"], "strata": sa,
                  "abstract_precision": abs_prec, "judged_n": n, "judged_tp": tp,
                  "abstract_precision_ci": [lo, hi], "weighted": wa,
                  "range": [wa_lo, wa_hi]},
        "net_change": wa - wb,
        "prose_only": {"weighted": w_prose,
                       "agree": prose_agree, "body": prose_body,
                       "abstract": prose_abs},
        "measured_strata": {
            "agree": {"precision": agree_p, "n": agree_n, "tp": agree_tp, "ci": agree_ci},
            "body_only": {"precision": body_p, "n": body_n, "tp": body_tp, "ci": body_ci},
        },
    }, indent=2) + "\n")
    print(f"weighted precision {100*wb:.1f}% -> {100*wa:.1f}%  ({100*(wa-wb):+.1f} points)")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
