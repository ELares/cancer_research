#!/usr/bin/env python3
"""Two modality descriptors, two different recalls, and a ratio that divided them

WHAT THIS WITHDRAWS
-------------------
`analysis/manuscript-vs-census.md` recomputes the manuscript's Section 8.2
PDT:SDT ratio on the census, finds a larger figure, and concludes the claim
"survives, UNDERSTATED by the manuscript". That verdict is the document's
headline and it is repeated in CLAUDE.md. Every figure below is rendered from
the measurement; none is written into this docstring.

It divides two MeSH descriptor counts, `Photochemotherapy` and
`Ultrasonic Therapy`. A count ratio is only a ratio of the underlying concepts
if both descriptors capture their concept equally well. Nobody measured that.

Measured here on the same ferroptosis-indexed census articles, with one text
rule applied identically to both arms, the two descriptors turn out comparably
PRECISE and markedly different in RECALL. So the descriptor ratio is
substantially a measurement of how consistently MeSH indexers apply two terms,
and only partly of how much literature exists.

WHAT THE DATA SUPPORTS, AND WHAT AN EARLIER DRAFT OVER-CLAIMED
----------------------------------------------------------------
That earlier draft said the symmetric ratio "reproduces the manuscript's ratio
to within a fraction of a percent". WITHDRAWN. On counts this size the symmetric ratio carries a
wide interval, the difference from the manuscript's figure sits well inside
sampling noise, and a reviewer's bootstrap put the probability of exceeding it
near one half. Quoting a coin flip to that precision is the same
error class this page exists to correct.

Two claims ARE well supported, and the page is now built on them:

  the two recalls have disjoint intervals, so the descriptors genuinely
  differ; and

  the PAIRED inflation -- descriptor ratio over symmetric ratio, computed on
  the same articles so shared noise cancels -- has an interval excluding 1.

Those withdraw the understatement verdict, because the evidence for it
disappears under a symmetric rule. They do NOT establish that the census
agrees with the manuscript, and this page does not say so. Every figure is
interpolated at render time; none is written here, because a docstring cannot
be kept fresh.

WHY THE EXISTING SYMMETRY CHECK MISSED IT
-------------------------------------------
`manuscript-vs-census.md` does check for asymmetry, and reports:

    "the over-estimation is symmetric and the ratio survives it"

That is a check on PRECISION -- what share of each descriptor's records are
on-modality -- and precision really is symmetric. RECALL was never checked,
and that is where the asymmetry lives. The same document's sensitivity
analysis varies several DESCRIPTOR SETS and reports that every variant exceeds
the manuscript's; every variant is built from descriptors, so all of them
inherit the same recall gap. Varying the descriptor set cannot detect a
descriptor-versus-text recall problem.

AND IT INVERTS A CAVEAT CARRIED ACROSS THE REPO
------------------------------------------------
Several documents state that `Ultrasonic Therapy` is BROADER than sonodynamic
therapy, so the SDT count is an OVER-estimate and every ratio against it is a
LOWER bound. The breadth is real but small; the recall shortfall is larger and
runs the other way, so the count is an UNDER-estimate and ratios against it
are UPPER bounds. Both quantities are rendered from the measurement rather
than stated here -- an earlier draft hand-wrote the number of affected
files, and got it wrong.

SCOPE OF THE WORD "RECALL". The denominator is what each record's own title
and abstract say, not ground truth. A record can carry a descriptor correctly
while its abstract never uses the word, so this measures TEXT AGREEMENT and is
labelled as such wherever it is used to compare the two arms.

WHY THE TWO DESCRIPTORS DIFFER, WHICH IS ITS OWN FINDING
----------------------------------------------------------
Photodynamic therapy has a dedicated MeSH descriptor. Sonodynamic therapy has
none, and is indexed under a broader physical-modality term when it is indexed
at all. That is the mechanism `analysis/atlas-untagged-partner.md` describes
for the frozen corpus -- a modality with no name of its own is recorded under
someone else's -- appearing again one layer up, in MeSH rather than in this
project's taxonomy.

WHAT THIS DOES NOT CLAIM
-------------------------
Not that the manuscript is wrong, and equally not that it is confirmed. A
symmetric measurement cannot DISTINGUISH the census from the manuscript's
figure, which is different from and weaker than agreement. What is withdrawn
is the "understated" verdict, because the evidence for it was an asymmetric
comparison. An earlier draft of this very section then said the manuscript's
ratio was "reproduced almost exactly", which is the withdrawn over-precision
restated a few lines below its own withdrawal.

Nor that the text rule is ground truth. It has its own errors in both
directions; the point is that it is applied IDENTICALLY to both arms, which
the descriptor comparison is not.

Usage:
    python scripts/atlas_descriptor_recall.py
    python scripts/atlas_descriptor_recall.py --render-only
"""

import argparse
import gzip
import json
import math
import re
import statistics
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = PROJECT_ROOT / "corpus" / "atlas"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-descriptor-recall.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-descriptor-recall.json"

SUBJECT = "ferroptosis"

# One shape of text rule for every arm: the full word, its hyphenation, and the
# conventional acronym. Asymmetry between arms here would reproduce the defect
# this analysis exists to correct, so the shape is held fixed and only the stem
# changes.
ARMS = {
    "PDT": {
        "label": "photodynamic therapy",
        "descriptors": {"photochemotherapy"},
        "text": r"photodynamic|photo-dynamic|\bPDT\b",
    },
    "SDT": {
        "label": "sonodynamic therapy",
        "descriptors": {"ultrasonic therapy"},
        "text": r"sonodynamic|sono-dynamic|\bSDT\b",
    },
}
# The pair the manuscript's Section 8.2 states, and its stated value.
RATIO_PAIR = ("PDT", "SDT")
MANUSCRIPT_RATIO = 2.93
# The interval level, interpolated wherever it is printed. A "95%" typed into
# prose is still a number a sentence can outlive, and this file's own subject
# is exactly that failure.
CONF = 0.95
# DERIVED from CONF, never typed beside it. Hardcoding Z let a
# mutation print a 68% interval under a "95% CI" label, because the
# label was interpolated and the arithmetic was not.
Z = statistics.NormalDist().inv_cdf(1 - (1 - CONF) / 2)


def scan() -> dict:
    pats = {k: re.compile(v["text"], re.I) for k, v in ARMS.items()}
    # the descriptor name travels WITH the counts, so a consumer quoting this
    # artifact names the right descriptor instead of falling back to the arm
    # key and printing "`PDT` recalls 80.2% of PDT papers"
    stat = {k: {"text": 0, "descriptor": 0, "both": 0,
                "descriptors": sorted(v["descriptors"]),
                "label": v["label"],
                "descriptor_not_text": []} for k, v in ARMS.items()}
    n = 0
    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                mesh = {m.lower() for m in (r.get("mesh") or [])}
                if SUBJECT not in mesh:
                    continue
                n += 1
                blob = (r.get("title") or "") + " " + (r.get("abstract") or "")
                for k, arm in ARMS.items():
                    t = bool(pats[k].search(blob))
                    d = bool(mesh & arm["descriptors"])
                    stat[k]["text"] += t
                    stat[k]["descriptor"] += d
                    if t and d:
                        stat[k]["both"] += 1
                    elif d and not t and len(stat[k]["descriptor_not_text"]) < 12:
                        stat[k]["descriptor_not_text"].append(
                            {"pmid": r.get("pmid"),
                             "title": (r.get("title") or "")[:110]})

    for k, s in stat.items():
        s["recall"] = s["both"] / s["text"] if s["text"] else None
        s["precision"] = s["both"] / s["descriptor"] if s["descriptor"] else None

    a, b = RATIO_PAIR
    by_desc = (stat[a]["descriptor"] / stat[b]["descriptor"]
               if stat[b]["descriptor"] else None)
    by_text = stat[a]["text"] / stat[b]["text"] if stat[b]["text"] else None
    recalls = [s["recall"] for s in stat.values() if s["recall"] is not None]

    # UNCERTAINTY, because the first version of this analysis quoted the
    # symmetric ratio as reproducing the manuscript to a fraction of a percent
    # of 182 and 63. A reviewer's bootstrap put the 95% interval at roughly
    # [2.2, 3.9] and P(ratio > manuscript) at 0.46 -- the verdict was a coin
    # flip quoted to three significant figures. Closed-form Poisson
    # log-ratio intervals are used so this stays deterministic and offline.
    def _logratio_ci(x, y):
        """Interval for x/y treating both as Poisson counts."""
        if not x or not y:
            return None
        se = math.sqrt(1 / x + 1 / y)
        pt = math.log(x / y)
        return [math.exp(pt - Z * se), math.exp(pt + Z * se)]

    def _wilson(k, n):
        if not n:
            return None
        p, z = k / n, Z
        d = 1 + z * z / n
        c = p + z * z / (2 * n)
        m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
        return [(c - m) / d, (c + m) / d]

    for k, s in stat.items():
        s["recall_ci"] = _wilson(s["both"], s["text"])
        s["precision_ci"] = _wilson(s["both"], s["descriptor"])

    text_ci = _logratio_ci(stat[a]["text"], stat[b]["text"])
    desc_ci = _logratio_ci(stat[a]["descriptor"], stat[b]["descriptor"])
    # the comparison that is actually decisive: how much the descriptor route
    # inflates the ratio relative to the symmetric one
    inflation = by_desc / by_text if (by_desc and by_text) else None
    # PAIRED, not independent. An earlier version summed 1/x over all four
    # counts, treating each arm's descriptor and text totals as independent
    # Poissons. They are not -- they SHARE `both`, so
    #     Var(log(D/T)) = 1/D + 1/T - 2*both/(D*T)
    # The independent form gave an SE 1.73x too large. It is CONSERVATIVE, so
    # nothing published on it was false, but it put the lower bound at 1.02
    # rather than 1.25 and, under two defensible text rules, flipped the
    # significance flag to False and DISCARDED the one claim this analysis
    # positively supports. An over-wide interval is not the safe direction
    # when a flag keys off it. Verified against a 4,000-resample
    # article-level bootstrap by a reviewer.
    var = 0.0
    for k in (a, b):
        D, T, B = stat[k]["descriptor"], stat[k]["text"], stat[k]["both"]
        if not (D and T):
            var = None
            break
        var += 1 / D + 1 / T - 2 * B / (D * T)
    infl_se = math.sqrt(var) if var and var > 0 else None
    infl_ci = ([math.exp(math.log(inflation) - Z * infl_se),
                math.exp(math.log(inflation) + Z * infl_se)]
               if (inflation and infl_se) else None)
    # and whether the recall gap itself excludes parity
    ra = (max(recalls) / min(recalls)) if len(recalls) > 1 and min(recalls) else None
    return {
        "subject": SUBJECT,
        "subject_articles": n,
        "arms": stat,
        "ratio_pair": list(RATIO_PAIR),
        "ratio_by_descriptor": by_desc,
        "ratio_by_descriptor_ci": desc_ci,
        "ratio_by_text": by_text,
        "ratio_by_text_ci": text_ci,
        "manuscript_ratio": MANUSCRIPT_RATIO,
        "recall_asymmetry": ra,
        # the decisive quantity: the DIFFERENCE between the two routes, which
        # is a paired comparison over the same articles and is far better
        # determined than either ratio on its own
        "descriptor_inflation": inflation,
        "descriptor_inflation_ci": infl_ci,
        "symmetric_ratio_covers_manuscript": bool(
            text_ci and text_ci[0] <= MANUSCRIPT_RATIO <= text_ci[1]),
        "descriptor_route_significantly_inflates": bool(
            infl_ci and infl_ci[0] > 1.0),
    }


def render(d: dict) -> str:
    a, b = d["ratio_pair"]
    A, B = d["arms"][a], d["arms"][b]
    L = ["# Two modality descriptors, two different recalls", ""]
    L += ["*Generated by `scripts/atlas_descriptor_recall.py`. Every figure is "
          "recomputed.*", ""]

    L += [f"Over the {d['subject_articles']:,} census articles indexed "
          f"`{d['subject'].title()}`, measured with one text rule applied "
          f"identically to both arms:", ""]
    L += ["| arm | descriptor | text says so | descriptor | both | recall | precision |",
          "|---|---|--:|--:|--:|--:|--:|"]
    for k in (a, b):
        s, arm = d["arms"][k], ARMS[k]
        L.append(
            f"| {k} | `{sorted(arm['descriptors'])[0].title()}` | "
            f"{s['text']:,} | {s['descriptor']:,} | {s['both']:,} | "
            f"**{100*s['recall']:.1f}%** | {100*s['precision']:.1f}% |")
    L += [""]

    if d["recall_asymmetry"]:
        L += [f"The two descriptors are comparably **precise** and differ by a "
              f"factor of **{d['recall_asymmetry']:.2f}** in **recall**.", ""]

    L += ["## What that does to the ratio the manuscript states", ""]
    def ci(v):
        return f"[{v[0]:.2f}, {v[1]:.2f}]" if v else "n/a"

    L += [f"| how the two are counted | ratio | {100*CONF:g}% CI |",
          "|---|--:|--:|"]
    L += [f"| by descriptor (what `manuscript-vs-census.md` reports) | "
          f"**{d['ratio_by_descriptor']:.2f}:1** | "
          f"{ci(d['ratio_by_descriptor_ci'])} |"]
    L += [f"| by one text rule applied to both | **{d['ratio_by_text']:.2f}:1** "
          f"| {ci(d['ratio_by_text_ci'])} |"]
    L += [f"| as the manuscript states it | {d['manuscript_ratio']:.2f}:1 | |",
          ""]

    L += ["**Both ratios are loosely determined and their intervals overlap, "
          "so neither point estimate should be read alone.** An earlier draft "
          "of this page said the symmetric ratio \"reproduces the manuscript's "
          "ratio to within a fraction of a percent\" -- a precision the counts "
          "come nowhere near supporting, on a difference well inside the "
          "sampling noise. That sentence is withdrawn.", ""]

    L += ["### What the data does support", ""]
    if d.get("symmetric_ratio_covers_manuscript"):
        L += [f"The manuscript's {d['manuscript_ratio']:.2f}:1 sits **inside** "
              f"the symmetric interval {ci(d['ratio_by_text_ci'])}, so a "
              f"symmetric measurement cannot distinguish the census from the "
              f"manuscript. The evidence for UNDERSTATEMENT disappears -- "
              f"which is what withdraws that verdict. It does not establish "
              f"agreement, and this page does not claim it.", ""]
    if d.get("descriptor_route_significantly_inflates"):
        L += [f"The decisive quantity is the **paired** comparison of the two "
              f"routes over the same articles: the descriptor route inflates "
              f"the ratio by **{d['descriptor_inflation']:.2f}x**, "
              f"{100*CONF:g}% CI "
              f"{ci(d['descriptor_inflation_ci'])}, which excludes 1. That is "
              f"far better determined than either ratio alone, because the "
              f"two routes share their articles and much of their noise "
              f"cancels.", ""]
    L += ["So the robust claims are that the two descriptors differ in recall "
          "and that the descriptor route inflates this particular ratio. The "
          "size of the true ratio is not settled here.", ""]

    L += ["## Why the existing symmetry check did not catch this", ""]
    L += ["`manuscript-vs-census.md` does test for asymmetry and concludes "
          "that \"the over-estimation is symmetric\" from a small gap in the "
          "share of each descriptor's records that are on-modality. That is a "
          "test of **precision**, and precision genuinely is symmetric. "
          "Recall was never measured, and that is where the asymmetry is.", ""]
    L += ["The same document varies several alternative descriptor sets and "
          "reports that every variant exceeds the manuscript's ratio. Every "
          "variant is built from descriptors, so all of them inherit the same "
          "recall gap. A sensitivity analysis over descriptor sets cannot "
          "detect a descriptor-versus-text recall problem.", ""]

    L += ["## The caveat this inverts", ""]
    L += [f"Several files state that `Ultrasonic Therapy` is broader than "
          f"sonodynamic therapy, so the count is an OVER-estimate and ratios "
          f"against it are LOWER bounds. The breadth is real and small: "
          f"precision {100*B['precision']:.1f}%, i.e. "
          f"{B['descriptor'] - B['both']} of {B['descriptor']} records. "
          f"Recall is {100*B['recall']:.1f}%, which is far larger and runs the "
          f"other way, so the count is an **under**-estimate and ratios "
          f"against it are **upper** bounds.", ""]
    if B["descriptor_not_text"]:
        L += ["The records carrying the descriptor whose text does not use "
              "this rule's words are listed below. They are an upper bound on "
              "descriptor breadth, NOT a measurement of it: at least one is "
              "matched by the wider vocabulary in "
              "`manuscript_vs_census.py` (it says *sonosensitizer*), so it is "
              "a miss of this rule rather than a wrong descriptor. An earlier "
              "draft called these \"the whole of the over-count\", which "
              "conflated the two.", ""]
        for r in B["descriptor_not_text"]:
            L.append(f"* {r['pmid']} — {r['title']}")
        L += [""]

    L += ["## Why the two descriptors differ, which is its own finding", ""]
    L += ["Photodynamic therapy has a dedicated MeSH descriptor. Sonodynamic "
          "therapy has none, and is indexed under a broader physical-modality "
          "term when it is indexed at all. That is the mechanism "
          "`analysis/atlas-untagged-partner.md` measured inside this project's "
          "own taxonomy -- a modality with no name of its own is recorded "
          "under someone else's -- appearing again one layer up, in MeSH.", ""]

    L += ["## What this does not claim", ""]
    L += ["* Not that the manuscript is wrong, and equally **not that it is "
          "confirmed**. A symmetric measurement cannot tell the census apart "
          "from the manuscript's figure, which is weaker than agreement. An "
          "earlier draft of this bullet said the ratio was \"reproduced "
          "almost exactly\" -- the withdrawn over-precision, restated inside "
          "the section that exists to disclaim it, and a reviewer found it "
          "there.",
          "* Not that the text rule is ground truth. It has errors in both "
          "directions. The point is that it is applied IDENTICALLY to both "
          "arms, which the descriptor comparison is not.",
          "* Not that either descriptor is a bad descriptor. They are "
          "comparably precise; they are applied at different rates.",
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
        if not d["subject_articles"]:
            raise SystemExit(
                f"no article carries `{SUBJECT}`, which is not a finding -- it "
                "is what a descriptor-case mismatch looks like.")
        for k, s in d["arms"].items():
            if not s["text"] or not s["descriptor"]:
                raise SystemExit(
                    f"{k} matched nothing on one axis, so no recall or "
                    "precision can be computed; that is a broken rule, not a "
                    "finding.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for k, s in d["arms"].items():
        print(f"  {k:4s} recall {100*s['recall']:.1f}%  "
              f"precision {100*s['precision']:.1f}%")
    print(f"  ratio by descriptor {d['ratio_by_descriptor']:.2f}:1  "
          f"by text {d['ratio_by_text']:.2f}:1  "
          f"manuscript {d['manuscript_ratio']:.2f}:1")


if __name__ == "__main__":
    main()
