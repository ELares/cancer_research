#!/usr/bin/env python3
"""Two modality descriptors, two different recalls, and a ratio that divided them

WHAT THIS WITHDRAWS
-------------------
`analysis/manuscript-vs-census.md` reports the manuscript's Section 8.2
PDT:SDT ratio as 2.93:1, recomputes it on the census as **4.75:1**, and
concludes the claim "survives, UNDERSTATED by the manuscript". That verdict is
the document's headline and it is repeated in CLAUDE.md.

It divides two MeSH descriptor counts, `Photochemotherapy` and
`Ultrasonic Therapy`. A count ratio is only a ratio of the underlying concepts
if both descriptors capture their concept equally well. Nobody measured that.

Measured here, on the same ferroptosis-indexed census articles with the same
generosity of text rule on both sides:

    Photochemotherapy   recall 80.2%   precision 96.1%
    Ultrasonic Therapy  recall 46.0%   precision 90.6%

The descriptors are comparably PRECISE and differ nearly twofold in RECALL. So
the 4.75:1 is substantially a measurement of how consistently MeSH indexers
apply two descriptors, and only partly of how much literature exists.

Applying one rule to both sides gives **2.89:1** against the manuscript's
2.93:1 -- a 1.4% difference. The honest conclusion is that the census
REPRODUCES the manuscript's ratio, not that the manuscript understated it by
62%.

WHY THE EXISTING SYMMETRY CHECK MISSED IT
-------------------------------------------
`manuscript-vs-census.md` does check for asymmetry, and reports:

    "The gap is 3.3 points, so the over-estimation is symmetric and the
     ratio survives it"

That is a check on PRECISION -- what share of each descriptor's records are
on-modality -- and precision really is symmetric (90.8% vs 87.5%). Recall was
never checked, and that is where the asymmetry lives. The same document's
sensitivity analysis varies four DESCRIPTOR SETS and reports "every variant
exceeds the manuscript's"; every variant is built from descriptors, so all
four inherit the same recall gap. Varying the descriptor set cannot detect a
descriptor-versus-text recall problem.

AND IT INVERTS A CAVEAT CARRIED IN SIX FILES
----------------------------------------------
Five documents state that `Ultrasonic Therapy` is BROADER than sonodynamic
therapy, so the SDT count is an OVER-estimate and every ratio against it is a
LOWER bound. Precision 90.6% says the breadth is real but small -- 3 of 32
records, and reading them they are ultrasound-activated redox, enzyodynamic
and piezocatalytic papers, sonodynamic-adjacent rather than physiotherapy.
Recall 46.0% is far larger and runs the other way. The count is an
UNDER-estimate by roughly a factor of two, and ratios against it are UPPER
bounds.

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
Not that the manuscript is wrong. Its 2.93:1 is reproduced almost exactly by
an independent instrument, which is a stronger result for it than a number
that disagreed. What is withdrawn is the "understated by 62%" verdict built on
an asymmetric comparison.

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
import re
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
    return {
        "subject": SUBJECT,
        "subject_articles": n,
        "arms": stat,
        "ratio_pair": list(RATIO_PAIR),
        "ratio_by_descriptor": by_desc,
        "ratio_by_text": by_text,
        "manuscript_ratio": MANUSCRIPT_RATIO,
        "recall_asymmetry": (max(recalls) / min(recalls)) if len(recalls) > 1
                            and min(recalls) else None,
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
    L += ["| how the two are counted | ratio |", "|---|--:|"]
    L += [f"| by descriptor (what `manuscript-vs-census.md` reports) | "
          f"**{d['ratio_by_descriptor']:.2f}:1** |"]
    L += [f"| by one text rule applied to both | **{d['ratio_by_text']:.2f}:1** |"]
    L += [f"| as the manuscript states it | {d['manuscript_ratio']:.2f}:1 |", ""]

    diff = abs(d["ratio_by_text"] - d["manuscript_ratio"]) / d["manuscript_ratio"]
    over = (d["ratio_by_descriptor"] - d["manuscript_ratio"]) / d["manuscript_ratio"]
    L += [f"So the census does not exceed the manuscript by "
          f"{100*over:.0f}%. Measured symmetrically it **reproduces** the "
          f"manuscript's ratio to within {100*diff:.1f}%.", ""]
    L += ["That is a stronger result for the manuscript than the one being "
          "withdrawn: an independent instrument over the whole indexed cancer "
          "literature lands on the figure it published. What does not survive "
          "is the claim that the manuscript **understated** its own case, "
          "which rested on comparing two descriptors that do not capture "
          "their concepts equally well.", ""]

    L += ["## Why the existing symmetry check did not catch this", ""]
    L += ["`manuscript-vs-census.md` does test for asymmetry and concludes "
          "\"the gap is 3.3 points, so the over-estimation is symmetric\". "
          "That is a test of **precision** -- what share of each descriptor's "
          "records are on-modality -- and precision genuinely is symmetric. "
          "Recall was never measured, and that is where the asymmetry is.", ""]
    L += ["The same document varies four alternative descriptor sets and "
          "reports that every variant exceeds the manuscript's ratio. Every "
          "variant is built from descriptors, so all four inherit the same "
          "recall gap. A sensitivity analysis over descriptor sets cannot "
          "detect a descriptor-versus-text recall problem.", ""]

    L += ["## The caveat this inverts", ""]
    L += [f"Five files state that `Ultrasonic Therapy` is broader than "
          f"sonodynamic therapy, so the count is an OVER-estimate and ratios "
          f"against it are LOWER bounds. The breadth is real and small: "
          f"precision {100*B['precision']:.1f}%, i.e. "
          f"{B['descriptor'] - B['both']} of {B['descriptor']} records. "
          f"Recall is {100*B['recall']:.1f}%, which is far larger and runs the "
          f"other way, so the count is an **under**-estimate and ratios "
          f"against it are **upper** bounds.", ""]
    if B["descriptor_not_text"]:
        L += ["The records carrying the descriptor whose text does not say "
              "sonodynamic -- the whole of the over-count -- are:", ""]
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
    L += ["* Not that the manuscript is wrong. Its ratio is reproduced almost "
          "exactly by an independent instrument.",
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
