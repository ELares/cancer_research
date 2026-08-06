#!/usr/bin/env python3
"""Atlas: find claims the cancer literature contradicts itself on (#ATLAS).

WHY
---
This project's stated ethos is to be honest about what is contested -- the
manuscript flags its own hypoxia leg as "genuinely contested" after a manual
five-lens review of a handful of papers. That review does not scale, and it only
covers what the author thought to check.

The relation graph makes contestedness measurable across the whole corpus. Two
shapes are detectable:

  DIRECTION  the same entity pair is asserted as `positive_correlate` in some
             papers and `negative_correlate` in others;
  VALENCE    a chemical is asserted to `treat` a disease and to `cause` it.

Neither is automatically an error. A gene can be protective in one tissue and
harmful in another, and a chemotherapeutic genuinely does both treat cancer and
cause secondary malignancy. That is the point: these are the pairs where a
single citation is least trustworthy, and where the project's "cite one paper
you read" pattern is most likely to pick a side by accident.

SCORING
-------
Ranked by the WEAKER side's count, not the total, so a 50-vs-1 split (a settled
claim with one outlier) ranks below a 12-vs-9 split (genuinely divided). The
balance ratio is reported so a reader can see which they are looking at.

LIMITS
------
`positive_correlate`/`negative_correlate` are extraction outputs, not curated
judgements; the extractor scores ~79.6 F1 on BioRED, so some conflicts are
extraction error rather than scientific disagreement.

Two of those failure modes have since been MEASURED
(`scripts/atlas_contradiction_quality.py`), with opposite answers:

  * A single paper extracted as asserting both directions would be inconsistency
    rather than disagreement. It happens to 1 paper in 115,024. This mode is
    effectively absent, and the conflicts really are between studies.
  * Merging two entities under one identifier merges two literatures, which will
    disagree. Pairs involving an identifier measured as a SENSE COLLISION are
    **1.45x** more likely to be flagged contradictory (95% CI 1.37-1.53). That
    survives stratifying by assertion count, so it is not the popularity artifact
    it could have been, and it rises with assertion count -- the direction
    conflation predicts. Check any conflict involving a blocklisted symbol for
    conflation before reading it as a scientific dispute. Directionality is also not
preserved: the graph does not record which entity is subject. And no context is
attached -- a relation true in one cell line and false in another appears here
as a contradiction. Treat the output as a QUEUE FOR READING, not a verdict.

Usage:
    python scripts/atlas_contradictions.py
    python scripts/atlas_contradictions.py --focus GPX4 SLC7A11 ACSL4 --top 40
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_graph import load_index, resolve, resolve_reason  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-contradictions.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-contradictions.json"

MIN_WEAK = 3      # the weaker side needs this much support to count as a dispute
MIN_TOTAL = 8     # and the pair needs this much total attention


def scan(idx: dict, focus_ids=None):
    direction, valence = [], []
    for key, preds in idx["edges"].items():
        if focus_ids and not (key[0] in focus_ids or key[1] in focus_ids):
            continue
        pos, neg = preds.get("positive_correlate", 0), preds.get("negative_correlate", 0)
        treat, cause = preds.get("treat", 0), preds.get("cause", 0)
        total = sum(preds.values())
        names = (idx["canon"].get(key[0], key[0]), idx["canon"].get(key[1], key[1]))

        if pos >= MIN_WEAK and neg >= MIN_WEAK and total >= MIN_TOTAL:
            weak, strong = min(pos, neg), max(pos, neg)
            direction.append(dict(a=key[0], b=key[1], a_name=names[0], b_name=names[1],
                                  positive=pos, negative=neg, total=total,
                                  weaker=weak, balance=weak / strong,
                                  pmids=idx["pmids"].get(key, [])[:6]))
        if treat >= MIN_WEAK and cause >= MIN_WEAK and total >= MIN_TOTAL:
            weak, strong = min(treat, cause), max(treat, cause)
            valence.append(dict(a=key[0], b=key[1], a_name=names[0], b_name=names[1],
                                treat=treat, cause=cause, total=total,
                                weaker=weak, balance=weak / strong,
                                pmids=idx["pmids"].get(key, [])[:6]))
    direction.sort(key=lambda r: (-r["weaker"], -r["balance"]))
    valence.sort(key=lambda r: (-r["weaker"], -r["balance"]))
    return direction, valence


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--focus", nargs="*", default=None,
                    help="restrict to pairs involving these entities")
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()

    idx = load_index(atlas_root())
    focus_ids = None
    if args.focus:
        focus_ids = set()
        for f in args.focus:
            r = resolve(idx, f)
            if r:
                focus_ids.add(r)
            else:
                # a blocked sense collision is not the same failure as an
                # unknown symbol, and the caller needs to know which it was
                print(f"  warning: could not resolve {f!r}. "
                      f"{resolve_reason(idx, f)}", file=sys.stderr)

    direction, valence = scan(idx, focus_ids)
    # The TOTALS, not just the top 500. They were previously printed to stdout
    # and nowhere else, so the figures quoted downstream had been copied from a
    # terminal: uncheckable, undetectable when they drifted, and by the time
    # anyone re-ran this they were wrong by about 50% (4,667 -> 7,068 direction,
    # 6,764 -> 9,094 valence) because the disambiguation corrections had grown
    # and merged more entities. The truncated lists are byte-identical across
    # that change -- the top conflicts are stable while the tail grows -- so
    # nothing in the artifact could have revealed it.
    RAW.write_text(json.dumps({
        "direction_conflicts": len(direction),
        "valence_conflicts": len(valence),
        "listed": 500,
        "direction": direction[:500], "valence": valence[:500],
    }, indent=1), encoding="utf-8")

    L = [
        "# Where the cancer literature contradicts itself (#ATLAS)", "",
        "Generated by `scripts/atlas_contradictions.py` over the atlas relation graph.", "",
        "## What this detects", "",
        "**Direction conflicts** -- the same entity pair asserted `positive_correlate` in",
        "some papers and `negative_correlate` in others. **Valence conflicts** -- a",
        "chemical asserted both to `treat` a disease and to `cause` it.", "",
        "Neither is automatically an error. A gene can be protective in one tissue and",
        "harmful in another; a chemotherapeutic genuinely does both treat cancer and cause",
        "secondary malignancy. The point is that these are the pairs where citing a single",
        "paper is least safe -- exactly the pattern the simulation modules use.", "",
        f"**{len(direction):,} direction conflicts** and **{len(valence):,} valence",
        f"conflicts** were found; the tables below list the top 500 of each.", "",
        "Ranked by the WEAKER side's count, so a 50-vs-1 split (settled, one outlier) ranks",
        "below a 12-vs-9 split (genuinely divided). `balance` is weaker/stronger.", "",
        "## Limits -- this is a reading queue, not a verdict", "",
        "The predicates are extraction outputs, not curated judgements, and the extractor",
        "scores ~79.6 F1 on BioRED, so some conflicts are extraction error. The graph does",
        "not record which entity is the subject, so direction of effect is lost. No context",
        "is attached, so a relation true in one cell line and false in another appears here",
        "as a contradiction.", "",
        f"Thresholds: weaker side >= {MIN_WEAK}, pair total >= {MIN_TOTAL}.", "",
        f"## Direction conflicts ({len(direction)} found)", "",
        "| pair | positive | negative | balance | total | example PMIDs |",
        "|---|---|---|---|---|---|",
    ]
    for r in direction[:args.top]:
        L.append(f"| {r['a_name']} — {r['b_name']} | {r['positive']} | {r['negative']} | "
                 f"{r['balance']:.2f} | {r['total']} | {', '.join(r['pmids'][:4])} |")

    L += ["", f"## Valence conflicts: treats and causes ({len(valence)} found)", "",
          "| pair | treat | cause | balance | total | example PMIDs |",
          "|---|---|---|---|---|---|"]
    for r in valence[:args.top]:
        L.append(f"| {r['a_name']} — {r['b_name']} | {r['treat']} | {r['cause']} | "
                 f"{r['balance']:.2f} | {r['total']} | {', '.join(r['pmids'][:4])} |")

    L += ["", "## How to use this", "",
          "Before a simulation module encodes a mechanism on the strength of one paper,",
          "check whether its entity pair appears here. If it does, the module is picking a",
          "side in a live disagreement, and the module docs should say so rather than",
          "citing the single paper that happened to be read.", ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"direction conflicts: {len(direction)}   valence conflicts: {len(valence)}")
    for r in direction[:10]:
        print(f"  {r['a_name'][:26]:<28}{r['b_name'][:26]:<28} +{r['positive']:<5}-{r['negative']:<5}"
              f"bal={r['balance']:.2f}")


if __name__ == "__main__":
    main()
