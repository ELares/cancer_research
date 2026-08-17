#!/usr/bin/env python3
"""Where does the thesis sit among ALL modalities, not the three it names? (#725)

WHY
---
`atlas_thesis_position.py` measures where this project's thesis sits by
intersecting MeSH `Ferroptosis` with a hand-written list of six modalities. Its
honest headline -- the sonodynamic leg rests on roughly thirty papers in the
entire indexed cancer literature -- is one of the best results in the repo.

But a hand-written list cannot surface a leg nobody thought to look for. The
table enumerates the modalities the thesis already believes in, so a large
ferroptosis intersection with something outside that list is structurally
invisible.

This ranks the intersection against every modality in the committed
`analysis/modality-partitions.json` -- a descriptor universe built for #724 by
five independent panels and reviewed for symmetry -- rather than a list written
by the same person making the claim.

WHAT THE ANSWER IS FOR
----------------------
Two things, and the second is why it is worth doing.

  RANK. Where do the thesis's own legs sit among all modalities that intersect
  ferroptosis? A leg being small is not news -- the repo already says so. A leg
  being small while something unexamined is large is news.

  PRECEDENT. If a modality the thesis does not discuss has a much larger
  ferroptosis literature, that is either strong support for the mechanism or a
  novelty claim that needs qualifying, and the manuscript should say which.

A WITHDRAWN CLAIM THIS REPLACES. An earlier version of #725 asserted that
ferroptosis x radiotherapy is "roughly an order of magnitude larger than the
sonodynamic leg". That reproduced under no descriptor set: the bare
`Radiotherapy` descriptor gives a count BELOW the sonodynamic leg, and the
apparent 10x came from comparing a wide text-stem count against a narrow
descriptor count -- asymmetric rules, the same error corrected in #722. This
script compares descriptor to descriptor throughout, and reports the whole
ranking so no single pair can be cherry-picked.

Usage:
    python scripts/atlas_thesis_rank.py
"""

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = PROJECT_ROOT / "corpus" / "atlas"
PARTITIONS = PROJECT_ROOT / "analysis" / "modality-partitions.json"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-thesis-rank.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-thesis-rank.json"

FERROPTOSIS = "ferroptosis"

# The legs the thesis names, so their rank can be located rather than assumed.
THESIS_LEGS = {
    "sonodynamic therapy": "ultrasonic therapy",
    "photodynamic therapy": "photochemotherapy",
    "focused ultrasound": "high-intensity focused ultrasound ablation",
    "drug resistance": "drug resistance, neoplasm",
}

# Modalities with no usable MeSH descriptor. Named rather than silently absent,
# which is the treatment `mechanism_recall.py` already gives unmeasurable
# mechanisms -- reporting them as a zero would be a different claim.
NOT_MEASURABLE = [
    "tumour-treating fields (no descriptor; `Electric Stimulation Therapy` is broader)",
    "cold atmospheric plasma (no descriptor)",
    "sonodynamic therapy specifically. `Ultrasonic Therapy` IS broader, but "
    "that is the smaller effect: measured, its precision is 90.6% (3 records "
    "of 32) while its recall is 46.0% (34 ferroptosis-SDT papers carry no "
    "such descriptor), so the count is an UNDER-estimate by roughly twofold "
    "and every ratio against it is an UPPER bound. An earlier version of this "
    "list stated the opposite; see analysis/atlas-descriptor-recall.md",
]


def modality_universe() -> set:
    """Every modality descriptor the committed partitions name."""
    d = json.loads(PARTITIONS.read_text())
    out = set()
    for spec in d.values():
        for side in ("pharmacological", "physical"):
            for x in spec.get(side, []):
                out.add(x.strip().lower())
    if not out:
        raise SystemExit(f"no descriptors in {PARTITIONS}")
    return out


def scan(universe: set) -> dict:
    inter = Counter()
    ferro_total = 0
    census = 0
    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                census += 1
                mesh = {m.lower() for m in (r.get("mesh") or [])}
                if FERROPTOSIS not in mesh:
                    continue
                ferro_total += 1
                for m in mesh & universe:
                    inter[m] += 1
    return {"census": census, "ferroptosis_total": ferro_total,
            "universe_size": len(universe),
            "intersections": [[k, v] for k, v in inter.most_common()]}


def render(d: dict) -> str:
    rows = d["intersections"]
    rank = {k: i + 1 for i, (k, _v) in enumerate(rows)}
    counts = dict(rows)
    L = ["# Where the thesis sits among all modalities", ""]
    L += ["*Generated by `scripts/atlas_thesis_rank.py`. The modality universe "
          "is the committed `analysis/modality-partitions.json`, built for #724 "
          "by five independent panels -- not a list written by the person "
          "making the claim.*", ""]

    L += [f"**{d['ferroptosis_total']:,} census articles carry MeSH "
          f"`Ferroptosis`.** Their intersection with each of "
          f"{d['universe_size']:,} modality descriptors, ranked:", ""]
    L += ["| rank | modality descriptor | ferroptosis articles | share of ferroptosis |",
          "|--:|---|--:|--:|"]
    for i, (k, v) in enumerate(rows[:20], 1):
        L.append(f"| {i} | {k} | {v:,} | "
                 f"{100*v/max(d['ferroptosis_total'],1):.2f}% |")
    L += [""]

    L += ["## Where the thesis's own legs land", ""]
    L += ["| leg | descriptor | count | rank of "
          f"{len(rows)} |", "|---|---|--:|--:|"]
    for leg, desc in THESIS_LEGS.items():
        c = counts.get(desc, 0)
        r = rank.get(desc)
        L.append(f"| {leg} | {desc} | {c:,} | "
                 f"{r if r else 'not in the universe'} |")
    L += [""]

    top = rows[0] if rows else ("none", 0)
    sdt = counts.get(THESIS_LEGS["sonodynamic therapy"], 0)
    if sdt:
        L += [f"The largest ferroptosis-modality intersection is **{top[0]}** at "
              f"{top[1]:,} articles, which is **{top[1]/sdt:.1f}x** the "
              f"sonodynamic leg the simulation half rests on.", ""]

    L += ["## Not measurable this way, and named rather than shown as zero", ""]
    for n in NOT_MEASURABLE:
        L.append(f"* {n}")
    L += [""]

    L += ["## What this does and does not establish", ""]
    L += ["* A large intersection is ATTENTION, not endorsement. It says a "
          "literature exists in which both concepts are indexed together, not "
          "that the combination works.",
          "* It does not replace the thesis-position analysis. That one asks "
          "whether the thesis's legs are thin; this one asks whether anything "
          "outside them is thick, which a hand-written list cannot answer.",
          "* Descriptor-to-descriptor throughout. An earlier version of this "
          "question compared a wide text-stem count against a narrow descriptor "
          "count and produced a ratio that reproduced under no definition.",
          "* The whole ranking is published so no single pair can be picked out "
          "to support a conclusion chosen first.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan(modality_universe())
        if d["ferroptosis_total"] == 0:
            raise SystemExit(
                "no ferroptosis articles found, which is not a finding -- it is "
                "what a descriptor-case mismatch looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for k, v in d["intersections"][:8]:
        print(f"  {v:>6,}  {k}")


if __name__ == "__main__":
    main()
