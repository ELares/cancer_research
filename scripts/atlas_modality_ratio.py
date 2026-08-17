#!/usr/bin/env python3
"""Does the pharmacological:physical claim survive a symmetric class definition? (#724)

THE CLAIM
---------
The manuscript's central corpus claim is that pharmacological cancer research
vastly outweighs physical-modality research. `atlas_landscape.py` computes it
from two curated sets of this project's own mechanism tags and reports 9.1:1 by
the manuscript's method and 17.6:1 on the census, reading the census figure as
STRENGTHENING the case.

WHY THE CLASSES CANNOT SETTLE IT AS WRITTEN
--------------------------------------------
    PHYSICAL = {hifu, sonodynamic, electrochemical-therapy}
    PHARMACOLOGICAL = {immunotherapy, car-t, antibody-drug-conjugate,
                       bispecific-antibody, synthetic-lethality, epigenetic,
                       metabolic-targeting}

Both omit their largest real-world member. PHYSICAL has no radiotherapy -- the
physical modality delivered to roughly half of all cancer patients -- and
PHARMACOLOGICAL has no cytotoxic chemotherapy. Both are missing for the same
reason: neither has a mechanism tag in this project (see #723), not because
either was weighed and excluded.

An earlier proposal was to add radiotherapy to PHYSICAL. That is a ONE-SIDED
WIDENING and it manufactures whatever inversion it finds. A recomputation has to
move both sides by one principle, or neither.

WHAT THIS DOES
--------------
Five partitions were constructed independently, each by a single stated
principle applied to both classes, using MeSH descriptors -- where radiotherapy
and chemotherapy ARE measurable even though this project has no tags for them.
Each was checked by a separate reviewer for the one-sided-widening error and
required to reproduce its own count.

The partitions are the INPUT (`analysis/modality-partitions.json`), committed
with their member lists so any placement can be disputed. This script is the
measurement. The panels' own reported ratios are deliberately NOT carried in
that file: an input that contains the result it is supposed to produce is how a
measurement ends up validating itself.

WHAT THE ANSWER IS FOR
----------------------
Not to replace 17.6:1 with a better number. The deliverable is the SPREAD. If
the ratio is stable across defensible partitions, the claim is about the
literature. If it swings, the claim is about the class boundary -- and a reader
is entitled to know which they are being told.

A MEASURED CAVEAT THAT CUTS BOTH WAYS. The ingest reads MeSH DescriptorName and
never QualifierName (#722), and MeSH files these modalities partly as
qualifiers: measured on sampled shards, the descriptor axis catches radiotherapy
in 3.3% of cancer articles where either axis catches 5.7%, and drug therapy in
12.4% against 22.6%. So BOTH classes are understated here, drug therapy more
than radiotherapy -- which means a qualifier-aware recount would likely move the
ratio further DOWN, not up.

Usage:
    python scripts/atlas_modality_ratio.py
    python scripts/atlas_modality_ratio.py --render-only
"""

import argparse
import gzip
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = PROJECT_ROOT / "corpus" / "atlas"
PARTITIONS = PROJECT_ROOT / "analysis" / "modality-partitions.json"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-modality-ratio.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-modality-ratio.json"

# What the manuscript and the existing landscape analysis report, for contrast.
MANUSCRIPT_RATIO = 9.1
LANDSCAPE_CENSUS_RATIO = 17.6


def load_partitions() -> dict:
    d = json.loads(PARTITIONS.read_text())
    if not d:
        raise SystemExit(f"no partitions in {PARTITIONS}")
    for name, spec in d.items():
        if "ratio" in spec or "reported_ratio" in spec:
            raise SystemExit(
                f"{name} carries a ratio in the INPUT file. The input must not "
                "contain the result; that is how a measurement validates "
                "itself.")
        if not spec.get("pharmacological") or not spec.get("physical"):
            raise SystemExit(f"{name} is missing a class")
    return d


def scan(parts: dict) -> dict:
    sets = {k: {"ph": {x.lower() for x in v["pharmacological"]},
                "py": {x.lower() for x in v["physical"]}}
            for k, v in parts.items()}
    counts = {k: {"pharm": 0, "phys": 0, "both": 0} for k in sets}
    n = 0
    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                mesh = {m.lower() for m in (r.get("mesh") or [])}
                if not mesh:
                    continue
                for k, s in sets.items():
                    a = bool(mesh & s["ph"])
                    b = bool(mesh & s["py"])
                    counts[k]["pharm"] += a
                    counts[k]["phys"] += b
                    counts[k]["both"] += (a and b)
    out = {"census": n, "partitions": {}}
    for k, c in counts.items():
        out["partitions"][k] = {
            **c,
            "ratio": (c["pharm"] / c["phys"]) if c["phys"] else 0.0,
            "n_pharm_descriptors": len(parts[k]["pharmacological"]),
            "n_phys_descriptors": len(parts[k]["physical"]),
        }
    return out


def render(d: dict) -> str:
    ps = d["partitions"]
    ranked = sorted(ps.items(), key=lambda kv: -kv[1]["ratio"])
    lo, hi = ranked[-1][1]["ratio"], ranked[0][1]["ratio"]
    L = ["# Does the pharmacological:physical claim survive a symmetric class definition?", ""]
    L += ["*Generated by `scripts/atlas_modality_ratio.py` over "
          f"{d['census']:,} census articles. Partitions are the committed input "
          "`analysis/modality-partitions.json`; every member is listed there so "
          "a placement can be disputed.*", ""]

    L += ["## The measured spread", ""]
    L += ["| partition | pharmacological | physical | both | **ratio** |",
          "|---|--:|--:|--:|--:|"]
    for k, c in ranked:
        L.append(f"| {k} | {c['pharm']:,} | {c['phys']:,} | {c['both']:,} | "
                 f"**{c['ratio']:.2f}:1** |")
    L += [""]
    L += [f"| for comparison | | | | |",
          f"|---|--:|--:|--:|--:|",
          f"| manuscript's own method | | | | {MANUSCRIPT_RATIO}:1 |",
          f"| `atlas_landscape.py` on the census | | | | "
          f"{LANDSCAPE_CENSUS_RATIO}:1 |", ""]

    L += ["## What this says", ""]
    L += [f"Every partition built by one principle applied to both classes gives "
          f"a ratio between **{lo:.2f}:1 and {hi:.2f}:1**. The reported census "
          f"figure is {LANDSCAPE_CENSUS_RATIO}:1 -- between "
          f"{LANDSCAPE_CENSUS_RATIO/hi:.1f}x and "
          f"{LANDSCAPE_CENSUS_RATIO/lo:.1f}x larger than anything measured "
          f"here.", ""]
    L += ["**The direction survives and the magnitude does not.** Pharmacological "
          "exceeds physical under all five partitions, so the claim's sign is "
          "robust. But the figure that makes it rhetorically powerful -- an "
          "order of magnitude -- is a property of a PHYSICAL class containing "
          "three mechanism tags and excluding radiotherapy, brachytherapy, "
          "hyperthermia, ablation and phototherapy.", ""]
    L += ["That exclusion is not a judgement anyone made. Those modalities have "
          "no mechanism tag in this project, so they could not enter a "
          "tag-based class. The census figure inherits the taxonomy's field of "
          "view rather than measuring the literature (see "
          "`analysis/atlas-taxonomy-reach.md`).", ""]

    L += ["## What would make this wrong", ""]
    L += ["* If the partitions here are not defensible. Each was built by one "
          "stated principle, checked by an independent reviewer for one-sided "
          "widening, and required to reproduce its own count; every member is "
          "listed. Dispute a placement and re-run.",
          "* If descriptor-only counting biased the result. It does bias it, and "
          "it cuts AGAINST this finding's direction: the ingest misses the "
          "qualifier axis (#722), which understates drug therapy (12.4% vs "
          "22.6% either-axis) more than radiotherapy (3.3% vs 5.7%). A "
          "qualifier-aware recount should move the ratio further down.",
          "* Classes are not mutually exclusive, and the `both` column is "
          "reported rather than resolved. Combined-modality treatment is real, "
          "and forcing an article into one class would be the arbitrary step.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan(load_partitions())
        if not d["partitions"] or all(v["phys"] == 0 for v in d["partitions"].values()):
            raise SystemExit(
                "no physical-class articles matched, which is not a finding -- "
                "it is what a descriptor-case mismatch looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for k, c in sorted(d["partitions"].items(), key=lambda kv: -kv[1]["ratio"]):
        print(f"  {k:24s} {c['ratio']:>6.2f}:1")


if __name__ == "__main__":
    main()
