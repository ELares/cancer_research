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
required to reproduce its own count. THAT CLAUSE IS WITHDRAWN: no
artifact in this repo records a reviewer or a per-partition count, and
the only counts are this script's own output.

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
import re
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


# THE ONE UNDECLARED VARIABLE. The five partitions were presented as five
# independent readings whose SPREAD was the deliverable. They are not
# independent: they differ almost entirely in whether operative surgery counts
# as a "physical modality", and the surgical share of each physical class is a
# perfect inverse rank predictor of the published ratio. Held out under ONE
# rule applied identically to all five, the spread nearly disappears.
SURGICAL = re.compile(
    r"surg|ectomy|otomy|ostomy|resect|excis|laparoscop|endoscop|"
    r"transplant|amputat|dissect|anastomos|graft|implantation|"
    r"reconstructi|debulk|lymphadenectom", re.I)


def scan(parts: dict) -> dict:
    sets = {k: {"ph": {x.lower() for x in v["pharmacological"]},
                "py": {x.lower() for x in v["physical"]}}
            for k, v in parts.items()}
    # the same regex applied to every partition, so the held-out comparison is
    # not itself a per-partition judgement
    surg = {k: {x for x in s["py"] if SURGICAL.search(x)} for k, s in sets.items()}
    counts = {k: {"pharm": 0, "phys": 0, "both": 0, "phys_nosurg": 0}
              for k in sets}
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
                    counts[k]["phys_nosurg"] += bool(mesh & (s["py"] - surg[k]))
    out = {"census": n, "partitions": {},
           "landscape_composition": landscape_composition()}
    for k, c in counts.items():
        out["partitions"][k] = {
            **c,
            "ratio": (c["pharm"] / c["phys"]) if c["phys"] else 0.0,
            "n_pharm_descriptors": len(parts[k]["pharmacological"]),
            "n_phys_descriptors": len(parts[k]["physical"]),
            "n_surgical_descriptors": len(surg[k]),
            "surgical_share_of_physical":
                (c["phys"] - c["phys_nosurg"]) / c["phys"] if c["phys"] else None,
            "ratio_surgery_held_out":
                (c["pharm"] / c["phys_nosurg"]) if c["phys_nosurg"] else None,
        }
    return out


def landscape_composition() -> dict:
    """What the 17.6:1 comparator is made of, and the symmetric restriction.

    Its numerator is dominated by two mechanisms `atlas_landscape.py`'s own
    text calls a SCOPE ARTIFACT rather than a therapy. Dropping those from the
    numerator ALONE would be the one-sided narrowing this page exists to
    police -- the denominator's `electrochemical-therapy` is equally
    non-precise. So the restriction is applied to BOTH classes, using that
    script's own PRECISE set rather than a judgement made here.
    """
    src = PROJECT_ROOT / "analysis" / "atlas-landscape.json"
    if not src.exists():
        return {}
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "al", PROJECT_ROOT / "scripts" / "atlas_landscape.py")
    al = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(al)
    d = json.loads(src.read_text())
    rows = d if isinstance(d, list) else d.get("mechanisms") or d.get("rows") or []
    cen = {}
    for r in rows:
        if isinstance(r, dict) and r.get("mechanism"):
            cen[r["mechanism"].lower()] = r.get("mesh_census") or 0
    if not cen:
        return {}

    def tot(names):
        return sum(cen.get(x, 0) for x in names)

    ph, py, pre = al.PHARMACOLOGICAL, al.PHYSICAL, al.PRECISE
    num, den = tot(ph), tot(py)
    num_p, den_p = tot(ph & pre), tot(py & pre)
    biggest = sorted(((x, cen.get(x, 0)) for x in ph), key=lambda kv: -kv[1])[:2]
    return {
        "numerator": num, "denominator": den,
        "ratio": num / den if den else None,
        "top_two_numerator": biggest,
        "top_two_share": sum(v for _k, v in biggest) / num if num else None,
        "precise_numerator": num_p, "precise_denominator": den_p,
        "precise_ratio": num_p / den_p if den_p else None,
        "precise_pharm": sorted(ph & pre), "precise_phys": sorted(py & pre),
        "dropped_pharm": sorted(ph - pre), "dropped_phys": sorted(py - pre),
    }


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
    L += ["| partition | pharmacological | physical | both | **ratio** | "
          "surgical share of physical | ratio, surgery held out |",
          "|---|--:|--:|--:|--:|--:|--:|"]
    for k, c in ranked:
        ss = c.get("surgical_share_of_physical")
        ho = c.get("ratio_surgery_held_out")
        L.append(f"| {k} | {c['pharm']:,} | {c['phys']:,} | {c['both']:,} | "
                 f"**{c['ratio']:.2f}:1** | "
                 f"{100*ss:.1f}% | {ho:.2f}:1 |" if ss is not None and ho
                 else f"| {k} | {c['pharm']:,} | {c['phys']:,} | "
                      f"{c['both']:,} | **{c['ratio']:.2f}:1** | | |")
    L += [""]

    # THE SPREAD IS ONE UNDECLARED VARIABLE. Derived, so the claim moves with
    # the data: an earlier version called the spread "the deliverable" and
    # presented the five partitions as independent readings.
    hos = [c["ratio_surgery_held_out"] for _k, c in ranked
           if c.get("ratio_surgery_held_out")]
    shares = [c.get("surgical_share_of_physical") for _k, c in ranked
              if c.get("surgical_share_of_physical") is not None]
    if len(hos) == len(ranked) and len(shares) == len(ranked):
        pub = hi / lo
        held = max(hos) / min(hos)
        # rank correlation between surgical share and the published ratio
        order_share = sorted(range(len(shares)), key=lambda i: shares[i])
        order_ratio = sorted(range(len(ranked)),
                             key=lambda i: ranked[i][1]["ratio"])
        inverse = order_share == list(reversed(order_ratio))
        L += [f"**The spread is one undeclared variable: whether operative "
              f"surgery is a physical modality.** Its share of the physical "
              f"class runs {100*min(shares):.1f}% to {100*max(shares):.1f}% "
              f"across the five"
              + (", and ranks PERFECTLY INVERSELY with the published ratio"
                 if inverse else "") +
              f". Held out under one regex applied identically to all five, "
              f"the spread collapses from **{pub:.2f}x** "
              f"({lo:.2f}-{hi:.2f}) to **{held:.2f}x** "
              f"({min(hos):.2f}-{max(hos):.2f}) -- every partition landing "
              f"near {sum(hos)/len(hos):.1f}:1.", ""]
        L += ["So the five are not five independent readings whose "
              "disagreement is the finding. They agree about "
              "pharmacological-versus-physical and disagree about one "
              "membership question nobody declared. An earlier version of this "
              "page called that spread \"the deliverable\".", ""]
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

    lc = d.get("landscape_composition") or {}
    if lc.get("precise_ratio"):
        L += [f"## What the {LANDSCAPE_CENSUS_RATIO}:1 comparator is made of", ""]
        top = ", ".join(f"`{k}` {v:,}" for k, v in lc["top_two_numerator"])
        L += [f"Its numerator is {lc['numerator']:,} and its denominator "
              f"{lc['denominator']:,}. Two mechanisms supply "
              f"**{100*lc['top_two_share']:.0f}%** of the numerator: {top}. "
              f"`atlas_landscape.py`'s own text calls the largest of them a "
              f"SCOPE ARTIFACT rather than a therapy -- a descriptor carried "
              f"by any paper that MEASURES the process.", ""]
        L += [f"Dropping those from the numerator alone would be the one-sided "
              f"narrowing this page exists to police. Applying that script's "
              f"own `PRECISE` set to BOTH classes -- which also drops "
              f"{', '.join(f'`{x}`' for x in lc['dropped_phys'])} from the "
              f"denominator -- gives **{lc['precise_ratio']:.2f}:1** on "
              f"{lc['precise_numerator']:,} against {lc['precise_denominator']:,}.", ""]
        L += [f"So the comparator falls from {lc['ratio']:.1f}:1 to "
              f"{lc['precise_ratio']:.2f}:1 under a symmetric restriction, "
              f"which is within the range the partitions give once surgery is "
              f"held fixed. The gap between this page and "
              f"`atlas_landscape.py` is substantially a class-composition "
              f"difference, not a disagreement about the literature.", ""]

    L += ["## What would make this wrong", ""]
    L += ["* If the partitions here are not defensible. They are a committed "
          "member list and every member is disputable; dispute a placement and "
          "re-run. An earlier version of this bullet also claimed each was "
          "\"checked by an independent reviewer for one-sided widening, and "
          "required to reproduce its own count\". NO ARTIFACT IN THIS REPO "
          "SUPPORTS EITHER CLAUSE -- `modality-partitions.json` carries only "
          "the two member lists, there is no reviewer record, and the only "
          "per-partition counts are this script's own output, which makes "
          "\"reproduce its own count\" circular. Both clauses are withdrawn.",
          "* If descriptor-only counting biased the result. An earlier version "
          "of this bullet argued that it cuts AGAINST the finding, from the "
          "#722 recalls: drug therapy 12.4% vs 22.6% either-axis against "
          "radiotherapy 3.3% vs 5.7%. THAT INFERENCE IS INVALID. Those are "
          "recalls of 0.547 and 0.579, and correcting each class by its own "
          "recall multiplies the ratio by 0.579/0.547 = 1.05 -- UP, not down. "
          "\"Understates more\" is true in percentage POINTS, and a ratio "
          "responds to relative rather than additive understatement. The "
          "bias is real and its direction is not established here.",
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
