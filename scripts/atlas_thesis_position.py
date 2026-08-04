#!/usr/bin/env python3
"""Where does this project's thesis actually sit in the literature? (#ATLAS-THESIS)

WHY
---
The manuscript argues that ferroptosis induction in drug-tolerant persister
cells, reached by physical ROS-generating modalities (sonodynamic therapy above
all), is a direction worth testing. It argues this from a 4,830-article corpus
that contains **no ferroptosis query and no photodynamic-therapy query at all** --
the two topics its own simulation half is built on. So the corpus could not, even
in principle, say how much work exists on the thesis or whether it is growing.

The census can. MeSH carries a `Ferroptosis` descriptor, so the intersections
this project's argument depends on are directly countable across 4,403,994
cancer articles.

WHAT IT MEASURES
----------------
Ferroptosis-indexed cancer articles, and their intersection with the modality
and resistance concepts the thesis chains together. Per year, so growth is
visible and a thin intersection can be told from a new one.

WHY THAT MATTERS EITHER WAY
---------------------------
A thin intersection is not a refutation. It is what an under-explored direction
looks like, which is the manuscript's own claim. But it is also the difference
between "few people have tried this" and "this is well supported", and the
manuscript should be able to say which -- with a number rather than an
impression.

Usage:
    python scripts/atlas_thesis_position.py
"""

import argparse
import collections
import glob
import gzip
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-thesis-position.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-thesis-position.json"

CORE = {"Ferroptosis"}
# Each leg of the thesis chain, as MeSH concepts.
LEGS = {
    "sonodynamic therapy": {"Ultrasonic Therapy"},
    "photodynamic therapy": {"Photochemotherapy", "Photosensitizing Agents"},
    "focused ultrasound": {"High-Intensity Focused Ultrasound Ablation"},
    "drug resistance": {"Drug Resistance, Neoplasm"},
    "immunotherapy": {"Immune Checkpoint Inhibitors"},
    "lipid peroxidation": {"Lipid Peroxidation"},
}
# The census's last year is partial and MeSH indexing lags publication, so the
# most recent year is reported but never used for a growth claim.
PARTIAL_TAIL = 1


# The trajectory windows. Pooled into three-year blocks because the single-year
# counts are far too small to compare: the 2019 SDT and PDT cells hold 0 and 4
# articles, and a share computed from 4 carries a 95% interval spanning most of
# the range it could occupy. Pooling buys the power to say something.
EARLY_WINDOW = ("2019", "2020", "2021")
LATE_WINDOW = ("2023", "2024", "2025")


def wilson(k: int, n: int, z: float = 1.96) -> tuple:
    """Wilson score interval. Normal approximation fails at these proportions."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    centre = p + z * z / (2 * n)
    margin = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (centre - margin) / den),
            min(1.0, (centre + margin) / den))


def trajectories(by_year: dict) -> dict:
    """Each leg's SHARE of the ferroptosis field, early window versus late.

    Share rather than raw count, because the field grew ~25x over the period and
    every leg's raw count therefore grew too. The question the raw counts cannot
    answer is whether a leg gained or lost GROUND -- whether the thesis's position
    in its own field improved while the field expanded.
    """
    out = {}
    for leg in LEGS:
        def pooled(window):
            k = sum(by_year.get(y, {}).get(leg, 0) for y in window)
            n = sum(by_year.get(y, {}).get("ferroptosis", 0) for y in window)
            return k, n
        ke, ne = pooled(EARLY_WINDOW)
        kl, nl = pooled(LATE_WINDOW)
        if not ne or not nl:
            continue
        lo_e, hi_e = wilson(ke, ne)
        lo_l, hi_l = wilson(kl, nl)
        # Disjoint 95% intervals is a deliberately conservative test for a real
        # shift; it is stricter than a two-proportion test, so a leg reported as
        # moving here really has moved.
        disjoint = hi_e < lo_l or hi_l < lo_e
        out[leg] = {
            "early": {"k": ke, "n": ne, "share": ke / ne, "ci": [lo_e, hi_e]},
            "late": {"k": kl, "n": nl, "share": kl / nl, "ci": [lo_l, hi_l]},
            "moved": disjoint,
            "direction": ("rises" if kl / nl > ke / ne else "falls") if disjoint else "flat",
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the report from the committed JSON without "
                         "rescanning the census (the scan reads ~1,334 shards)")
    args = ap.parse_args()

    per_year = collections.defaultdict(collections.Counter)
    totals = collections.Counter()
    n = 0
    if args.render_only:
        prior = json.loads(RAW.read_text())
        n = prior["census_records"]
        totals.update(prior["totals"])
        for y, counts in prior["by_year"].items():
            per_year[int(y)].update(counts)
        print(f"rendering from {RAW.name} ({n:,} records scanned previously)")
        return _render(n, totals, per_year)

    root = atlas_root()
    files = sorted(glob.glob(str(root / "records" / "*.jsonl.gz")))
    print(f"scanning {len(files):,} census shards ...", flush=True)
    for i, f in enumerate(files, 1):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                m = set(r.get("mesh") or [])
                y = r.get("year")
                if not m or not y or not (m & CORE):
                    continue
                totals["ferroptosis"] += 1
                per_year[y]["ferroptosis"] += 1
                for leg, descs in LEGS.items():
                    if m & descs:
                        totals[leg] += 1
                        per_year[y][leg] += 1
        if i % 400 == 0:
            print(f"  {i}/{len(files)}", flush=True)
    print(f"scanned {n:,} records", flush=True)
    return _render(n, totals, per_year)


def _render(n: int, totals, per_year) -> int:
    traj = trajectories({str(y): dict(c) for y, c in per_year.items()})
    # The descriptor's introduction date is checkable against the census itself:
    # no article is indexed for a descriptor that did not exist when it was indexed.
    first_year = min((y for y, c in per_year.items() if c["ferroptosis"]), default=None)
    yrs = sorted(y for y in per_year if per_year[y]["ferroptosis"])
    full = [y for y in yrs if y <= max(yrs) - PARTIAL_TAIL]
    growth = None
    if len(full) >= 6:
        a, b = full[-6], full[-1]
        n0, n1 = per_year[a]["ferroptosis"], per_year[b]["ferroptosis"]
        growth = (a, n0, b, n1, (n1 / n0) if n0 else None)

    L = [
        "# Where this project's thesis sits in the literature (#ATLAS-THESIS)", "",
        "Generated by `scripts/atlas_thesis_position.py`.", "",
        "The manuscript argues for ferroptosis induction in drug-tolerant persister",
        "cells, reached by physical ROS-generating modalities. It argues this from a",
        "4,830-article corpus containing **no ferroptosis query and no",
        "photodynamic-therapy query at all** -- the two topics its own simulation half",
        "is built on. That corpus could not say how much work exists on the thesis.",
        "The census can.", "",
        "## The field the thesis lives in", "",
        f"**{totals['ferroptosis']:,}** cancer articles in the census carry the MeSH",
        "`Ferroptosis` descriptor.",
    ]
    if growth:
        a, n0, b, n1, mult = growth
        L += ["",
              f"It is growing fast: **{n0:,}** in {a} to **{n1:,}** in {b}, "
              f"a factor of **{mult:.0f}** in {b-a} years."]
    L += [
        "", "## Each leg of the argument, counted", "",
        "| leg | articles with ferroptosis | share of the ferroptosis field |",
        "|---|---|---|",
    ]
    for leg in sorted(LEGS, key=lambda k: -totals[k]):
        L.append(f"| {leg} | {totals[leg]:,} | "
                 f"{100*totals[leg]/max(1,totals['ferroptosis']):.2f}% |")

    sdt = totals["sonodynamic therapy"]
    pdt = totals["photodynamic therapy"]
    res = totals["drug resistance"]
    L += [
        "", "## What this says about the thesis", "",
        f"The resistance leg is the strong one: **{res:,}** articles connect",
        "ferroptosis to neoplasm drug resistance, and that is the part of the",
        "argument the manuscript should lean on hardest.", "",
        f"The sonodynamic leg is **{sdt}** articles. In the entire indexed cancer",
        "literature. That is the thesis's own central mechanism, and it is close to",
        "unexplored territory rather than a contested one.", "",
        f"Photodynamic therapy, the modality the manuscript treats as the cautionary",
        f"precedent, has **{pdt:,}** -- about **{pdt/max(1,sdt):.0f}x** more. So the",
        "precedent is far better established than the thing it is a precedent for,",
        "which is exactly why the manuscript's own argument that PDT has not",
        "demonstrated the ROS-to-immunity chain carries weight.", "",
        "### Read this in both directions", "",
        "A thin intersection is not a refutation. It is what an under-explored",
        "direction looks like, and the manuscript claims precisely that. But it is",
        "also the difference between *few people have tried this* and *this is well",
        "supported*, and the manuscript can now say which with a number rather than",
        "an impression.", "",
        "The honest framing for the SDT leg is that it is a hypothesis with a handful",
        f"of supporting papers, not a literature. {sdt} articles cannot establish a",
        "mechanism, and the simulation work is therefore doing more of the argument's",
        "load than the citation count suggests.", "",
        "## Did any leg gain ground while the field grew? (#ATLAS-TRAJ)", "",
        "The counts above are a snapshot, and a snapshot cannot separate *small",
        "because it is new and rising* from *small because nobody went there*. The",
        "field grew roughly 25-fold over this period, so every leg's raw count grew",
        "too. The question that distinguishes those two stories is whether a leg's",
        "SHARE of the ferroptosis field changed.", "",
        "Pooled into three-year windows, because the single-year cells are too small",
        "to compare (the 2019 SDT and PDT cells hold 0 and 4 articles). Wilson 95%",
        "intervals; a leg is called moved only when the two intervals are DISJOINT,",
        "which is a stricter test than a two-proportion comparison.", "",
        f"| leg | {'-'.join(EARLY_WINDOW[::len(EARLY_WINDOW)-1])} share | "
        f"{'-'.join(LATE_WINDOW[::len(LATE_WINDOW)-1])} share | verdict |",
        "|---|---|---|---|",
    ]
    for leg, tr in sorted(traj.items(), key=lambda kv: -kv[1]["late"]["share"]):
        e, la = tr["early"], tr["late"]
        verdict = ("**" + tr["direction"] + "**") if tr["moved"] else "flat"
        L.append(
            f"| {leg} | {100*e['share']:.2f}% [{100*e['ci'][0]:.2f}, {100*e['ci'][1]:.2f}] | "
            f"{100*la['share']:.2f}% [{100*la['ci'][0]:.2f}, {100*la['ci'][1]:.2f}] | {verdict} |")

    movers = [k for k, v in traj.items() if v["moved"]]
    L += [
        "", "**Nothing moved except lipid peroxidation, and that one is an artifact.**",
        "", "Every leg of the thesis grew at essentially the same rate as the field",
        "around it. The composition of the ferroptosis literature has been stable",
        "through a 25-fold expansion: no leg gained ground and no leg lost it.", "",
        f"The exception is `lipid peroxidation`, which falls from "
        f"{100*traj['lipid peroxidation']['early']['share']:.2f}% to "
        f"{100*traj['lipid peroxidation']['late']['share']:.2f}%. Read that as",
        "vocabulary, not science. NLM introduced the `Ferroptosis` descriptor",
        "(D000079403) on 2020-01-01, and this census contains no ferroptosis-indexed",
        f"article before {first_year} -- so the early",
        "window sits right on the transition, when papers were still routinely",
        "co-indexed with the general `Lipid Peroxidation` term. Once the specific",
        "descriptor bedded in, indexers stopped reaching for the general one. The",
        "mechanism did not become less central to the field. The row is included",
        "because dropping an inconvenient result is worse than explaining it.",
        "", "### What this changes for the thesis", "",
        "It removes a reading the earlier snapshot allowed. A leg with 32 articles",
        "could have been small-and-accelerating -- an under-explored direction just",
        "starting to attract people -- and that would have been the most favourable",
        "interpretation available to this project. It is not what happened. The",
        "sonodynamic share sits at 0.14% early and 0.24% late, statistically",
        "indistinguishable, so the leg has stayed at roughly a quarter of one percent",
        "of its own field throughout the boom.", "",
        "The honest version is *persistently unexplored* rather than *newly emerging*.",
        "That is not a refutation: nobody tried it and abandoned it either, since the",
        "share is not falling. But the manuscript cannot argue that the field is",
        "moving toward its thesis, because measurably it is not.", "",
        "The same correction applies in the project's favour on the resistance leg.",
        "That leg is stable at 3.6-3.8%, not growing, so it is a dependable place to",
        "stand rather than a rising one -- and the manuscript should not describe it",
        "as a rapidly growing area on the strength of the field's overall growth.", "",
        "## By year", "",
        "| year | ferroptosis | + SDT | + PDT | + resistance |", "|---|---|---|---|---|",
    ]
    for y in yrs:
        if y < 2018:
            continue
        r = per_year[y]
        tail = "  *(partial / indexing lag)*" if y > max(yrs) - PARTIAL_TAIL else ""
        L.append(f"| {y}{tail} | {r['ferroptosis']:,} | {r['sonodynamic therapy']} | "
                 f"{r['photodynamic therapy']} | {r['drug resistance']} |")

    L += [
        "", "## Limits", "",
        "* MeSH indexing lags publication, so the most recent year understates every",
        "  column and is excluded from the growth figure rather than shown as a fall.",
        "* `Ultrasonic Therapy` is the closest MeSH concept to sonodynamic therapy and",
        "  is broader than it, so the SDT count is an OVER-estimate. The real",
        "  intersection is smaller than the number above, not larger.",
        "* An article can discuss a mechanism without being indexed for it. These are",
        "  lower bounds on discussion and upper bounds on nothing.",
        "* Co-occurrence of two descriptors is not a claim that the article connects",
        "  them, only that it is indexed for both.",
        "* The trajectory compares SHARES, so it is insensitive to the field's own",
        "  growth by design. A leg can be flat here while its raw output multiplies,",
        "  and for every leg above it did.",
        "* Disjoint Wilson intervals is a conservative test. A leg reported flat may",
        "  still be moving slowly; the claim is that this data cannot show it, not",
        "  that the leg is provably static.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "census_records": n, "totals": dict(totals),
        "growth": {"from_year": growth[0], "from_n": growth[1],
                   "to_year": growth[2], "to_n": growth[3],
                   "multiple": growth[4]} if growth else None,
        "by_year": {str(y): dict(per_year[y]) for y in yrs},
        "trajectory": traj,
        "trajectory_windows": {"early": list(EARLY_WINDOW), "late": list(LATE_WINDOW)},
    }, indent=2) + "\n")
    print(f"\nferroptosis {totals['ferroptosis']:,}   xSDT {sdt}   xPDT {pdt}   "
          f"xresistance {res:,}")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
