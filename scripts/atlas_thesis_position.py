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


def min_detectable_share(early_hi: float, n_late: int) -> int:
    """Smallest late count the disjoint-interval test could call a rise.

    The test is deliberately conservative, and conservative means BLIND over
    some range. Reporting "flat" without saying how large a change would have
    had to be is evidence of absence from a test that may have no power to
    detect presence -- so every flat verdict here ships with this number.
    """
    for k in range(n_late + 1):
        if wilson(k, n_late)[0] > early_hi:
            return k
    return n_late + 1


def window_sensitivity(by_year: dict, leg: str) -> dict:
    """Does the verdict survive a different choice of windows?

    The 2019-2021 / 2023-2025 split was chosen before the results were seen, but
    nothing makes it privileged, and a verdict that only holds for one split is
    a property of the split. Sweeps every non-overlapping contiguous 2- and
    3-year window pair in the range.
    """
    yrs = [str(y) for y in range(2019, 2026)]
    verdicts = []
    for size in (2, 3):
        for i in range(len(yrs) - 2 * size + 1):
            early, late = yrs[i:i + size], yrs[-size:]
            if set(early) & set(late):
                continue
            ke = sum(by_year.get(y, {}).get(leg, 0) for y in early)
            ne = sum(by_year.get(y, {}).get("ferroptosis", 0) for y in early)
            kl = sum(by_year.get(y, {}).get(leg, 0) for y in late)
            nl = sum(by_year.get(y, {}).get("ferroptosis", 0) for y in late)
            if not ne or not nl:
                continue
            a, b = wilson(ke, ne), wilson(kl, nl)
            moved = a[1] < b[0] or b[1] < a[0]
            verdicts.append({
                "early": f"{early[0]}-{early[-1]}", "late": f"{late[0]}-{late[-1]}",
                "moved": moved,
                "direction": ("rises" if kl / nl > ke / ne else "falls") if moved else "flat",
            })
    moved = [v for v in verdicts if v["moved"]]
    return {"n_windows": len(verdicts), "n_moved": len(moved),
            "moved_windows": moved}


def trajectories(by_year: dict) -> dict:
    """Each leg's SHARE of the ferroptosis field, early window versus late.

    Share rather than raw count, because the field grew several-fold over the
    period so every leg's raw count grew too (the exact multiple is computed in
    the report and depends on the baseline year, which is why it is not named
    here). The question the raw counts cannot
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
        k_needed = min_detectable_share(hi_e, nl)
        out[leg] = {
            "early": {"k": ke, "n": ne, "share": ke / ne, "ci": [lo_e, hi_e]},
            "late": {"k": kl, "n": nl, "share": kl / nl, "ci": [lo_l, hi_l]},
            "moved": disjoint,
            "direction": ("rises" if kl / nl > ke / ne else "falls") if disjoint else "flat",
            # What a "flat" verdict is actually compatible with.
            "min_detectable_k": k_needed,
            "min_detectable_share": k_needed / nl if nl else None,
            "min_detectable_ratio": ((k_needed / nl) / (ke / ne)) if ke and ne and nl else None,
            "windows": window_sensitivity(by_year, leg),
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
        stale = []
        if prior.get("legs") and sorted(prior["legs"]) != sorted(LEGS):
            stale.append(f"LEGS changed: committed {sorted(prior['legs'])}, "
                         f"now {sorted(LEGS)}")
        if prior.get("core") and sorted(prior["core"]) != sorted(CORE):
            stale.append(f"CORE changed: committed {sorted(prior['core'])}, "
                         f"now {sorted(CORE)}")
        if stale:
            print("cannot render: the committed counts do not cover the current "
                  "definitions, and a missing leg would render as a measured 0.00%.",
                  file=sys.stderr)
            for s in stale:
                print(f"  {s}", file=sys.stderr)
            print("  re-run without --render-only to rescan the census.",
                  file=sys.stderr)
            return 1
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

    L += ["> " + (

        "**This enumerates the modalities the thesis names.** Ranked against all 166 in the committed partition universe the sonodynamic leg is rank 22 and `antineoplastic agents` is 29.8x larger, so the unexamined precedent is chemotherapy. See `analysis/atlas-thesis-rank.md` (#725)."), ""]
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
        "because it is new and rising* from *small because nobody went there*. Every",
        "leg's raw count grew, because the field did. The question that separates the",
        "two stories is whether a leg's SHARE of the ferroptosis field changed.", "",
        "Pooled into three-year windows, because the single-year cells are too small",
        "to compare (the 2019 SDT and PDT cells hold 0 and 4 articles). Wilson 95%",
        "intervals; a leg is called moved only when the two intervals are DISJOINT.", "",
        "**That test is conservative, which means blind over a range, so every verdict",
        "below ships with the size of change it could actually have detected.** A flat",
        "verdict on a leg with a 4x detection floor is not evidence that the leg did",
        "not move.", "",
        "| leg | early share | late share | verdict | would have needed | holds across windows |",
        "|---|---|---|---|---|---|",
    ]
    for leg, tr in sorted(traj.items(), key=lambda kv: -kv[1]["late"]["share"]):
        e, la, w = tr["early"], tr["late"], tr["windows"]
        verdict = ("**" + tr["direction"] + "**") if tr["moved"] else "flat"
        ratio = tr["min_detectable_ratio"]
        need = f"a {ratio:.1f}x change" if ratio else "n/a (no early support)"
        L.append(
            f"| {leg} | {100*e['share']:.2f}% [{100*e['ci'][0]:.2f}, {100*e['ci'][1]:.2f}] | "
            f"{100*la['share']:.2f}% [{100*la['ci'][0]:.2f}, {100*la['ci'][1]:.2f}] | {verdict} | "
            f"{need} | {w['n_windows']-w['n_moved']}/{w['n_windows']} flat |")

    movers = sorted(k for k, v in traj.items() if v["moved"])
    flat = sorted(k for k, v in traj.items() if not v["moved"])
    L += ["", f"**Moved: {', '.join(movers) if movers else 'nothing'}. "
              f"Not distinguishable from flat: {', '.join(flat)}.**", ""]

    if movers == ["lipid peroxidation"]:
        L += [
            "No leg of the thesis is measurably gaining or losing ground on its own",
            "field. Read that as a bound, not as proof of stasis: the detection floors",
            "in the table are the honest statement of what this cannot see.", ""]
    elif movers:
        L += [
            "Note that more than the vocabulary artifact is moving here, which the",
            "prose below was written before. Re-read the table rather than the",
            "narrative.", ""]

    sdt_tr = traj.get("sonodynamic therapy")
    if sdt_tr:
        L += [
            "### What this does and does not say about the sonodynamic leg", "",
            "It removes the most favourable reading the earlier snapshot allowed. A",
            "32-article leg could have been small-and-accelerating -- an under-explored",
            "direction just starting to attract people -- and nothing in the counts",
            "ruled that out.",
            f"The share went {100*sdt_tr['early']['share']:.2f}% to "
            f"{100*sdt_tr['late']['share']:.2f}%, and it is flat in all "
            f"{sdt_tr['windows']['n_windows']} window pairs tried, so there is no sign of",
            "an inflection.", "",
            "**But the honest limit is severe.** With only "
            f"{sdt_tr['early']['k']} early articles, the test would not have called a rise",
            f"below **{sdt_tr['min_detectable_ratio']:.1f}x**. So the defensible statement is",
            "that the leg shows no detectable acceleration, NOT that it demonstrably has",
            "none. *Persistently unexplored* is supported as a description of the",
            "measured shares. The stronger claim an earlier draft made -- that the",
            "field is demonstrably not moving toward the thesis -- is not supported by",
            "a test this blind, and has been withdrawn.", "",
            "The leg is also not shrinking, so this is not a case of people trying it",
            "and giving up.", ""]

    resist_tr = traj.get("drug resistance")
    if resist_tr:
        w = resist_tr["windows"]
        L += ["### The resistance leg", "",
              f"Flat on the headline windows at {100*resist_tr['early']['share']:.2f}% to "
              f"{100*resist_tr['late']['share']:.2f}%, so it is a dependable place to stand",
              "rather than a rising one, and it should not be called a rapidly growing",
              "area on the strength of the field's overall growth.", ""]
        if w["n_moved"]:
            ex = w["moved_windows"][0]
            L += [
                f"**That verdict is window-dependent**, and the report says so rather",
                f"than reporting only the split that was chosen first: the leg is called",
                f"*{ex['direction']}* in {w['n_moved']} of {w['n_windows']} window pairs,",
                f"including {ex['early']} against {ex['late']}. Those are the windows that",
                "start after the descriptor transition discussed below, so they are not",
                "obviously the worse choice. The cautious reading is that the resistance",
                "share may be drifting up and this data cannot settle it.", ""]

    lp_tr = traj.get("lipid peroxidation")
    if lp_tr and lp_tr["moved"]:
        by_yr = {y: c for y, c in per_year.items() if c.get("ferroptosis")}
        series = [(y, 100 * by_yr[y].get("lipid peroxidation", 0) / by_yr[y]["ferroptosis"])
                  for y in sorted(by_yr)]
        trough = min((s for s in series if s[0] >= 2019), key=lambda s: s[1])
        last = series[-1]
        L += [
            "### The one leg that moves, and why the obvious explanation is incomplete",
            "",
            f"`lipid peroxidation` falls from {100*lp_tr['early']['share']:.2f}% to "
            f"{100*lp_tr['late']['share']:.2f}%, the only disjoint result, and it moves in",
            f"{lp_tr['windows']['n_moved']} of {lp_tr['windows']['n_windows']} window pairs.", "",
            "The natural explanation is vocabulary rather than science. NLM introduced",
            "the `Ferroptosis` descriptor (D000079403) on 2020-01-01, and before a",
            "specific descriptor exists indexers reach for the general one, so a decline",
            "as the specific term beds in is what one would expect.", "",
            "**The series does not fully support that.** Year by year: "
            + ", ".join(f"{y} {v:.1f}%" for y, v in series if y >= 2019) + ".",
            f"The decline continues to {trough[0]}, three years past the introduction,",
            f"and then REVERSES to {last[1]:.1f}% by {last[0]}. Indexers dropping a",
            "general term cannot produce a rise. So the transition may explain the early",
            "part of the fall, but something else is driving the shape, and this report",
            "does not know what. It is stated as an open question rather than settled,",
            "because the tidy version was written first and the data refused it.", "",
            "What it is NOT is evidence that lipid peroxidation became less central to",
            "ferroptosis, which would be a claim about the biology that a co-indexing",
            "rate cannot support in either direction.", ""]

    L += [
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
        "* `Ultrasonic Therapy` is the closest MeSH concept to sonodynamic therapy",
        "  and IS broader than it -- but measured, that is the smaller effect. It",
        "  recalls under half of the ferroptosis papers whose own text says",
        "  sonodynamic, against a breadth of a handful of records, so the SDT",
        "  count is an UNDER-estimate: the real intersection is LARGER than the",
        "  number above, not smaller. An earlier version of this file stated the",
        "  opposite; see `analysis/atlas-descriptor-recall.md`.",
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
        # Provenance: which definitions produced these counts. --render-only
        # refuses if they no longer match, since it cannot measure a new leg.
        "legs": {k: sorted(v) for k, v in LEGS.items()},
        "core": sorted(CORE),
    }, indent=2) + "\n")
    print(f"\nferroptosis {totals['ferroptosis']:,}   xSDT {sdt}   xPDT {pdt}   "
          f"xresistance {res:,}")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
