#!/usr/bin/env python3
"""Every modality's kill is reported; none of them has its selectivity. (#728)

WHAT THIS IS AND, MORE IMPORTANTLY, WHAT IT IS NOT
---------------------------------------------------
Every therapy in this project's corpus fails clinically for the same reason: it
harms normal tissue before it clears the tumour. Selectivity is the whole
problem. The engine reports kill for every modality it models and selectivity
for none, and the word "toxicity" appears once in the ~48,900-word manuscript.

A model that computes kill without harm will systematically prefer whichever
modality kills hardest, and that bias runs toward the thesis: ferroptosis
inducers are attractive precisely on the claim that normal cells resist them.

SO THIS IS NOT A THERAPEUTIC INDEX, AND SAYING OTHERWISE WOULD BE THE ERROR
---------------------------------------------------------------------------
The engine has a `Stromal` phenotype and it is tempting to divide by it. It
represents CANCER-ASSOCIATED FIBROBLASTS: a tumour-resident, metabolically
altered cell type recruited by the tumour, whose parameters were chosen to model
shielding. It is not healthy tissue at a distant site. Reading a
tumour-over-stromal ratio as a therapeutic index would be a category error, and
an earlier version of #728 proposed exactly that.

What this measures is the CAF-versus-tumour contrast, which is a real and
reportable property of the model, and which the engine has always been able to
produce and has never published.

WHAT A REAL THERAPEUTIC INDEX WOULD NEED, stated so the gap stays visible: a
normal-tissue phenotype distinct from Stromal, with a parameter set traceable to
the ferroptosis-resistance literature (ACSL4-low, MUFA-high, GSH-high). No such
phenotype exists, so the project cannot currently check its own most
load-bearing selectivity assumption.

Control is included as a baseline. A contrast computed against a phenotype that
dies on its own is not a contrast about treatment.

Usage:
    python scripts/engine_selectivity.py
    python scripts/engine_selectivity.py --n 20000
"""

import argparse
import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_MD = PROJECT_ROOT / "analysis" / "engine-selectivity.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "engine-selectivity.json"

TUMOUR = ["Glycolytic", "OXPHOS", "Persister", "PersisterNrf2"]
NON_TUMOUR = "Stromal"
TREATMENTS = ["Control", "RSL3", "SDT", "PDT"]
SEED = 20260817
# INDEPENDENT SEEDS MUST BE SPACED BY 2n. `sim_batch` draws cell i from
# `seed + 2i` and its simulation RNG from `seed + 2i + 1`, so one run consumes
# the whole span `seed .. seed + 2n - 1` (lib.rs, sim_batch). Two runs whose
# seeds differ by less than 2n therefore SHARE almost every cell -- seed and
# seed+2 differ by one cell in n -- and anyone "checking robustness" by
# nudging the seed gets a false confirmation. Nothing said so.
N_SEEDS = 8


RESISTANCE_AXES = {
    # axis -> True if a LOWER value means more ferroptosis-resistant
    "iron": True, "basal_ros": True, "lipid_unsat": True,
    "gsh": False, "gpx4": False, "fsp1": False, "nrf2": False,
}


def phenotype_params() -> dict:
    """Each phenotype's mean parameters, parsed from cell.rs.

    The withdrawn fingerprint claim was given a derived control while the
    SURVIVING claim -- that `Stromal` encodes the resistance assumption --
    stayed bare prose in a page whose thesis is derive-don't-assert. This
    makes it checkable.
    """
    src = (PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src" / "cell.rs")
    if not src.exists():
        return {}
    text = src.read_text(errors="ignore")
    out = {}
    for m in re.finditer(r"Phenotype::(\w+)\s*=>\s*Cell\s*\{(.*?)\n\s*\}", text, re.S):
        name, body = m.group(1), m.group(2)
        vals = {}
        for f in RESISTANCE_AXES:
            mm = re.search(rf"\b{f}:\s*norm\(rng,\s*([0-9.]+)", body)
            if mm:
                vals[f] = float(mm.group(1))
        if vals:
            out[name] = vals
    return out


def resistance_rank(params: dict, target: str) -> dict:
    """On how many axes is `target` the most ferroptosis-resistant phenotype?"""
    if target not in params:
        return {}
    axes, wins, detail = 0, 0, {}
    for f, lower_is_resistant in RESISTANCE_AXES.items():
        vals = {k: v[f] for k, v in params.items() if f in v}
        if target not in vals or len(vals) < 2:
            continue
        axes += 1
        best = min(vals.values()) if lower_is_resistant else max(vals.values())
        won = abs(vals[target] - best) < 1e-9
        wins += won
        detail[f] = {"value": vals[target], "most_resistant": won}
    return {"axes": axes, "wins": wins, "detail": detail}


def run(n: int) -> dict:
    import ferroptosis_core as fc
    pp = phenotype_params()
    out = {"n_per_cell_type": n, "seed": SEED, "seed_stride": 2 * n,
           "phenotype_params": pp,
           "stromal_resistance": resistance_rank(pp, NON_TUMOUR),
           "n_seeds": N_SEEDS, "treatments": {}}
    for t in TREATMENTS:
        row = {}
        for ph in TUMOUR + [NON_TUMOUR]:
            r = fc.sim_batch(ph, t, n, SEED)
            row[ph] = {"death_rate": r["death_rate"],
                       "ci_low": r["ci_low"], "ci_high": r["ci_high"]}
        # the contrast, computed against the WORST-killed tumour phenotype and
        # the best, so a single flattering phenotype cannot carry it
        tum = [row[p]["death_rate"] for p in TUMOUR]
        caf = row[NON_TUMOUR]["death_rate"]
        row["_contrast"] = {
            "tumour_max": max(tum), "tumour_min": min(tum), "caf": caf,
            "ratio_best_case": (max(tum) / caf) if caf > 0 else None,
            "ratio_worst_case": (min(tum) / caf) if caf > 0 else None,
        }
        out["treatments"][t] = row

    # ACROSS GENUINELY DISJOINT SEEDS, spaced by 2n so no two runs share a
    # cell. The published ratios carried no interval, and a re-run at a
    # neighbouring seed would have confirmed them by construction.
    spread = {}
    for t in TREATMENTS:
        best, worst = [], []
        for k in range(N_SEEDS):
            # k starts at 1 so the SPREAD is disjoint from the point
            # estimate above. An earlier version started at 0, so the
            # published point was a MEMBER of the sample it was being
            # compared against, and "sits near the bottom of the range"
            # described a 37.5%-likely outcome for any member.
            s = SEED + (k + 1) * 2 * n
            tum, caf = [], None
            for ph in TUMOUR + [NON_TUMOUR]:
                r = fc.sim_batch(ph, t, n, s)
                if ph == NON_TUMOUR:
                    caf = r["death_rate"]
                else:
                    tum.append(r["death_rate"])
            if caf and caf > 0:
                best.append(max(tum) / caf)
                worst.append(min(tum) / caf)
        if best:
            spread[t] = {"best_lo": min(best), "best_hi": max(best),
                         "worst_lo": min(worst), "worst_hi": max(worst),
                         "n_seeds": len(best)}
    out["seed_spread"] = spread
    return out


def render(d: dict) -> str:
    L = ["# What each modality kills, and what it spares", ""]
    L += ["*Generated by `scripts/engine_selectivity.py`. "
          f"{d['n_per_cell_type']:,} cells per phenotype per treatment, seed "
          f"{d['seed']}.*", ""]

    L += ["The engine reports kill for every modality it models and selectivity "
          "for none. This is the contrast it has always been able to produce.", ""]

    L += ["| treatment | " + " | ".join(TUMOUR) + f" | {NON_TUMOUR} (CAF) |",
          "|---|" + "--:|" * (len(TUMOUR) + 1)]
    for t in TREATMENTS:
        row = d["treatments"][t]
        cells = " | ".join(f"{100*row[p]['death_rate']:.2f}%" for p in TUMOUR)
        L.append(f"| {t} | {cells} | **{100*row[NON_TUMOUR]['death_rate']:.2f}%** |")
    L += [""]

    L += ["## The contrast, best and worst case", ""]
    L += ["A ratio taken against the most-killed tumour phenotype flatters the "
          "modality; against the least-killed it does the opposite. Both are "
          "given, because a single phenotype should not carry the claim.", ""]
    L += ["| treatment | best-case tumour:CAF | worst-case tumour:CAF |",
          "|---|--:|--:|"]
    for t in TREATMENTS:
        c = d["treatments"][t]["_contrast"]
        b = f"{c['ratio_best_case']:.1f}x" if c["ratio_best_case"] else "n/a"
        w = f"{c['ratio_worst_case']:.1f}x" if c["ratio_worst_case"] else "n/a"
        L.append(f"| {t} | {b} | {w} |")
    L += [""]

    sp = d.get("seed_spread") or {}
    if sp:
        L += [f"Those are single-seed point estimates. Across "
              f"**{d['n_seeds']} genuinely disjoint seeds**, spaced by "
              f"{d['seed_stride']:,} so that no two runs share a cell:", ""]
        L += ["| treatment | best-case range | worst-case range |",
              "|---|--:|--:|"]
        for tt in TREATMENTS:
            if tt not in sp:
                continue
            s = sp[tt]
            L.append(f"| {tt} | {s['best_lo']:.2f}x - {s['best_hi']:.2f}x | "
                     f"{s['worst_lo']:.2f}x - {s['worst_hi']:.2f}x |")
        L += [""]
        L += [f"**Spacing matters and nothing said so.** `sim_batch` draws "
              f"cell *i* from `seed + 2i` and its simulation RNG from "
              f"`seed + 2i + 1`, so one run consumes the whole span "
              f"`seed .. seed + {d['seed_stride'] - 1:,}` -- EVEN offsets for "
              f"the cells, ODD offsets for the simulation.", ""]
        L += [f"That makes the overlap depend on the PARITY of the gap, which "
              f"an earlier version of this paragraph got wrong by saying any "
              f"gap under {d['seed_stride']:,} shares almost every cell. An "
              f"EVEN gap *d* shares "
              f"{d['n_per_cell_type']:,} - d/2 cells, so a gap of 2 differs in "
              f"one cell and often returns the same count. An ODD gap shares "
              f"**none** -- one run's cell seeds land on the other's "
              f"simulation seeds -- so it is already an independent sample. "
              f"The safe rule is unchanged and is what this table uses: space "
              f"runs by {d['seed_stride']:,}. What is corrected is the reason, "
              f"and the claim that a nudged seed necessarily returns a "
              f"bit-identical count.", ""]
        # DERIVED. An earlier version wrote "the point estimates sit near the
        # bottom of both ranges" as prose, which survived being flipped to
        # "the top", and computed it against a sample the point was a member
        # of.
        for tt in TREATMENTS:
            if tt not in sp:
                continue
            pt = d["treatments"][tt]["_contrast"].get("ratio_best_case")
            s0 = sp[tt]
            span = s0["best_hi"] - s0["best_lo"]
            if pt is None or span <= 0:
                continue
            frac = (pt - s0["best_lo"]) / span
            where = "bottom" if frac < 0.5 else "top"
            L += [f"The single-seed point for {tt} ({pt:.2f}x) sits at "
                  f"{100*frac:.0f}% of that range, i.e. near the {where} of "
                  f"it. The point's own seed is excluded from the spread, so "
                  f"the two are independent -- an earlier version drew the "
                  f"spread from a set containing the point.", ""]
            break

    ctrl = d["treatments"]["Control"]["_contrast"]
    L += [f"Control is the baseline: with no treatment the CAF phenotype dies at "
          f"{100*ctrl['caf']:.2f}% and tumour phenotypes at "
          f"{100*ctrl['tumour_min']:.2f}-{100*ctrl['tumour_max']:.2f}%. A "
          f"contrast that does not exceed this is not about treatment.", ""]

    # Two properties the table exposes that the project has never reported.
    sdt, pdt = d["treatments"]["SDT"], d["treatments"]["PDT"]
    identical = all(abs(sdt[p]["death_rate"] - pdt[p]["death_rate"]) < 1e-12
                    for p in sdt if not p.startswith("_"))
    rsl3_caf = d["treatments"]["RSL3"][NON_TUMOUR]["death_rate"]

    if identical:
        L += ["## SDT and PDT are the same modality here", ""]
        L += ["Their death rates are **bit-identical across every phenotype**, "
              "because `sdt_ros` and `pdt_ros` share a default value and the "
              "single-cell path differs in nothing else. The two modalities are "
              "distinguished only by their DEPTH PHYSICS, in a different "
              "module.", ""]
        L += ["So any single-cell comparison of SDT against PDT is comparing a "
              "modality with itself. That is not a defect -- it is what the "
              "model says -- but it is worth knowing before reading a "
              "single-cell contrast between them as evidence about two "
              "different therapies.", ""]

    # IS THE ZERO UNIQUE TO THE NON-TUMOUR PHENOTYPE? Derived, because the
    # earlier version read it as a fingerprint of `Stromal` having been
    # parameterised to produce it -- and a control printed in the same row of
    # the same table refutes that: RSL3 also kills exactly zero of a TUMOUR
    # phenotype. The exactness is a property of the RSL3 path against a
    # resistant parameterisation, not evidence about how Stromal was chosen.
    tumour_zeros = sorted(
        ph for ph, row in d["treatments"]["RSL3"].items()
        if not ph.startswith("_") and ph != NON_TUMOUR
        and row["death_rate"] == 0.0)

    if rsl3_caf == 0.0:
        L += ["## The zero denominator, and what it does and does not show", ""]
        L += [f"RSL3 kills **exactly zero** of {d['n_per_cell_type']:,} "
              f"non-tumour cells. Not a small number -- zero, with the interval "
              f"running to "
              f"{100*d['treatments']['RSL3'][NON_TUMOUR]['ci_high']:.3f}%.", ""]
        if tumour_zeros:
            L += [f"**An earlier version of this page read that as a "
                  f"fingerprint** -- \"what a parameter set chosen to produce "
                  f"it looks like\". A control in the same row refutes it: "
                  f"RSL3 also kills exactly zero "
                  f"{', '.join(f'`{x}`' for x in tumour_zeros)} cells, and "
                  f"{'that is a TUMOUR phenotype' if len(tumour_zeros) == 1 else 'those are TUMOUR phenotypes'}. "
                  f"The exactness is a property of the RSL3 path against a "
                  f"resistant parameterisation, not evidence about how "
                  f"`Stromal` was chosen. That inference is withdrawn.", ""]
        rr = d.get("stromal_resistance") or {}
        if rr.get("axes"):
            won = [k for k, v in rr["detail"].items() if v["most_resistant"]]
            L += [f"**The surviving criticism, now measured rather than "
                  f"asserted.** Parsed from `cell.rs`, `{NON_TUMOUR}` is the "
                  f"most ferroptosis-resistant phenotype on "
                  f"**{rr['wins']} of {rr['axes']}** parameter axes "
                  f"({', '.join(f'`{x}`' for x in sorted(won))}). So the "
                  f"resistance really is encoded in the parameters -- that "
                  f"part of the original criticism stands, and only the "
                  f"inference FROM THE EXACT ZERO is withdrawn.", ""]
        L += ["What survives, and it is the part that matters: a ratio with a "
              "zero denominator is undefined, so no selectivity figure can be "
              "computed here at all. And the project's load-bearing assumption "
              "-- that normal cells resist ferroptosis inducers -- is encoded "
              "in the `Stromal` parameters, so a ratio computed against them "
              "would restate the assumption rather than test it. That does not "
              "need the fingerprint argument, and is why it was never load "
              "bearing.", ""]
        L += ["This is the concrete reason the missing normal-tissue phenotype "
              "matters. Until one exists whose parameters come from the "
              "literature rather than from the claim, the model cannot "
              "disagree with the project about selectivity.", ""]

    L += ["## This is not a therapeutic index", ""]
    L += ["`Stromal` models **cancer-associated fibroblasts** -- a "
          "tumour-resident, metabolically altered cell type recruited by the "
          "tumour, whose parameters were chosen to model shielding. It is not "
          "healthy tissue at a distant site, and dividing by it does not give a "
          "therapeutic index. An earlier version of #728 proposed exactly that, "
          "and it is a category error.", ""]
    L += ["A real therapeutic index needs a normal-tissue phenotype distinct "
          "from `Stromal`, with a parameter set traceable to the "
          "ferroptosis-resistance literature (ACSL4-low, MUFA-high, GSH-high). "
          "No such phenotype exists.", ""]
    L += ["That absence is the finding this table exists to make concrete: the "
          "project's most load-bearing selectivity assumption -- that normal "
          "cells resist ferroptosis inducers -- cannot currently be checked "
          "inside the model that assumes it.", ""]

    L += ["## Other limits", ""]
    L += ["* Uncalibrated. These are the engine's default parameters; the "
          "CALIBRATION_STATUS accounting applies unchanged.",
          "* Single-cell, 2D context. The spatial binaries add shielding "
          "geometry that changes the contrast, and are not run here.",
          "* Wilson intervals are in the JSON. A ratio of two rates with "
          "overlapping intervals is not a precise quantity, and the "
          "best/worst-case pair is given instead of a single figure for that "
          "reason.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = run(args.n)
        ctrl = d["treatments"]["Control"]["_contrast"]
        if ctrl["tumour_max"] > 0.5:
            raise SystemExit(
                f"untreated tumour death is {100*ctrl['tumour_max']:.1f}%, which "
                "is not a baseline -- every contrast below it would be measuring "
                "a population that dies on its own.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for t in TREATMENTS:
        c = d["treatments"][t]["_contrast"]
        print(f"  {t:8s} tumour {100*c['tumour_min']:6.2f}-{100*c['tumour_max']:6.2f}%  "
              f"CAF {100*c['caf']:6.2f}%")


if __name__ == "__main__":
    main()
