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


def run(n: int) -> dict:
    import ferroptosis_core as fc
    out = {"n_per_cell_type": n, "seed": SEED, "seed_stride": 2 * n,
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
            s = SEED + k * 2 * n
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
              f"cell *i* from `seed + 2i`, so one run consumes the whole span "
              f"`seed .. seed + {d['seed_stride'] - 1:,}`. Two runs whose "
              f"seeds differ by less than {d['seed_stride']:,} share almost "
              f"every cell -- `seed` and `seed + 2` differ by one cell in "
              f"{d['n_per_cell_type']:,} and return bit-identical counts -- so "
              f"anyone checking robustness by nudging the seed gets a false "
              f"confirmation. The point estimates above sit near the bottom of "
              f"both ranges.", ""]

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
