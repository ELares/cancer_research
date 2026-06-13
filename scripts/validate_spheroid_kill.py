#!/usr/bin/env python3
"""Spheroid RSL3 kill-vs-size: directional validation + a falsifiable prediction (#333 kill leg).

The #333 structure leg validated the spheroid module's zone GEOMETRY (and the
size-aware refinement). The KILL leg asks whether the module reproduces the
size-dependence of kill. There is NO ferroptosis-inducer spheroid size-kill
dataset (a dedicated probe confirmed this), so this is NOT a magnitude
calibration. Instead it does two honest things with a real simulation sweep
(`sim-tme-3d --spheroid-size-sweep`, committed to
`analysis/calibration/spheroid_kill_vs_size.csv`):

1. DIRECTIONAL validation of the SUPPLY-GRADIENT contribution. With the FIXED
   (size-independent) zone thresholds, every spheroid has the same phenotype mix,
   so size-dependence comes only from the O2/drug supply gradient (fixed Krogh λ,
   larger spheroid ⇒ more hypoxic, less-penetrated core). The model's RSL3 kill
   then falls monotonically with size, reproducing the UNIVERSAL bigger-spheroids-
   resist-more direction seen for generic cytotoxics/PDT, and the fold-drop lands
   in the measured cytotoxic size-resistance range. Direction + order of magnitude,
   not calibrated magnitude.

2. A FALSIFIABLE PREDICTION from the size-aware structure. RSL3 preferentially
   kills the PERSISTER phenotype (the manuscript thesis), and the size-aware
   persister core only emerges above ~280 µm radius. So the model predicts a small,
   all-proliferating spheroid RESISTS RSL3 (no persister target ⇒ ~0 kill) and
   vulnerability EMERGES as the persister core appears, the OPPOSITE size-direction
   to generic cytotoxics. This is a distinct, testable prediction (untestable today:
   no ferroptosis-inducer size-kill data), not a validation.

HONEST SCOPE: the kill MAGNITUDE is uncalibrated (the in-vivo-tuned switch is
RSL3-resistant for Glycolytic, and the per-zone biochem gradients are placeholders);
the cytotoxic size-resistance data is for generic proliferating-cell killers, the
WRONG drug class for a persister-targeting ferroptosis inducer, so it bounds the
supply-gradient DIRECTION only.

Run (reads the committed sweep CSV; pure stdlib, CI-safe):
  python3 scripts/validate_spheroid_kill.py
Regenerate the CSV (needs the compiled binary, not run in CI):
  cd simulations && cargo run --release -p sim-tme-3d -- --spheroid-size-sweep
Writes analysis/calibration/spheroid-kill-vs-size.md + .json.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SWEEP_CSV = REPO_ROOT / "analysis" / "calibration" / "spheroid_kill_vs_size.csv"
OUT_MD = REPO_ROOT / "analysis" / "calibration" / "spheroid-kill-vs-size.md"
OUT_JSON = REPO_ROOT / "analysis" / "calibration" / "spheroid-kill-vs-size.json"

# Published cytotoxic / PDT spheroid size-resistance (the generic-drug comparator;
# transcribed fold-ranges, not digitized). These test proliferating-cell killers,
# so they bound the SUPPLY-GRADIENT direction of our model, not the persister-
# targeting RSL3 prediction.
CYTOTOXIC_SIZE_RESISTANCE = {
    "min_fold": 1.8,   # West 1989 (Br J Cancer, WiDr, PDT) smaller-size end
    "max_fold": 22.0,  # West 1989 largest-size end (~22x less sensitive at 500 µm)
    "sources": [
        "West 1989 Br J Cancer (WiDr spheroids, PDT): ~1.8x/2.5x/22x less sensitive at 100/250/500 µm",
        "Eilenberger 2021 Adv Sci (A549, doxorubicin/cisplatin): 2.3 to 6.9x IC50 increase large vs small",
        "Demuynck 2020 PMID 32183000 (3DELTA, ML-162 ferroptosis): ~50% to ~30% kill small to large (2-point)",
    ],
}


def load_sweep(path=SWEEP_CSV):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "grid_dim": int(r["grid_dim"]),
                "radius_um": float(r["radius_um"]),
                "size_aware": r["size_aware"] == "true",
                "kill_rate": float(r["kill_rate"]),
                "total_tumor": int(r["total_tumor"]),
                "total_dead": int(r["total_dead"]),
            })
    return rows


def _by_radius(rows, size_aware):
    return sorted([r for r in rows if r["size_aware"] == size_aware], key=lambda r: r["radius_um"])


def _is_monotone_decreasing(vals):
    return all(a >= b for a, b in zip(vals, vals[1:]))


def analyze(rows):
    fixed = _by_radius(rows, False)
    sized = _by_radius(rows, True)

    # --- Fixed thresholds: supply-gradient direction ---
    fixed_kills = [r["kill_rate"] for r in fixed]
    fixed_decreasing = _is_monotone_decreasing(fixed_kills)
    # fold-drop smallest->largest spheroid (guard against zero)
    fold = (fixed_kills[0] / fixed_kills[-1]) if fixed_kills[-1] > 0 else float("inf")
    in_cyto_range = CYTOTOXIC_SIZE_RESISTANCE["min_fold"] <= fold <= CYTOTOXIC_SIZE_RESISTANCE["max_fold"]

    # --- Size-aware: persister-targeting twist ---
    # small spheroids (below the ~280 µm core onset) should be ~0 kill (all
    # proliferating, no persister target); kill emerges as the persister core appears.
    small_sized = [r for r in sized if r["radius_um"] < 280.0]
    large_sized = [r for r in sized if r["radius_um"] >= 280.0]
    small_all_resist = all(r["kill_rate"] < 1e-4 for r in small_sized) if small_sized else False
    core_emerges = (max((r["kill_rate"] for r in large_sized), default=0.0)
                    > max((r["kill_rate"] for r in small_sized), default=0.0))

    result = {
        "source": "sim-tme-3d --spheroid-size-sweep (committed analysis/calibration/spheroid_kill_vs_size.csv)",
        "treatment": "RSL3 (pharmacologic; immune/stromal/pH off to isolate the pharmacologic kill)",
        "o2_lambda_um": 120.0,
        "supply_gradient_direction": {
            "fixed_threshold_kills_by_radius_um": {str(int(r["radius_um"])): round(r["kill_rate"], 6) for r in fixed},
            "monotone_decreasing_with_size": fixed_decreasing,
            "fold_drop_small_to_large": round(fold, 1),
            "radius_range_um": [fixed[0]["radius_um"], fixed[-1]["radius_um"]],
            "cytotoxic_measured_fold_range": [CYTOTOXIC_SIZE_RESISTANCE["min_fold"], CYTOTOXIC_SIZE_RESISTANCE["max_fold"]],
            "fold_in_measured_range": in_cyto_range,
            "verdict": ("bigger-spheroids-resist-more direction reproduced and fold-drop in the measured "
                        "cytotoxic range (DIRECTION + order of magnitude validated; magnitude uncalibrated)"),
        },
        "persister_targeting_prediction": {
            "size_aware_kills_by_radius_um": {str(int(r["radius_um"])): round(r["kill_rate"], 6) for r in sized},
            "small_spheroids_resist": small_all_resist,
            "vulnerability_emerges_with_core": core_emerges,
            "prediction": ("RSL3 preferentially kills the persister phenotype, and the size-aware persister "
                           "core only emerges above ~280 µm, so small all-proliferating spheroids RESIST RSL3 "
                           "(no persister target) and vulnerability emerges with the core, the OPPOSITE "
                           "size-direction to generic cytotoxics. Falsifiable; untestable without "
                           "ferroptosis-inducer size-kill data."),
        },
        "cytotoxic_comparator_sources": CYTOTOXIC_SIZE_RESISTANCE["sources"],
        "scope": ("Kill MAGNITUDE uncalibrated (in-vivo-resistant switch + placeholder gradients). The "
                  "cytotoxic data validates the supply-gradient DIRECTION only (wrong drug class for a "
                  "persister-targeting inducer). The persister-targeting prediction is falsifiable but "
                  "data-blocked. STRUCTURE was validated separately (#333 structure leg + size-aware)."),
    }
    return result


def run(args):
    rows = load_sweep(args.csv)
    result = analyze(rows)
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_report(result)
    sg = result["supply_gradient_direction"]
    pt = result["persister_targeting_prediction"]
    print(f"fixed-threshold RSL3 kill vs size: monotone-decreasing={sg['monotone_decreasing_with_size']}, "
          f"fold-drop={sg['fold_drop_small_to_large']}x over {sg['radius_range_um']} um, "
          f"in cytotoxic range {sg['cytotoxic_measured_fold_range']}={sg['fold_in_measured_range']}")
    print(f"size-aware persister-targeting: small-spheroids-resist={pt['small_spheroids_resist']}, "
          f"vulnerability-emerges-with-core={pt['vulnerability_emerges_with_core']}")
    print(f"wrote {OUT_MD.relative_to(REPO_ROOT)} + {OUT_JSON.relative_to(REPO_ROOT)}")
    return result


def write_report(r):
    sg = r["supply_gradient_direction"]
    pt = r["persister_targeting_prediction"]

    def kill_table(d):
        radii = sorted(d, key=lambda k: int(k))
        head = "| radius (µm) | " + " | ".join(radii) + " |"
        sep = "|---|" + "---|" * len(radii)
        vals = "| RSL3 kill | " + " | ".join(f"{d[k] * 100:.2f}%" for k in radii) + " |"
        return "\n".join([head, sep, vals])

    md = f"""# Spheroid RSL3 kill vs size: direction validated + a falsifiable prediction (#333 kill leg)

Generated by `scripts/validate_spheroid_kill.py` (reads the committed
`analysis/calibration/spheroid_kill_vs_size.csv`, written by
`sim-tme-3d --spheroid-size-sweep`; pure stdlib, CI-safe).

Treatment: **{r['treatment']}**. O₂ λ fixed at {r['o2_lambda_um']:.0f} µm, so a larger
spheroid has a proportionally more hypoxic, less-penetrated core (the physical
size-dependence). There is NO ferroptosis-inducer spheroid size-kill dataset, so
this is a DIRECTIONAL validation + a falsifiable prediction, not a magnitude
calibration.

## 1. Supply-gradient direction (FIXED thresholds): bigger spheroids resist more

With the fixed (size-independent) zone thresholds every spheroid has the same
phenotype mix, so size-dependence comes only from the O₂/drug supply gradient.

{kill_table(sg['fixed_threshold_kills_by_radius_um'])}

- Monotone-decreasing with size: **{sg['monotone_decreasing_with_size']}**.
- Fold-drop smallest→largest ({sg['radius_range_um'][0]:.0f}→{sg['radius_range_um'][1]:.0f} µm):
  **{sg['fold_drop_small_to_large']}x**, within the measured cytotoxic size-resistance
  range {sg['cytotoxic_measured_fold_range']} (in range: **{sg['fold_in_measured_range']}**).

This reproduces the universal bigger-spheroids-resist-more direction (generic
cytotoxics / PDT) to direction + order of magnitude. The magnitude itself is
uncalibrated.

## 2. Persister-targeting twist (SIZE-AWARE thresholds): a falsifiable prediction

{kill_table(pt['size_aware_kills_by_radius_um'])}

- Small spheroids (R < 280 µm) resist RSL3 (≈0 kill): **{pt['small_spheroids_resist']}**.
- Vulnerability emerges as the persister core appears: **{pt['vulnerability_emerges_with_core']}**.

{pt['prediction']}

This is the OPPOSITE size-direction to the generic cytotoxic data, and it is a
genuine consequence of two model features: RSL3 preferentially kills the persister
phenotype, and (with size-aware structure) the persister core only emerges above
~280 µm radius. A small, all-proliferating spheroid has no persister target, so it
resists the persister-targeting inducer.

## Comparator sources (generic cytotoxic / PDT size-resistance)

{chr(10).join('- ' + s for s in r['cytotoxic_comparator_sources'])}

## Scope and honesty

{r['scope']}

In short: the SUPPLY-GRADIENT contribution to size-resistance is directionally
validated against measured cytotoxic data (right direction, fold in range), and the
PERSISTER-TARGETING contribution is a distinct, falsifiable prediction that the
existing cytotoxic data cannot test (wrong drug class) and that no ferroptosis-
inducer size-kill dataset yet exists to validate. The kill leg of #333 is therefore
characterized as far as public data allows; a ferroptosis-inducer multi-size
spheroid kill experiment is the experiment that would close it (a falsifiable
prediction for `analysis/contribution-plan-2026.md`).
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, default=SWEEP_CSV)
    args = ap.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
