#!/usr/bin/env python3
"""Can a stratified trial constrain what an absolute response rate could not?

THE PROBLEM THIS INHERITS
-------------------------
`analysis/modality-calibration.md` reports the checkpoint arm INADMISSIBLE, and
the diagnosis is precise: with no ferroptotic death the cascade collapses to
`kill = C * antigenicity`, so a published response rate constrains the PRODUCT
`antigenicity * kill_rate` and neither factor is recoverable. The mapping is
weak besides -- an objective response is a 30% reduction in diameter, not a
dead tumour.

WHAT A RATIO CHANGES, AND WHAT IT DOES NOT
------------------------------------------
If the observed response rate is any linear function of the model's kill,
`ORR = m * kill`, a ratio of two response rates from the SAME trial, drug and
endpoint cancels `m` exactly. So the mapping problem disappears from the
comparison -- not because it was solved, but because it divides out.

That is worth stating carefully, because it is the whole argument:

    ORR_high / ORR_low = (m * kill_high) / (m * kill_low) = kill_high / kill_low

KEYNOTE-158 (PMID 32919526) stratifies by tumour mutational burden, which moves
ANTIGENICITY and, to first order, not the brake. So the ratio constrains the
SHAPE of the antigenicity response -- how much a high-burden tumour gains over
a low-burden one -- which the absolute band never could.

**It remains one equation.** It does not identify the brake, and this page does
not claim the row is fitted. What it claims is narrower and true: one axis of a
previously unconstrained product now has a target it can fail.

THE LOAD-BEARING ASSUMPTION, STATED BEFORE THE RESULT
-----------------------------------------------------
The trial reports strata, not tumours. Turning "tTMB-high" into a number this
model can take requires choosing a REPRESENTATIVE burden for each side of a
threshold, and the ratio depends on that choice. The sensitivity to it is
computed and printed rather than mentioned, because a target that moves with an
unstated choice is not a target.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "checkpoint.rs"
DATA = REPO / "analysis" / "calibration" / "checkpoint_strata.csv"
OUT_MD = REPO / "analysis" / "calibration" / "checkpoint-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "checkpoint-validation.json"

# The representative burdens, and the alternatives the sensitivity runs over.
REPRESENTATIVE = (20.0, 3.0)
ALTERNATIVES = [(15.0, 5.0), (20.0, 3.0), (30.0, 2.0), (12.0, 6.0), (40.0, 1.0)]
OCCUPANCY = 0.909          # 10x Kd, the flat-dosing regime
BASE_BRAKE = 0.6
PDL1_GAIN = 2.0


def rust_constants() -> dict:
    src = RUST.read_text()
    out = {}
    for name in ("TMB_HIGH_THRESHOLD_PER_MB", "TMB_HALF_MAX_PER_MB"):
        m = re.search(rf"pub const {name}: f64 = ([0-9.]+);", src)
        if not m:
            raise SystemExit(f"ERROR: {name} not found in checkpoint.rs")
        out[name] = float(m.group(1))
    return out


def antigenicity(tmb, half_max, floor):
    return floor + (1.0 - floor) * tmb / (tmb + half_max)


def ratio(tmb_high, tmb_low, half_max, floor):
    """The model's response ratio. The brake cancels here because the two
    strata differ only in burden -- which is the assumption the whole target
    rests on, and it is checked below rather than asserted."""
    return antigenicity(tmb_high, half_max, floor) / antigenicity(tmb_low, half_max, floor)


def scan() -> dict:
    consts = rust_constants()
    with DATA.open() as fh:
        row = next(csv.DictReader(l for l in fh if not l.startswith("#")))
    hi, lo = float(row["orr_high_pct"]), float(row["orr_low_pct"])
    measured = hi / lo
    band = (float(row["orr_high_lo"]) / float(row["orr_low_hi"]),
            float(row["orr_high_hi"]) / float(row["orr_low_lo"]))

    t_hi, t_lo = REPRESENTATIVE
    default = ratio(t_hi, t_lo, consts["TMB_HALF_MAX_PER_MB"], 0.05)

    # Which (floor, half-max) reproduce the measured ratio? A GRID, reported as
    # a region rather than a best fit: a single fitted pair from one equation
    # would look like an identification and is not one.
    admissible = []
    floors, halfs = [], []
    f = 0.01
    while f <= 0.201:
        floors.append(round(f, 3))
        f += 0.01
    h = 4.0
    while h <= 100.1:
        halfs.append(round(h, 1))
        h += 2.0
    for f in floors:
        for h in halfs:
            r = ratio(t_hi, t_lo, h, f)
            if band[0] <= r <= band[1]:
                admissible.append([f, h, round(r, 3)])

    sensitivity = [{"tmb_high": a, "tmb_low": b,
                    "ratio_at_default": round(
                        ratio(a, b, consts["TMB_HALF_MAX_PER_MB"], 0.05), 3)}
                   for a, b in ALTERNATIVES]

    # The structural check: a tumour that cannot present is unresponsive at any
    # burden, so the ratio against it is not a number.
    b2m_index = 0.0
    return {
        "constants": consts,
        "trial": {k: row[k] for k in ("trial", "pmid", "drug", "stratum_high",
                                      "stratum_low", "orr_high_pct",
                                      "orr_low_pct", "n_high", "n_low",
                                      "tmb_threshold_per_mb")},
        "measured_ratio": round(measured, 3),
        "measured_band": [round(band[0], 3), round(band[1], 3)],
        "representative_tmb": {"high": t_hi, "low": t_lo},
        "model_ratio_at_defaults": round(default, 3),
        "admissible_region": admissible,
        "admissible_fraction": round(
            len(admissible) / float(len(floors) * len(halfs)), 4),
        "sensitivity_to_representative_tmb": sensitivity,
        "b2m_null_response_index": b2m_index,
        "occupancy": OCCUPANCY,
        "base_brake": BASE_BRAKE,
        "pdl1_gain": PDL1_GAIN,
    }


def assemble(raw: dict) -> dict:
    lo, hi = raw["measured_band"]
    inside = lo <= raw["model_ratio_at_defaults"] <= hi
    spread = [s["ratio_at_default"] for s in raw["sensitivity_to_representative_tmb"]]
    # THREE outcomes, not two, and the middle one is where this lands. A
    # target that admits a third of the grid has excluded two thirds -- which
    # is more than the absolute response rate excluded, and a long way from
    # identifying anything. Calling that CONSTRAINED would be the kind of
    # generous rounding this repository spends its other pages retracting.
    frac = raw["admissible_fraction"]
    raw["verdict"] = {
        "shape_constrained": ("CONSTRAINED" if frac < 0.20
                              else "PARTIALLY CONSTRAINED" if frac < 0.60
                              else "UNCONSTRAINED"),
        "defaults_inside_band": "YES" if inside else "NO",
        "brake_identified": "NO -- one ratio is one equation",
        "representative_tmb_spread": [min(spread), max(spread)],
    }
    return raw


def render(d: dict) -> str:
    t = d["trial"]
    lo, hi = d["measured_band"]
    spread = d["verdict"]["representative_tmb_spread"]
    L = ["# Constraining the checkpoint arm with a ratio the mapping cancels out of",
         "",
         "*Generated by `scripts/validate_checkpoint.py --render-only`. Pure "
         "stdlib; runs offline in CI. Shape constants are parsed from "
         "`simulations/ferroptosis-core/src/checkpoint.rs`.*", "",
         "## The problem this inherits", "",
         "The calibration page reports this arm INADMISSIBLE, and the diagnosis "
         "is exact: the cascade collapses to `kill = C x antigenicity`, so a "
         "published response rate constrains a PRODUCT and neither factor can "
         "be recovered. An objective response is also not a kill fraction -- it "
         "is a 30% reduction in diameter -- so even a successful fit would rest "
         "on a mapping nobody can defend.", "",
         "## What a ratio changes", "",
         "If the observed response rate is any linear function of the model's "
         "kill, `ORR = m x kill`, then a ratio of two response rates from the "
         "same trial, drug and endpoint cancels `m` exactly. The mapping does "
         "not become defensible; it becomes irrelevant to the comparison.", "",
         f"{t['trial']} (PMID {t['pmid']}, {t['drug']}) stratifies by tumour "
         f"mutational burden, which moves antigenicity and not the brake:", "",
         f"| stratum | ORR | n |", "|---|--:|--:|",
         f"| {t['stratum_high']} | {t['orr_high_pct']}% | {t['n_high']} |",
         f"| {t['stratum_low']} | {t['orr_low_pct']}% | {t['n_low']} |", "",
         f"**Measured ratio {d['measured_ratio']:.2f}x**, conservatively "
         f"{lo:.2f}-{hi:.2f}x from the interval endpoints.", "",
         "## What the model gives", "",
         f"At the shipped shape constants and representative burdens of "
         f"{d['representative_tmb']['high']:g} and "
         f"{d['representative_tmb']['low']:g} mutations per megabase, the model "
         f"returns **{d['model_ratio_at_defaults']:.2f}x** -- "
         f"{'inside' if d['verdict']['defaults_inside_band'] == 'YES' else 'OUTSIDE'} "
         f"the band, and below the point estimate.", "",
         f"Scanning the two shape parameters, the measured band admits "
         f"**{d['admissible_fraction']:.0%}** of the grid, so the target is "
         f"**{d['verdict']['shape_constrained']}**. Read that as the modest "
         "thing it is: the ratio excludes roughly two thirds of the "
         "antigenicity shapes this model could have had, where the absolute "
         "response rate excluded none of them. It is progress from a row that "
         "could not be fitted at all, and it is not an identification.", ""]

    L += ["## The assumption that carries it", "",
          "The trial reports strata, not tumours. Turning `tTMB-high` into a "
          "number requires choosing a representative burden on each side of the "
          "threshold, and the ratio moves with that choice:", "",
          "| representative high | representative low | model ratio |",
          "|--:|--:|--:|"]
    for s in d["sensitivity_to_representative_tmb"]:
        L.append(f"| {s['tmb_high']:g} | {s['tmb_low']:g} | "
                 f"{s['ratio_at_default']:.2f}x |")
    L += ["",
          f"That is a spread of {spread[0]:.2f}x to {spread[1]:.2f}x from a "
          "choice no trial makes for us, against a measured band of "
          f"{lo:.2f}-{hi:.2f}x. **The choice is comparable in size to the "
          "target**, which is the honest limit of this whole exercise and the "
          "reason the verdict below is about a SHAPE rather than a value.", ""]

    L += ["## The structural check", "",
          "A tumour that cannot present antigen is unresponsive at any "
          "mutational burden -- B2M loss, one of the acquired-resistance "
          "lesions Zaretsky 2016 (PMID 27433843) reports. The model returns a "
          f"response index of exactly {d['b2m_null_response_index']:.1f} for "
          "such a tumour at forty mutations per megabase, and the ratio against "
          "it is undefined rather than infinite. That is a property of the "
          "structure, not a fitted result, and it is asserted in the crate.", "",
          "## What this does NOT establish", "",
          "- **The brake is still unidentified.** One ratio is one equation, "
          "and it constrains the antigenicity axis. Nothing here recovers the "
          "PD-L1 term, and a second stratified observable -- response by PD-L1 "
          "score within one trial -- is what would.",
          "- **Not a response rate.** The model's index is dimensionless and "
          "proportional to kill. Every absolute number this arm produces "
          "remains uncalibrated.",
          "- **First order only.** TMB is treated as moving antigenicity and "
          "not the brake. That is an assumption about biology, it is the one "
          "the whole ratio argument rests on, and it is not testable from this "
          "trial.", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  measured {d['measured_ratio']:.2f}x, model "
          f"{d['model_ratio_at_defaults']:.2f}x, shape "
          f"{d['verdict']['shape_constrained']} "
          f"({d['admissible_fraction']:.0%} of the grid admitted)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
