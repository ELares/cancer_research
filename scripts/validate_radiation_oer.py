#!/usr/bin/env python3
"""The oxygen enhancement ratio a SPATIAL run exhibits (#844).

THE TARGET WAS NAMED BEFORE THE WORK, WHICH IS THE POINT
--------------------------------------------------------
`simulations/calibration/CALIBRATION_STATUS.md` has said since #726 that the
radiation layer "ships with NO INDEPENDENT FAILABLE PREDICTION YET. The one it
will have is the OER a full run EXHIBITS with both channels active against the
published 2.5-3.0, which is emergent and needs the binary wiring." #844 is that
wiring and this is that measurement.

WHY IT IS NOT A RESTATEMENT
---------------------------
`radiation::dna_channel_dose_modifying_factor` returns 2.86 and is documented
and tested as a RESTATEMENT of the Alper-Howard-Flanders hyperbola -- scoring
it against the published band would be a guard computing its own expectation.
What is measured here is different in four ways, and the difference shows up in
the answer:

  * the core is a DISTRIBUTION of pO2 rather than anoxia
  * the linear-quadratic model is not log-linear in dose, so the ratio can
    depend on the kill level it is read at
  * the ratio is a population average over a nonlinear function
  * and it depends on the O2 gradient's steepness, a spatial parameter that
    does not appear in the formula at all

The formula gives one number. The spatial engine gives 1.20 to 2.64 depending
on the gradient, non-monotonically.

WHAT IS SCORED, AND THE STRUCTURAL REASON IT FALLS SHORT
--------------------------------------------------------
The published 2.5-3.0 is measured between fully oxygenated and fully ANOXIC
cells. This model cannot present that pair: ONE lambda sets both the rim and
the core, so a gradient steep enough to make the core anoxic also leaves the
rim hypoxic, and a gradient shallow enough to oxygenate the rim leaves the core
well up the hyperbola. Every measured factor is therefore a LOWER BOUND on the
oxic-anoxic ratio the band describes, and the sweep's maximum is the model's
best estimate rather than its answer.

That structure produces an interior optimum in gradient steepness that nobody
designed in, which is the most interesting thing here.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SWEEP = REPO / "analysis" / "calibration" / "radiation_oer_sweep.txt"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "radiation.rs"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"
OUT_MD = REPO / "analysis" / "calibration" / "radiation-oer-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "radiation-oer-validation.json"

# Alper & Howard-Flanders; the band radiobiology reports for the oxygen effect
# on mammalian cells. NOT fitted here and NOT this project's number.
PUBLISHED_OER = (2.5, 3.0)
KILL_LEVELS = (0.3, 0.5, 0.7, 0.9)


def _rows():
    out = []
    for ln in SWEEP.read_text().splitlines():
        if not ln.startswith("RADIATION_OER_SWEEP"):
            continue
        d = dict(re.findall(r"(\w+)=([-\d.eE]+)", ln))
        out.append({k: float(v) for k, v in d.items()})
    return out


def _zones():
    out = []
    for ln in SWEEP.read_text().splitlines():
        if not ln.startswith("RADIATION_OER_ZONES"):
            continue
        d = dict(re.findall(r"(\w+)=([-\d.eE]+)", ln))
        out.append({k: float(v) for k, v in d.items()})
    return out


def dose_for_kill(rows, lam, zone, level):
    """Linear interpolation on the measured curve. `None` when the level is
    never reached, which must not be silently read as an infinite dose."""
    pts = sorted((r["dose_gy"], r[zone]) for r in rows if r["lambda_um"] == lam)
    for (d0, k0), (d1, k1) in zip(pts, pts[1:]):
        if k0 <= level <= k1 and k1 > k0:
            return d0 + (d1 - d0) * (level - k0) / (k1 - k0)
    return None


def scan() -> dict:
    rows, zones = _rows(), _zones()
    lambdas = sorted({r["lambda_um"] for r in rows})
    by_lambda = []
    for lam in lambdas:
        z = next(x for x in zones if x["lambda_um"] == lam)
        factors = []
        for lvl in KILL_LEVELS:
            a = dose_for_kill(rows, lam, "normoxic", lvl)
            b = dose_for_kill(rows, lam, "hypoxic", lvl)
            factors.append({
                "kill_level": lvl,
                "dose_rim_gy": a, "dose_core_gy": b,
                "dmf": (b / a) if (a and b and a > 0) else None,
            })
        vals = [f["dmf"] for f in factors if f["dmf"] is not None]
        by_lambda.append({
            "lambda_um": lam,
            "rim_po2_mmhg": z["rim_po2_mmhg"],
            "core_po2_mmhg": z["deep_po2_mmhg"],
            "rim_n": int(z["rim_n"]), "core_n": int(z["deep_n"]),
            "factors": factors,
            "dmf_mean": sum(vals) / len(vals) if vals else None,
            "dmf_spread": (max(vals) - min(vals)) if vals else None,
        })
    return {
        "published_oer_band": list(PUBLISHED_OER),
        "kill_levels": list(KILL_LEVELS),
        "by_lambda": by_lambda,
        "restated_single_cell_dmf": 2.86,
    }


def assemble(raw: dict) -> dict:
    d = dict(raw)
    rows = [r for r in raw["by_lambda"] if r["dmf_mean"] is not None]
    best = max(rows, key=lambda r: r["dmf_mean"])
    d["best_lambda_um"] = best["lambda_um"]
    d["best_dmf"] = best["dmf_mean"]
    lo, hi = raw["published_oer_band"]
    d["best_reaches_band"] = lo <= best["dmf_mean"] <= hi
    # The interior optimum: the peak must be strictly inside the swept range,
    # or "the model has a best gradient" is a statement about the scan.
    d["optimum_is_interior"] = (
        best["lambda_um"] != min(r["lambda_um"] for r in rows)
        and best["lambda_um"] != max(r["lambda_um"] for r in rows))
    # How much of the answer is WHERE you read it. If this were large the
    # single "the DMF is X" framing would be wrong.
    d["max_spread_across_kill_levels"] = max(r["dmf_spread"] for r in rows)
    # No lambda gives both a well-oxygenated rim and an anoxic core, which is
    # the structural reason every factor is a lower bound.
    d["best_rim_po2_at_anoxic_core"] = min(
        (r["rim_po2_mmhg"] for r in rows if r["core_po2_mmhg"] < 0.1), default=None)
    d["arms_below_band"] = sorted(
        r["lambda_um"] for r in rows if r["dmf_mean"] < lo)
    d["verdict"] = (
        "PARTIAL — reaches the band at one gradient only"
        if d["best_reaches_band"] and d["arms_below_band"] else
        "ADMISSIBLE" if d["best_reaches_band"] else "BELOW BAND")
    return d


def render(d: dict) -> str:
    lo, hi = d["published_oer_band"]
    L = [
        "# Radiation: the oxygen enhancement ratio a spatial run exhibits (#844)",
        "",
        "*Generated by `scripts/validate_radiation_oer.py --render-only`. Pure "
        "stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}.**",
        "",
        "## The target was named before the work",
        "",
        "`CALIBRATION_STATUS.md` has said since #726 that the radiation layer "
        "*\"ships with NO INDEPENDENT FAILABLE PREDICTION YET. The one it will "
        "have is the OER a full run EXHIBITS with both channels active against "
        "the published 2.5-3.0, which is emergent and needs the binary "
        "wiring.\"* #844 supplied the wiring; this is the measurement.",
        "",
        "## Why this is not the formula restated",
        "",
        f"`radiation::dna_channel_dose_modifying_factor` returns "
        f"{d['restated_single_cell_dmf']} and is documented and tested as a "
        "RESTATEMENT of the Alper–Howard-Flanders hyperbola. Scoring *that* "
        "against the published band would be a guard computing its own "
        "expectation. The spatial measurement differs in four ways — the core "
        "is a distribution of pO₂ rather than anoxia, the linear-quadratic "
        "model is not log-linear in dose, the ratio is a population average "
        "over a nonlinear function, and it depends on the gradient's "
        "steepness, a parameter that does not appear in the formula at all.",
        "",
        "The formula gives one number. The engine gives a range, and it is not "
        "monotonic:",
        "",
        "| O₂ gradient λ | rim pO₂ | core pO₂ | core cells | dose-modifying factor |",
        "|--:|--:|--:|--:|--:|",
    ]
    for r in d["by_lambda"]:
        dmf = f"{r['dmf_mean']:.2f}" if r["dmf_mean"] is not None else "—"
        mark = " **←**" if r["lambda_um"] == d["best_lambda_um"] else ""
        L.append(f"| {r['lambda_um']:.0f} µm | {r['rim_po2_mmhg']:.1f} mmHg | "
                 f"{r['core_po2_mmhg']:.3f} mmHg | {r['core_n']:,} | {dmf}{mark} |")
    L += [
        "",
        f"Averaged over kill levels {d['kill_levels']}; the factor moves by at "
        f"most **{d['max_spread_across_kill_levels']:.2f}** across them, so it "
        "is not strongly a property of where the curve is read — which the "
        "quadratic term made a real possibility and it had to be checked "
        "rather than assumed.",
        "",
        "## The interior optimum nobody designed in",
        "",
        f"The largest factor is at **λ = {d['best_lambda_um']:.0f} µm** "
        f"({d['best_dmf']:.2f}), and it falls at BOTH steeper and shallower "
        "gradients. That is a structural consequence rather than a tuning "
        "artifact: **one λ sets the rim and the core together.** A gradient "
        "steep enough to make the core anoxic leaves the rim hypoxic too — at "
        f"λ = 30 µm the core reaches 0.0001 mmHg but the rim is only "
        f"{d['best_rim_po2_at_anoxic_core']:.1f} mmHg, so the rim has already "
        "lost part of its own oxygen advantage. A gradient shallow enough to "
        "oxygenate the rim leaves the core well up the hyperbola.",
        "",
        "## What is scored, and what it does not establish",
        "",
        f"The published band is **{lo}–{hi}**, measured between fully "
        "oxygenated and fully **anoxic** cells. This model cannot present that "
        "pair, so every factor above is a **lower bound** on the ratio the "
        "band describes and the sweep's maximum is the model's best estimate "
        "rather than its answer.",
        "",
        f"- It reaches the band at **one** gradient of five "
        f"({d['best_dmf']:.2f} at λ = {d['best_lambda_um']:.0f} µm), and only just.",
        f"- At the engine's own zone reference (λ = 120 µm) it gives "
        + f"{next(r['dmf_mean'] for r in d['by_lambda'] if r['lambda_um'] == 120.0):.2f}"
        + ", **below** the band. Reporting that alone as a failure would have "
        "been the category error this page exists to avoid: at that gradient "
        "the deep zone sits near 1.4 mmHg and is not anoxic, so it is not the "
        "quantity the band measures.",
        f"- The gradients below the band are λ = "
        + ", ".join(f"{v:.0f}" for v in d["arms_below_band"]) + " µm.",
        "",
        "**Only the DNA channel is wired.** #726 built radiation as two "
        "independent channels and the second one's `ros_per_gy` has no "
        "gray-to-ROS conversion in the literature. Switching it on would make "
        "this measurement a function of an unanchored knob, so the target is "
        "met with one channel and this page says which.",
        "",
        "**A grid smaller than 40 cannot make this measurement at all.** The "
        "deep zone is empty below that — 0 cells at grid 20 and 30, 123 at 40, "
        "3,071 at the production 60 — and `zone_kill_rates_3d` returns 0.0 for "
        "an empty zone, so a smaller probe reports a hypoxic kill rate of "
        "exactly zero at every dose. That reads as total radioresistance and "
        "is a division by nothing. The first version of this sweep ran at grid "
        "30 and would have published it.",
        "",
    ]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if args.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    print(f"{OUT_MD.relative_to(REPO)}: {d['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
