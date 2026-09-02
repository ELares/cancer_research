#!/usr/bin/env python3
"""How fast to deliver the light, and why the answer is not 'as fast as possible'.

THE CLAIM BEING TESTED
----------------------
PDT is dosed as a total fluence in J/cm2, and for most of this engine's life
the fluence RATE did not exist as a quantity at all: `pdt_intensity_at_depth`
samples the photosensitizer once, at the moment illumination begins, and light
is otherwise treated as a scalar. That treats a 10-minute and a 2-hour
illumination delivering the same joules as the same treatment.

They are not, and PDT is the one arm in this chapter where the field says so
loudly. A Type II photosensitizer works by handing energy to ground-state O2,
so illumination CONSUMES ITS OWN SUBSTRATE at the site it is treating.
Henderson & Busch, "Fluence rate as a modulator of PDT mechanisms" (Lasers
Surg Med 2006, PMID 16615136) report oxygen depletion within seconds of
illumination at 75 mW/cm2 in murine tumours, and the fluence-rate effect --
the same total light delivered more slowly producing a greater response -- is
the practical consequence.

WHAT THE LAYER ADDS, AND WHY THE OPTIMUM IS INTERIOR
----------------------------------------------------
Two effects oppose across the choice of rate:

  * FASTER depletes oxygen at the site, cutting the per-photon singlet-oxygen
    yield -- so the same joules do less
  * SLOWER runs the illumination further into the sensitizer's own clearance,
    and drug that has left cannot be excited

The second term is not new machinery. It is `Photosensitizer::yield_at`, the
module's existing pharmacokinetics, INTEGRATED over the illumination window
rather than sampled at its start. Sampling once is precisely the approximation
that makes a long illumination look free.

WHAT IS ANCHORED AND WHAT IS NOT
--------------------------------
The DIRECTION is anchored: an oxygen-blind model has no reason to prefer any
rate, and the measured fluence-rate effect says slower is better over the
clinical range. The MAGNITUDE is not. The optimum's position scales with
`phi_crit`, the fluence rate at which photochemical consumption matches
perfusive resupply, and nothing in this repository measures that quantity --
so the milliwatt-per-square-centimetre figure is a restatement of an
assumption, and this script measures how strongly, rather than hiding it.

THE PREDICTION THAT IS THE MODEL'S OWN
--------------------------------------
The optimal rate falls as the sensitizer's half-life rises. That is a
statement about DRUGS rather than about doses: a fast-clearing sensitizer
should be illuminated fast, a slow-clearing one slowly, and the two should
differ by roughly the ratio of their half-lives raised to a fractional power.
Registered as P21; no dataset here tests it.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "photosensitizer_pk.rs"
OXY = REPO / "simulations" / "ferroptosis-core" / "src" / "oxygen.rs"
OUT_MD = REPO / "analysis" / "calibration" / "pdt-fluence-rate-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "pdt-fluence-rate-validation.json"

TOTAL_FLUENCE_J = 150.0
SCAN_LO, SCAN_HI = 5.0, 400.0
SCAN_N = 4000
PHI_CRIT_DEFAULT = 50.0


def oxygen_enhancement_ratio(p_mmhg):
    p = max(p_mmhg, 0.0)
    return (3.0 * p + 3.0) / (p + 3.0)


def oer_relative_efficacy(o2_supply, p_full):
    s = min(max(o2_supply, 0.0), 1.0)
    return min(max(oxygen_enhancement_ratio(s * p_full) / oxygen_enhancement_ratio(p_full), 0.0), 1.0)


def o2_fraction(rate, phi_crit):
    if rate <= 0:
        return 1.0
    return 0.0 if phi_crit <= 0 else 1.0 / (1.0 + rate / phi_crit)


def delivered(t_half_h, rate, phi_crit, p_full=40.0, steps=64):
    """The crate's `delivered_singlet_oxygen`, re-implemented."""
    if rate <= 0:
        return 0.0
    duration_h = TOTAL_FLUENCE_J * 1000.0 / rate / 3600.0
    g = oer_relative_efficacy(o2_fraction(rate, phi_crit), p_full)
    k = math.log(2.0) / t_half_h
    dt = duration_h / steps
    drug_hours = sum(math.exp(-k * ((i + 0.5) * dt)) * dt for i in range(steps))
    return g * rate * drug_hours * 3.6


def best_rate(t_half_h, phi_crit):
    best = (SCAN_LO, -1.0)
    for i in range(SCAN_N):
        f = i / (SCAN_N - 1)
        r = SCAN_LO * (SCAN_HI / SCAN_LO) ** f
        d = delivered(t_half_h, r, phi_crit)
        if d > best[1]:
            best = (r, d)
    return best


def _rust_reference_po2():
    m = re.search(r"OER_REFERENCE_PO2_MMHG:\s*f64\s*=\s*([0-9.]+)", OXY.read_text())
    return float(m.group(1))


def scan() -> dict:
    p_full = _rust_reference_po2()

    half_lives = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 48.0]
    by_half_life = []
    for th in half_lives:
        r, d = best_rate(th, PHI_CRIT_DEFAULT)
        by_half_life.append({
            "t_half_h": th, "optimal_mw_cm2": r,
            "interior": SCAN_LO * 1.001 < r < SCAN_HI * 0.999,
            "duration_h": TOTAL_FLUENCE_J * 1000.0 / r / 3600.0,
            "delivered": d,
        })

    # How much of the answer is the uncalibrated parameter?
    by_phi_crit = []
    for pc in (10.0, 25.0, 50.0, 100.0, 200.0):
        r, _ = best_rate(1.0, pc)
        by_phi_crit.append({"phi_crit": pc, "optimal_mw_cm2": r})

    # The cost of ignoring the rate: what an oxygen-blind model over-reports
    # when it delivers the same joules at a clinical high rate.
    blind = []
    for rate in (25.0, 50.0, 100.0, 200.0, 400.0):
        o2 = o2_fraction(rate, PHI_CRIT_DEFAULT)
        blind.append({"rate_mw_cm2": rate, "o2_fraction": o2,
                      "yield_factor": oer_relative_efficacy(o2, p_full)})

    return {
        "oer_reference_po2_mmhg_from_oxygen_rs": p_full,
        "total_fluence_j_cm2": TOTAL_FLUENCE_J,
        "phi_crit_default_mw_cm2": PHI_CRIT_DEFAULT,
        "scan_range_mw_cm2": [SCAN_LO, SCAN_HI],
        "by_half_life": by_half_life,
        "by_phi_crit": by_phi_crit,
        "oxygen_blind_gap": blind,
    }


def assemble(raw: dict) -> dict:
    d = dict(raw)
    hl = raw["by_half_life"]
    d["all_optima_interior"] = all(r["interior"] for r in hl)
    d["optimum_falls_with_half_life"] = all(
        hl[i]["optimal_mw_cm2"] > hl[i + 1]["optimal_mw_cm2"] for i in range(len(hl) - 1))
    # How much of the position is the uncalibrated knob: the ratio of the
    # optimum's span to the parameter's span. Near 1 would mean the answer IS
    # the parameter; near 0 that the parameter barely matters.
    pcs = raw["by_phi_crit"]
    d["phi_crit_span"] = pcs[-1]["phi_crit"] / pcs[0]["phi_crit"]
    d["optimum_span_over_phi_crit"] = pcs[-1]["optimal_mw_cm2"] / pcs[0]["optimal_mw_cm2"]
    d["optimum_scaling_exponent"] = (
        math.log(d["optimum_span_over_phi_crit"]) / math.log(d["phi_crit_span"]))
    d["verdict"] = ("DIRECTION ANCHORED, MAGNITUDE UNCONSTRAINED"
                    if d["all_optima_interior"] and d["optimum_falls_with_half_life"]
                    else "UNRESOLVED")
    return d


def render(d: dict) -> str:
    L = [
        "# PDT fluence rate: validation (#831)",
        "",
        "*Generated by `scripts/validate_pdt_fluence_rate.py --render-only`. "
        "Pure stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}.**",
        "",
        "## What was missing",
        "",
        "PDT is dosed in J/cm², and until this layer the fluence *rate* did "
        "not exist in the engine as a quantity: `pdt_intensity_at_depth` "
        "samples the photosensitizer once, when illumination starts. That "
        "makes a ten-minute and a two-hour delivery of the same joules the "
        "same treatment, and PDT is the arm where the field says loudest that "
        "they are not.",
        "",
        "## The two opposing effects",
        "",
        "A Type II sensitizer consumes ground-state O₂ to make singlet oxygen, "
        "so illumination depletes its own substrate. Henderson & Busch (Lasers "
        "Surg Med 2006, PMID 16615136) report depletion within *seconds* at 75 "
        "mW/cm² in murine tumours. Against that, a slower delivery runs "
        "further into the sensitizer's clearance — and drug that has left "
        "cannot be excited.",
        "",
        "The oxygen term reuses `oxygen::oer_relative_efficacy`, the "
        "Alper–Howard-Flanders hyperbola #726 adopted after the engine's "
        "linear O₂ form was measured against eighty years of radiobiology and "
        "found wrong below 10 mmHg. Radiation and PDT are both oxygen-"
        "dependent arms here and are not allowed two different O₂ curves.",
        "",
        f"What the oxygen term costs at a fixed reference of "
        f"{d['oer_reference_po2_mmhg_from_oxygen_rs']:.0f} mmHg and "
        f"φ_crit = {d['phi_crit_default_mw_cm2']:.0f} mW/cm²:",
        "",
        "| fluence rate | site O₂, relative | per-photon yield |",
        "|--:|--:|--:|",
    ]
    for r in d["oxygen_blind_gap"]:
        L.append(f"| {r['rate_mw_cm2']:.0f} mW/cm² | {r['o2_fraction']:.2f} | "
                 f"{r['yield_factor']:.3f} |")
    L += [
        "",
        "## The optimum, and what sets it",
        "",
        f"Holding the total at {d['total_fluence_j_cm2']:.0f} J/cm² and varying "
        "only the rate:",
        "",
        "| sensitizer half-life | optimal rate | illumination time | interior |",
        "|--:|--:|--:|---|",
    ]
    for r in d["by_half_life"]:
        L.append(f"| {r['t_half_h']:.2f} h | {r['optimal_mw_cm2']:.0f} mW/cm² | "
                 f"{r['duration_h']:.2f} h | {'yes' if r['interior'] else 'EDGE'} |")
    L += [
        "",
        "Every optimum is **interior** — strictly inside the scanned range, "
        "not sitting on its floor. That distinction is load-bearing: a scan "
        "returning its own lower bound means the model is monotonic and the "
        "real limit is something it does not carry, which is the degenerate "
        "row the oncolytic and ablation sections each shipped once before "
        "marking it.",
        "",
        "**The model's own prediction is the first column, not the second.** "
        "The optimal rate falls monotonically as the sensitizer clears more "
        "slowly, which is a statement about *drugs* rather than doses: a "
        "fast-clearing sensitizer should be illuminated fast and a "
        "slow-clearing one slowly. Registered as P21.",
        "",
        "## How much of the answer is the uncalibrated knob",
        "",
        "`φ_crit` — the rate at which photochemical consumption matches "
        "perfusive resupply — is measured nowhere in this repository. Asking "
        "what it costs, rather than caveating it:",
        "",
        "| φ_crit | optimal rate |",
        "|--:|--:|",
    ]
    for r in d["by_phi_crit"]:
        L.append(f"| {r['phi_crit']:.0f} mW/cm² | {r['optimal_mw_cm2']:.0f} mW/cm² |")
    L += [
        "",
        f"A **{d['phi_crit_span']:.0f}×** span in the parameter moves the "
        f"optimum **{d['optimum_span_over_phi_crit']:.2f}×** — an exponent of "
        f"**{d['optimum_scaling_exponent']:.2f}**, so the optimum rides on "
        "roughly the square root of a quantity nobody has measured. That is "
        "why the verdict separates direction from magnitude: the milliwatt "
        "figure is a restatement of an assumption and the *ordering* is not.",
        "",
        "## What is NOT claimed",
        "",
        "- **No fitted rate.** No clinical fluence-rate dose-finding dataset "
        "is read here, and the numbers above would move with `φ_crit`.",
        "- **The oxygen model is quasi-steady and site-local.** Real "
        "photochemical depletion is transient, heterogeneous, and coupled to "
        "the vasculature this layer does not see.",
        "- **The clinic bounds the slow end and the model does not.** A "
        "two-hour illumination is a scheduling fact as much as a physical "
        "one, and where the model's optimum sits below what a treatment room "
        "will tolerate, the binding constraint is not in the model.",
        "",
    ]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    # `assemble` runs on BOTH paths deliberately. Reading the stored JSON and
    # rendering it directly would re-print derived fields -- the verdicts --
    # from storage instead of recomputing them, which makes every guard that
    # reads a derived field inert. `tests/test_artifact_freshness.py` checks
    # exactly this, and caught it here.
    d = assemble(json.loads(OUT_JSON.read_text()) if args.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    print(f"{OUT_MD.relative_to(REPO)}: {d['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
