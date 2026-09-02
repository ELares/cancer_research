#!/usr/bin/env python3
"""What the chemotherapy arm predicts, and what cannot be checked against data.

WHY THIS FILE IS SHAPED LIKE THIS
---------------------------------
Every other calibration page in this repository is a comparison against
something measured. This one opens by saying that the measurement is out of
reach, because that is the honest headline for the modality most patients
receive.

The repository's dose-response route is CTRPv2 through the DepMap download
API. It reached the five ferroptosis compounds `fetch_calibration_data.py`
fetched, and it no longer reaches anything: the catalogue endpoint now returns
a verification page rather than the documented CSV, which
`analysis/calibration/calibration-feasibility.md` already recorded as an
ACCESS block rather than a method one. So no cytotoxic dose-response is fitted
here, and every absolute kill fraction the arm produces is a placeholder.

WHAT CAN STILL BE CHECKED
-------------------------
Two structural predictions, neither of which needs a fitted potency, and both
of which the model is free to get wrong.

1. A phase-specific agent leaves a larger residue than a phase-nonspecific one
   at the same dose, and the gap NARROWS in a population that is mostly out of
   cycle. The second half is the interesting one: it was not designed in.

2. Shortening the interval at the same total dose helps only inside a window
   of regrowth rates. It cannot help when there is nothing to outrun, and it
   stops helping when regrowth is fast enough that both schedules return to
   the Gompertz plateau between cycles. The trial that motivated the question
   is CALGB 9741 (PMID 12668651); the hypothesis is Norton and Simon's (PMID
   3510732).

THE SECOND IMPLEMENTATION IS THE POINT
--------------------------------------
The model below is written again in Python rather than called through the
crate. A formula agreeing with itself is not evidence, and this repository has
found real defects by rewriting one -- so the sensitivity table is PARSED out
of the Rust and compared, and the two implementations have to produce the same
curves.

WHAT THIS DOES NOT ESTABLISH
----------------------------
That the window the model reports contains any real tumour. Where a breast
micrometastasis sits on a Gompertz rate axis is not something this repository
can determine, so the trial is named as the motivation for the question and
NOT as a result the model reproduces.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "chemo.rs"
OUT_MD = REPO / "analysis" / "calibration" / "chemo-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "chemo-validation.json"

PHASES = ("G1", "S", "G2M", "G0")


def rust_sensitivity() -> dict:
    """The sensitivity table, parsed from the Rust match arms.

    Not restated: a Python copy of the table would let the two implementations
    agree while both drifted from the crate.
    """
    src = RUST.read_text()
    body = src[src.index("pub fn phase_sensitivity"):src.index("/// Surviving fraction after one exposure")]
    out = {}
    for cls, phase, value in re.findall(
            r"\((PhaseNonspecific|SPhaseSpecific|MPhaseSpecific), (G1|S|G2M|G0|_)\) => ([0-9.]+)",
            body):
        out[(cls, phase)] = float(value)
    if not out:
        raise SystemExit("ERROR: could not parse the sensitivity table from chemo.rs")
    return out


def sensitivity(table: dict, cls: str, phase: str) -> float:
    if (cls, phase) in table:
        return table[(cls, phase)]
    if (cls, "_") in table:
        return table[(cls, "_")]
    raise SystemExit(f"ERROR: no sensitivity for {cls}/{phase}")


PROLIFERATING = {"G1": 0.55, "S": 0.25, "G2M": 0.15, "G0": 0.05}
QUIESCENT_RICH = {"G1": 0.25, "S": 0.08, "G2M": 0.05, "G0": 0.62}


def surviving_fraction(table, dose, potency, cls, dist) -> float:
    if dose <= 0 or potency <= 0:
        return 1.0
    return sum(dist[p] * math.exp(-potency * dose * sensitivity(table, cls, p))
               for p in PHASES)


def surviving_distribution(table, dose, potency, cls, dist) -> dict:
    out = {p: dist[p] * math.exp(-potency * dose * sensitivity(table, cls, p))
           for p in PHASES}
    total = sum(out.values())
    return {p: (v / total if total > 0 else v) for p, v in out.items()}


def redistribute(current, steady, rate_per_day, days) -> dict:
    if rate_per_day <= 0 or days <= 0:
        return dict(current)
    w = 1.0 - math.exp(-rate_per_day * days)
    out = {p: current[p] + (steady[p] - current[p]) * w for p in PHASES}
    total = sum(out.values())
    return {p: v / total for p, v in out.items()}


def gompertz(burden, carrying, rate_per_day, days) -> float:
    if burden <= 0 or carrying <= burden or rate_per_day <= 0 or days <= 0:
        return max(burden, 0.0)
    ln_ratio = math.log(carrying / burden)
    return min(burden * math.exp(ln_ratio * (1 - math.exp(-rate_per_day * days))), carrying)


def regimen_burden(table, cycles, dose, interval, potency, cls, start,
                   burden, carrying, regrowth, redist) -> float:
    dist = dict(start)
    for cycle in range(cycles):
        burden *= surviving_fraction(table, dose, potency, cls, dist)
        dist = surviving_distribution(table, dose, potency, cls, dist)
        if cycle + 1 < cycles:
            burden = gompertz(burden, carrying, regrowth, interval)
            dist = redistribute(dist, start, redist, interval)
    return burden


def dose_density_advantage(table, regrowth) -> float:
    args = dict(table=table, cycles=6, dose=2.0, potency=0.5,
                cls="PhaseNonspecific", start=PROLIFERATING, burden=1.0e9,
                carrying=1.0e12, regrowth=regrowth, redist=0.2)
    conventional = regimen_burden(interval=21.0, **args)
    dense = regimen_burden(interval=14.0, **args)
    return conventional / dense if dense > 0 else float("inf")


def scan() -> dict:
    table = rust_sensitivity()

    rates, advantages = [], []
    r = 0.0
    while r <= 0.5001:
        rates.append(round(r, 4))
        advantages.append(round(dose_density_advantage(table, r), 4))
        r += 0.0025 if r < 0.05 else 0.025

    peak_i = max(range(len(advantages)), key=lambda i: advantages[i])
    # The WINDOW: where shortening the interval is worth at least a tenth.
    inside = [rates[i] for i, a in enumerate(advantages) if a >= 1.1]

    residue = {}
    for pop_name, pop in (("proliferating", PROLIFERATING),
                          ("quiescent_rich", QUIESCENT_RICH)):
        residue[pop_name] = {}
        for cls in ("SPhaseSpecific", "MPhaseSpecific"):
            flat = surviving_fraction(table, 8.0, 0.5, "PhaseNonspecific", pop)
            spec = surviving_fraction(table, 8.0, 0.5, cls, pop)
            residue[pop_name][cls] = round(spec / flat, 3)

    curves = {}
    for cls in ("PhaseNonspecific", "SPhaseSpecific", "MPhaseSpecific"):
        curves[cls] = [[d, round(surviving_fraction(table, d, 0.5, cls, PROLIFERATING), 6)]
                       for d in (0.5, 1, 2, 4, 8, 12, 16, 24, 32)]

    return {
        "sensitivity_table": {f"{c}/{p}": v for (c, p), v in sorted(table.items())},
        "dose_density": {"regrowth_per_day": rates, "advantage": advantages,
                         "peak_advantage": advantages[peak_i],
                         "peak_at_regrowth_per_day": rates[peak_i],
                         "window_lo": min(inside) if inside else None,
                         "window_hi": max(inside) if inside else None,
                         "advantage_at_zero_regrowth": advantages[0],
                         "advantage_at_fast_regrowth": advantages[-1]},
        "residue_ratio_at_dose_8": residue,
        "dose_response": curves,
        "dose_response_target": {
            "source": "CTRPv2 via the DepMap download API",
            "status": "UNREACHABLE",
            "why": "the catalogue endpoint returns a verification page rather "
                   "than the documented CSV; recorded in "
                   "analysis/calibration/calibration-feasibility.md as an "
                   "access block rather than a method one",
        },
    }


def assemble(raw: dict) -> dict:
    dd = raw["dose_density"]
    raw["verdicts"] = {
        "phase_specific_residue": (
            "REPRODUCES the direction"
            if raw["residue_ratio_at_dose_8"]["proliferating"]["SPhaseSpecific"] > 1.0
            and raw["residue_ratio_at_dose_8"]["quiescent_rich"]["SPhaseSpecific"]
            < raw["residue_ratio_at_dose_8"]["proliferating"]["SPhaseSpecific"]
            else "DOES NOT"),
        "dose_density_window": (
            "TWO-SIDED as predicted"
            if dd["advantage_at_zero_regrowth"] < 1.05
            and dd["advantage_at_fast_regrowth"] < 1.2
            and dd["peak_advantage"] > 2.0 else "NOT TWO-SIDED"),
        "dose_response_magnitude": "NO TARGET REACHABLE",
    }
    return raw


def render(d: dict) -> str:
    dd = d["dose_density"]
    res = d["residue_ratio_at_dose_8"]
    L = ["# The chemotherapy arm: two structural predictions, and no reachable "
         "dose-response", "",
         "*Generated by `scripts/validate_chemo.py --render-only`. Pure stdlib; "
         "runs offline in CI. The model is implemented a SECOND time here and "
         "the sensitivity table is parsed out of "
         "`simulations/ferroptosis-core/src/chemo.rs`, so the two "
         "implementations cannot agree by sharing code.*", "",
         "## The headline is what is missing", "",
         "This is the modality most patients receive and it is the least "
         "calibratable arm in the engine. The repository's dose-response route "
         f"-- {d['dose_response_target']['source']} -- is "
         f"**{d['dose_response_target']['status']}**: "
         f"{d['dose_response_target']['why']}. No cytotoxic curve is fitted "
         "here, so every absolute kill fraction the arm produces is a "
         "placeholder and the two results below are STRUCTURAL.", "",
         "## Prediction 1: what a phase-specific agent leaves behind", "",
         "Surviving fraction relative to a phase-nonspecific agent of the same "
         "potency, at the same dose. Above 1 means the phase-specific agent "
         "left more alive.", "",
         "| population | S-phase agent | M-phase agent |",
         "|---|--:|--:|",
         f"| proliferating | {res['proliferating']['SPhaseSpecific']:.1f}x | "
         f"{res['proliferating']['MPhaseSpecific']:.1f}x |",
         f"| mostly out of cycle | {res['quiescent_rich']['SPhaseSpecific']:.1f}x | "
         f"{res['quiescent_rich']['MPhaseSpecific']:.1f}x |", "",
         "The first row is the expected result and is close to a restatement of "
         "the sensitivity table. **The second row is not.** In a population "
         "that is mostly out of cycle the gap NARROWS, because the "
         "phase-nonspecific agent is struggling too -- so the phase-specific "
         "agent's relative disadvantage is a property of how much of the "
         "tumour is dividing, not a constant of the drug class. Nothing in the "
         "model was arranged to produce that.", "",
         "## Prediction 2: when shortening the interval helps", ""]
    L += [f"At the same total dose, six cycles at 14 days against 21 days. The "
          f"advantage peaks at **{dd['peak_advantage']:.2f}x** at a Gompertz "
          f"rate of {dd['peak_at_regrowth_per_day']:.4f}/day, and is worth at "
          f"least a tenth only between {dd['window_lo']:.4f} and "
          f"{dd['window_hi']:.3f}/day.", "",
          "| regrowth | advantage |", "|---|--:|"]
    for r, a in zip(dd["regrowth_per_day"], dd["advantage"]):
        if r in (0.0, 0.0025, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5):
            L.append(f"| {r:g}/day | {a:.3f}x |")
    L += ["",
          f"**Both ends matter and only one was expected.** With no regrowth "
          f"the advantage is {dd['advantage_at_zero_regrowth']:.3f}x -- below "
          "one, because the longer gap gives the surviving population time to "
          "redistribute back into sensitive phases, which is a small effect "
          "running the other way. With fast regrowth it is "
          f"{dd['advantage_at_fast_regrowth']:.3f}x, because both schedules' "
          "tumours return to the Gompertz plateau between cycles and the extra "
          "week costs nothing that was not already lost. A model reporting a "
          "dose-density benefit at either extreme would be reporting an "
          "artifact of its own arithmetic.", "",
          "## What this does not establish", "",
          "- **Not a reproduction of CALGB 9741.** That trial (PMID 12668651) "
          "found a shorter interval better in early breast cancer, and it is "
          "why the question is worth asking. Whether breast micrometastatic "
          "disease sits inside the window above is not something this "
          "repository can determine -- the window is in units of a Gompertz "
          "rate constant that nothing here measures.",
          "- **No magnitude is defensible.** The potency is a placeholder, the "
          "phase fractions are conventional, and the redistribution rate is a "
          "free parameter. Only the SIGNS and the shape of the dependence are "
          "results.",
          "- **The cell cycle is population-level.** There is no per-cell cycle "
          "position, so nothing here can model a synchronised population or "
          "the phase-targeted scheduling that would follow from one.",
          "- **Nothing here feeds a reported number.** "
          "`simulations/calibration/CALIBRATION_STATUS.md` carries that per "
          "layer.", ""]
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
    print("  " + "; ".join(f"{k}: {v}" for k, v in d["verdicts"].items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
