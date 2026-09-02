#!/usr/bin/env python3
"""What the CAR-T arm predicts, and why its anchor cannot be used as one.

THE ANCHOR, AND WHY IT IS NOT A TARGET
--------------------------------------
Tisagenlecleucel produced an overall remission rate of 81% within three months
in relapsed or refractory B-cell ALL (ELIANA, Maude 2018, PMID 29385370).
Solid-tumour CAR-T has produced nothing comparable. The contrast is real, it is
the reason this arm exists, and it CANNOT be used to fit this model.

The reason is the one the checkpoint arm ran into: a remission is not a kill
fraction. There the fix was a ratio between two strata of ONE trial, where the
unknown mapping cancels. Here the two settings are different trials, different
diseases and different endpoints, so nothing cancels -- the mapping constant is
not shared, and a ratio between them measures the difference between two
protocols as much as between two biologies.

**So this page does not fit the arm to the ELIANA number, and does not compute
a blood-versus-solid ratio and call it a check.** Saying that plainly is worth
more than a fit that would look like one.

WHAT CAN BE PREDICTED WITHOUT A FIT
-----------------------------------
A discrimination between two failures that look identical in an outcome table.

A tumour can fail to respond because too few effectors reach it -- a delivery
failure, which every barrier in this module already describes -- or because the
cells that reach it cannot engage what they find, the antigen density being
below the threshold a CAR needs to trigger lysis. Those are different in KIND,
and they behave differently under the one intervention a clinician can most
easily make: escalating the dose rescues the first and does nothing for the
second.

That is an experimental prediction rather than a comparison to a published
number, and it is the honest thing this arm can offer.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "adoptive.rs"
OUT_MD = REPO / "analysis" / "calibration" / "adoptive-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "adoptive-validation.json"

ELIANA = {"trial": "ELIANA", "pmid": "29385370", "product": "tisagenlecleucel",
          "setting": "relapsed/refractory B-cell ALL",
          "remission_rate_pct": 81, "window": "3 months"}


def rust_constants() -> dict:
    src = RUST.read_text()
    m = re.search(r"pub const ANTIGEN_DENSITY_THRESHOLD: f64 = ([0-9.]+);", src)
    if not m:
        raise SystemExit("ERROR: ANTIGEN_DENSITY_THRESHOLD not found")
    defaults = {}
    body = src[src.index("impl Default for ExpansionKinetics"):]
    for field in ("growth_per_day", "contraction_per_day", "max_fold",
                  "memory_fraction"):
        fm = re.search(rf"{field}: ([0-9.]+),", body)
        if not fm:
            raise SystemExit(f"ERROR: ExpansionKinetics::{field} not found")
        defaults[field] = float(fm.group(1))
    return {"ANTIGEN_DENSITY_THRESHOLD": float(m.group(1)), **defaults}


def engagement(density, threshold, steepness=3.0):
    x = (max(density, 0.0) / threshold) ** steepness
    return x / (1.0 + x)


def scan() -> dict:
    c = rust_constants()
    t = c["ANTIGEN_DENSITY_THRESHOLD"]
    curve = [[d, round(engagement(d, t), 5)]
             for d in (50, 100, 200, 400, 700, 1000, 1500, 2500, 5000, 10000)]

    # The discrimination, computed here independently of the crate: a
    # delivery-limited failure has headroom under the density cap and a
    # density-limited one does not.
    def gain(density, delivery_fraction, tumour=2.0e4, infused=1.0e5):
        cap = engagement(density, t) * tumour
        base = min(infused * delivery_fraction, tumour, cap)
        high = min(infused * 10.0 * delivery_fraction, tumour, cap)
        return high / base if base > 0 else 1.0

    discrimination = [
        {"case": "delivery-limited (sparse infiltration, antigen fine)",
         "density": 5000, "delivery_fraction": 0.002,
         "gain": round(gain(5000, 0.002), 3)},
        {"case": "density-limited (infiltration fine, antigen sparse)",
         "density": 200, "delivery_fraction": 0.002,
         "gain": round(gain(200, 0.002), 3)},
        {"case": "neither (already at the tumour ceiling)",
         "density": 5000, "delivery_fraction": 1.0,
         "gain": round(gain(5000, 1.0), 3)},
    ]

    # Expansion is antigen-driven, so peak expansion tracks burden.
    def peak_fold(drive, days=28.0):
        net = c["growth_per_day"] * drive - c["contraction_per_day"] * (1 - drive)
        best = 1.0
        d = 0.0
        while d <= days:
            best = max(best, min(math.exp(net * d), c["max_fold"]))
            d += 0.5
        return best

    expansion = [{"antigen_drive": round(x, 2), "peak_fold": round(peak_fold(x), 1)}
                 for x in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0)]

    return {
        "constants": c,
        "anchor": ELIANA,
        "engagement_curve": curve,
        "discrimination": discrimination,
        "expansion_vs_drive": expansion,
    }


def assemble(raw: dict) -> dict:
    d = {x["case"].split(" (")[0]: x["gain"] for x in raw["discrimination"]}
    raw["verdict"] = {
        "anchor_used_as_target": "NO -- a remission is not a kill fraction and "
                                 "no ratio cancels across two trials",
        "discrimination_holds": ("YES" if d["delivery-limited"] > 5.0
                                 and d["density-limited"] < 1.1 else "NO"),
        "delivery_gain": d["delivery-limited"],
        "density_gain": d["density-limited"],
    }
    return raw


def render(d: dict) -> str:
    a = d["anchor"]
    v = d["verdict"]
    L = ["# The CAR-T arm: a discrimination, not a fit", "",
         "*Generated by `scripts/validate_adoptive.py --render-only`. Pure "
         "stdlib; runs offline in CI. Constants are parsed from "
         "`simulations/ferroptosis-core/src/adoptive.rs`.*", "",
         "## The anchor, and why it is not a target", "",
         f"{a['product'].title()} produced an overall remission rate of "
         f"**{a['remission_rate_pct']}%** within {a['window']} in "
         f"{a['setting']} ({a['trial']}, PMID {a['pmid']}). Solid-tumour CAR-T "
         "has produced nothing comparable. That contrast is why this arm "
         "exists and it cannot be used to fit this model.", "",
         "A remission is not a kill fraction -- the same objection that made "
         "the checkpoint arm's response rate unusable. There the fix was a "
         "RATIO between two strata of one trial, where the unknown mapping "
         "cancels. **It does not transfer here.** Blood and solid CAR-T are "
         "different trials, different diseases and different endpoints, so no "
         "constant is shared and a ratio between them would measure the "
         "difference between two protocols as much as between two biologies. "
         "This page therefore fits nothing.", "",
         "## What can be predicted without a fit", "",
         "Two tumours can fail for reasons that look identical in an outcome "
         "table: too few effectors arriving, or effectors that arrive and "
         "cannot engage what they find. The second is a THRESHOLD -- a CAR "
         "needs a minimum target density to trigger lysis at all -- and a "
         "threshold is not a barrier. They separate under dose escalation.", "",
         "| case | antigen density | gain at 10x dose |", "|---|--:|--:|"]
    for row in d["discrimination"]:
        L.append(f"| {row['case']} | {row['density']:,} | {row['gain']:.2f}x |")
    L += ["",
          f"**{v['delivery_gain']:.1f}x against {v['density_gain']:.2f}x, with "
          "the same barriers and the same poor outcome.** That is an "
          "experiment: escalate the dose and the two failures diverge. It is "
          "not a comparison to a published number, and this arm does not have "
          "one it can use.", "",
          "## Expansion tracks the antigen it consumes", "",
          "Infused cells are not the cells that do the work. Expansion is "
          "driven by antigen encounter and is therefore self-limiting: a "
          "product that clears its target stops expanding.", "",
          "| antigen drive | peak fold expansion |", "|--:|--:|"]
    for row in d["expansion_vs_drive"]:
        L.append(f"| {row['antigen_drive']:.2f} | {row['peak_fold']:,.0f}x |")
    L += ["",
          "That is the shape behind an observation the field reports and this "
          "model was not fitted to: peak expansion tracks tumour burden. It is "
          "reported as a SHAPE. The rate constants behind it are placeholders "
          "and no cellular-kinetics dataset in this repository constrains "
          "them.", "",
          "## What this does not establish", "",
          "- **Nothing is fitted.** Not the density threshold, which varies by "
          "orders of magnitude with construct and target; not the expansion "
          "rates; not the cytokine scale. Every absolute number here is a "
          "placeholder.",
          "- **The discrimination is a prediction about EXPERIMENTS**, not a "
          "reproduction of one. No published series has, to this project's "
          "knowledge, escalated dose in matched density-stratified cohorts.",
          "- **The threshold is sharp in this model and graded in reality.** "
          "Antigen density varies within a tumour, so a real failure is a "
          "mixture of both modes rather than one of them.", ""]
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
    print(f"  discrimination {d['verdict']['discrimination_holds']}: "
          f"{d['verdict']['delivery_gain']:.1f}x vs "
          f"{d['verdict']['density_gain']:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
