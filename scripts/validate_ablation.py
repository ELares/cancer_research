#!/usr/bin/env python3
"""Where a thermal ablation fails, and why a coverage fraction cannot say.

WHY THIS ARM NEEDED SOMETHING OTHER THAN A FIT
----------------------------------------------
`analysis/modality-calibration.md` reports this arm UNCONSTRAINED: its target
is a temperature threshold, a threshold is nearly binary, and almost the whole
scanned parameter range reproduces it. That is not a failure of the fitting --
it is what happens when the observable carries almost no information about the
parameter.

The information is somewhere else. Ablation's clinical failure is SPATIAL:
tissue beside a large vessel survives, because flowing blood carries heat away
and that tissue never reaches a lethal thermal dose. Perivascular local
progression after thermal ablation is a recognised problem (PMID 35114665).

WHAT THE MODEL PREDICTS
-----------------------
A surviving sleeve with a radius, rather than a coverage percentage. And a
DISCRIMINATION: irreversible electroporation is not thermal, so a heat sink
does not apply to it and the sleeve is exactly zero. The two modalities fail
in different PLACES, which is checkable in a way a fraction is not -- and it
is the documented reason electroporation is reached for near vessels.

WHAT IS ANCHORED AND WHAT IS NOT
--------------------------------
The DIRECTION is anchored: perivascular tumours progress locally more often
after thermal ablation. The sleeve's SIZE is not anchored to anything. The
cooling length is a placeholder that depends on vessel calibre and flow, and
this layer represents neither.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "ablation.rs"
OUT_MD = REPO / "analysis" / "calibration" / "ablation-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "ablation-validation.json"

BODY_C = 37.0
MINUTES = 5.0
CEM43_THRESHOLD = 240.0
COOLING_LENGTH_MM = 2.0


def cem43(temperature_c, minutes):
    """The crate's own thermal-dose form, re-implemented."""
    if temperature_c >= 43.0:
        return minutes * (0.5 ** (43.0 - temperature_c))
    return minutes * (0.25 ** (43.0 - temperature_c))


def sleeve_radius(applicator_c, cooling_length=COOLING_LENGTH_MM, max_mm=20.0):
    step = max_mm / 500.0
    d = 0.0
    while d <= max_mm:
        cooling = 1.0 - math.exp(-d / cooling_length)
        t = BODY_C + max(applicator_c - BODY_C, 0.0) * cooling
        if cem43(t, MINUTES) >= CEM43_THRESHOLD:
            return round(d, 3)
        d += step
    return max_mm


def scan() -> dict:
    src = RUST.read_text()
    assert "pub fn electroporation_failure_radius_mm" in src, (
        "the crate no longer carries the electroporation contrast")
    temps = [45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0]
    rows = []
    for t in temps:
        # A TOTAL FAILURE is not a sleeve, and reporting it as one would be the
        # same error the oncolytic page had to fix: an applicator that cannot
        # reach a lethal dose even infinitely far from a vessel fails
        # EVERYWHERE, and the scan limit it returns reads as a very wide
        # perivascular sleeve while meaning something else entirely.
        unimpeded_lethal = cem43(t, MINUTES) >= CEM43_THRESHOLD
        rows.append({"applicator_c": t,
                     "thermal_sleeve_mm": sleeve_radius(t),
                     "electroporation_sleeve_mm": 0.0,
                     "total_failure": not unimpeded_lethal})
    lengths = [0.5, 1.0, 2.0, 4.0, 8.0]
    sensitivity = [{"cooling_length_mm": l,
                    "sleeve_at_60c_mm": sleeve_radius(60.0, l)} for l in lengths]
    return {
        "body_temperature_c": BODY_C, "minutes": MINUTES,
        "cem43_threshold": CEM43_THRESHOLD,
        "cooling_length_mm": COOLING_LENGTH_MM,
        "by_applicator_temperature": rows,
        "sensitivity_to_cooling_length": sensitivity,
        "anchor": {"pmid": "35114665",
                   "claim": "local tumour progression remains a challenge in "
                            "perivascular hepatocellular carcinoma after "
                            "radiofrequency ablation",
                   "kind": "DIRECTION only -- no sleeve radius is published here"},
    }


def assemble(raw: dict) -> dict:
    rows = [r for r in raw["by_applicator_temperature"] if not r["total_failure"]]
    raw["n_total_failure"] = sum(1 for r in raw["by_applicator_temperature"]
                                 if r["total_failure"])
    sleeves = [r["thermal_sleeve_mm"] for r in rows]
    raw["verdict"] = {
        "sleeve_shrinks_with_temperature":
            "YES" if all(a >= b for a, b in zip(sleeves, sleeves[1:])) else "NO",
        "electroporation_contrast": "YES -- exactly zero at every temperature",
        "sleeve_size_anchored": "NO",
        "sleeve_range_mm": [min(sleeves), max(sleeves)],
    }
    return raw


def render(d: dict) -> str:
    v = d["verdict"]
    L = ["# Where a thermal ablation fails, and why a coverage fraction cannot say",
         "",
         "*Generated by `scripts/validate_ablation.py --render-only`. Pure "
         "stdlib; runs offline in CI. The thermal-dose form is re-implemented "
         "here rather than imported.*", "",
         "## Why this arm needed something other than a fit", "",
         "The calibration page reports this arm UNCONSTRAINED, and that is not "
         "a failure of the fitting: its target is a temperature threshold, a "
         "threshold is nearly binary, and almost the whole scanned range "
         "reproduces it. The observable carries almost no information about "
         "the parameter.", "",
         "The information is somewhere else. Ablation's clinical failure is "
         "SPATIAL -- tissue beside a large vessel survives because flowing "
         "blood carries the heat away.", "",
         "## The surviving sleeve", "",
         "| applicator temperature | thermal sleeve | electroporation sleeve |",
         "|--:|--:|--:|"]
    for r in d["by_applicator_temperature"]:
        thermal = ("fails everywhere" if r["total_failure"]
                   else f"{r['thermal_sleeve_mm']:.2f} mm")
        L.append(f"| {r['applicator_c']:.0f} °C | {thermal} | "
                 f"{r['electroporation_sleeve_mm']:.2f} mm |")
    L += ["",
          f"The sleeve shrinks as the applicator gets hotter "
          f"({v['sleeve_range_mm'][1]:.2f} mm down to "
          f"{v['sleeve_range_mm'][0]:.2f} mm across the range above) and never "
          "closes. It is a RADIUS, not a percentage: the survivors are in a "
          "specific place and a clinician can look there, which is what makes "
          "it checkable where a coverage fraction is not.", "",
          "**Electroporation's sleeve is exactly zero**, and that is the "
          "point rather than a stub. It kills by permeabilising membranes "
          "with an electric field and deposits little heat, so a heat sink "
          "removes nothing that matters to it -- which is the documented "
          "reason it is reached for near vessels. The two modalities fail in "
          "different PLACES.", "",
          "## What moves the answer", "",
          "| cooling length | sleeve at 60 °C |", "|--:|--:|"]
    for r in d["sensitivity_to_cooling_length"]:
        L.append(f"| {r['cooling_length_mm']:.1f} mm | "
                 f"{r['sleeve_at_60c_mm']:.2f} mm |")
    L += ["",
          "The cooling length is a placeholder standing in for vessel calibre "
          "and flow rate, neither of which this layer represents, and the "
          "sleeve scales with it almost proportionally. So the SIZE of the "
          "sleeve is a restatement of that placeholder.", "",
          "## What is anchored and what is not", "",
          f"- **Direction only.** PMID {d['anchor']['pmid']}: "
          f"{d['anchor']['claim']}. That supports the existence of the effect "
          "and not its magnitude -- no sleeve radius is published there, and "
          "none is fitted here.",
          "- **The contrast is structural.** Electroporation's zero follows "
          "from its being non-thermal, not from a measurement, and a change "
          "giving it a thermal dependence breaks a test in the crate rather "
          "than passing quietly.",
          f"- **A total failure is not a sleeve, and it is marked.** "
          f"{d['n_total_failure']} of "
          f"{len(d['by_applicator_temperature'])} rows have an applicator too "
          "cool to reach a lethal dose even far from any vessel, so the scan "
          "limit they return would read as a very wide perivascular sleeve "
          "while meaning the ablation failed everywhere. Excluded from the "
          "shrinking claim.",
          "- **The threshold is the module's existing CEM43 criterion**, which "
          "is a thermal DOSE and not a temperature: a brief exposure at a "
          "lethal temperature is not lethal, and the guard for that exists "
          "because a mutation replacing the dose with a bare temperature "
          "survived every other test.", ""]
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
    print(f"  sleeve {d['verdict']['sleeve_range_mm'][0]}-"
          f"{d['verdict']['sleeve_range_mm'][1]} mm; electroporation zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
