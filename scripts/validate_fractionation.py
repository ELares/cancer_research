#!/usr/bin/env python3
"""Check the fractionation layer against trials it was not fitted to.

WHAT IS BEING CHECKED, AND WHY IT CAN FAIL
------------------------------------------
The engine's radiation arm can now express a SCHEDULE -- fractions, the size of
each, the gap between them, the calendar the course runs over. That machinery
has an unusual property for this repository: it can be checked against
randomised trials, because radiotherapy trials publish their two schedules and
report whether the outcomes differed.

The check runs BACKWARDS, which is what makes it a test rather than a fit. Two
schedules reported as not differing imply, under the linear-quadratic model,
the alpha/beta ratio at which they deliver the same effect:

    D1 (1 + d1/x) = D2 (1 + d2/x)   =>   x = (D2 d2 - D1 d1) / (D1 - D2)

Nothing about that expression knows what tissue it is describing. So the value
it returns can be compared against alpha/beta estimates the radiobiology
literature derives by other means, from other data. Two numbers out of a trial
PROTOCOL predict a third quantity nobody put in.

THE HONEST HALF: ONE LEG PASSES AND ONE FAILS
---------------------------------------------
Prostate passes and breast does not, and the failure is the more interesting
of the two. It is reported at the same weight as the pass, and the arithmetic
that produces it is the same arithmetic -- there is no separate treatment for
the leg that disagrees.

WHAT THIS IS NOT
----------------
Not a trial analysis, and it cannot become one. A non-inferiority result says
the arms did not differ by more than a pre-specified margin; it does not say
they were equal, so the implied alpha/beta inherits the width of that margin
and this page states the point estimate as an ANCHOR rather than as an
estimate with a confidence interval. Reconstructing an interval would need the
trials' own outcome distributions, which are not in this repository.

Nor does it license any clinical statement. The layer is not used by any
number the manuscript's quantitative chapters report, which
`simulations/calibration/CALIBRATION_STATUS.md` records per layer.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "radiation.rs"
DATA = REPO / "analysis" / "calibration" / "fractionation_trials.csv"
OUT_MD = REPO / "analysis" / "calibration" / "fractionation-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "fractionation-validation.json"


def rust_constants() -> dict:
    """The constants, parsed from the Rust rather than restated here.

    A Python copy of a Rust number is a second source of truth that goes stale
    silently; this repository has fixed that defect in three other validation
    scripts and the guard below asserts the parse actually found something.
    """
    src = RUST.read_text()
    out = {}
    for name in ("ALPHA_BETA_TUMOUR_GY", "ALPHA_BETA_LATE_GY",
                 "ALPHA_HEAD_NECK_PER_GY", "ALPHA_GBM_PARAMETERISATION_PER_GY",
                 "REPOP_DOUBLING_DAYS", "REPOP_KICKOFF_DAYS",
                 "SUBLETHAL_REPAIR_HALF_TIME_H"):
        m = re.search(rf"pub const {name}: f64 = ([0-9.]+);", src)
        if not m:
            raise SystemExit(f"ERROR: {name} not found in {RUST.name}")
        out[name] = float(m.group(1))
    m = re.search(r"pub const D_PROLIF_PUBLISHED_GY_PER_DAY: \(f64, f64\) = "
                  r"\(([0-9.]+), ([0-9.]+)\);", src)
    if not m:
        raise SystemExit("ERROR: D_PROLIF_PUBLISHED_GY_PER_DAY not found")
    out["D_PROLIF_PUBLISHED_GY_PER_DAY"] = (float(m.group(1)), float(m.group(2)))
    return out


def isoeffect_alpha_beta(n1, d1, n2, d2):
    """The same expression the Rust implements, re-derived here.

    Deliberately NOT imported from the crate: an independent implementation of
    a four-line formula is worth more than a shared one, because the failure
    this catches is a sign error in the Rust, and importing it would make the
    check compare the code to itself.
    """
    total1, total2 = n1 * d1, n2 * d2
    if abs(total1 - total2) < 1e-9:
        return None
    return (total2 * d2 - total1 * d1) / (total1 - total2)


def repopulation_dose_per_day(alpha, doubling_days):
    import math
    return math.log(2.0) / (alpha * doubling_days)


def scan() -> dict:
    consts = rust_constants()
    rows = []
    with DATA.open() as fh:
        for row in csv.DictReader(l for l in fh if not l.startswith("#")):
            n1 = int(row["arm_a_fractions"])
            d1 = float(row["arm_a_dose_per_fraction_gy"])
            n2 = int(row["arm_b_fractions"])
            d2 = float(row["arm_b_dose_per_fraction_gy"])
            implied = isoeffect_alpha_beta(n1, d1, n2, d2)
            lo = float(row["published_ab_lo_gy"])
            hi = float(row["published_ab_hi_gy"])
            rows.append({
                "site": row["site"],
                "trial": row["trial"],
                "pmid": row["pmid"],
                "arm_a": f"{n1} x {d1:g} Gy",
                "arm_b": f"{n2} x {d2:g} Gy",
                "arm_a_total_gy": round(n1 * d1, 2),
                "arm_b_total_gy": round(n2 * d2, 2),
                "weeks_a": float(row["arm_a_weeks"]),
                "weeks_b": float(row["arm_b_weeks"]),
                "outcome": row["outcome"],
                "implied_alpha_beta_gy": round(implied, 3) if implied else None,
                "published_lo_gy": lo,
                "published_hi_gy": hi,
                "published_source_pmid": row["published_ab_source_pmid"],
                "verdict": ("REPRODUCES" if implied is not None and lo <= implied <= hi
                            else "DISAGREES"),
            })

    d_prolif = repopulation_dose_per_day(
        consts["ALPHA_HEAD_NECK_PER_GY"], consts["REPOP_DOUBLING_DAYS"])
    d_prolif_gbm = repopulation_dose_per_day(
        consts["ALPHA_GBM_PARAMETERISATION_PER_GY"], consts["REPOP_DOUBLING_DAYS"])
    lo, hi = consts["D_PROLIF_PUBLISHED_GY_PER_DAY"]
    return {
        "constants": consts,
        "isoeffect": rows,
        "repopulation": {
            "d_prolif_head_and_neck_gy_per_day": round(d_prolif, 3),
            "d_prolif_with_gbm_alpha_gy_per_day": round(d_prolif_gbm, 3),
            "published_lo": lo, "published_hi": hi,
            "verdict": "REPRODUCES" if lo <= d_prolif <= hi else "DISAGREES",
            "gbm_verdict": "REPRODUCES" if lo <= d_prolif_gbm <= hi else "DISAGREES",
        },
    }


def assemble(raw: dict) -> dict:
    raw["n_reproduces"] = sum(1 for r in raw["isoeffect"] if r["verdict"] == "REPRODUCES")
    raw["n_legs"] = len(raw["isoeffect"])
    return raw


def render(d: dict) -> str:
    L = ["# Checking the fractionation layer against trials it was not fitted to",
         "",
         "*Generated by `scripts/validate_fractionation.py --render-only`. Pure "
         "stdlib; runs offline in CI. Constants are parsed from "
         "`simulations/ferroptosis-core/src/radiation.rs`, not restated here.*",
         "",
         "The engine can now express a schedule, and schedules are the one part "
         "of this project with randomised trials behind them. The check runs "
         "BACKWARDS: two schedules a trial reported as not differing imply, "
         "under the linear-quadratic model, the α/β at which they are "
         "equivalent. That implied value is then compared against α/β estimates "
         "the literature derives from other data entirely.", "",
         "## The isoeffect inversion", "",
         "| site | trial | arms | total dose | implied α/β | published α/β | verdict |",
         "|---|---|---|---|--:|---|---|"]
    for r in d["isoeffect"]:
        L.append(
            f"| {r['site']} | {r['trial']} (PMID {r['pmid']}) | "
            f"{r['arm_a']} vs {r['arm_b']} | "
            f"{r['arm_a_total_gy']:g} vs {r['arm_b_total_gy']:g} Gy | "
            f"{r['implied_alpha_beta_gy']:.2f} Gy | "
            f"{r['published_lo_gy']}--{r['published_hi_gy']} Gy | "
            f"**{r['verdict']}** |")

    prostate = next(r for r in d["isoeffect"] if r["site"] == "prostate")
    breast = next(r for r in d["isoeffect"] if r["site"] == "breast")
    L += ["",
          f"**Prostate reproduces.** CHHiP's two schedules imply "
          f"α/β = {prostate['implied_alpha_beta_gy']:.2f} Gy, inside the band "
          f"the literature estimates for prostate by other means. Nothing in "
          f"that computation knows it is about prostate: the inputs are "
          f"{prostate['arm_a']} and {prostate['arm_b']}, and the low ratio -- "
          f"the observation that motivated hypofractionating prostate cancer in "
          f"the first place -- falls out of the arithmetic.", "",
          f"**Breast does not, and that is the more useful row.** START-B's "
          f"schedules imply α/β = {breast['implied_alpha_beta_gy']:.2f} Gy "
          f"against a published {breast['published_lo_gy']}--"
          f"{breast['published_hi_gy']} Gy. The reason is visible in the total "
          f"dose column: the 15-fraction arm delivers "
          f"{breast['arm_a_total_gy'] - breast['arm_b_total_gy']:.2f} Gy LESS "
          f"than the 25-fraction arm and was still not inferior, so no positive "
          f"α/β makes the two EQD2-identical. Under the linear-quadratic model "
          f"alone the trial's result is not reproducible -- which is a "
          f"statement about what EQD2 leaves out (the arms also differ by two "
          f"weeks of elapsed time) and not about the trial.", ""]

    rp = d["repopulation"]
    L += ["## The repopulation prediction", "",
          "The dose each further day of treatment costs once accelerated "
          "repopulation has begun is `ln2 / (α·T_p)`. α comes from a survival "
          "parameterisation and `T_p` from a repopulation study; neither was "
          "chosen to make this number land anywhere.", "",
          "| α used | D_prolif | published band | verdict |",
          "|---|--:|---|---|",
          f"| head and neck ({d['constants']['ALPHA_HEAD_NECK_PER_GY']}/Gy) | "
          f"{rp['d_prolif_head_and_neck_gy_per_day']:.2f} Gy/day | "
          f"{rp['published_lo']}--{rp['published_hi']} Gy/day | "
          f"**{rp['verdict']}** |",
          f"| glioblastoma "
          f"({d['constants']['ALPHA_GBM_PARAMETERISATION_PER_GY']}/Gy) | "
          f"{rp['d_prolif_with_gbm_alpha_gy_per_day']:.2f} Gy/day | "
          f"same band | **{rp['gbm_verdict']}** |", "",
          "The second row is the honest limit of the first. The prediction is "
          "about head and neck cancer, and transplanting it to another tumour's "
          "α takes it outside the band it was checked against. A guard in the "
          "crate asserts that, so the caveat cannot quietly stop being true.",
          ""]

    L += ["## What this does not establish", "",
          "- **Not a trial analysis.** A non-inferiority result says the arms "
          "did not differ by more than a margin, not that they were equal. The "
          "implied α/β inherits the width of that margin, and no interval is "
          "reported here because reconstructing one needs outcome "
          "distributions this repository does not hold.",
          "- **Not a validation of the engine's kill magnitudes.** What is "
          "checked is the SCHEDULE arithmetic. Every absolute surviving "
          "fraction the layer produces still rests on α and β taken from a "
          "parameterisation.",
          "- **Not in use.** No number the manuscript's quantitative chapters "
          "report comes from this layer; "
          "`simulations/calibration/CALIBRATION_STATUS.md` carries that per "
          "layer.",
          "- **The repair half-time is uncalibrated.** Published values span "
          "roughly 0.5 to 4 hours and differ between tissues; the layer's "
          "incomplete-repair correction is direction-only, and the six-hour "
          "interval it makes an argument about is a clinical convention rather "
          "than a fitted result.", ""]
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
    print(f"  isoeffect: {d['n_reproduces']} of {d['n_legs']} reproduce; "
          f"repopulation {d['repopulation']['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
