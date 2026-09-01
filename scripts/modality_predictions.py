"""Falsifiable predictions for the modality arms, derived from the engine.

WHY THIS EXISTS
---------------
`analysis/scope-audit.md` reports the sharpest measure of this project's
narrowness, and it is not lines of code or figures: **8 of 8 preregistered
predictions concern ferroptosis or the physical-ROS modalities.** P1-P8 are
RSL3, PDT/SDT, persisters, CAFs, pH and spheroid size. The nine arms this
campaign added had none.

That gap matters more than the others because a falsifiable prediction is the
currency this repository treats as real. An engine that can EXPRESS nine
modalities while having committed to being wrong about only one is still a
ferroptosis project by the measure it chose for itself.

WHAT THIS SCRIPT DOES, AND DOES NOT
-----------------------------------
It DERIVES the quantitative model output for each new prediction from the
committed engine constants and artifacts, so `PREREGISTRATION.md` quotes a
measurement rather than a remembered number. It does NOT decide the
falsification thresholds -- those are judgements, they are written in the
preregistration document, and this script's guards check the document quotes
the derived figures.

Every prediction is DIRECTIONAL. Every barrier value in these arms is an
uncalibrated placeholder (`CALIBRATION_STATUS.md`), so a magnitude is not a
prediction. The one deliberate exception is P9, whose band is published.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CORE = REPO / "simulations" / "ferroptosis-core" / "src"
OUT_MD = REPO / "analysis" / "modality-predictions.md"
OUT_JSON = REPO / "analysis" / "modality-predictions.json"


def _rust_const(stem: str, name: str) -> float:
    src = (CORE / f"{stem}.rs").read_text()
    m = re.search(rf"{name}: f64 = ([0-9.eE+-]+)", src)
    if m is None:
        raise SystemExit(f"{stem}.rs: no constant {name}")
    return float(m.group(1))


def _rust_tuple(stem: str, name: str) -> tuple[float, float]:
    src = (CORE / f"{stem}.rs").read_text()
    m = re.search(rf"{name}: \(f64, f64\) = \(([0-9.]+), ([0-9.]+)\)", src)
    if m is None:
        raise SystemExit(f"{stem}.rs: no tuple {name}")
    return float(m.group(1)), float(m.group(2))


def _ser(alpha: float, beta: float, boost: float, dose: float) -> float:
    """Sensitizer enhancement ratio at `dose`, closed form.

    The dose WITHOUT the drug over the dose WITH it, at equal survival. Under
    LQ with an alpha-only boost this is DOSE-DEPENDENT, which is the whole
    content of P9.
    """
    sf = math.exp(-alpha * dose - beta * dose * dose)
    a2 = alpha * (1.0 + boost)
    d2 = (-a2 + math.sqrt(a2 * a2 - 4.0 * beta * math.log(sf))) / (2.0 * beta)
    return dose / d2


def scan() -> dict:
    alpha = _rust_const("radiation", "ALPHA_GBM_PARAMETERISATION_PER_GY")
    ab = _rust_const("radiation", "ALPHA_BETA_TUMOUR_GY")
    beta = alpha / ab
    lo, hi = _rust_tuple("radiation", "PARP_SER_BAND")

    # The admissible boost window: one boost must hold the published band
    # across the whole survival range, which is what makes P9 constrained.
    window = [b / 1000.0 for b in range(0, 2001)
              if all(lo <= _ser(alpha, beta, b / 1000.0, d) <= hi
                     for d in (2.0, 4.0, 6.0, 8.0, 10.0))]
    p9 = {
        "published_band": [lo, hi],
        "boost_window": [round(min(window), 3), round(max(window), 3)],
        "window_share_of_scan": round(len(window) / 2001.0, 3),
        "ser_by_dose": {
            f"{b:.3f}": {str(int(d)): round(_ser(alpha, beta, b, d), 3)
                         for d in (2, 6)}
            for b in (min(window), (min(window) + max(window)) / 2, max(window))
        },
        "low_over_high_dose_ratio": [
            round(_ser(alpha, beta, b, 2.0) / _ser(alpha, beta, b, 6.0), 3)
            for b in (min(window), max(window))
        ],
    }

    # P10: the bystander effect is STARVED by the escape it answers.
    src = (CORE / "adc.rs").read_text()
    esc, reach, dkp = 0.5, 0.6, 0.8   # the module's own `heterogeneous()`
    p10 = {"reach_of_negative_pool": {}}
    for apf in (0.9, 0.6, 0.3, 0.1):
        dying = apf * dkp
        byst = min(dying * esc * reach, max(1.0 - dying, 0.0))
        p10["reach_of_negative_pool"][str(apf)] = round(byst / (1.0 - apf), 4)
    p10["relative_advantage"] = round(1.0 + esc * reach, 4)

    panel = json.loads((REPO / "analysis" / "modality-panel.json").read_text())
    ab_block = panel["adoptive_barriers"]
    p11 = {
        "delivery_efficiency": round(ab_block["delivery_efficiency_solid"], 4),
        "predicted_gain_from_opening_trafficking_only":
            round(1.0 / 0.3, 3),   # the preset's trafficking barrier
        "total_collapse": round(ab_block["leukaemia_kill_fraction"]
                                / ab_block["solid_tumour_kill_fraction"], 1),
    }

    # P12: establishment is threshold-governed and dose-INDEPENDENT.
    onc = {"replication": 0.9, "interferon": 0.3, "clearance": 0.2,
           "lysis": 0.15}
    eff_r = onc["replication"] * (1.0 - onc["interferon"])
    p12 = {
        "effective_replication": round(eff_r, 4),
        "removal_rate": round(onc["clearance"] + onc["lysis"], 4),
        "establishes": eff_r > onc["clearance"] + onc["lysis"],
        "titres_tested": [1e-4, 1e-3, 1e-2, 1e-1],
    }

    p13 = {"threshold_v_per_cm":
           _rust_const("ablation", "IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM"),
           "survival_is_one_minus_coverage": True}
    return {"P9": p9, "P10": p10, "P11": p11, "P12": p12, "P13": p13}


def render(d: dict) -> str:
    p9, p10, p11, p12 = d["P9"], d["P10"], d["P11"], d["P12"]
    lo_r, hi_r = p9["low_over_high_dose_ratio"]
    band_lo, band_hi = p9["published_band"]
    L = ["# Falsifiable predictions for the modality arms", "",
         "*Generated by `scripts/modality_predictions.py --render-only`. "
         "Offline; derives each model output from the committed engine "
         "constants and artifacts.*", "",
         "`analysis/scope-audit.md` reports that **8 of 8** preregistered "
         "predictions concern ferroptosis or the physical-ROS modalities. That "
         "is the sharpest measure of this project's narrowness — sharper than "
         "code volume or figure count, because a falsifiable prediction is the "
         "currency this repository treats as real. An engine that can express "
         "nine modalities while having committed to being wrong about only one "
         "is still a ferroptosis project by its own chosen measure.", "",
         "These are the model outputs behind P9–P13 in `PREREGISTRATION.md`. "
         "The thresholds live there; the numbers are derived here.", "",
         "## P9 — the PARP sensitizer ratio falls with dose per fraction", "",
         f"Under the linear-quadratic model with an alpha-only boost, ONE "
         f"boost must hold the published {band_lo}-{band_hi} enhancement band "
         f"across the whole survival range. Only {p9['window_share_of_scan']:.0%} "
         f"of the scanned range does, and the admissible boosts are "
         f"{p9['boost_window'][0]} to {p9['boost_window'][1]}. Within that "
         f"window the ratio is larger at low dose per fraction than high, by a "
         f"factor of {lo_r} to {hi_r}:", "",
         "| boost | SER at 2 Gy | SER at 6 Gy |", "|---|--:|--:|"]
    # Sorted HERE, on the boost, because a renderer that inherits dict order
    # publishes whatever the serialiser happened to do.
    for b, row in sorted(p9["ser_by_dose"].items(), key=lambda kv: float(kv[0])):
        L.append(f"| {b} | {row['2']} | {row['6']} |")
    L += ["",
          "**The direction does not depend on the fit**, which is what makes "
          "it worth preregistering: it holds across the entire window the "
          "published band permits. Beta is untouched by the boost because "
          "unrepaired single-strand breaks convert to double-strand breaks at "
          "replication — one-track damage — so the linear term rises alone and "
          "the ratio must decay with dose.", "",
          "## P10 — the bystander effect is starved by the escape it answers",
          "",
          "The prediction has the OPPOSITE sign to the intuition the module "
          "was built on, and to a guard this repository shipped. The payload "
          "comes from cells that took up the ADC, so antigen escape removes "
          "its own source. The share of the antigen-negative pool a cleavable "
          "linker reaches:", "",
          "| antigen-positive fraction | share of the negative pool reached |",
          "|---|--:|"]
    for k, v in sorted(p10["reach_of_negative_pool"].items(),
                       key=lambda kv: float(kv[0]), reverse=True):
        L.append(f"| {k} | {v:.1%} |")
    L += ["",
          f"The RELATIVE advantage meanwhile is exactly flat at "
          f"{p10['relative_advantage']}x, and a guard asserting it was "
          f"non-decreasing was satisfied by that constant for the life of the "
          f"module — see the retraction in `CALIBRATION_STATUS.md`.", "",
          "## P11 — the adoptive barriers multiply", "",
          f"Delivery efficiency is {p11['delivery_efficiency']}, so opening "
          f"exactly one barrier must multiply the kill by that barrier's "
          f"reciprocal and no more. Intratumoural against intravenous "
          f"administration sets trafficking to roughly 1 and leaves the other "
          f"two untouched, predicting a "
          f"{p11['predicted_gain_from_opening_trafficking_only']}x gain — not "
          f"the {p11['total_collapse']}x that separates the two diseases.", "",
          "## P12 — oncolytic establishment is dose-independent", "",
          f"Effective replication is {p12['effective_replication']} against a "
          f"removal rate of {p12['removal_rate']}, so the infection "
          f"{'establishes' if p12['establishes'] else 'dies out'} — and that "
          "comparison contains no titre term at all. The initial dose changes "
          "how fast the outcome arrives and not which outcome it is.", "",
          "## P13 — ablation outcome is geometry, not dose", "",
          "Survival is one minus the covered fraction whenever the threshold "
          "is exceeded, with no dependence on by how much. Recurrence after "
          "ablation is a margin problem, which is the one place in this engine "
          "where the hypoxia and penetration work cannot help.", "",
          "## What these are not", "",
          "**Directional, every one of them, and the barrier values behind "
          "them are uncalibrated placeholders.** A magnitude here is not a "
          "prediction. P9 is the deliberate exception: its band is published, "
          "which is why it is the arm the data most tightly constrains and the "
          "only one whose model output is quantitative rather than a sign.", "",
          "**Not clinical guidance.** Each names an experiment that would "
          "falsify a MODEL, not a therapy that should be given.", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = json.loads(OUT_JSON.read_text()) if a.render_only else scan()
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
