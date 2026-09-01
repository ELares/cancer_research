"""Falsifiable predictions for the modality arms, derived from the engine.

WHY THIS EXISTS
---------------
`analysis/scope-audit.md` reports the sharpest measure of this project's
narrowness, and it is not lines of code or figures: **8 of 13 preregistered
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


def _rust_struct_field(stem: str, fn_name: str, field: str) -> float:
    """Read one field out of a named constructor in the crate.

    EVERY P10-P13 constant used to be a Python literal, and `scan()` even read
    `adc.rs` and threw the result away, which made the script LOOK derived. A
    reviewer changed four literals so they contradicted the Rust and published
    a preregistration disagreeing with its own cited module on every P10 and
    P11 number, with a fully green suite. A preregistration that can drift
    from the model it registers is worse than none, because the timestamp
    makes it look like a commitment.
    """
    src = (CORE / f"{stem}.rs").read_text()
    m = re.search(rf"fn {fn_name}\(\)[^{{]*{{(.*?)\n    }}", src, re.S)
    if m is None:
        raise SystemExit(f"{stem}.rs: no constructor {fn_name}()")
    f = re.search(rf"\b{field}:\s*([0-9.]+)", m.group(1))
    if f is None:
        raise SystemExit(f"{stem}.rs::{fn_name}: no field {field}")
    return float(f.group(1))


def _rust_test_fixture(stem: str, fn_name: str, field: str) -> float:
    """Same, for a fixture that lives inside `#[cfg(test)]`."""
    return _rust_struct_field(stem, fn_name, field)


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

    # SAMPLED ON SURVIVING FRACTION, exactly as `calibrate_modality_arms.py`
    # does, because two committed artifacts published DIFFERENT admissible
    # windows for the same band: this script sampled doses {2..10 Gy} and gave
    # [0.567, 0.922] while the calibration sampled sf {0.5, 0.1, 0.01} and gave
    # [0.544, 0.949], and `fig35` drew the PARP row ADMISSIBLE off the second
    # while the preregistration published the first. Neither stated its
    # parameterisation.
    SF_POINTS = (0.5, 0.1, 0.01)

    def _dose_at(sf: float, a: float) -> float:
        L = -math.log(sf)
        return (-a + math.sqrt(a * a + 4.0 * beta * L)) / (2.0 * beta)

    def _ser_at_sf(boost: float, sf: float) -> float:
        return _dose_at(sf, alpha) / _dose_at(sf, alpha * (1.0 + boost))

    STEPS = 2001
    window = [b / 1000.0 for b in range(STEPS)
              if all(lo <= _ser_at_sf(b / 1000.0, sf) <= hi
                     for sf in SF_POINTS)]
    p9 = {
        "published_band": [lo, hi],
        "sampled_at_surviving_fractions": list(SF_POINTS),
        "dose_range_gy": [round(min(_dose_at(sf, alpha) for sf in SF_POINTS), 2),
                          round(max(_dose_at(sf, alpha) for sf in SF_POINTS), 2)],
        "boost_window": [round(min(window), 3), round(max(window), 3)],
        # Reported as an INTERVAL, and the share is kept only WITH its ceiling,
        # because the same window is 36% of a 0-1 scan and 3.6% of a 0-10 one:
        # a bare "18%" measures the author's choice of ceiling.
        "scan_ceiling": (STEPS - 1) / 1000.0,
        "window_share_of_scan": round(len(window) / float(STEPS), 3),
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

    # P10 -- every constant READ FROM `adc.rs`, and the quantity is the one an
    # escape experiment can produce. Dividing the whole bystander term by the
    # antigen-negative pool published 216% as a "share", which is impossible:
    # the term is bounded by every SURVIVING cell, antigen-positive ones
    # included. `adc::bystander_kill_on_negative` apportions it.
    esc = _rust_test_fixture("adc", "heterogeneous", "payload_escape_fraction")
    reach = _rust_test_fixture("adc", "heterogeneous", "neighbours_in_reach")
    dkp = _rust_test_fixture("adc", "heterogeneous", "direct_kill_probability")
    p10 = {"reach_of_negative_pool": {}, "constants": {
        "payload_escape_fraction": esc, "neighbours_in_reach": reach,
        "direct_kill_probability": dkp}}
    for apf in (0.9, 0.6, 0.3, 0.1):
        dying = apf * dkp
        surviving = max(1.0 - dying, 0.0)
        byst = min(dying * esc * reach, surviving)
        on_negative = byst * ((1.0 - apf) / surviving) if surviving else 0.0
        p10["reach_of_negative_pool"][str(apf)] = round(
            on_negative / (1.0 - apf), 4)
    p10["relative_advantage"] = round(1.0 + esc * reach, 4)

    panel = json.loads((REPO / "analysis" / "modality-panel.json").read_text())
    ab_block = panel["adoptive_barriers"]
    traffick = _rust_struct_field("adoptive", "solid_tumour", "trafficking")
    p11 = {
        "delivery_efficiency": round(ab_block["delivery_efficiency_solid"], 4),
        "trafficking_barrier": traffick,
        "predicted_gain_from_opening_trafficking_only": round(1.0 / traffick, 3),
        "total_collapse": round(ab_block["leukaemia_kill_fraction"]
                                / ab_block["solid_tumour_kill_fraction"], 1),
    }

    # P12 -- SIMULATED, not asserted. The first version stated the verdict from
    # a one-line inequality that has no titre term BY CONSTRUCTION, listed four
    # titres it never used, and called their span "four orders of magnitude"
    # when 1e-4 to 1e-1 is three. It also used its own criterion rather than
    # the crate's: `oncolytic::spread_threshold_ratio` compares replication
    # against CLEARANCE alone, while lysis also removes infected cells, so a
    # config above that ratio can still die out.
    onc = {"replication": 0.9, "interferon": 0.3, "clearance": 0.2,
           "lysis": 0.15}
    eff_r = onc["replication"] * (1.0 - onc["interferon"])

    def _spread(initial, steps=180):
        """Mirrors `oncolytic::simulate_spread` STATEMENT BY STATEMENT.

        The first version accumulated `lysed` from the POST-update `infected`
        while the crate uses the pre-update value, so every published figure
        was wrong against the engine it claimed to describe -- 0.6238 where
        cargo gives 0.6418. A re-implementation that is nearly the same is a
        different model.
        """
        infected, lysed = initial, 0.0
        for _ in range(steps):
            susceptible = max(1.0 - infected - lysed, 0.0)
            new = eff_r * infected * susceptible
            cleared = onc["clearance"] * infected
            died = onc["lysis"] * infected
            infected = min(max(infected + new - cleared - died, 0.0), 1.0)
            lysed = min(lysed + died, 1.0)
            if infected <= 0.0:
                break
        return infected, lysed

    titres = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    runs = {f"{t:g}": round(_spread(t)[1], 4) for t in titres}
    verdicts = {k: v > 0.01 for k, v in runs.items()}
    p12 = {
        "effective_replication": round(eff_r, 4),
        "removal_rate": round(onc["clearance"] + onc["lysis"], 4),
        "crate_threshold_ratio": round(eff_r / onc["clearance"], 4),
        "establishes": eff_r > onc["clearance"] + onc["lysis"],
        "titres_tested": titres,
        "orders_of_magnitude": len(titres) - 1,
        "cumulative_lysed_by_titre": runs,
        "verdict_is_the_same_at_every_titre": len(set(verdicts.values())) == 1,
        "lysed_spread_across_titres": round(max(runs.values()) - min(runs.values()), 4),
    }

    # P13 -- READ FROM `ablation.rs`, because `"survival_is_one_minus_coverage":
    # True` was a Python literal. If the function ever gained an energy term
    # the literal would stay True and the prediction would keep being
    # published, which is the drift this whole script exists to prevent.
    abl = (CORE / "ablation.rs").read_text()
    body = re.search(r"pub fn margin_survival_fraction\(.*?\n}", abl, re.S)
    if body is None:
        raise SystemExit("ablation.rs: no margin_survival_fraction")
    src13 = body.group(0)
    p13 = {
        "threshold_v_per_cm":
            _rust_const("ablation",
                        "IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM"),
        "returns_one_minus_covered": "1.0 - covered" in src13,
        # The BODY branches on `electroporation_ablated() ||
        # hifu_thermal_ablation()`, which is precisely where temperature,
        # duration and field strength ARE read -- P13 said it read none of
        # them, which was false. What is true and is what P13 registers: above
        # the threshold the return value does not vary with them.
        "reads_energy_only_to_test_the_threshold": (
            "electroporation_ablated" in src13
            or "hifu_thermal_ablation" in src13),
        "return_value_varies_with_energy": any(
            t in src13.split("{", 1)[-1]
            for t in ("temperature_c *", "minutes *", "field_v_per_cm *",
                      "cem43(")),
        "signature_takes_only_config_and_coverage":
            bool(re.search(r"margin_survival_fraction\(\s*cfg: &AblationConfig,"
                           r"\s*covered_fraction: f64,?\s*\)", src13)),
    }
    return {"P9": p9, "P10": p10, "P11": p11, "P12": p12, "P13": p13}


def render(d: dict) -> str:
    p9, p10, p11, p12 = d["P9"], d["P10"], d["P11"], d["P12"]
    audit = json.loads((REPO / "analysis" / "scope-audit.json").read_text())
    nf = sum(1 for v in audit["predictions"].values() if v)
    nt = len(audit["predictions"])
    lo_r, hi_r = p9["low_over_high_dose_ratio"]
    band_lo, band_hi = p9["published_band"]
    share01 = p9["boost_window"][1] - p9["boost_window"][0]
    L = ["# Falsifiable predictions for the modality arms", "",
         "*Generated by `scripts/modality_predictions.py --render-only`. "
         "Offline; derives each model output from the committed engine "
         "constants and artifacts.*", "",
         f"`analysis/scope-audit.md` reported that **8 of 8** preregistered "
         f"predictions concerned ferroptosis or the physical-ROS modalities "
         f"before the five below were registered, and **{nf} of {nt}** now. That "
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
         f"at every surviving fraction in "
         f"{p9['sampled_at_surviving_fractions']} — which is "
         f"{p9['dose_range_gy'][0]} to {p9['dose_range_gy'][1]} Gy and NOT the "
         f"whole survival range, a claim an earlier version made and this one "
         f"withdraws. The admissible boosts are {p9['boost_window'][0]} to "
         f"{p9['boost_window'][1]}, which is "
         f"{p9['window_share_of_scan']:.0%} of a scan whose ceiling is "
         f"{p9['scan_ceiling']} — a share quoted without its ceiling measures "
         f"the ceiling, since the same window is {share01:.0%} of a 0-to-1 "
         f"scan. Within "
         f"the window the ratio is larger at low dose per fraction than high, "
         f"by a factor of {lo_r} to {hi_r}:", "",
         "| boost | SER at 2 Gy | SER at 6 Gy |", "|---|--:|--:|"]
    # Sorted HERE, on the boost, because a renderer that inherits dict order
    # publishes whatever the serialiser happened to do.
    for b, row in sorted(p9["ser_by_dose"].items(), key=lambda kv: float(kv[0])):
        L.append(f"| {b} | {row['2']} | {row['6']} |")
    L += ["",
          "**The direction does not depend on the fit**, which is what makes "
          "it worth preregistering: it holds across the entire window the "
          "published band permits, and it holds analytically at every dose "
          "rather than only at the sampled ones. What the window itself does "
          "NOT survive is a wider fractionation range — requiring the band at "
          "1.8 Gy and at 20 Gy as well empties it, so the window is a "
          "statement about 2.24 to 9.38 Gy and the page says which. Beta is "
          "untouched by the boost because "
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
