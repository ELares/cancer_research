#!/usr/bin/env python3
"""Fit each non-ferroptosis arm to its published target, and say what is left.

WHY THIS EXISTS
---------------
`analysis/modality-coverage.md` moved this project from thirteen of sixteen
taxonomy mechanisms absent to none, and `analysis/modality-panel.md` made every
arm produce a number. Both then said the same honest thing: **most of those
numbers are functions of parameters nobody fitted.** A mechanism the engine can
name, can apply, and has not calibrated is a long way from one it can answer a
question with, and `CALIBRATION_STATUS.md` records that per row as
`used in any reported number: N`.

This is the page that attacks that. Each arm has a NAMED target -- a published
number, cited on its constant, inside this project's own frozen corpus -- and
this fits the arm's one free parameter to it and reports what the fit costs.

WHAT A FIT IS AND IS NOT, because the distinction is the whole content
---------------------------------------------------------------------
Three outcomes are possible per arm and they are NOT the same result:

  ADMISSIBLE     a parameter value reproduces the target, and the range that
                 does is narrow enough that the data has said something.
  UNCONSTRAINED  a value exists, but so does most of the search range: the
                 target is satisfiable and uninformative, which is a fact
                 about the TARGET, not the model.
  INADMISSIBLE   no value in the range reproduces it. That falsifies the
                 FORM, and it is the outcome worth wanting, because it is the
                 only one that can teach the model something.

Reported per arm, with the admissible WIDTH as a fraction of the range
scanned, because "we fitted it" means nothing without that number.

WHAT THIS CANNOT DO
-------------------
**An ORR is not a kill fraction.** A RECIST objective response in a patient and
a dead cell in a lattice are different quantities, and no amount of fitting
makes them the same one. Where a target is clinical the fit is reported as
DIRECTION-ANCHORED: it pins the parameter to a value that reproduces the right
magnitude under an explicit and stated mapping, and the mapping is the weakest
link. Each row says which kind it is.

**One parameter per arm.** Fitting two would make most of these targets
underdetermined, and an underdetermined fit that reports a point estimate is
worse than no fit. Where an arm has more than one free parameter the others are
held at their documented defaults and the row says so.

Offline: reads the crate constants and the committed corpus citations. No
corpus scan, no simulation run.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CORE = REPO / "simulations" / "ferroptosis-core" / "src"

OUT_MD = REPO / "analysis" / "modality-calibration.md"
OUT_JSON = REPO / "analysis" / "modality-calibration.json"

# How wide a fitted window may be before the target has told us nothing. A
# window covering more than this share of the scanned range is reported
# UNCONSTRAINED rather than fitted -- the value exists, and so does most of
# the range, which is a fact about the target.
UNCONSTRAINED_WIDTH = 0.5


def _rust_const(path: Path, name: str) -> float:
    m = re.search(rf"\b{re.escape(name)}: f64 = ([0-9.]+)", path.read_text())
    if not m:
        raise SystemExit(f"{name} not found in {path}")
    return float(m.group(1))


def _rust_tuple(path: Path, name: str) -> tuple:
    m = re.search(rf"\b{re.escape(name)}: \(f64, f64\) = \(([0-9.]+), ([0-9.]+)\)",
                  path.read_text())
    if not m:
        raise SystemExit(f"{name} not found in {path}")
    return float(m.group(1)), float(m.group(2))


# ── the arms ─────────────────────────────────────────────────────────────

def _fit(predict, target_lo, target_hi, lo, hi, steps=4000):
    """Every parameter in [lo, hi] whose prediction lands in the target band."""
    ok = []
    for i in range(steps + 1):
        x = lo + (hi - lo) * i / steps
        y = predict(x)
        if y is None or not math.isfinite(y):
            continue
        if target_lo <= y <= target_hi:
            ok.append(x)
    if not ok:
        return None
    return {"lo": min(ok), "hi": max(ok),
            "width_fraction": (max(ok) - min(ok)) / (hi - lo)}


def _parp_arm() -> dict:
    """SER band -> alpha boost. A MODEL-INTERNAL target: the observable and
    the model output are the same quantity, so no mapping is needed."""
    rad = CORE / "radiation.rs"
    alpha = _rust_const(rad, "ALPHA_GBM_PARAMETERISATION_PER_GY")
    ratio = _rust_const(rad, "ALPHA_BETA_TUMOUR_GY")
    beta = alpha / ratio
    lo_band, hi_band = _rust_tuple(rad, "PARP_SER_BAND")

    def dose(a, b, L):
        return (-a + math.sqrt(a * a + 4 * b * L)) / (2 * b)

    def ser(boost, sf):
        L = -math.log(sf)
        return dose(alpha, beta, L) / dose(alpha * (1 + boost), beta, L)

    # ONE boost must satisfy the band across the whole survival range, which
    # is what makes this a joint constraint rather than a single equation.
    def worst_case(boost):
        vals = [ser(boost, sf) for sf in (0.5, 0.1, 0.01)]
        return vals[0] if not (lo_band <= min(vals) and max(vals) <= hi_band) else vals[1]

    def admissible(boost):
        vals = [ser(boost, sf) for sf in (0.5, 0.1, 0.01)]
        return min(vals) if (min(vals) >= lo_band and max(vals) <= hi_band) else None

    fit = _fit(admissible, lo_band, hi_band, 0.0, 2.0)
    return {
        "arm": "Synthetic lethality (PARP)",
        "parameter": "RadiationConfig::parp_alpha_boost",
        "target": f"sensitizer enhancement ratio {lo_band}-{hi_band}",
        "source": "PMID 35205750 (glioma cells)",
        "target_kind": "model-internal (the observable IS the model output)",
        "range_scanned": [0.0, 2.0],
        "fit": fit,
        "mapping": None,
    }


def _radiation_arm() -> dict:
    """SF2 -> alpha. A ROUND TRIP, and the report has already withdrawn the
    claim that it is a prediction."""
    rad = CORE / "radiation.rs"
    ratio = _rust_const(rad, "ALPHA_BETA_TUMOUR_GY")
    sf2 = 0.5445  # the published input, quoted on the constant's doc comment

    def predict(alpha):
        beta = alpha / ratio
        return math.exp(-(alpha * 2.0 + beta * 4.0))

    fit = _fit(predict, sf2 - 1e-4, sf2 + 1e-4, 0.0, 2.0, steps=200_000)
    return {
        "arm": "Radiation (DNA channel)",
        "parameter": "RadiationConfig::alpha_per_gy",
        "target": f"SF2 = {sf2} at alpha/beta = {ratio}",
        "source": "PMID 32307022 (glioblastoma, EUD-model input)",
        "target_kind": "round trip -- the paper DERIVED its third number from "
                       "the first two through this same model, so reproducing "
                       "it checks the implementation and not the biology",
        "range_scanned": [0.0, 2.0],
        "fit": fit,
        "mapping": None,
    }


def _immunotherapy_arm() -> dict:
    """Published anti-PD-1 ORR -> baseline antigenicity.

    This one comes out INADMISSIBLE, and the diagnosis is the result.
    """
    params = (CORE / "params.rs").read_text()

    def cfg(field):
        m = re.search(rf"\b{field}: ([0-9.]+),", params)
        return float(m.group(1)) if m else None

    dc_rate = cfg("dc_maturation_rate")
    prime = cfg("tcell_priming_rate")
    brake = cfg("pd1_brake")
    anti = cfg("anti_pd1_efficacy")
    default_kill = cfg("tcell_kill_rate")
    panel_kill = 0.02  # the panel's rate, chosen to stay off the survivor cap
    eff_brake = brake * (1.0 - anti)

    # `immune_cascade` with no ferroptotic death reduces to a PRODUCT:
    #   kill_fraction = dc_rate * prime * kill_rate * (1 - eff_brake) * a
    # so the target constrains the product and not either factor.
    def coefficient(kill_rate):
        return dc_rate * prime * kill_rate * (1.0 - eff_brake)

    def predict(antigenicity):
        return min(coefficient(panel_kill) * antigenicity, 1.0)

    # `baseline_antigenicity` is a FRACTION of the tumour presenting, so 1.0
    # is its own upper bound -- scanning past it would fit a number that
    # cannot mean anything.
    fit = _fit(predict, 0.20, 0.30, 0.0, 1.0)
    ceiling = predict(1.0)
    needed = 0.20 / coefficient(panel_kill)
    return {
        "arm": "Checkpoint blockade",
        "parameter": "ImmuneParams::baseline_antigenicity",
        "target": "anti-PD-1 monotherapy ORR 20-30%",
        "source": "PMID 31877341 (melanoma 30%, RCC 30%, NSCLC ~20%)",
        "target_kind": "clinical -- DIRECTION-ANCHORED only",
        "range_scanned": [0.0, 1.0],
        "fit": fit,
        "ceiling_at_full_presentation": ceiling,
        "antigenicity_required": needed,
        "diagnosis": (
            f"At the panel's `tcell_kill_rate = {panel_kill}`, the whole "
            f"cascade collapses to `kill = {coefficient(panel_kill):.4f} * "
            "antigenicity`, so even at FULL presentation "
            f"({ceiling:.1%}) it cannot reach the 20% floor -- it would need "
            f"an antigenicity of {needed:.2f}, and that parameter is a "
            "fraction bounded by 1. Raising the kill rate to the engine "
            f"default ({default_kill}) instead saturates the survivor cap "
            "immediately. **The target constrains the PRODUCT "
            "`antigenicity x kill_rate`, and neither factor is identifiable "
            "from it alone.** That is the same finding "
            "`analysis/identifiability-report.md` reports for the ferroptosis "
            "engine's rate constants, arriving independently in the immune "
            "layer, and it is why this row is INADMISSIBLE rather than "
            "fitted to a number that would look precise and mean nothing."),
        "mapping": "An objective response in a patient is mapped onto a KILL "
                   "FRACTION in a lattice. Those are different quantities and "
                   "no fit makes them the same one: a partial response is a "
                   "30% diameter reduction, roughly a 66% volume reduction, "
                   "and a responding patient is not a dead tumour. Even had "
                   "the fit succeeded, this mapping would be the weakest link "
                   "in the row.",
    }


def _oncolytic_arm() -> dict:
    """T-VEC durable response, treated vs control -> immunogenicity.

    The RATIO is the target rather than the treated arm alone, because 16%
    quoted without its 2% control invites the wrong baseline -- which is why
    the constant stores both.
    """
    immune = CORE / "immune.rs"
    treated, control = _rust_tuple(immune, "TVEC_DURABLE_RESPONSE")
    ratio = treated / control
    params = (CORE / "params.rs").read_text()
    kd = float(re.search(r"\bdc_activation_kd: ([0-9.]+),", params).group(1))

    def predict(immunogenicity):
        # The virus's ADVANTAGE over baseline is the ICD quality it adds, so
        # the modelled ratio is its DC activation against the baseline's.
        q = immunogenicity * 10.0     # the panel's LP-scale conversion
        base_q = 0.5                  # spontaneous death quality, documented
        if base_q <= 0:
            return None
        return (q / (q + kd)) / (base_q / (base_q + kd))

    fit = _fit(predict, ratio * 0.8, ratio * 1.2, 0.0, 1.0)
    return {
        "arm": "Oncolytic virus",
        "parameter": "oncolytic_lysis immunogenicity",
        "target": f"durable-response ratio {ratio:.1f}x "
                  f"({treated:.0%} treated vs {control:.0%} control)",
        "source": "PMID 27298410 (T-VEC, OPTiM)",
        "target_kind": "clinical RATIO -- less mapping-sensitive than a rate, "
                       "because a ratio of response rates and a ratio of DC "
                       "activations are both dimensionless",
        "range_scanned": [0.0, 1.0],
        "fit": fit,
        "mapping": "A durable-response RATIO is mapped onto a ratio of "
                   "DC-activation fractions. Both are dimensionless, which "
                   "makes this less fragile than the ORR row -- but it "
                   "assumes the clinical advantage is entirely immunological, "
                   "and T-VEC also lyses cells directly.",
    }


def _ablation_arm() -> dict:
    """Two published numbers that must reconcile, which is a stronger check
    than either alone.

    The CEM43 threshold for coagulation is conventionally 240 cumulative
    equivalent minutes. Separately, PMID 22180520 states that "increase in
    tissue temperature beyond 57-60 C leads to denaturation of cell proteins
    sufficient to cause coagulative necrosis" -- a TEMPERATURE, with no time
    attached, which reads as effectively immediate.

    Those are different kinds of statement about the same event, so the model
    has to make them agree: the exposure its threshold implies at 57 C must be
    SHORT. If it came out at hours the two numbers would be describing
    different phenomena and one of them would be wrong for this model.
    """
    rs = (CORE / "ablation.rs").read_text()
    ref = _rust_const(CORE / "ablation.rs", "CEM43_REFERENCE_C")
    m = re.search(r"fn default_cem43_threshold\(\) -> f64 \{\s*([0-9.]+)", rs)
    if not m:
        raise SystemExit("the CEM43 threshold default is gone from ablation.rs")
    default_threshold = float(m.group(1))
    necrosis_c = 57.0  # the low end of the published band

    def exposure_seconds(threshold):
        # cem43 = t * 0.5^(43 - T)  =>  t = threshold / 0.5^(43 - T)
        return threshold / (0.5 ** (ref - necrosis_c)) * 60.0

    # ADMISSIBLE if the implied exposure at the published necrosis temperature
    # is between a tenth of a second and a minute -- short enough to read as
    # "beyond 57 C causes necrosis" and long enough not to be instantaneous.
    fit = _fit(exposure_seconds, 0.1, 60.0, 1.0, 2000.0)
    return {
        "arm": "Ablation (thermal)",
        "parameter": "AblationConfig::cem43_threshold",
        "target": f"exposure at {necrosis_c:.0f} C must be seconds, not hours",
        "source": "PMID 22180520 (57-60 C -> coagulative necrosis) against the "
                  "CEM43 = 240 convention",
        "target_kind": "consistency between TWO published numbers of different "
                       "kinds -- a temperature with no time, and a thermal dose "
                       "with no temperature",
        "range_scanned": [1.0, 2000.0],
        "fit": fit,
        "default_threshold": default_threshold,
        "implied_exposure_s": exposure_seconds(default_threshold),
        "mapping": None,
    }


def _strip_rust_comments(src: str) -> str:
    """Reuse `modality_coverage`'s scanner rather than writing a second one.

    A regex here would be a second implementation of a thing this repository
    already got wrong once (nested block comments, `//` inside a string), and
    two scanners that drift is its own defect class.
    """
    spec = importlib.util.spec_from_file_location(
        "modality_coverage", REPO / "scripts" / "modality_coverage.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.strip_rust_comments(src)


def _rust_default_barriers() -> dict:
    """`AdoptiveBarriers::default()`, parsed from the crate.

    A drift guard of the same shape `validate_spheroid_structure.py` and
    `validate_trigger_wave.py` already use: the Python fit and the Rust module
    must not be able to disagree about the setting the published band was
    measured in.
    """
    # Comments are STRIPPED first, using the same scanner the coverage
    # generators use. Without it the field regex reads values out of prose: two
    # explanatory comment lines inside the `Default` impl mentioning
    # `activation: 0.5` and `antigen_positive_fraction: 0.8` -- rustfmt-stable,
    # `cargo test` green -- made this return the solid-tumour values and flipped
    # the published CAR-T verdict from ADMISSIBLE to UNCONSTRAINED. The dict
    # comprehension below is last-wins, so a comment AFTER the real field silently
    # replaced it. Every other rustfmt-legal reformatting (struct-update syntax,
    # a const, `1.0f64`, `1e0`, delegation) already failed loudly; this one did not.
    src = _strip_rust_comments((CORE / "adoptive.rs").read_text())
    body = re.search(r"impl Default for AdoptiveBarriers \{.*?\n    \}",
                     src, re.S)
    if body is None:
        raise SystemExit("adoptive.rs: no Default impl to read barriers from")
    out = {k: float(v) for k, v in
           re.findall(r"(\w+): ([0-9.]+),", body.group(0))}
    need = {"trafficking", "infiltration", "activation",
            "exhaustion_rate", "antigen_positive_fraction"}
    missing = need - set(out)
    if missing:
        raise SystemExit(f"adoptive.rs Default is missing {sorted(missing)}")
    return out


def _cart_arm() -> dict:
    """Published B-ALL complete-remission band -> effector supply.

    A CLINICAL target, and the source carries its own caveat which is quoted
    beside the constant in `immune.rs` and repeated here.
    """
    params = (CORE / "params.rs").read_text()
    immune_rs = (CORE / "immune.rs").read_text()
    lo, hi = _rust_tuple(CORE / "immune.rs", "CART_B_ALL_CR_BAND")

    def cfg(field):
        m = re.search(rf"\b{field}: ([0-9.]+),", params)
        return float(m.group(1)) if m else None

    brake, anti = cfg("pd1_brake"), cfg("anti_pd1_efficacy")
    kill = 0.02
    n = 20_000
    eff_brake = brake * (1.0 - anti)

    # The barriers this fit runs at, READ FROM `adoptive.rs::Default` rather
    # than assumed. For one review round the module documented this default as
    # load-bearing for this fit while the fit reimplemented the kill in Python
    # and never opened the file -- a claim about one artifact written in
    # another, which is this repository's most-repeated defect and was found
    # here by a reviewer tracing callers instead of grepping for the type.
    # Now a solid-tumour default would move the fitted band, and
    # `test_modality_calibration.py` fails if these stop being all-open.
    barriers = _rust_default_barriers()
    delivery = (barriers["trafficking"] * barriers["infiltration"]
                * barriers["activation"])
    ceiling = barriers["antigen_positive_fraction"]

    def predict(effectors):
        # `adoptive_transfer_kills` composed with `adoptive::` -- the leukaemia
        # setting (no suppression, every barrier open), which is what the band
        # is measured in. `barrier_limited_kills` is the antigen cap.
        raw = effectors * delivery * kill * (1.0 - eff_brake)
        return min(raw, float(n) * ceiling) / n

    fit = _fit(predict, lo, hi, 0.0, 5_000_000.0, steps=20_000)
    return {
        "arm": "CAR-T (adoptive transfer)",
        "parameter": "adoptive_transfer_kills effector count",
        "target": f"B-ALL complete remission {lo:.0%}-{hi:.0%}",
        "source": "PMID 32607912 (CD19 CAR-T)",
        "target_kind": "clinical -- DIRECTION-ANCHORED only",
        "range_scanned": [0.0, 5_000_000.0],
        "fit": fit,
        "mapping": "A complete remission is mapped onto a KILL FRACTION. The "
                   "same source says the result has NOT transferred to solid "
                   "tumours and that 30-50% of responses were not durable, so "
                   "a fit to the headline band describes the indication these "
                   "therapies were approved for and not the setting this "
                   "engine simulates. The fit is to the leukaemia case with "
                   "suppression at zero, which is the only setting the band "
                   "is measured in.",
    }


def _adc_arm() -> dict:
    """No fittable published endpoint in this corpus, reported as such.

    This is the fourth outcome and the one a calibration page is most tempted
    to hide. The ADC module's mechanism is corpus-anchored on BOTH arms of its
    own comparison, which is why the layer landed -- but an anchor for a
    MECHANISM is not a target for a PARAMETER, and nothing in the frozen
    corpus gives a number `payload_escape_fraction` or `neighbours_in_reach`
    could be fitted to.
    """
    return {
        "arm": "ADC bystander effect",
        "parameter": "AdcConfig::payload_escape_fraction",
        "target": None,
        "source": "PMID 31930187 anchors the MECHANISM, not a rate",
        "target_kind": "no fittable endpoint in this corpus",
        "range_scanned": None,
        "fit": None,
        "mapping": None,
        "no_target_reason": (
            "The corpus establishes that a cleavable linker causes bystander "
            "killing and a non-cleavable one does not -- a qualitative "
            "contrast with both arms named, which is what let the layer land "
            "under the layer-freeze policy. It gives no escape fraction, no "
            "diffusion radius and no bystander kill rate. Reporting a fit "
            "here would mean inventing a target, and an invented target that "
            "a flexible form then satisfies is worse than no target at all: "
            "it looks like calibration and constrains nothing."),
    }


def scan() -> dict:
    return {"arms": [_radiation_arm(), _parp_arm(), _immunotherapy_arm(),
                     _oncolytic_arm(), _ablation_arm(), _cart_arm(),
                     _adc_arm()],
            "unconstrained_width": UNCONSTRAINED_WIDTH}


def _verdict(a: dict, cap: float) -> str:
    # A row with no TARGET is a different outcome from one whose target no
    # value satisfies, and collapsing them would hide the more interesting
    # admission: that this project has no number to fit.
    if a.get("target") is None:
        return "NO TARGET"
    f = a["fit"]
    if f is None:
        return "INADMISSIBLE"
    if f["width_fraction"] > cap:
        return "UNCONSTRAINED"
    return "ADMISSIBLE"


def assemble(raw: dict) -> dict:
    cap = raw["unconstrained_width"]
    arms = [dict(a, verdict=_verdict(a, cap)) for a in raw["arms"]]
    counts = {v: sum(1 for a in arms if a["verdict"] == v)
              for v in ("ADMISSIBLE", "UNCONSTRAINED", "INADMISSIBLE",
                        "NO TARGET")}
    return dict(raw, arms=arms, verdict_counts=counts,
                n_clinical=sum(1 for a in arms if a["mapping"]))


def render(d: dict) -> str:
    arms = d["arms"]
    c = d["verdict_counts"]
    L = ["# Fitting the non-ferroptosis arms to their published targets", "",
         "*Generated by `scripts/calibrate_modality_arms.py --render-only`. "
         "Offline; reads the crate constants. No simulation run.*", "",
         "Every arm this project added carries a NAMED target — a published "
         "number, cited on its constant, inside the frozen corpus. Naming one "
         "is what `CONTRIBUTING.md`'s layer-freeze policy requires. FITTING to "
         "it is what turns a layer from present into calibrated, and this is "
         "where that is attempted and where it is reported honestly when it "
         "does not work.", "",
         "| arm | parameter | target | verdict | fitted range | width |",
         "|---|---|---|---|---|--:|"]
    for a in arms:
        f = a["fit"]
        rng = (f"{f['lo']:.4g} – {f['hi']:.4g}" if f else "—")
        tgt = a["target"] or "*none in this corpus*"
        wid = (f"{f['width_fraction'] * 100:.0f}%" if f else "—")
        L.append(f"| {a['arm']} | `{a['parameter']}` | {tgt} | "
                 f"**{a['verdict']}** | {rng} | {wid} |")
    L += ["",
          f"**{c['ADMISSIBLE']} admissible, {c['UNCONSTRAINED']} "
          f"unconstrained, {c['INADMISSIBLE']} inadmissible, "
          f"{c['NO TARGET']} with no target at all.** The width "
          "column is the point: a fit that admits most of the search range "
          "has been given a target that cannot discriminate, and reporting it "
          "as \"calibrated\" would be the same error as reporting a "
          "p-value without an effect size.", ""]

    L += ["## Four outcomes, and they are not the same result", "",
          "**ADMISSIBLE** — a value reproduces the target and the range that "
          "does is narrow. The data has said something.", "",
          "**UNCONSTRAINED** — a value exists and so does most of the search "
          "range. That is a fact about the TARGET, not the model, and it is "
          "reported rather than dressed up: a satisfiable target that "
          "excludes nothing has not calibrated anything.", "",
          "**INADMISSIBLE** — nothing in the range reproduces it, which "
          "falsifies the FORM. **This is the outcome worth wanting**, because "
          "it is the only one that can teach the model something.", "",
          "**NO TARGET** — the corpus anchors the MECHANISM and gives no "
          "number the parameter could be fitted to. This is the outcome a "
          "calibration page is most tempted to hide, and the temptation is "
          "specific: inventing a target that a flexible form then satisfies "
          "looks exactly like calibration and constrains nothing. A row that "
          "says it has no number is more use than a row that has made one "
          "up.", ""]

    none = [a for a in arms if a["verdict"] == "NO TARGET"]
    if none:
        L += ["## The rows with nothing to fit to", ""]
        for a in none:
            L += [f"**{a['arm']}** — {a.get('no_target_reason', 'no reason recorded')}",
                  ""]

    unc = [a for a in arms if a["verdict"] == "UNCONSTRAINED"]
    if unc:
        L += ["## The rows whose target excludes almost nothing", ""]
        for a in unc:
            f = a["fit"]
            L += [f"**{a['arm']}** — {f['width_fraction'] * 100:.0f}% of the "
                  "searched range satisfies the target, so the target is "
                  "consistent with the model and constrains it barely at all. "
                  "That is a fact about the target and it is reported rather "
                  "than counted as a fit.", ""]

    bad = [a for a in arms if a["verdict"] == "INADMISSIBLE"]
    if bad:
        L += ["## The inadmissible row is the informative one", ""]
        for a in bad:
            L += [f"**{a['arm']}** — {a.get('diagnosis', 'no diagnosis recorded')}",
                  ""]
        L += ["A row like this is worth more than the three that fitted. The "
              "fitted ones show that a flexible form can be made to hit a "
              "number, which is weak evidence; this one says something about "
              "the model that was not put in by hand.", ""]

    clinical = [a for a in arms if a["mapping"]]
    if clinical:
        L += ["## The mapping is the weak link, and it is stated per row", "",
              f"{len(clinical)} of the {len(arms)} targets are CLINICAL, and "
              "no amount of fitting makes a clinical endpoint and a lattice "
              "kill fraction the same quantity.", ""]
        for a in clinical:
            L += [f"**{a['arm']}** — {a['mapping']}", ""]

    L += ["## What this does not say", "",
          "**A fitted parameter is not a validated model.** Every fit here "
          "moves ONE parameter with the rest held at documented defaults, "
          "because fitting two would leave most of these targets "
          "underdetermined and an underdetermined fit reporting a point "
          "estimate is worse than no fit at all. What a row shows is that the "
          "form CAN reproduce the target, not that it does so for the right "
          "reason.", "",
          "**Hitting a target is weak evidence when nothing could miss it.** "
          "The radiation row is the clearest case and the report already "
          "withdrew the stronger claim: the paper DERIVED its third number "
          "from the first two through this same model, so reproducing it "
          "checks the implementation and not the biology.", "",
          "**None of these fits is wired into a binary.** They constrain "
          "parameters; they do not change any committed number. "
          "`simulations/calibration/CALIBRATION_STATUS.md` remains the "
          "authority on what is used where, and these rows are still "
          "`used in any reported number: N`.", ""]
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
    for arm in d["arms"]:
        f = arm["fit"]
        w = f"{f['width_fraction'] * 100:5.1f}%" if f else "   --"
        print(f"  {arm['arm']:28s} {arm['verdict']:14s} width {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
