#!/usr/bin/env python3
"""Anchor the tumor-PK model to measured ferroptosis-inducer pharmacokinetics (#334).

The repo's `tumor_pk` ODE (plasma -> vascular -> interstitial) and its per-tumor
presets are currently order-of-magnitude ESTIMATES, not anchored to any measured
drug (see `simulations/calibration/CALIBRATION_STATUS.md`, "RSL3 pharmacokinetics:
Uncalibrated"). This script anchors the plasma + tumor disposition to the only
PUBLIC, ferroptosis-specific in-vivo PK dataset that reports a PAIRED plasma AND
tumor concentration-time course: imidazole ketone erastin (IKE), a system-xc-
inhibitor developed for in-vivo stability (Zhang et al. 2019, Cell Chem Biol,
PMID 30799221). Sorafenib (a clinical kinase inhibitor that also induces
ferroptosis via system-xc-) provides a second, human-scale anchor from a published
population-PK model (Jain et al. 2011, PMID 21392074).

WHAT IS ANCHORED
----------------
A minimal two-compartment plasma -> tumor model:

  plasma:  Cp(t) = scale * (exp(-ke*t) - exp(-ka*t))     # 1-cmt, first-order absorption
  tumor:   dCt/dt = k_pt * Cp(t) - k_te * Ct,  Ct(0)=0   # tumor tissue compartment

fit to the IKE summary NCA targets (Tmax, Cmax, terminal half-life, AUC) for each
compartment. The tumor transfer/elimination rates are DERIVED from two purely
measured quantities, then the tumor Tmax and Cmax are left as out-of-fit
PREDICTIONS validated against the measured values:

  * k_te = ln2 / t_half_tumor                  (measured tumor terminal half-life)
  * k_pt = Kp * k_te,  Kp = AUC_tumor/AUC_plasma   (measured tissue:plasma partition)

The headline anchor is the measured partition Kp ~ 0.90 and the measured plasma ->
tumor delay (plasma Tmax 1.35 h -> tumor Tmax 3.30 h), both of which the model's
uncalibrated `partition_coeff` (0.15 to 0.5 in the presets) and transfer rates can
now be checked against.

HONESTY / SCOPE
---------------
- This anchors the PK model STRUCTURE and its partition/delay to measured data for
  two ferroptosis-relevant drugs. It does NOT recalibrate every per-tumor-type
  preset (breast/pancreatic/GBM/...): no public per-tumor-type measured PK exists
  for a ferroptosis inducer, so that stays a documented gap.
- The canonical tool compounds the biochem layer uses (RSL3, ML162, ML210) have
  no usable published in-vivo PK (erastin's poor metabolic stability is exactly why
  IKE was engineered). IKE is therefore the in-vivo system-xc- proxy; RSL3/ML PK is
  flagged absent rather than fabricated.
- A single 1-compartment absorption model cannot match all four plasma summary
  metrics simultaneously (the real disposition is multi-exponential); the residual
  is reported, not hidden.
- This is an in-vivo MOUSE anchor (IKE) plus a human popPK anchor (sorafenib); it is
  not a whole-body multi-organ PBPK build.

Run (pure Python + scipy; no compiled extension; CI-safe):
  python3 scripts/calibrate_pk.py
Writes analysis/calibration/pk-calibration.md + .json.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

from scipy.optimize import brentq

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_CSV = REPO_ROOT / "analysis" / "calibration" / "pk_measured_data.csv"
OUT_MD = REPO_ROOT / "analysis" / "calibration" / "pk-calibration.md"
OUT_JSON = REPO_ROOT / "analysis" / "calibration" / "pk-calibration.json"

LN2 = math.log(2.0)

# IKE distribution study: the paired plasma+tumor cohort we fit (NCG SUDHL6 xenograft).
IKE_POP = "NCG SUDHL6 xenograft"


# ---------------------------------------------------------------------------
# Pure PK model (no scipy; unit-tested directly)
# ---------------------------------------------------------------------------

def plasma_conc(t, ka, ke, scale):
    """One-compartment first-order absorption concentration at time t (t >= 0)."""
    return scale * (math.exp(-ke * t) - math.exp(-ka * t))


def plasma_tmax(ka, ke):
    """Analytical time of peak for the absorption model (ka != ke)."""
    return math.log(ka / ke) / (ka - ke)


def plasma_auc(ka, ke, scale):
    """Analytical AUC(0..inf) of the absorption model."""
    return scale * (1.0 / ke - 1.0 / ka)


def tumor_conc(t, ka, ke, scale, k_pt, k_te):
    """Analytical tumor concentration: linear ODE driven by the plasma curve.

    dCt/dt = k_pt*Cp(t) - k_te*Ct, Ct(0)=0, with Cp = scale*(exp(-ke t)-exp(-ka t)).
    Falls back to a numerical solve if any rate pair is near-degenerate.
    """
    if min(abs(k_te - ke), abs(k_te - ka)) < 1e-6:
        return tumor_conc_numeric(t, ka, ke, scale, k_pt, k_te)
    term_e = (math.exp(-ke * t) - math.exp(-k_te * t)) / (k_te - ke)
    term_a = (math.exp(-ka * t) - math.exp(-k_te * t)) / (k_te - ka)
    return k_pt * scale * (term_e - term_a)


def tumor_conc_numeric(t_end, ka, ke, scale, k_pt, k_te, dt=0.001):
    """Fixed-step RK4 of the tumor ODE to t_end (cross-check for tumor_conc)."""
    def deriv(tt, ct):
        return k_pt * plasma_conc(tt, ka, ke, scale) - k_te * ct
    n = max(1, int(round(t_end / dt)))
    h = t_end / n
    ct = 0.0
    tt = 0.0
    for _ in range(n):
        k1 = deriv(tt, ct)
        k2 = deriv(tt + h / 2, ct + h / 2 * k1)
        k3 = deriv(tt + h / 2, ct + h / 2 * k2)
        k4 = deriv(tt + h, ct + h * k3)
        ct += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        tt += h
    return ct


def tumor_auc(scale, ka, ke, k_pt, k_te):
    """AUC(0..inf) of the tumor compartment = (k_pt/k_te) * plasma AUC (mass balance)."""
    return (k_pt / k_te) * plasma_auc(ka, ke, scale)


def argmax_on_grid(fn, t_hi, n=4000):
    """Time and value of the maximum of fn(t) on a fine grid over [0, t_hi]."""
    best_t, best_v = 0.0, fn(0.0)
    for i in range(1, n + 1):
        t = t_hi * i / n
        v = fn(t)
        if v > best_v:
            best_v, best_t = v, t
    return best_t, best_v


def terminal_half_life(fn, t0, t1):
    """Apparent terminal half-life from two late-time points on a curve."""
    c0, c1 = fn(t0), fn(t1)
    if c0 <= 0 or c1 <= 0 or c1 >= c0:
        return float("nan")
    k = (math.log(c0) - math.log(c1)) / (t1 - t0)
    return LN2 / k


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

def ka_from_tmax(tmax, ke):
    """Solve ln(ka/ke)/(ka-ke) = tmax for ka > ke (the absorption root).

    plasma_tmax decreases from its supremum 1/ke (as ka->ke+) toward 0 (as
    ka->inf), so a 1-compartment first-order-absorption model can only realize a
    peak time tmax < 1/ke. If tmax >= 1/ke the equation has no ka > ke root, which
    we raise on rather than silently returning a wrong ka (a model not actually
    fitted to the requested peak time).
    """
    f = lambda ka: plasma_tmax(ka, ke) - tmax
    lo, hi = ke * (1.0 + 1e-4), ke * 1.0e4
    if f(lo) * f(hi) > 0:
        raise ValueError(
            f"tmax={tmax:.4f} h is not achievable for a 1-compartment "
            f"first-order-absorption model with ke={ke:.4f}/h: the peak time must "
            f"be below 1/ke={1.0 / ke:.4f} h."
        )
    return brentq(f, lo, hi, xtol=1e-9)


def one_cmt_auc_over_cmax_floor(tmax):
    """Minimum achievable AUC/Cmax for ANY 1-compartment first-order-absorption curve
    at a fixed Tmax. Feasibility requires Tmax < 1/ke. The normalized exposure
    AUC/(Cmax*Tmax) is minimized at `e` in the equal-rate limit (the absorption rate
    constant ka approaching the elimination rate constant ke) and grows without bound
    as ka -> inf (the instantaneous-absorption/bolus limit), so AUC/Cmax >= `e * Tmax`
    with equality reached as ka -> ke. A measured AUC/Cmax below this floor is proof
    that the disposition is multi-compartment (a fast distribution phase)."""
    return math.e * tmax


def fit_plasma(tmax, cmax, thalf, auc):
    """Anchor (ka, ke, scale) EXACTLY to the three shape-defining plasma metrics
    (Tmax, Cmax, terminal half-life): `ke = ln2/thalf`, `ka` solved from `Tmax`,
    `scale` from `Cmax`. AUC is then a PREDICTION; for IKE the measured AUC/Cmax
    (2.11 h) is below the hard 1-compartment floor e*Tmax (~3.67 h), so a 1-cmt model
    necessarily OVER-predicts AUC. That over-prediction is reported as the evidence
    of an omitted fast distribution phase, not papered over. Returns (params,
    predicted metrics, relative residuals); Tmax/Cmax/thalf are ~0 by construction.
    """
    ke = LN2 / thalf
    ka = ka_from_tmax(tmax, ke)
    if ka <= ke:
        ka = ke * 1.0001
    tm = plasma_tmax(ka, ke)
    scale = cmax / max(plasma_conc(tm, ka, ke, 1.0), 1e-12)
    params = {"ka": ka, "ke": ke, "scale": scale}
    pred = {
        "Tmax": plasma_tmax(ka, ke),
        "Cmax": plasma_conc(plasma_tmax(ka, ke), ka, ke, scale),
        "thalf": LN2 / ke,
        "AUC": plasma_auc(ka, ke, scale),
    }
    rel = {k: pred[k] / m - 1.0 for k, m in (("Tmax", tmax), ("Cmax", cmax), ("thalf", thalf), ("AUC", auc))}
    return params, pred, rel


def derive_tumor(plasma_params, kp_measured, thalf_tumor):
    """Tumor rates from purely measured quantities (partition + terminal half-life).

    k_te from the measured tumor terminal half-life; k_pt from the measured
    tissue:plasma AUC ratio (mass balance: AUC_t/AUC_p = k_pt/k_te). Tumor Tmax and
    Cmax then follow as PREDICTIONS to validate against measured values.
    """
    ka, ke, scale = plasma_params["ka"], plasma_params["ke"], plasma_params["scale"]
    k_te = LN2 / thalf_tumor
    k_pt = kp_measured * k_te
    fn = lambda t: tumor_conc(t, ka, ke, scale, k_pt, k_te)
    t_hi = max(24.0, 8.0 * thalf_tumor)
    tmax_pred, cmax_pred = argmax_on_grid(fn, t_hi)
    th_pred = terminal_half_life(fn, 0.7 * t_hi, 0.9 * t_hi)
    return {
        "k_pt": k_pt,
        "k_te": k_te,
        "pred": {"Tmax": tmax_pred, "Cmax": cmax_pred, "thalf": th_pred,
                 "AUC": tumor_auc(scale, ka, ke, k_pt, k_te)},
    }


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_targets(path=DATA_CSV):
    """Nested dict: data[drug][population][route][compartment][parameter] = value."""
    data = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            (data.setdefault(r["drug"], {}).setdefault(r["population"], {})
                 .setdefault(r["route"], {}).setdefault(r["compartment"], {}))[r["parameter"]] = float(r["value"])
    return data


def _nca(d):
    return d["Tmax"], d["Cmax"], d["thalf"], d["AUC"]


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run(args):
    data = load_targets(args.data)
    ike = data["IKE"][IKE_POP]["IP"]
    p_tmax, p_cmax, p_thalf, p_auc = _nca(ike["plasma"])
    t_tmax, t_cmax, t_thalf, t_auc = _nca(ike["tumor"])
    kp_measured = t_auc / p_auc

    # Closed-form multi-compartment finding: a 1-cmt absorption curve has a hard
    # floor AUC/Cmax >= e*Tmax. A measured ratio below it proves a fast distribution
    # phase (multi-compartment disposition the simple model omits).
    auc_over_cmax = p_auc / p_cmax
    floor = one_cmt_auc_over_cmax_floor(p_tmax)
    multicompartment = {
        "plasma_auc_over_cmax_h": round(auc_over_cmax, 3),
        "one_cmt_floor_e_tmax_h": round(floor, 3),
        "below_one_cmt_floor": auc_over_cmax < floor,
        "shortfall_fraction": round(1.0 - auc_over_cmax / floor, 3),
    }

    plasma_params, plasma_pred, plasma_rel = fit_plasma(p_tmax, p_cmax, p_thalf, p_auc)
    tumor = derive_tumor(plasma_params, kp_measured, t_thalf)

    # Out-of-fit tumor prediction errors (Tmax/Cmax were NOT fit).
    tumor_meas = {"Tmax": t_tmax, "Cmax": t_cmax, "thalf": t_thalf, "AUC": t_auc}
    tumor_rel = {k: tumor["pred"][k] / tumor_meas[k] - 1.0 for k in tumor_meas}

    # Exported plasma + tumor curves (normalized to plasma Cmax) for a measured-anchored
    # DoseSchedule::FromPk series the spatial sims can consume.
    ka, ke, scale = plasma_params["ka"], plasma_params["ke"], plasma_params["scale"]
    norm = plasma_pred["Cmax"]
    grid = [round(0.25 * i, 2) for i in range(0, 49)]  # 0..12 h, 0.25 h step
    plasma_curve = [round(plasma_conc(t, ka, ke, scale) / norm, 5) for t in grid]
    tumor_curve = [round(tumor_conc(t, ka, ke, scale, tumor["k_pt"], tumor["k_te"]) / norm, 5) for t in grid]

    # Sorafenib human-scale anchor (forward sanity check from published popPK).
    sor = data["sorafenib"]["solid-tumor popPK"]["PO"]["central"]
    cl_f, v_f = sor["CL_F"], sor["V_F"]
    ke_sor = cl_f / v_f
    tau = 12.0  # 400 mg BID
    dose_mg = 400.0
    cavg_ss_mg_l = dose_mg / (cl_f * tau)  # = Dose/(CL/F * tau), apparent (per F)
    sorafenib = {
        "CL_F_L_per_h": cl_f,
        "V_F_L": v_f,
        "derived_terminal_thalf_h": round(LN2 / ke_sor, 2),
        "cavg_ss_mg_per_L_400mg_BID": round(cavg_ss_mg_l, 3),
        "clinical_cavg_ss_range_mg_per_L": [3.0, 5.0],
        "in_clinical_range": 3.0 <= cavg_ss_mg_l <= 5.0,
    }

    result = {
        "drug_primary": "IKE",
        "drug_primary_pmid": 30799221,
        "population": IKE_POP,
        "route": "IP",
        "dose_mg_per_kg": 50,
        "measured": {
            "plasma": {"Tmax": p_tmax, "Cmax": p_cmax, "thalf": p_thalf, "AUC": p_auc},
            "tumor": tumor_meas,
            "partition_kp_auc_ratio": round(kp_measured, 4),
            "plasma_to_tumor_tmax_delay_h": round(t_tmax - p_tmax, 3),
        },
        "multicompartment_finding": multicompartment,
        "plasma_fit": {
            "params": {k: round(v, 6) for k, v in plasma_params.items()},
            "predicted": {k: round(v, 3) for k, v in plasma_pred.items()},
            "rel_residual": {k: round(v, 4) for k, v in plasma_rel.items()},
            "max_abs_rel_residual": round(max(abs(v) for v in plasma_rel.values()), 4),
        },
        "tumor_derived": {
            "k_pt_per_h": round(tumor["k_pt"], 5),
            "k_te_per_h": round(tumor["k_te"], 5),
            "predicted": {k: round(v, 3) for k, v in tumor["pred"].items()},
            "rel_residual_out_of_fit": {k: round(v, 4) for k, v in tumor_rel.items()},
            "tmax_pred_h": round(tumor["pred"]["Tmax"], 3),
            "tmax_meas_h": t_tmax,
        },
        "sorafenib_anchor": sorafenib,
        "frompk_series": {
            "time_h": grid,
            "plasma_norm": plasma_curve,
            "tumor_norm": tumor_curve,
            "note": ("normalized to fitted plasma Cmax; 1-cmt shape (Cmax/Tmax/terminal-slope "
                     "accurate, terminal over-sustained because the fast distribution phase is "
                     "omitted); usable as DoseSchedule::FromPk (#239). Rescale to measured AUC "
                     "for exposure-sensitive uses."),
        },
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_report(result)
    print(f"IKE measured Kp (tumor/plasma AUC) = {result['measured']['partition_kp_auc_ratio']}")
    print(f"plasma AUC/Cmax {multicompartment['plasma_auc_over_cmax_h']} h < 1-cmt floor "
          f"{multicompartment['one_cmt_floor_e_tmax_h']} h => multi-compartment "
          f"(below floor = {multicompartment['below_one_cmt_floor']})")
    print(f"plasma 1-cmt AUC over-prediction = {result['plasma_fit']['rel_residual']['AUC']:+.0%} "
          f"(Tmax/Cmax/thalf exact)")
    print(f"tumor Tmax predicted {result['tumor_derived']['tmax_pred_h']} h vs measured {t_tmax} h "
          f"(out-of-fit rel err {result['tumor_derived']['rel_residual_out_of_fit']['Tmax']})")
    print(f"sorafenib Cavg,ss(400 BID) = {sorafenib['cavg_ss_mg_per_L_400mg_BID']} mg/L "
          f"(clinical 3 to 5; in range = {sorafenib['in_clinical_range']})")
    print(f"wrote {OUT_MD.relative_to(REPO_ROOT)} + {OUT_JSON.relative_to(REPO_ROOT)}")
    return result


def write_report(r):
    m = r["measured"]
    mc = r["multicompartment_finding"]
    pf = r["plasma_fit"]
    td = r["tumor_derived"]
    sa = r["sorafenib_anchor"]

    def row(name, meas, pred, rel):
        return f"| {name} | {meas} | {pred} | {rel:+.1%} |"

    plasma_rows = "\n".join(
        row(k, m["plasma"][k], pf["predicted"][k], pf["rel_residual"][k])
        for k in ("Tmax", "Cmax", "thalf", "AUC")
    )
    tumor_rows = "\n".join(
        row(k + (" (predicted)" if k in ("Tmax", "Cmax") else " (derived)"),
            m["tumor"][k], td["predicted"][k], td["rel_residual_out_of_fit"][k])
        for k in ("Tmax", "Cmax", "thalf", "AUC")
    )

    md = f"""# Tumor-PK anchored to measured ferroptosis-inducer PK (#334)

Generated by `scripts/calibrate_pk.py` (pure Python + scipy; runs in CI). Target
data: `analysis/calibration/pk_measured_data.csv` (see `calibration-targets-pk.md`).

## What this anchors

The `tumor_pk` ODE and its per-tumor presets were order-of-magnitude estimates,
not tied to any measured drug (`CALIBRATION_STATUS.md`: "RSL3 pharmacokinetics,
Uncalibrated"). This leg anchors the plasma + tumor disposition to imidazole
ketone erastin (**IKE**, a system-xc- ferroptosis inducer engineered for in-vivo
stability; Zhang 2019, PMID 30799221), the only public ferroptosis-specific
dataset with a PAIRED plasma + tumor concentration-time course, and adds a
human-scale anchor from a published sorafenib population-PK model (Jain 2011,
PMID 21392074).

## Headline measured anchors (IKE, IP 50 mg/kg, {r['population']})

- **Tissue:plasma partition Kp = {m['partition_kp_auc_ratio']}** (tumor AUC / plasma AUC).
  The presets' `partition_coeff` is currently 0.15 to 0.5 (estimated); the measured
  value for this system-xc- inducer is ~0.90, i.e. tumor exposure nearly matches
  plasma over the dosing window.
- **Plasma to tumor delay = {m['plasma_to_tumor_tmax_delay_h']} h** (plasma Tmax
  {m['plasma']['Tmax']} h, tumor Tmax {m['tumor']['Tmax']} h). The tumor peaks later
  and clears slower (tumor terminal half-life {m['tumor']['thalf']} h vs plasma
  {m['plasma']['thalf']} h), the signature of a distributional tissue compartment.

## A 1-compartment plasma model is provably insufficient (structural finding)

Before fitting, an exact arithmetic check on the IKE plasma NCA: a 1-compartment
first-order-absorption curve has a hard lower bound on its total exposure per peak,
`AUC/Cmax >= e * Tmax` (the floor is reached in the equal-rate ka->ke limit, not the instantaneous-absorption limit).
For IKE the measured ratio is **{mc['plasma_auc_over_cmax_h']} h**, which is
**{mc['shortfall_fraction']:.0%} below** the floor `e * Tmax = {mc['one_cmt_floor_e_tmax_h']} h`.
That is impossible for any single compartment, so IKE has a fast distribution phase
that clears most drug shortly after the peak: the disposition is genuinely
multi-compartment. This is itself a result. It explains why the repo's `tumor_pk`
is a two-compartment (vascular to interstitial) ODE rather than a single
compartment, and it bounds what the summary NCA alone can identify (a unique multi-compartment plasma
model is NOT identifiable from four summary numbers, so we anchor the robust,
identifiable quantities and report the rest as flagged residuals rather than
over-fitting).

## Plasma curve anchored to Tmax, Cmax, and terminal half-life

The plasma curve is anchored EXACTLY to its three shape-defining metrics
(`ke = ln2/thalf`, `ka` from Tmax, `scale` from Cmax): `ka = {pf['params']['ka']:.4f}/h`,
`ke = {pf['params']['ke']:.4f}/h`, `scale = {pf['params']['scale']:.1f}` ng/mL.

| plasma metric | measured | model | rel. error |
|---|---|---|---|
{plasma_rows}

Tmax, Cmax, and terminal half-life match by construction. The AUC is then a
PREDICTION, and the 1-compartment curve necessarily over-predicts it (by
{pf['rel_residual']['AUC']:+.0%}) because it has no fast distribution phase to remove
the early exposure: the same finding as above, now quantified. The curve is
therefore accurate in peak level, peak time, and terminal slope, and over-sustained
in the middle, which is the documented limitation of the exported `frompk_series`.

## Tumor compartment (rates from measured partition + half-life, Tmax/Cmax predicted)

The tumor transfer and elimination rates are DERIVED from two purely measured
quantities, so the tumor **Tmax and Cmax are genuine out-of-fit predictions**:

- `k_te = ln2 / t_half_tumor = {td['k_te_per_h']}/h` (measured tumor half-life).
- `k_pt = Kp * k_te = {td['k_pt_per_h']}/h` (measured tissue:plasma AUC ratio).

| tumor metric | measured | model | rel. error |
|---|---|---|---|
{tumor_rows}

The model captures the DIRECTION of the plasma to tumor delay (predicted tumor Tmax
{td['tmax_pred_h']} h vs measured {td['tmax_meas_h']} h) and the slower tumor
clearance (tumor terminal half-life matches by construction). The Tmax magnitude is
over-predicted ({td['rel_residual_out_of_fit']['Tmax']:+.0%}). The tumor AUC inherits
the plasma limitation EXACTLY: by the mass-balance identity tumor_AUC = Kp * plasma_AUC,
the {pf['rel_residual']['AUC']:+.0%} plasma AUC over-prediction propagates to the same
{td['rel_residual_out_of_fit']['AUC']:+.0%} tumor AUC over-prediction, and the tumor
Tmax is pushed late because the over-sustained plasma input feeds the tumor too long.
Both are the same omitted fast distribution phase, not independent errors. The robust
results are the measured partition, the slower tumor clearance, and the delay
direction, not the Tmax to the hour.

## Sorafenib human-scale anchor (forward sanity check)

From the published popPK (Jain 2011): `CL/F = {sa['CL_F_L_per_h']} L/h`,
`V/F = {sa['V_F_L']} L`, implying a terminal half-life of
~{sa['derived_terminal_thalf_h']} h. A 400 mg twice-daily schedule gives a predicted
steady-state average concentration `Cavg,ss = Dose/(CL/F . tau) =
{sa['cavg_ss_mg_per_L_400mg_BID']} mg/L`, within the clinically reported ~3 to 5
mg/L range (in range: {sa['in_clinical_range']}). This is a forward consistency
check, not a fit: the published human PK parameters produce a plausible clinical
exposure, anchoring the model's PK scale at human dose.

## How this connects to the simulation

- The fitted plasma and tumor curves are exported in the JSON as a normalized
  `frompk_series` (time in hours, concentrations normalized to plasma Cmax), which
  is directly consumable as a measured-anchored `DoseSchedule::FromPk` series (the
  #239 plasma-to-spatial bridge). For the first time the spatial sims can be driven
  by a drug exposure profile tied to measured data, not a placeholder shape.
- The measured `Kp ~ 0.90` and the ~2 h plasma to tumor delay are reference anchors
  for the `tumor_pk` presets' `partition_coeff` and transfer rates.

## Caveats (what this is and is NOT)

1. **Two drugs, not every preset.** This anchors the PK structure + partition +
   delay to IKE (mouse) and sorafenib (human). It does NOT recalibrate the per-tumor
   presets (breast/pancreatic/GBM/melanoma/sarcoma): no public per-tumor-type
   measured PK exists for a ferroptosis inducer, so those stay documented estimates.
2. **RSL3/ML tool-compound PK is genuinely absent.** The biochem layer's RSL3 /
   ML162 / ML210 have no usable in-vivo PK (erastin's poor stability is why IKE was
   engineered). IKE is the in-vivo system-xc- proxy; the tool-compound PK gap is
   flagged, not fabricated.
3. **Mouse + summary NCA.** IKE is a mouse anchor and the targets are published
   summary NCA metrics (the per-timepoint curve is in the cited supplement), not
   digitized points. Sorafenib is human but a forward check, not a per-timepoint fit.
4. **Not whole-body PBPK.** The coupling is plasma + a single tumor tissue
   compartment with a measured partition, not a multi-organ physiological model.
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=DATA_CSV)
    args = ap.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
