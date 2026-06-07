#!/usr/bin/env python3
"""Validate the Krogh drug-penetration model against measured penetration data (#335).

The `drug_transport` module models drug concentration vs distance from the nearest
vessel as the exponential Krogh-cylinder steady state `C(r) = C0 * exp(-r/lambda)`,
with `lambda = sqrt(D / k_total)` (k_total = uptake_rate + metabolism_rate). The
penetration length per drug has been an estimate (CALIBRATION_STATUS.md: "Tissue drug
penetration, Uncalibrated"); the only external check was a self-consistency unit test
asserting lambda is in a literature range. This script anchors that check to PUBLISHED
MEASURED penetration metrics.

WHAT IS VALIDATED
-----------------
1. The exponential FORM. Primeau 2005 and Tannock 2002 both report drug concentration
   declining exponentially with distance from the vessel, so the model's
   C(r) = C0*exp(-r/lambda) functional form is the right one for this in-vivo geometry.
2. The penetration LENGTH for the doxorubicin transport reference. The model computes
   lambda = 50 um (half-distance lambda*ln2 = 34.7 um). We compare that half-distance to
   the measured doxorubicin half-distances: Tannock 2002 (25 to 75 um, time-resolved),
   Primeau 2005 (40 to 50 um), and the Minchinton & Tannock 2006 review range (40 to
   80 um). The model is at the conservative (shorter) end but within or adjacent to the
   measured spread.

WHAT IS NOT VALIDATED (honest scope)
------------------------------------
- The RSL3-like lambda (100 um). No ferroptosis inducer has published spatial
  penetration data (erastin/RSL3/ML/IKE), so this is an unvalidated extrapolation from
  the small-molecule transport class; flagged, not fabricated.
- The "binding-site barrier" half of #335. The model has NO binding-site-barrier
  mechanism: uptake is LINEAR (a constant `uptake_rate` folded into lambda), so
  penetration depth is dose-INDEPENDENT. A true binding-site barrier (Tannock) is
  SATURABLE: perivascular target binding sequesters drug until it saturates at higher
  dose, making penetration depth dose-dependent. A data-availability review concluded
  this mechanism is NOT warranted for the small molecules this model targets, and the
  mechanism is therefore deliberately NOT added (rather than added as unvalidatable
  complexity): (1) the strong dose-dependent binding-site barrier is fundamentally an
  ANTIBODY phenomenon (Fujimori/Weinstein 1990 PMID 2362198; Saga 1995 PMID 7568060
  dose-titration; Thurber/Wittrup 2008 PMID 18541331 derive penetration ~ sqrt(dose)
  for antibodies); (2) for small molecules the effect is physically WEAK, because a
  deep, slow-saturating binding sink barely saturates at achievable doses (the
  El-Kareh/Secomb tumor-cord model predicts little dose-through deepening for
  doxorubicin); and (3) there is NO extractable small-molecule dose-resolved
  penetration-depth dataset to validate such a term against. So #335's binding-site
  clause is resolved by evidence (not warranted + not validatable for this drug class),
  not by building a mechanism the data cannot support.
- Spheroid-specific data. The cleanest exponential penetration data is in-vivo tumor
  (drug vs distance-from-vessel), the same Krogh geometry the model assumes. Small-
  molecule drugs often fully penetrate spheroids (a poor exponential test), so the
  honest spheroid contrast is the slow/surface-limited drugs (platinum, paclitaxel),
  reported qualitatively.

A drift-guard re-reads the Rust drug_transport.rs presets so the validated lambda values
cannot silently diverge from the model they claim to validate.

Run (pure Python stdlib; runs in CI):
  python3 scripts/validate_penetration.py
Writes analysis/calibration/penetration-validation.md + .json.
"""

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_CSV = REPO_ROOT / "analysis" / "calibration" / "penetration_measured_data.csv"
RUST_SRC = REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "drug_transport.rs"
OUT_MD = REPO_ROOT / "analysis" / "calibration" / "penetration-validation.md"
OUT_JSON = REPO_ROOT / "analysis" / "calibration" / "penetration-validation.json"

LN2 = math.log(2.0)

# Model drug-transport presets, mirroring drug_transport.rs (verified by the drift-guard
# against the Rust source so these cannot silently diverge from the model).
MODEL_DRUGS = {
    "RSL3-like": {"d_cm2_s": 5.0e-7, "uptake": 0.004, "metabolism": 0.001},
    "Doxorubicin-transport": {"d_cm2_s": 3.0e-7, "uptake": 0.01, "metabolism": 0.002},
}


def penetration_length_um(d_cm2_s, uptake, metabolism):
    """Replicates ferroptosis-core drug_transport::penetration_length_um exactly:
    lambda = sqrt(D[um^2/s] / k_total), D converted cm^2/s -> um^2/s by 1e8."""
    k_total = uptake + metabolism
    if k_total <= 0:
        return float("inf")
    return math.sqrt((d_cm2_s * 1e8) / k_total)


def half_distance_um(lambda_um):
    """Distance at which C = 0.5*C0 for an exponential profile: lambda*ln2."""
    return lambda_um * LN2


def model_penetration():
    """Model lambda + half-distance for each shipped drug-transport preset."""
    out = {}
    for name, p in MODEL_DRUGS.items():
        lam = penetration_length_um(p["d_cm2_s"], p["uptake"], p["metabolism"])
        out[name] = {"lambda_um": lam, "half_distance_um": half_distance_um(lam)}
    return out


# ---------------------------------------------------------------------------
# Drift-guard: the Rust source is the source of truth for the preset constants
# ---------------------------------------------------------------------------

def _parse_rust_preset(src, fn_name):
    """Extract diffusion_coeff_cm2_s / uptake_rate / metabolism_rate from a Rust
    preset fn body (e.g. rsl3_like, doxorubicin_transport_reference)."""
    m = re.search(r"fn\s+" + re.escape(fn_name) + r"\s*\(\s*\)\s*->\s*DrugParams\s*\{(.*?)\n\}", src, re.DOTALL)
    if not m:
        return None
    body = m.group(1)

    def field(name):
        fm = re.search(re.escape(name) + r"\s*:\s*([0-9.eE+-]+)", body)
        return float(fm.group(1)) if fm else None

    return {
        "d_cm2_s": field("diffusion_coeff_cm2_s"),
        "uptake": field("uptake_rate"),
        "metabolism": field("metabolism_rate"),
    }


def drift_guard():
    """Confirm the Python-encoded MODEL_DRUGS match the Rust presets (recomputed
    lambda must agree). Returns a per-drug dict; raises on mismatch."""
    src = RUST_SRC.read_text(encoding="utf-8")
    mapping = {"RSL3-like": "rsl3_like", "Doxorubicin-transport": "doxorubicin_transport_reference"}
    report = {}
    for py_name, rust_fn in mapping.items():
        rust = _parse_rust_preset(src, rust_fn)
        if rust is None or any(v is None for v in rust.values()):
            raise ValueError(f"could not parse Rust preset {rust_fn} from {RUST_SRC.name}")
        py = MODEL_DRUGS[py_name]
        for k in ("d_cm2_s", "uptake", "metabolism"):
            if abs(rust[k] - py[k]) > 1e-12:
                raise ValueError(
                    f"{py_name}.{k} drift: Python {py[k]} vs Rust {rust[k]} in {RUST_SRC.name}. "
                    f"Update MODEL_DRUGS in validate_penetration.py to match the Rust preset."
                )
        lam_rust = penetration_length_um(rust["d_cm2_s"], rust["uptake"], rust["metabolism"])
        report[py_name] = {"rust_lambda_um": lam_rust, "matches_python": True}
    return report


# ---------------------------------------------------------------------------
# Measured data
# ---------------------------------------------------------------------------

def load_measured(path=DATA_CSV):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def doxorubicin_targets(rows):
    """Quantitative doxorubicin half-distance / penetration ranges (um) with source."""
    out = []
    for r in rows:
        if r["drug"] != "doxorubicin" or not r["value_low_um"]:
            continue
        out.append({
            "system": r["system"],
            "measurement": r["measurement"],
            "low_um": float(r["value_low_um"]),
            "high_um": float(r["value_high_um"]),
            "exponential_form": r["exponential_form"] == "yes",
            "source": r["source"],
            "pmid": r["source_pmid"],
        })
    return out


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run(args):
    drift = drift_guard()
    model = model_penetration()
    rows = load_measured(args.data)
    dox_targets = doxorubicin_targets(rows)

    dox_half = model["Doxorubicin-transport"]["half_distance_um"]

    # Compare the model doxorubicin half-distance to each measured range.
    comparisons = []
    for t in dox_targets:
        within = t["low_um"] <= dox_half <= t["high_um"]
        # signed gap vs the nearest edge (0 if within)
        if within:
            gap = 0.0
        elif dox_half < t["low_um"]:
            gap = dox_half / t["low_um"] - 1.0  # negative => model shorter
        else:
            gap = dox_half / t["high_um"] - 1.0  # positive => model deeper
        comparisons.append({
            "source": t["source"], "pmid": t["pmid"], "system": t["system"],
            "measured_range_um": [t["low_um"], t["high_um"]],
            "model_half_distance_um": round(dox_half, 1),
            "within_measured_range": within,
            "rel_gap_to_nearest_edge": round(gap, 3),
            "reports_exponential_form": t["exponential_form"],
        })

    exp_supported = any(c["reports_exponential_form"] for c in comparisons)
    within_any = any(c["within_measured_range"] for c in comparisons)

    qualitative = [
        {"drug": r["drug"], "system": r["system"], "measurement": r["measurement"],
         "source": r["source"], "pmid": r["source_pmid"], "note": r["note"]}
        for r in rows if r["drug"] in ("cisplatin", "paclitaxel")
    ]

    result = {
        "model": {k: {"lambda_um": round(v["lambda_um"], 1),
                      "half_distance_um": round(v["half_distance_um"], 1)} for k, v in model.items()},
        "drift_guard": {k: {"rust_lambda_um": round(v["rust_lambda_um"], 1)} for k, v in drift.items()},
        "exponential_form_supported": exp_supported,
        "doxorubicin_half_distance_within_a_measured_range": within_any,
        "doxorubicin_comparisons": comparisons,
        "rsl3_penetration_validated": False,
        "rsl3_note": ("No ferroptosis inducer (erastin/RSL3/ML/IKE) has published spatial "
                      "penetration data; the RSL3-like lambda=100 um is an unvalidated "
                      "extrapolation from the small-molecule transport class."),
        "binding_site_barrier_modeled": False,
        "binding_site_barrier_note": ("The model has no binding-site barrier: uptake is LINEAR "
                                      "(folded into lambda) so penetration depth is dose-INDEPENDENT. "
                                      "A true binding-site barrier is saturable (dose-dependent depth). "
                                      "A data-availability review concluded this mechanism is NOT warranted "
                                      "for the small molecules this model targets and is deliberately NOT "
                                      "added: the strong dose-dependent barrier is an ANTIBODY phenomenon "
                                      "(Fujimori 1990, Saga 1995, Thurber/Wittrup 2008 penetration~sqrt(dose)), "
                                      "it is physically WEAK for small molecules (a deep binding sink barely "
                                      "saturates at achievable doses, El-Kareh/Secomb), and NO extractable "
                                      "small-molecule dose-resolved penetration-depth data exists to validate "
                                      "it. #335 binding-site clause resolved by evidence, not by building it."),
        "binding_site_barrier_decision": "not_added_not_warranted_for_small_molecules",
        "qualitative_contrasts": qualitative,
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_report(result)
    print(f"drift-guard OK (Rust == Python presets)")
    print(f"model doxorubicin half-distance = {dox_half:.1f} um (lambda {model['Doxorubicin-transport']['lambda_um']:.0f} um)")
    print(f"exponential form supported = {exp_supported}; within a measured range = {within_any}")
    for c in comparisons:
        flag = "within" if c["within_measured_range"] else f"gap {c['rel_gap_to_nearest_edge']:+.0%}"
        print(f"  vs {c['source']} {c['measured_range_um']} um: {flag}")
    print(f"RSL3 penetration validated = False (no ferroptosis-inducer data); binding-site barrier modeled = False")
    print(f"wrote {OUT_MD.relative_to(REPO_ROOT)} + {OUT_JSON.relative_to(REPO_ROOT)}")
    return result


def write_report(r):
    m = r["model"]

    def comp_row(c):
        flag = "within range" if c["within_measured_range"] else f"{c['rel_gap_to_nearest_edge']:+.0%} vs nearest edge"
        return (f"| {c['source']} (PMID {c['pmid']}) | {c['system']} | "
                f"{c['measured_range_um'][0]:.0f} to {c['measured_range_um'][1]:.0f} | "
                f"{c['model_half_distance_um']:.1f} | {flag} |")

    comp_rows = "\n".join(comp_row(c) for c in r["doxorubicin_comparisons"])
    qual_rows = "\n".join(
        f"| {q['drug']} | {q['system']} | {q['note']} |" for q in r["qualitative_contrasts"]
    )

    md = f"""# Krogh drug-penetration validation vs measured data (#335)

Generated by `scripts/validate_penetration.py` (pure Python stdlib; runs in CI).
Target data: `analysis/calibration/penetration_measured_data.csv`.

## What this validates

The `drug_transport` Krogh model `C(r) = C0 * exp(-r/lambda)` with
`lambda = sqrt(D / k_total)`. The model lambda is drug-dependent (not tissue-dependent):

| model drug | lambda (um) | half-distance (um) |
|---|---|---|
| RSL3-like | {m['RSL3-like']['lambda_um']:.0f} | {m['RSL3-like']['half_distance_um']:.1f} |
| Doxorubicin-transport | {m['Doxorubicin-transport']['lambda_um']:.0f} | {m['Doxorubicin-transport']['half_distance_um']:.1f} |

A drift-guard re-reads the Rust `drug_transport.rs` presets and confirms these lambda
values are the ones the model actually uses (so the validation cannot silently diverge
from the model).

## 1. Exponential form: supported

Primeau 2005 and Tannock 2002 both report drug concentration declining EXPONENTIALLY
with distance from the vessel, so the model's `C(r) = C0 * exp(-r/lambda)` functional
form is the right one for this in-vivo geometry. Exponential form supported:
**{r['exponential_form_supported']}**.

## 2. Doxorubicin penetration length: within / adjacent to measured

The doxorubicin transport reference is the model preset with published comparators.
Its half-distance ({m['Doxorubicin-transport']['half_distance_um']:.1f} um) compares to
the measured doxorubicin half-distances as:

| measured source | system | measured half-distance (um) | model (um) | result |
|---|---|---|---|---|
{comp_rows}

The model sits at the conservative (shorter-penetration) end: within Tannock 2002's
time-resolved 25 to 75 um spread, and modestly below the Primeau 2005 / Minchinton
review ranges (40 to 80 um). Within a measured range:
**{r['doxorubicin_half_distance_within_a_measured_range']}**. The model is therefore the
right order of magnitude and right functional form, slightly conservative on absolute
depth versus the deeper in-vivo datasets.

## What this does NOT validate (honest scope)

### RSL3 penetration length is unvalidated
{r['rsl3_note']} So the validation above transfers to the model's RSL3 penetration only
by analogy of molecular weight and transport class, not by direct measurement.

### The binding-site barrier: deliberately not added (resolved by evidence)
{r['binding_site_barrier_note']}

Issue #335 asks to validate the binding-site barrier. The honest, evidence-led
finding from a data-availability review is that this mechanism is NOT warranted for
the small molecules this model targets, so it is deliberately NOT added rather than
introduced as unvalidatable complexity:

1. The strong dose-dependent binding-site barrier is fundamentally an ANTIBODY
   phenomenon. It was discovered with antibodies (Fujimori and Weinstein 1990,
   PMID 2362198), the canonical dose-titration that shows penetration depth rising
   with dose is an antibody study (Saga 1995, PMID 7568060), and the closed-form
   signature (penetration distance scaling as the square root of dose) is derived
   for antibodies (Thurber and Wittrup 2008, PMID 18541331).
2. For small molecules the effect is physically WEAK: a deep, slowly-saturating
   binding sink barely saturates at achievable doses, so a saturable-uptake tumor-cord
   model predicts little dose-through deepening for doxorubicin (El-Kareh and Secomb).
3. There is NO extractable small-molecule dose-resolved penetration-depth dataset to
   validate such a term against (dose-dependence is asserted only qualitatively, e.g.
   gemcitabine in multicellular layers, with no numeric depth-vs-dose profiles).

So the binding-site clause of #335 is resolved by evidence (not warranted and not
validatable for this drug class), not by adding a mechanism the data cannot support.
Should this model ever be extended to antibodies or ADCs, a saturable Michaelis-Menten
uptake term validated against the antibody dose-titration data above would be the
appropriate addition.

### Spheroid-specific data and slow-penetrator contrasts
The cleanest exponential penetration data is in-vivo tumor (drug vs distance-from-vessel),
the same Krogh geometry the model assumes. Small-molecule drugs often fully penetrate
spheroids (a poor exponential test), so the honest spheroid contrast is the slow /
surface-limited drugs:

| drug | system | behavior (qualitative contrast) |
|---|---|---|
{qual_rows}

These are poorly-penetrating contrasts (the model's lambda is a passive-diffusion length,
so it does not by itself reproduce surface-limited or saturable behavior); they bound the
model's applicability rather than calibrate its lambda.

## Bottom line

The Krogh exponential penetration FORM is validated against measured in-vivo data, and
the doxorubicin penetration LENGTH is the right order of magnitude (conservative end of
the measured spread). The ferroptosis-inducer penetration length remains unvalidated for
lack of any ferroptosis-inducer penetration data, and the dose-dependent binding-site
barrier is deliberately NOT added because it is an antibody phenomenon, physically weak
for the small molecules this model targets, and unvalidatable for this drug class from
public data. The penetration-gradient comparisons in the manuscript stay directional,
not calibrated magnitudes.
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
