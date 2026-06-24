#!/usr/bin/env python3
"""Validate the model's single source-independent exo-ROS kill threshold against
measured in-vivo singlet-oxygen necrosis thresholds (#464).

The simulation feeds SDT and PDT exogenous ROS into the SAME bistable lipid-
peroxidation death switch, with one `death_threshold` and `sdt_ros == pdt_ros`
(both 5.0). That encodes a strong modeling assumption: the dose-to-kill is a
property of the cell/tissue, NOT of the particular ROS source. Zhu et al. 2015
(PMID 25927018; PMC4410434, "In-vivo singlet oxygen threshold doses for PDT")
measured the reacted-singlet-oxygen necrosis threshold [1O2]rx in vivo for three
clinically distinct photosensitizers and found it approximately PHOTOSENSITIZER-
INDEPENDENT:

  Photofrin  ~0.56 mM
  BPD        ~0.72 mM
  mTHPC      ~0.40 mM

i.e. one reacted-ROS dose-to-kill (~0.5 mM, within a factor of 1.8 across three
chemically different sensitizers). That is real-world support for the model's
design choice; this script makes that validation explicit and drift-guards the
two Rust constants it depends on (sdt_ros == pdt_ros).

SCOPE / HONESTY. The model is dimensionless, so this validates the FORM (a single
source-independent kill threshold) and the ORDER OF MAGNITUDE, not a unit-rigorous
calibration of `death_threshold`. The separate question of the O2-DEPENDENCE
functional form is only partly addressed: the model's `oxygen::o2_dependent_exo_factor`
(#336) is a LINEAR `(1 - dep) + dep*o2_supply`, whereas Zhu's macroscopic
singlet-oxygen model makes the singlet-oxygen quantum yield a saturating
(Michaelis-type) function of [O2]; the linear form is a first-order approximation,
and the precise oxygen-quenching constant lives in the full text (not the abstract)
so it is flagged here rather than fabricated.

Pure stdlib; runs in CI. Committed result: analysis/calibration/pdt-threshold-validation.{md,json}.
"""
import json
import re
import statistics as st
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PARAMS_RS = REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "params.rs"
OUT_DIR = REPO_ROOT / "analysis" / "calibration"

# Measured in-vivo reacted-singlet-oxygen necrosis thresholds [1O2]rx (mM),
# Zhu et al. 2015 (PMID 25927018 / PMC4410434).
MEASURED_THRESHOLDS_MM = {
    "Photofrin": 0.56,
    "BPD": 0.72,
    "mTHPC": 0.40,
}


def _rust_default_value(src: str, field: str) -> float:
    """Parse a numeric field from the params.rs Default impl (drift guard)."""
    m = re.search(re.escape(field) + r"\s*:\s*([0-9.]+)", src)
    if not m:
        raise SystemExit(f"could not parse {field} from {PARAMS_RS.name}")
    return float(m.group(1))


def validate() -> dict:
    src = PARAMS_RS.read_text(encoding="utf-8")
    sdt_ros = _rust_default_value(src, "sdt_ros")
    pdt_ros = _rust_default_value(src, "pdt_ros")
    death_threshold = _rust_default_value(src, "death_threshold")

    vals = list(MEASURED_THRESHOLDS_MM.values())
    lo, hi = min(vals), max(vals)
    ratio = hi / lo
    cv = st.pstdev(vals) / st.mean(vals)

    # The model encodes a single source-independent kill threshold IFF the two
    # exogenous-ROS sources share a peak (sdt_ros == pdt_ros) feeding one
    # death_threshold. Zhu 2015 supports source-independence (ratio < ~2).
    model_source_independent = sdt_ros == pdt_ros
    measured_source_independent = ratio < 2.0

    return {
        "source": "Zhu et al. 2015, In-vivo singlet oxygen threshold doses for PDT (PMID 25927018, PMC4410434)",
        "measured_threshold_mm": MEASURED_THRESHOLDS_MM,
        "measured_threshold_mean_mm": round(st.mean(vals), 3),
        "measured_max_over_min_ratio": round(ratio, 3),
        "measured_coefficient_of_variation": round(cv, 3),
        "measured_is_photosensitizer_independent": measured_source_independent,
        "model_sdt_ros": sdt_ros,
        "model_pdt_ros": pdt_ros,
        "model_death_threshold": death_threshold,
        "model_uses_single_source_independent_threshold": model_source_independent,
        "validation": (
            "PASS: the model uses one death_threshold with sdt_ros == pdt_ros (a single "
            "source-independent kill threshold), and the measured reacted-singlet-oxygen "
            "necrosis threshold is photosensitizer-independent within a factor of "
            f"{round(ratio, 2)} (~0.5 mM), so the design choice is supported."
            if (model_source_independent and measured_source_independent)
            else "FAIL: model/measured source-independence disagree."
        ),
        "o2_dependence_note": (
            "The model's oxygen::o2_dependent_exo_factor (#336) is LINEAR in O2; "
            "Zhu's macroscopic singlet-oxygen model makes the quantum yield a saturating "
            "(Michaelis-type) function of [O2]. The linear form is a first-order "
            "approximation; the precise oxygen-quenching constant is in the full text "
            "(not the abstract) and is flagged here, not fabricated."
        ),
        "scope": (
            "Validates the FORM (single source-independent kill threshold) + order of "
            "magnitude, NOT a unit-rigorous calibration of the dimensionless death_threshold."
        ),
    }


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="print, write nothing")
    args = ap.parse_args()

    result = validate()
    print(result["validation"])
    print(
        f"  measured: {result['measured_threshold_mm']} mM "
        f"(mean {result['measured_threshold_mean_mm']}, max/min {result['measured_max_over_min_ratio']}x)"
    )
    print(
        f"  model: sdt_ros={result['model_sdt_ros']} == pdt_ros={result['model_pdt_ros']}, "
        f"death_threshold={result['model_death_threshold']}"
    )
    if args.check:
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "pdt-threshold-validation.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    print(f"wrote {OUT_DIR / 'pdt-threshold-validation.json'}")


if __name__ == "__main__":
    main()
