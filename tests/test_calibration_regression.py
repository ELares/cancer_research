"""Calibration-regression gate (#499).

The `sim-tme-3d` byte-identity SHA gate proves a software invariant (the default
24-condition matrix is unchanged), not biology. This test guards the thing that
actually matters scientifically: that the legs anchored to independent published
data still reproduce their committed held-out metrics. If a change silently moves
a calibrated leg off its anchor, this fails.

Pure-Python stdlib (no compiled extension), so it runs in CI. It (1) asserts the
committed calibration result artifacts still hold their anchored metrics within
tolerance, including the legs whose full fit needs the compiled extension (so the
committed numbers cannot be edited away unnoticed), and (2) live-re-runs the
trigger-wave validator (which is pure-Python) as a computed cross-check.

Anchors guarded:
- CTRPv2 GPX4-inhibitor kill-switch fit (#330): ML162 fit RMSE, ML210 held-out RMSE.
- CTRPv2 System Xc-/erastin fit (#502): erastin fit RMSE + mechanism specificity.
- Tumor-PK partition vs IKE (#334): tissue:plasma Kp.
- Ferroptotic trigger-wave speed vs Co 2024 (#482): baseline 5.52 um/min.
- Krogh drug-penetration lengths vs measured (#335).
- Spheroid zone geometry vs Browning 2021 (#333): size-aware boundary error.
- PDT/SDT source-independent kill threshold vs Zhu 2015 (#464).
"""

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CALIB = REPO_ROOT / "analysis" / "calibration"
sys.path.insert(0, str(REPO_ROOT / "scripts"))


def _load(name: str) -> dict:
    return json.loads((CALIB / name).read_text(encoding="utf-8"))


def test_kill_switch_fit_and_heldout_rmse_hold(_=None):
    """#330: the in-vitro RSL3 kill switch fit (ML162) and held-out validation
    (ML210) must stay anchored, and the default in-vivo switch must stay clearly
    worse (the 'calibration matters' invariant)."""
    d = _load("kill-switch-calibration.json")
    assert d["fit_compound"] == "ML162" and d["heldout_compound"] == "ML210"
    assert d["fit_rmse"] <= 0.06, f"ML162 fit RMSE drifted: {d['fit_rmse']}"
    assert d["heldout_rmse"] <= 0.08, f"ML210 held-out RMSE drifted: {d['heldout_rmse']}"
    # The uncalibrated in-vivo default must remain much worse than the fit, or
    # the calibration is no longer doing anything.
    assert d["default_uncalibrated_rmse"] >= 5.0 * d["fit_rmse"], (
        f"calibration no longer separates from the default: "
        f"fit={d['fit_rmse']} default={d['default_uncalibrated_rmse']}"
    )


def test_erastin_system_xc_fit_holds():
    """#502: the core's second data-anchored inducer mechanism (System Xc-/erastin)
    must stay fit to the CTRPv2 erastin curve, and erastin must still RAISE death
    monotonically (a broken mechanism would zero the increment). The erastin curve
    is flat-then-steep so the fit is a PARTIAL anchor (RMSE ~0.10, looser than the
    GPX4i leg), reflected in the tolerance."""
    d = _load("erastin-calibration.json")
    assert d["compound"] == "ERASTIN"
    assert d["fit_rmse"] <= 0.12, f"erastin fit RMSE drifted: {d['fit_rmse']}"
    # Erastin must drive a substantial dose-dependent kill via System Xc- (the
    # mechanism works), well above the cascade's Control baseline.
    assert d["erastin_increment_top_dose"] >= 0.3, (
        f"erastin top-dose kill increment collapsed: {d['erastin_increment_top_dose']}"
    )
    for k in ("lp_propagation", "lp_rate", "k_erastin", "hill"):
        assert k in d["calibrated_params"], f"missing calibrated param {k}"


def test_tumor_pk_partition_holds():
    """#334: the measured IKE tissue:plasma partition Kp ~= 0.90 anchor."""
    d = _load("pk-calibration.json")
    kp = d["measured"]["partition_kp_auc_ratio"]
    assert 0.85 <= kp <= 0.95, f"tumor-PK partition Kp drifted: {kp}"


def test_trigger_wave_baseline_reproduced_live():
    """#482: re-run the (pure-Python) trigger-wave validator and confirm the
    model baseline still lands on the measured 5.52 um/min and all checks pass."""
    import validate_trigger_wave as vt

    result = vt.validate()
    assert result["all_passed"], result["checks"]
    base = result["model_analytical"]["baseline"]
    assert abs(base - 5.52) < 0.6, f"trigger-wave baseline drifted: {base}"
    # The committed artifact must agree with the live run.
    committed = _load("trigger-wave-validation.json")["model_analytical"]["baseline"]
    assert abs(committed - base) < 1e-6, "committed trigger-wave result is stale"


def test_penetration_lengths_hold():
    """#335: the doxorubicin transport-reference length stays within the
    measured band and the drift-guarded lambda presets are unchanged."""
    d = _load("penetration-validation.json")
    assert d["model"]["Doxorubicin-transport"]["lambda_um"] == 50.0
    assert d["model"]["RSL3-like"]["lambda_um"] == 100.0
    assert d["drift_guard"]["Doxorubicin-transport"]["rust_lambda_um"] == 50.0
    assert d["doxorubicin_half_distance_within_a_measured_range"] is True


def test_spheroid_zone_geometry_holds():
    """#333: the size-aware zone thresholds keep their large improvement over the
    fixed thresholds on the Browning 2021 bins."""
    d = _load("spheroid-structure-validation.json")
    sar = d["size_aware_refinement"]
    assert sar["size_aware_mean_abs_err"] < 0.05
    assert sar["size_aware_mean_abs_err"] < sar["fixed_mean_abs_err"]


def test_pdt_threshold_design_supported():
    """#464: the single source-independent kill-threshold design stays supported
    by the measured photosensitizer-independent singlet-O2 necrosis threshold."""
    d = _load("pdt-threshold-validation.json")
    assert d["model_uses_single_source_independent_threshold"] is True
    assert d["measured_max_over_min_ratio"] <= 2.0
    assert d["validation"].startswith("PASS")


def test_all_anchored_artifacts_present():
    """Guards that no anchored calibration artifact is deleted (the gate would
    otherwise silently pass on a missing leg)."""
    for name in (
        "kill-switch-calibration.json",
        "erastin-calibration.json",
        "pk-calibration.json",
        "trigger-wave-validation.json",
        "penetration-validation.json",
        "spheroid-structure-validation.json",
        "pdt-threshold-validation.json",
    ):
        assert (CALIB / name).exists(), f"missing anchored artifact: {name}"
