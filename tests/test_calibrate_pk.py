"""Tests for the tumor-PK measured-data anchoring (#334).

The PK fit is pure Python + scipy (no compiled `ferroptosis_core` extension), so
unlike the kill-switch calibration these run in CI. They cover the analytical PK
model, the closed-form multi-compartment floor finding, the exact plasma anchoring,
the out-of-fit tumor delay prediction, the data file, and the end-to-end run.
"""

import json
import math
import sys
from types import SimpleNamespace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import calibrate_pk as cp  # noqa: E402


# --- analytical PK model ---

def test_plasma_starts_zero_peaks_at_tmax():
    ka, ke, scale = 1.3, 0.38, 12000.0
    assert cp.plasma_conc(0.0, ka, ke, scale) == 0.0
    tm = cp.plasma_tmax(ka, ke)
    peak = cp.plasma_conc(tm, ka, ke, scale)
    # the analytical Tmax is the true argmax of the curve
    for dt in (-0.2, -0.05, 0.05, 0.2):
        assert cp.plasma_conc(tm + dt, ka, ke, scale) <= peak + 1e-9


def test_plasma_auc_matches_numeric_integral():
    ka, ke, scale = 1.3, 0.38, 12000.0
    analytic = cp.plasma_auc(ka, ke, scale)
    # trapezoid over a long horizon (>> 1/ke) approximates AUC(0..inf)
    n, t_hi = 200000, 200.0
    h = t_hi / n
    s = 0.5 * (cp.plasma_conc(0.0, ka, ke, scale) + cp.plasma_conc(t_hi, ka, ke, scale))
    for i in range(1, n):
        s += cp.plasma_conc(i * h, ka, ke, scale)
    numeric = s * h
    assert abs(numeric - analytic) / analytic < 1e-3


def test_ka_from_tmax_roundtrip():
    ke = 0.38
    for tmax in (0.8, 1.35, 2.0):
        ka = cp.ka_from_tmax(tmax, ke)
        assert ka > ke
        assert abs(cp.plasma_tmax(ka, ke) - tmax) < 1e-6


def test_ka_from_tmax_raises_on_unachievable_tmax():
    # a 1-cmt absorption model cannot peak later than 1/ke; asking for it must raise
    # rather than silently return a wrong ka (ke=0.38 -> 1/ke=2.63 h)
    import pytest
    ke = 0.38
    with pytest.raises(ValueError):
        cp.ka_from_tmax(1.0 / ke + 0.5, ke)


def test_tumor_analytical_matches_rk4():
    ka, ke, scale, k_pt, k_te = 1.28, 0.38, 12000.0, 0.179, 0.198
    for t in (1.0, 3.0, 6.0, 10.0):
        a = cp.tumor_conc(t, ka, ke, scale, k_pt, k_te)
        n = cp.tumor_conc_numeric(t, ka, ke, scale, k_pt, k_te, dt=0.0005)
        assert abs(a - n) / max(abs(a), 1.0) < 2e-3


def test_tumor_auc_is_mass_balance():
    ka, ke, scale, k_pt, k_te = 1.28, 0.38, 12000.0, 0.179, 0.198
    expected = (k_pt / k_te) * cp.plasma_auc(ka, ke, scale)
    assert abs(cp.tumor_auc(scale, ka, ke, k_pt, k_te) - expected) < 1e-6


# --- the closed-form multi-compartment finding ---

def test_one_cmt_floor_is_e_times_tmax():
    for tmax in (0.5, 1.35, 2.0):
        assert abs(cp.one_cmt_auc_over_cmax_floor(tmax) - math.e * tmax) < 1e-12


def test_ike_plasma_is_below_one_cmt_floor():
    # measured IKE plasma: AUC/Cmax must lie below e*Tmax (proof of a fast phase)
    auc_over_cmax = 10926.0 / 5185.0
    floor = cp.one_cmt_auc_over_cmax_floor(1.35)
    assert auc_over_cmax < floor
    assert auc_over_cmax == 10926.0 / 5185.0  # exact arithmetic, no rounding drift


# --- exact plasma anchoring + tumor prediction ---

def test_fit_plasma_exact_on_shape_metrics_auc_overpredicts():
    params, pred, rel = cp.fit_plasma(1.35, 5185.0, 1.83, 10926.0)
    # Tmax, Cmax, terminal half-life are exact by construction
    assert abs(rel["Tmax"]) < 1e-3
    assert abs(rel["Cmax"]) < 1e-3
    assert abs(rel["thalf"]) < 1e-3
    # AUC is a prediction and a 1-cmt model necessarily over-predicts it here
    assert rel["AUC"] > 0.5
    assert params["ka"] > params["ke"]


def test_derive_tumor_captures_delay_direction():
    params, _, _ = cp.fit_plasma(1.35, 5185.0, 1.83, 10926.0)
    tumor = cp.derive_tumor(params, kp_measured=0.9022, thalf_tumor=3.50)
    # tumor peaks LATER than plasma (delay direction) and clears slower
    assert tumor["pred"]["Tmax"] > 1.35
    assert abs(cp.LN2 / tumor["k_te"] - 3.50) < 0.1  # terminal half-life by construction
    # k_pt = Kp * k_te
    assert abs(tumor["k_pt"] - 0.9022 * tumor["k_te"]) < 1e-9


# --- the committed data file ---

def test_load_targets_has_ike_pair_and_sorafenib():
    data = cp.load_targets()
    ike = data["IKE"][cp.IKE_POP]["IP"]
    assert ike["plasma"]["Cmax"] == 5185.0 and ike["plasma"]["Tmax"] == 1.35
    assert ike["tumor"]["Cmax"] == 2516.0 and ike["tumor"]["thalf"] == 3.50
    sor = data["sorafenib"]["solid-tumor popPK"]["PO"]["central"]
    assert sor["CL_F"] == 8.13 and sor["V_F"] == 213.0


# --- end-to-end run (CI-safe; deterministic) ---

def test_run_writes_outputs_and_headline_anchors(tmp_path):
    result = cp.run(SimpleNamespace(data=cp.DATA_CSV))
    assert cp.OUT_MD.exists() and cp.OUT_JSON.exists()
    # headline measured partition ~ 0.90
    assert abs(result["measured"]["partition_kp_auc_ratio"] - 0.9022) < 1e-3
    # the multi-compartment finding fires
    mc = result["multicompartment_finding"]
    assert mc["below_one_cmt_floor"] is True
    assert mc["plasma_auc_over_cmax_h"] < mc["one_cmt_floor_e_tmax_h"]
    # sorafenib forward check lands in the clinical window
    assert result["sorafenib_anchor"]["in_clinical_range"] is True
    # the FromPk export has the documented shape
    fp = result["frompk_series"]
    assert len(fp["time_h"]) == len(fp["plasma_norm"]) == len(fp["tumor_norm"])
    assert fp["plasma_norm"][0] == 0.0  # starts at zero


def test_run_is_deterministic():
    a = cp.run(SimpleNamespace(data=cp.DATA_CSV))
    b = cp.run(SimpleNamespace(data=cp.DATA_CSV))
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
