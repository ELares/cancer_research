"""Tests for the Krogh drug-penetration validation (#335).

Pure Python stdlib (no compiled extension), so these run in CI. They cover the
lambda formula, the half-distance conversion, the drift-guard against the Rust
drug_transport.rs presets, the doxorubicin comparison outcome, and the honest
scope flags (RSL3 unvalidated, binding-site barrier not modeled).
"""

import json
import math
import sys
from types import SimpleNamespace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import validate_penetration as vp  # noqa: E402


def test_lambda_formula_matches_rust_values():
    # RSL3-like: D=5e-7, k=0.005 -> lambda = sqrt(50/0.005) = 100 um
    assert abs(vp.penetration_length_um(5.0e-7, 0.004, 0.001) - 100.0) < 1e-9
    # Doxorubicin: D=3e-7, k=0.012 -> lambda = sqrt(30/0.012) = 50 um
    assert abs(vp.penetration_length_um(3.0e-7, 0.01, 0.002) - 50.0) < 1e-9


def test_half_distance_is_lambda_ln2():
    assert abs(vp.half_distance_um(100.0) - 100.0 * math.log(2)) < 1e-9
    assert abs(vp.half_distance_um(50.0) - 34.657) < 1e-2


def test_model_penetration_values():
    m = vp.model_penetration()
    assert abs(m["RSL3-like"]["lambda_um"] - 100.0) < 1e-6
    assert abs(m["Doxorubicin-transport"]["lambda_um"] - 50.0) < 1e-6
    assert abs(m["Doxorubicin-transport"]["half_distance_um"] - 34.657) < 1e-2


def test_drift_guard_matches_rust_source():
    # The Rust presets must equal the Python-encoded MODEL_DRUGS (recomputed lambda agrees)
    report = vp.drift_guard()
    assert abs(report["RSL3-like"]["rust_lambda_um"] - 100.0) < 1e-6
    assert abs(report["Doxorubicin-transport"]["rust_lambda_um"] - 50.0) < 1e-6
    assert all(v["matches_python"] for v in report.values())


def test_measured_data_has_doxorubicin_anchors():
    rows = vp.load_measured()
    dox = vp.doxorubicin_targets(rows)
    sources = {t["source"].split()[0] for t in dox}
    assert "Primeau" in sources and "Tannock" in sources
    # Tannock 2002 range brackets the model half-distance (34.7 um); Primeau (40-50) does not
    tannock = next(t for t in dox if t["source"].startswith("Tannock"))
    assert tannock["low_um"] <= 34.7 <= tannock["high_um"]
    primeau = next(t for t in dox if t["source"].startswith("Primeau"))
    assert 34.7 < primeau["low_um"]


def test_run_outcome_and_honest_scope_flags():
    result = vp.run(SimpleNamespace(data=vp.DATA_CSV))
    assert vp.OUT_MD.exists() and vp.OUT_JSON.exists()
    # exponential form validated, model within at least one measured range
    assert result["exponential_form_supported"] is True
    assert result["doxorubicin_half_distance_within_a_measured_range"] is True
    # honest scope: RSL3 unvalidated, binding-site barrier not a model feature
    assert result["rsl3_penetration_validated"] is False
    assert result["binding_site_barrier_modeled"] is False
    # the Tannock comparison is within-range; Primeau is a small negative (model shorter) gap
    primeau = next(c for c in result["doxorubicin_comparisons"] if c["source"].startswith("Primeau"))
    assert primeau["within_measured_range"] is False
    assert -0.5 < primeau["rel_gap_to_nearest_edge"] < 0.0


def test_run_is_deterministic():
    a = vp.run(SimpleNamespace(data=vp.DATA_CSV))
    b = vp.run(SimpleNamespace(data=vp.DATA_CSV))
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
