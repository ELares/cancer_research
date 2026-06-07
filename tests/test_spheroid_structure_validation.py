"""Tests for the #333 spheroid zone-geometry validation vs Browning 2021.

The fetch needs network and is not run in CI; these cover the pure binning /
valid-range logic, a drift guard tying the Python model thresholds to the Rust
`spheroid` defaults, and the committed validation result.
"""

import json
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import validate_spheroid_structure as vss  # noqa: E402

VAL_JSON = REPO_ROOT / "analysis" / "calibration" / "spheroid-structure-validation.json"
SPHEROID_RS = REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "spheroid.rs"


# --------------------------------------------------------------------------
# Pure logic
# --------------------------------------------------------------------------


def test_summarize_bins_and_core_fraction():
    rows = (
        [{"cell_line": "x", "R": 150.0, "phi": 0.1, "eta": 0.0}] * 4
        + [{"cell_line": "x", "R": 350.0, "phi": 0.9, "eta": 0.65}] * 6
    )
    s = vss.summarize(rows)
    small = next(b for b in s if b["r_lo_um"] == 0)
    big = next(b for b in s if b["r_lo_um"] == 300)
    assert small["n"] == 4 and small["eta_median"] == 0.0 and small["frac_with_necrotic_core"] == 0.0
    assert big["n"] == 6 and big["eta_median"] == pytest.approx(0.65) and big["frac_with_necrotic_core"] == 1.0


def test_valid_size_range_picks_matching_large_bin():
    summary = [
        {"r_lo_um": 0, "model_phi_abs_err": 0.76, "model_eta_abs_err": 0.73, "frac_with_necrotic_core": 0.0},
        {"r_lo_um": 200, "model_phi_abs_err": 0.31, "model_eta_abs_err": 0.73, "frac_with_necrotic_core": 0.45},
        {"r_lo_um": 300, "model_phi_abs_err": 0.03, "model_eta_abs_err": 0.10, "frac_with_necrotic_core": 1.0},
    ]
    assert vss.valid_size_range(summary) == 300


def test_valid_size_range_none_when_no_bin_matches():
    summary = [
        {"r_lo_um": 0, "model_phi_abs_err": 0.7, "model_eta_abs_err": 0.7, "frac_with_necrotic_core": 0.0},
    ]
    assert vss.valid_size_range(summary) is None


def test_model_thresholds_match_rust_spheroid_defaults():
    """The Python validation's model thresholds must stay tied to the Rust
    `spheroid` volume-fraction defaults (glycolytic_frac, oxphos_frac); the radius
    boundaries are their cube roots. Catches drift if the Rust defaults are retuned."""
    src = SPHEROID_RS.read_text()
    gly = float(re.search(r"glycolytic_frac:\s*([0-9.]+)", src).group(1))
    oxp = float(re.search(r"oxphos_frac:\s*([0-9.]+)", src).group(1))
    assert vss.MODEL_RIM_BOUNDARY == pytest.approx(gly ** (1 / 3), abs=1e-6)
    assert vss.MODEL_CORE_BOUNDARY == pytest.approx(oxp ** (1 / 3), abs=1e-6)


# --------------------------------------------------------------------------
# Committed validation result
# --------------------------------------------------------------------------


def test_committed_validation_structure_and_finding():
    r = json.loads(VAL_JSON.read_text())
    assert r["n_spheroids"] > 900
    assert r["valid_radius_um_min"] is not None
    assert r["valid_diameter_um_min"] == 2 * r["valid_radius_um_min"]
    bins = {b["r_lo_um"]: b for b in r["size_bins"]}
    # the load-bearing finding: small spheroids have ~no necrotic core, large ones do
    assert bins[0]["frac_with_necrotic_core"] < 0.1
    assert bins[300]["frac_with_necrotic_core"] > 0.9
    # the fixed model core boundary is a poor fit for small spheroids (eta ~ 0 vs 0.73)
    assert bins[0]["eta_median"] < 0.1
    assert bins[300]["eta_median"] > 0.4


def test_valid_range_is_large_spheroids_only():
    r = json.loads(VAL_JSON.read_text())
    # the fixed limiting-structure thresholds should only validate for larger spheroids
    assert r["valid_radius_um_min"] >= 250


# --------------------------------------------------------------------------
# Size-aware refinement (#333)
# --------------------------------------------------------------------------


def test_size_aware_boundaries_ramp_from_zero_to_fixed_limit():
    # Small spheroid: below both onsets -> no core, no thinned rim (both 0).
    phi_s, eta_s = vss.size_aware_boundaries(150.0)
    assert phi_s == 0.0 and eta_s == 0.0
    # Large spheroid: above both full radii -> reduces exactly to the fixed limits.
    phi_l, eta_l = vss.size_aware_boundaries(600.0)
    assert phi_l == pytest.approx(vss.MODEL_RIM_BOUNDARY, abs=1e-9)
    assert eta_l == pytest.approx(vss.MODEL_CORE_BOUNDARY, abs=1e-9)


def test_size_aware_reduces_small_spheroid_error_on_committed_bins():
    r = json.loads(VAL_JSON.read_text())
    sa = vss.evaluate_size_aware(r["size_bins"])
    # Overall the size-aware model fits the bin medians better than the fixed one.
    assert sa["improves"] is True
    assert sa["size_aware_mean_abs_err"] < sa["fixed_mean_abs_err"]
    per_bin = {x["r_lo_um"]: x for x in sa["per_bin"]}
    # The gain is concentrated in the small bins (the core over-prediction the
    # fixed model makes); the large bin is essentially unchanged (reduces to fixed).
    small = per_bin[0]
    assert small["size_aware_eta_abs_err"] < small["fixed_eta_abs_err"]
    assert small["size_aware_phi_abs_err"] < small["fixed_phi_abs_err"]
    big = per_bin[300]
    assert big["size_aware_eta_abs_err"] == pytest.approx(big["fixed_eta_abs_err"], abs=0.12)


def test_committed_result_has_size_aware_refinement():
    r = json.loads(VAL_JSON.read_text())
    sa = r["size_aware_refinement"]
    assert sa["improves"] is True
    assert sa["size_aware_mean_abs_err"] < sa["fixed_mean_abs_err"]


def test_size_aware_params_match_rust_constants():
    """Drift guard: the Python SIZE_AWARE ramp params must equal the Rust
    spheroid::SizeAwareZones::literature() values."""
    rust = vss.size_aware_rust_constants()
    assert rust == vss.SIZE_AWARE
