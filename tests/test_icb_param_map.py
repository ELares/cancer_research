"""Tests for the Talkington & Kearsley 2025 ICB parameter mapping (#472).

Pure-Python (no compiled extension, runs in CI). Guards (1) the verified
published constants against drift and (2) the conversion math the doc relies on.
"""

import math
from pathlib import Path

import pytest

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import icb_param_map as ipm  # noqa: E402


def test_published_constants_match_the_paper():
    """The verified Table 1 values (DOI 10.1002/cso2.70007) must not silently
    drift. These are the load-bearing numbers the mapping doc cites."""
    assert ipm.TALKINGTON_PARAMS["kill_rate"].value == pytest.approx(1.101e-7)
    assert ipm.TALKINGTON_PARAMS["exhaustion_threshold"].value == pytest.approx(1.0e4)
    assert ipm.TALKINGTON_PARAMS["t_stimulation"].value == pytest.approx(0.1245)
    assert ipm.TALKINGTON_PARAMS["t_death"].value == pytest.approx(0.0412)
    # The kill rate is mouse-literature-cited, not fit here.
    assert ipm.TALKINGTON_PARAMS["kill_rate"].provenance == "cited"
    # I is the swept control knob.
    assert ipm.TALKINGTON_PARAMS["icb_efficiency"].provenance == "varied"


def test_kill_rate_to_per_step_probability_form():
    """p = 1 - exp(-r2 * density * dt): monotone in each arg, in [0,1), and
    exactly 0 when any factor is 0."""
    r2 = 1.101e-7
    assert ipm.kill_rate_to_per_step_probability(r2, 0.0, 0.25) == 0.0
    assert ipm.kill_rate_to_per_step_probability(r2, 1.0e6, 0.0) == 0.0
    p_small = ipm.kill_rate_to_per_step_probability(r2, 1.0e5, 0.25)
    p_big = ipm.kill_rate_to_per_step_probability(r2, 1.0e7, 0.25)
    assert 0.0 < p_small < p_big < 1.0
    # Matches the closed form exactly.
    assert ipm.kill_rate_to_per_step_probability(r2, 1.0e6, 0.25) == pytest.approx(
        1.0 - math.exp(-r2 * 1.0e6 * 0.25)
    )


def test_kill_rate_rejects_negative_inputs():
    with pytest.raises(ValueError):
        ipm.kill_rate_to_per_step_probability(1e-7, -1.0, 0.25)
    with pytest.raises(ValueError):
        ipm.kill_rate_to_per_step_probability(1e-7, 1.0, -0.25)


def test_icb_efficiency_maps_to_residual_and_clamps():
    """I=0 perfect blockade -> residual 0; I=1 no blockade -> residual 1; the
    80-90%-blockade optimum (I~0.1-0.2) is a low residual; out-of-range clamps."""
    assert ipm.icb_efficiency_to_checkpoint_residual(0.0) == 0.0
    assert ipm.icb_efficiency_to_checkpoint_residual(1.0) == 1.0
    assert ipm.icb_efficiency_to_checkpoint_residual(0.15) == pytest.approx(0.15)
    assert ipm.icb_efficiency_to_checkpoint_residual(-0.5) == 0.0
    assert ipm.icb_efficiency_to_checkpoint_residual(2.0) == 1.0


def test_exhaustion_threshold_to_rate():
    """exhaustion_rate ~ 1/n; with n=1e4 that is 1e-4 per interaction."""
    assert ipm.exhaustion_threshold_to_rate(1.0e4) == pytest.approx(1.0e-4)
    # Larger threshold -> slower exhaustion (smaller rate).
    assert ipm.exhaustion_threshold_to_rate(2.0e4) < ipm.exhaustion_threshold_to_rate(1.0e4)
    with pytest.raises(ValueError):
        ipm.exhaustion_threshold_to_rate(0.0)


def test_mapping_table_flags_repo_specific_layers():
    """The mapping must honestly mark the spatial/DAMP/IFN-gamma layers as
    repo-specific (no external anchor), not silently imply they were calibrated."""
    notes = " ".join(note for _, _, note in ipm.MAPPING).lower()
    assert "repo-specific" in notes
    # The dimensional kill rate must be flagged as not directly transplantable.
    assert "not transplantable" in notes or "double-count" in notes
