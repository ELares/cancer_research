"""Guards for the #464 PDT/SDT exo-ROS kill-threshold validation (Zhu 2015).

Validates the committed result and re-runs the drift guard so the validated claim
(the model's single source-independent death threshold is supported by the measured
photosensitizer-independence of the singlet-oxygen kill dose) cannot silently rot.
Pure stdlib; runs in CI.
"""
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULT = REPO_ROOT / "analysis" / "calibration" / "pdt-threshold-validation.json"
sys.path.insert(0, str(REPO_ROOT / "scripts"))


def test_committed_result_passes():
    r = json.loads(RESULT.read_text())
    assert r["model_uses_single_source_independent_threshold"] is True
    assert r["measured_is_photosensitizer_independent"] is True
    assert r["validation"].startswith("PASS")
    # Measured threshold ~0.5 mM, photosensitizer-independent within a factor of 2.
    assert 0.3 <= r["measured_threshold_mean_mm"] <= 0.8
    assert r["measured_max_over_min_ratio"] < 2.0


def test_drift_guard_model_constants_match():
    # Re-run the validation live against params.rs so a future change to sdt_ros /
    # pdt_ros / death_threshold that broke source-independence would fail here.
    import validate_pdt_threshold as v

    live = v.validate()
    # The two exogenous-ROS sources must stay equal (the source-independence the
    # Zhu-2015 data supports).
    assert live["model_sdt_ros"] == live["model_pdt_ros"]
    assert live["model_uses_single_source_independent_threshold"] is True
    # And the committed artifact must match the live parse.
    committed = json.loads(RESULT.read_text())
    assert live["model_sdt_ros"] == committed["model_sdt_ros"]
    assert live["model_death_threshold"] == committed["model_death_threshold"]


def test_o2_dependence_honesty_flagged():
    # The doc/result must keep the honest caveat that the linear O2-dependence is a
    # first-order approximation of Zhu's Michaelis form (constant not fabricated).
    r = json.loads(RESULT.read_text())
    note = r["o2_dependence_note"].lower()
    assert "linear" in note and "michaelis" in note
    assert "fabricated" in note
