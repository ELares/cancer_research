"""Guards for the #408 reaction-diffusion vs BioFVM external cross-check.

The cross-check itself needs a local BioFVM build + cargo and is NOT run in CI;
this validates the committed result (analysis/calibration/rd-biofvm-crosscheck.json)
so the agreement claim cannot silently rot.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CROSSCHECK = REPO_ROOT / "analysis" / "calibration" / "rd-biofvm-crosscheck.json"


def _load():
    return json.loads(CROSSCHECK.read_text())


def test_structure():
    r = _load()
    assert r["single_source_dt_sweep"]
    assert "two_source" in r
    for s in r["single_source_dt_sweep"]:
        assert {"dt", "shape_log_pearson_r", "rust_decay_length_um",
                "biofvm_decay_length_um", "median_ratio_rust_over_biofvm"} <= set(s)


def test_shape_agrees_everywhere():
    """The load-bearing finding: the two independent solvers produce the SAME field
    shape (log-field Pearson r > 0.99) in every run, single- and two-source."""
    r = _load()
    for s in r["single_source_dt_sweep"]:
        assert s["shape_log_pearson_r"] > 0.99
    assert r["two_source"]["shape_log_pearson_r"] > 0.99


def test_decay_length_matches():
    """Both solvers reproduce the same effective decay length (the physics matches)."""
    r = _load()
    for s in r["single_source_dt_sweep"]:
        assert abs(s["rust_decay_length_um"] - s["biofvm_decay_length_um"]) < 3.0


def test_magnitude_ratio_converges_to_one_as_dt_falls():
    """The only difference is BioFVM's LOD operator-splitting error: the Rust/BioFVM
    magnitude ratio moves toward 1 as BioFVM's dt shrinks, so the residual is a
    controllable numerical artifact, not a solver disagreement."""
    sweep = sorted(_load()["single_source_dt_sweep"], key=lambda s: -s["dt"])  # large dt -> small
    ratios = [s["median_ratio_rust_over_biofvm"] for s in sweep]
    # monotone improvement toward 1 as dt falls
    assert all(abs(b - 1.0) <= abs(a - 1.0) + 1e-9 for a, b in zip(ratios, ratios[1:]))
    # at the finest dt the solvers agree to within ~10%
    assert abs(ratios[-1] - 1.0) < 0.12


def test_two_source_multivessel_agreement():
    """The agreement holds in the multi-vessel case the nearest-vessel proxy averages
    away (the regime #343 motivates)."""
    t = _load()["two_source"]
    assert t["shape_log_pearson_r"] > 0.99
    assert abs(t["rust_decay_length_um"] - t["biofvm_decay_length_um"]) < 3.0
