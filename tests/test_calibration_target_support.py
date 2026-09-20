"""Integration guards from screened CTRPv2 ranges to fitted result artifacts.

These checks reconstruct targets from the committed data. They do not rerun
the simulation or mistake stored fit metrics for independent biological tests.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
CALIB = REPO / "analysis" / "calibration"
CURVES = CALIB / "ctrpv2_ferroptosis_curves.csv"
sys.path.insert(0, str(REPO / "scripts"))

import calibrate_erastin as ce  # noqa: E402
import calibrate_kill_switch as ck  # noqa: E402
from ctrp_dose_support import select_supported_cohort  # noqa: E402


def _curve(model, compound="ERASTIN", **overrides):
    return {"ModelID": model, "SampleID": f"sample-{model}",
            "CompoundName": compound, "MinimumDose": "0.1",
            "MaximumDose": "30", "DoseUnit": "uM",
            "LowerAsymptote": "0", "UpperAsymptote": "1",
            "EC50": "1", "Slope": "-2", **overrides}


@pytest.mark.parametrize("fit", [ck, ce], ids=["gpx4-point-fit", "erastin-point-fit"])
def test_unsupported_point_fit_fails_before_simulation_or_output(fit, tmp_path, monkeypatch):
    path = tmp_path / "curves.csv"
    rows = [_curve(name, compound=name) for name in ("ML162", "ML210", "ERASTIN")]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    def unexpected_simulation(*args, **kwargs):
        pytest.fail("unsupported dose grid reached a simulation or parameter search")

    monkeypatch.setattr(fit, "DOSE_GRID_UM", (0.1, 100.0))
    monkeypatch.setattr(ck, "_fc", unexpected_simulation)
    monkeypatch.setattr(fit, "grid_search", unexpected_simulation)
    outputs = [tmp_path / "result.json", tmp_path / "result.md"]
    for key, output in zip(("OUT_JSON", "OUT_MD"), outputs):
        output.write_text("existing result must survive\n")
        monkeypatch.setattr(fit, key, output)

    with pytest.raises(ValueError, match="100 uM: 0"):
        fit.run(SimpleNamespace(curves=path))
    assert all(p.read_text() == "existing result must survive\n" for p in outputs)


def test_empirical_target_uses_one_cohort_without_partial_support_outliers():
    shared = _curve("whole-grid")
    rows = [
        _curve("low-only", MaximumDose="1", LowerAsymptote="100", UpperAsymptote="100"),
        shared,
        _curve("high-only", MinimumDose="1", LowerAsymptote="-100", UpperAsymptote="-100"),
    ]
    doses = (0.1, 1, 10)
    medians, support = ck.empirical_target(rows, doses)
    # Hand-computed logistic for the single common row: 1 / (1 + dose**2).
    assert medians == pytest.approx([100 / 101, 0.5, 1 / 101])
    assert support["n_retained"] == 1 and support["n_excluded"] == 2
    assert [item["n_supported"] for item in support["per_dose_support"]] == [2, 3, 2]
    assert ck.empirical_median_viability(rows, doses) == medians
    assert ck.empirical_target([shared], doses)[0] == medians


@pytest.mark.parametrize("field,value", [
    ("EC50", "nan"), ("EC50", "inf"), ("EC50", "0"), ("EC50", "-1"),
    ("LowerAsymptote", "nan"), ("UpperAsymptote", "inf"), ("Slope", "nan"),
])
def test_invalid_curve_coefficient_fails_before_fitting(field, value, tmp_path, monkeypatch):
    row = _curve("bad-coefficient", **{field: value})
    path = tmp_path / "invalid-coefficient.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)

    def unexpected_fit(*args, **kwargs):
        pytest.fail("invalid target coefficient reached the fit")

    monkeypatch.setattr(ce, "DOSE_GRID_UM", (0.1, 1.0))
    monkeypatch.setattr(ce, "grid_search", unexpected_fit)
    monkeypatch.setattr(ck, "_fc", unexpected_fit)
    output = tmp_path / "result.json"
    monkeypatch.setattr(ce, "OUT_JSON", output)
    with pytest.raises(ValueError, match="finite|EC50"):
        ce.run(SimpleNamespace(curves=path))
    assert not output.exists()


def test_finite_curve_values_cannot_overflow_the_target_median():
    rows = [_curve(str(i), LowerAsymptote="1e308", UpperAsymptote="1e308")
            for i in range(2)]
    with pytest.raises(ValueError, match="nonfinite median"):
        ck.empirical_target(rows, (0.1, 1))


def test_fitted_viability_above_one_is_preserved():
    # Assay normalization can put a fitted upper asymptote above one. Finiteness
    # validation must not convert it to a silently clipped target.
    medians, _ = ck.empirical_target([_curve("valid", UpperAsymptote="1.1")], (0.1, 1))
    assert medians == pytest.approx([1.1 / 1.01, 0.55])
    assert medians[0] > 1


@pytest.mark.parametrize("dose,ec50,expected", [
    (1e-300, 1e300, 1.0), (1e300, 1e-300, 0.0),
])
def test_finite_extreme_dose_ratios_reach_the_correct_asymptote(dose, ec50, expected):
    assert ck.ctrp_viability(dose, 0, 1, ec50, -2) == expected


def _stored_targets(artifact_name, result):
    """Map each compound to its recorded grid and any stored target array."""
    if artifact_name == "kill-switch-calibration.json":
        doses = result["dose_grid_um"]
        curves = result["curves"]
        return {"ML162": (doses, curves["empirical_fit"]),
                "ML210": (doses, curves["empirical_heldout"]),
                "ERASTIN": (doses, curves["empirical_cross"])}
    if artifact_name == "erastin-calibration.json":
        return {"ERASTIN": (result["dose_grid_um"], result["curves"]["empirical"])}
    if artifact_name == "abc-posterior.json":
        heldout = result["posterior_predictive_heldout"]
        # The single-inducer artifact stores only the held-out target array.
        # Training support is nevertheless checked against the full source CSV.
        return {"ML162": (heldout["dose_um"], None),
                "ML210": (heldout["dose_um"], heldout["empirical_heldout"])}
    curves = result["curves"]
    return {"ML162": (curves["rsl3_doses_um"], curves["empirical_rsl3_ml162"]),
            "ERASTIN": (curves["erastin_doses_um"], curves["empirical_erastin"]),
            "ML210": (curves["rsl3_doses_um"],
                      result["heldout_posterior_predictive"]["empirical"])}


@pytest.mark.parametrize("artifact_name", [
    "kill-switch-calibration.json", "erastin-calibration.json",
    "abc-posterior.json", "joint-posterior.json",
])
def test_committed_artifact_targets_reconstruct_from_supported_source(artifact_name):
    result = json.loads((CALIB / artifact_name).read_text())
    assert "target_source" in result, f"{artifact_name}: regenerate with source provenance"
    source = result["target_source"]
    assert source["file"] == CURVES.name
    assert source["sha256"] == hashlib.sha256(CURVES.read_bytes()).hexdigest()
    assert "fitted" in source["target_kind"]
    targets = _stored_targets(artifact_name, result)
    assert set(result["target_support"]) == set(targets)
    rows = ck.load_curves(CURVES)

    for compound, (doses, stored_medians) in targets.items():
        expected_grid = (ce.DOSE_GRID_UM if compound == "ERASTIN" and
                         artifact_name != "kill-switch-calibration.json"
                         else ck.DOSE_GRID_UM)
        assert doses == list(expected_grid), f"{artifact_name}: stale {compound} grid"
        retained, support = select_supported_cohort(rows[compound], doses)
        assert result["target_support"][compound] == support
        assert len(retained) == support["n_retained"] > 0
        if stored_medians is not None:
            reconstructed, _ = ck.empirical_target(rows[compound], doses)
            # Generators store four decimals; allow only their rounding error.
            assert stored_medians == pytest.approx(reconstructed, rel=0, abs=0.000051), (
                f"{artifact_name}: {compound} target no longer matches supported curves")
