"""Report claims must follow the target definition and the recorded inference."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CAL = REPO / "analysis" / "calibration"
sys.path.insert(0, str(REPO / "scripts"))


def _module(name):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _artifact(name):
    return json.loads((CAL / f"{name}.json").read_text())


def _report(module, result, tmp_path, monkeypatch):
    # Older committed fixtures predate dose-support metadata. Its construction
    # has separate tests; these checks exercise the interpretation of the fit.
    result.setdefault("target_support", {})
    output = tmp_path / "report.md"
    monkeypatch.setattr(module, "OUT_MD", output)
    module.write_report(result)
    return output.read_text()


@pytest.mark.parametrize("underpowered", [False, True])
def test_joint_report_states_tolerance_and_sample_size_null(
        tmp_path, monkeypatch, underpowered):
    module = _module("abc_joint_posterior")
    result = _artifact("joint-posterior")
    result["underpowered"] = underpowered
    result["n_accepted"] = result["min_posterior"] - int(underpowered)
    text = _report(module, result, tmp_path, monkeypatch)
    assert f"**{result['reference_distance']}**" in text
    assert f"**{result['tolerance_factor']}**" in text
    assert "fraction is an outcome, not a quota" in text
    assert f"**{result['uninformative_null_p5_width']}**" in text
    assert "former fixed 0.6 threshold is no longer used" in text
    assert ("**Underpowered:**" in text) == underpowered
    assert "accepting the closest" not in text


def test_joint_report_does_not_assert_separation_when_interval_overlaps(tmp_path, monkeypatch):
    module = _module("abc_joint_posterior")
    result = _artifact("joint-posterior")
    result["underpowered"] = False
    result["n_accepted"] = result["min_posterior"]
    result["disjunction_with_invivo_priors"]["lp_rate"][
        "entire_95pct_posterior_above_invivo_max"] = False
    text = _report(module, result, tmp_path, monkeypatch)
    assert "does not lie entirely above" in text
    assert "Both reported 95% intervals lie above" not in text


@pytest.mark.parametrize("name,artifact", [
    ("abc_joint_posterior", "joint-posterior"),
    ("abc_posterior", "abc-posterior"),
])
def test_parameter_draw_bands_do_not_claim_experimental_coverage(
        tmp_path, monkeypatch, name, artifact):
    result = _artifact(artifact)
    if name == "abc_joint_posterior":
        result.update(underpowered=False, n_accepted=result["min_posterior"])
    text = _report(_module(name), result, tmp_path, monkeypatch)
    assert "parameter-draw band" in text
    assert "experimental outcomes" in text
    assert "not independent" in text
    assert "cell lines" in text


def test_single_abc_report_does_not_praise_zero_heldout_coverage(tmp_path, monkeypatch):
    module = _module("abc_posterior")
    result = _artifact("abc-posterior")
    n = len(result["posterior_predictive_heldout"]["dose_um"])
    result["heldout_coverage_strict"] = f"0/{n}"
    result["heldout_coverage_tolerant"] = f"0/{n}"
    text = _report(module, result, tmp_path, monkeypatch)
    assert f"**0/{n}**" in text
    assert "in the right place" not in text
    assert "epsilon is an output" in text


def test_point_fit_reports_population_target_without_claiming_no_overfit(tmp_path, monkeypatch):
    module = _module("calibrate_kill_switch")
    result = _artifact("kill-switch-calibration")
    result["heldout_rmse"] = 2 * result["default_uncalibrated_rmse"]
    text = _report(module, result, tmp_path, monkeypatch)
    assert str(result["heldout_rmse"]) in text
    assert "pointwise median of fitted viability curves" in text
    assert "target need not describe any individual cell line" in text
    assert "does not rule out overfitting" in text
    assert "only ~1.4x" not in text


def test_erastin_report_does_not_infer_biological_ceiling_or_target_specificity(
        tmp_path, monkeypatch):
    module = _module("calibrate_erastin")
    text = _report(module, _artifact("erastin-calibration"), tmp_path, monkeypatch)
    assert "not experimental evidence of target specificity" in text
    assert "does not measure a high-dose biological ceiling" in text
    assert "pointwise medians" in text
    assert "below the measured erastin ceiling" not in text


@pytest.mark.parametrize("current_reference", [0.1, 0.5])
def test_diagnostic_keeps_historical_reference_on_its_original_objective(current_reference):
    module = _module("abc_acceptance_diagnostic")
    result = _artifact("abc-acceptance-diagnostic")
    result["acceptance_rule"] = "tolerance"
    result["committed_distance"] = current_reference
    text = module.render(result)
    history = text.split("## Historical run: original targets")[1].split(
        "## Current-target check")[0]
    assert str(module.HISTORICAL["reference_distance"]) in history
    assert "not directly comparable" in history
    assert "went from **losing**" not in text
    assert "| | before | after |" not in text


def test_diagnostic_does_not_treat_median_vector_as_an_accepted_draw():
    module = _module("abc_acceptance_diagnostic")
    result = _artifact("abc-acceptance-diagnostic")
    result["acceptance_rule"] = "tolerance"
    result["posterior_median_distance"] = result["committed_distance"] + 0.1
    text = module.render(result)
    assert "Status: ACCEPTANCE RULE FIXED" in text
    assert "does not fit better than the reference" in text
    assert "median need not be an accepted vector" in text
    assert "did not test\n  the corrected posterior" in text


def test_information_report_labels_legacy_flag_and_limits_width_inference():
    module = _module("abc_posterior_information")
    text = module.render(_artifact("abc-information-content"))
    assert "no longer the joint generator's rule" in text
    assert "legacy 0.6 would flag" in text
    assert "width test does not rule out shifts" in text
    assert "are flagged *unconstrained* by the generator" not in text
    assert "Their reported credible intervals are the prior's" not in text
    assert "did not test the corrected posterior" in text


@pytest.mark.parametrize("flag,count", [(True, 20), (False, 4), (True, 4)])
def test_joint_withholds_inference_when_either_sampling_guard_fails(
        tmp_path, monkeypatch, flag, count):
    module = _module("abc_joint_posterior")
    result = _artifact("joint-posterior")
    result.update(underpowered=flag, n_accepted=count, min_posterior=20)
    text = _report(module, result, tmp_path, monkeypatch)
    assert "**Underpowered:**" in text
    assert "claims are withheld" in text
    assert "| parameter | 2.5%" not in text
    assert "Both reported 95% intervals lie above" not in text
    # Plot dispatch must not access draws or construct intervals when invalid.
    plotted = []
    monkeypatch.setattr(module, "_plot_underpowered", lambda r: plotted.append(r))
    module._plot(result, None, None, None, None)
    assert plotted == [result]
