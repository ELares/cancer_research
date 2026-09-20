"""Keep acceptance-rule compliance separate from posterior sample adequacy.

The original fixed-quota defect is historical. Current probes must be reported
as measured, and a favorable four-draw median must not hide an underpowered run.
"""

import importlib.util
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
CAL = REPO_ROOT / "analysis" / "calibration"
DIAG_JSON = CAL / "abc-acceptance-diagnostic.json"
DIAG_MD = CAL / "abc-acceptance-diagnostic.md"
JOINT = CAL / "joint-posterior.json"
GENERATOR = REPO_ROOT / "scripts" / "abc_joint_posterior.py"


def diag() -> dict:
    return json.loads(DIAG_JSON.read_text())


def renderer():
    spec = importlib.util.spec_from_file_location(
        "abc_acceptance_diagnostic", REPO_ROOT / "scripts" / "abc_acceptance_diagnostic.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.render


def test_the_rule_fix_does_not_hide_the_recorded_sampling_shortfall():
    d = diag()
    joint = json.loads(JOINT.read_text())
    assert d["acceptance_rule"] == joint["acceptance_rule"] == "tolerance"
    assert d["n_accepted_now"] == joint["n_accepted"]
    assert d["min_posterior"] == joint["min_posterior"]
    assert d["underpowered"] == joint["underpowered"]
    assert d["sampling"]["draws_of_record"] == joint["n_draws"]
    assert d["reported_epsilon"] == joint["epsilon_joint_distance"]
    assert d["target_source"] == joint["target_source"]
    assert d["target_dose_grids_um"] == {
        "ML162": joint["curves"]["rsl3_doses_um"],
        "ERASTIN": joint["curves"]["erastin_doses_um"],
    }
    text = DIAG_MD.read_text()
    assert "Status: ACCEPTANCE RULE FIXED" in text
    if d["underpowered"] or d["n_accepted_now"] < d["min_posterior"]:
        assert "POSTERIOR UNDERPOWERED" in text
        assert "not reliable\nposterior estimates" in text
    assert renderer()(d) == text


def test_boundary_probe_conclusion_follows_the_measured_distances():
    t = diag()["prior_truncation_test"]
    at_bound = t["k_erastin"]["3.0"]
    outside = [v for k, v in t["k_erastin"].items() if float(k) < 3.0]
    assert outside, "the k_erastin truncation test no longer probes outside the bound"
    text = DIAG_MD.read_text()
    expected = ("tested lower `k_erastin` values do not improve" if all(
        v >= at_bound for v in outside) else
        "`k_erastin` value outside the prior improves")
    assert expected in text
    assert "do not search for a global optimum outside the prior" in text


def test_hill_probe_does_not_claim_inertness_throughout_the_prior():
    hl = diag()["prior_truncation_test"]["hill"]
    text = DIAG_MD.read_text()
    if len(set(hl.values())) == 1:
        assert "unchanged at the probed values 6, 8 and 10" in text
        assert "does not establish that the parameter is inert throughout its prior" in text
    else:
        assert "distances differ across the tested values" in text


def test_fresh_search_reports_its_rate_without_augmenting_the_posterior():
    s_ = diag()["sampling"]
    rate = s_["n_beating_committed"] / s_["draws"]
    assert 0 <= s_["n_beating_committed"] <= s_["draws"]
    assert 0 <= s_["n_inside_epsilon"] <= s_["draws"]
    assert s_["frac_inside_epsilon"] == round(s_["n_inside_epsilon"] / s_["draws"], 4)
    assert s_["draws"] >= 10000, (
        "the rate is measured on too few draws to distinguish 'rare' from "
        "'absent' -- the earlier 300-draw sample reported zero and was wrong")
    text = DIAG_MD.read_text()
    assert f"about {rate:.1e} per draw" in text
    assert "not added to the joint posterior" in text


def test_epsilon_is_now_anchored_rather_than_floating():
    """Epsilon must sit NEAR the achievable distance, not far above it.

    This asserted the gap was LARGE, which was the defect's signature. Under the
    fix epsilon is the reference times the tolerance factor, so the same quantity
    is now small and positive -- the assertion had to invert with the thing it
    measures.
    """
    d = diag()
    gap = d["epsilon_excess_over_committed"]
    assert gap is not None, "the epsilon-to-reference gap is no longer recorded"
    assert 0 < gap <= 0.15, (
        f"epsilon is {100*gap:.0f}% above the achievable distance; anchored "
        "acceptance should keep it within the tolerance factor")


def test_the_report_says_what_it_does_not_overturn():
    """It must not read as 'the posterior is worthless'."""
    txt = DIAG_MD.read_text().lower()
    assert "does not claim" in txt or "does not overturn" in txt, (
        "the report no longer bounds what the fix does and does not establish")
    assert "not a validated one" in txt, (
        "the report no longer says a better-fitting posterior is not a validated "
        "one -- the overclaim this whole arc is most exposed to")
    assert "inadmissible" in txt, (
        "the report no longer points at the in-vivo regime where substituting "
        "these values remains inadmissible")


def test_the_acceptance_rule_is_a_tolerance_not_a_quantile():
    """The fix itself. This guard used to assert the OPPOSITE.

    It was written as "the rule is still a fixed fraction -- if that ever changes
    the document must be reworded", the right shape for a guard on a known-live
    defect. Fixing the rule fired it, the document was reworded, and it now pins
    the fix.
    """
    src = GENERATOR.read_text()
    assert "n_accept = max(10, int(args.n_draws * ACCEPT_FRAC))" not in src, (
        "the generator is back to a fixed-fraction acceptance")
    assert "distances <= eps" in src, "acceptance is no longer a tolerance"
    assert "REFERENCE_VECTOR" in src and "TOLERANCE_FACTOR" in src, (
        "the tolerance is no longer anchored to a reachable reference distance")
    assert "underpowered" in src, (
        "the run no longer reports a shortfall; padding back to a quota is the "
        "behaviour that produced the original defect")


@pytest.mark.parametrize("median_distance", [0.05, 0.3])
@pytest.mark.parametrize("reported_underpowered", [True, False])
def test_four_draw_median_cannot_override_a_sampling_shortfall(
        median_distance, reported_underpowered):
    d = diag()
    d.update(acceptance_rule="tolerance", n_accepted_now=4, min_posterior=20,
             underpowered=reported_underpowered, committed_distance=0.2,
             posterior_median_distance=median_distance)
    text = renderer()(d)
    assert "Status: ACCEPTANCE RULE FIXED; POSTERIOR UNDERPOWERED" in text
    assert "4 accepted draws; recorded minimum\n20" in text
    assert "not reliable\nposterior estimates" in text
    assert "favorable median-vector distance cannot resolve" in text


def test_recorded_shortfall_is_not_overridden_by_count_alone():
    d = diag()
    d.update(acceptance_rule="tolerance", n_accepted_now=25, min_posterior=20,
             underpowered=True)
    assert "POSTERIOR UNDERPOWERED" in renderer()(d)


def test_meeting_minimum_is_not_presented_as_validating_the_posterior():
    d = diag()
    d.update(acceptance_rule="tolerance", n_accepted_now=25, min_posterior=20,
             underpowered=False)
    text = renderer()(d)
    assert "POSTERIOR UNDERPOWERED" not in text
    assert "accepted count meets the recorded minimum of 20" in text
    assert "minimum alone does not establish Monte Carlo stability" in text


@pytest.mark.parametrize("missing", ["underpowered", "min_posterior"])
def test_missing_metadata_cannot_establish_sampling_adequacy(missing):
    d = diag()
    d.update(acceptance_rule="tolerance", n_accepted_now=25, min_posterior=20,
             underpowered=False)
    d.pop(missing)
    assert "SAMPLING ADEQUACY UNKNOWN" in renderer()(d)
