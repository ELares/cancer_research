"""The "unconstrained" flag must be judged against a null, not a constant.

WHAT WAS WRONG
--------------
`abc_joint_posterior.py` flagged a parameter unconstrained when its 95% posterior
interval still spanned at least 0.6 of the prior width. That threshold ignores
how many draws were accepted, and that is precisely what sets how narrow an
UNINFORMATIVE posterior looks: with 30 accepted draws, samples taken from the
prior and nothing else still span a median of ~0.90 of the prior width, because
30 points rarely reach the corners.

So 0.6 sat far below anything noise produces, and the flag fired on WELL
determined parameters. It labelled `lp_propagation`, `lp_rate` and `gpx4_rate`
unconstrained while they sit at the 0th percentile of that null — and those are
the cascade parameters the manuscript quotes credible intervals for. The
mislabel had propagated into `joint-posterior.md` ("the parameters the in-vitro
dose-response panel does not identify") and into CLAUDE.md ("the LP cascade ...
stay loosely constrained"), both asserting the opposite of the truth.

WHY THE NULL IS THE RIGHT REFERENCE
-----------------------------------
"Posterior width as a fraction of prior width" has no fixed scale. Its
uninformative value depends only on the number of accepted draws, so the honest
comparison is against that distribution rather than against a number someone
picked. These guards pin the null itself, so a change to the sampler cannot
quietly move what counts as informed.
"""

import json
import math
import random
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import abc_posterior_information as info  # noqa: E402
import identifiability_report as ident  # noqa: E402

CAL = REPO_ROOT / "analysis" / "calibration"
INFO_JSON = CAL / "abc-information-content.json"
INFO_MD = CAL / "abc-information-content.md"
JOINT_MD = CAL / "joint-posterior.md"
GENERATOR = REPO_ROOT / "scripts" / "abc_joint_posterior.py"


def reports() -> list:
    return json.loads(INFO_JSON.read_text())


def test_the_null_is_what_an_uninformative_posterior_looks_like():
    """Independent re-derivation, so the reported null is not self-certifying.

    Computed here with a different seed and a different generator than the
    analysis uses. If the two disagreed, the null would be an artifact of one
    RNG rather than a property of the statistic.
    """
    r = next(r for r in reports() if r.get("parameters"))
    n = r["n_accepted"]
    rng = random.Random(12345)
    widths = []
    for _ in range(4000):
        u = sorted(rng.random() for _ in range(n))
        k_hi = (len(u) - 1) * 0.975
        k_lo = (len(u) - 1) * 0.025
        hi = u[int(k_hi)] + (u[min(int(k_hi) + 1, len(u) - 1)] - u[int(k_hi)]) * (k_hi - int(k_hi))
        lo = u[int(k_lo)] + (u[min(int(k_lo) + 1, len(u) - 1)] - u[int(k_lo)]) * (k_lo - int(k_lo))
        widths.append(hi - lo)
    widths.sort()
    median = widths[len(widths) // 2]
    assert abs(median - r["null_median_width"]) < 0.02, (
        f"independent null median {median:.3f} disagrees with the reported "
        f"{r['null_median_width']}")
    # The property that makes 0.6 wrong: noise does NOT look like 0.6.
    assert median > 0.8, (
        f"an uninformative posterior of {n} draws spans {median:.3f} of the "
        "prior; if this ever fell near 0.6 the original threshold would have "
        "been defensible and this analysis would need revisiting")


def test_the_legacy_threshold_mislabels_informed_parameters():
    """Test the retired rule's defect independently of a run's validity."""
    null = sorted(info.null_widths(30, replicates=4000))
    width = 0.65
    assert width >= info.LEGACY_THRESHOLD
    assert width < info._pct(null, 5)


def test_joint_information_is_withheld_when_the_current_run_is_underpowered():
    """A sample-size diagnostic cannot license intervals from four draws."""
    joint = [r for r in reports() if "joint-posterior" in r["artifact"]][0]
    source = json.loads((CAL / "joint-posterior.json").read_text())
    assert source["min_posterior"] == 20
    if source["underpowered"] or source["n_accepted"] < source["min_posterior"]:
        assert joint["underpowered"] is True
        assert "underpowered" in joint["unassessable"]
        assert "parameters" not in joint and "null_p5_width" not in joint
    else:
        assert joint.get("parameters") and not joint.get("unassessable")


def test_assessable_parameter_flags_follow_the_null_criterion():
    for report in reports():
        if report.get("parameters"):
            for parameter in report["parameters"].values():
                assert parameter["informed"] == (parameter["null_percentile"] <= 5.0)


def _posterior_fixture(path, n=4, underpowered=True):
    path.write_text(json.dumps({
        "n_draws": 40000, "n_accepted": n, "min_posterior": 20,
        "underpowered": underpowered,
        "posterior": {
            "narrow": {"posterior_width_frac_of_prior": 0.1},
            "wide": {"posterior_width_frac_of_prior": 0.99},
        },
    }))


@pytest.mark.parametrize("n,flag", [(4, True), (4, False), (20, True)])
def test_underpowered_status_or_count_blocks_information_before_null(n, flag, tmp_path, monkeypatch):
    path = tmp_path / "joint-posterior.json"
    _posterior_fixture(path, n, flag)
    monkeypatch.setattr(info, "null_widths", lambda *a, **k: pytest.fail("invalid posterior reached null test"))
    result = info.assess(path)
    assert result["underpowered"] is True and result["n_accepted"] == n
    assert "parameters" not in result and "null_median_width" not in result
    rendered = info.render([result])
    assert "Not assessed: underpowered" in rendered
    assert "parameters are informed" not in rendered


def test_minimum_accepted_count_permits_informative_and_uninformative_results(tmp_path):
    path = tmp_path / "joint-posterior.json"
    _posterior_fixture(path, n=20, underpowered=False)
    result = info.assess(path)
    assert "unassessable" not in result
    assert result["parameters"]["narrow"]["informed"] is True
    assert result["parameters"]["wide"]["informed"] is False


@pytest.mark.parametrize("current", ["underpowered", "missing"])
def test_identifiability_does_not_reuse_stale_joint_information(current, tmp_path, monkeypatch):
    cal = tmp_path / "analysis" / "calibration"
    cal.mkdir(parents=True)
    (cal / "abc-information-content.json").write_text(json.dumps([{
        "artifact": "analysis/calibration/joint-posterior.json",
        "parameters": {"old_parameter": {"informed": True}},
    }]))
    if current == "underpowered":
        _posterior_fixture(cal / "joint-posterior.json")
    monkeypatch.setattr(ident, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(ident, "OUT_MD", tmp_path / "identifiability.md")
    facts = ident._fitted_cascade_facts()
    assert "informed" not in facts and "uninformed" not in facts
    assert facts["joint_information"]["unassessable"]
    report = ident.build()
    assert "parameters" not in report["joint_posterior_information"]
    ident.write_report(report)
    text = ident.OUT_MD.read_text()
    assert "Current joint inference" in text and "NOT ASSESSED" in text
    assert "1 of its 1 parameters" not in text
    assert "including the whole LP cascade" not in text


def test_the_generator_no_longer_uses_a_bare_constant():
    src = GENERATOR.read_text()
    assert 'post[n]["posterior_width_frac_of_prior"] >= 0.6' not in src, (
        "the generator is back to the fixed 0.6 threshold")
    assert "null_p5" in src, (
        "the generator does not compute a null-calibrated threshold")


def test_the_history_survives_regeneration():
    """The corrections were written INTO a generated file and were erased.

    joint-posterior.md is produced by abc_joint_posterior.py. Two correction
    notes were added directly to the markdown; the next run of the generator
    overwrote both without warning. The history now lives in the generator's own
    template, so regenerating reproduces it -- which is what this checks.
    """
    txt = JOINT_MD.read_text()
    source_result = json.loads((CAL / "joint-posterior.json").read_text())
    underpowered = (source_result["underpowered"] or
                    source_result["n_accepted"] < source_result["min_posterior"])
    history_heading = ("Acceptance and history" if underpowered else
                       "History of this run's acceptance rule")
    assert history_heading in txt, (
        "the acceptance-rule history is missing from joint-posterior.md; if it "
        "was written into the markdown directly it will be erased again on the "
        "next regeneration -- it belongs in the generator")
    src = (REPO_ROOT / "scripts" / "abc_joint_posterior.py").read_text()
    assert history_heading in src, (
        "the history is in the generated document but not in the generator, so "
        "it will not survive the next run")


def test_the_analysis_states_what_it_does_not_measure():
    txt = INFO_MD.read_text().lower()
    assert "not whether the model" in txt or "not whether the" in txt, (
        "the report no longer distinguishes 'the data moved this parameter' "
        "from 'the fit is right'")


# --- the identifiability report's closing section --------------------------

IDENT_MD = REPO_ROOT / "analysis" / "identifiability-report.md"
IDENT_GEN = REPO_ROOT / "scripts" / "identifiability_report.py"


def test_the_closing_section_separates_current_inference_from_historical_substitution():
    txt = IDENT_MD.read_text()
    tail = txt[txt.index("## What would make a headline point-estimable"):]
    assert "Current joint inference" in tail
    assert "Historical substitution experiment" in tail
    assert "not an admissibility result for the corrected joint" in tail
    assert "INADMISSIBLE" in tail, (
        "the closing section does not record the demonstrated reason the "
        "substitution route is closed")


def test_the_closing_figures_are_derived_from_the_artifacts():
    """Every number in that section must come from a committed artifact.

    They were hand-typed beside the two files that compute them, which is the
    shape this repository keeps rediscovering: the artifact moves, the sentence
    does not, and the stale figure reads as freshly checked.
    """
    import re
    txt = IDENT_MD.read_text()
    tail = txt[txt.index("## What would make a headline point-estimable"):]

    h = json.loads((CAL.parent / "headline-at-fitted-cascade.json").read_text())
    worst = max(("ctrpv2_point", "posterior_median"),
                key=lambda k: h[k]["admissibility"]["worst_rate"])
    a = h[worst]["admissibility"]
    assert f"{a['worst_rate']*100:.2f}%" in tail, (
        "the inadmissible untreated-death rate in the prose is not the one in "
        "headline-at-fitted-cascade.json")
    assert f"{h['default']['admissibility']['worst_rate']*100:.2f}%" in tail

    joint = [r for r in reports() if "joint-posterior" in r["artifact"]][0]
    if joint.get("parameters"):
        informed = sum(1 for v in joint["parameters"].values() if v["informed"])
        assert f"{informed} of its {len(joint['parameters'])} parameters" in tail
    else:
        assert joint["unassessable"] in tail
        assert "NOT ASSESSED" in tail

    # The generator must CALL the helper, not merely define it. Asserting the
    # name appears anywhere passed a mutation that deleted the call site, because
    # the function's own `def` line contains the name -- a guard satisfied by the
    # thing it is checking for.
    src = IDENT_GEN.read_text()
    body = src[src.index("def write_report("):]
    assert "_fitted_cascade_facts()" in body, (
        "write_report no longer calls _fitted_cascade_facts(); the closing "
        "figures are not being derived even if the helper still exists")
    for literal in ("99.96%", "5 of its 7 parameters are informed"):
        assert literal not in src, (
            f"{literal!r} is hardcoded in the generator again")


def test_the_disjunction_framing_carries_its_qualification():
    """The in-vivo PRCC ranges are +/-50% of the defaults under question.

    So the disjunction restates the falsification rather than supplying
    independent grounds to discount the fit. Without that, "DISJOINT" reads as
    stronger evidence than it is.
    """
    txt = IDENT_MD.read_text()
    assert "DISJOINT" in txt, "the disjunction is no longer reported at all"
    # NOT an `or`. The first version accepted either the explanation or the bare
    # "+/-50%", and the bare number appears elsewhere in the document -- so
    # deleting the sentence that does the explaining still passed.
    assert "restates the falsification" in txt, (
        "the report asserts DISJOINT without the circularity qualification that "
        "analysis/calibration/in-vivo-prior-provenance.md establishes: the "
        "in-vivo PRCC ranges are +/-50% bands around the defaults under "
        "question, so the disjunction restates the falsification rather than "
        "supplying independent grounds to discount the fit")
