"""Guards for P9 to P13 and `analysis/modality-predictions.{md,json}`.

A preregistration is worth exactly what its numbers are worth. These check the
document quotes what the engine derives, that each prediction is FALSIFIABLE
rather than decorative, and that the two directions the campaign got wrong
stay recorded.
"""
import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/modality-predictions.md"
JSON_ = REPO / "analysis/modality-predictions.json"
PREREG = REPO / "PREREGISTRATION.md"
NEW = ("P9", "P10", "P11", "P12", "P13")


def _load():
    spec = importlib.util.spec_from_file_location(
        "modality_predictions", REPO / "scripts/modality_predictions.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MP = _load()


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


def test_the_model_outputs_are_recomputed_not_stored(d):
    """A hand-edited prediction is a prediction about nothing."""
    live = MP.scan()
    assert live == d, "the stored model outputs drifted from the engine"


def test_every_new_prediction_is_registered_with_a_threshold():
    """A prediction with no falsification threshold cannot be scored.

    The document's own format requires both halves, and a preregistration that
    states only the model output is a description.
    """
    text = PREREG.read_text()
    for p in NEW:
        m = re.search(rf"^\*\*{p}\. (.+?)\*\*\n(.*?)(?=^\*\*P|\Z)",
                      text, re.M | re.S)
        assert m, f"{p} is not registered in PREREGISTRATION.md"
        body = m.group(2)
        assert "*Quantitative model output:*" in body, f"{p} has no model output"
        assert "*Falsification threshold:*" in body, (
            f"{p} states no falsification threshold, so nothing can score it")
        thr = body.split("*Falsification threshold:*", 1)[1]
        assert len(thr.split()) >= 12, (
            f"{p}'s threshold is too short to name an outcome: {thr[:80]}")


def test_the_registered_numbers_are_the_derived_ones(d):
    """P9 is the one quantitative prediction, so its figures must be live."""
    text = PREREG.read_text()
    p9 = d["P9"]
    lo_r, hi_r = p9["low_over_high_dose_ratio"]
    for frag in (f"{p9['published_band'][0]} to {p9['published_band'][1]} band",
                 f"{p9['boost_window'][0]} to {p9['boost_window'][1]}",
                 f"{lo_r} to {hi_r}"):
        assert frag in text, f"PREREGISTRATION.md does not quote {frag!r}"
    # and the direction it rests on must actually hold in the derived data
    assert lo_r > 1.0 and hi_r > 1.0, (
        "the model no longer predicts a larger SER at low dose; P9 must be "
        "rewritten rather than left standing")


def test_p10_is_registered_in_the_corrected_direction(d):
    """The campaign's own prior belief was wrong, and that is the point.

    The ADC module was built expecting the bystander advantage to GROW as
    antigen is lost and shipped a guard asserting it; the guard was satisfied
    by a constant. If this ever flips back, the retraction has to be revisited
    rather than the prediction quietly restated.
    """
    reach = d["P10"]["reach_of_negative_pool"]
    keys = sorted(reach, key=float, reverse=True)
    vals = [reach[k] for k in keys]
    assert vals == sorted(vals, reverse=True), (
        f"the bystander reach no longer FALLS as antigen is lost ({reach}); "
        "P10 is registered on that decline")
    text = PREREG.read_text()
    assert "STARVED by the antigen escape" in text
    assert "OPPOSITE sign to the intuition" in text, (
        "the preregistration no longer records that this contradicts the "
        "project's own prior belief")
    status = (REPO / "simulations/calibration/CALIBRATION_STATUS.md").read_text()
    assert "THAT GUARD WAS VACUOUS AND ITS CLAIM IS RETRACTED" in status


def test_the_scope_audit_reports_the_widened_denominator():
    """The measure this work exists to move."""
    audit = json.loads((REPO / "analysis/scope-audit.json").read_text())
    preds = audit["predictions"]
    assert len(preds) >= 13, f"only {len(preds)} predictions found"
    ferro = sum(1 for v in preds.values() if v)
    assert ferro < len(preds), (
        "every preregistered prediction still concerns ferroptosis or the "
        "physical-ROS modalities; the modality arms remain unfalsifiable")
    for p in NEW:
        assert preds.get(p) is False, (
            f"{p} is classified as a ferroptosis prediction, which would mean "
            "the audit cannot see the widening")


def test_the_page_refuses_to_present_placeholders_as_predictions():
    md = MD.read_text()
    for frag in ("Directional, every one of them",
                 "uncalibrated placeholders",
                 "A magnitude here is not a prediction",
                 "Not clinical guidance"):
        assert frag in md, f"the page no longer says: {frag}"
