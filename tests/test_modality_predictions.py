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


def test_every_constant_is_read_from_the_crate_not_from_this_script():
    """The finding that made all the others possible.

    Every P10-P13 constant used to be a Python literal, and `scan()` even read
    `adc.rs` and threw the result away, which made the script LOOK derived. A
    reviewer changed four literals so they contradicted the Rust, regenerated,
    and published a preregistration disagreeing with its own cited module on
    every P10 and P11 number -- with all six guards and the freshness gate
    green. This re-reads the crate INDEPENDENTLY of the script's helpers.
    """
    core = REPO / "simulations/ferroptosis-core/src"
    d = json.loads(JSON_.read_text())

    adc = core / "adc.rs"
    fixture = re.search(r"fn heterogeneous\(\) -> AdcConfig \{(.*?)\n    \}",
                        adc.read_text(), re.S)
    assert fixture, "adc.rs no longer defines the heterogeneous() fixture"
    for field, key in (("payload_escape_fraction", "payload_escape_fraction"),
                       ("neighbours_in_reach", "neighbours_in_reach"),
                       ("direct_kill_probability", "direct_kill_probability")):
        m = re.search(rf"\b{field}:\s*([0-9.]+)", fixture.group(1))
        assert m, f"adc.rs::heterogeneous has no {field}"
        assert float(m.group(1)) == d["P10"]["constants"][key], (
            f"P10 uses {d['P10']['constants'][key]} for {field}; adc.rs says "
            f"{m.group(1)}")

    adoptive = (core / "adoptive.rs").read_text()
    m = re.search(r"fn solid_tumour\(\) -> Self \{(.*?)\n    \}", adoptive, re.S)
    assert m, "adoptive.rs no longer defines solid_tumour()"
    t = re.search(r"\btrafficking:\s*([0-9.]+)", m.group(1))
    assert float(t.group(1)) == d["P11"]["trafficking_barrier"], (
        f"P11 uses {d['P11']['trafficking_barrier']} for trafficking; "
        f"adoptive.rs says {t.group(1)}")

    # P13 must be READ, not a literal True.
    abl = (core / "ablation.rs").read_text()
    body = re.search(r"pub fn margin_survival_fraction\(.*?\n\}", abl, re.S)
    assert body, "ablation.rs no longer defines margin_survival_fraction"
    assert d["P13"]["returns_one_minus_covered"] == ("1.0 - covered" in body.group(0))
    # P13 registers that the RETURN VALUE does not vary with energy, not that
    # the body never reads it -- the body reads temperature, duration and field
    # strength to test the threshold, and an earlier version of this claim
    # denied that.
    assert d["P13"]["return_value_varies_with_energy"] is False, (
        "margin_survival_fraction's return now varies with an energy "
        "quantity, so P13's geometry-not-dose claim must be rewritten rather "
        "than left standing")
    assert d["P13"]["reads_energy_only_to_test_the_threshold"] is True, (
        "the threshold test is gone from margin_survival_fraction, so the "
        "'above threshold' scope P13 registers no longer applies")


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
        # BOUNDED AT THE NEXT FIELD OR THE NEXT HEADING. P13 is last, so its
        # slice ran on into the Honesty clause and a threshold reading
        # literally "no." passed on 77 borrowed words.
        thr = re.split(r"\n- \*|\n#{2,3} ", body.split(
            "*Falsification threshold:*", 1)[1])[0]
        assert len(thr.split()) >= 12, (
            f"{p}'s threshold is too short to name an outcome: {thr!r}")
        assert any(w in thr.lower() for w in
                   ("or ", "exceeds", "below", "above", "fails", "flat",
                    "rising", "at or")), (
            f"{p}'s threshold names no comparison, so no result can be scored "
            f"against it: {thr!r}")


def test_the_registered_numbers_are_the_derived_ones(d):
    """P9 is the one quantitative prediction, so its figures must be live."""
    text = PREREG.read_text()
    p9 = d["P9"]
    lo_r, hi_r = p9["low_over_high_dose_ratio"]
    # INSIDE P9's OWN BLOCK. Checking the document as a whole let a reviewer
    # gut P9's model output and satisfy every fragment from an unrelated
    # appendix sentence about a thermostat and a coffee machine.
    m = re.search(r"^\*\*P9\. (.+?)\*\*\n(.*?)(?=^\*\*P|^### )",
                  text, re.M | re.S)
    assert m, "P9 is not registered"
    block = m.group(2)
    for frag in (f"{p9['published_band'][0]} to {p9['published_band'][1]} band",
                 f"{p9['boost_window'][0]} to {p9['boost_window'][1]}",
                 f"{lo_r} to {hi_r}",
                 f"{p9['dose_range_gy'][0]} to {p9['dose_range_gy'][1]} Gy"):
        assert frag in block, (
            f"P9's own entry does not quote {frag!r}; it reads: {block[:200]}")
    # and the falsifier must score the BAND, not only the ratio: a measurement
    # violating the band at every dose can satisfy a ratio-only threshold.
    thr = block.split("*Falsification threshold:*", 1)[1]
    assert "outside" in thr and f"{p9['published_band'][0]}" in thr, (
        "P9's threshold scores only the ratio, so the arm's whole calibration "
        "could be refuted while the prediction stands")
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
    # STRICT, and steep. `vals == sorted(vals, reverse=True)` accepts a
    # CONSTANT -- which is the exact vacuity this prediction exists to retract,
    # reproduced in the guard written to replace it.
    assert all(a > b for a, b in zip(vals, vals[1:])), (
        f"the bystander reach no longer strictly FALLS as antigen is lost "
        f"({reach}); P10 is registered on that decline and a flat series is "
        "the shape the retracted guard could not distinguish from a real one")
    assert vals[0] > vals[-1] * 5.0, (
        f"the decline is shallow ({vals[0]:.4f} to {vals[-1]:.4f}); P10 claims "
        "the mechanism is starved, not mildly reduced")
    # and it must be a SHARE, which the retracted version was not: it
    # published 216%.
    for k, v in reach.items():
        assert 0.0 <= v <= 1.0, (
            f"the reach at {k} antigen-positive is {v}, which is not a share "
            "of a pool and cannot be produced by the registered experiment")
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
    # THE CLASSIFIER READS THE HEADLINE ONLY, so the widening this work exists
    # to demonstrate is decided by five sentences the author wrote. Pin each
    # new prediction to the ARM it is about, so renaming one cannot silently
    # reverse the measure.
    text = PREREG.read_text()
    for p, must in (("P9", "PARP"), ("P10", "ADC"), ("P11", "adoptive"),
                    ("P12", "oncolytic"), ("P13", "ablation")):
        m = re.search(rf"^\*\*{p}\. (.+?)\*\*", text, re.M)
        assert m and must.lower() in m.group(1).lower(), (
            f"{p}'s headline no longer names {must}; the audit classifies on "
            f"that sentence alone, so this is what decides the 8-of-13 measure")


def test_the_page_refuses_to_present_placeholders_as_predictions():
    md = MD.read_text()
    for frag in ("Directional, every one of them",
                 "uncalibrated placeholders",
                 "A magnitude here is not a prediction",
                 "Not clinical guidance"):
        assert frag in md, f"the page no longer says: {frag}"
