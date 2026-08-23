"""Guards for the oxygen-effect functional form (#726).

`analysis/oxygen-form-check.md` found that the engine's linear O2 dependence is
not the oxygen enhancement ratio, and that the two disagree most below 10 mmHg
-- the regime carrying section 7.1 and prediction P4, the leg this project
calls its weakest. The measured form now ships alongside the linear one.

THE PROPERTY THAT PROTECTS EVERYTHING ELSE is that selecting the form cannot
move an unconfigured run. Both forms return exactly 1.0 when the O2 dependence
is off, so the production matrix SHA is untouched; the Rust side asserts equal
bit patterns and these guards keep the Python side honest about it.

THE PROPERTY THAT MAKES IT WORTH HAVING is that it is not a renaming. If the
two forms ever agree across the hypoxic range, the layer has stopped being a
different answer and should be retired rather than kept.

OFFLINE: pure stdlib plus the committed artifacts and the Rust source as text.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/calibration/oer-form-validation.json"
MD = REPO / "analysis/calibration/oer-form-validation.md"
OXY = REPO / "simulations/ferroptosis-core/src/oxygen.rs"
STATUS = REPO / "simulations/calibration/CALIBRATION_STATUS.md"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_published_anchors_hold(d):
    """m(0)=1, half the rise at 3 mmHg, asymptote 3. The relation has no free
    parameter, so any of these failing means the formula changed."""
    assert d["anchors"], "no anchors checked"
    for a in d["anchors"]:
        assert a["ok"], f"{a['name']}: {a['computed']} vs {a['published']}"


def test_the_validation_reads_the_rust_rather_than_reimplementing_it(d):
    """A validation that only reimplements the formula checks itself.

    AND SO DOES A GUARD THAT READS ONLY THE ARTIFACT. The first version of
    this test asserted the stored `formula_matches_published` flag, so
    changing the Rust formula without regenerating left it green -- the
    generator-guard defect `tests/test_generator_guards_pin_source.py` exists
    for, found here by mutation rather than by the suite. The formula is
    re-parsed from source at test time now.
    """
    src = OXY.read_text()
    body = re.search(
        r"pub fn oxygen_enhancement_ratio\(p_mmhg: f64\) -> f64 \{(.*?)\n\}",
        src, re.S)
    assert body, "oxygen_enhancement_ratio not found in oxygen.rs"
    live = " ".join(body.group(1).split())
    assert live == d["expected_formula"], (
        f"the Rust formula is {live!r}, expected {d['expected_formula']!r}")
    assert d["formula_matches_published"], (
        "the committed artifact records a formula mismatch")
    assert d["rust"]["formula"] == live, (
        "the artifact's stored formula and the live source disagree, so the "
        "validation was not regenerated after the Rust changed")
    m = re.search(r"pub const OER_REFERENCE_PO2_MMHG: f64 = ([\d.]+);", src)
    assert m and float(m.group(1)) == d["rust"]["reference_po2_mmhg"], (
        "the committed artifact and the Rust source disagree on the reference "
        "pO2, so the drift guard has drifted")
    assert d["rust"]["has_oer_exo_factor"]
    assert d["rust"]["has_relative_efficacy"]


def test_the_forms_agree_at_full_supply(d):
    """Normalising at full supply is what makes an A/B measure SHAPE.

    If they disagreed at the endpoint too, switching forms would move the
    normoxic baseline as well, and a difference in kill could not be
    attributed to the hypoxic response.
    """
    assert d["endpoints_agree"], (
        "the two forms disagree at full O2 supply, so an A/B between them "
        "confounds shape with scale")
    last = d["rows"][-1]
    assert last["pO2_mmhg"] == d["rust"]["reference_po2_mmhg"]


def test_it_is_a_different_answer_not_a_renaming(d):
    """#726 originally claimed the engine ALREADY contained the OER. It did
    not, and if it ever does, this layer stops earning its place."""
    assert d["anoxic_linear"] == 0.0, (
        "the linear form no longer sends an anoxic cell to zero, so the "
        "motivating disagreement has changed and this analysis needs redoing")
    assert 0.3 < d["anoxic_oer"] < 0.4
    assert d["worst_gap"] > 0.5, (
        f"the largest gap between the forms is only {d['worst_gap']}, so they "
        "have converged and the second form is no longer a different answer")
    assert d["worst_gap_pO2"] <= 10, (
        "the forms now diverge most OUTSIDE the sub-10 mmHg range, which is "
        "the range this layer was justified by")


def test_the_report_scopes_itself_to_the_form(d):
    """A form validated against data sits beside a fraction that is not, and
    the two are easy to conflate into a calibrated layer."""
    md = " ".join(MD.read_text().split())
    assert "does NOT calibrate the Type I/Type II fraction" in md
    assert "drift guard rather than a fit" in md
    assert "says nothing about whether the exogenous-ROS mechanism is right" in md


def test_the_layer_freeze_row_exists_and_states_both_halves():
    """CONTRIBUTING.md requires a CALIBRATION_STATUS row naming the target."""
    s = " ".join(STATUS.read_text().split())
    assert "Oxygen-effect FUNCTIONAL FORM (OER)" in s, (
        "the new axis has no CALIBRATION_STATUS row, which the layer-freeze "
        "policy requires in the same PR")
    row = s[s.index("Oxygen-effect FUNCTIONAL FORM (OER)"):]
    row = row[:row.index("Used in any reported number") + 60]
    assert "the Type II fraction it scales remains uncalibrated" in row, (
        "the row does not separate the calibrated form from the uncalibrated "
        "fraction, so the layer reads as fully calibrated")
    assert "Used in any reported number: **N**" in row
    assert "0.0% (linear) -> 7.9% (OER)" in row, (
        "the row asserts the change matters without the measured consequence")


def test_the_rust_guards_the_inert_case_with_bit_equality():
    """A float comparison with a tolerance would let the flag move the
    production matrix by an amount smaller than the tolerance."""
    src = OXY.read_text()
    main = (REPO / "simulations/sim-tme-3d/src/main.rs").read_text()
    assert "fn oer_exo_factor_is_bit_identical_to_the_linear_form_at_zero_dependence" in src
    assert "to_bits()" in src.split("oer_exo_factor_is_bit_identical")[1][:900], (
        "the identity test does not compare bit patterns")
    fn = "the_oer_form_raises_hypoxic_kill_and_is_inert_when_the_dependence_is_off"
    assert fn in main, "sim-tme-3d has no A/B for the form"
    body = main.split(fn)[1][:3000]
    assert "to_bits()" in body, (
        "the spatial inert-case check does not compare bit patterns")


def test_the_layer_has_a_caller():
    """A correct layer with no caller is not progress -- a defect this repo has
    recorded three times. The form must be selectable from the binary."""
    main = (REPO / "simulations/sim-tme-3d/src/main.rs").read_text()
    assert "sdt_o2_use_oer" in main
    assert "oer_exo_factor(" in main, (
        "sim-tme-3d imports the flag but never calls the function")
    assert 'name: "sdt-oer"' in main, "no preset exercises the form"
