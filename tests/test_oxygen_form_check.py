"""Guards for the O2 functional-form comparison (#726).

THE CLAIM
---------
#726 rested on the sentence that `oxygen::o2_dependent_exo_factor` "is the
oxygen enhancement ratio under another name". It is not: the engine is linear in
O2 and the Alper-Howard-Flanders OER is a saturating hyperbola. At 7.5 mmHg the
engine gives 0.125 where the OER gives 0.810, and at anoxia 0.000 against 0.333.

That withdraws the issue's premise and strengthens its conclusion: the engine's
least-certain parameter would genuinely improve by adopting the measured form,
and the OER supplies the calibration target the layer-freeze policy demands.

WHAT MUST NOT DRIFT
-------------------
1. THE ENGINE ARM MUST TRACK THE ENGINE. The comparison is only meaningful if
   one side really is the shipped function. A guard checks the reproduction
   against oxygen.rs, so a change there fails here rather than silently making
   this analysis describe a function nobody runs.

2. THE TWO FORMS MUST STAY DIFFERENT FUNCTIONS. If they ever agree everywhere,
   that is a bug in one of them, not a finding.

3. IT MUST NOT SILENTLY BECOME AN ENGINE CHANGE. Adopting the hyperbola moves
   every result that depends on it, including committed byte-identity gates.
   This analysis quantifies; changing is a separate reviewed act.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "oxygen_form_check.py"
MD = REPO_ROOT / "analysis" / "oxygen-form-check.md"
JSON_OUT = REPO_ROOT / "analysis" / "oxygen-form-check.json"
OXYGEN_RS = REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "oxygen.rs"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("ofc", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_engine_arm_reproduces_the_shipped_function():
    """A comparison against a function nobody runs proves nothing."""
    rs = OXYGEN_RS.read_text()
    assert "(1.0 - d + d * s).clamp(0.0, 1.0)" in rs, (
        "oxygen.rs no longer computes the linear form this analysis reproduces; "
        "the comparison now describes a function the engine does not run")
    m = _mod()
    # spot values, computed from the Rust expression by hand
    assert abs(m.engine_linear(0.0, 1.0) - 0.0) < 1e-12
    assert abs(m.engine_linear(1.0, 1.0) - 1.0) < 1e-12
    assert abs(m.engine_linear(0.5, 0.4) - (1 - 0.4 + 0.4 * 0.5)) < 1e-12
    # and the clamps must be there, since the Rust clamps
    assert m.engine_linear(2.0, 1.0) == 1.0
    assert m.engine_linear(-1.0, 1.0) == 0.0


def test_the_two_forms_are_genuinely_different():
    """Agreement everywhere would be a bug in one of them."""
    d = _doc()
    assert d["max_abs_diff"] > 0.05, (
        f"the forms differ by at most {d['max_abs_diff']:.4f}; either the "
        "hyperbola has been flattened or the linear arm is no longer linear")
    src = SCRIPT.read_text()
    assert 'if d["max_abs_diff"] == 0.0:' in src, (
        "the generator no longer refuses to publish when the forms agree")


def test_the_hyperbola_is_the_published_one():
    """K and m_max are canonical values, not free parameters fitted here."""
    m = _mod()
    assert 1.0 <= m.OER_K_MMHG <= 5.0, (
        f"K = {m.OER_K_MMHG} mmHg is outside the reported Alper-Howard-Flanders "
        "range; a fitted K would make this a model of this model")
    assert 2.0 <= m.OER_M_MAX <= 3.5, (
        f"m_max = {m.OER_M_MAX} is outside the range reported for sparsely "
        "ionising radiation")
    # the defining property: half-maximal at K
    lo, hi = m.oer_hyperbola(0.0), m.oer_hyperbola(1e6)
    mid = m.oer_hyperbola(m.OER_K_MMHG)
    assert abs(mid - (lo + hi) / 2) < 0.02, (
        "the hyperbola is not half-way between anoxic and oxic at pO2 = K, so "
        "it is not the Alper-Howard-Flanders form")


def test_the_disagreement_is_worst_where_the_contested_claim_lives():
    """The finding is not that they differ but WHERE."""
    d, md = _doc(), MD.read_text()
    assert d["max_abs_diff_at"] <= 20.0, (
        f"maximum disagreement is at {d['max_abs_diff_at']} mmHg, outside the "
        "hypoxic range; the argument that this bears on P4 would not hold")
    assert f"{d['max_abs_diff']:.3f}" in md, (
        "the report does not state its own maximum disagreement")
    assert "P4" in md, (
        "the report no longer names the prediction this bears on, which is what "
        "makes it a finding rather than a curiosity")


def test_it_does_not_change_the_engine():
    """Adopting the hyperbola moves committed byte-identity gates."""
    md, src = MD.read_text(), SCRIPT.read_text()
    assert "does not change the engine" in md, (
        "the report no longer disclaims modifying the engine")
    # concretely: this script must never write into simulations/
    assert "simulations" not in src.split("OUT_JSON")[1], (
        "the generator writes into simulations/; this analysis quantifies a "
        "difference and must not apply one")


def test_the_withdrawn_premise_is_recorded():
    """The issue's own sentence was wrong and the correction must stay visible."""
    md = MD.read_text()
    assert "**It is not.**" in md or "It is not" in md, (
        "the report no longer states that the OER-equivalence claim is false")
    assert "under another name" in md, (
        "the withdrawn wording is no longer quoted, so a reader cannot tell "
        "what was corrected")


def test_the_mapping_assumption_is_disclosed():
    """pO2 to o2_supply is an assumption this comparison makes."""
    md = MD.read_text()
    assert "assumption this comparison makes" in md, (
        "the report no longer discloses that the pO2-to-supply mapping is its "
        "own assumption, which sets the size of the gap")
