"""Guards for the PDT fluence-rate validation.

The failure mode this page invites is specific: it produces a number in
mW/cm² that looks like a recommendation, and that number rides on `phi_crit`,
which nothing in this repository measures. A caveat sentence is not enough --
the repository has shipped a placeholder as a finding before -- so the guards
require the page to keep MEASURING the dependence rather than mentioning it,
and require the interior-vs-edge distinction to stay live, since an edge
answer reported as an optimum states a scan bound as a result.
"""
import json
import math
import subprocess
import sys
from importlib import util
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_pdt_fluence_rate.py"
MD = REPO / "analysis" / "calibration" / "pdt-fluence-rate-validation.md"
JSON = REPO / "analysis" / "calibration" / "pdt-fluence-rate-validation.json"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "photosensitizer_pk.rs"
OXY = REPO / "simulations" / "ferroptosis-core" / "src" / "oxygen.rs"


def _d():
    return json.loads(JSON.read_text())


def _mod():
    spec = util.spec_from_file_location("vp", SCRIPT)
    m = util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_committed_page_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, (
            "pdt-fluence-rate-validation.md is stale; run the generator")
        assert JSON.read_text() == before_json, "the JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_every_reported_optimum_is_interior():
    """An edge answer means the model is monotonic over the scanned range and
    the binding limit is outside it. Reporting one as an optimum publishes a
    scan parameter as a finding, which is the defect the oncolytic and
    ablation sections each shipped once."""
    d = _d()
    assert d["all_optima_interior"]
    lo, hi = d["scan_range_mw_cm2"]
    for row in d["by_half_life"]:
        assert row["interior"], f"t_half {row['t_half_h']} sits on a scan edge"
        assert lo < row["optimal_mw_cm2"] < hi
        assert row["delivered"] > 0


def test_the_prediction_is_a_direction_and_it_is_monotone():
    """P21: the optimal rate falls as the sensitizer clears more slowly. If
    this ever stops being monotone the prediction has to be rewritten, not
    the threshold moved."""
    d = _d()
    assert d["optimum_falls_with_half_life"]
    rows = d["by_half_life"]
    assert len(rows) >= 5
    for a, b in zip(rows, rows[1:]):
        assert b["t_half_h"] > a["t_half_h"]
        assert b["optimal_mw_cm2"] < a["optimal_mw_cm2"]
    # And the effect must be large enough to be a claim rather than noise.
    assert rows[0]["optimal_mw_cm2"] / rows[-1]["optimal_mw_cm2"] > 5.0


def test_the_uncalibrated_knob_is_measured_and_not_merely_caveated():
    d = _d()
    assert d["phi_crit_span"] > 1.0
    assert d["optimum_span_over_phi_crit"] > 1.0
    exp = d["optimum_scaling_exponent"]
    assert 0.0 < exp < 1.0, (
        f"the optimum's dependence on phi_crit is {exp}; outside (0,1) the "
        "page's 'roughly the square root' reading is wrong")
    md = MD.read_text()
    assert "How much of the answer is the uncalibrated knob" in md
    assert "φ_crit" in md
    # The verdict must keep the two apart.
    assert d["verdict"] == "DIRECTION ANCHORED, MAGNITUDE UNCONSTRAINED"
    assert "MAGNITUDE UNCONSTRAINED" in md


def test_the_oxygen_term_is_the_engines_measured_shape_not_a_second_one():
    """Radiation and PDT are both oxygen-dependent arms in this engine. Two
    different O2 curves would let them drift apart for no reason other than
    which module was written first, which #726 already had to fix once."""
    src = RUST.read_text()
    assert "crate::oxygen::oer_relative_efficacy" in src, (
        "the fluence-rate yield no longer uses the AHF hyperbola")
    assert "o2_dependent_exo_factor" not in src.split("mod tests")[0], (
        "the linear O2 form has reappeared in the PDT production path")
    ref = _d()["oer_reference_po2_mmhg_from_oxygen_rs"]
    assert f"OER_REFERENCE_PO2_MMHG: f64 = {ref}" in OXY.read_text()


def test_the_oxygen_penalty_grows_with_rate_and_is_not_flat():
    rows = _d()["oxygen_blind_gap"]
    assert len(rows) >= 4
    for a, b in zip(rows, rows[1:]):
        assert b["rate_mw_cm2"] > a["rate_mw_cm2"]
        assert b["o2_fraction"] < a["o2_fraction"]
        assert b["yield_factor"] < a["yield_factor"]
    # If the penalty were negligible there would be no fluence-rate effect to
    # model at all, so the span has to be large enough to matter.
    assert rows[0]["yield_factor"] / rows[-1]["yield_factor"] > 1.1


def test_the_python_reimplementation_matches_the_rust_contract_at_its_fixed_point():
    """phi_crit is DEFINED as the rate where consumption matches resupply, so
    half the oxygen must remain there. A drift in either implementation turns
    the parameter into an arbitrary scale."""
    m = _mod()
    assert abs(m.o2_fraction(50.0, 50.0) - 0.5) < 1e-12
    assert m.o2_fraction(0.0, 50.0) == 1.0
    assert m.o2_fraction(10.0, 0.0) == 0.0
    assert "1.0 / (1.0 + phi / phi_crit_mw_cm2)" in RUST.read_text(), (
        "the Rust quasi-steady form has changed; the Python one must too")
    # The AHF hyperbola, checked by hand at a point neither file computes.
    assert abs(m.oxygen_enhancement_ratio(0.0) - 1.0) < 1e-12
    assert abs(m.oxygen_enhancement_ratio(3.0) - 2.0) < 1e-12


def test_the_page_states_what_it_does_not_claim():
    md = MD.read_text()
    for phrase in ("No fitted rate", "quasi-steady", "the binding constraint "
                   "is not in the model"):
        assert phrase in md, f"the page no longer states its limit: {phrase!r}"
    assert "PMID 16615136" in md, (
        "the fluence-rate anchor is no longer cited, so the direction is "
        "unsupported")


def test_the_integral_is_below_the_start_sample():
    """The defect the layer fixes, checked in Python too: under a decaying
    sensitizer, integrating the drug over the illumination must give LESS
    than sampling it once at the start -- and more so for a longer one."""
    m = _mod()
    for rate in (20.0, 200.0):
        hours = m.TOTAL_FLUENCE_J * 1000.0 / rate / 3600.0
        g = m.oer_relative_efficacy(m.o2_fraction(rate, 50.0), 40.0)
        naive = g * rate * 1.0 * hours * 3.6
        assert m.delivered(1.0, rate, 50.0) < naive
