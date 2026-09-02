"""Guards for the chemotherapy rim-bias decomposition (#844).

A decomposition is only worth having if its null case is genuinely null, so the
first guard is on the control: with neither term active the kill must be flat
with depth, or the other three columns are decomposing the zoning rather than
the biology. The second is on the discriminator: the cycle term must behave
differently by agent class, or it is measuring the drug field twice under
another name.
"""
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_chemo_decomposition.py"
MD = REPO / "analysis" / "calibration" / "chemo-decomposition-validation.md"
JSON_ = REPO / "analysis" / "calibration" / "chemo-decomposition-validation.json"
SWEEP = REPO / "analysis" / "calibration" / "chemo_decomposition_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"


def _d():
    return json.loads(JSON_.read_text())


def test_the_committed_page_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON_.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "the page is stale"
        assert JSON_.read_text() == before_json, "the JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON_.write_text(before_json)


def test_the_null_case_is_actually_null():
    """The load-bearing control.

    With neither delivery gradient nor cycle coupling the kill must be flat
    with depth. If it were not, the other three columns would be decomposing
    however the zones happen to be drawn rather than the two named terms.
    """
    d = _d()
    assert d["control_is_flat"]
    for c in d["classes"]:
        assert c["control_flat"] < d["flat_tolerance"], (
            f"{c['class']} already has a {c['control_flat']:.1%} depth gradient "
            "with both terms OFF, so the decomposition has a term it did not "
            "name")
        assert c["neither"]["core"] > 0.05, (
            f"{c['class']}'s null case barely kills anything in the core, so "
            "the ratios below divide by nearly nothing")


def test_the_cycle_term_is_class_dependent_or_it_is_not_the_cell_cycle():
    """The discriminator.

    A phase-nonspecific alkylator damages DNA whatever the cell is doing, so a
    quiescent core should cost it comparatively little. If every class lost the
    same share, the 'cycle' arm would be measuring the drug field twice.
    """
    d = _d()
    assert d["cycle_term_is_class_dependent"]
    assert d["cycle_cost_ratio"] > 1.3, (
        f"the phase-nonspecific agent keeps only {d['cycle_cost_ratio']:.2f}x "
        "as much core kill as the phase-specific ones under quiescence, which "
        "is too little separation to call it a cell-cycle effect")
    nonspec = next(c for c in d["classes"] if c["class"] == "PhaseNonspecific")
    assert nonspec["core_retained_cycle_only"] > 0.5
    for c in d["classes"]:
        if c["class"] != "PhaseNonspecific":
            assert c["core_retained_cycle_only"] < 0.5


def test_delivery_dominates_and_the_page_says_it_depends_on_parameters():
    d = _d()
    assert d["delivery_dominates_everywhere"]
    for c in d["classes"]:
        assert c["core_retained_delivery_only"] < 0.05
    md = MD.read_text()
    assert "Which term dominates depends on the parameters" in md, (
        "the page no longer says the dominance is parameter-dependent, which "
        "is the caveat that keeps it from reading as a general result")
    assert "longer penetration length would shift the balance" in md


def test_all_four_cells_of_the_two_by_two_are_present_per_class():
    d = _d()
    assert len(d["classes"]) == 3
    for c in d["classes"]:
        for key in ("neither", "cycle_only", "delivery_only", "both"):
            assert key in c
            for zone in ("overall", "rim", "core"):
                assert 0.0 <= c[key][zone] <= 1.0


def test_the_sweep_is_the_binarys_own_output():
    txt = SWEEP.read_text()
    for l in txt.splitlines():
        assert l.startswith("CHEMO_DECOMP"), f"stray line: {l[:60]!r}"
    src = BIN.read_text()
    assert '"--chemo-decomposition-sweep"' in src
    assert "fn run_chemo_decomposition_sweep" in src
    # The quiescent core has to come from somewhere: without the spheroid
    # layer every cell has the same phenotype and the cycle arm reads nothing.
    assert "spheroid: Some(SpheroidConfig::literature())" in src
    assert "the cycle arm has nothing to read" in src


def test_the_page_states_what_it_does_not_establish():
    md = MD.read_text()
    for phrase in ("NO TARGET", "cycle is read from PHENOTYPE, not simulated",
                   "conventional"):
        assert phrase in md, f"the page no longer states its limit: {phrase!r}"
    src = (REPO / "simulations" / "ferroptosis-core" / "src" / "chemo.rs").read_text()
    assert "conventional rather than measured" in src, (
        "chemo.rs no longer records that its phase sensitivities are "
        "conventional, which is what the page's last caveat rests on")
