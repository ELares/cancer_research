"""The step->wall-clock binding, pinned to the source it is read from.

`PREREGISTRATION.md` states every prediction in hours or days and the model
runs in steps, so nothing time-stated can be scored until the conversion is
written down. `parameter_provenance.md` now records it. These guards keep that
record tied to the code rather than to somebody's memory of the code.

THE FAILURE MODE THIS PREVENTS is the one the issue itself fell into: counting
quantities that measure different things as competing bindings. A PDE
integrator timestep, a defence half-life in days and a per-step PK sample rate
are three different clocks, and pooling them manufactures a disagreement.

OFFLINE: reads only committed source and artifacts.
"""
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOC = REPO / "simulations/calibration/parameter_provenance.md"
AUDIT = REPO / "analysis/engine-time-audit.json"
CORE = REPO / "simulations/ferroptosis-core/src"


def _doc():
    return " ".join(DOC.read_text().split())


def test_the_binding_matches_the_only_one_the_code_carries():
    """Derived from the audit, not typed into the prose."""
    d = json.loads(AUDIT.read_text())
    prod = d["production_minutes_per_step"]
    assert prod == [1.0], (
        f"the engine now carries {prod} minutes/step reaching production, so "
        "the reconciliation in parameter_provenance.md is stale")
    wall = [b for b in d["step_bindings"] if b["kind"] == "wall-clock"]
    assert wall and {b["module"] for b in wall} == {"tumor_pk.rs"}, (
        "a second module now prices a step in wall-clock time, so 'exactly "
        "one binding' is false")
    doc = _doc()
    assert "180 steps = 3.0 hours" in doc
    assert "**1 minute**" in doc


def test_the_two_axes_are_kept_apart():
    """Comparing the inner loop against the recovery axis is the asymmetric
    comparison `engine_time_audit.py` retracted."""
    doc = _doc()
    assert "two separate time axes" in doc
    assert "GSH 1, GPX4 3, NRF2 5, FSP1 7" in doc
    # And those half-lives must still be what cell.rs says.
    src = (CORE / "cell.rs").read_text()
    for field, days in (("gsh", 1.0), ("gpx4", 3.0), ("nrf2", 5.0), ("fsp1", 7.0)):
        m = re.search(rf"{field}_half_recovery_days: ([\d.]+)", src)
        assert m and float(m.group(1)) == days, (
            f"{field}_half_recovery_days is {m.group(1) if m else '?'} in "
            f"cell.rs, not {days} as the doc records")


def test_the_integrator_timestep_is_not_counted_as_a_binding():
    doc = _doc()
    assert "Not a step binding" in doc
    assert "CFL-constrained integrator timestep" in doc
    d = json.loads(AUDIT.read_text())
    tw = [b for b in d["step_bindings"] if b["module"] == "trigger_wave.rs"]
    assert tw and all(b["kind"] == "solver-timestep" for b in tw), (
        "trigger_wave's dt_min is no longer classified as a solver timestep")


def test_the_cost_of_the_binding_is_recorded():
    """A binding that makes an existing anchor unreachable has to say so."""
    doc = _doc()
    assert "cannot represent that accumulation at all" in doc, (
        "the doc adopts 1 min/step without recording that the 48-72h MUFA "
        "window does not fit inside a 3.0-hour run")
    assert "48-72 hour" in doc
    # The claim about params.rs must still hold.
    src = (CORE / "params.rs").read_text()
    assert "48-72h" in src or "48 to 72" in src


def test_it_does_not_attribute_a_figure_to_the_manuscript_that_is_not_there():
    """#727 credits the manuscript with ~16 min/step. It states no duration."""
    doc = _doc()
    assert "states no per-step duration" in doc
    txt = " ".join((REPO / "article/drafts/v1.md").read_text().split())
    assert "180 steps within a single treatment window" in txt
    assert "16 min/step" not in txt and "16 minutes per step" not in txt, (
        "the manuscript now states a per-step duration, so the doc's claim "
        "that it does not is stale")


def test_the_production_matrix_is_marked_dimensionless():
    """The narrowness is the point: the default matrix never enters the PK
    solver, so calling it 3.0 hours would be a claim about a path it does not
    take."""
    doc = _doc()
    assert "no wall-clock reading at" in doc and "dimensionless" in doc
    main = (REPO / "simulations/sim-tme-3d/src/main.rs").read_text()
    assert "DoseSchedule::Constant" in main


def test_the_status_line_does_not_overclaim():
    doc = _doc()
    assert "Assumed" in doc
    assert "No figure in this book is computed from a step duration" in doc
    for over in ("calibrated against", "measured at 1 minute", "validated"):
        assert over not in doc.split("### Status")[1]
