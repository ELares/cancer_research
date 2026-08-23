"""The step-to-wall-clock readings, pinned to the sources they are read from.

Prediction P3 is stated in days and the model runs in steps, so scoring it
needs a conversion. `parameter_provenance.md` records that the engine does not
have one conversion but TWO, 16x apart, in different binaries.

WHAT THE FIRST VERSION OF THIS FILE GOT WRONG, kept here because the shape
recurs. It asserted there was exactly one binding and guarded that with a grep
for the literal string "16 min/step" -- a quotient no source or document would
ever contain -- so it was inert against precisely the reading it named. The
audit behind it globbed the library only, so the claim was measured over a
scope narrower than the claim itself.

These guards therefore check the MEASUREMENT, not the prose: the audit must
still scan the binaries, must still find both readings, and the doc must still
name both.

OFFLINE: reads only committed source and artifacts.
"""
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOC = REPO / "simulations/calibration/parameter_provenance.md"
AUDIT = REPO / "analysis/engine-time-audit.json"
CORE = REPO / "simulations/ferroptosis-core/src"
SIMS = REPO / "simulations"


def _doc():
    return " ".join(DOC.read_text().split())


def _audit():
    return json.loads(AUDIT.read_text())


def test_the_audit_scans_the_binaries_not_only_the_library():
    """The scope error that made the wrong claim measurable-and-green."""
    src = (REPO / "scripts/engine_time_audit.py").read_text()
    assert "_rust_sources" in src
    import sys
    sys.path.insert(0, str(REPO / "scripts"))
    import engine_time_audit as m

    files = m._rust_sources()
    scanned = {p.parent.parent.name for p in files}
    # EVERY sim crate, not "at least one": narrowing the glob to a single
    # binary re-commits the original scope error for the other eleven and an
    # `any()` check cannot see it.
    on_disk = {d.name for d in (REPO / "simulations").glob("sim-*")
               if (d / "src").is_dir()}
    assert on_disk <= scanned, (
        f"the audit does not scan {sorted(on_disk - scanned)}, so a binding "
        "living there is invisible to every claim it makes")
    assert "ferroptosis-core" in scanned
    # And the file list must be deterministic: the artifact is byte-compared
    # in CI on two filesystems.
    assert files == m._rust_sources()
    assert files == sorted(files, key=lambda q: (q.parent.parent.name, q.name))


def test_both_competing_readings_are_still_found():
    """Derived from the audit, not typed into the prose."""
    d = _audit()
    explicit = {b["minutes_per_step"] for b in d["step_bindings"]
                if b["kind"] == "wall-clock"}
    implied = {w["minutes_per_step"] for w in d.get("implied_windows", [])}
    assert explicit == {1.0}, f"explicit wall-clock readings moved: {explicit}"
    assert implied == {16.0}, (
        f"implied-window readings moved: {implied}. If a binary changed its "
        "stated validity window or its step count, the reconciliation is stale")
    doc = _doc()
    assert "**1.0**" in doc and "**16.0**" in doc
    assert "16x apart" in doc


def test_the_implied_window_is_still_what_sim_tme_says():
    """A doc-comment scope claim plus a loop length. Checked at the source."""
    main = (SIMS / "sim-tme/src/main.rs").read_text()
    assert re.search(r"const N_STEPS:\s*\w+\s*=\s*180", main), (
        "sim-tme's loop length changed, so 16 min/step no longer follows")
    assert "resident T cell phase (0-48h)" in main
    w = _audit()["implied_windows"][0]
    assert w["binary"] == "sim-tme" and w["n_steps"] == 180
    assert w["minutes_per_step"] == round(w["window_hours"] * 60 / w["n_steps"], 4)


def test_the_doc_does_not_adopt_either_reading_as_correct():
    """Neither is measured, so picking one would be a decision dressed as a
    finding. The doc says which applies WHERE instead."""
    doc = _doc()
    assert "Neither is adopted as correct" in doc
    assert "Anything combining the two is unconvertible until one is measured" in doc, (
        "the doc no longer states that a quantity spanning both readings "
        "cannot be converted, which is the practical consequence of adopting "
        "neither")
    for over in ("the step duration is 1 minute", "we adopt 1 min/step",
                 "the binding is 1 min/step"):
        assert over not in doc.lower()


def test_the_classification_is_not_keyed_on_spelling():
    """`kind` used to be `solver-timestep if named-field else wall-clock`, so
    any second binding spelled `dt_<unit>` was excluded by construction."""
    src = (REPO / "scripts/engine_time_audit.py").read_text()
    assert "SOLVER_MODULES" in src
    # Narrowly: the KIND assignment must not key on the spelling. The same
    # test elsewhere legitimately decides how to READ a magnitude.
    kind_line = [l for l in src.splitlines() if l.strip().startswith("kind = ")]
    assert kind_line and "SOLVER_MODULES" in kind_line[0], (
        f"kind is assigned as {kind_line!r}, not from SOLVER_MODULES")
    assert not any('how ==' in l for l in kind_line)
    d = _audit()
    tw = [b for b in d["step_bindings"] if b["module"] == "trigger_wave.rs"]
    assert tw and all(b["kind"] == "solver-timestep" for b in tw)
    assert all(abs(b["minutes_per_step"] - 0.02) < 1e-9 for b in tw), (
        "trigger_wave's dt_min value moved; `_default_for` used to ignore its "
        "stem argument and could report another module's literal here")


def test_the_two_axes_and_the_shared_inner_loop_are_both_stated():
    doc = _doc()
    assert "GSH 1, GPX4 3, NRF2 5, FSP1 7" in doc
    assert "P3 is scored across BOTH axes" in doc, (
        "the doc claims P3 is scored on the recovery axis alone; sim-window "
        "sets initial state and then runs the same 180-step loop")
    src = (CORE / "cell.rs").read_text()
    for field, days in (("gsh", 1.0), ("gpx4", 3.0), ("nrf2", 5.0), ("fsp1", 7.0)):
        m = re.search(rf"{field}_half_recovery_days: ([\d.]+)", src)
        assert m and float(m.group(1)) == days


def test_the_costs_are_recorded_on_both_sides():
    """A reconciliation that only lists what the value it prefers makes
    unreachable is an argument, not an accounting."""
    doc = _doc()
    assert "cannot represent that accumulation at all" in doc
    assert "borderline reachable" in doc, (
        "the doc records what 1 min/step silences without recording that the "
        "other reading relieves it, which is an argument for the other side")
    assert "clinically ordinary" in doc


def test_the_siblings_do_not_contradict_it():
    """Two pages said this mapping did not exist and the decision was open."""
    cal = " ".join((SIMS / "calibration/CALIBRATION_STATUS.md").read_text().split())
    # Checked as AGREEMENT, not as a banned phrase: rewriting the line to
    # claim a single settled mapping is a flat contradiction and passed a
    # substring ban.
    assert "a step-to-time mapping the suite lacks" not in cal
    assert "16x apart" in cal or "TWO of" in cal, (
        "CALIBRATION_STATUS does not record that two readings compete")
    for settled in ("single settled step-to-time", "the suite's step duration is",
                    "decided in #727a"):
        assert settled not in cal, (
            f"CALIBRATION_STATUS claims the mapping is settled ({settled!r}); "
            "the reconciliation adopts neither reading")
    audit_md = " ".join((REPO / "analysis/engine-time-audit.md").read_text().split())
    assert "It does not choose a step duration" in audit_md, (
        "engine-time-audit.md should still decline to CHOOSE -- the doc does "
        "not choose either -- but it must point at the reconciliation")
    assert "parameter_provenance" in audit_md


# ---------------------------------------------------------------------------
# Guards that bind the DOCUMENT to the measurement.
#
# A reviewer broke four independent doc mutations while every test stayed
# green: swapping which subsystem gets which conversion, changing 3.0 h to
# 99 h, falsifying both source attributions, and writing "16 min/step x 180 =
# 6 minutes". The guards derived the two numbers from JSON and then checked
# literals and fixed phrases, so nothing bound the table rows, the arithmetic,
# or the rule the document exists to state.
# ---------------------------------------------------------------------------

def _table_rows():
    """The reading table, parsed out of the doc."""
    rows = []
    for line in DOC.read_text().splitlines():
        if line.startswith("| ") and "**" in line and "|" in line[2:]:
            cells = [c.strip() for c in line.strip("|").split("|")]
            if len(cells) == 4:
                rows.append(cells)
    return rows


def test_every_row_of_the_reading_table_matches_the_measurement():
    """Row by row: the number, its arithmetic, and its cited sources."""
    d = _audit()
    rows = _table_rows()
    assert len(rows) == 2, f"expected two readings, parsed {len(rows)}"
    got = {}
    for label, source, mins, span in rows:
        v = float(mins.strip("* "))
        got[v] = (source, span)
        # 180 steps x minutes must equal the span the row claims.
        m = re.match(r"([\d.]+)\s*(h|min)", span.strip())
        assert m, f"unparseable span {span!r}"
        hours = float(m.group(1)) * (1 / 60 if m.group(2) == "min" else 1)
        assert abs(180 * v / 60 - hours) < 0.05, (
            f"row {label!r}: 180 steps x {v} min is "
            f"{180 * v / 60:.2f} h, not {hours} h")
    declared = {b["minutes_per_step"] for b in d["step_bindings"]
                if b["kind"] == "wall-clock"}
    implied = {w["minutes_per_step"] for w in d["implied_windows"]}
    assert set(got) == declared | implied, (
        f"the table lists {sorted(got)} and the audit measures "
        f"{sorted(declared | implied)}")
    # And each row must cite a file the audit actually attributes it to.
    for v, (source, _) in got.items():
        srcs = ([b["module"] for b in d["step_bindings"]
                 if b.get("minutes_per_step") == v]
                + [f"{w['binary']}/{w['module']}" for w in d["implied_windows"]
                   if w["minutes_per_step"] == v])
        assert any(x.split("/")[-1].split(".")[0] in source for x in srcs), (
            f"the {v} min/step row cites {source!r}, and the audit attributes "
            f"that reading to {srcs}")


def test_the_subsystem_rule_is_stated_and_not_a_per_binary_one():
    """The rule the document exists to give. An earlier draft split by BINARY
    and was wrong, because `sim-tme-3d` carries both readings."""
    doc = _doc()
    assert "not cleanly separable by binary" in doc
    assert "1 min/step** applies to the PK trajectory" in doc
    assert "16 min/step** applies to the immune cascade" in doc
    # And the reason must be checkable: sim-tme-3d really does state both.
    d = _audit()
    bins = {w["binary"] for w in d["implied_windows"]}
    assert "sim-tme-3d" in bins, (
        "sim-tme-3d no longer states an immune window, so the doc's reason for "
        "refusing a per-binary rule is stale")
    conv = {b["module"] for b in d["step_bindings"] if b["kind"] == "wall-clock"}
    assert any("tumor_pk" in c for c in conv)


def test_the_retraction_is_recorded_not_deleted():
    doc = _doc()
    assert "What this section previously claimed, and why it was wrong" in doc
    assert "a quotient no document would ever write" in doc
    assert "A claim is only as wide as the sweep behind it" in doc


def test_the_superseded_headline_is_gone_from_every_site():
    """It survived in four files after the correction, including the
    artifact's own opening line and this guard file's docstring."""
    sites = [REPO / "analysis/engine-time-audit.md",
             REPO / "scripts/engine_time_audit.py",
             REPO / "tests/test_engine_time_audit.py",
             REPO / "CLAUDE.md"]
    for f in sites:
        flat = " ".join(f.read_text().split())
        assert "EXACTLY ONE module prices a step" not in flat, (
            f"{f.name} still carries the retracted one-binding headline")
    md = " ".join((REPO / "analysis/engine-time-audit.md").read_text().split())
    assert "A second reading is IMPLIED and never declared" in md, (
        "the artifact reports only the declared binding, so its headline "
        "contradicts the reconciliation it points at")
