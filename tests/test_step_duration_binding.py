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
    # And a withdrawn cost must stay withdrawn rather than quietly returning:
    # the persister anchor was cited as priced when `params.rs` says it is
    # explicitly NOT mapped.
    assert "NOT a third anchor" in doc
    assert "Turning a stated non-mapping into a pricing" in doc


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
    """The reading table, parsed out of the doc.

    Every four-cell row between the header and the blank line -- NOT only the
    bolded ones. Requiring `"**"` let an unbolded row be added invisibly, so
    the two-row assertion above it could not see a third reading appear.
    """
    lines = DOC.read_text().splitlines()
    try:
        i = next(n for n, l in enumerate(lines)
                 if l.startswith("| reading | source |"))
    except StopIteration:
        return []
    rows = []
    for line in lines[i + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) == 4:
            rows.append(cells)
    return rows


def test_every_row_of_the_reading_table_matches_the_measurement():
    """Row by row: the LABEL, the number, its arithmetic, and its sources.

    The first version bound only the numbers, and a reviewer inverted the
    document's central distinction with every test green -- calling the
    declared reading "implied" and vice versa -- then re-attributed the
    16 min/step reading to two unrelated binaries, because the source check
    reduced `sim-tme/README.md` to the token `README`, which is a substring of
    almost any path. It also accepted a third row and never read a line number
    the doc cites.
    """
    d = _audit()
    rows = _table_rows()
    assert len(rows) == 2, f"expected two readings, parsed {len(rows)}"

    declared = {b["minutes_per_step"] for b in d["step_bindings"]
                if b["kind"] == "wall-clock"}
    implied = {w["minutes_per_step"] for w in d["implied_windows"]}
    assert declared and implied and not (declared & implied)

    seen = {}
    for label, source, mins, span in rows:
        v = float(mins.strip("* "))
        # THE LABEL IS THE CLAIM. A reading the audit measures as declared may
        # not be presented as implied, or the whole finding inverts.
        lab = label.lower()
        if v in declared:
            assert "declar" in lab and "implied" not in lab, (
                f"{v} min/step is a DECLARED binding and the table labels it "
                f"{label!r}")
        elif v in implied:
            assert "implied" in lab and "declar" not in lab.replace("never declar", ""), (
                f"{v} min/step is only IMPLIED and the table labels it {label!r}")
        else:
            raise AssertionError(f"the table lists {v} min/step, which the "
                                 f"audit does not measure")
        m = re.match(r"([\d.]+)\s*(h|min)", span.strip())
        assert m, f"unparseable span {span!r}"
        hours = float(m.group(1)) * (1 / 60 if m.group(2) == "min" else 1)
        assert abs(180 * v / 60 - hours) < 0.05, (
            f"row {label!r}: 180 steps x {v} min is {180 * v / 60:.2f} h, "
            f"not {hours} h")
        seen[v] = source

    assert set(seen) == declared | implied

    # Sources compared on the FULL path, and every cited line number checked
    # against the audit. A token match let `sim-tme/README.md` be satisfied by
    # `sim-tumor-pk/src/main.rs`.
    for v, source in seen.items():
        paths = ([b["module"] for b in d["step_bindings"]
                  if b.get("minutes_per_step") == v and b["kind"] == "wall-clock"]
                 + [f"{w['binary']}/{w['module']}" for w in d["implied_windows"]
                    if w["minutes_per_step"] == v])
        assert paths
        assert any(x in source for x in paths), (
            f"the {v} min/step row cites {source!r}; the audit attributes that "
            f"reading to {sorted(set(paths))}")
        # EVERY cited file:line must exist, and the text quoted beside it must
        # actually be on that line. The previous version checked only that a
        # line number appeared in the audit's implied-window list -- so a row
        # could quote words the cited file does not contain, which is the "a
        # claim in file A about file B" defect this whole section retracts.
        cites = re.findall(r"`([^`]+?):(\d+)`\s*--\s*\"([^\"]+)\"", source)
        assert cites, f"the {v} min/step row cites no file:line with a quote"
        for path_s, lineno, quote in cites:
            f = REPO / "simulations" / path_s if not (REPO / path_s).exists() \
                else REPO / path_s
            for cand in (REPO / path_s, REPO / "simulations" / path_s,
                         REPO / "simulations/ferroptosis-core/src" / path_s):
                if cand.exists():
                    f = cand
                    break
            assert f.exists(), f"cited file {path_s} does not exist"
            src_lines = f.read_text().splitlines()
            n = int(lineno)
            assert 1 <= n <= len(src_lines), (
                f"{path_s} has {len(src_lines)} lines; the row cites {n}")
            line = src_lines[n - 1]
            assert quote in line, (
                f"the {v} min/step row quotes {quote!r} at {path_s}:{n}, "
                f"which reads {line.strip()!r}")


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


# ---------------------------------------------------------------------------
# THE PROSE, derived from source the way the table already is.
#
# Guarding the table and stopping there left eleven prose mutations passing a
# green suite: the 60-step delay read as "6 hours and ten minutes", the front
# speed as "~36 minute", the integrator gap as "5000x", the retraction as
# "The claim was TRUE", and an appended "in practice, adopt 1 min/step
# everywhere" -- the last two inverting the document's two central stances.
#
# Every number the prose states about the engine is recomputed here from the
# file that determines it. A hand-written figure in a document about
# hand-written figures being wrong is the defect this whole PR is retracting.
# ---------------------------------------------------------------------------

CORE = REPO / "simulations/ferroptosis-core/src"


def _num_near(doc: str, phrase: str, pattern: str) -> str:
    """The first match of `pattern` in the sentence containing `phrase`."""
    i = doc.index(phrase)
    window = doc[max(0, i - 200):i + 300]
    m = re.search(pattern, window)
    assert m, f"no {pattern!r} near {phrase!r}"
    return m.group(1)


def test_the_chapter_7_attribution_is_the_binary_that_produced_it():
    """Which binary produced the published immune numbers is load-bearing:
    it is why the implied reading matters rather than being a curiosity."""
    doc = _doc()
    d = _audit()
    bins = sorted({w["binary"] for w in d["implied_windows"]})
    m = re.search(r"`(sim-[\w-]+)` produced this book's", doc)
    assert m, "the doc no longer attributes the published immune numbers"
    named = m.group(1)
    assert named in bins, f"{named} states no immune window; audit has {bins}"
    # The 2D Chapter 7 numbers come from the 2D binary, not the 3D one.
    ch7 = (REPO / "article/drafts/v1.md").read_text()
    assert "104:1 in 2D" in ch7
    assert named == "sim-tme", (
        f"the doc credits {named}; the 104:1 figure is the 2D result and "
        "sim-tme is the 2D binary")


def test_the_immune_delay_figures_are_what_the_code_implies():
    """"reads as 16 hours under it and one hour under the other"."""
    doc = _doc()
    main = (SIMS / "sim-tme/src/main.rs").read_text()
    m = re.search(r"let immune_start_step\s*=\s*(\d+)", main)
    assert m, "sim-tme no longer declares an immune activation delay"
    steps = int(m.group(1))
    d = _audit()
    implied = max(w["minutes_per_step"] for w in d["implied_windows"])
    declared = max(b["minutes_per_step"] for b in d["step_bindings"]
                   if b["kind"] == "wall-clock")
    hours = steps * implied / 60
    assert f"{steps}-step immune activation delay" in doc, (
        f"the doc no longer cites the {steps}-step delay it computes from")
    assert f"reads as {int(hours)} hours" in doc, (
        f"{steps} steps x {implied} min is {hours:g} hours; the doc says "
        f"otherwise")
    other = int(steps * declared / 60)
    words = {1: "one", 2: "two", 3: "three"}
    assert (f"{other} hour under the other" in doc
            or f"{words.get(other, other)} hour under the other" in doc), (
        f"{steps} steps x {declared} min is {other} hour(s); the doc says "
        "otherwise")


def test_the_integrator_gap_is_recomputed():
    """"manufactures a 50x disagreement"."""
    doc = _doc()
    d = _audit()
    solver = min(b["minutes_per_step"] for b in d["step_bindings"]
                 if b["kind"] == "solver-timestep")
    declared = max(b["minutes_per_step"] for b in d["step_bindings"]
                   if b["kind"] == "wall-clock")
    assert f"{int(declared / solver)}x disagreement" in doc, (
        f"{declared} / {solver} is {declared / solver:g}x; the doc says "
        "otherwise")


def test_the_derived_front_speed_candidate_is_recomputed():
    """"5.52 um/min ... at the 20 um cell pitch ... a ~3.6 minute crossing"."""
    doc = _doc()
    tw = (CORE / "trigger_wave.rs").read_text()
    m = re.search(r"measured a baseline front speed of (\d+\.\d+)", tw)
    assert m, "trigger_wave no longer states a measured BASELINE front speed"
    speed = float(m.group(1))
    pitch = None
    for b in ("sim-tme", "sim-tme-3d"):
        mm = re.search(r"const CELL_SIZE_UM:\s*f64\s*=\s*(\d+(?:\.\d+)?)",
                       (SIMS / b / "src/main.rs").read_text())
        if mm:
            pitch = float(mm.group(1))
            break
    assert pitch, "no binary declares a cell pitch"
    crossing = pitch / speed
    assert f"{speed} um/min" in doc
    assert f"{pitch:g} um cell pitch" in doc
    assert f"~{crossing:.1f} minute cell crossing" in doc, (
        f"{pitch} um / {speed} um/min is {crossing:.2f} min; the doc says "
        "otherwise")
    # And it must be marked as DERIVED, not as a third reading in the source.
    assert "DERIVED, not stated anywhere in the engine's text" in doc
    assert "not a third entry in the table above" in doc
    assert "a third entry in the table above" not in doc.replace(
        "not a third entry in the table above", ""), (
        "the derived candidate is presented as a third reading in the table, "
        "which lists what the code SAYS")
    # The table has exactly the readings the audit measures, and this is not
    # one of them.
    assert len(_table_rows()) == 2


def test_the_manuscript_adjacency_claim_is_checked_against_the_manuscript():
    """"the bullet IMMEDIATELY above the 180-step sentence"."""
    doc = _doc()
    assert "bullet IMMEDIATELY above" in doc
    md = (REPO / "article/drafts/v1.md").read_text().splitlines()
    bullets = [n for n, l in enumerate(md) if l.startswith("**")
               and 1195 <= n <= 1225]
    window = next((n for n in bullets if "0-48 hour" in md[n]), None)
    steps = next((n for n in bullets if "180 steps within a single" in md[n]), None)
    assert window is not None and steps is not None, (
        "Section 8.4 no longer carries both bullets the doc cites")
    later = [n for n in bullets if n > window]
    assert later and later[0] == steps, (
        f"the 0-48h bullet (line {window + 1}) is not immediately above the "
        f"180-step bullet (line {steps + 1}); intervening bullets at "
        f"{[n + 1 for n in later if n < steps]}")


def test_neither_stance_can_be_inverted_in_prose():
    """The two things the document exists to say.

    Banning three literal spellings of "adopt one" let a fourth walk past, and
    nothing at all stopped the retraction being flipped to "The claim was TRUE".
    """
    doc = _doc()
    # Stance 1: adopt neither.
    assert "Neither is adopted as correct" in doc
    assert not re.search(
        r"adopt(?:ing)?\s+(?:the\s+)?1\s*min/step\s+(?:everywhere|as correct|"
        r"throughout)|the\s+PK\s+declaration\s+is\s+the\s+correct\s+one|"
        r"in practice,?\s+adopt", doc, re.I), (
        "the doc adopts one reading in prose while its heading says it adopts "
        "neither")
    # Stance 2: the retraction retracts.
    i = doc.index("What this section previously claimed")
    tail = doc[i:i + 1200]
    assert "The claim was false" in tail, (
        "the retraction no longer says the claim was false")
    assert "was TRUE" not in tail and "was correct" not in tail


def test_the_outer_axis_figures_are_derived_from_sim_window():
    """"0 to 28 days" and "minimum timepoint spacing is 6 hours".

    Two engine numbers the prose stated that nothing recomputed: a reviewer
    moved them to 280 days and 60 hours against a green suite. Both are fixed
    by `sim-window`'s own timepoint list.
    """
    doc = _doc()
    main = (SIMS / "sim-window/src/main.rs").read_text()
    m = re.search(r"timepoints_hours:\s*Vec<f64>\s*=\s*vec!\[(.*?)\];",
                  main, re.S)
    assert m, "sim-window no longer declares timepoints_hours"
    body = re.sub(r"//[^\n]*", "", m.group(1))
    hrs = sorted({int(float(x)) for x in re.findall(r"[\d.]+", body)})
    assert len(hrs) >= 5, hrs
    span_days = max(hrs) // 24
    gaps = [b - a for a, b in zip(hrs, hrs[1:])]
    assert f"over 0 to {span_days} days" in doc, (
        f"sim-window spans {max(hrs)} h = {span_days} days; the doc says "
        "otherwise")
    assert f"spacing is {min(gaps)} hours" in doc, (
        f"the smallest gap in {hrs} is {min(gaps)} h; the doc says otherwise")


def test_the_prose_verbs_are_not_inverted():
    """Meaning-level mutations that change no number.

    A reviewer flipped `compete`->`agree`, `coarser`->`finer`,
    `minimum`->`maximum`, `uncalibrated`->`calibrated`, `no caller`->`callers`
    and `agrees with neither`->`with both`, all against a green suite. Each
    pairs a claim with the source that settles it.
    """
    doc = _doc()
    d = _audit()
    declared = max(b["minutes_per_step"] for b in d["step_bindings"]
                   if b["kind"] == "wall-clock")
    implied = max(w["minutes_per_step"] for w in d["implied_windows"])
    assert declared != implied
    assert "readings compete for" in doc and "readings agree on" not in doc, (
        "the doc says two readings AGREE while measuring them "
        f"{implied / declared:g}x apart")
    # trigger_wave really has no caller outside the library re-export.
    callers = [q for q in (SIMS).glob("sim-*/src/*.rs")
               if "trigger_wave" in q.read_text()]
    assert not callers, f"trigger_wave is now called by {callers}"
    assert "no caller in any binary" in doc
    # scd_mufa_rate really is documented uncalibrated, and the doc must not
    # say the opposite.
    params = (CORE / "params.rs").read_text()
    assert "uncalibrated" in params, (
        "params.rs no longer describes any rate as uncalibrated")
    assert "as uncalibrated against" in doc, (
        "the doc no longer records scd_mufa_rate as uncalibrated")
    assert "as calibrated against" not in doc
    # the derived candidate agrees with neither reading
    assert "It agrees with neither reading" in doc
    assert "agrees with both" not in doc
    # The two the docstring named and the body did not check. `6 hours` is the
    # MINIMUM gap (the maximum is 336), so "maximum ... is 6 hours" is flatly
    # false, and the coarser/finer pair is the sentence's whole point.
    assert "coarser than the whole inner assay" in doc and \
        "finer than the whole inner assay" not in doc
    assert "minimum\ntimepoint spacing" in DOC.read_text() or \
        "minimum timepoint spacing" in doc
    assert "maximum timepoint spacing" not in doc


def test_the_detectors_stated_reach_matches_what_it_does():
    """The Status paragraph describes the detector's shape, so it must.

    An unbounded "will now find a third reading if one appears" was the same
    over-claim this section retracts: the detector takes ONE window per line,
    needs the scope verb and the window on one line, and walks `sim-*` only.
    `sim-tme/README.md:138` already carries a second window it drops.
    """
    doc = _doc()
    assert "ONE window per line" in doc
    assert "if one appears" not in doc, (
        "the reach claim is unbounded again")
    # The dropped second window really is on that line.
    line = (SIMS / "sim-tme/README.md").read_text().splitlines()[137]
    assert "0-48h" in line and "1-7 days" in line, (
        "README:138 no longer carries two windows, so the example is stale")
    import sys
    sys.path.insert(0, str(REPO / "scripts"))
    import engine_time_audit as m
    hits = [w for w in m.find_implied_windows()
            if w["module"] == "README.md" and w["binary"] == "sim-tme"]
    assert len(hits) == 1 and hits[0]["window_hours"] == 48.0, (
        "the detector now reports more than one window for that line, so the "
        "stated limitation is stale")


def test_the_citation_is_pinned_to_the_module_it_comes_from():
    doc = _doc()
    tw = (CORE / "trigger_wave.rs").read_text()
    m = re.search(r"PMID[:\s]*(\d{7,8})", tw)
    assert m, "trigger_wave no longer carries a PMID"
    assert f"PMID {m.group(1)}" in doc, (
        f"the doc cites a PMID trigger_wave.rs does not ({m.group(1)})")
