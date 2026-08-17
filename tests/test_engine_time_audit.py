"""Guards for the engine time-units audit (#727).

WHAT THIS DOCUMENT NOW CLAIMS
------------------------------
Six declarations bind one engine step to a real duration, in two modules, and
they DISAGREE by a factor of fifty (`tumor_pk` 1 min/step, `trigger_wave`
0.02 min/step). Only `tumor_pk`'s reaches a binary, and `sim-tumor-pk` runs it
for the same 180 steps the production matrix uses -- so the production run spans
three hours, which is shorter than P3's own 24-hour falsification threshold.

WHY THE PREVIOUS GUARDS DID NOT CATCH THE PREVIOUS HEADLINE
-------------------------------------------------------------
The previous headline was "Nothing in the engine states it", and it was FALSE:
`tumor_pk.rs:354` says "Time points in minutes (one per simulation step)". Three
guards stood over that sentence and all three passed.

  `test_the_absence_is_searched_for_not_assumed` asserted
  `declares_step_duration is False` -- it PINNED THE WRONG ANSWER, so the only
  way to fail it was to fix the bug.

  Its "prove the search can succeed" probe re-declared the SAME step-first
  regex and ran it against `"One step is 15 minutes of simulated time."` -- a
  string written to match it. A probe shaped to fit the pattern cannot discover
  that the pattern misses the phrasing the source actually uses. That is a
  guard computing its own expectation.

  `test_the_timescale_span_claim_is_supported` checked the 30-minute and 14-21
  day constants EXIST, which they do. It never checked they were COMPOSED,
  which was the actual claim and was false.

So the guards below are built the other way round: the fixtures come from the
engine source rather than from the regex, and the headline is exercised in both
directions.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "engine_time_audit.py"
MD = REPO_ROOT / "analysis" / "engine-time-audit.md"
JSON_OUT = REPO_ROOT / "analysis" / "engine-time-audit.json"
SRC = REPO_ROOT / "simulations" / "ferroptosis-core" / "src"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("eta", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_detector_sees_the_phrasing_that_fooled_it():
    """The regression test for the actual bug, using the REAL source line.

    The fixture is read out of tumor_pk.rs rather than written here, so it
    cannot be shaped to fit the pattern -- which is exactly how the previous
    version of this guard passed over a live counter-example.
    """
    m = _mod()
    rs = (SRC / "tumor_pk.rs").read_text(errors="ignore")
    real = [ln for ln in rs.split("\n")
            if "one per simulation step" in ln.lower()]
    assert real, (
        "tumor_pk.rs no longer contains the unit-first step declaration this "
        "audit was corrected to detect; re-check the finding against the "
        "current source rather than deleting this test")

    found = m.find_step_bindings()
    assert found, "the detector finds no step binding anywhere in the engine"
    tp = [b for b in found if b["module"] == "tumor_pk.rs"]
    assert tp, (
        "the detector no longer sees tumor_pk.rs's step-to-minute binding -- "
        "the precise blindness that made the previous headline false")
    assert any(b["minutes_per_step"] == 1.0 for b in tp), (
        "tumor_pk's binding is no longer read as one minute per step")


def test_the_detector_is_not_limited_to_one_word_order():
    """Step-first AND unit-first, or the absence is a fact about the regex."""
    m = _mod()
    hows = {h for h, _ in m.BINDINGS}
    assert {"step-first", "unit-first"} <= hows, (
        "the binding patterns no longer cover both word orders; a declaration "
        "written the other way round would read as an absence")
    # exercised, not merely declared
    step_first = [p for h, p in m.BINDINGS if h == "step-first"][0]
    unit_first = [p for h, p in m.BINDINGS if h == "unit-first"][0]
    assert step_first.search("One step is 15 minutes of simulated time.")
    assert unit_first.search("Time points in minutes (one per simulation step).")
    # and neither pattern may match the other's sentence, which would mean one
    # of them is loose enough to fire on anything
    assert not step_first.search("Time points in minutes (one per simulation step).")


def test_doc_comment_runs_are_joined_before_matching():
    """The binding phrase routinely spans two `///` lines."""
    m = _mod()
    sample = ("    /// Returns concentration time-courses at 1-minute\n"
              "    /// resolution (one value per simulation step)\n"
              "    pub fn solve() {}\n")
    joined = [t for _, t in m.logical_lines(sample)]
    assert any("1-minute resolution (one value per simulation step)" in t
               for t in joined), (
        "consecutive doc-comment lines are no longer joined, so a declaration "
        "split across two lines is invisible -- the second half of the bug")


def test_the_headline_can_say_the_opposite():
    """A headline that cannot flip is decoration, not a finding.

    The previous one was a hardcoded 'Nothing in the engine states it' sitting
    above a table that listed the module which states it.
    """
    m = _mod()
    d = _doc()
    absent = m.render({**d, "step_bindings": [], "n_step_bindings": 0,
                       "distinct_minutes_per_step": [],
                       "production_minutes_per_step": []})
    assert "Nothing in the engine states it" in absent, (
        "with no bindings the report no longer reports an absence")
    present = m.render(d)
    assert "Nothing in the engine states it" not in present, (
        "the report still claims nothing states a step duration while "
        f"{d['n_step_bindings']} declarations are in its own artifact")
    assert str(d["n_step_bindings"]) in present


def test_solver_timesteps_are_not_counted_as_wall_clock_bindings():
    """A CFL-constrained integrator dt is not a claim about a step's worth.

    An earlier draft pooled trigger_wave's `dt_min` (a PDE stability step, with
    the CFL assertion `dt < h^2/(2D)` on the line beside it) with tumor_pk's
    step-to-minute alignment and reported "a factor of 50 apart, unreconciled".
    Unlike objects; two solvers having different dt is not a disagreement.
    """
    d = _doc()
    kinds = {b["module"]: b["kind"] for b in d["step_bindings"]}
    assert kinds.get("trigger_wave.rs") == "solver-timestep", (
        "trigger_wave's dt_min is classified as a wall-clock step binding "
        "again; it is a numerical stability parameter")
    assert kinds.get("tumor_pk.rs") == "wall-clock"
    for c in d["wall_clock_conventions"]:
        assert c["module"] != "trigger_wave.rs"
    md = MD.read_text()
    assert "factor of 50 apart, unreconciled" not in md, (
        "the withdrawn 50x disagreement between a PDE dt and a step alignment "
        "is back")


def test_conventions_are_counted_not_regex_hits():
    """`\\bdt_min\\b` matches a field three times; that is one convention."""
    d = _doc()
    assert d["n_wall_clock_conventions"] <= d["n_step_bindings"], (
        "more conventions than matches, which is impossible")
    # the concrete case: one dt_min field, matched three times
    tw = [b for b in d["step_bindings"] if b["module"] == "trigger_wave.rs"]
    assert len(tw) >= 3, (
        "trigger_wave's dt_min no longer produces multiple matches, so this "
        "test no longer exercises the hits-versus-conventions distinction")
    assert len({b["minutes_per_step"] for b in tw}) == 1, (
        "the three trigger_wave matches disagree on the duration; they are "
        "the same field and must collapse to one convention")
    md = MD.read_text()
    assert f"**{d['n_wall_clock_conventions']} module" in md, (
        "the headline count is not the convention count from the artifact")


def test_every_binary_step_count_is_read_regardless_of_int_type():
    """The 'production binary' was being selected by a type annotation.

    `const N_STEPS: usize` matched only sim-tumor-pk; sim-tme-3d, sim-tme and
    sim-combo-mech declare `u32` and were invisible, so changing sim-tme-3d's
    step count left the report unchanged and a retyped sim-combo-mech could
    become "the production matrix" for a binary that never touches the module.
    """
    d = _doc()
    seen = {s["binary"]: s["steps"] for s in d["step_counts"]}
    for path in (REPO_ROOT / "simulations").glob("sim-*/src/main.rs"):
        m = re.search(r"const\s+N_STEPS\s*:\s*\w+\s*=\s*(\d+)", path.read_text())
        if m:
            binary = path.parts[-3]
            assert binary in seen, (
                f"{binary} declares N_STEPS and the audit does not see it")
            assert seen[binary] == int(m.group(1)), (
                f"{binary}: audit says {seen[binary]} steps, source says "
                f"{m.group(1)}")


def test_a_span_is_only_reported_for_a_binary_that_consumes_the_binding():
    """Pricing a run needs the binary to actually use the pricing module."""
    d, md = _doc(), MD.read_text()
    # Derived from the ARTIFACT's own data, never by calling _priced_runs --
    # using the function under test to compute the expectation let a mutation
    # that hardcoded minutes-per-step move both sides together and survive.
    conv = {c["module"]: c["minutes_per_step"] for c in d["wall_clock_conventions"]}
    assert conv, "no wall-clock convention, but the report shows spans"
    priced = 0
    for sc in d["step_counts"]:
        mods = sc["consumes_binding_module"]
        if not mods:
            continue
        src = "".join(p.read_text(errors="ignore") for p in
                      (REPO_ROOT / "simulations" / sc["binary"] / "src").glob("*.rs"))
        for stem in mods:
            assert re.search(rf"\b{re.escape(stem)}\b", src), (
                f"{sc['binary']} is priced by {stem} and does not reference it")
        mps = next(v for k, v in conv.items() if Path(k).stem in mods)
        hours = sc["steps"] * mps / 60
        assert f"| `{sc['binary']}` | {sc['steps']} | {mps:g} | {hours:.1f} h |" in md, (
            f"{sc['binary']} should be priced at {hours:.1f} h "
            f"({sc['steps']} steps x {mps:g} min) and the report does not say so")
        priced += 1
    assert priced, "no binary consumes a pricing module, but spans are shown"


def test_every_reported_caller_actually_references_the_module():
    """The caller column was unverified: fabricating it left the suite green."""
    d = _doc()
    checked = 0
    for mod, callers in d["real_time_module_callers"].items():
        stem = mod[:-3] if mod.endswith(".rs") else mod
        for binary in callers:
            src_dir = REPO_ROOT / "simulations" / binary / "src"
            assert src_dir.is_dir(), f"{binary} is listed as a caller and does not exist"
            src = "".join(p.read_text(errors="ignore") for p in src_dir.glob("*.rs"))
            assert re.search(rf"\b{re.escape(stem)}\b", src), (
                f"{binary} is reported as calling {mod} and does not reference it")
            checked += 1
    assert checked > 0, "no caller relationships were checked"
    for b in d["step_bindings"]:
        stem = b["module"][:-3]
        for binary in b["callers"]:
            src = "".join(p.read_text(errors="ignore") for p in
                          (REPO_ROOT / "simulations" / binary / "src").glob("*.rs"))
            assert re.search(rf"\b{re.escape(stem)}\b", src), (
                f"{b['module']}'s binding claims caller {binary}, which does "
                "not reference it")


def test_the_p3_representation_finding_is_in_the_rendered_report():
    """The whole conclusion could be deleted from the page with tests green.

    The previous guard checked the JSON field and PREREGISTRATION.md and never
    that the RENDERED document said anything about P3 -- so deleting the entire
    consequence from render() left 11/11 passing.
    """
    d, md = _doc(), MD.read_text()
    mod = d.get("p3_modelled")
    assert mod, "the P3-representation measurement is gone from the artifact"
    assert "P3 is represented" in md, (
        "the rendered report no longer carries the P3 finding at all")
    for k, v in mod["recovery_days"].items():
        assert f"`{k}` {v:g}d" in md, (
            f"{k}'s recovery half-time is measured but not rendered")
    if mod.get("sweep"):
        assert f"{mod['sweep']['max_hours']/24:g} days" in md


def test_the_false_p3_impossibility_claim_cannot_return():
    """P3 IS represented; the previous draft said it could not be."""
    d, md = _doc(), MD.read_text()
    if d.get("p3_modelled"):
        assert "8x shorter than the threshold" not in md
        assert "cannot represent either outcome of its own most directly " \
               "testable prediction -- not because" not in md, (
            "the withdrawn claim that the model cannot represent P3 is back, "
            "while cell.rs carries its recovery half-times in days")


def test_the_p3_order_verdict_is_derived_from_both_sources():
    """The contradiction must come from the files, not from prose."""
    d, md = _doc(), MD.read_text()
    o = d.get("p3_order")
    assert o, "the P3 ordering comparison is gone"
    cell = (REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "cell.rs").read_text()
    for name, days in d["p3_modelled"]["recovery_days"].items():
        assert re.search(rf"{name}_half_recovery_days\s*:\s*{days:g}\b", cell), (
            f"{name}'s {days} days is not what cell.rs actually declares")
    preg = (REPO_ROOT / "PREREGISTRATION.md").read_text()
    blk = re.search(r"\*\*P3\..*?(?=\n\*\*P4\.)", preg, re.S)
    assert blk, "P3's block is no longer locatable in the preregistration"
    for nm in o["stated_first"]:
        assert re.search(rf"\b{nm}\b", blk.group(0), re.I), (
            f"{nm} is reported as stated-first but is not in P3's text")

    # RE-DERIVE the verdict from the two raw sources. Branching on the
    # artifact's own `agrees` field is a guard computing its own expectation:
    # forcing `agrees = True` in the generator, or flipping fsp1's default in
    # cell.rs, both left this test green because it simply followed whichever
    # branch the artifact claimed.
    rates = {m.group(1).lower(): float(m.group(2)) for m in
             re.finditer(r"\b(\w+)_half_recovery_days\s*:\s*(\d+(?:\.\d+)?)", cell)}
    grp = re.search(r"\(([^)]*?)\s+first,\s*([^)]*?)\s+later\)", blk.group(0), re.I)
    assert grp, "P3 no longer states a recovery ordering"

    def names(s):
        return {w.lower() for w in re.findall(r"[A-Za-z][A-Za-z0-9]+", s)
                if w.lower() not in ("and", "or")}

    first = {k: rates[k] for k in names(grp.group(1)) if k in rates}
    later = {k: rates[k] for k in names(grp.group(2)) if k in rates}
    assert first and later, "P3's named defenses do not match cell.rs's fields"
    truth = max(first.values()) < min(later.values())
    assert o["agrees"] == truth, (
        f"the artifact says agrees={o['agrees']} but cell.rs "
        f"({rates}) against P3's stated grouping gives {truth}")
    if not truth:
        assert "They **contradict**" in md
        assert f"`{o['violator']}` is stated to recover early" in md
        assert o["violator"] == max(first, key=first.get)
    else:
        assert "They **contradict**" not in md


def test_both_reachable_and_orphan_timescales_are_printed():
    """Printing only the orphans is how the counter-example was dropped."""
    d, md = _doc(), MD.read_text()
    fields = d["orphan_timescale_fields"]
    reach = [(m, f) for m, fs in fields.items() for f, r in fs.items() if r]
    assert reach, "no real-time field is reachable, which contradicts sim-window"
    for mod, f in reach:
        assert f"`{mod}` `{f}`" in md, (
            f"{mod}.{f} is reachable and measured but not rendered -- the "
            "renderer is dropping counter-examples again")


def test_the_core_loop_sentence_is_derived():
    """A fixed sentence would survive biochem.rs acquiring a binding."""
    m, d = _mod(), _doc()
    assert m._core_prices(d) is False, (
        "biochem.rs now prices a step; the report's sentence about the core "
        "loop must be re-read rather than trusted")
    faked = {**d, "wall_clock_conventions": [
        {"module": "biochem.rs", "minutes_per_step": 60.0}]}
    assert "NOW declares one" in m.render(faked), (
        "the core-loop sentence does not change when biochem.rs acquires a "
        "wall-clock binding, so it is decoration")


def test_it_does_not_choose_or_reconcile_a_step_duration():
    """Choosing one moves every calibrated layer and the byte-identity gates."""
    md, src = MD.read_text(), SCRIPT.read_text()
    assert "does not choose a step duration" in md
    assert not re.search(r"\bSTEP_MINUTES\b|\bstep_duration\s*=\s*\d", src), (
        "the audit has started asserting a step duration, which is a modelling "
        "decision and not this script's to make")
    # it must not write into the engine either
    assert "ferroptosis-core" not in src.split("OUT_JSON")[1], (
        "the generator reaches into the engine source after defining its "
        "outputs; this audit measures and must not modify")


def test_the_false_absence_claim_cannot_return():
    """The retracted headline, pinned by IDENTIFIER not by substring.

    The generator quotes the old sentence in order to withdraw it, so a bare
    substring ban would trip on the retraction itself -- the shape that has now
    caught guards in this repo four separate times.
    """
    md = MD.read_text()
    d = _doc()
    if d["n_step_bindings"] > 0:
        assert "Nothing in the engine states it" not in md, (
            "the report asserts nothing states a step duration while its own "
            "artifact lists declarations that do")
        assert "three orders of magnitude" not in md, (
            "the withdrawn span sentence is back; it was built on a "
            "composition that was never measured")


def test_the_composition_claim_is_measured_at_field_level():
    """Module-level reachability is too coarse and passed over a false claim."""
    d = _doc()
    fields = d.get("orphan_timescale_fields")
    assert fields, "the field-level composition check is gone"
    orphans = {m: [f for f, r in fs.items() if not r] for m, fs in fields.items()}
    orphans = {m: v for m, v in orphans.items() if v}
    if orphans:
        md = MD.read_text()
        tot = sum(len(v) for v in orphans.values())
        assert f"**{tot} real-time configuration fields are referenced" in md, (
            "the report does not state how many declared timescales are "
            "referenced nowhere, which is what refutes the composition claim")
        for mod in orphans:
            assert f"`{mod}`" in md


def test_an_empty_scan_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'if not d["modules"]:' in src
    assert "is not a finding" in src and "raise SystemExit" in src
