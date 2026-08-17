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


def test_the_conflict_is_derived_from_the_numbers():
    """Two durations disagreeing is the finding; asserting it would not be."""
    d = _doc()
    vals = d["distinct_minutes_per_step"]
    assert len(vals) >= 2, (
        "the engine no longer declares two different step durations. If they "
        "were reconciled, that is the outcome this audit exists to prompt -- "
        "update the report rather than deleting the test")
    md = MD.read_text()
    ratio = max(vals) / min(vals)
    assert f"factor of {ratio:.0f} apart" in md, (
        "the report does not state the measured ratio between the conflicting "
        "step durations")
    src = SCRIPT.read_text()
    assert 'if len(d["distinct_minutes_per_step"]) > 1:' in src, (
        "the disagreement is no longer conditional on the measurement, so the "
        "sentence would print whether or not they disagree")


def test_the_production_span_is_read_from_source_not_written():
    """180 steps and 1 min/step are facts about files, not prose."""
    d = _doc()
    prod = d["production"]
    assert prod and prod["steps"] > 0
    main = (REPO_ROOT / "simulations" / prod["binary"] / "src" / "main.rs")
    assert re.search(rf"const\s+N_STEPS\s*:\s*usize\s*=\s*{prod['steps']}\b",
                     main.read_text()), (
        f"{prod['binary']} no longer runs {prod['steps']} steps; the span "
        "figure in the report is stale")
    span = prod["steps"] * min(d["production_minutes_per_step"]) / 60
    assert f"{span:.1f} hours" in MD.read_text(), (
        "the report does not state the span its own artifact implies")


def test_the_p3_consequence_is_parsed_from_the_preregistration():
    """The threshold must come from PREREGISTRATION.md, not from memory."""
    d = _doc()
    p3 = d.get("p3") or {}
    assert "threshold_hours" in p3, (
        "P3's falsification threshold is no longer parseable from "
        "PREREGISTRATION.md; the consequence paragraph rests on it")
    preg = (REPO_ROOT / "PREREGISTRATION.md").read_text()
    blk = re.search(r"\*\*P3\..*?(?=\n\*\*P4\.)", preg, re.S)
    assert blk, "P3's block is no longer locatable in the preregistration"
    assert re.search(rf"within\s+{p3['threshold_hours']:g}\s*hours",
                     blk.group(0), re.I), (
        "the parsed threshold is not the one the preregistration states")


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
