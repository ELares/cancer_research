"""Guards for the engine time-units audit (#727).

THE CLAIM
---------
Nothing in the engine states what one step is worth in real time, while modules
carrying real time -- a 30-minute drug half-life, a 24-48 hour photosensitizer,
a 14-21 day senescence programme -- are composed into the same 180-step run.
The consequence is that PREREGISTRATION.md's P3, stated in days, cannot be
scored against the model that produced it.

WHAT WOULD MAKE THIS WRONG
--------------------------
1. IF SOMETHING DOES DECLARE A STEP DURATION. The whole claim is an absence, and
   an absence is the easiest thing to assert carelessly. The audit SEARCHES for
   a declaration rather than assuming there is none, and this pins that the
   search is real by checking it can find one when one exists.

2. IF THE REAL-TIME MODULES ARE NOT ACTUALLY COMPOSED. Two systems in one
   codebase are only a problem if a single run uses both. The claim rests on
   composition, not co-existence.

3. IF THE ISSUE'S ORIGINAL WORDING CREPT BACK. It said the engine has "four
   conflicting step-duration bindings". Inspection found none: there are modules
   with their own units that were never reconciled. "Definitions disagree"
   implies four choices were made; nobody made any. The report must keep saying
   the absence is the finding.
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


def test_the_absence_is_searched_for_not_assumed():
    """An absence asserted is worth nothing; an absence searched for is a result."""
    m = _mod()
    assert _doc()["declares_step_duration"] is False, (
        "something now declares a step duration -- which is the outcome this "
        "audit exists to prompt. Update the report: the finding has been fixed")
    # and prove the search can succeed, or 'none found' means nothing
    probe = re.compile(
        r"(one|a|each|per)\s+step\s+(is|=|represents?|corresponds? to|equals?)\s*"
        r"[^.\n]{0,40}\b(min|minute|hour|h|day|sec)", re.I)
    assert probe.search("One step is 15 minutes of simulated time."), (
        "the declaration pattern cannot match an obvious declaration, so "
        "'none found' is a statement about the regex rather than the engine")
    # and the SEARCH must still traverse real files. The artifact does not move
    # when the scan is edited, so a scan neutered to look at nothing would keep
    # reporting declares_step_duration=False and pass everything above.
    assert 'for p in list(SRC.glob("*.rs"))' in SCRIPT.read_text(), (
        "the declaration search no longer walks the engine source; "
        "'none found' would then be a statement about an empty loop")


def test_both_time_systems_are_actually_present():
    """Co-existence is the premise; without both sets the claim is empty."""
    d = _doc()
    assert d["real_time_only"], "no module carries real time units"
    assert d["per_step_only"], "no module carries per-step rates"
    assert d["modules_total"] > 20, "the module scan found suspiciously few files"
    # the specific ones the report names must really carry real time
    for expected in ("tumor_pk.rs",):
        assert expected in d["real_time_only"] + d["carry_both"], (
            f"{expected} no longer carries a real-time constant, but the "
            "report's argument names it")


def test_the_real_time_constants_are_real():
    """Spot-check against the source, so the scan is not matching noise."""
    d = _doc()
    for name in d["real_time_only"] + d["carry_both"]:
        ex = d["modules"][name]["examples"]
        assert ex, f"{name} is classified real-time with no example"
        src = (SRC / name).read_text(errors="ignore")
        assert ex[0]["context"][:40] in src, (
            f"{name}: the quoted context is not in the file, so the audit is "
            "reporting something the source does not say")


def test_the_timescale_span_claim_is_supported():
    """Three orders of magnitude is a strong sentence and needs its evidence."""
    md = MD.read_text()
    assert "three orders of magnitude" in md
    joined = " ".join(
        e["context"] for v in _doc()["modules"].values() for e in v["examples"])
    assert re.search(r"\b30\s*min", joined, re.I), (
        "no 30-minute constant found; the span claim's short end is missing")
    assert re.search(r"1[0-9]\s*to\s*21\s*day|14 to 21 day", joined, re.I), (
        "no multi-week constant found; the span claim's long end is missing")


def test_it_does_not_invent_a_step_duration():
    """Choosing one belongs to whoever owns the calibrated layers."""
    md, src = MD.read_text(), SCRIPT.read_text()
    assert "does not choose a step duration" in md, (
        "the report no longer disclaims choosing a duration")
    assert not re.search(r"\bSTEP_MINUTES\b|\bstep_duration\s*=", src), (
        "the audit has started asserting a step duration, which is a modelling "
        "decision and not this script's to make")


def test_the_withdrawn_four_bindings_claim_does_not_return():
    """It said four definitions conflict. There are none, and that is different."""
    md = MD.read_text()
    # NOT a bare substring check. The generator's docstring QUOTES the phrase
    # in order to withdraw it, so forbidding the string outright fails on the
    # correction itself -- the same shape as an earlier guard today that
    # forbade "not a change here" and then tripped on the sentence retracting
    # it. What must not return is the CLAIMING form, so the quote is required
    # to stay attached to its withdrawal.
    src = SCRIPT.read_text()
    if "four conflicting" in src.lower():
        assert "an earlier version of this issue said" in src.lower(), (
            "the 'four conflicting step-duration bindings' phrase appears in "
            "the generator without the withdrawal that makes it historical")
    assert "four conflicting" not in md.lower(), (
        "the withdrawn 'four conflicting step-duration bindings' framing is "
        "back; inspection found no definitions at all, and 'definitions "
        "disagree' implies choices nobody made")
    assert "the absence is the finding" in md.lower() or \
        "absence is the finding" in SCRIPT.read_text().lower(), (
        "the report no longer distinguishes an absence from a disagreement")


def test_the_prediction_consequence_is_named_specifically():
    """'Units matter' is a truism; naming the unfalsifiable prediction is not."""
    md = MD.read_text()
    assert "P3" in md, (
        "the report no longer names which preregistered prediction this makes "
        "unscoreable, which is what turns a units complaint into a finding")
    preg = (REPO_ROOT / "PREREGISTRATION.md").read_text()
    assert re.search(r"3 to 7 days", preg), (
        "P3 no longer states a day-scale window, so the report's example is "
        "stale and the claim needs re-checking against the current text")


def test_an_empty_scan_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'if not d["modules"]:' in src
    assert "is not a finding" in src and "raise SystemExit" in src
