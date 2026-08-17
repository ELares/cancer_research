"""A guard that reads only a generated artifact cannot fail when the generator changes.

THE PATTERN, FOUND SIX TIMES IN ONE DAY
----------------------------------------
Analyses here follow one shape: a script under `scripts/` computes something and
writes a committed report under `analysis/`. The natural guard reads the report
and asserts the numbers and sentences it should contain.

That guard cannot catch a change to the script. The report is a committed file:
it does not move when the generator is edited, only when someone re-runs it. So
a mutation that deletes the central finding from the renderer, hardcodes a
result the code used to derive, drops a filter, or neuters a refusal will leave
every artifact-reading assertion passing.

Measured instances today, each caught only by a mutation sweep:
  * the sharpest-case sentence hardcoded to a chosen modality
  * the check-tag exclusion deleted from a descriptor profile
  * an empty-result refusal replaced with `if False:`
  * a precision-instrument caveat dropped from a renderer
  * partner-tag ranking reverted to a dict that sort_keys reorders
  * the direction-vs-magnitude finding removed from a report

Every one was fixed by adding an assertion against the generator's SOURCE
alongside the one against its output. This makes that structural rather than
remembered.

THE RULE
--------
A test module that knows about BOTH a generator script and its committed
artifact must assert against the script's source at least once. Knowing about
both is what makes it this kind of guard; a test that only reads an artifact
(a pure data check) is out of scope, and so is one that only reads source.

WHAT THIS DELIBERATELY DOES NOT DO
-----------------------------------
It does not check that the source assertion is a GOOD one. A test could satisfy
this by asserting the script contains the word "import". The rule buys the habit,
not the judgement -- and the habit is what was missing six times.
"""

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS = REPO_ROOT / "tests"

# The rule is BLOCKING only for guards written against the analyses added in
# this campaign, where the pattern was verified by mutation. Applying it to the
# whole suite would flag dozens of files at roughly half precision, and the repo
# already learned -- with `audit_heading_only_assertions.py` -- that a gate at
# that precision trains people to add exemptions instead of fixing findings.
#
# The advisory sweep over everything else is `scripts/audit_generator_guards.py`.
# A new analysis should join this list, not the advisory one.
IN_SCOPE = {
    "test_taxonomy_reach.py",
    "test_scope_audit.py",
    "test_ingest_sensitivity.py",
    "test_untagged_partner.py",
    "test_modality_ratio.py",
    "test_engine_time_audit.py",
    "test_recent_window.py",
    "test_baseline_updates.py",
    "test_thesis_rank.py",
    "test_engine_selectivity.py",
    "test_site_coverage.py",
}


def _paths_named(tree: ast.AST) -> tuple:
    """(names a scripts/ path, names an analysis/ artifact) for a module."""
    script = artifact = False
    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
            continue
        v = node.value
        if v.endswith(".py") and "test" not in v:
            script = True
        if v.endswith((".md", ".json", ".tsv.gz")):
            artifact = True
        if v in ("scripts",):
            script = True
        if v in ("analysis",):
            artifact = True
    return script, artifact


def _asserts_on_source(src: str) -> bool:
    """Does the module ever assert against a script's text?

    Deliberately textual: the shapes in use are `SCRIPT.read_text()`,
    `_src()` and `src = SCRIPT.read_text()`, and an AST rule strict enough to
    catch all three is more brittle than the thing it guards.
    """
    return any(marker in src for marker in (
        "SCRIPT.read_text()", "_src()", "SCRIPT).read_text()",
        "script.read_text()", "SRC.read_text()"))


def test_every_generator_guard_pins_its_generator_source():
    offenders = []
    for p in sorted(TESTS.glob("test_*.py")):
        if p.name not in IN_SCOPE:
            continue
        src = p.read_text(encoding="utf-8")
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        names_script, names_artifact = _paths_named(tree)
        if names_script and names_artifact and not _asserts_on_source(src):
            offenders.append(p.name)
    assert not offenders, (
        "these guards read a generated artifact and know about its generator, "
        "but never assert against the generator's source:\n  "
        + "\n  ".join(offenders)
        + "\n\nA committed artifact does not change when its generator is "
          "edited, so such a guard cannot catch a renderer that stops emitting "
          "a finding, a derived value replaced by a literal, or a refusal "
          "neutered to `if False:`. Add one assertion against SCRIPT.read_text() "
          "for whatever the artifact assertion is really about. This pattern "
          "was found six times in a single day, every time by a mutation sweep "
          "rather than by the suite.")


def test_the_detector_would_catch_a_violation():
    """A scan that returns nothing because it is broken looks like a clean suite.

    This repo has that lesson written down, so the detector is exercised against
    a synthetic module that exhibits exactly the pattern.
    """
    bad = ast.parse(
        'SCRIPT = "scripts/x.py"\nMD = "analysis/x.md"\n'
        'def test_a():\n    assert "finding" in open(MD).read()\n')
    s, a = _paths_named(bad)
    assert s and a, "the detector does not recognise a script+artifact module"
    assert not _asserts_on_source(
        'def test_a():\n    assert "finding" in open(MD).read()\n'), (
        "the source-assertion detector fires on a module that has none")
    assert _asserts_on_source(
        'def test_a():\n    assert "x" in SCRIPT.read_text()\n'), (
        "the source-assertion detector misses the standard shape")


def test_the_scope_list_only_holds_real_files():
    """A stale entry silently narrows the rule."""
    missing = [n for n in IN_SCOPE if not (TESTS / n).exists()]
    assert not missing, (
        f"in-scope tests no longer exist: {missing}. Remove them, or the list "
        "is protecting nothing.")
    assert len(IN_SCOPE) >= 6, (
        "the in-scope set has shrunk below the analyses this rule was written "
        "from; narrowing it is how a rule stops applying to anything")
