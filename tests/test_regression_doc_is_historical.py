"""A historical comparison must not drift with a live artifact.

`analysis/comention-regression.md` measures ONE change: the #617 filter swap.
Its "before" side was pinned as a constant; its "after" side was read live from
`analysis/atlas-comention-audit.json`. That was correct exactly until the audit
was regenerated for a DIFFERENT build -- which happened when the authority
filter shipped (#646).

From then on, re-running the generator recomputed a historical regression
against a build it does not describe, reporting precision RISING 10.9 points
where the finding is that it FELL 9. Nothing failed. The document rewrote its
own conclusion, and it would have shipped the moment anyone regenerated it.

These guards pin the property: the comparison's endpoints are constants, the
document still reports the direction it was written to report, and the generator
notices when the live artifact has moved on rather than silently consuming it.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

DOC = REPO_ROOT / "analysis" / "comention-regression.md"
AUDIT = REPO_ROOT / "analysis" / "atlas-comention-audit.json"


def test_both_endpoints_of_the_comparison_are_pinned():
    import comention_regression as cr

    for name in ("PRE_REBUILD_AUDIT", "POST_617_AUDIT"):
        d = getattr(cr, name, None)
        assert isinstance(d, dict), f"{name} is not a pinned constant"
        for k in ("mentions", "pubtator_agree", "in_abstract", "body_only"):
            assert isinstance(d.get(k), int), f"{name}[{k!r}] is not a fixed count"
    assert cr.PRE_REBUILD_AUDIT["mentions"] != cr.POST_617_AUDIT["mentions"], (
        "the two endpoints describe the same build")


def test_the_generator_does_not_consume_the_live_audit_as_its_after_side():
    """The defect itself: reading whatever the artifact currently says.

    Checked on the source, because the failure is silent at run time -- the
    generator produced a confident, wrong document and exited 0.
    """
    src = (REPO_ROOT / "scripts" / "comention_regression.py").read_text()
    assert "after = json.loads(AUDIT.read_text())" not in src, (
        "the after side is read live again; regenerating this document will "
        "recompute the #617 regression against whatever build the audit "
        "currently describes")


def test_the_live_audit_has_in_fact_moved_on():
    """Proof the hazard is real here and not hypothetical.

    If this ever fails because the two agree again, the guard above still
    stands; this one exists to document that they diverged.
    """
    if not AUDIT.exists():
        return
    import comention_regression as cr

    live = json.loads(AUDIT.read_text())
    assert live["mentions"] != cr.POST_617_AUDIT["mentions"], (
        "the live audit matches the pinned post-#617 build; if that is genuine "
        "the pin is harmless, but check that the audit was not reverted")


def test_the_document_still_reports_the_regression_it_measures():
    txt = DOC.read_text()
    m = re.search(r"fell from ([\d.]+)% to ([\d.]+)%, about ([\d.]+) points", txt)
    assert m, "the headline no longer states a fall"
    before, after, drop = float(m.group(1)), float(m.group(2)), float(m.group(3))
    assert after < before, (
        f"the document says precision went {before}% -> {after}%, which is a "
        "RISE; #617 made the layer worse and this document exists to say so")
    assert abs((before - after) - drop) < 0.15, (
        f"the stated drop of {drop} points does not match {before} - {after}")
