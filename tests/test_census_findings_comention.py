"""The census findings page must report the MEASURED co-mention precision.

Two failure modes this guards:

1. Reporting the corroboration BOUND (44.3%-91.8%) once a hand-judged
   measurement exists. The bound is what the layer had before anyone read its
   output; keeping it overstates what is known, and here it also understates how
   bad the layer is.
2. Omitting the regression from "What the census did NOT support". A repair this
   project made and justified, which measurement then showed made things worse,
   is the single most load-bearing entry that section can carry -- and the one
   most likely to quietly disappear.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DOC = REPO_ROOT / "analysis" / "census-findings.md"
REG = REPO_ROOT / "analysis" / "comention-regression.json"


def test_the_measured_precision_supersedes_the_bound():
    """The page must quote a MEASUREMENT, and it must be of the layer as SHIPPED.

    This guard originally pinned the figure to `comention-regression.json`,
    which was correct while that document described the running layer. Once the
    authority filter was promoted (#628) that artifact described a
    configuration nothing runs, and a bound on a build nobody uses is not a
    bound -- so the guard is pinned to whichever artifact describes the shipped
    layer, preferring the authority result when it exists.
    """
    txt = DOC.read_text()
    auth = REPO_ROOT / "analysis" / "comention-authority-result.json"
    expected = (json.loads(auth.read_text())["weighted"] if auth.exists()
                else json.loads(REG.read_text())["after"]["weighted"])
    m = re.search(r"co-mention precision: \*\*(\d+)% measured\*\*", txt)
    assert m, "the page still reports a bound rather than the measurement"
    assert int(m.group(1)) == round(100 * expected), (
        f"page says {m.group(1)}%, the shipped layer measures {100*expected:.1f}%")


def test_the_negative_result_is_not_left_at_its_first_half():
    """A closed arc reported only up to the failure is its own kind of stale.

    'A repair made it worse' was true and complete when written. It stopped
    being complete when the real fix shipped, and a findings page that keeps
    only the half that went wrong is as misleading as one that keeps only the
    half that went right -- which is the failure mode this page exists to name.
    """
    auth = REPO_ROOT / "analysis" / "comention-authority-result.json"
    if not auth.exists():
        return
    d = json.loads(auth.read_text())
    section = _section("## What the census did NOT support")
    assert "made it worse" in section, "the failure itself must still be reported"
    # Anchored to its own SENTENCE, not to the bare number. The number alone
    # also appears in "Every layer now carries a bound" further down the page,
    # so a slice-to-EOF assertion on it passed even with this whole paragraph
    # deleted -- and passed with the figure rewritten to contradict the artifact.
    assert f"the fix is measured at {round(100*d['weighted'])}%" in section, (
        "the page reports the regression without reporting that it was then fixed")
    # The lesson, not just the recovery: what the fix cost is part of the finding.
    assert f"{round(100*d['recall_cost']['tp_lost_share'])}% of true matches" in section


def _section(heading: str) -> str:
    """One section of the page, bounded at the next horizontal rule.

    Slicing to EOF makes every assertion below it satisfiable by text from a
    LATER section, which is how the guard above was briefly vacuous.
    """
    txt = DOC.read_text()
    start = txt.index(heading)
    end = txt.find("\n---\n", start + len(heading))
    return txt[start:end if end != -1 else len(txt)]


def test_the_regression_is_listed_under_what_was_not_supported():
    txt = DOC.read_text()
    start = txt.index("## What the census did NOT support")
    section = txt[start:]
    assert "made it worse" in section, (
        "the co-mention regression is not reported among the negative results")
    # It must name the mechanism, not just the number -- the transferable part.
    assert "moved the pressure to the channel just opened" in section


def test_the_page_does_not_still_claim_the_layer_improved():
    """The pre-measurement framing said the filters repaired precision."""
    txt = DOC.read_text()
    reg = json.loads(REG.read_text())
    assert reg["net_change"] < 0
    assert "co-mention precision: 44.3% to 91.8%" not in txt
