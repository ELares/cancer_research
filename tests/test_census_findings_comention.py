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
    txt = DOC.read_text()
    reg = json.loads(REG.read_text())
    m = re.search(r"co-mention precision: \*\*(\d+)% measured\*\*", txt)
    assert m, "the page still reports a bound rather than the measurement"
    assert int(m.group(1)) == round(100 * reg["after"]["weighted"]), (
        f"page says {m.group(1)}%, the measurement is "
        f"{100*reg['after']['weighted']:.1f}%")


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
