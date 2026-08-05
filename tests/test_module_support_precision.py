"""The co-mention column in atlas-module-support must carry its error bound.

That column is load-bearing: the document uses it to argue a zero in the
relation column is an extraction failure rather than absence of evidence. The
layer's precision is read from its own artifact rather than stated here -- it
has moved from ~47% to 88% across two rebuilds, and a figure written into this
docstring would be the very drift these tests exist to catch. The argument is
made most often on single-digit counts, where whatever the error rate is can
account for a large share of the figure, so presenting the counts without the
bound overstates the evidence.
"""

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))
DOC = REPO_ROOT / "analysis" / "atlas-module-support.md"
REG = REPO_ROOT / "analysis" / "comention-regression.json"


def test_the_column_is_populated_or_declared_unbuilt():
    txt = DOC.read_text()
    built = "co-mention column is measured at roughly" in txt
    unbuilt = "co-mention layer is not built" in txt
    assert built ^ unbuilt, "the doc must state exactly one of the two states"


def test_the_stated_precision_tracks_the_layer_as_built():
    """It must quote the precision of the layer this column came from.

    The authority filter is on by default since #628, so the figure is the
    filtered layer's. Reading the unfiltered one would have this document
    quoting 42% for a layer measured at 88% and understating its own evidence
    by half.
    """
    txt = DOC.read_text()
    if "co-mention layer is not built" in txt:
        return
    m = re.search(r"roughly (\d+)% precision", txt)
    assert m, "the co-mention precision is not stated"
    auth = REPO_ROOT / "analysis" / "comention-authority-result.json"
    measured = (json.loads(auth.read_text())["weighted"] if auth.exists()
                else json.loads(REG.read_text())["after"]["weighted"])
    assert abs(int(m.group(1)) / 100 - measured) < 0.02, (
        f"doc says {m.group(1)}%, the layer as built measures {100*measured:.1f}%")


def test_small_counts_are_flagged_where_they_carry_the_argument():
    """A handful of co-mentions is a handful of chances for the error rate."""
    txt = DOC.read_text()
    if "co-mention layer is not built" in txt:
        return
    import atlas_module_support as ms

    thin = ms._thin_threshold()
    for line in txt.splitlines():
        m = re.search(r"\(but \*\*([\d,]+)\*\* full-text co-mentions", line)
        if m and int(m.group(1).replace(",", "")) < thin:
            assert "still thin against" in line, (
                f"a single-digit-to-tens count is offered as evidence without "
                f"the caveat: {line[:120]}")
