"""The contradiction totals must live in the artifact, not in a terminal.

`atlas_contradictions.py` printed its two headline counts to stdout and stored
only the truncated top-500 lists. Those counts were quoted downstream, having
been copied out of a terminal by hand -- so nothing could check them, and
nothing noticed when they drifted about 50% (4,667 -> 7,068 direction,
6,764 -> 9,094 valence) as the disambiguation corrections grew and merged more
entities.

The committed lists are byte-identical across that change, because the biggest
conflicts are stable while the tail grows. So no drift guard over the artifact
could have caught it either. The only fix is for the number to exist.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "analysis" / "atlas-contradictions.json"
DOC = REPO_ROOT / "analysis" / "atlas-contradictions.md"
INDEX = REPO_ROOT / "CLAUDE.md"


def test_the_totals_are_in_the_artifact():
    d = json.loads(RAW.read_text())
    for k in ("direction_conflicts", "valence_conflicts"):
        assert isinstance(d.get(k), int) and d[k] > 0, (
            f"{k} is not recorded; it exists only wherever someone last ran the "
            "script, which is how the previous figures went 50% stale unnoticed")
    assert d["direction_conflicts"] >= len(d["direction"]), (
        "the total is smaller than the list it summarises")
    assert d["valence_conflicts"] >= len(d["valence"])


def test_the_document_states_the_totals_it_lists_from():
    d, txt = json.loads(RAW.read_text()), " ".join(DOC.read_text().split())
    assert f"{d['direction_conflicts']:,} direction conflicts" in txt
    assert f"{d['valence_conflicts']:,} valence" in txt
    assert f"top {d['listed']}" in txt, (
        "the document does not say the tables are truncated, so a reader may "
        "take the listed rows for the whole result")


def test_the_index_quotes_the_artifact():
    """CLAUDE.md is where the stale pair sat for months."""
    d = json.loads(RAW.read_text())
    txt = INDEX.read_text()
    m = re.search(r"asserts in BOTH directions \(([\d,]+)\)", txt)
    assert m, "the index no longer states the direction total"
    assert int(m.group(1).replace(",", "")) == d["direction_conflicts"], (
        f"the index says {m.group(1)}, the artifact says "
        f"{d['direction_conflicts']:,}")
    m2 = re.search(r"cause a disease \(([\d,]+)\)", txt)
    assert m2 and int(m2.group(1).replace(",", "")) == d["valence_conflicts"]
