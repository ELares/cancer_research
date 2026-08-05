"""Guards for the before/after pair comparison (#628).

The comparison exists to answer one question the offline prediction could not:
did the MeSH-only filter leak into gene pairs? A leak there is a correctness
bug, not a tuning choice, so it is asserted rather than reported.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

RAW = REPO_ROOT / "analysis" / "comention-rebuild-compare.json"
DOC = REPO_ROOT / "analysis" / "comention-rebuild-compare.md"


def test_the_namespace_classifier_is_right():
    """MeSH, OMIM and bare gene ids must be told apart, since the whole rule
    turns on it."""
    import comention_rebuild_compare as c

    assert c.ns("MESH:D001943") == "mesh"
    assert c.ns("OMIM:612348") == "omim"
    assert c.ns("2879") == "gene"
    assert c.ns("51062") == "gene"


@pytest.mark.skipif(not RAW.exists(), reason="rebuild comparison not run yet")
def test_the_gene_pair_movement_is_split_and_the_confound_declared():
    """The invariant as first written was wrong, and the data said so.

    "Gene-gene pairs must be unchanged" is not a property a MeSH-only filter
    has: removing a MeSH alias frees the tokens it consumed, so a shorter gene
    alias can match at the same position and CREATE pairs. The net figure was
    negative, which a leak check reading it as a loss would have called a pass
    in the wrong direction.

    What genuinely needs explaining is the LOSS, and it is not attributable to
    the filter, because the two builds also differ by a rebuilt index. The
    document must declare that rather than assign the loss to the rule.
    """
    import json

    d = json.loads(RAW.read_text())
    assert "gene_gene_lost" in d and "gene_gene_gained" in d, (
        "the net figure hides two opposite movements and must be split")
    assert d["gene_gene_gained"] > 0, "unmasking produced no new gene pairs"
    txt = DOC.read_text()
    if d.get("confounded"):
        assert "TWO changes, not one" in txt, (
            "the builds differ by more than the filter and the document does "
            "not say so")
        assert "recorded as unresolved rather than attributed" in txt
    else:
        assert d["gene_gene_lost"] == 0, (
            f"{d['gene_gene_lost']} gene-gene pairs lost with no confound to "
            "explain them; a MeSH-only filter cannot do that")


@pytest.mark.skipif(not RAW.exists(), reason="rebuild comparison not run yet")
def test_the_comparison_reports_gained_pairs_not_only_lost():
    """The filter is not purely subtractive: removing a long alias unmasks a
    shorter surviving one, because the matcher consumes tokens. A comparison
    that only counted losses would misreport the change."""
    import json

    d = json.loads(RAW.read_text())
    assert "gained" in d and "lost" in d
    assert d["rebuilt_pairs"] == d["baseline_pairs"] - d["lost"] + d["gained"], (
        "the pair arithmetic does not close; lost/gained are not a partition")
