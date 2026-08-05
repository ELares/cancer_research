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
def test_the_filter_did_not_leak_into_gene_pairs():
    """The invariant. A MeSH-only rule that loses gene-gene pairs is broken,
    and the layer's only consumer is 19/20 gene-gene."""
    import json

    d = json.loads(RAW.read_text())
    assert d["gene_gene_leak"] == 0, (
        f"{d['gene_gene_leak']} gene-gene pairs lost to a MeSH-only filter; this "
        "is a leak, not a design choice")


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
