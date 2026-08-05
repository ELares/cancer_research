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
        assert "every number below carries the difference" in txt, (
            "the builds differ by more than the filter and the document does "
            "not say so")
        assert "recorded as unresolved rather than attributed" in txt
    else:
        # Clean A/B. A MeSH-only rule should lose no gene pairs at all; the
        # measured residue is 21 of 3.3M, which is real but below the resolution
        # of any probe run here. It must stay negligible, and it must be
        # explained rather than reported.
        assert d["gene_gene_lost"] < 0.0001 * d["by_namespace"]["gene-gene"]["baseline"], (
            f"{d['gene_gene_lost']} gene-gene pairs lost in a CLEAN A/B; a "
            "MeSH-only filter cannot do that at this scale")
        assert "clean A/B proves it" in txt
        assert "unexplained residue" in txt, (
            "the residual loss is reported without being accounted for")


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


def test_the_control_alias_map_is_read_from_its_log_not_assumed(tmp_path):
    """Assuming the control matched the current code is what hid a confound.

    A rebuilt graph index changed the alias map by 853 forms between two builds,
    and because nothing read what each build ACTUALLY used, a 40,050-pair
    gene-gene loss looked like the filter's doing for several hours.
    """
    import comention_rebuild_compare as c

    log = tmp_path / "run.log"
    log.write_text(
        "loading the entity alias map ...\n"
        "  739,383 aliases -> 44,287 usable (6.0% kept after disambiguation)\n"
        "  [1/28] shard: 1 docs, 2 sentences, 3 pairs so far, 0.1s\n")
    assert c._forms_from_log(log) == 44287
    assert c._forms_from_log(tmp_path / "absent.log") is None


@pytest.mark.skipif(not RAW.exists(), reason="rebuild comparison not run yet")
def test_the_comparison_declares_whether_it_is_clean():
    """A reader must not have to work out whether the two builds are comparable."""
    import json

    d = json.loads(RAW.read_text())
    txt = DOC.read_text()
    assert "Is this a clean A/B?" in txt
    if d.get("confounded"):
        assert "every number below carries the difference" in txt
    else:
        assert "**Yes.**" in txt
        assert d["control_alias_forms"] == d["flagoff_alias_forms"]
