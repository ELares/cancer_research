"""Guards for the corpus-scale alias-name check (#628).

The finding runs against a layer this repo built and against a rule this repo
proposed, so it is the kind of result that quietly softens. Everything here is
pinned to a recomputed quantity or to the committed authority table.
"""

import csv
import gzip
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

RAW = REPO_ROOT / "analysis" / "comention-name-check.json"
DOC = REPO_ROOT / "analysis" / "comention-name-check.md"
TABLE = REPO_ROOT / "analysis" / "comention" / "authority-labels.tsv.gz"


def _raw():
    return json.loads(RAW.read_text())


def test_the_authority_table_is_committed_and_covers_both_namespaces():
    """The blocker this closed: gene identifiers had no labels at all, and half
    the map is genes."""
    assert TABLE.exists()
    ids = []
    with gzip.open(TABLE, "rt") as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                ids.append(line.split("\t", 1)[0])
    genes = [i for i in ids if not i.startswith(("MESH:", "OMIM:"))]
    assert len(ids) > 20000, f"only {len(ids)} identifiers"
    assert len(genes) > 10000, (
        f"only {len(genes)} gene identifiers; the gene half is what was blocked")
    assert _raw()["identifiers_with_labels"] == len(ids)


def test_the_table_carries_aliases_not_just_preferred_names():
    """Comparing against a single preferred label rejected `xCT` for SLC7A11.

    Using every name an authority lists is what dropped the rule's true-positive
    cost from 60% to a third of that.
    """
    import build_label_source as bls

    t = bls.load_table()
    assert "xCT" in t.get("23657", []), "SLC7A11 aliases are missing"
    assert "PHGPx" in t.get("2879", []), "GPX4 aliases are missing"
    # And the collision this repo already documented must be visible in it.
    assert "FSP1" in t.get("51062", []), "ATL1's FSP1 alias is missing"


def test_most_of_the_layers_volume_does_not_sit_on_a_name():
    """The headline. It must be recomputable and must stay a minority claim."""
    d = _raw()
    m, f = d["mentions"], d["forms"]
    assert abs(m["share"] - m["name"] / (m["name"] + m["other"])) < 1e-9
    assert abs(f["share"] - f["name"] / (f["name"] + f["other"])) < 1e-9
    assert 0.4 < m["share"] < 0.7, m["share"]
    # The mention share must stay BELOW the form share -- non-name forms being
    # used more is the unfavourable direction, and an earlier draft asserted the
    # opposite before the numbers were checked.
    assert m["share"] < f["share"], (
        "non-name forms are no longer over-used; the explanation in the document "
        "is stale")


def test_the_discriminator_is_recommended_only_where_it_works():
    """It works on MeSH and does nothing useful on genes."""
    d = _raw()["discriminator"]
    assert d["mesh"]["fp_removed"] > 0.8, d["mesh"]["fp_removed"]
    assert d["mesh"]["kept_precision"] > d["mesh"]["base"]
    # Genes are already precise, which is why the rule has nothing to gain.
    assert d["gene"]["base"] > d["mesh"]["base"] * 2, (
        "the gene subset is no longer much cleaner than MeSH; the reason given "
        "for not applying the rule there is stale")
    assert d["gene"]["fp_removed"] < 0.5
    txt = DOC.read_text()
    assert "apply\nit to MeSH identifiers only" in txt or \
           "MeSH identifiers only" in txt
    assert "weakly supported" in txt, "the small gene sample is not caveated"


def test_the_apoptosis_finding_is_stated_and_true():
    """The alias map has no route to MeSH Apoptosis; every `apoptosis`
    co-mention is filed under a cortical-malformation descriptor."""
    import build_label_source as bls

    t = bls.load_table()
    assert "MESH:D017209" not in t, (
        "Apoptosis is now reachable; the finding in the document is stale")
    assert t.get("MESH:D065703") == ["Malformations of Cortical Development, Group I"]
    assert "no route to MeSH `Apoptosis`" in DOC.read_text()


def test_names_are_not_claimed_to_be_senses():
    """`FSP1` is a listed alias of ATL1, so a name match can still be the wrong
    gene. The document must not present this as a correctness measure."""
    txt = DOC.read_text()
    assert "Names, not senses" in txt
    assert "never whether it is the right" in txt
