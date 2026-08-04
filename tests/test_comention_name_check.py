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


def test_the_unreachable_class_is_measured_and_reported_as_near_empty():
    """A hypothesis raised and then refuted by counting.

    `Apoptosis` being locked out by a cortical-malformation descriptor looked
    like the tip of a damaging class. It is not: most unreachable identifiers
    are unreachable because nobody writes their name, or because a
    near-identical entity holds it. The negative bounds how far the Apoptosis
    finding should be read, so it must not quietly drop out.
    """
    u = _raw()["unreachable"]
    assert u["annotated"] > 200000
    assert u["reachable"] < u["annotated"]
    n = u["examined"]
    assert n > 2000, f"only {n} examined"
    # The benign explanations must dominate; if the different-entity class ever
    # grows past a fifth, the document's framing is stale and must be re-read.
    benign = u["name_unused"] + u["name_to_same"]
    assert benign / n > 0.8, f"benign causes are only {100*benign/n:.0f}%"
    assert u["name_to_different"] / n < 0.2
    txt = DOC.read_text()
    assert "close to empty" in txt
    assert "outlier rather than a symptom" in txt, (
        "the bound on the Apoptosis finding is missing")


def test_the_cross_namespace_cases_are_not_called_corruption():
    """`Prostate-Specific Antigen` losing to `KLK3` is redundancy, not capture.

    Counting those as an entity being locked out by a competitor would overstate
    the finding; they are the same biological entity in two vocabularies.
    """
    txt = DOC.read_text()
    assert "not what it looks like either" in txt
    assert "cross-namespace redundancy" in txt
    assert "are not competitors" in txt


def test_the_insertion_point_cost_is_reported_not_just_the_diagnostic():
    """A filter inside `build_alias_map` sees the FORM, not the sentence.

    The document's headline rows score canonical-form OR matched-span, which is
    a fair diagnostic and not what a filter could use. Quoting only that
    understated the true-positive cost by half again, and it is the number a
    reader would build against.
    """
    d = _raw()["at_insertion_point"]
    strict, norm = d["strict"], d["normalised"]
    assert strict["tp_removed"] > _raw()["discriminator"]["mesh"]["tp_removed"], (
        "the form-only cost is no longer higher than the diagnostic; if that is "
        "real the document's warning is stale")
    # Cleaner but more expensive: it removes every span-bearing false positive.
    assert strict["fp_removed"] > 0.99 and strict["kept_precision"] > 0.95
    flat = " ".join(DOC.read_text().split())
    assert "form only, as a filter would" in flat
    assert "favourable number is the one this document originally quoted" in flat


def test_normalising_recovers_the_cost_at_no_precision_loss():
    """The named change that makes the rule buildable."""
    d = _raw()["at_insertion_point"]
    strict, norm = d["strict"], d["normalised"]
    assert norm["tp_removed"] < strict["tp_removed"], "normalising no longer helps"
    assert norm["kept_precision"] >= strict["kept_precision"] - 1e-9
    assert norm["fp_removed"] >= strict["fp_removed"] - 1e-9
    assert norm["kept"] > strict["kept"]


def test_the_cancer_vocabulary_cost_is_measured():
    """Precision on judged mentions cannot see this, and it is decisive.

    The strict rule leaves most of MeSH tree C04 -- the cancer definition the
    census is built on -- unreachable. A rule that doubles precision by deleting
    the vocabulary it exists to index is not an improvement.
    """
    c = _raw()["c04_cost"]
    assert c["strict"] and c["normalised"], "the C04 cost is not computed"
    assert c["strict"]["mass_retained"] < 0.4, (
        "the strict rule no longer destroys the cancer vocabulary; re-read")
    assert c["normalised"]["mass_retained"] > c["strict"]["mass_retained"] * 1.5
    assert c["normalised"]["descriptors_killed"] < c["strict"]["descriptors_killed"]
    # The generator wraps prose, so collapse whitespace before matching a phrase
    # that can straddle a line break.
    flat = " ".join(DOC.read_text().split())
    assert "should not be built" in flat, (
        "the document does not warn against building the strict form")


def test_the_table_carries_mesh_entry_terms():
    """A preferred label is not how the literature writes a concept.

    MeSH says `Breast Neoplasms`; papers say `breast cancer`. Without entry
    terms the rule leaves 260 of 670 cancer descriptors unreachable.
    """
    import build_label_source as bls

    t = bls.load_table()
    breast = t.get("MESH:D001943", [])
    assert len(breast) > 3, f"Breast Neoplasms has only {len(breast)} names"
    assert breast[0] == "Breast Neoplasms", "the preferred label must come first"
    joined = " | ".join(breast).lower()
    assert "cancer" in joined, "no cancer-worded entry term for Breast Neoplasms"
    # Across the table, entry terms must be a substantial addition.
    mesh_names = sum(len(v) for k, v in t.items() if k.startswith("MESH:"))
    mesh_ids = sum(1 for k in t if k.startswith("MESH:"))
    assert mesh_names / mesh_ids > 1.2, (
        f"only {mesh_names/mesh_ids:.2f} names per MeSH identifier; entry terms "
        "are missing")


def test_stopwords_are_dropped_because_mesh_stores_inverted_forms():
    """`Cancer of Breast` must match `breast cancer` as a bag."""
    import comention_name_check as c

    assert c.norm_bag("Cancer of Breast") == c.norm_bag("breast cancer")
    assert c.norm_bag("Breast Neoplasms") == c.norm_bag("breast cancer")
    # But it must not collapse genuinely different concepts.
    assert c.norm_bag("lung cancer") != c.norm_bag("breast cancer")


def test_the_recommended_rule_beats_the_strict_one_on_both_axes():
    """Precision is held while the two costs both fall. If either reverses, the
    recommendation in the document is stale."""
    d = _raw()
    i, c = d["at_insertion_point"], d["c04_cost"]
    assert i["normalised"]["tp_removed"] < i["strict"]["tp_removed"]
    assert i["normalised"]["kept_precision"] >= i["strict"]["kept_precision"] - 1e-9
    assert i["normalised"]["fp_removed"] >= i["strict"]["fp_removed"] - 1e-9
    assert c["normalised"]["mass_retained"] > c["strict"]["mass_retained"]
    assert c["normalised"]["mass_retained"] > 0.65, (
        f"cancer vocabulary retention fell to "
        f"{100*c['normalised']['mass_retained']:.1f}%; the rule is no longer safe "
        "to build")
