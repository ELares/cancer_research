"""Guards for the combination-by-target analysis.

WHAT IT MEASURES
----------------
224,146 `cotreat` rows over 40,878 drug pairs, joined to what each drug is tied
to in the same paper. Two layers: by gene, and by specific variant.

THE FINDING IS A NEGATIVE ONE, AND THE GUARDS EXIST TO KEEP IT NEGATIVE
-----------------------------------------------------------------------
The analysis set out to answer "which resistance mutations have a combination
answer and which do not". At variant level that produces a clean, plausible and
WRONG table: PIK3CA H1047R shows zero, which reads as a therapeutic gap.

Alpelisib plus fulvestrant is FDA-approved for exactly that population, and 59
papers in this corpus assert the co-treatment. Of those, 45 annotate PIK3CA as
a GENE and none annotates a variant. The combination literature writes
"PIK3CA-mutated", not "H1047R", so a variant-level absence measures whether a
paper happened to name a substitution.

WHY THE CONTROL IS GUARDED RATHER THAN THE CONCLUSION
------------------------------------------------------
A document can be edited to remove a caveat. What cannot be edited away is a
measurement that fails: these assert that the control STILL FAILS at variant
level and STILL PASSES at gene level. If the variant join ever starts finding
the control, the negative result is no longer true and the document has to be
rewritten rather than the guard relaxed.

The risk this closes is a real one for this repository, which has twice shipped
a measured gap that turned out to be an artifact: a 24% performance deficit
that was a configuration mistake, and a replication "collapse" that was a
censoring artifact of the author's own window.
"""

import gzip
import importlib.util
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JSON = REPO_ROOT / "analysis" / "atlas-combination-gaps.json"
DOC = REPO_ROOT / "analysis" / "atlas-combination-gaps.md"
TSV = REPO_ROOT / "analysis" / "atlas-combination-gaps.tsv.gz"
SCRIPT = REPO_ROOT / "scripts" / "atlas_combination_gaps.py"


def d() -> dict:
    return json.loads(JSON.read_text())


def flat() -> str:
    return " ".join(DOC.read_text().split())


def mod():
    spec = importlib.util.spec_from_file_location("atlas_combination_gaps", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def rows() -> list:
    with gzip.open(TSV, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh]


# --- the control is the analysis ------------------------------------------

def test_the_control_still_fails_at_variant_level():
    """The negative result. If this passes, the document is wrong."""
    c = d()["control"]
    assert c["papers_asserting_the_combination"] > 0, (
        "the control combination is no longer found in the corpus at all, so "
        "it cannot serve as a control and the analysis needs a new one")
    assert not c["variant_join_passes"], (
        f"the variant-level join now finds the control "
        f"({c['found_by_the_variant_join']} papers). That is good news and it "
        "INVALIDATES this document, whose whole finding is that a variant-level "
        "absence carries no information. Rewrite the analysis rather than "
        "relaxing this guard")
    assert c["of_those_annotating_any_variant"] == 0, (
        "some control paper now carries a variant entity, so the explanation "
        "the document gives for the failure no longer holds")


def test_the_control_passes_at_gene_level():
    """The positive half: the gene join is not merely a looser sieve."""
    c = d()["control"]
    assert c["gene_join_passes"], (
        "the gene-level join no longer finds the control either, so the "
        "document's 'what works instead' section is not supported and the "
        "recommendation to join at gene level should be withdrawn")
    assert c["of_those_annotating_the_gene"] > c["of_those_annotating_any_variant"], (
        "gene annotation is no longer commoner than variant annotation in the "
        "control papers, which is the mechanism the document names")


def test_the_control_is_pinned_by_identifier_not_by_name():
    """A label change must not be able to switch the control off silently."""
    m = mod()
    assert m.CONTROL["drugs"] == ("MESH:C585539", "MESH:D000077267")
    assert m.CONTROL["gene"] == "5290"
    c = d()["control"]
    assert c["drug_ids"] == list(m.CONTROL["drugs"]), (
        "the shipped control identifiers are not the ones the script defines")


def test_the_document_states_the_control_before_the_tables():
    """A reader must meet the invalidation before the numbers it invalidates."""
    txt = DOC.read_text()
    ctrl = txt.index("## The positive control")
    for heading in ("## What works instead", "## The variant layer"):
        assert heading in txt, f"{heading} is missing"
        assert ctrl < txt.index(heading), (
            f"'{heading}' appears above the positive control, so a reader "
            "meets the tables before learning what they can and cannot say")


# --- the join itself -------------------------------------------------------

def test_the_join_requires_both_drugs_tied_to_the_target():
    """Same-paper co-occurrence alone would make the table meaningless."""
    m = mod()
    cotreat = {"p1": {frozenset(("A", "B"))}, "p2": {frozenset(("A", "B"))}}
    drug_target = {
        "p1": {"A": {"g1", "g2"}, "B": {"g1"}},   # only g1 is tied to BOTH
        "p2": {"A": {"g3"}, "B": {"g4"}},         # no shared target at all
    }
    out = m.join(cotreat, drug_target)
    assert set(out) == {(frozenset(("A", "B")), "g1")}, (
        "the join no longer requires BOTH drugs to be tied to the target; a "
        "target mentioned beside only one of them now enters the table")
    assert out[(frozenset(("A", "B")), "g1")] == {"p1"}


def test_rows_and_paper_pair_combinations_are_separate_measures():
    """A paper asserting one pair twice is two rows and one combination.

    The first version labelled the second as the first, which is the same
    defect the variant-map PR shipped and had to correct.
    """
    r = d()
    assert r["cotreat_rows"] > r["paper_pair_combinations"], (
        "co-treatment rows no longer exceed paper-pair combinations, so either "
        "no paper repeats a pair or the two are being counted by the same "
        "expression again")


def test_the_overlap_is_small_enough_to_justify_the_finding():
    r = d()
    assert r["overlap_share_of_cotreatment_papers"] < 10.0, (
        f"the variant and co-treatment layers now overlap on "
        f"{r['overlap_share_of_cotreatment_papers']}% of co-treatment papers; "
        "the document argues an absence is uninformative BECAUSE the overlap "
        "is tiny, so a large overlap means that argument needs revisiting")
    # Both the table row AND the sentence, checked separately. A single
    # substring test over the whole document passes when only one of them is
    # wrong, because the other still supplies a match.
    share = f"{r['overlap_share_of_cotreatment_papers']}%"
    txt = flat()
    assert f"share of the co-treatment literature | {share}" in txt, (
        f"the overlap row in the summary table is not {share}")
    assert f"The join has {share} of the" in txt, (
        f"the overlap share in the prose sentence is not {share}; it and the "
        "table row must both come from the JSON, and checking only that the "
        "number appears somewhere lets either one drift alone")


def test_the_gene_layer_recovers_known_regimens():
    """The check that it measures something real rather than plausible noise.

    These are named trial regimens, verified independently of this corpus. If
    the join stops recovering them it has broken, whatever its row count says.
    """
    by_pair = {}
    for r in rows():
        if r["layer"] != "gene":
            continue
        by_pair.setdefault((r["drug_a"], r["drug_b"]), set()).add(r["target"])
    def find(a, b, gene):
        for (x, y), genes in by_pair.items():
            if {x.lower(), y.lower()} == {a, b} and gene in genes:
                return True
        return False
    for a, b, gene, trial in (
            ("trastuzumab", "pertuzumab", "ERBB2", "CLEOPATRA"),
            ("dabrafenib", "trametinib", "BRAF", "COMBI-d"),
            ("vemurafenib", "cobimetinib", "BRAF", "coBRIM"),
            ("cetuximab", "encorafenib", "BRAF", "BEACON")):
        assert find(a, b, gene), (
            f"the gene layer no longer recovers {a} + {b} against {gene} "
            f"({trial}); that is a real regimen and its absence means the join "
            "has broken rather than that the field changed")


def test_the_caveats_that_bound_the_gene_layer_are_stated():
    """A gene is not a target, and PubTator's assignment is not adjudicated."""
    txt = flat()
    assert "A gene is not a target" in txt, (
        "the report no longer distinguishes a gene the paper ties both drugs "
        "to from the drug's target; KLK3 is PSA, a response biomarker")
    assert "KLK3" in txt and "MAP2K7" in txt, (
        "the two worked examples of that limit are gone, so the caveat is "
        "asserted without an instance a reader can check")
    assert "does not apply" in txt, (
        "the report no longer states that the atlas_ambiguity blocklist does "
        "NOT guard this path, which a reader would otherwise assume it does")


def test_cotreat_is_not_presented_as_efficacy():
    txt = flat()
    assert "not that the combination worked" in txt, (
        "the report no longer states that `cotreat` records co-administration "
        "rather than benefit, which is the most likely misreading of every "
        "row in it")
