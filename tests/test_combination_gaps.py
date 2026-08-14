"""Guards for the co-treatment layer analysis.

WHAT IT MEASURES
----------------
224,146 `cotreat` rows joined to what each drug is tied to in the same paper,
at gene and at variant level. Two things had to be established before the join
could be read at all, and both became findings.

1. `cotreat` IS NOT CO-ADMINISTRATION. Alectinib and crizotinib are sequential
   ALK inhibitors, never given together, and carry 86 `cotreat` rows. The
   corpus's own `compare` predicate supplies a discriminator, and the ratio
   separates a panel where the truth is known independently.

2. AN ABSENCE IS INFORMATIVE ONLY ABOUT SUBSTITUTION-KEYED REGIMENS. The first
   version claimed a variant-level absence carries NO information, generalising
   from one failed control. Measured against ten approved regimens the variant
   join recovers eight, and the two misses are both keyed to a gene-level
   biomarker rather than a substitution.

WHAT THESE GUARD, AND WHY THE PANELS ARE THE THING PINNED
----------------------------------------------------------
A document can be edited to drop a caveat; a panel that stops separating
cannot. So these assert the measured splits (8 of 8 substitution-keyed, 0 of 2
gene-keyed; every co-administered pair above every sequential one) rather than
the sentences describing them.

THE FAILURE THAT MADE THE PANEL RESOLUTION A HARD ERROR
--------------------------------------------------------
The first version of the panels hardcoded MeSH identifiers from memory and got
several wrong: `MESH:C571179` is eravacycline, not vemurafenib. Every wrong one
shipped as a row reading MISSED, indistinguishable from a real negative, and
dropped measured recovery from 8/10 to 2/10 while looking like a result. Names
are now resolved through the authority table and an unresolvable one raises.
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


# --- 1. cotreat semantics --------------------------------------------------

def test_the_panel_shows_cotreat_is_not_co_administration():
    """The finding: pairs that are never given together carry cotreat rows."""
    sem = d()["cotreat_semantics"]
    seq = [x for x in sem["panel"] if not x["co_administered"]]
    assert seq, "the sequential half of the panel is gone"
    assert all(x["cotreat"] > 0 for x in seq), (
        "no sequential pair carries a cotreat row any more; if PubTator's "
        "predicate was corrected upstream this document's premise is void and "
        "it should be rewritten rather than the guard relaxed")
    assert sem["pairs_where_compare_outnumbers_cotreat"] > 0


def test_the_ratio_separates_the_panel_with_no_overlap():
    """A discriminator that overlaps cannot license the flag column."""
    sem = d()["cotreat_semantics"]
    assert sem["ratio_separates_the_panel"], (
        f"the cotreat:compare ratio no longer separates: lowest "
        f"co-administered {sem['lowest_co_administered_ratio']}, highest "
        f"sequential {sem['highest_sequential_ratio']}. The flag column and "
        "the prose that justifies it both have to go")
    assert sem["lowest_co_administered_ratio"] > sem["highest_sequential_ratio"]
    # The threshold must sit INSIDE the measured gap, or it is a number
    # someone picked rather than one the panel supports.
    m = mod()
    assert sem["highest_sequential_ratio"] < m.CO_ADMIN_RATIO \
        <= sem["lowest_co_administered_ratio"], (
        f"the {m.CO_ADMIN_RATIO} threshold is outside the gap the panel "
        f"measures ({sem['highest_sequential_ratio']} to "
        f"{sem['lowest_co_administered_ratio']}), so it is not supported by it")


def test_the_ratio_is_reported_never_used_to_filter():
    """A nine-pair panel cannot support a classifier, and the doc says so."""
    r = d()
    flagged = sum(1 for x in rows() if x["reads_as_co_administration"] == "0")
    assert flagged > 0, (
        "no row is flagged compare-dominated any more; either the flag stopped "
        "being computed or the rows were filtered out, and filtering on a "
        "nine-pair heuristic is what the document promises not to do")
    assert r["gene_layer"]["triples"] > r["gene_layer"]["reading_as_co_administration"], (
        "every gene row now reads as co-administration, so the flag is inert")
    assert "not a validated classifier" in flat(), (
        "the report no longer states the ratio's status, so a reader may take "
        "the flag as a verdict")


# --- 2. what the variant join recovers -------------------------------------

def test_every_panel_regimen_resolves_to_exactly_one_drug_pair():
    """An unresolvable name ships as a false MISSED.

    This is not hypothetical: hardcoded identifiers put eravacycline where
    vemurafenib belonged and fluticasone where panitumumab belonged, which read
    as four extra misses and took measured recovery from 8/10 to 2/10.
    """
    m = mod()
    lab = m.load_labels()
    assert lab, "the authority table is missing, so no panel entry can resolve"
    for name in {n for p in m.RECOVERY_PANEL for n in p[:2]} | \
                {n for p in m.SEMANTIC_PANEL for n in p[:2]}:
        ident = m.resolve_drug(name, lab)   # raises if it does not resolve
        assert lab[ident].lower() == name.lower()
    # And the FAILURE path, which no panel name can exercise because they all
    # resolve. Without this the raising behaviour is untested: a version that
    # returned an empty identifier instead of raising passed the loop above
    # unchanged, because the difference is invisible on names that resolve.
    import pytest
    for missing in ("a drug that does not exist", ""):
        with pytest.raises(SystemExit):
            m.resolve_drug(missing, lab)


def test_the_recovery_split_is_the_finding():
    """8 of 8 substitution-keyed, 0 of 2 gene-keyed. Pinned as a split."""
    rec = d()["recovery"]
    assert rec["substitution_keyed_total"] > 0 and rec["gene_keyed_total"] > 0, (
        "the panel no longer contains both kinds of regimen, so it cannot "
        "measure the split that is the whole finding")
    assert rec["substitution_keyed_recovered"] == rec["substitution_keyed_total"], (
        f"only {rec['substitution_keyed_recovered']} of "
        f"{rec['substitution_keyed_total']} substitution-keyed regimens are "
        "recovered; the document claims the join works for these, so a miss "
        "means either the join broke or a panel entry is wrong")
    assert rec["gene_keyed_recovered"] == 0, (
        f"{rec['gene_keyed_recovered']} gene-keyed regimen(s) are now "
        "recovered at variant level. That is good news and it invalidates the "
        "mechanism this document gives for the failure. Rewrite it rather than "
        "relaxing this guard")


def test_the_document_does_not_reinstate_the_retracted_claim():
    """The first version said an absence carries NO information. It was wrong.

    Its own top variant row, dabrafenib + trametinib against BRAF V600E at 181
    papers, refutes it.
    """
    import re
    txt = flat()
    # PROXIMITY, not absence. The document states the retracted claim in order
    # to retract it, so a bare `not in` fires on the retraction itself -- the
    # trap this repository records as pinning a guard to a substring rather
    # than to what it means. Every occurrence must sit beside its withdrawal.
    hits = list(re.finditer(r"carries NO information", txt))
    for m in hits:
        window = txt[max(0, m.start() - 300): m.end() + 300]
        assert "earlier version" in window and "too strong" in window, (
            "the retracted categorical claim appears without its withdrawal "
            f"beside it; measured recovery is {d()['recovery']['recovered']} "
            f"of {d()['recovery']['of']}")
    assert len(hits) == 1, (
        f"the retracted claim appears {len(hits)} times; it belongs in exactly "
        "one place, the sentence that withdraws it")
    assert "evidence about substitution-keyed regimens" in txt, (
        "the narrowed claim is gone, so the document states the limitation "
        "without saying what it is limited to")


def test_the_control_explains_the_failure_mode():
    """The control is kept because it makes the mechanism visible."""
    c = d()["control"]
    assert c["papers_asserting_the_combination"] > 0, (
        "the control combination is no longer in the corpus")
    assert c["of_those_annotating_any_variant"] == 0, (
        "some control paper now carries a variant entity, so the explanation "
        "the document gives for the miss no longer holds")
    assert c["of_those_annotating_the_gene"] > c["of_those_annotating_any_variant"], (
        "gene annotation is no longer commoner than variant annotation in the "
        "control papers, which is the mechanism the document names")
    assert c["found_by_the_gene_join"] > 0, (
        "the gene join no longer finds the control either, so reporting it "
        "beside the variant layer is not supported")


# --- the join itself -------------------------------------------------------

def test_the_join_requires_both_drugs_tied_to_the_target():
    m = mod()
    cotreat = {"p1": {frozenset(("A", "B"))}, "p2": {frozenset(("A", "B"))}}
    drug_target = {"p1": {"A": {"g1", "g2"}, "B": {"g1"}},
                   "p2": {"A": {"g3"}, "B": {"g4"}}}
    out = m.join(cotreat, drug_target)
    assert set(out) == {(frozenset(("A", "B")), "g1")}, (
        "the join no longer requires BOTH drugs to be tied to the target")
    assert out[(frozenset(("A", "B")), "g1")] == {"p1"}


def test_the_shipped_layers_match_the_shipped_table():
    """The summary counts must BE the table's, not a parallel computation."""
    r, rs = d(), rows()
    gene = [x for x in rs if x["layer"] == "gene"]
    var = [x for x in rs if x["layer"] == "variant"]
    assert r["gene_layer"]["triples"] == len(gene), (
        f"the JSON reports {r['gene_layer']['triples']:,} gene triples and the "
        f"table ships {len(gene):,}")
    assert r["variant_layer"]["triples"] == len(var)
    assert r["gene_layer"]["single_paper"] == sum(
        1 for x in gene if x["papers"] == "1")


def test_rows_and_paper_pair_combinations_are_separate_measures():
    r = d()
    assert r["cotreat_rows"] > r["paper_pair_combinations"], (
        "co-treatment rows no longer exceed paper-pair combinations, so the "
        "two are being counted by the same expression again")


def test_the_overlap_figure_appears_in_both_places_it_is_quoted():
    r = d()
    assert r["overlap_share_of_cotreatment_papers"] < 10.0
    share = f"{r['overlap_share_of_cotreatment_papers']}%"
    assert f"({share})" in flat(), f"the overlap row is not {share}"


# --- what the document must keep saying ------------------------------------

def test_the_provenance_claim_is_narrowed_to_what_is_true():
    """"Nothing had read them" was FALSE and had to be narrowed to "joined".

    `atlas-retraction-exposure.md` already ships the 224,146 cotreat count and
    `atlas-emergence.md` already ships per-pair cotreat counts. This is the
    second time in two PRs that a "nothing had read this" claim needed
    narrowing, so it is checked against the files rather than trusted.
    """
    txt = flat()
    assert "not new" in txt, (
        "the report implies the cotreat counts are new here; two committed "
        "reports already carry them")
    others = [p.name for p in (REPO_ROOT / "analysis").glob("atlas-*.md")
              if p.name != "atlas-combination-gaps.md"
              and "cotreat" in p.read_text(errors="ignore")]
    assert others, (
        "no other committed report mentions cotreat any more, so the "
        "'not new here' framing describes something no longer true")
    assert any(o.replace(".md", "") in txt for o in others), (
        f"the report says the counts are not new without naming where they "
        f"already appear; they are in {others}")


def test_the_caveats_that_bound_the_gene_layer_are_stated_in_prose():
    """Checked in the PROSE section, not the whole document.

    An earlier version asserted `"KLK3" in txt`, which the generated TABLE
    satisfies -- so both worked-example sentences could be deleted and the
    guard still passed.
    """
    txt = DOC.read_text()
    start = txt.index("## What neither layer can say")
    section = " ".join(txt[start:].split())
    assert "A gene is not a target" in section
    for token in ("KLK3", "MAP2K7", "SLTM", "Neptunium"):
        assert token in section, (
            f"{token} is no longer in the caveats section; the collision "
            "classes are asserted without an instance a reader can check")
    assert "does not apply" in section, (
        "the section no longer states that the atlas_ambiguity blocklist does "
        "NOT guard this path")
    assert "not measured at scale" in section, (
        "the section no longer says the gene layer's precision is unmeasured, "
        "so the sampled bound reads as a rate")


def test_cotreat_is_not_presented_as_efficacy():
    txt = DOC.read_text()
    section = " ".join(txt[txt.index("## What neither layer can say"):].split())
    assert "No outcome" in section and "worked" in section, (
        "the report no longer states that nothing here says a combination "
        "worked, which is the most likely misreading of every row in it")
