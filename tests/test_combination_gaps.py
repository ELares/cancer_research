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

2. WHAT THE VARIANT JOIN RECOVERS IS SET BY PARTNER-AGENT CLASS. The first
   version claimed a variant-level absence carries NO information, generalising
   from one failed control. The second claimed the recoveries split by
   BIOMARKER class -- and that was wrong too, because two panel labels were the
   author's to assign and both were assigned in the direction that made the
   split clean. Corrected and extended to 14 regimens, biomarker class does not
   separate (7/12 against 1/2) and partner-agent class does (7/7 against 1/7):
   the extractor does not tie chemotherapy, a hypomethylating agent, a
   checkpoint antibody or an endocrine agent to a substitution.

WHAT THESE GUARD, AND WHY THE PANELS ARE THE THING PINNED
----------------------------------------------------------
A document can be edited to drop a caveat; a panel that stops separating
cannot. So these assert the measured COMPARISON between the two explanations,
not a panel outcome -- an earlier guard pinned "8 of 8 substitution-keyed",
which made adding a real approved regimen a test failure and meant the panel
could only grow in directions that preserved the headline.

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
import re
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
    # The nine-pair gap is very wide, so sitting inside it is a weak test. An
    # independent 36-pair sweep put the real boundary at 2.98 (sorafenib then
    # regorafenib, sequential in HCC), and a threshold below that flagged
    # thousands of rows the wider panel says are not co-administered.
    assert m.CO_ADMIN_RATIO >= 3.0, (
        f"the threshold is {m.CO_ADMIN_RATIO}, below the 2.98 boundary a wider "
        "panel measured; the nine-pair gap admits it but a wider one does not")


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


def test_partner_agent_class_explains_the_recoveries_better_than_biomarker_class():
    """The finding, pinned as a COMPARISON rather than as a panel outcome.

    An earlier guard asserted `substitution_keyed_recovered == total`, which
    pinned the panel's COMPOSITION: adding ivosidenib + azacitidine, a real
    approved regimen, turned the suite red. A guard that makes honest extension
    a failure is a guard that protects a number rather than a claim.

    This one compares two explanations of the same panel, so extending the
    panel can change the verdict -- which is what a measurement should do.
    """
    rec = d()["recovery"]
    bio, par = rec["by_biomarker_class"], rec["by_partner_class"]
    for name, x in (("biomarker", bio), ("partner", par)):
        assert x["yes_total"] and x["no_total"], (
            f"the panel no longer contains both sides of the {name} split, so "
            "it cannot compare anything")
    assert rec["partner_class_explains_better"], (
        f"biomarker class now explains the recoveries at least as well "
        f"(separation {bio['separation']:+.2f}) as partner-agent class "
        f"({par['separation']:+.2f}). The document asserts the opposite and "
        "must be rewritten rather than this guard relaxed")
    assert par["separation"] > 0.5, (
        f"partner-agent separation has fallen to {par['separation']:+.2f}; the "
        "document presents it as the explanation, which a weak separation "
        "would not support")
    # EVERY both-targeted entry must recover, because that is the claim. This
    # DOES pin the panel, and deliberately: unlike the guard it replaced, a
    # both-targeted regimen that misses does not merely change a count, it
    # contradicts the stated mechanism and the document should be revisited.
    # It is also the only guard that catches a panel drug swapped for a
    # different REAL drug, which resolves cleanly and would otherwise ship as
    # a silent false MISSED.
    missed = [r["trial"] for r in d()["recovery"]["panel"]
              if r["both_drugs_are_targeted"] and not r["recovered"]]
    assert not missed, (
        f"these targeted-plus-targeted regimens are not recovered: {missed}. "
        "Either the join broke, a panel drug is the wrong drug, or the "
        "partner-agent explanation does not hold and the document is wrong")


def test_the_two_relabelled_regimens_keep_their_corrected_labels():
    """Both were mine to assign and I assigned both to flatter my own claim.

    FLAURA2 enrols exon 19 deletion OR L858R -- the SAME biomarker as
    MARIPOSA -- so it cannot be gene-keyed while MARIPOSA is substitution-keyed.
    CAPItello-291's indication is PIK3CA/AKT1/PTEN-altered, a three-gene list.
    """
    m = mod()
    by_trial = {p[4]: p for p in m.RECOVERY_PANEL}
    mariposa, flaura = by_trial["MARIPOSA"], by_trial["FLAURA2"]
    assert flaura[3] == mariposa[3], (
        f"FLAURA2 is labelled {flaura[3]!r} and MARIPOSA {mariposa[3]!r}; they "
        "enrol the same biomarker and cannot carry different labels")
    assert by_trial["CAPItello-291"][3] == "", (
        "CAPItello-291 is labelled with a substitution again; its indication "
        "is a three-gene alteration list")


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
    assert "evidence about targeted-plus-targeted regimens" in txt, (
        "the narrowed claim is gone, so the document states the limitation "
        "without saying what it is limited to")
    # The SECOND retracted claim, also stated only to withdraw it.
    for m2 in re.finditer(r"split by BIOMARKER", txt):
        w = txt[max(0, m2.start() - 200): m2.end() + 400]
        assert "wrong too" in w, (
            "the biomarker-class explanation appears without its withdrawal; "
            "it was measured not to separate")


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
    # And the NUMBER must be the one that file actually ships. The claim is
    # about another artifact but was filled from this run's value, so inflating
    # this count made the document assert something false ABOUT A SIBLING
    # REPORT with every guard green.
    n = f"{d()['cotreat_rows']:,}"
    sibling = (REPO_ROOT / "analysis" / "atlas-retraction-exposure.md")
    assert n in sibling.read_text(), (
        f"the report claims atlas-retraction-exposure.md ships {n} cotreat "
        "rows; that file does not contain that number")


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
    # Pinned to the SUBSTANCE, not the bold lead-in. The previous version
    # asserted only the lead, so the measured claim beneath it could be
    # deleted with the guard still green -- the heading-only shape this same
    # test's docstring says it was written to eliminate.
    assert "not measured, at all" in section, (
        "the section no longer says the gene layer's precision is unmeasured")
    assert "No sample has been judged" in section, (
        "the substance of the precision caveat is gone; a bold lead-in alone "
        "lets a reader assume a figure exists somewhere")
    assert "never ran" in section, (
        "the section no longer records that an earlier version reported a "
        "sample this repository did not run, which is the reason the caveat "
        "is worded this strongly")


def test_cotreat_is_not_presented_as_efficacy():
    txt = DOC.read_text()
    section = " ".join(txt[txt.index("## What neither layer can say"):].split())
    assert "No outcome" in section and "worked" in section, (
        "the report no longer states that nothing here says a combination "
        "worked, which is the most likely misreading of every row in it")
