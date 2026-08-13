"""Guards for the drug-by-variant map.

WHAT IT MEASURES
----------------
`relations.tsv.gz` carries `ProteinMutation`/`DNAMutation`/`SNP` entities whose
identifier holds the variant inline. Those entities already flow through this
repository as opaque node labels; what nothing here had ever done is PARSE the
fields. The map is the drug-by-variant slice of them, reconciled and resolved.

THE FOUR THINGS THAT CAN GO WRONG, and why each gets a guard
-------------------------------------------------------------
1. IT READS AS A CLINICAL RESOURCE. It is machine-extracted co-assertion with
   no directionality: `associate` covers "confers resistance" and "predicts
   response" identically. A reader who takes a row as clinical guidance is the
   most probable way this document causes harm, so the caveat has to reach them
   before the table does.

2. RECONCILIATION OVER-MERGES. This has already happened once. The first rule
   tested agreement among PROTEIN spellings only and swept the rest onto the
   winner, so under rs1801131 a 17-row `p.E429A` captured `c.1298A>C` (710)
   together with `c.1286A>C` and `c.1298A>T` -- a different position and a
   different allele. It over-merged 232 rsids and 7,764 rows while the document
   presented it as the safety property. The worse version of the same mistake
   is available: rs121913529 covers KRAS G12D, G12V AND G12A, and a merged row
   would still look entirely plausible.

3. RECONCILIATION UNDER-MERGES. The un-reconciled map read JAK2 V617F at a
   seventh of its support, because PubTator emits `p.V61F` -- the canonical MPN
   driver with a position digit dropped -- for 86% of that variant's rows.

4. THE ARTIFACT DESCRIBES A DECISION THE CODE DID NOT MAKE. Also already
   happened: a form that is a digit-drop twin of two different canonicals had
   its correction decided by dict-write order, and the JSON shipped a verdict
   that was not applied.

WHY THE JAK2 GUARD KEYS ON THE RSID AND NOT ON THE COUNTS
---------------------------------------------------------
Counts cannot adjudicate a digit-drop twin, and reading the majority as correct
is wrong in both directions: TP53 `p.R72P` and ERBB2 `p.I655V` are real rs-backed
polymorphisms whose longer twins are the typos. These check that the report gives
the independent evidence, not merely that it reached the right answer.

WHY THERE IS A NO-HANDWRITTEN-FIGURES GUARD
-------------------------------------------
The first version hand-wrote its worked examples and shipped "rs77375493 on 620
of 624 rows" when the truth is 624 of 624 -- in a paragraph whose neighbouring
figures were freshly computed, which is what made the stale one credible.
"""

import gzip
import importlib.util
import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JSON = REPO_ROOT / "analysis" / "atlas-variant-drug-map.json"
DOC = REPO_ROOT / "analysis" / "atlas-variant-drug-map.md"
TSV = REPO_ROOT / "analysis" / "atlas-variant-drug-map.tsv.gz"
SCRIPT = REPO_ROOT / "scripts" / "atlas_variant_drug_map.py"


def d() -> dict:
    return json.loads(JSON.read_text())


def flat() -> str:
    """Document text with whitespace collapsed.

    Content assertions must not depend on where a line happens to wrap; two
    guards in this repository have already been broken by a reflow.
    """
    return " ".join(DOC.read_text().split())


def mod():
    spec = importlib.util.spec_from_file_location("atlas_variant_drug_map", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def rows() -> list:
    with gzip.open(TSV, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh]


# --- 1. it must not read as a clinical resource ----------------------------

def test_the_curated_databases_are_named_before_any_drug_variant_row():
    """A reader meets the caveat first or the caveat is not doing its job."""
    txt = DOC.read_text()
    for db in ("CIViC", "OncoKB", "COSMIC"):
        assert db in txt, f"{db} is no longer named; it is the correct source"
    # Anchored on the table HEADER, not on whichever drug happens to rank first.
    # `txt.index("| osimertinib")` raised a bare ValueError the moment the top
    # row changed, discarding every failure message below it.
    header = "| drug | gene | variant | papers |"
    assert header in txt, "the drug-variant table header has changed shape"
    first_row = txt.index(header)
    for db in ("CIViC", "OncoKB", "COSMIC"):
        assert txt.index(db) < first_row, (
            f"{db} is named only BELOW the drug-variant table, so a reader "
            "meets the rows before learning what the curated sources are")


def test_the_absent_directionality_is_stated_not_implied():
    """`associate` covers resistance AND sensitivity, which is the whole risk.

    NOT an `or` over paraphrases. Three guards in this repository have been
    satisfied by a different fix because an `or` let one arm stand in for the
    other, so this pins the substantive claim itself.
    """
    txt = flat()
    assert "carries **no direction**" in txt, (
        "the report no longer states that the predicate carries no direction")
    assert "confers RESISTANCE" in txt and "predicts RESPONSE" in txt, (
        "the report no longer spells out that resistance and response land in "
        "the same bucket, which is the specific thing a reader must not "
        "misread a high count as resolving")


def test_the_single_paper_share_is_reported():
    r = d()
    share = r["single_paper_pairs"] / r["distinct_drug_variant_pairs"]
    assert share > 0.5, (
        f"single-paper pairs are now {share:.1%} of the map; the report frames "
        "this as the dominant shape, so a drop below half means that framing "
        "needs rewriting rather than the guard relaxing")
    assert f"{100 * share:.0f}% of pairs rest on a single paper" in flat(), (
        "the single-paper share in the prose is not the one in the JSON")


# --- 2. reconciliation must not over-merge ---------------------------------

def test_the_agreement_test_covers_every_spelling_not_just_protein_forms():
    """The defect that shipped: 232 rsids, 7,764 rows, silently over-merged."""
    m = mod()
    import collections
    fix, tally, refused = m.resolve_rsids({
        # THE REGRESSION CASE, from the corpus. One protein spelling at 17 rows
        # against coding spellings at a different position (1286) and a
        # different allele (A>T). The superseded rule saw one protein triple,
        # concluded "no disagreement", and captured all 742 rows onto p.E429A.
        ("1801131", "4524"): collections.Counter({
            "c.1298A>C": 710, "p.E429A": 17, "c.-1298A>C": 3,
            "c.1286A>C": 3, "c.1298A>T": 3, "c.1298C>T": 2}),
    })
    assert ("1801131", "4524") not in fix, (
        "a 17-row protein spelling captured 742 rows spanning a different "
        "position and a different allele; the agreement test is back to "
        "looking at protein forms only")
    assert tally.get("classes_disagree_refused") == 1
    assert refused[0]["old_rule_would_have_collapsed"], (
        "the record of what the superseded rule would have done is no longer "
        "computed, so the size of that defect stops being a derived number")
    # The report shows the offending spellings; a truncated set once flagged a
    # genuinely disagreeing spelling as clean.
    assert "c.1286A>C" in refused[0]["offending"]
    assert "c.1298C>T" in refused[0]["offending"], (
        "a disagreeing spelling is missing from the offending set, so the "
        "report would display it with no flag and read as though it were fine")
    assert "p.E429A" not in refused[0]["offending"], (
        "the lone protein spelling is flagged as offending; it agrees with "
        "itself and is the form the old rule wrongly merged ONTO")


def test_the_multi_allelic_refusal_is_still_protecting_something():
    """If this reaches zero the rule has stopped distinguishing anything."""
    t = d()["rsid_resolution"]
    assert t.get("classes_disagree_refused", 0) > 0, (
        "no rsID is refused for disagreeing spellings any more, so either the "
        "corpus changed shape or the rule collapsed into an unconditional "
        "merge, which would silently fuse KRAS G12C with G12D")
    assert t.get("one_change_many_spellings", 0) > 0, (
        "no rsID is collapsed either, so the rule is inert and the map is back "
        "to counting `p.T790M` and `c.2369C>T` as two variants")


def test_kras_codon_twelve_substitutions_stay_distinct():
    """The load-bearing safety property, checked on the shipped table.

    rs121913529 covers G12D, G12V and G12A; rs121913530 covers G12C, G12R and
    G12S. Any rule keyed on the rsID would fuse them, and the fused row would
    still look entirely plausible.
    """
    kras = {r["variant"] for r in rows() if r["gene"] == "KRAS"}
    for v in ("p.G12C", "p.G12D", "p.G12V"):
        assert v in kras, (
            f"KRAS {v} is absent from the map; it shares an rsID with the other "
            "codon-12 substitutions, so its disappearance means they merged")
    by_var = {}
    for r in rows():
        if r["gene"] == "KRAS":
            by_var.setdefault(r["variant"], set()).add(r["drug"])
    assert by_var["p.G12C"] != by_var["p.G12D"], (
        "KRAS G12C and G12D now carry an identical drug set, which is what a "
        "merge on the shared rsID would produce")


def test_the_resolver_refuses_a_multi_allelic_site_and_collapses_a_respelling():
    m = mod()
    import collections
    fix, tally, _ = m.resolve_rsids({
        # one change, protein and coding spellings -> collapse to the protein form
        ("121434569", "1956"): collections.Counter({"p.T790M": 4146, "c.2369C>T": 9}),
        # three substitutions at one codon -> refuse
        ("121913529", "3845"): collections.Counter({"p.G12D": 1830, "p.G12V": 873,
                                                    "p.G12A": 126}),
        # a truncated/HTML-escaped spelling cannot be checked -> refuse
        ("113488022", "673"): collections.Counter({"p.V600E": 12488, "c.1799T&gt": 7}),
        # a lone spelling -> pass through unchanged
        ("77375493", "3717"): collections.Counter({"p.V617F": 620}),
    })
    assert fix[("121434569", "1956")] == "p.T790M", (
        "the coding and protein spellings of one change did not collapse onto "
        "the protein form")
    assert ("121913529", "3845") not in fix, (
        "a multi-allelic rsID was resolved to a single HGVS; this is the merge "
        "that would fuse KRAS G12D with G12V")
    assert ("113488022", "673") not in fix, (
        "an rsID carrying an unparsable spelling was collapsed; the rule can "
        "say nothing about a string it cannot compare, so it must refuse")
    assert fix[("77375493", "3717")] == "p.V617F"
    assert tally["classes_disagree_refused"] == 1
    assert tally["unparsable_spelling_refused"] == 1
    assert tally["one_change_many_spellings"] == 1


def test_a_minority_protein_spelling_cannot_capture_a_majority():
    """The property in one line, on synthetic input, with no corpus dependency."""
    m = mod()
    import collections
    fix, _, _ = m.resolve_rsids({
        ("1", "1"): collections.Counter({"p.A1B": 1, "c.100A>C": 999,
                                         "c.200A>C": 1}),
    })
    assert fix == {}, (
        "a 1-row protein spelling captured 1,000 rows whose coding spellings "
        "name two different positions")


# --- 3. reconciliation must not under-merge --------------------------------

def test_the_jak2_correction_is_applied():
    m = mod()
    import collections
    gene_hgvs = {"3717": {
        "p.V61F": [3747, collections.Counter()],
        "p.V617F": [622, collections.Counter({"77375493": 620})],
    }}
    fix, twins = m.adjudicate_twins(gene_hgvs)
    assert fix.get(("3717", "p.V61F")) == "p.V617F", (
        "the digit-drop twin is no longer corrected, so JAK2 V617F reads at a "
        "fraction of its support")
    assert len(twins) == 1 and twins[0]["moves_majority"], (
        "the correction is no longer flagged as moving the majority onto the "
        "minority string, which is the case a reader most needs to see")


def test_the_rule_abstains_when_the_rsid_cannot_adjudicate():
    """Abstention is the point: counts alone are not evidence."""
    m = mod()
    import collections
    gene_hgvs = {"9999": {
        "p.A11B": [500, collections.Counter({"1": 1})],
        "p.A111B": [3, collections.Counter({"2": 1})],   # both carry rsIDs
    }, "8888": {
        "p.C22D": [500, collections.Counter()],
        "p.C222D": [3, collections.Counter()],           # neither does
    }}
    fix, twins = m.adjudicate_twins(gene_hgvs)
    assert fix == {}, (
        "a twin was corrected although the rsIDs cannot adjudicate it; the "
        "rule would then be deciding on counts, which is wrong in both "
        "directions (TP53 p.R72P is the short form and correct, EGFR p.T790M "
        "is the long form and correct)")
    assert all(t["verdict"].startswith("abstain:") for t in twins)


def test_a_form_one_digit_from_two_canonicals_abstains():
    """KRAS `p.G1V` is one digit from both `p.G12V` and `p.G13V`.

    Those are different substitutions. The first version let dict-write order
    pick, silently deciding by ASCII sort and reporting in the JSON a verdict
    it had not applied.
    """
    m = mod()
    import collections
    none, rs = collections.Counter(), collections.Counter({"7": 1})
    fix, twins = m.adjudicate_twins({"3845": {
        "p.G1V": [2, none], "p.G12V": [873, rs], "p.G13V": [14, rs],
    }})
    assert ("3845", "p.G1V") not in fix, (
        "an ambiguous digit-drop form was assigned a canonical; G12V and G13V "
        "are different substitutions and this is a guess")
    ambiguous = [t for t in twins if t["shorter"] == "p.G1V"]
    assert ambiguous and all("abstain:" in t["verdict"] for t in ambiguous), (
        "the ambiguous pairs do not record an abstention verdict")


def test_the_artifact_reports_only_decisions_the_code_actually_made():
    """The defect class: a verdict in the JSON that was never applied."""
    tw = d()["twins"]
    assert tw["corrections_applied"] == tw["adjudicated"], (
        f"the report claims {tw['adjudicated']} adjudicated twins but only "
        f"{tw['corrections_applied']} corrections were applied, so a reader "
        "auditing a verdict from the artifact is auditing a decision the code "
        "did not make")


def test_the_jak2_verdict_is_justified_by_the_rsid_not_by_the_counts():
    """The report must give the independent evidence, not just the answer."""
    txt = flat()
    assert "rs77375493" in txt, (
        "the report no longer cites the rsID that adjudicates JAK2 V617F, so "
        "it is asserting the correction rather than evidencing it")
    i = txt.index("rs77375493")
    window = txt[max(0, i - 700): i + 700]
    assert "carries no rsid at all" in window.lower(), (
        "the report cites the canonical form's rsID without stating that the "
        "malformed twin carries none, which is the half that makes it evidence")
    assert "Polycythemia Vera" in window, (
        "the corroborating disease profile is no longer stated beside the "
        "rsID; the two together are what make the verdict independent of the "
        "row counts, which point the other way")


def test_the_majority_moving_rule_reports_how_thin_its_evidence_is():
    """Most of these rest on a single rsID-bearing row. Listing without saying
    so presents four equally-supported corrections when only one is."""
    movers = d()["twins"]["moving_the_majority"]
    assert movers, "no twin is reported as moving the majority any more"
    assert all("canonical_rsid_rows" in t for t in movers), (
        "the majority-moving table no longer carries the row count its rsID "
        "evidence rests on")
    assert min(t["canonical_rsid_rows"] for t in movers) < 5, (
        "every majority-moving correction now rests on substantial evidence; "
        "if that is real the prose calling most of them thin must be rewritten")
    assert "rsid-bearing rows" in flat(), (
        "the evidence column is missing from the shipped table")


# --- 4. no figure may be hand-written --------------------------------------

# Figures that must appear in the document AND be derivable from the JSON.
# Keyed by the JSON path, valued by how it is written in prose.
def _derived_figures(r: dict) -> list:
    tw = r["twins"]
    f, ref = tw["featured"], r["rsid_refusals"]["featured"]
    old = r["rsid_refusals"]["the_old_rule_would_have_collapsed"]
    out = [
        (f"{tw['found']} twin pairs", "twins.found"),
        (f"{r['rsid_rows_given_an_hgvs']:,}", "rsid_rows_given_an_hgvs"),
        (f"{old['rsids']:,} of the rsids refused", "old_rule.rsids"),
        (f"{old['rows']:,} rows", "old_rule.rows"),
    ]
    if f:
        out += [(f"rs{f['canonical_rsid']} on {f['canonical_rsid_rows']:,} of its "
                 f"{f['canonical_rows']:,} rows", "twins.featured rsid support")]
    if ref:
        out += [(f"rs{ref['rsid']}", "rsid_refusals.featured rsid"),
                (f"{ref['rows']:,} rows", "rsid_refusals.featured rows")]
    return out


def test_every_headline_figure_in_the_prose_comes_from_the_artifact():
    r, txt = d(), flat()
    missing = [(s, where) for s, where in _derived_figures(r) if s not in txt]
    assert not missing, (
        "these figures are in the JSON but not written as stated in the prose, "
        "so the prose is carrying a hand-typed number: "
        + "; ".join(f"{w} -> {s!r}" for s, w in missing))


def test_the_generator_does_not_hardcode_the_worked_example_figures():
    """The first version hardcoded "620 of 624" and it was wrong.

    Checked against the generator source, because a figure that is correct
    today and hand-written is the shape this repository keeps rediscovering:
    the artifact moves, the sentence does not, and the stale number reads as
    freshly checked -- the more so when the figures beside it ARE derived.
    """
    src = SCRIPT.read_text()
    r = d()
    f = r["twins"]["featured"]
    forbidden = ["620 of 624", "12,488", "13,069"]
    if f:
        forbidden += [f"{f['canonical_rsid_rows']:,} of its"]
    present = [lit for lit in forbidden if lit in src]
    assert not present, (
        f"the generator hardcodes {present}; these must be derived from the "
        "data so the example cannot go stale against the artifact beside it")


# --- provenance and readability --------------------------------------------

def test_the_unread_claim_is_narrowed_to_what_is_true():
    """"Absent from every other analysis" was FALSE and its guard checked
    a different file entirely.

    `analysis/atlas-emergence.md` already ships a drug-variant row carrying the
    raw identifier, for the same adagrasib/KRAS-G12C pair this map reports. What
    is true is narrower: nothing else PARSES the inline fields.
    """
    others = [p for p in (REPO_ROOT / "analysis").rglob("*.md")
              if p.name != "atlas-variant-drug-map.md"
              and "CorrespondingGene" in p.read_text(errors="ignore")]
    txt = flat()
    if others:
        assert "not absent" in txt, (
            f"{len(others)} other committed report(s) already carry variant "
            f"identifiers (e.g. {others[0].name}), so the document must not "
            "imply the entities are absent from this repository")
        assert "atlas-emergence.md" in txt, (
            "the report claims the entities appear elsewhere without naming "
            "where, so a reader cannot check it")
    parsers = [p.name for p in (REPO_ROOT / "scripts").glob("*.py")
               if "CorrespondingGene" in p.read_text(errors="ignore")]
    assert set(parsers) <= {"atlas_variant_drug_map.py", "atlas_discovery.py"}, (
        f"{parsers} now reference the inline variant fields; the claim that "
        "nothing else parses them needs re-checking against whatever changed")


def test_a_bare_rsid_in_the_variant_column_carries_its_prefix():
    """`113488022` is indistinguishable from a codon position."""
    bad = [r for r in rows() if re.fullmatch(r"\d+", r["variant"])]
    assert not bad, (
        f"{len(bad)} rows show a bare numeric variant, e.g. "
        f"{bad[0]['drug']}/{bad[0]['gene']}/{bad[0]['variant']}; an rsID in "
        "that column reads as a residue position")


def test_the_obsolete_prefix_is_not_shown_beside_an_approved_drug():
    """NLM prefixes retired DESCRIPTORS; the drug is not retired.

    "[OBSOLETE] avapritinib" reads as a withdrawn medicine. Avapritinib is
    approved, so displaying the prefix states something false about the drug.
    """
    bad = [r["drug"] for r in rows() if r["drug"].startswith("[OBSOLETE]")]
    assert not bad, (
        f"{len(bad)} rows show NLM's retired-descriptor prefix as part of the "
        f"drug name, e.g. {bad[0]!r}; that reads as a claim about the drug")


def test_rows_and_occurrences_are_not_the_same_measure():
    """They were reported as one number under two labels, and the second was
    4,300 short: 4,300 rows carry a variant on BOTH sides."""
    r = d()
    assert r["variant_entity_occurrences"] > r["relation_rows_touching_a_variant"], (
        "variant occurrences no longer exceed variant-touching rows, so either "
        "no relation has a variant on both sides or the two are being counted "
        "by the same expression again")


def test_the_map_never_fuses_genes_under_a_gene_less_variant():
    """`p.G12C` with no CorrespondingGene could be KRAS, NRAS or HRAS.

    The document says these are left unresolved. An earlier count keyed them all
    on ("", spelling), which silently made one variant out of three genes'.
    """
    r = d()
    assert r["variants_in_the_map_with_no_gene"] > 0, (
        "no gene-less variants remain in the map; if they were dropped rather "
        "than kept separate the document's account of them is wrong")
    # The reported count must BE the shipped table's, or the summary can claim
    # a separation the table does not have.
    assert r["distinct_variants_in_the_map"] == len(
        {(x["gene_id"], x["variant"]) for x in rows()}), (
        "the reported distinct-variant count is not the number of distinct "
        "(gene, variant) keys in the shipped table, so the summary and the "
        "table disagree about how variants were keyed")
    gene_less = {x["variant"] for x in rows() if not x["gene_id"]}
    with_gene = {x["variant"] for x in rows() if x["gene_id"]}
    assert gene_less & with_gene, (
        "no spelling appears both with and without a gene, so this guard is "
        "vacuous on the current corpus and the fusion risk it checks is gone")
    # The real property: those shared spellings are DISTINCT keys, not merged.
    for v in sorted(gene_less & with_gene)[:5]:
        genes = {x["gene_id"] for x in rows() if x["variant"] == v}
        assert len(genes) > 1 and "" in genes, (
            f"{v} appears under one gene key only, so the gene-less rows were "
            "folded into a gene rather than kept separate")
