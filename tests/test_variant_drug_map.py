"""Guards for the drug-by-variant map.

WHAT IT MEASURES
----------------
`relations.tsv.gz` carries `ProteinMutation`/`DNAMutation`/`SNP` entities that
no analysis in this repository had ever read -- `atlas_discovery.py:145` skips
them by design. The map is the drug-by-variant slice of them, resolved to
readable names.

THE THREE THINGS THAT CAN GO WRONG, and why each gets a guard
--------------------------------------------------------------
1. IT READS AS A CLINICAL RESOURCE. It is machine-extracted co-assertion with
   no directionality: `associate` covers "confers resistance" and "predicts
   response" identically. A reader who takes a row as clinical guidance is the
   most probable way this document causes harm, so the caveat has to reach them
   before the table does.

2. CANONICALIZATION OVER-MERGES. Collapsing an rsid to one HGVS is the obvious
   fix for identifier fragmentation and it is WRONG: rs121913529 covers KRAS
   G12D, G12V and G12A -- one multi-allelic site, three substitutions, three
   different drug programs. Over-merging would destroy exactly the distinction
   the map exists to show, and would do it silently, since the merged row still
   looks like a plausible result.

3. CANONICALIZATION UNDER-MERGES. The un-canonicalized map read JAK2 V617F at a
   seventh of its support, because PubTator emits `p.V61F` -- the canonical MPN
   driver with a position digit dropped -- for 86% of that variant's rows.

WHY THE JAK2 GUARD KEYS ON THE RSID AND NOT ON THE COUNTS
---------------------------------------------------------
Counts cannot adjudicate a digit-drop twin and reading the majority as correct
is wrong in both directions: TP53 `p.R72P` and ERBB2 `p.I655V` are real rs-backed
polymorphisms whose longer twins are the typos. The rsid is independent evidence,
so these check that the report gives it, not merely that it reached the right
answer.
"""

import gzip
import importlib.util
import json
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
    first_row = txt.index("| osimertinib")
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


# --- 2. canonicalization must not over-merge -------------------------------

def test_the_multi_allelic_refusal_is_still_protecting_something():
    """If this reaches zero the rule has stopped distinguishing anything."""
    t = d()["rsid_resolution"]
    assert t.get("multi_allelic_refused", 0) > 0, (
        "no rsID is being refused as multi-allelic any more, so either the "
        "corpus changed shape or the rule collapsed into an unconditional "
        "merge -- which would silently fuse KRAS G12C with G12D")
    assert t.get("one_change_many_spellings", 0) > 0, (
        "no rsID is being collapsed either, so the rule is now inert and the "
        "map is back to counting `p.T790M` and `c.2369C>T` as two variants")


def test_kras_codon_twelve_substitutions_stay_distinct():
    """The load-bearing safety property, checked on the shipped table.

    rs121913529 covers G12D, G12V and G12A; rs121913530 covers G12C, G12R and
    G12S. Any rule that keyed on the rsID would fuse them, and the fused row
    would still look entirely plausible.
    """
    kras = {r["variant"] for r in rows() if r["gene"] == "KRAS"}
    for v in ("p.G12C", "p.G12D", "p.G12V"):
        assert v in kras, (
            f"KRAS {v} is absent from the map; it shares an rsID with the other "
            "codon-12 substitutions, so its disappearance means they were "
            "merged")
    # And they must carry DIFFERENT drugs, which is what the distinction is for.
    by_var = {}
    for r in rows():
        if r["gene"] == "KRAS":
            by_var.setdefault(r["variant"], set()).add(r["drug"])
    assert by_var["p.G12C"] != by_var["p.G12D"], (
        "KRAS G12C and G12D now carry an identical drug set, which is what a "
        "merge on the shared rsID would produce")


def test_the_resolver_refuses_a_multi_allelic_site_and_collapses_a_respelling():
    """Unit test of the rule itself, on synthetic input.

    Built from cases the corpus supplies, but exercised directly so the rule is
    pinned independently of whatever the current dump happens to contain.
    """
    m = mod()
    import collections
    fix, tally = m.resolve_rsids({
        # one change, protein and coding spellings -> collapse to the protein form
        ("121434569", "1956"): collections.Counter({"p.T790M": 4146, "c.2369C>T": 9}),
        # three substitutions at one codon -> refuse
        ("121913529", "3845"): collections.Counter({"p.G12D": 1830, "p.G12V": 873,
                                                    "p.G12A": 126}),
        # a lone form -> pass through unchanged
        ("77375493", "3717"): collections.Counter({"p.V617F": 620}),
    })
    assert fix[("121434569", "1956")] == "p.T790M", (
        "the coding and protein spellings of one change did not collapse onto "
        "the protein form")
    assert ("121913529", "3845") not in fix, (
        "a multi-allelic rsID was resolved to a single HGVS; this is the merge "
        "that would fuse KRAS G12D with G12V")
    assert fix[("77375493", "3717")] == "p.V617F"
    assert tally["multi_allelic_refused"] == 1
    assert tally["one_change_many_spellings"] == 1


# --- 3. canonicalization must not under-merge ------------------------------

def test_the_jak2_correction_is_applied():
    m = mod()
    import collections
    gene_hgvs = {"3717": {
        "p.V61F": [3747, collections.Counter()],
        "p.V617F": [622, collections.Counter({"77375493": 620})],
    }}
    fix, _, twins, _ = m.build_canonical_maps(gene_hgvs, {})
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
    fix, _, twins, _ = m.build_canonical_maps(gene_hgvs, {})
    assert fix == {}, (
        "a twin was corrected although the rsIDs cannot adjudicate it; the "
        "rule would then be deciding on counts, which is wrong in both "
        "directions (TP53 p.R72P is the short form and correct, EGFR p.T790M "
        "is the long form and correct)")
    assert {t["verdict"] for t in twins} == {"abstain"}


def test_the_jak2_verdict_is_justified_by_the_rsid_not_by_the_counts():
    """The report must give the independent evidence, not just the answer."""
    txt = flat()
    assert "rs77375493" in txt, (
        "the report no longer cites the rsID that adjudicates JAK2 V617F, so "
        "it is asserting the correction rather than evidencing it")
    i = txt.index("rs77375493")
    window = txt[max(0, i - 600): i + 600]
    assert "no rsid at all" in window.lower(), (
        "the report cites the canonical form's rsID without stating that the "
        "malformed twin carries none, which is the half that makes it evidence")
    assert "Polycythemia Vera" in window, (
        "the corroborating disease profile is no longer stated beside the "
        "rsID; the two together are what make the verdict independent of the "
        "3,747-to-622 count that points the other way")


def test_the_counts_that_would_mislead_are_present_so_the_reader_can_see_them():
    """The majority-moving table must show the losing side's volume."""
    movers = d()["twins"]["moving_the_majority"]
    jak2 = [t for t in movers if t["gene"] == "3717"]
    assert jak2, "the JAK2 twin is no longer reported as moving the majority"
    t = jak2[0]
    assert t["shorter_rows"] > t["longer_rows"], (
        "the malformed JAK2 string no longer outnumbers the canonical one; if "
        "the upstream extractor was fixed this whole correction should be "
        "revisited rather than kept")


# --- provenance -------------------------------------------------------------

def test_the_skipped_discovery_filter_is_still_the_reason_this_was_unread():
    """The claim that nothing else reads this layer must stay checkable."""
    src = (REPO_ROOT / "scripts" / "atlas_discovery.py").read_text()
    assert 'startswith(("RS#", "HGVS", "CorrespondingGene"))' in src, (
        "atlas_discovery.py no longer skips variant identifiers, so the claim "
        "that this layer was unread needs re-checking against whatever now "
        "consumes it")


def test_the_obsolete_prefix_is_not_shown_beside_an_approved_drug():
    """NLM prefixes retired DESCRIPTORS; the drug is not retired.

    "[OBSOLETE] avapritinib" reads as a withdrawn medicine. Avapritinib is
    approved, so displaying the prefix unchanged states something false about
    the drug rather than about the descriptor.
    """
    bad = [r["drug"] for r in rows() if r["drug"].startswith("[OBSOLETE]")]
    assert not bad, (
        f"{len(bad)} rows show NLM's retired-descriptor prefix as part of the "
        f"drug name, e.g. {bad[0]!r}; that reads as a claim about the drug")
