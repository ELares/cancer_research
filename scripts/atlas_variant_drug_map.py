#!/usr/bin/env python3
"""Which drugs the cancer literature discusses alongside which mutations.

WHY THIS DATA WAS SITTING UNREAD
--------------------------------
`corpus/atlas/relations/relations.tsv.gz` carries entity types this project
never ingested: `ProteinMutation`, `DNAMutation`, `SNP`. Each variant entity
arrives annotated inline -- `RS#:121434568;HGVS:p.L858R;CorrespondingGene:1956`
is EGFR L858R -- and 228,501 relation rows touch one.

The only code that ever saw them discards them. `scripts/atlas_discovery.py:145`
skips any identifier beginning `RS#`, `HGVS` or `CorrespondingGene`, for a good
local reason: that function asks PubMed whether two terms co-occur in title or
abstract, and an HGVS string is not a searchable term. But the effect is that
the variant layer is absent from every other analysis in this repository.

WHAT THIS BUILDS
----------------
The drug-by-variant slice: relations with a `Chemical` on one side and a variant
on the other, resolved to readable names through the authority table already
committed for the co-mention work (`analysis/comention/authority-labels.tsv.gz`,
NLM for MeSH and NCBI for genes). Resistance mutations are the therapeutically
interesting case and surface by construction: a drug discussed alongside a
specific substitution is usually reported to work, or to stop working, against
it.

WHAT THIS IS NOT, STATED FIRST BECAUSE IT DECIDES HOW TO READ EVERYTHING ELSE
-----------------------------------------------------------------------------
CIViC, OncoKB and COSMIC curate drug-variant relationships with clinical
evidence levels, expert review and directionality. This has none of those. It
measures ATTENTION, not clinical actionability. `associate` carries no
direction, so a paper reporting that a mutation CONFERS RESISTANCE and one
reporting that it PREDICTS RESPONSE land in the same bucket. For any clinical
question the curated databases are the correct source.

THE IDENTIFIER IS NOT THE VARIANT, WHICH IS WHY THIS SCRIPT HAS TWO PASSES
--------------------------------------------------------------------------
The first version keyed a pair on whatever the entity string carried, and that
splits one variant across several keys. Two distinct defects, both measured
here rather than assumed, and both correctable without guessing:

1. THE RSID IS NOT THE VARIANT EITHER. The same change appears as
   `RS#+HGVS+CorrespondingGene`, as `HGVS+CorrespondingGene`, as
   `RS#+CorrespondingGene` with no HGVS, and at both protein and coding level
   (`p.T790M` and `c.2369C>T` are one variant under rs121434569).

   But an rsid CANNOT simply be collapsed to one HGVS: rs121913529 covers KRAS
   `p.G12D`, `p.G12V` AND `p.G12A` -- one multi-allelic codon-12 site, three
   different substitutions, targeted by different drugs. Merging on rsid would
   destroy exactly the distinction this map exists to show.

   The rule that separates them needs no threshold: parse each protein form to
   its (origin, position, new residue) triple, and collapse an rsid's forms
   only when every parsable one agrees. That merges representation variety and
   refuses multi-allelic sites. Measured over 2,983 rsids mapping to several
   HGVS, 1,857 are one substitution written several ways and 620 are genuinely
   multi-allelic. BRAF rs113488022 stays split, correctly: alongside 12,488
   `p.V600E` it carries `p.V600G` and `p.E600V`, which the rule cannot tell
   apart from a real second allele.

2. A DIGIT DROPPED FROM THE POSITION. `p.V61F` on gene 3717 (JAK2) is the
   canonical MPN driver V617F with a digit missing, and it carries 86% of that
   variant's volume -- so the shipped map read JAK2 V617F at a seventh of its
   real support. THE COUNTS CANNOT ADJUDICATE THIS and the first attempt to let
   them was wrong in both directions: TP53 `p.R72P` and ERBB2 `p.I655V` are
   real rs-backed polymorphisms whose longer twins are the typos, while EGFR
   `p.T790M` is the real one and `p.T90M` the typo.

   THE RSID ADJUDICATES. dbSNP assignment is independent of the string that
   was extracted, so when exactly one side of a digit-drop twin carries rsids
   that side is canonical -- and where neither or both do, this abstains rather
   than guessing, the same choice `atlas_graph.resolve` makes for a contested
   surface form. On JAK2 the verdict is independent of any count: `p.V617F`
   carries rs77375493 on 620 of 624 rows, `p.V61F` carries none at all, and
   both show the identical Polycythemia Vera / Myelofibrosis / Essential
   Thrombocythemia profile in the same rank order.

Usage:
    python scripts/atlas_variant_drug_map.py
    python scripts/atlas_variant_drug_map.py --gene EGFR
"""

import argparse
import collections
import gzip
import io
import itertools
import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RELATIONS = PROJECT_ROOT / "corpus" / "atlas" / "relations" / "relations.tsv.gz"
LABELS = PROJECT_ROOT / "analysis" / "comention" / "authority-labels.tsv.gz"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-variant-drug-map.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-variant-drug-map.json"
OUT_TSV = PROJECT_ROOT / "analysis" / "atlas-variant-drug-map.tsv.gz"

VARIANT_TYPES = ("ProteinMutation", "DNAMutation", "SNP")

# A single-substitution protein change: origin residue, position, new residue.
# Only this shape can be checked for a dropped position digit; indels, frame
# shifts and coding-DNA changes are left alone.
SUBSTITUTION = re.compile(r"^p\.([A-Z])(\d+)([A-Z*])$")

# NLM prefixes a retired descriptor's label. Displaying it unchanged reads as a
# statement about the DRUG -- "[OBSOLETE] avapritinib" suggests a withdrawn
# medicine, and avapritinib is approved. The descriptor was retired, not the
# drug, so the prefix is stripped for display and counted.
OBSOLETE = "[OBSOLETE] "


def load_labels() -> dict:
    """identifier -> primary authority name (first field before the pipe)."""
    lab = {}
    if not LABELS.exists():
        return lab
    with gzip.open(LABELS, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                lab[parts[0]] = parts[1].split("|")[0]
    return lab


def parse_variant(ident: str) -> dict:
    """Split a variant identifier into its `KEY:value;` fields."""
    return dict(f.split(":", 1) for f in ident.split(";") if ":" in f)


def iter_variants(path: Path):
    """Yield (pmid, predicate, variant_fields, other_type, other_id) per row."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, pred, a, b = parts[0], parts[1], parts[2], parts[3]
            ta, _, ia = a.partition("|")
            tb, _, ib = b.partition("|")
            a_var, b_var = ta in VARIANT_TYPES, tb in VARIANT_TYPES
            if not (a_var or b_var):
                continue
            var, other = ((ia, (tb, ib)) if a_var else (ib, (ta, ia)))
            yield pmid, pred, parse_variant(var), other[0], other[1]


def digit_drop_twins(gene_hgvs: dict) -> list:
    """Every (gene, shorter, longer) pair differing by one deleted position digit.

    Returns the full set with each side's occurrence count and rsid support, so
    the adjudication and the abstentions are both reportable.
    """
    out = []
    for gene, forms in gene_hgvs.items():
        parsed = {}
        for h, (n, rsids) in forms.items():
            m = SUBSTITUTION.match(h)
            if m:
                parsed[h] = (m.group(1), m.group(2), m.group(3), n, rsids)
        for x, y in itertools.combinations(sorted(parsed), 2):
            px, py = parsed[x], parsed[y]
            if px[0] != py[0] or px[2] != py[2]:
                continue
            (short, ps), (long_, pl) = sorted(((x, px), (y, py)),
                                              key=lambda t: len(t[1][1]))
            if len(pl[1]) != len(ps[1]) + 1:
                continue
            if not any(pl[1][:k] + pl[1][k + 1:] == ps[1] for k in range(len(pl[1]))):
                continue
            out.append({
                "gene": gene,
                "shorter": short, "shorter_rows": ps[3], "shorter_rsids": len(ps[4]),
                "longer": long_, "longer_rows": pl[3], "longer_rsids": len(pl[4]),
            })
    return out


def resolve_rsids(rs_hgvs: dict) -> tuple:
    """(rsid, gene) -> the one HGVS its forms all denote, where they do.

    Collapses representation variety (`p.T790M` and `c.2369C>T`) and REFUSES
    multi-allelic sites (rs121913529 = KRAS G12D and G12V and G12A). The test
    is that every parsable protein form agrees on its substitution triple, so
    nothing here rests on a dominance threshold.
    """
    fix, tally = {}, collections.Counter()
    for key, forms in rs_hgvs.items():
        if len(forms) == 1:
            fix[key] = next(iter(forms))
            tally["single_form"] += 1
            continue
        triples = collections.defaultdict(collections.Counter)
        for h, n in forms.items():
            m = SUBSTITUTION.match(h)
            if m:
                triples[m.groups()][h] += n
        if len(triples) == 1:
            # Every protein form denotes one change; the cDNA spellings under
            # the same rsid denote it too. Canonical = the commonest protein
            # form, so the key a reader recognises wins.
            fix[key] = next(iter(triples.values())).most_common(1)[0][0]
            tally["one_change_many_spellings"] += 1
        elif len(triples) > 1:
            tally["multi_allelic_refused"] += 1
        else:
            tally["no_parsable_substitution_refused"] += 1
    return fix, dict(tally)


def build_canonical_maps(gene_hgvs: dict, rs_hgvs: dict) -> tuple:
    """(twin corrections, rsid->HGVS resolutions, twin table, rsid tally).

    The twin rule fires only when exactly one side carries rsids. Corrections
    chain -- p.V2617F and p.V61F both reach p.V617F -- so each is followed to a
    fixed point.
    """
    twins = digit_drop_twins(gene_hgvs)
    raw = {}
    for t in twins:
        has_s, has_l = t["shorter_rsids"] > 0, t["longer_rsids"] > 0
        if has_s == has_l:
            t["verdict"] = "abstain"
            continue
        canon, wrong = ((t["shorter"], t["longer"]) if has_s
                        else (t["longer"], t["shorter"]))
        t["verdict"] = "canonical:" + canon
        t["rows_moved"] = t["longer_rows"] if has_s else t["shorter_rows"]
        # The rule is load-bearing exactly when it moves the MAJORITY onto the
        # minority string. Flagged rather than suppressed by a threshold: a
        # hand-picked cutoff would hide the cases worth a reader's attention.
        t["moves_majority"] = t["rows_moved"] > (
            t["shorter_rows"] + t["longer_rows"] - t["rows_moved"])
        raw[(t["gene"], wrong)] = canon

    fixed = {}
    for (gene, wrong), canon in raw.items():
        seen = {wrong}
        while (gene, canon) in raw and canon not in seen:
            seen.add(canon)
            canon = raw[(gene, canon)]
        fixed[(gene, wrong)] = canon

    rs_fix, rs_tally = resolve_rsids(rs_hgvs)
    return fixed, rs_fix, twins, rs_tally


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default=None,
                    help="narrow the report's TOP-PAIRS view to one gene symbol")
    ap.add_argument("--min-papers", type=int, default=3)
    args = ap.parse_args()

    lab = load_labels()
    if not lab:
        print(f"no authority labels at {LABELS}; names will be identifiers only")

    # --- pass 1: learn what the identifiers mean -------------------------------
    gene_hgvs = collections.defaultdict(
        lambda: collections.defaultdict(lambda: [0, collections.Counter()]))
    rs_hgvs = collections.defaultdict(collections.Counter)
    variant_rows = no_gene = 0
    for _, _, f, _, _ in iter_variants(RELATIONS):
        variant_rows += 1
        gene, hgvs, rs = (f.get("CorrespondingGene", ""), f.get("HGVS", ""),
                          f.get("RS#", ""))
        if not gene:
            no_gene += 1
        if gene and hgvs:
            rec = gene_hgvs[gene][hgvs]
            rec[0] += 1
            if rs:
                rec[1][rs] += 1
            if rs:
                rs_hgvs[(rs, gene)][hgvs] += 1

    twin_fix, rs_fix, twins, rs_tally = build_canonical_maps(gene_hgvs, rs_hgvs)

    # --- pass 2: the drug-by-variant table, on canonical keys ------------------
    pairs = collections.defaultdict(
        lambda: {"pmids": set(), "preds": collections.Counter()})
    rows = chem_variant = rs_filled = rs_respelled = twin_corrected = 0
    variants = set()
    genes = collections.Counter()

    for pmid, pred, f, other_t, other_i in iter_variants(RELATIONS):
        rows += 1
        gene, hgvs, rs = (f.get("CorrespondingGene", ""), f.get("HGVS", ""),
                          f.get("RS#", ""))
        # The rsid resolution applies to EVERY row carrying it, not only the
        # ones missing an HGVS: `c.2369C>T` and `p.T790M` are the same change
        # and were being counted as two variants.
        if rs and (rs, gene) in rs_fix:
            canon = rs_fix[(rs, gene)]
            if not hgvs:
                rs_filled += 1
            elif hgvs != canon:
                rs_respelled += 1
            hgvs = canon
        if hgvs and (gene, hgvs) in twin_fix:
            hgvs = twin_fix[(gene, hgvs)]
            twin_corrected += 1
        if gene:
            genes[gene] += 1
        variants.add((gene, hgvs or rs))
        if other_t != "Chemical":
            continue
        chem_variant += 1
        key = (other_i, gene, hgvs or rs)
        pairs[key]["pmids"].add(pmid)
        pairs[key]["preds"][pred] += 1

    def name(ident, prefix=""):
        n = lab.get(ident)
        if n is None:
            return f"{prefix}{ident}" if ident else "?"
        return n[len(OBSOLETE):] if n.startswith(OBSOLETE) else n

    table = []
    for (chem, gene, var), v in pairs.items():
        table.append({
            "drug": name(chem), "drug_id": chem,
            "gene": name(gene, "gene:"), "gene_id": gene,
            "variant": var,
            "papers": len(v["pmids"]),
            "assertions": sum(v["preds"].values()),
            "predicates": dict(v["preds"].most_common()),
        })
    table.sort(key=lambda r: (-r["papers"], r["drug"], r["gene"], r["variant"]))

    # mtime=0 and no embedded filename, so an unchanged input regenerates a
    # BYTE-IDENTICAL artifact. gzip.open() stamps the current time into the
    # header, which makes every rerun a spurious diff and defeats any
    # byte-identity check on a committed file.
    with open(OUT_TSV, "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as fh:
        fh.write("drug\tdrug_id\tgene\tgene_id\tvariant\tpapers\tassertions\tpredicates\n")
        for r in table:
            fh.write(f"{r['drug']}\t{r['drug_id']}\t{r['gene']}\t{r['gene_id']}\t"
                     f"{r['variant']}\t{r['papers']}\t{r['assertions']}\t"
                     f"{';'.join(f'{k}={n}' for k, n in r['predicates'].items())}\n")

    view = [r for r in table if r["gene"] == args.gene] if args.gene else table
    applied = [t for t in twins if t["verdict"] != "abstain"]
    res = {
        "relation_rows_touching_a_variant": rows,
        "variant_entity_occurrences": variant_rows,
        "occurrences_with_no_gene": no_gene,
        "chemical_variant_rows": chem_variant,
        "rsid_rows_given_an_hgvs": rs_filled,
        "rsid_rows_respelled_to_one_form": rs_respelled,
        "rsid_resolution": rs_tally,
        "digit_drop_rows_corrected": twin_corrected,
        "distinct_variants": len(variants),
        "distinct_drug_variant_pairs": len(table),
        "pairs_at_min_papers": sum(1 for r in table if r["papers"] >= args.min_papers),
        "single_paper_pairs": sum(1 for r in table if r["papers"] == 1),
        "twins": {
            "found": len(twins),
            "adjudicated": len(applied),
            "abstained": len(twins) - len(applied),
            "moving_the_majority": [t for t in applied if t.get("moves_majority")],
            "largest": sorted(twins, key=lambda t: -(t["shorter_rows"] + t["longer_rows"]))[:12],
        },
        "top_genes_by_variant_rows": {name(g, "gene:"): n for g, n in genes.most_common(15)},
        "top_pairs": view[:60],
        "gene_filter": args.gene,
        "min_papers": args.min_papers,
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(res), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}\nwrote {OUT_TSV}")
    print(f"  {chem_variant:,} chemical-variant rows -> {len(table):,} pairs")
    print(f"  canonicalized: {rs_filled:,} rsid rows given an HGVS, "
          f"{rs_respelled:,} respelled, {twin_corrected:,} digit-drop "
          f"({len(applied)} of {len(twins)} twins)")
    return 0


def render(r: dict) -> str:
    tw = r["twins"]
    rsr = r["rsid_resolution"]
    single_pct = 100.0 * r["single_paper_pairs"] / max(r["distinct_drug_variant_pairs"], 1)
    nogene_pct = 100.0 * r["occurrences_with_no_gene"] / max(r["variant_entity_occurrences"], 1)
    L = [
        "# Which drugs the literature discusses alongside which mutations", "",
        "Generated by `scripts/atlas_variant_drug_map.py`.", "",
        "## Read this first", "",
        "**CIViC, OncoKB and COSMIC are the correct sources for any clinical",
        "question about a drug and a variant.** They curate directionality,",
        "evidence level and expert review. This has none of that. It is",
        "machine-extracted co-assertion over the whole cancer census, and it",
        "measures **attention, not actionability**.", "",
        "In particular the `associate` predicate carries **no direction**: a paper",
        "reporting that a mutation confers RESISTANCE and one reporting that it",
        "predicts RESPONSE land in the same bucket, and nothing here separates",
        "them. A high paper count means the pair is discussed, not that the drug",
        "works.", "",
        "## Why this data was unread", "",
        "`relations.tsv.gz` carries `ProteinMutation`, `DNAMutation` and `SNP`",
        "entities annotated inline: `RS#:121434568;HGVS:p.L858R;"
        "CorrespondingGene:1956` is EGFR L858R. The only code that ever saw them",
        "discards them -- `atlas_discovery.py:145` skips any identifier beginning",
        "`RS#`, `HGVS` or `CorrespondingGene`, for the good local reason that it",
        "asks PubMed whether two terms co-occur and an HGVS string is not a",
        "searchable term. The effect was that the variant layer was absent from",
        "every other analysis here.", "",
        "## The identifier is not the variant", "",
        "Keying a pair on whatever string the entity carried splits one variant",
        "across several keys. Two defects, both measured rather than assumed:", "",
        "**The rsid is not the variant either.** One change arrives with an",
        "rsid and without, and at protein and coding level -- `p.T790M` and",
        "`c.2369C>T` are one variant under rs121434569. But an rsid cannot",
        "simply be collapsed to one HGVS: **rs121913529 covers KRAS `p.G12D`,",
        "`p.G12V` and `p.G12A`** -- one multi-allelic codon-12 site, three",
        "substitutions, different drugs. Merging on rsid would destroy exactly",
        "the distinction this map exists to show.", "",
        "The rule that separates them needs no threshold: parse each protein",
        "form to its (origin, position, new residue) triple and collapse an",
        "rsid's forms only when every parsable one agrees.", "",
        "| rsIDs mapping to several HGVS | |", "|---|--:|",
        f"| one substitution, several spellings: collapsed | "
        f"{rsr.get('one_change_many_spellings', 0):,} |",
        f"| genuinely multi-allelic: refused | "
        f"{rsr.get('multi_allelic_refused', 0):,} |",
        f"| no parsable substitution: refused | "
        f"{rsr.get('no_parsable_substitution_refused', 0):,} |", "",
        f"That gave **{r['rsid_rows_given_an_hgvs']:,}** rows an HGVS they",
        f"lacked and respelled **{r['rsid_rows_respelled_to_one_form']:,}** onto",
        "a single form. BRAF rs113488022 stays split, correctly: alongside",
        "12,488 `p.V600E` it carries `p.V600G` and `p.E600V`, which the rule",
        "cannot distinguish from a real second allele.", "",
        f"**A dropped position digit.** {tw['found']} twin pairs differ by one",
        "deleted digit. Counts cannot adjudicate them and reading the majority as",
        "correct is wrong in both directions: TP53 `p.R72P` and ERBB2 `p.I655V`",
        "are real rs-backed polymorphisms whose longer twins are typos, while",
        "EGFR `p.T790M` is real and `p.T90M` the typo.", "",
        "**The rsid adjudicates**, because dbSNP assignment is independent of the",
        "string that was extracted. Where exactly one side carries rsids that side",
        f"is canonical ({tw['adjudicated']} pairs, "
        f"{r['digit_drop_rows_corrected']:,} rows); where neither or both do this",
        f"abstains ({tw['abstained']} pairs) rather than guessing, the same choice",
        "`atlas_graph.resolve` makes for a contested surface form.", "",
    ]
    if tw["moving_the_majority"]:
        L += ["The rule is load-bearing exactly where it moves the MAJORITY onto",
              "the minority string. Those cases are listed rather than suppressed",
              "by a threshold, because a hand-picked cutoff would hide precisely",
              "the rows worth checking:", "",
              "| gene | non-canonical | rows | canonical | rows | rsids |",
              "|---|---|--:|---|--:|--:|"]
        for t in tw["moving_the_majority"]:
            s_can = t["verdict"] == "canonical:" + t["shorter"]
            bad, good = ((t["longer"], t["shorter"]) if s_can
                         else (t["shorter"], t["longer"]))
            badn, goodn = ((t["longer_rows"], t["shorter_rows"]) if s_can
                           else (t["shorter_rows"], t["longer_rows"]))
            rsn = t["shorter_rsids"] if s_can else t["longer_rsids"]
            L.append(f"| {t['gene']} | `{bad}` | {badn:,} | `{good}` | {goodn:,} | {rsn} |")
        L += ["",
              "The JAK2 row is the one that matters. Gene 3717 `p.V61F` is the",
              "canonical myeloproliferative-neoplasm driver V617F with a digit",
              "missing, and it carried 86% of that variant's volume, so the",
              "uncorrected map read JAK2 V617F at a seventh of its real support.",
              "The verdict does not rest on counts: `p.V617F` carries rs77375493",
              "on 620 of 624 rows and `p.V61F` carries no rsid at all, and both",
              "show the identical Polycythemia Vera / Primary Myelofibrosis /",
              "Essential Thrombocythemia profile in the same rank order.", ""]
    L += ["## What is in it", "",
          "| | count |", "|---|--:|",
          f"| relation rows touching a variant | {r['relation_rows_touching_a_variant']:,} |",
          f"| variant entity occurrences | {r['variant_entity_occurrences']:,} |",
          f"| ...carrying no gene, left unresolved | {r['occurrences_with_no_gene']:,} ({nogene_pct:.1f}%) |",
          f"| **chemical-to-variant rows** | **{r['chemical_variant_rows']:,}** |",
          f"| distinct variants | {r['distinct_variants']:,} |",
          f"| distinct (drug, gene, variant) pairs | {r['distinct_drug_variant_pairs']:,} |",
          f"| pairs with >= {r['min_papers']} papers | {r['pairs_at_min_papers']:,} |",
          f"| pairs resting on ONE paper | {r['single_paper_pairs']:,} ({single_pct:.1f}%) |",
          "",
          f"An HGVS string with no `CorrespondingGene` ({nogene_pct:.1f}% of",
          "occurrences) is genuinely ambiguous -- `p.G12C` alone could be KRAS,",
          "NRAS or HRAS -- so those are left unresolved rather than assigned to",
          "the most likely gene.", "",
          f"**{single_pct:.0f}% of pairs rest on a single paper.** The retraction",
          "analysis found the same shape across the whole graph (70.2%), and it",
          "is the first thing to know before reading any row below as evidence.", "",
          "## Genes carrying the most variant relations", "",
          "| gene | variant-touching rows |", "|---|--:|"]
    for g, n in list(r["top_genes_by_variant_rows"].items())[:15]:
        L.append(f"| {g} | {n:,} |")
    L += ["", "## The most-discussed drug-variant pairs", "",
          "Ranked by asserting papers. `papers` counts distinct PMIDs; the",
          "predicate counts are ASSERTIONS, so they exceed the paper count when",
          "one paper states a pair more than once. Again: this ranks how much a",
          "pair is *written about*, and `associate` does not say which direction.", "",
          "| drug | gene | variant | papers | predicates (assertions) |",
          "|---|---|---|--:|---|"]
    for row in r["top_pairs"][:40]:
        preds = ", ".join(f"`{k}` {n}" for k, n in list(row["predicates"].items())[:3])
        L.append(f"| {row['drug']} | {row['gene']} | `{row['variant']}` | "
                 f"{row['papers']} | {preds} |")
    L += ["", "## What this cannot say", "",
          "* **No direction.** `associate` covers resistance and sensitivity alike.",
          "* **Not curated.** No evidence level, no expert review, no clinical",
          "  interpretation. Use CIViC/OncoKB/COSMIC for those.",
          "* **Extractor error dominates.** PubTator's own error rate is larger",
          "  than most differences between rows here, and the digit-drop finding",
          "  above is a direct instance of it.",
          "* **Canonicalization is incomplete.** It corrects the twins an rsid can",
          f"  adjudicate; the {tw['abstained']} abstentions stay fragmented, as do",
          "  variants written in forms this does not compare (indels, frame",
          "  shifts, coding-DNA changes).",
          "* **Attention is not importance.** A well-studied pair outranks a real",
          "  but rarely-written-about one, exactly as `atlas_model_gaps.py` warns",
          "  for its own ranking.",
          "* **Full table** is `analysis/atlas-variant-drug-map.tsv.gz`; the rows",
          "  above are a view, not the result.", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
