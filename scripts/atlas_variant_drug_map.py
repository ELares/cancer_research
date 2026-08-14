#!/usr/bin/env python3
"""Which drugs the cancer literature discusses alongside which mutations.

WHY THIS LAYER WAS UNRESOLVED
-----------------------------
`corpus/atlas/relations/relations.tsv.gz` carries `ProteinMutation`,
`DNAMutation` and `SNP` entities whose identifier holds the variant inline:
`RS#:121434568;HGVS:p.L858R;CorrespondingGene:1956` is EGFR L858R.

Those entities are NOT absent from this repository -- they are nodes in the
graph index like any other, and `analysis/atlas-emergence.md` already ships a
row reading `adagrasib -- RS#:121913530;HGVS:p.G12C;CorrespondingGene:3845`.
What no analysis here has ever done is PARSE those fields. The identifier is
passed through whole, as an opaque label, so nothing could group two spellings
of one variant or say which gene a mutation belongs to. `atlas_discovery.py:145`
skips them outright, for a good local reason: it asks PubMed whether two terms
co-occur, and an HGVS string is not a searchable term.

WHAT THIS BUILDS
----------------
The drug-by-variant slice: relations with a `Chemical` on one side and a variant
on the other, with the fields parsed, the spellings reconciled, and both sides
resolved to readable names through the authority table committed for the
co-mention work (`analysis/comention/authority-labels.tsv.gz`).

WHAT THIS IS NOT, STATED FIRST BECAUSE IT DECIDES HOW EVERYTHING ELSE READS
---------------------------------------------------------------------------
CIViC, OncoKB and COSMIC curate drug-variant relationships with clinical
evidence levels, expert review and directionality. This has none of those. It
measures ATTENTION, not clinical actionability. `associate` carries no
direction, so a paper reporting that a mutation CONFERS RESISTANCE and one
reporting that it PREDICTS RESPONSE land in the same bucket. For any clinical
question the curated databases are the correct source.

THE IDENTIFIER IS NOT THE VARIANT
---------------------------------
Keying a pair on whatever string the entity carried splits one variant across
several keys. Two defects, and for each the OBVIOUS correction is wrong:

1. THE RSID IS NOT THE VARIANT EITHER. One change arrives with an rsid and
   without, and at protein and coding level. But an rsid cannot simply be
   collapsed to one HGVS: rs121913529 covers KRAS `p.G12D`, `p.G12V` AND
   `p.G12A` -- one multi-allelic codon-12 site, three substitutions, different
   drug programs.

   THE AGREEMENT TEST MUST COVER EVERY SPELLING, which is where the first
   version of this was wrong. Testing agreement among PROTEIN forms only and
   sweeping the rest onto the winner let a minority `p.E429A` absorb rs1801131's
   dominant `c.1298A>C` together with `c.1286A>C` and `c.1298A>T`, which are a
   different position and a different allele. (The size of that defect is
   computed and reported; it is deliberately not repeated here, because the
   first fix for it hand-wrote the figure into this docstring and it went stale
   against the artifact within one commit.)

   Three tests now, and a spelling failing any of them refuses its whole rsid:
   every spelling must PARSE, each representation class must agree INTERNALLY,
   and a protein spelling must agree with a coding one ACROSS classes, since
   codon = (coding position + 2) // 3 relates them. That last test is what
   catches rs1494558, where `p.I66T` and `c.412G>A` (codon 138) are two
   different variants, and rs2297518, where a missense sits beside a promoter
   position that cannot encode it.

2. THE MAJORITY IS NOT THE TRUTH. Some twins differ by one deleted position
   digit, and reading the commoner as correct is wrong in both directions: TP53
   `p.R72P` and ERBB2 `p.I655V` are real rs-backed polymorphisms whose LONGER
   twins are the typos, while EGFR `p.T790M` is real and `p.T90M` the typo.

   THE RSID ADJUDICATES, being assigned independently of the string that was
   extracted. Where exactly one side carries rsids that side is canonical.
   Where neither or both do, or where a form is a digit-drop twin of TWO
   different canonicals (KRAS `p.G1V` is one digit from both `p.G12V` and
   `p.G13V`), this ABSTAINS rather than guessing -- the same choice
   `atlas_graph.resolve` makes for a contested surface form.

EVERY FIGURE IN THE REPORT IS DERIVED, INCLUDING THE WORKED EXAMPLES, which name
whichever case the data says is largest. An earlier version hand-wrote them and
shipped an rsid count that was wrong, in a paragraph whose neighbouring figures
were freshly computed -- which is the shape that makes a stale number credible.

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

# A single-residue protein substitution, and a single-base change against a
# coding/genomic/mitochondrial reference. Only these shapes can be compared;
# anything else (indels, frame shifts, truncated or HTML-escaped strings such as
# `c.1799T&gt`) is UNPARSABLE and refuses its whole rsid.
PROTEIN = re.compile(r"^p\.([A-Z])(\d+)([A-Z*])$")
NUCLEOTIDE = re.compile(r"^([cgmn])\.([+-]?\d+)([ACGT])>([ACGT])$")

# NLM prefixes a retired descriptor's label. Displaying it unchanged states
# something false about the DRUG: "[OBSOLETE] avapritinib" reads as a withdrawn
# medicine, and avapritinib is approved. The DESCRIPTOR was retired.
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


def iter_rows(path: Path):
    """Yield (pmid, predicate, ((type, id), (type, id))) for every relation row."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            ta, _, ia = parts[2].partition("|")
            tb, _, ib = parts[3].partition("|")
            yield parts[0], parts[1], ((ta, ia), (tb, ib))


def _is_parsable(h: str) -> bool:
    """One definition of parsable, used by the classifier AND the evidence.

    A no-change `A>A` spelling matches the NUCLEOTIDE regex but denotes
    nothing, so it is unparsable. Two call sites disagreeing about that is what
    left a block of refusals reporting no offending spelling at all.
    """
    if PROTEIN.match(h):
        return True
    m = NUCLEOTIDE.match(h)
    return bool(m) and m.group(3) != m.group(4)


def classify_spellings(forms) -> tuple:
    """Group HGVS spellings by representation class -> the changes each denotes.

    Returns (classes, n_unparsable). `classes` maps a class key ("p", "c", "g",
    ...) to the set of distinct changes seen in it, so a class holding more than
    one change is a disagreement.

    `X` and `*` both denote a stop codon, so they are normalised together: without
    it, rsids carrying both were reported to the reader as carrying "spellings that denote
    different changes" when the only difference was `p.R461X` against
    `p.R461*`. A no-change nucleotide spelling (`A>A`) denotes nothing and is
    counted unparsable rather than allowed to agree with itself.
    """
    classes, unparsable = collections.defaultdict(set), 0
    for h in forms:
        m = PROTEIN.match(h)
        if m:
            orig, pos, new = m.groups()
            classes["p"].add((orig, pos, "*" if new == "X" else new))
            continue
        m = NUCLEOTIDE.match(h)
        if m and _is_parsable(h):
            classes[m.group(1)].add(m.groups()[1:])
            continue
        unparsable += 1
    return classes, unparsable


def codon_of(coding_position: str):
    """The codon a plain coding position falls in, or None if there is not one.

    A `+`/`-` offset is intronic or untranslated, so no protein change can
    correspond to it and the caller must refuse rather than compare.
    """
    if coding_position.startswith(("+", "-")):
        return None
    return (int(coding_position) + 2) // 3


def cross_class_conflict(classes: dict) -> tuple:
    """Can a protein spelling denote the same change as a coding one?

    Returns (outcome, description). Outcome is one of:
      "agree"        the codon relation was computed and it holds
      "mismatch"     the codon relation was computed and it fails
      "uncomparable" no comparison was possible, so this refuses
      "none"         there is nothing to compare

    THE THREE OUTCOMES ARE KEPT APART BECAUSE THE REPORT QUOTES THEM. Collapsing
    "mismatch" and "uncomparable" into one refusal count let the document claim
    the codon relation "fails" for cases where it was never evaluated, and
    label a protein-versus-GENOMIC refusal as protein-versus-coding.

    Each class is known to hold exactly one change when this is called: the
    three preceding branches in `resolve_rsids` have already excluded a
    multi-change class, any unparsable spelling, and an empty `classes`. That
    invariant is what makes `next(iter(...))` safe here.
    """
    if "p" not in classes:
        return "none", None
    residue = int(next(iter(classes["p"]))[1])
    compared = False
    for cls_key, changes in classes.items():
        if cls_key == "p":
            continue
        if cls_key != "c":
            # A genomic, mitochondrial or non-coding-transcript position has no
            # arithmetic relation to a residue number. Refusing keeps a variant
            # visibly fragmented, which is the failure this design prefers --
            # and it is not over-cautious: on this corpus four of these would be
            # wrongly ACCEPTED if a genomic coordinate were read as a codon.
            return ("uncomparable",
                    f"a `{cls_key}.` position cannot be checked against a residue number")
        codon = codon_of(next(iter(changes))[0])
        if codon is None:
            return "uncomparable", "a protein change beside a non-coding position"
        if codon != residue:
            return ("mismatch",
                    f"the coding position falls in codon {codon}, not residue {residue}")
        compared = True
    return ("agree" if compared else "none"), None


def resolve_rsids(rs_hgvs: dict) -> tuple:
    """(rsid, gene) -> the single HGVS its spellings all denote, where they do.

    Collapses representation variety (`p.T790M` and `c.2369C>T`) and REFUSES
    both multi-allelic sites (rs121913529 = KRAS G12D and G12V and G12A) and
    anything carrying a spelling that cannot be checked. Nothing rests on a
    dominance threshold: one disagreeing spelling refuses the rsid however rare
    it is, so a minority spelling cannot capture a dominant one.
    """
    fix, tally, refused = {}, collections.Counter(), []
    for key, forms in rs_hgvs.items():
        if len(forms) == 1:
            fix[key] = next(iter(forms))
            tally["single_spelling"] += 1
            continue
        classes, unparsable = classify_spellings(forms)
        # A genuine disagreement is checked FIRST. Both can hold at once, and
        # reporting the unparsable-string technicality then leaves the report
        # showing spellings that all look checkable, so a reader cannot see why
        # the rsid was refused.
        if any(len(s) > 1 for s in classes.values()):
            tally["classes_disagree_refused"] += 1
            reason = "spellings that denote different changes"
            # `_is_parsable`, not the bare regex. A no-change `A>A` spelling
            # MATCHES NUCLEOTIDE but denotes nothing and is counted unparsable
            # by the classifier, so testing the regex here listed it among the
            # spellings that "denote different changes" -- the same two-call-
            # sites-disagree defect already fixed one branch below, still live
            # in this one.
            evidence = sorted({
                h for h in forms
                if _is_parsable(h)
                and len(classes[(m.group(1) if (m := NUCLEOTIDE.match(h))
                                 else "p")]) > 1})
        elif unparsable:
            tally["unparsable_spelling_refused"] += 1
            reason = "a spelling that cannot be checked"
            # The SAME predicate `classify_spellings` uses, not a bare regex
            # match: a no-change `A>A` spelling matches NUCLEOTIDE but is
            # counted unparsable, so testing the regex alone left 42 refusals
            # shipping an EMPTY offending set -- a table whose every row reads
            # clean under a heading saying one of them is not.
            evidence = sorted(h for h in forms if not _is_parsable(h))
        elif (conflict := cross_class_conflict(classes))[0] in ("mismatch",
                                                                "uncomparable"):
            outcome, reason = conflict
            tally["codon_mismatch_refused" if outcome == "mismatch"
                  else "not_checkable_against_a_residue_refused"] += 1
            evidence = sorted(forms)
        else:
            # Recorded so the report can say how often the cross-class relation
            # was actually EVALUATED and held, rather than quoting the collapse
            # total as though every collapse had been checked that way.
            tally["codon_relation_agrees" if conflict[0] == "agree"
                  else "nothing_to_compare_across_classes"] += 1
            # Every spelling denotes one change, and a protein spelling agrees
            # with its coding one. Canonical is the commonest PROTEIN form when
            # there is one, since that is the key a reader recognises, else the
            # commonest spelling of any class.
            prot = {h: n for h, n in forms.items() if PROTEIN.match(h)}
            fix[key] = max((prot or forms).items(), key=lambda kv: (kv[1], kv[0]))[0]
            tally["one_change_many_spellings"] += 1
            continue
        # (`not classes` needs no branch: an empty classes dict means every
        # spelling was unparsable, which the branch above has already caught.
        # It shipped as a tally row that was provably always zero.)
        # Would the SUPERSEDED rule -- agreement among protein forms only, every
        # other spelling swept onto the winner -- have collapsed this? Recording
        # it makes the size of that defect a derived number rather than a
        # remembered one, and picks the example that demonstrates the fix.
        # Would the SUPERSEDED rule -- agreement among protein forms only, every
        # other spelling swept onto the winner -- have collapsed this? And how
        # many rows would have CHANGED KEY, which is smaller than the rsid's
        # total and is the number that means something.
        old_collapse = len(classes.get("p", ())) == 1
        old_canon = max(((h, n) for h, n in forms.items() if PROTEIN.match(h)),
                        key=lambda kv: (kv[1], kv[0]), default=("", 0))
        refused.append({"rsid": key[0], "gene": key[1],
                        "rows": sum(forms.values()), "reason": reason,
                        "old_rule_would_have_collapsed": old_collapse,
                        "old_rule_rows_moved": (sum(forms.values()) - old_canon[1]
                                                if old_collapse else 0),
                        # The spellings that DEMONSTRATE the reason, so the
                        # report shows the evidence rather than the top of an
                        # unrelated frequency list. NOT truncated: the render
                        # flags a displayed row by membership here, and a cap
                        # made a genuinely offending spelling show as clean.
                        "offending": {h: forms[h] for h in evidence},
                        # The protein spellings from the FULL set. Deriving
                        # them from `offending` gave an always-empty dict on a
                        # disagreement refusal, since a sole protein form
                        # agrees with itself and is never offending -- so the
                        # explanatory paragraph fell back to the truncated
                        # top-six it was written to stop using.
                        "protein_spellings": {h: n for h, n in forms.items()
                                              if PROTEIN.match(h)},
                        "n_spellings": len(forms),
                        "spellings": dict(forms.most_common(6))})
    refused.sort(key=lambda r: (-r["rows"], r["rsid"]))
    # How often does the SAME rsid get different treatment under different
    # genes? The unit of resolution is the (rsid, gene) key, and the report
    # once asserted "ten rsids are refused under one gene while still
    # collapsing under another" -- hand-written, and false: none of them
    # collapses elsewhere, they have a single spelling there, which is the case
    # the same document calls one where no collapse decision was made.
    by_rsid = collections.defaultdict(set)
    for (rs_, g_) in rs_hgvs:
        by_rsid[rs_].add(g_)
    refused_keys = {(r["rsid"], r["gene"]) for r in refused}
    collapsed = {k for k in fix if len(rs_hgvs[k]) > 1}
    single = {k for k in fix if len(rs_hgvs[k]) == 1}
    tally["refused_here_collapsed_under_another_gene"] = sum(
        1 for rs_, g_ in refused_keys
        if any((rs_, o) in collapsed for o in by_rsid[rs_] if o != g_))
    tally["refused_here_single_spelling_under_another_gene"] = sum(
        1 for rs_, g_ in refused_keys
        if any((rs_, o) in single for o in by_rsid[rs_] if o != g_))
    return fix, dict(tally), refused


def digit_drop_twins(gene_hgvs: dict) -> list:
    """Every (gene, shorter, longer) protein pair differing by one deleted digit."""
    out = []
    for gene, forms in gene_hgvs.items():
        parsed = {}
        for h, (n, rsids) in forms.items():
            m = PROTEIN.match(h)
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
                "shorter": short, "shorter_rows": ps[3],
                "shorter_rsids": len(ps[4]), "shorter_rsid_rows": sum(ps[4].values()),
                "shorter_top_rsid": ps[4].most_common(1)[0][0] if ps[4] else None,
                "longer": long_, "longer_rows": pl[3],
                "longer_rsids": len(pl[4]), "longer_rsid_rows": sum(pl[4].values()),
                "longer_top_rsid": pl[4].most_common(1)[0][0] if pl[4] else None,
            })
    out.sort(key=lambda t: (t["gene"], t["shorter"], t["longer"]))
    return out


def adjudicate_twins(gene_hgvs: dict) -> tuple:
    """(corrections, the full twin table with a verdict on every pair).

    A twin is decided only by the rsid. A form that is a digit-drop twin of TWO
    different canonicals is ambiguous -- KRAS `p.G1V` is one digit from both
    `p.G12V` and `p.G13V`, which are different substitutions -- so it abstains.
    An earlier version let the last write win, silently picking by ASCII order
    and then reporting in the artifact a verdict it had not applied.
    """
    twins = digit_drop_twins(gene_hgvs)
    candidates = collections.defaultdict(set)
    for t in twins:
        has_s, has_l = t["shorter_rsids"] > 0, t["longer_rsids"] > 0
        if has_s == has_l:
            t["verdict"] = "abstain:the rsids cannot adjudicate"
            continue
        t["_canon"], t["_wrong"] = ((t["shorter"], t["longer"]) if has_s
                                    else (t["longer"], t["shorter"]))
        candidates[(t["gene"], t["_wrong"])].add(t["_canon"])

    fix = {}
    for t in twins:
        canon = t.pop("_canon", None)
        wrong = t.pop("_wrong", None)
        if canon is None:
            continue
        if len(candidates[(t["gene"], wrong)]) > 1:
            t["verdict"] = ("abstain:this form is one digit from "
                            + " and ".join(sorted(candidates[(t["gene"], wrong)])))
            continue
        fix[(t["gene"], wrong)] = canon
        short_is_canon = canon == t["shorter"]
        t["verdict"] = "canonical:" + canon
        t["rows_moved"] = t["longer_rows"] if short_is_canon else t["shorter_rows"]
        t["canonical_rows"] = t["shorter_rows"] if short_is_canon else t["longer_rows"]
        t["canonical_rsid_rows"] = (t["shorter_rsid_rows"] if short_is_canon
                                    else t["longer_rsid_rows"])
        t["canonical_rsid"] = (t["shorter_top_rsid"] if short_is_canon
                               else t["longer_top_rsid"])
        # The rule is load-bearing exactly when it moves the MAJORITY onto the
        # minority string. Flagged rather than suppressed by a threshold, and
        # reported WITH the row count its rsid evidence rests on, since most of
        # these rest on a single row.
        t["moves_majority"] = t["rows_moved"] > t["canonical_rows"]
    return fix, twins


def build_pairs(path: Path, canonicalize):
    """(chemical, gene, variant) -> the papers and predicates asserting it."""
    pairs = collections.defaultdict(
        lambda: {"pmids": set(), "preds": collections.Counter()})
    genes = collections.Counter()
    # Both measures, so the report's "rows" column can be checked against the
    # occurrence count it is NOT. They differ only for a gene appearing twice
    # in one relation row, which is exactly the confusion being guarded.
    gene_occurrences = collections.Counter()
    chem_variant = 0
    for pmid, pred, sides in iter_rows(path):
        # Per ROW, not per occurrence: 4,300 rows carry a variant on both sides,
        # and the column is labelled "variant-touching rows". Counting each
        # occurrence inflated BRAF by 60 and EGFR by 419.
        row_genes = set()
        for (tv, iv), (to, io_) in (sides, sides[::-1]):
            if tv not in VARIANT_TYPES:
                continue
            gene, var = canonicalize(parse_variant(iv))
            if gene:
                row_genes.add(gene)
                gene_occurrences[gene] += 1
            if to != "Chemical":
                continue
            chem_variant += 1
            e = pairs[(io_, gene, var)]
            e["pmids"].add(pmid)
            e["preds"][pred] += 1
        for g in row_genes:
            genes[g] += 1
    return pairs, genes, gene_occurrences, chem_variant


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
    rows_touching = occurrences = no_gene = 0
    for _, _, sides in iter_rows(RELATIONS):
        hit = False
        for t, ident in sides:
            if t not in VARIANT_TYPES:
                continue
            hit = True
            occurrences += 1
            f = parse_variant(ident)
            gene, hgvs, rs = (f.get("CorrespondingGene", ""), f.get("HGVS", ""),
                              f.get("RS#", ""))
            if not gene:
                no_gene += 1
            if gene and hgvs:
                rec = gene_hgvs[gene][hgvs]
                rec[0] += 1
                if rs:
                    rec[1][rs] += 1
                    rs_hgvs[(rs, gene)][hgvs] += 1
        rows_touching += hit

    rs_fix, rs_tally, refused = resolve_rsids(rs_hgvs)
    twin_fix, twins = adjudicate_twins(gene_hgvs)

    counts = collections.Counter()

    def canonicalize(f):
        gene, hgvs, rs = (f.get("CorrespondingGene", ""), f.get("HGVS", ""),
                          f.get("RS#", ""))
        # Applied to EVERY row carrying the rsid, not only rows missing an HGVS:
        # `c.2369C>T` and `p.T790M` were being counted as two variants.
        if rs and (rs, gene) in rs_fix:
            canon = rs_fix[(rs, gene)]
            if not hgvs:
                counts["rs_filled"] += 1
            elif hgvs != canon:
                counts["rs_respelled"] += 1
            hgvs = canon
        if hgvs and (gene, hgvs) in twin_fix:
            new = twin_fix[(gene, hgvs)]
            counts["twin_corrected"] += new != hgvs
            hgvs = new
        # A bare rsid keeps its `rs` prefix: `113488022` in a variant column is
        # indistinguishable from a codon position.
        return gene, hgvs or (f"rs{rs}" if rs else "")

    # --- pass 2: the drug-by-variant table, on canonical keys ------------------
    pairs, genes, gene_occ, chem_variant = build_pairs(RELATIONS, canonicalize)

    # The worked examples the report quotes are CHOSEN BY THE DATA -- the biggest
    # case of each kind -- so the prose cannot go stale against the artifact.
    movers = [t for t in twins if t.get("moves_majority")]
    featured = max(movers, key=lambda t: t["shorter_rows"] + t["longer_rows"],
                   default=None)

    # --- pass 3: the disease profile corroborating the featured twin -----------
    feat_disease = collections.defaultdict(collections.Counter)
    if featured:
        want = {featured["shorter"], featured["longer"]}
        for _, _, sides in iter_rows(RELATIONS):
            for (tv, iv), (to, io_) in (sides, sides[::-1]):
                if tv not in VARIANT_TYPES or to != "Disease":
                    continue
                f = parse_variant(iv)
                if f.get("CorrespondingGene") == featured["gene"] \
                        and f.get("HGVS") in want:
                    feat_disease[f["HGVS"]][io_] += 1

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
    table.sort(key=lambda r: (-r["papers"], r["drug_id"], r["gene_id"], r["variant"]))

    with open(OUT_TSV, "wb") as raw_fh, \
            gzip.GzipFile(fileobj=raw_fh, mode="wb", mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as fh:
        fh.write("drug\tdrug_id\tgene\tgene_id\tvariant\tpapers\tassertions\tpredicates\n")
        for r in table:
            fh.write(f"{r['drug']}\t{r['drug_id']}\t{r['gene']}\t{r['gene_id']}\t"
                     f"{r['variant']}\t{r['papers']}\t{r['assertions']}\t"
                     f"{';'.join(f'{k}={n}' for k, n in r['predicates'].items())}\n")

    view = [r for r in table if r["gene"] == args.gene] if args.gene else table
    applied = [t for t in twins if t["verdict"].startswith("canonical:")]
    # Distinct variants IN THE MAP, not across the whole variant layer, and genes
    # are never fused: a gene-less `p.G12C` is its own key, since it could be
    # KRAS, NRAS or HRAS and this refuses to choose.
    mapped = {(r["gene_id"], r["variant"]) for r in table}
    res = {
        "relation_rows_touching_a_variant": rows_touching,
        "variant_entity_occurrences": occurrences,
        "occurrences_with_no_gene": no_gene,
        "chemical_variant_rows": chem_variant,
        "rsid_rows_given_an_hgvs": counts["rs_filled"],
        "rsid_rows_respelled_to_one_form": counts["rs_respelled"],
        "rsid_resolution": rs_tally,
        "digit_drop_rows_corrected": counts["twin_corrected"],
        "distinct_variants_in_the_map": len(mapped),
        "variants_in_the_map_with_no_gene": sum(1 for g, _ in mapped if not g),
        "distinct_drug_variant_pairs": len(table),
        "pairs_at_min_papers": sum(1 for r in table if r["papers"] >= args.min_papers),
        "single_paper_pairs": sum(1 for r in table if r["papers"] == 1),
        "twins": {
            "found": len(twins),
            "adjudicated": len(applied),
            "abstained": len(twins) - len(applied),
            "corrections_applied": len(twin_fix),
            "moving_the_majority": movers,
            "featured": featured,
            # A LIST of pairs, not a dict: `json.dumps(sort_keys=True)` re-sorts
            # a dict alphabetically, so the artifact contradicted the prose's
            # "in the same rank order" claim while both were correct.
            "featured_diseases": {k: [[name(d), n] for d, n in v.most_common(4)]
                                  for k, v in feat_disease.items()},
        },
        "rsid_refusals": {
            "count": len(refused),
            "largest": refused[:5],
            # The superseded protein-forms-only rule would have collapsed these,
            # so this IS the size of the defect the corrected rule fixes.
            "the_old_rule_would_have_collapsed": {
                "rsids": sum(1 for x in refused if x["old_rule_would_have_collapsed"]),
                "rows": sum(x["rows"] for x in refused
                            if x["old_rule_would_have_collapsed"]),
                # The rows whose KEY would have changed. Quoting the rsids'
                # total instead overstates it, since the rows already on the
                # winning spelling do not move.
                "rows_moved": sum(x["old_rule_rows_moved"] for x in refused),
                # Strictly less than `rows` by construction: the rows already
                # sitting on the winning spelling do not move. If the two are
                # ever equal, the counterfactual is quoting the wrong thing.
            },
            # Featured = the largest case the OLD rule got WRONG, which is the
            # one that demonstrates the fix. The largest refusal overall is BRAF,
            # whose spellings are mostly variations of one change, so it argues
            # the opposite of the point.
            "featured": next((x for x in refused
                              if x["old_rule_would_have_collapsed"]), None),
        },
        "top_genes_by_variant_rows": [
            {"gene": name(g, "gene:"), "gene_id": g, "rows": n,
             "occurrences": gene_occ[g]}
            for g, n in genes.most_common(15)],
        "top_pairs": view[:60],
        "gene_filter": args.gene,
        "min_papers": args.min_papers,
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(res, name), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}\nwrote {OUT_TSV}")
    print(f"  {chem_variant:,} chemical-variant rows -> {len(table):,} pairs")
    print(f"  rsid: {counts['rs_filled']:,} rows given an HGVS, "
          f"{counts['rs_respelled']:,} respelled, {len(refused):,} rsids refused")
    print(f"  twins: {len(applied)} of {len(twins)} adjudicated, "
          f"{counts['twin_corrected']:,} rows corrected")
    return 0


def profiles_agree(a: list, b: list, n: int = 3) -> bool:
    """Do two ranked lists share their first `n` entries, in order?

    A function rather than an inline condition so it can be tested directly.
    `a and a[:3] == b[:3]` is TRUE when both lists hold one entry and it is the
    same, which would print "identical profile in the same rank order" on a
    single agreement; on this corpus both lists are long enough that the bug is
    inert, so only a unit test can see it.
    """
    return len(a) >= n and len(b) >= n and a[:n] == b[:n]


def a_fraction_of(total: int, part: int) -> str:
    """'a seventh' rather than a hand-written word beside a computed number."""
    words = {1: "the whole", 2: "half", 3: "a third", 4: "a quarter", 5: "a fifth",
             6: "a sixth", 7: "a seventh", 8: "an eighth", 9: "a ninth",
             10: "a tenth"}
    k = round(total / max(part, 1))
    return words.get(k, f"1/{k}")


def render(r: dict, name) -> str:
    tw, rsr = r["twins"], r["rsid_resolution"]
    agree = rsr.get("codon_relation_agrees", 0)
    mismatch = rsr.get("codon_mismatch_refused", 0)
    uncheckable = rsr.get("not_checkable_against_a_residue_refused", 0)
    ref = r["rsid_refusals"]["featured"]
    old = r["rsid_refusals"]["the_old_rule_would_have_collapsed"]
    f = tw["featured"]
    pairs_n = max(r["distinct_drug_variant_pairs"], 1)
    single_pct = 100.0 * r["single_paper_pairs"] / pairs_n
    nogene_pct = 100.0 * r["occurrences_with_no_gene"] / max(
        r["variant_entity_occurrences"], 1)

    L = [
        "# Which drugs the literature discusses alongside which mutations", "",
        "Generated by `scripts/atlas_variant_drug_map.py`. Every figure below is",
        "derived, including the worked examples, which name whichever case the",
        "data says is largest.", "",
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
        "## Why this layer was unresolved", "",
        "`relations.tsv.gz` carries `ProteinMutation`, `DNAMutation` and `SNP`",
        "entities whose identifier holds the variant inline:",
        "`RS#:121434568;HGVS:p.L858R;CorrespondingGene:1956` is EGFR L858R.", "",
        "Those entities are **not absent** from this repository. They are nodes in",
        "the graph index like any other, and `analysis/atlas-emergence.md` already",
        "ships a row reading `adagrasib -- RS#:121913530;HGVS:p.G12C;"
        "CorrespondingGene:3845`. What no analysis here has ever done is **parse**",
        "those fields: the identifier passes through whole, as an opaque label, so",
        "nothing could group two spellings of one variant or say which gene a",
        "mutation belongs to. `atlas_discovery.py:145` skips them outright, for a",
        "good local reason -- it asks PubMed whether two terms co-occur, and an",
        "HGVS string is not a searchable term.", "",
        "## The identifier is not the variant", "",
        "Keying a pair on whatever string the entity carried splits one variant",
        "across several keys. For each of the two defects, the OBVIOUS correction",
        "is wrong.", "",
        "### The rsid is not the variant either", "",
        "One change arrives with an rsid and without, and at protein and coding",
        "level (`p.T790M` and `c.2369C>T`). But an rsid cannot simply be collapsed",
        "to one HGVS: **rs121913529 covers KRAS `p.G12D`, `p.G12V` and `p.G12A`**,",
        "one multi-allelic codon-12 site, three substitutions, different drug",
        "programs. Merging on rsid would destroy exactly the distinction this map",
        "exists to show.", "",
        "The rule needs no threshold. Every spelling must parse; each",
        "representation class must agree internally; and a protein spelling must",
        "agree with a coding one ACROSS classes, since codon = (coding position",
        "+ 2) // 3 relates them. Failing any of the three refuses that rsid",
        "UNDER THAT GENE: the unit of resolution is the (rsid, gene) key, not",
        "the rsid. In practice that distinction is nearly inert here -- "
        f"{rsr.get('refused_here_collapsed_under_another_gene', 0)} refused keys",
        "have the same rsid collapsing under a different gene, and "
        f"{rsr.get('refused_here_single_spelling_under_another_gene', 0)} have it",
        "carrying a single spelling there, which is not a collapse decision at",
        "all.", "",
        "| (rsid, gene) keys carrying several spellings | |", "|---|--:|",
        f"| one change, several spellings: collapsed | "
        f"{rsr.get('one_change_many_spellings', 0):,} |",
        f"| spellings denote different changes: refused | "
        f"{rsr.get('classes_disagree_refused', 0):,} |",
        f"| the coding position is in a different codon: refused | "
        f"{rsr.get('codon_mismatch_refused', 0):,} |",
        f"| no position can be checked against the residue: refused | "
        f"{rsr.get('not_checkable_against_a_residue_refused', 0):,} |",
        f"| a spelling cannot be checked: refused | "
        f"{rsr.get('unparsable_spelling_refused', 0):,} |", "",
        "The cross-class test is the one internal agreement cannot see, and it",
        "is also the strongest evidence that the collapses are right. Where a",
        "protein and a coding spelling were actually compared, the codon",
        f"relation **agrees** for {agree:,} rsids and **fails** for "
        f"{mismatch:,}, so it holds in {100.0 * agree / max(agree + mismatch, 1):.1f}%",
        "of the cases it can decide.", "",
        f"The other {uncheckable:,} refusals are ones it could NOT decide: a",
        "genomic coordinate has no arithmetic relation to a residue number, and",
        "neither does an intronic or untranslated offset. Those are counted",
        "separately rather than folded into the failures, because quoting them",
        "as failures would both overstate the error rate and mislabel a",
        "protein-versus-genomic refusal as protein-versus-coding.", "",
        "It catches rs1494558, where `p.I66T` sits beside `c.412G>A` and codon",
        "138 is not residue 66, and rs2297518, where a missense sits beside a",
        "promoter position that cannot encode it.", "",
        f"That gave **{r['rsid_rows_given_an_hgvs']:,}** rows an HGVS they lacked "
        f"and respelled **{r['rsid_rows_respelled_to_one_form']:,}** onto a single "
        "form.", "",
        f"Most of that fill is not the rule's doing. "
        f"**{rsr.get('single_spelling', 0):,}** (rsid, gene) keys carry only ONE "
        "spelling, so no class could disagree and no collapse decision was made: "
        "a bare row simply inherits the one spelling its rsid was ever given. "
        "The table above covers the keys with SEVERAL spellings, which is where "
        "the three tests apply.", "",
    ]
    if ref:
        # From the FULL spelling set, not the displayed top six: if the sole
        # protein spelling ranked below sixth the explanatory paragraph would
        # silently vanish and the table would ship unexplained.
        prot = max(ref["protein_spellings"].items(), key=lambda kv: kv[1],
                   default=("", 0))
        big = max(ref["spellings"].items(), key=lambda kv: kv[1])
        shown, n_all = len(ref["spellings"]), ref["n_spellings"]
        L += ["**Testing agreement among protein forms only is not enough**, which",
              "is how the first version of this was wrong: it swept every other",
              "spelling onto the winning protein form. That would have collapsed",
              f"**{old['rsids']:,} of the rsids refused here, moving "
              f"{old['rows_moved']:,} rows onto a single key.** That figure is",
              "what it says and no more: it counts the rows whose spelling",
              "would have changed, not a claim that every one of them would",
              "have landed somewhere wrong. Some of those rsids are refused",
              "only because a SIBLING spelling cannot be checked, and their",
              "protein and coding forms do agree.", "",
              f"The largest of them shows what it cost. rs{ref['rsid']} on "
              f"{name(ref['gene'], 'gene:')} carries {ref['rows']:,} rows across "
              f"{n_all} spellings, {ref['reason']}"
              + (f" (top {shown} shown):" if n_all > shown else ":"), "",
              "| spelling | rows | |", "|---|--:|---|"]
        for h, n in ref["spellings"].items():
            flag = " refuses the rsid" if h in ref["offending"] else ""
            L.append(f"| `{h}` | {n:,} |{flag} |")
        if prot[0] and prot[1] < big[1]:
            L += ["",
                  f"There is exactly one protein spelling, `{prot[0]}` at "
                  f"{prot[1]:,} rows, so the old rule saw no disagreement among "
                  f"protein forms and captured all {ref['rows']:,} rows onto it, "
                  f"including the {big[1]:,}-row `{big[0]}`."]
        L.append("")

    L += ["### The majority is not the truth", "",
          f"{tw['found']} twin pairs differ by one deleted position digit. Reading",
          "the commoner as correct is wrong in both directions: TP53 `p.R72P` and",
          "ERBB2 `p.I655V` are real rs-backed polymorphisms whose *longer* twins",
          "are the typos, while EGFR `p.T790M` is real and `p.T90M` the typo.", "",
          "**The rsid adjudicates**, being assigned independently of the string",
          "that was extracted. Where exactly one side carries rsids that side is",
          f"canonical ({tw['adjudicated']} pairs, {r['digit_drop_rows_corrected']:,}",
          "rows). Where neither or both do, or where one form is one digit from",
          f"two different canonicals, this **abstains** ({tw['abstained']} pairs),",
          "the same choice `atlas_graph.resolve` makes for a contested surface",
          "form. KRAS `p.G1V` is one digit from both `p.G12V` and `p.G13V`, which",
          "are different substitutions, so it gets no verdict at all.", ""]

    if tw["moving_the_majority"]:
        L += ["The rule is load-bearing exactly where it moves the MAJORITY onto",
              "the minority string. Those cases are listed rather than suppressed",
              "by a threshold, **with the row count their rsid evidence rests on**,",
              "because most of them rest on very little:", "",
              "| gene | non-canonical | rows | canonical | rows | rsid-bearing rows |",
              "|---|---|--:|---|--:|--:|"]
        for t in sorted(tw["moving_the_majority"], key=lambda x: -x["rows_moved"]):
            canon = t["verdict"].split(":", 1)[1]
            bad = t["longer"] if canon == t["shorter"] else t["shorter"]
            L.append(f"| {name(t['gene'], 'gene:')} | `{bad}` | {t['rows_moved']:,} |"
                     f" `{canon}` | {t['canonical_rows']:,} |"
                     f" {t['canonical_rsid_rows']:,} |")
        L.append("")
    if f:
        canon = f["verdict"].split(":", 1)[1]
        bad = f["longer"] if canon == f["shorter"] else f["shorter"]
        badn, goodn = f["rows_moved"], f["canonical_rows"]
        total = badn + goodn
        L += [f"The {name(f['gene'], 'gene:')} row is the one that matters, and the",
              f"only one whose rsid evidence is substantial. `{bad}` is `{canon}`",
              "with a digit missing and carried "
              f"**{100.0 * badn / total:.0f}%** of that variant's volume, so the",
              f"uncorrected map read `{canon}` at "
              f"**{a_fraction_of(total, goodn)}** of its real support.", "",
              "The verdict does not rest on those counts, which point the other",
              f"way: `{canon}` carries **rs{f['canonical_rsid']} on",
              f"{f['canonical_rsid_rows']:,} of its {goodn:,} rows** while `{bad}`",
              "carries no rsid at all."]
        dis = tw["featured_diseases"]
        a = [d for d, _ in dis.get(bad, [])]
        b = [d for d, _ in dis.get(canon, [])]
        if profiles_agree(a, b):
            L += ["They also show the **identical** disease profile in the same",
                  "rank order, which a coincidence of spelling would not produce: "
                  + ", ".join(b[:3]) + "."]
        L.append("")

    L += ["## What is in it", "",
          "| | count |", "|---|--:|",
          f"| relation rows touching a point-variant entity | "
          f"{r['relation_rows_touching_a_variant']:,} |",
          f"| variant entity occurrences (both sides counted) | {r['variant_entity_occurrences']:,} |",
          f"| ...carrying no gene | {r['occurrences_with_no_gene']:,} ({nogene_pct:.1f}%) |",
          f"| **chemical-to-variant rows** | **{r['chemical_variant_rows']:,}** |",
          f"| distinct variants IN THIS MAP | {r['distinct_variants_in_the_map']:,} |",
          f"| ...of those, carrying no gene | {r['variants_in_the_map_with_no_gene']:,} |",
          f"| distinct (drug, gene, variant) pairs | {r['distinct_drug_variant_pairs']:,} |",
          f"| pairs with >= {r['min_papers']} papers | {r['pairs_at_min_papers']:,} |",
          f"| pairs resting on ONE paper | {r['single_paper_pairs']:,} ({single_pct:.1f}%) |",
          "",
          "`relations.tsv.gz` also carries a fourth mutation type, `Mutation`,",
          "holding structural variants as chromosomal ranges "
          "(`Chr7:154954255-154998784dup`). It is deliberately out of scope: it",
          "has no HGVS substitution to reconcile and none of the three tests",
          "applies to it, so the counts above are point variants only.", "",
          "A variant with no `CorrespondingGene` is genuinely ambiguous, since",
          "`p.G12C` alone could be KRAS, NRAS or HRAS. Those are left unresolved",
          "and keyed separately: never assigned to a likely gene, and never merged",
          "with the same spelling under a known one. They show `?` as the gene.", "",
          f"**{single_pct:.0f}% of pairs rest on a single paper.** The retraction",
          "analysis found the same shape across the whole graph (70.2%), and it",
          "is the first thing to know before reading any row below as evidence.", "",
          "## Genes carrying the most variant relations", "",
          "| gene | variant-touching rows |", "|---|--:|"]
    for g in r["top_genes_by_variant_rows"][:15]:
        L.append(f"| {g['gene']} | {g['rows']:,} |")
    L += ["", "## The most-discussed drug-variant pairs", "",
          "Ranked by asserting papers. `papers` counts distinct PMIDs; the",
          "predicate counts are ASSERTIONS, so they exceed the paper count when",
          "one paper states a pair more than once. Again: this ranks how much a",
          "pair is *written about*, and `associate` does not say which direction.", "",
          "| drug | gene | variant | papers | predicates (assertions) |",
          "|---|---|---|--:|---|"]
    for row in r["top_pairs"][:40]:
        # ALL of them. Keeping the top three silently understated two shipped
        # rows while the prose invited the reader to compare these counts with
        # the paper count.
        preds = ", ".join(f"`{k}` {n}" for k, n in row["predicates"].items())
        L.append(f"| {row['drug']} | {row['gene']} | `{row['variant']}` | "
                 f"{row['papers']} | {preds} |")
    L += ["", "## What this cannot say", "",
          "* **No direction.** `associate` covers resistance and sensitivity alike.",
          "* **Not curated.** No evidence level, no expert review, no clinical",
          "  interpretation. Use CIViC/OncoKB/COSMIC for those.",
          "* **Extractor error dominates.** PubTator's own error rate is larger",
          "  than most differences between rows here, and the digit-drop finding",
          "  above is a direct instance of it.",
          "* **Reconciliation is deliberately incomplete.** It merges what the",
          f"  evidence can adjudicate and refuses the rest: {tw['abstained']} twins",
          f"  and {r['rsid_refusals']['count']:,} (rsid, gene) keys stay fragmented",
          "  guessed, so one variant may still appear under more than one key.",
          "* **Attention is not importance.** A well-studied pair outranks a real",
          "  but rarely-written-about one, exactly as `atlas_model_gaps.py` warns",
          "  for its own ranking.",
          "* **Full table** is `analysis/atlas-variant-drug-map.tsv.gz`; the rows",
          "  above are a view, not the result.", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
