#!/usr/bin/env python3
"""What the census's co-treatment layer can and cannot support.

WHAT THIS READS
---------------
`relations.tsv.gz` carries 224,146 `cotreat` rows over 40,878 drug pairs. Those
counts are NOT new here -- `atlas-retraction-exposure.md` already ships the
224,146 and `atlas-emergence.md` already ships per-pair cotreat counts. What
nothing had done is JOIN them to the gene and variant layers to ask what each
combination is studied against.

The question that motivated it was the therapeutically interesting one: which
resistance mutations have a combination answer and which do not. Answering it
required establishing two things first, and both turned out to be findings.

1. `cotreat` IS NOT CO-ADMINISTRATION
-------------------------------------
The predicate name promises that two drugs were given together. It does not
deliver that: alectinib and crizotinib are sequential ALK inhibitors, never
co-administered, and they carry 86 `cotreat` rows. Osimertinib and gefitinib,
likewise mutually exclusive, carry 28.

The corpus supplies its own discriminator, because it has a separate `compare`
predicate over the same pairs, and the RATIO separates cleanly on a panel where
the truth is known independently: every genuinely co-administered pair scores
above 3.7 and every sequential-or-compared pair below 1.3, with no overlap.
Every row below therefore carries its cotreat:compare ratio, and a row where
`compare` dominates is flagged rather than silently presented as a combination.

The ratio is a HEURISTIC on a nine-pair panel, not a validated classifier. It
is reported so a reader can discount rows, not so anything can be filtered out
automatically.

2. WHAT THE VARIANT JOIN RECOVERS IS SET BY PARTNER-AGENT CLASS
---------------------------------------------------------------
Two earlier explanations were measured and refuted. The first said a
variant-level absence carries NO information, generalising from one failed
control, and the shipped table's own top row refutes it. The second said the
recoveries split by whether the biomarker is a SUBSTITUTION -- also wrong,
because two panel labels were the author's to assign and both were assigned in
the direction that made the split clean.

Corrected and extended to 14 regimens: biomarker class does not separate
(7/12 against 1/2) and partner-agent class does (7/7 against 1/7). FLAURA2 is
labelled with the substitution it actually enrols, the same one as MARIPOSA,
and still misses -- because pemetrexed is a chemotherapy partner.

WHAT IS MEASURED AND WHAT IS ONLY OBSERVED. The join needs BOTH drugs tied to
the same variant IN ONE PAPER. Whether the extractor CAN tie a chemotherapy or
endocrine agent to a substitution at all is a separate question, and it is
measured here per missed regimen rather than asserted: for each drug, whether
it is tied to the target variant anywhere in the corpus, beside whether the
pair is tied to it in a single paper.

WHAT NEITHER LAYER CAN SAY
--------------------------
* `cotreat` carries no outcome. Nothing here says a combination worked.
* A GENE IS NOT A TARGET. Docetaxel + estramustine ties to KLK3, which is PSA,
  a response biomarker.
* This reads PubTator's OWN gene assignment, so the `atlas_ambiguity` blocklist
  -- which guards the alias-resolution path in `atlas_graph.resolve` -- does not
  apply. The collision examples in the report are DERIVED from the authority
  table at render time, because a hand-written one was wrong: an earlier version
  claimed `NP` reaches Neptunium "verified in the authority table itself", and
  Neptunium's alias list contains no `NP` at all.
* THE GENE LAYER'S PRECISION IS NOT MEASURED, AT ALL. No sample has been judged
  and no precision figure is claimed. An earlier version of this section
  reported a 40-row sample that this repository never ran.

Usage:
    python scripts/atlas_combination_gaps.py
"""

import collections
import gzip
import io
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RELATIONS = PROJECT_ROOT / "corpus" / "atlas" / "relations" / "relations.tsv.gz"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-combination-gaps.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-combination-gaps.json"
OUT_TSV = PROJECT_ROOT / "analysis" / "atlas-combination-gaps.tsv.gz"

VARIANT_TYPES = ("ProteinMutation", "DNAMutation", "SNP")

# 37 pairs whose real-world relationship is known INDEPENDENTLY of this corpus,
# used to locate the boundary the co-administration flag sits on. An earlier
# version cited "an independent 36-pair sweep" that this repository had never
# run, in three places, to justify a threshold that flags every shipped row.
# The sweep is now this panel and its boundary is derived below.
#
# NAMES, resolved through the authority table at runtime and asserted to
# resolve. An earlier version hardcoded MeSH identifiers from memory and got
# several wrong -- MESH:C571179 is eravacycline, not vemurafenib -- which
# produced a row reading MISSED. A false negative that looks like a finding is
# the defect class this whole document is about, so an unresolvable name is a
# hard failure rather than a miss.
SEMANTIC_PANEL = [
    # Co-administered: every one a standard regimen given together.
    ("trastuzumab", "pertuzumab", True), ("dabrafenib", "trametinib", True),
    ("encorafenib", "binimetinib", True), ("ipilimumab", "nivolumab", True),
    ("cisplatin", "etoposide", True), ("carboplatin", "paclitaxel", True),
    ("fluorouracil", "leucovorin", True), ("oxaliplatin", "fluorouracil", True),
    ("doxorubicin", "cyclophosphamide", True),
    ("rituximab", "cyclophosphamide", True), ("gemcitabine", "cisplatin", True),
    ("lenalidomide", "dexamethasone", True),
    ("bortezomib", "dexamethasone", True), ("venetoclax", "azacitidine", True),
    ("cetuximab", "irinotecan", True), ("pembrolizumab", "pemetrexed", True),
    ("olaparib", "bevacizumab", True), ("tucatinib", "trastuzumab", True),
    ("abiraterone", "prednisone", True), ("cytarabine", "daunorubicin", True),
    # Never co-administered: sequential lines, or two agents of one class.
    ("alectinib", "crizotinib", False), ("ceritinib", "crizotinib", False),
    ("osimertinib", "gefitinib", False), ("lorlatinib", "crizotinib", False),
    ("nivolumab", "pembrolizumab", False), ("cetuximab", "panitumumab", False),
    ("paclitaxel", "docetaxel", False), ("cisplatin", "carboplatin", False),
    ("anastrozole", "letrozole", False), ("azacitidine", "decitabine", False),
    ("everolimus", "temsirolimus", False), ("palbociclib", "ribociclib", False),
    ("sorafenib", "regorafenib", False), ("dasatinib", "nilotinib", False),
    ("sunitinib", "pazopanib", False), ("afatinib", "osimertinib", False),
    ("vemurafenib", "dabrafenib", False),
]

# Regimens with a molecular biomarker. ELEVEN are approved combinations; three
# (NEJ009, AG221-005, KRYSTAL-7) are major published trials without a combination
# approval, and all three sit in the not-both-targeted arm and all three MISS --
# so each widens the partner separation and narrows the biomarker one. Named
# rather than quietly dropped, and the approved-only figures are reported too.
# Used to measure what the
# variant-level join RECOVERS rather than asserting that it recovers nothing.
# An empty `variant` means the regimen's biomarker is stated at gene level,
# which is the failure mode the misses share.
# (drug a, drug b, gene, variant, trial, both drugs are targeted agents)
#
# `variant` empty means the regimen's biomarker is stated at gene level.
# TWO LABELS I ORIGINALLY GOT WRONG, both in the direction that flattered my
# first explanation: FLAURA2 enrols "exon 19 deletion OR L858R", the SAME
# biomarker as MARIPOSA, so it cannot be gene-keyed while MARIPOSA is
# substitution-keyed; and CAPItello-291's indication is PIK3CA/AKT1/PTEN-
# altered, a three-gene alteration list, so it is gene-keyed rather than
# E17K-keyed. Corrected, the biomarker split stops separating at all -- which
# is what led to the partner-agent column.
RECOVERY_PANEL = [
    ("dabrafenib", "trametinib", "673", "p.V600E", "COMBI-d", True),
    ("encorafenib", "cetuximab", "673", "p.V600E", "BEACON", True),
    ("encorafenib", "binimetinib", "673", "p.V600E", "COLUMBUS", True),
    ("vemurafenib", "cobimetinib", "673", "p.V600E", "coBRIM", True),
    ("amivantamab", "lazertinib", "1956", "p.L858R", "MARIPOSA", True),
    ("sotorasib", "panitumumab", "3845", "p.G12C", "CodeBreaK 300", True),
    ("adagrasib", "cetuximab", "3845", "p.G12C", "KRYSTAL-1", True),
    ("osimertinib", "pemetrexed", "1956", "p.L858R", "FLAURA2", False),
    ("gefitinib", "carboplatin", "1956", "p.L858R", "NEJ009", False),
    ("ivosidenib", "azacitidine", "3417", "p.R132H", "AGILE", False),
    ("enasidenib", "azacitidine", "3418", "p.R140Q", "AG221-005", False),
    ("adagrasib", "pembrolizumab", "3845", "p.G12C", "KRYSTAL-7", False),
    ("alpelisib", "fulvestrant", "5290", "", "SOLAR-1", False),
    ("capivasertib", "fulvestrant", "207", "", "CAPItello-291", False),
]

# The ratio above which a pair reads as co-administered. 3.0 sits inside the
# gap the committed 37-pair panel MEASURES, reported in `co_admin_boundary`.
# Reported, never used to filter: no panel this size can support a classifier.
CO_ADMIN_RATIO = 3.0


def resolve_drug(name: str, lab: dict) -> str:
    """Authority name -> identifier, or raise.

    Raising is the point. A panel entry that silently fails to resolve emits a
    row reading MISSED, which is indistinguishable from a real negative and is
    exactly the failure this document exists to warn about.
    """
    hits = [k for k, v in lab.items() if v.lower() == name.lower()]
    if len(hits) != 1:
        raise SystemExit(
            f"panel drug {name!r} resolves to {len(hits)} identifiers "
            f"({hits[:4]}); a panel entry must resolve to exactly one, because "
            "an unresolved one would ship as a false MISSED")
    return hits[0]


def derive_collisions(paths, wanted) -> list:
    """Confirm each claimed gene-symbol collision against the authority table.

    Hand-writing these is how `NP -> Neptunium` shipped under the words
    "verified in the authority table itself" when Neptunium's alias list holds
    no `NP` at all. Each candidate is now checked against the FULL alias list
    and one that does not hold is dropped rather than printed.
    """
    full = {}
    with gzip.open(paths, "rt") as fh:
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) >= 2:
                full[q[0]] = [x.strip() for x in q[1].split("|")]
    out = []
    for alias, ident, note in wanted:
        names = full.get(ident)
        if names and alias in names:
            out.append({"alias": alias, "resolves_to": names[0], "id": ident,
                        "note": note})
    return out


def load_labels() -> dict:
    path = PROJECT_ROOT / "analysis" / "comention" / "authority-labels.tsv.gz"
    lab = {}
    if not path.exists():
        return lab
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1]:
                lab[p[0]] = p[1].split("|")[0]
    return lab


def scan(path: Path) -> dict:
    """One pass over the relation file."""
    cotreat = collections.defaultdict(set)
    pair_pred = collections.defaultdict(collections.Counter)
    drug_gene = collections.defaultdict(lambda: collections.defaultdict(set))
    drug_var = collections.defaultdict(lambda: collections.defaultdict(set))
    has_var, genes_in = set(), collections.defaultdict(set)
    cotreat_rows = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, pred = parts[0], parts[1]
            ta, _, ia = parts[2].partition("|")
            tb, _, ib = parts[3].partition("|")
            if ta == "Chemical" and tb == "Chemical":
                pair_pred[frozenset((ia, ib))][pred] += 1
                if pred == "cotreat":
                    cotreat_rows += 1
                    cotreat[pmid].add(frozenset((ia, ib)))
            for (tx, ix), (ty, iy) in (((ta, ia), (tb, ib)), ((tb, ib), (ta, ia))):
                if tx in VARIANT_TYPES:
                    has_var.add(pmid)
                    if ty == "Chemical":
                        f = dict(k.split(":", 1) for k in ix.split(";") if ":" in k)
                        drug_var[pmid][iy].add(
                            (f.get("CorrespondingGene", ""), f.get("HGVS", "")))
                elif tx == "Gene":
                    genes_in[pmid].add(ix)
                    if ty == "Chemical":
                        drug_gene[pmid][iy].add(ix)
    return {"cotreat": cotreat, "cotreat_rows": cotreat_rows,
            "pair_pred": pair_pred, "drug_gene": drug_gene, "drug_var": drug_var,
            "has_var": has_var, "genes_in": genes_in}


def join(cotreat: dict, drug_target: dict) -> dict:
    """(drug pair, target) -> the papers tying BOTH drugs to that target."""
    out = collections.defaultdict(set)
    for pmid, pairs in cotreat.items():
        dt = drug_target.get(pmid)
        if not dt:
            continue
        for pair in pairs:
            a, b = tuple(pair)
            for target in dt.get(a, set()) & dt.get(b, set()):
                out[(pair, target)].add(pmid)
    return out


def ratio(pair_pred: dict, pair):
    """(cotreat, compare, ratio-or-None). None when there is nothing to divide by.

    `cotreat / max(compare, 1)` returns the raw cotreat COUNT when a pair has no
    `compare` rows, and the flag built on it then labelled pairs with no
    comparison evidence at all "compare-dominated". Most shipped rows are in
    that state, so the label misdescribed the majority of the table.
    """
    c = pair_pred.get(pair, {})
    ct, cm = c.get("cotreat", 0), c.get("compare", 0)
    return ct, cm, (round(ct / cm, 2) if cm else None)


def main() -> int:
    lab = load_labels()

    def nm(i, prefix=""):
        n = lab.get(i)
        if n is None:
            return f"{prefix}{i}" if i else "?"
        return n[len("[OBSOLETE] "):] if n.startswith("[OBSOLETE] ") else n

    s = scan(RELATIONS)
    pp = s["pair_pred"]
    by_gene = join(s["cotreat"], s["drug_gene"])
    by_variant = join(s["cotreat"], s["drug_var"])

    # --- 1. does `cotreat` mean co-administration? ---------------------------
    semantic = []
    for an, bn, truth in SEMANTIC_PANEL:
        a, b = resolve_drug(an, lab), resolve_drug(bn, lab)
        c = pp.get(frozenset((a, b)), {})
        ct, cm, rr = ratio(pp, frozenset((a, b)))
        if not ct:
            continue          # no cotreat rows: the pair says nothing here
        semantic.append({"pair": f"{an} + {bn}", "co_administered": truth,
                         "cotreat": ct, "compare": cm, "ratio": rr})
    # Only pairs with BOTH counts can locate a boundary; a null ratio has no
    # position on the axis.
    usable = [x for x in semantic if x["ratio"] is not None]
    lo_true = min((x["ratio"] for x in usable if x["co_administered"]), default=0)
    hi_false = max((x["ratio"] for x in usable if not x["co_administered"]),
                   default=0)

    # --- 2. what does the variant join recover? ------------------------------
    recovery = []
    for an, bn, gene, variant, trial, both_targeted in RECOVERY_PANEL:
        a, b = resolve_drug(an, lab), resolve_drug(bn, lab)
        pair = frozenset((a, b))
        n = sum(len(pm) for (pr, t), pm in by_variant.items()
                if pr == pair and t[0] == gene and (not variant or t[1] == variant))
        recovery.append({
            "trial": trial, "regimen": f"{nm(a)} + {nm(b)}",
            "drug_ids": sorted((a, b)),
            "gene_id": gene, "variant_key": variant,
            "gene": nm(gene, "gene:"), "variant": variant or "(gene-level)",
            "biomarker_is_a_substitution": bool(variant),
            "both_drugs_are_targeted": both_targeted,
            "papers": n, "recovered": n > 0})
    def split(key):
        yes = [r for r in recovery if r[key]]
        no = [r for r in recovery if not r[key]]
        return {"yes_recovered": sum(1 for r in yes if r["recovered"]),
                "yes_total": len(yes),
                "no_recovered": sum(1 for r in no if r["recovered"]),
                "no_total": len(no)}
    by_biomarker = split("biomarker_is_a_substitution")
    by_partner = split("both_drugs_are_targeted")
    # Which explanation fits? Separation = the gap between the two rates.
    def gap(x):
        return (x["yes_recovered"] / max(x["yes_total"], 1)
                - x["no_recovered"] / max(x["no_total"], 1))

    # WHY each miss misses, MEASURED. The document previously asserted that the
    # extractor cannot tie a chemotherapy or endocrine agent to a substitution.
    # Nothing tested that. What the join actually requires is both drugs tied to
    # one variant IN ONE PAPER, so the two are separated here: is each drug ever
    # tied to the target variant anywhere, and are both tied to it together?
    anywhere = collections.defaultdict(set)
    for pmid, dv in s["drug_var"].items():
        for drug, vs in dv.items():
            for v in vs:
                anywhere[(drug, v)].add(pmid)
    for r_ in recovery:
        if r_["recovered"]:
            continue
        a2, b2 = r_["drug_ids"]
        key = (r_["gene_id"], r_["variant_key"])
        r_["drug_a_tied_anywhere"] = len(anywhere.get((a2, key), ()))
        r_["drug_b_tied_anywhere"] = len(anywhere.get((b2, key), ()))
        r_["both_tied_but_never_together"] = bool(
            r_["drug_a_tied_anywhere"] and r_["drug_b_tied_anywhere"])

    # The control, kept because it is what made the failure MODE visible.
    ctrl_pair = frozenset((resolve_drug("alpelisib", lab),
                           resolve_drug("fulvestrant", lab)))
    ctrl_papers = {p for p, prs in s["cotreat"].items() if ctrl_pair in prs}
    control = {
        "regimen": "alpelisib + fulvestrant (SOLAR-1)",
        "papers_asserting_the_combination": len(ctrl_papers),
        "of_those_annotating_the_gene": sum(
            1 for p in ctrl_papers if "5290" in s["genes_in"][p]),
        "of_those_annotating_any_variant": sum(
            1 for p in ctrl_papers if p in s["has_var"]),
        "found_by_the_gene_join": len(by_gene.get((ctrl_pair, "5290"), ())),
    }

    def table(joined, is_variant):
        rows = []
        for (pr, t), pm in joined.items():
            a, b = sorted(pr)
            ct, cm, rr = ratio(pp, pr)
            rows.append({
                "drug_a": nm(a), "drug_b": nm(b), "drug_a_id": a, "drug_b_id": b,
                "target": nm(t[0] if is_variant else t, "gene:"),
                "target_id": t[0] if is_variant else t,
                "variant": t[1] if is_variant else "",
                "papers": len(pm), "cotreat": ct, "compare": cm,
                "cotreat_compare_ratio": rr,
                # THREE states, not two. A pair with no `compare` rows has no
                # comparison evidence either way and must not be called
                # compare-dominated.
                "verdict": ("no comparison evidence" if rr is None
                            else "reads as co-administration" if rr >= CO_ADMIN_RATIO
                            else "compare-dominated")})
        rows.sort(key=lambda x: (-x["papers"], x["drug_a_id"], x["drug_b_id"],
                                 x["target_id"], x["variant"]))
        return rows

    gene_rows, var_rows = table(by_gene, False), table(by_variant, True)

    with open(OUT_TSV, "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as fh:
        fh.write("layer\tdrug_a\tdrug_b\ttarget\ttarget_id\tvariant\tpapers"
                 "\tcotreat\tcompare\tcotreat_compare_ratio\tverdict\n")
        for layer, rows in (("gene", gene_rows), ("variant", var_rows)):
            for r in rows:
                rr = r["cotreat_compare_ratio"]
                fh.write(f"{layer}\t{r['drug_a']}\t{r['drug_b']}\t{r['target']}\t"
                         f"{r['target_id']}\t{r['variant']}\t{r['papers']}\t"
                         f"{r['cotreat']}\t{r['compare']}\t"
                         f"{'' if rr is None else rr}\t{r['verdict']}\n")

    cot_papers = set(s["cotreat"])
    overlap = len(cot_papers & s["has_var"])
    res = {
        "cotreat_rows": s["cotreat_rows"],
        "paper_pair_combinations": sum(len(v) for v in s["cotreat"].values()),
        "distinct_drug_pairs": len({p for v in s["cotreat"].values() for p in v}),
        "papers_asserting_a_cotreatment": len(cot_papers),
        "papers_asserting_a_variant_relation": len(s["has_var"]),
        "papers_asserting_both": overlap,
        "overlap_share_of_cotreatment_papers": round(
            100.0 * overlap / max(len(cot_papers), 1), 2),
        "cotreat_semantics": {
            "panel": semantic,
            "ratio_separates_the_panel": lo_true > hi_false,
            "lowest_co_administered_ratio": lo_true,
            "highest_sequential_ratio": hi_false,
            "pairs_where_compare_outnumbers_cotreat": sum(
                1 for c in pp.values()
                if c.get("cotreat", 0) and c.get("compare", 0) > c["cotreat"]),
        },
        "recovery": {
            "panel": recovery,
            "recovered": sum(1 for r in recovery if r["recovered"]),
            "of": len(recovery),
            "by_biomarker_class": {**by_biomarker,
                                   "separation": round(gap(by_biomarker), 3)},
            "by_partner_class": {**by_partner,
                                 "separation": round(gap(by_partner), 3)},
            "partner_class_explains_better":
                gap(by_partner) > gap(by_biomarker),
        },
        "collisions": derive_collisions(
            PROJECT_ROOT / "analysis" / "comention" / "authority-labels.tsv.gz",
            [("MEK", "5609", "MAP2K1 carries `MEK1`, so `MEK` alone lands on MKK7"),
             ("Met", "79811", "an alias of SLTM, not the MET receptor"),
             ("PGP", "283871", "phosphoglycolate phosphatase, not P-glycoprotein"),
             ("NO", "MESH:D009614", "nobelium, not nitric oxide"),
             ("PSA", "354", "KLK3, a response biomarker rather than a target"),
             ("NP", "MESH:D009405", "neptunium (claimed once and NOT true)")]),
        "co_admin_boundary": {
            "lowest_co_administered": lo_true,
            "highest_sequential": hi_false,
            "threshold": CO_ADMIN_RATIO,
            "threshold_inside_the_gap": hi_false < CO_ADMIN_RATIO <= lo_true,
            "panel_pairs": len(semantic),
        },
        "control": control,
        "gene_layer": {
            "triples": len(gene_rows),
            "single_paper": sum(1 for r in gene_rows if r["papers"] == 1),
            "at_least_three_papers": sum(1 for r in gene_rows if r["papers"] >= 3),
            "by_verdict": dict(collections.Counter(
                r["verdict"] for r in gene_rows)),
            "top": gene_rows[:40],
        },
        "variant_layer": {
            "triples": len(var_rows),
            "single_paper": sum(1 for r in var_rows if r["papers"] == 1),
            "top": var_rows[:25],
        },
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(res), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}\nwrote {OUT_TSV}")
    print(f"  gene layer {len(gene_rows):,}, variant layer {len(var_rows):,}")
    print(f"  cotreat:compare separates the panel: {lo_true > hi_false} "
          f"({lo_true} vs {hi_false})")
    print(f"  variant join recovers {res['recovery']['recovered']}/"
          f"{res['recovery']['of']} known regimens")
    return 0


def render(r: dict) -> str:
    sem, rec, c = r["cotreat_semantics"], r["recovery"], r["control"]
    g, v = r["gene_layer"], r["variant_layer"]
    alectinib = next(x["cotreat"] for x in sem["panel"] if "alectinib" in x["pair"])
    L = [
        "# What the co-treatment layer can and cannot support", "",
        "Generated by `scripts/atlas_combination_gaps.py`.", "",
        f"The census carries {r['cotreat_rows']:,} `cotreat` rows over "
        f"{r['distinct_drug_pairs']:,} drug pairs. Those counts are **not new",
        f"here**: `atlas-retraction-exposure.md` already ships the "
        f"{r['cotreat_rows']:,} and `atlas-emergence.md` already ships per-pair",
        "cotreat counts. What nothing had done is **join** them to the gene and",
        "variant layers to ask what each combination is studied against.", "",
        "The question that motivated it was the therapeutically interesting one:",
        "which resistance mutations have a combination answer and which do not.",
        "Answering it required establishing two things first, and both turned",
        "out to be findings.", "",
        "## 1. `cotreat` is not co-administration", "",
        "The predicate name promises that two drugs were given together. It does",
        "not deliver that. Alectinib and crizotinib are sequential ALK",
        f"inhibitors, never co-administered, and they carry {alectinib} `cotreat`",
        "rows.", "",
        "The corpus supplies its own discriminator, because it carries a separate",
        "`compare` predicate over the same pairs. On a panel where the truth is",
        "known independently of the corpus, the ratio separates:", "",
        "| pair | cotreat | compare | ratio | truth |", "|---|--:|--:|--:|---|"]
    for x in sem["panel"]:
        L.append(f"| {x['pair']} | {x['cotreat']:,} | {x['compare']:,} | "
                 f"{x['ratio']} | "
                 f"{'co-administered' if x['co_administered'] else 'sequential or compared'} |")
    L += ["",
          (f"Every co-administered pair scores at or above "
           f"**{sem['lowest_co_administered_ratio']}** and every sequential one "
           f"at or below **{sem['highest_sequential_ratio']}**, with no overlap."
           if sem["ratio_separates_the_panel"] else
           "**The ratio does NOT separate this panel**, so it cannot be used to "
           "discount rows and the tables below carry it for information only."),
          "",
          f"Across the whole corpus, "
          f"**{sem['pairs_where_compare_outnumbers_cotreat']:,}** pairs carrying",
          "a `cotreat` row have MORE `compare` rows than `cotreat` ones. Every",
          "row below therefore carries its ratio, and one where `compare`",
          "dominates is flagged rather than silently presented as a combination.", "",
          "**This is a heuristic on a nine-pair panel, not a validated",
          "classifier.** It is reported so a reader can discount rows. Nothing",
          "is filtered out automatically.", "",
          "## 2. What the variant join recovers, and why", "",
          "An earlier version of this analysis claimed a variant-level absence",
          "carries NO information, generalising from one control that failed.",
          "That is too strong, and the analysis's own top variant row refutes",
          "it. The second version claimed the recoveries split by BIOMARKER",
          "class -- substitution-keyed regimens recovered, gene-keyed ones not.",
          "That was wrong too, and it was wrong because two panel labels were",
          "mine to assign and I assigned them in the direction that made the",
          "split clean.", "",
          f"Measured over {rec['of']} regimens with a molecular biomarker "
          f"({rec['of'] - 3} approved combinations and 3 major published trials "
          "without one: NEJ009, AG221-005, KRYSTAL-7), with FLAURA2 relabelled",
          "to the substitution it actually enrols (the same one as MARIPOSA)",
          "and CAPItello-291 relabelled to the multi-gene alteration list it",
          "actually uses:", "",
          "| | regimen | gene | biomarker | both targeted | papers | trial |",
          "|---|---|---|---|---|--:|---|"]
    for x in rec["panel"]:
        L.append(f"| {'found' if x['recovered'] else '**MISSED**'} | "
                 f"{x['regimen']} | {x['gene']} | `{x['variant']}` | "
                 f"{'yes' if x['both_drugs_are_targeted'] else 'no'} | "
                 f"{x['papers']} | {x['trial']} |")
    bio, par = rec["by_biomarker_class"], rec["by_partner_class"]
    L += ["", "Two explanations, one of which does not survive contact with the",
          "panel:", "",
          "| explanation | holds | does not hold | separation |",
          "|---|--:|--:|--:|",
          f"| the biomarker is a substitution | {bio['yes_recovered']}/"
          f"{bio['yes_total']} | {bio['no_recovered']}/{bio['no_total']} | "
          f"{bio['separation']:+.2f} |",
          f"| **both drugs are targeted agents** | **{par['yes_recovered']}/"
          f"{par['yes_total']}** | **{par['no_recovered']}/{par['no_total']}** | "
          f"**{par['separation']:+.2f}** |", "",
          ("**Partner-agent class predicts the recoveries and biomarker class"
           " does not.**" if rec["partner_class_explains_better"] else
           "**Biomarker class predicts the recoveries at least as well.**"), "",
          "### But not for the reason an earlier version gave", "",
          "That version asserted the extractor *cannot* tie a chemotherapy, a",
          "hypomethylating agent, a checkpoint antibody or an endocrine agent to",
          "a substitution. Nothing measured that, and measuring it refutes it.",
          "For each missed regimen, how many papers tie each drug to the target",
          "variant ANYWHERE, beside the zero that tie both in one paper:", "",
          "| trial | targeted partner | other partner | both, never together |",
          "|---|--:|--:|---|"]
    for x in rec["panel"]:
        if x["recovered"]:
            continue
        A, B = sorted((x.get("drug_a_tied_anywhere", 0),
                       x.get("drug_b_tied_anywhere", 0)), reverse=True)
        L.append(f"| {x['trial']} | {A} | {B} | "
                 f"{'yes' if x.get('both_tied_but_never_together') else 'no'} |")
    L += ["",
          "Pemetrexed IS tied to EGFR L858R, azacitidine to IDH1 R132H,",
          "pembrolizumab to KRAS G12C. The extractor does it. What is thin is",
          "the VOLUME on one side: the targeted partner is tied to the variant",
          "in tens or hundreds of papers and the other in nought to a handful,",
          "so a single paper carrying BOTH links is rare rather than",
          "impossible.", "",
          "So the honest causal statement is weaker than the previous one and",
          "is what the join actually requires: **both drugs must be tied to the",
          "same variant in ONE paper, and for a targeted-plus-other regimen the",
          "second link is too rare for that to happen.** Partner class predicts",
          "the outcome; it is not the mechanism.", "",
          "So: **an absence is evidence about targeted-plus-targeted regimens",
          "and says nothing about a regimen with a chemotherapy, endocrine or",
          "immunotherapy partner.** The gene layer is reported beside the",
          "variant layer because it recovers several of those misses.", "",
          "## The two layers", "",
          "| | count |", "|---|--:|",
          f"| papers asserting a co-treatment | {r['papers_asserting_a_cotreatment']:,} |",
          f"| papers asserting a variant relation | "
          f"{r['papers_asserting_a_variant_relation']:,} |",
          f"| **both** | **{r['papers_asserting_both']:,}** "
          f"({r['overlap_share_of_cotreatment_papers']}%) |",
          f"| gene-level triples | {g['triples']:,} |",
          f"| ...resting on one paper | {g['single_paper']:,} |",
          *[f"| ...{k} | {v:,} |" for k, v in sorted(g["by_verdict"].items())],
          f"| variant-level triples | {v['triples']:,} |", "",
          "| drug | drug | gene | papers | cotreat | compare | ratio | verdict |",
          "|---|---|---|--:|--:|--:|--:|---|"]
    for row in g["top"][:25]:
        rr = row["cotreat_compare_ratio"]
        L.append(f"| {row['drug_a']} | {row['drug_b']} | {row['target']} | "
                 f"{row['papers']} | {row['cotreat']} | {row['compare']} | "
                 f"{'--' if rr is None else rr} | {row['verdict']} |")
    L += ["", "### Variant level", "",
          "| drug | drug | gene | variant | papers | ratio |",
          "|---|---|---|---|--:|--:|"]
    for row in v["top"][:12]:
        L.append(f"| {row['drug_a']} | {row['drug_b']} | {row['target']} | "
                 f"`{row['variant']}` | {row['papers']} | "
                 f"{row['cotreat_compare_ratio']} |")
    L += ["", "## What neither layer can say", "",
          "* **No outcome.** Nothing here says a combination worked.",
          "* **A gene is not a target.** Docetaxel + estramustine ties to KLK3,",
          "  which is PSA: a response biomarker in prostate cancer.",
          "* **Gene-symbol collisions pass straight through.** This reads",
          "  PubTator's own gene assignment, so the `atlas_ambiguity` blocklist,",
          "  which guards the alias-resolution path in `atlas_graph.resolve`,",
          "  does not apply. Each collision below is CHECKED against the full",
          "  alias list at render time and dropped if it does not hold, because",
          "  a hand-written one did not: an earlier version claimed `NP` reaches",
          "  Neptunium *\"verified in the authority table itself\"*, and",
          "  Neptunium's alias list contains no `NP`.",
          *[f"    * `{x['alias']}` resolves to **{x['resolves_to']}** "
            f"({x['note']})" for x in r["collisions"]],
          "* **The gene layer's precision is not measured, at all.** No sample",
          "  has been judged and no precision figure is claimed here. An earlier",
          "  version of this section reported a 40-row sample that this",
          f"  repository never ran. {g['single_paper']:,} of {g['triples']:,}",
          "  rows rest on a single paper, and the collision classes above are",
          "  the known error source, but the rate is unknown.",
          "* **Full table** is `analysis/atlas-combination-gaps.tsv.gz`, both",
          "  layers, with the ratio column.", "",
          "## The transferable part", "",
          "Two of the three things this document reports began as assumptions",
          "that survived until something checked them: that a predicate means",
          "what it is named, and that one failed control licenses a general",
          "claim. Both were wrong, and both were cheap to test against cases",
          "whose answers are known independently of the corpus. This repository",
          "has recorded the same shape before: `analysis/atlas-replication.md`'s",
          "replication 'collapse' was a censoring artifact of the author's own",
          "window, and `scripts/calibration_feasibility.py` found a prevalence",
          "prior that was the normal distribution being recovered rather than",
          "biology.", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
