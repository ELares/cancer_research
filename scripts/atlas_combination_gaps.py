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

2. AN ABSENCE IS INFORMATIVE ONLY ABOUT SUBSTITUTION-KEYED REGIMENS
-------------------------------------------------------------------
The first version of this analysis claimed a variant-level absence carries NO
information, generalising from one control that failed. That is too strong, and
the analysis's own top row refutes it: dabrafenib + trametinib against BRAF
V600E, 181 papers.

Measured against a panel of ten approved regimens, the variant join recovers
EIGHT. The two it misses are both the same shape: alpelisib + fulvestrant
(SOLAR-1) and osimertinib + chemotherapy (FLAURA2) are keyed to
"PIK3CA-mutated" and "EGFR-mutated", a GENE-level biomarker, not to a
substitution. The control makes the failure mode visible: of the 59 papers
asserting alpelisib + fulvestrant, 45 annotate PIK3CA as a gene and NOT ONE
annotates a variant.

So the honest statement is narrow and useful: an absence at variant level is
evidence about substitution-keyed regimens and says nothing about gene-keyed
ones, AND YOU CANNOT TELL WHICH YOU ARE LOOKING AT FROM THE ABSENCE ALONE. That
is why the gene-level join is reported beside it rather than instead of it.

WHAT NEITHER LAYER CAN SAY
--------------------------
* `cotreat` carries no outcome. Nothing here says a combination worked.
* A GENE IS NOT A TARGET. Docetaxel + estramustine ties to KLK3, which is PSA,
  a response biomarker.
* This reads PubTator's OWN gene assignment, so the `atlas_ambiguity` blocklist
  -- which guards the alias-resolution path in `atlas_graph.resolve` -- does not
  apply. Verified collisions in the authority table itself: `MEK` is an alias of
  MAP2K7 while MAP2K1 carries `MEK1`; `Met` is an alias of SLTM; `PGP` resolves
  to phosphoglycolate phosphatase rather than P-glycoprotein; `NP` reaches
  Neptunium. The gene layer's precision is NOT measured at scale here and a
  sampled estimate is reported as a bound, not a rate.
* Most gene-layer triples rest on one paper.

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

# NAMES, resolved through the authority table at runtime and asserted to
# resolve. An earlier version hardcoded MeSH identifiers from memory and got
# several wrong -- MESH:C571179 is eravacycline, not vemurafenib -- which
# produced a row reading MISSED. A false negative that looks like a finding is
# the defect class this whole document is about, so an unresolvable name is a
# hard failure rather than a miss.
SEMANTIC_PANEL = [
    ("trastuzumab", "pertuzumab", True),
    ("dabrafenib", "trametinib", True),
    ("encorafenib", "binimetinib", True),
    ("ipilimumab", "nivolumab", True),
    ("cisplatin", "etoposide", True),
    ("alectinib", "crizotinib", False),
    ("ceritinib", "crizotinib", False),
    ("osimertinib", "gefitinib", False),
    ("lorlatinib", "crizotinib", False),
]

# Approved regimens with a molecular biomarker, used to measure what the
# variant-level join RECOVERS rather than asserting that it recovers nothing.
# An empty `variant` means the regimen's biomarker is stated at gene level,
# which is the failure mode the misses share.
RECOVERY_PANEL = [
    ("dabrafenib", "trametinib", "673", "p.V600E", "COMBI-d"),
    ("encorafenib", "cetuximab", "673", "p.V600E", "BEACON"),
    ("encorafenib", "binimetinib", "673", "p.V600E", "COLUMBUS"),
    ("vemurafenib", "cobimetinib", "673", "p.V600E", "coBRIM"),
    ("amivantamab", "lazertinib", "1956", "p.L858R", "MARIPOSA"),
    ("sotorasib", "panitumumab", "3845", "p.G12C", "CodeBreaK 300"),
    ("adagrasib", "cetuximab", "3845", "p.G12C", "KRYSTAL-1"),
    ("capivasertib", "fulvestrant", "207", "p.E17K", "CAPItello-291"),
    ("alpelisib", "fulvestrant", "5290", "", "SOLAR-1"),
    ("osimertinib", "pemetrexed", "1956", "", "FLAURA2"),
]

# The ratio above which a pair reads as co-administered, taken from the gap in
# the panel (every true pair well above, every false pair well below).
# Reported, never used to filter: a nine-pair panel cannot support a classifier.
CO_ADMIN_RATIO = 2.0


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


def ratio(pair_pred: dict, pair) -> float:
    """cotreat rows per compare row. High means co-administration."""
    c = pair_pred.get(pair, {})
    return round(c.get("cotreat", 0) / max(c.get("compare", 0), 1), 2)


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
        semantic.append({"pair": f"{an} + {bn}", "co_administered": truth,
                         "cotreat": c.get("cotreat", 0),
                         "compare": c.get("compare", 0),
                         "ratio": ratio(pp, frozenset((a, b)))})
    lo_true = min((x["ratio"] for x in semantic if x["co_administered"]), default=0)
    hi_false = max((x["ratio"] for x in semantic if not x["co_administered"]),
                   default=0)

    # --- 2. what does the variant join recover? ------------------------------
    recovery = []
    for an, bn, gene, variant, trial in RECOVERY_PANEL:
        a, b = resolve_drug(an, lab), resolve_drug(bn, lab)
        pair = frozenset((a, b))
        n = sum(len(pm) for (pr, t), pm in by_variant.items()
                if pr == pair and t[0] == gene and (not variant or t[1] == variant))
        recovery.append({
            "trial": trial, "regimen": f"{nm(a)} + {nm(b)}",
            "drug_ids": sorted((a, b)),
            "gene": nm(gene, "gene:"), "variant": variant or "(gene-level)",
            "biomarker_is_a_substitution": bool(variant),
            "papers": n, "recovered": n > 0})
    subst = [r for r in recovery if r["biomarker_is_a_substitution"]]
    gene_keyed = [r for r in recovery if not r["biomarker_is_a_substitution"]]

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
            rr = ratio(pp, pr)
            rows.append({
                "drug_a": nm(a), "drug_b": nm(b), "drug_a_id": a, "drug_b_id": b,
                "target": nm(t[0] if is_variant else t, "gene:"),
                "target_id": t[0] if is_variant else t,
                "variant": t[1] if is_variant else "",
                "papers": len(pm), "cotreat_compare_ratio": rr,
                "reads_as_co_administration": rr >= CO_ADMIN_RATIO})
        rows.sort(key=lambda x: (-x["papers"], x["drug_a_id"], x["drug_b_id"],
                                 x["target_id"], x["variant"]))
        return rows

    gene_rows, var_rows = table(by_gene, False), table(by_variant, True)

    with open(OUT_TSV, "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as fh:
        fh.write("layer\tdrug_a\tdrug_b\ttarget\ttarget_id\tvariant\tpapers"
                 "\tcotreat_compare_ratio\treads_as_co_administration\n")
        for layer, rows in (("gene", gene_rows), ("variant", var_rows)):
            for r in rows:
                fh.write(f"{layer}\t{r['drug_a']}\t{r['drug_b']}\t{r['target']}\t"
                         f"{r['target_id']}\t{r['variant']}\t{r['papers']}\t"
                         f"{r['cotreat_compare_ratio']}\t"
                         f"{int(r['reads_as_co_administration'])}\n")

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
            "substitution_keyed_recovered": sum(1 for r in subst if r["recovered"]),
            "substitution_keyed_total": len(subst),
            "gene_keyed_recovered": sum(1 for r in gene_keyed if r["recovered"]),
            "gene_keyed_total": len(gene_keyed),
        },
        "control": control,
        "gene_layer": {
            "triples": len(gene_rows),
            "single_paper": sum(1 for r in gene_rows if r["papers"] == 1),
            "at_least_three_papers": sum(1 for r in gene_rows if r["papers"] >= 3),
            "reading_as_co_administration": sum(
                1 for r in gene_rows if r["reads_as_co_administration"]),
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
          "## 2. An absence is informative only about substitution-keyed regimens",
          "",
          "An earlier version of this analysis claimed a variant-level absence",
          "carries NO information, generalising from one control that failed.",
          "That is too strong, and the analysis's own top variant row refutes it.",
          "",
          f"Measured against {rec['of']} approved regimens, the variant join",
          f"recovers **{rec['recovered']}**:", "",
          "| | regimen | gene | biomarker | papers | trial |",
          "|---|---|---|---|--:|---|"]
    for x in rec["panel"]:
        L.append(f"| {'found' if x['recovered'] else '**MISSED**'} | "
                 f"{x['regimen']} | {x['gene']} | `{x['variant']}` | "
                 f"{x['papers']} | {x['trial']} |")
    L += ["",
          f"The split is the finding: **{rec['substitution_keyed_recovered']} of "
          f"{rec['substitution_keyed_total']}** regimens whose biomarker IS a",
          f"substitution are recovered, against **{rec['gene_keyed_recovered']} "
          f"of {rec['gene_keyed_total']}** whose biomarker is stated at gene",
          "level.", "",
          "The control makes that failure mode visible. Of the",
          f"{c['papers_asserting_the_combination']} papers asserting",
          f"{c['regimen']}, **{c['of_those_annotating_the_gene']}** annotate the",
          f"gene and **{c['of_those_annotating_any_variant']}** annotate any",
          "variant. The literature writes *PIK3CA-mutated*, not *H1047R*.", "",
          "So the honest statement is narrow: **an absence at variant level is",
          "evidence about substitution-keyed regimens and says nothing about",
          "gene-keyed ones, and you cannot tell which you are looking at from",
          "the absence alone.** That is why the gene layer is reported beside it",
          "rather than instead of it.", "",
          "## The two layers", "",
          "| | count |", "|---|--:|",
          f"| papers asserting a co-treatment | {r['papers_asserting_a_cotreatment']:,} |",
          f"| papers asserting a variant relation | "
          f"{r['papers_asserting_a_variant_relation']:,} |",
          f"| **both** | **{r['papers_asserting_both']:,}** "
          f"({r['overlap_share_of_cotreatment_papers']}%) |",
          f"| gene-level triples | {g['triples']:,} |",
          f"| ...resting on one paper | {g['single_paper']:,} |",
          f"| ...reading as co-administration | {g['reading_as_co_administration']:,} |",
          f"| variant-level triples | {v['triples']:,} |", "",
          "| drug | drug | gene | papers | ratio | |", "|---|---|---|--:|--:|---|"]
    for row in g["top"][:25]:
        flag = "" if row["reads_as_co_administration"] else " **compare-dominated**"
        L.append(f"| {row['drug_a']} | {row['drug_b']} | {row['target']} | "
                 f"{row['papers']} | {row['cotreat_compare_ratio']} |{flag} |")
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
          "  does not apply. Verified in the authority table itself: `MEK` is an",
          "  alias of **MAP2K7** while MAP2K1 carries `MEK1`; `Met` is an alias",
          "  of **SLTM**; `PGP` reaches phosphoglycolate phosphatase rather than",
          "  P-glycoprotein; `NP` reaches **Neptunium**.",
          "* **The gene layer's precision is not measured at scale.** A 40-row",
          "  sample judged from title and abstract put roughly a third of",
          "  single-paper rows in a clearly-wrong class, mostly from the",
          "  collisions above. Treat that as a bound on the tail rather than a",
          f"  rate, and note that {g['single_paper']:,} of {g['triples']:,} rows",
          "  are single-paper.",
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
