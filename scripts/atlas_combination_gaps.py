#!/usr/bin/env python3
"""Which drug combinations does the literature study, and against what?

THE QUESTION THIS SET OUT TO ANSWER, AND WHY IT CANNOT BE ASKED
---------------------------------------------------------------
`relations.tsv.gz` carries 224,146 `cotreat` rows over 40,878 distinct drug
pairs, and nothing in this repository had read them. Joined to the variant
layer (#ATLAS-VARIANT) the obvious question is the therapeutically interesting
one: WHICH RESISTANCE MUTATIONS HAVE A COMBINATION ANSWER AND WHICH DO NOT.

Asked at the variant level it produces a clean, plausible and WRONG table.
BRAF V600E shows hundreds of combination papers; PIK3CA H1047R shows zero. The
second reads as a therapeutic gap. It is not one: alpelisib plus fulvestrant is
FDA-approved for exactly that population on the strength of SOLAR-1, and this
same corpus holds 59 papers asserting that co-treatment.

THE POSITIVE CONTROL IS THE METHOD, not a footnote to it. Before believing any
absence, run the measurement against a case where the answer is known. Here it
fails outright: of those 59 papers, 45 annotate PIK3CA as a GENE and NOT ONE
annotates a variant entity. The combination literature writes "PIK3CA-mutated",
not "H1047R". So a variant-level gap measures whether a paper happened to name
a specific substitution, and reports the field's most successful biomarker-
directed combinations as missing.

The scale of the mismatch is the real finding: 132,520 papers assert a
co-treatment and 68,420 assert a variant relation, but only 2,605 do both, so
the join has 2.0% of the co-treatment literature to work with.

WHAT WORKS INSTEAD
------------------
The same join at GENE level. It passes the control (alpelisib + fulvestrant ->
PIK3CA, 21 papers) and its top rows are named trial regimens rather than
plausible-looking noise: trastuzumab + pertuzumab -> ERBB2 (CLEOPATRA),
dabrafenib + trametinib -> BRAF, vemurafenib + cobimetinib -> BRAF (coBRIM),
encorafenib + cetuximab -> BRAF (BEACON), ipilimumab + nivolumab -> PDCD1 and
CTLA4, crizotinib + alectinib -> ALK.

Both layers are reported. The variant-level one is real where it fires and is
NOT to be read as coverage; the gene-level one is the one that can carry a
question about what the field studies.

WHAT NEITHER CAN SAY
--------------------
* `cotreat` records that two drugs were given together, NOT that the
  combination worked, and carries no direction or outcome.
* A GENE IS NOT A TARGET. Docetaxel + estramustine ties to KLK3, which is PSA:
  a response biomarker in prostate cancer, not the drug's target. The join
  says the paper connects both drugs to that gene, nothing more.
* It reads PubTator's own gene assignment, so `atlas_ambiguity`'s blocklist --
  which guards the alias-resolution path in `atlas_graph.resolve` -- does not
  apply here. `dabrafenib + trametinib -> MAP2K7` is almost certainly "MEK"
  resolving to MKK7 rather than to MAP2K1/2, and is flagged, not fixed.
* Most triples rest on a single paper, the same shape found across the whole
  graph.

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

# The control. Alpelisib + fulvestrant is approved for PIK3CA-mutated advanced
# breast cancer (SOLAR-1), so a measurement that reports no combination for
# PIK3CA is measuring its own coverage, not the field's.
#
# It is chosen because the answer is known INDEPENDENTLY of this corpus, and
# pinned by identifier rather than by name so a label change cannot silently
# turn the control off.
CONTROL = {
    "name": "alpelisib + fulvestrant, for PIK3CA-mutated breast cancer (SOLAR-1)",
    "drugs": ("MESH:C585539", "MESH:D000077267"),
    "gene": "5290",
    "why": ("FDA-approved for this population, so an analysis that finds no "
            "combination for PIK3CA has measured its own coverage"),
}


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
    """One pass. Co-treatment pairs, and what each drug is tied to, per paper."""
    cotreat = collections.defaultdict(set)
    cotreat_rows = 0
    drug_gene = collections.defaultdict(lambda: collections.defaultdict(set))
    drug_var = collections.defaultdict(lambda: collections.defaultdict(set))
    has_var = set()
    genes_in = collections.defaultdict(set)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, pred = parts[0], parts[1]
            ta, _, ia = parts[2].partition("|")
            tb, _, ib = parts[3].partition("|")
            if pred == "cotreat" and ta == "Chemical" and tb == "Chemical":
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
            "drug_gene": drug_gene, "drug_var": drug_var,
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


def main() -> int:
    lab = load_labels()

    def nm(i, prefix=""):
        n = lab.get(i)
        if n is None:
            return f"{prefix}{i}" if i else "?"
        return n[len("[OBSOLETE] "):] if n.startswith("[OBSOLETE] ") else n

    s = scan(RELATIONS)
    by_gene = join(s["cotreat"], s["drug_gene"])
    by_variant = join(s["cotreat"], s["drug_var"])

    # --- the control, run on BOTH layers -------------------------------------
    a, b = CONTROL["drugs"]
    pair = frozenset((a, b))
    control_papers = {p for p, prs in s["cotreat"].items() if pair in prs}
    ctrl = {
        **{k: v for k, v in CONTROL.items() if k != "drugs"},
        "drug_ids": list(CONTROL["drugs"]),
        "papers_asserting_the_combination": len(control_papers),
        "of_those_annotating_the_gene": sum(
            1 for p in control_papers if CONTROL["gene"] in s["genes_in"][p]),
        "of_those_annotating_any_variant": sum(
            1 for p in control_papers if p in s["has_var"]),
        "found_by_the_gene_join": len(by_gene.get((pair, CONTROL["gene"]), ())),
        "found_by_the_variant_join": sum(
            len(v) for (pr, t), v in by_variant.items()
            if pr == pair and t[0] == CONTROL["gene"]),
    }
    ctrl["gene_join_passes"] = ctrl["found_by_the_gene_join"] > 0
    ctrl["variant_join_passes"] = ctrl["found_by_the_variant_join"] > 0

    cot_papers, var_papers = set(s["cotreat"]), s["has_var"]
    overlap = len(cot_papers & var_papers)

    def table(joined, label):
        rows = [{"drug_a": nm(sorted(pr)[0]), "drug_b": nm(sorted(pr)[1]),
                 "drug_a_id": sorted(pr)[0], "drug_b_id": sorted(pr)[1],
                 "target": nm(t if isinstance(t, str) else t[0], "gene:"),
                 "target_id": t if isinstance(t, str) else t[0],
                 "variant": "" if isinstance(t, str) else t[1],
                 "papers": len(pm)}
                for (pr, t), pm in joined.items()]
        rows.sort(key=lambda r: (-r["papers"], r["drug_a_id"], r["drug_b_id"],
                                 r["target_id"], r["variant"]))
        return rows

    gene_rows, var_rows = table(by_gene, "gene"), table(by_variant, "variant")

    with open(OUT_TSV, "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as fh:
        fh.write("layer\tdrug_a\tdrug_b\ttarget\ttarget_id\tvariant\tpapers\n")
        for layer, rows in (("gene", gene_rows), ("variant", var_rows)):
            for r in rows:
                fh.write(f"{layer}\t{r['drug_a']}\t{r['drug_b']}\t{r['target']}\t"
                         f"{r['target_id']}\t{r['variant']}\t{r['papers']}\n")

    res = {
        # Rows and paper-pair combinations are DIFFERENT measures: a paper
        # asserting the same pair twice is two rows and one combination. The
        # first version labelled the second as the first.
        "cotreat_rows": s["cotreat_rows"],
        "paper_pair_combinations": sum(len(v) for v in s["cotreat"].values()),
        "distinct_drug_pairs": len({p for v in s["cotreat"].values() for p in v}),
        "papers_asserting_a_cotreatment": len(cot_papers),
        "papers_asserting_a_variant_relation": len(var_papers),
        "papers_asserting_both": overlap,
        "overlap_share_of_cotreatment_papers": round(
            100.0 * overlap / max(len(cot_papers), 1), 2),
        "control": ctrl,
        "gene_layer": {
            "triples": len(gene_rows),
            "single_paper": sum(1 for r in gene_rows if r["papers"] == 1),
            "at_least_three_papers": sum(1 for r in gene_rows if r["papers"] >= 3),
            "top": gene_rows[:40],
        },
        "variant_layer": {
            "triples": len(var_rows),
            "single_paper": sum(1 for r in var_rows if r["papers"] == 1),
            "at_least_three_papers": sum(1 for r in var_rows if r["papers"] >= 3),
            "top": var_rows[:25],
        },
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(res), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}\nwrote {OUT_TSV}")
    print(f"  gene layer {len(gene_rows):,} triples, variant layer {len(var_rows):,}")
    print(f"  control: gene join {'PASSES' if ctrl['gene_join_passes'] else 'FAILS'}, "
          f"variant join {'PASSES' if ctrl['variant_join_passes'] else 'FAILS'}")
    return 0


def render(r: dict) -> str:
    c, g, v = r["control"], r["gene_layer"], r["variant_layer"]
    L = [
        "# Which drug combinations the literature studies, and against what", "",
        "Generated by `scripts/atlas_combination_gaps.py`.", "",
        "## The question this set out to answer, and why it cannot be asked", "",
        f"The census carries {r['cotreat_rows']:,} `cotreat` rows over "
        f"{r['distinct_drug_pairs']:,} distinct drug pairs, and nothing here had",
        "read them. Joined to the variant layer, the obvious question is the",
        "therapeutically interesting one: **which resistance mutations have a",
        "combination answer and which do not.**", "",
        "Asked at the variant level it produces a clean, plausible and wrong",
        "table. BRAF V600E shows hundreds of combination papers; PIK3CA H1047R",
        "shows zero, which reads as a therapeutic gap.", "",
        "## The positive control", "",
        f"**{c['name']}.**", "",
        f"Chosen because {c['why']}, and because the answer is known",
        "INDEPENDENTLY of this corpus. Before believing any absence, the",
        "measurement is run against a case whose answer is not in doubt.", "",
        "| | |", "|---|--:|",
        f"| papers asserting the co-treatment | {c['papers_asserting_the_combination']} |",
        f"| ...annotating the gene | {c['of_those_annotating_the_gene']} |",
        f"| ...annotating ANY variant | **{c['of_those_annotating_any_variant']}** |",
        f"| found by the variant-level join | "
        f"**{c['found_by_the_variant_join']}** "
        f"({'passes' if c['variant_join_passes'] else 'FAILS'}) |",
        f"| found by the gene-level join | {c['found_by_the_gene_join']} "
        f"({'passes' if c['gene_join_passes'] else 'FAILS'}) |", "",
        "The variant-level join reports nothing for a combination that is",
        "approved for exactly this population. The reason is visible in the row",
        "above it: the combination literature writes *PIK3CA-mutated*, not",
        "*H1047R*. **So a variant-level gap measures whether a paper happened to",
        "name a specific substitution, and reports the field's most successful",
        "biomarker-directed combinations as missing.**", "",
        "## The scale of the mismatch", "",
        "| | count |", "|---|--:|",
        f"| papers asserting a co-treatment | {r['papers_asserting_a_cotreatment']:,} |",
        f"| papers asserting a variant relation | "
        f"{r['papers_asserting_a_variant_relation']:,} |",
        f"| **papers asserting both** | **{r['papers_asserting_both']:,}** |",
        f"| share of the co-treatment literature | "
        f"{r['overlap_share_of_cotreatment_papers']}% |", "",
        f"The join has {r['overlap_share_of_cotreatment_papers']}% of the",
        "co-treatment literature to work with, which is why an absence in it",
        "carries no information about the field.", "",
        "## What works instead: the same join at gene level", "",
        f"{g['triples']:,} (combination, gene) triples where both drugs are tied",
        f"to the same gene in the same paper; {g['at_least_three_papers']:,} rest",
        f"on three or more papers and {g['single_paper']:,} on one.", "",
        "It passes the control, and its top rows are named trial regimens rather",
        "than plausible-looking noise, which is the check that it is measuring",
        "something real.", "",
        "| drug | drug | gene | papers |", "|---|---|---|--:|"]
    for row in g["top"][:25]:
        L.append(f"| {row['drug_a']} | {row['drug_b']} | {row['target']} | "
                 f"{row['papers']} |")
    L += ["", "## The variant layer, reported for what it is", "",
          f"{v['triples']:,} (combination, gene, variant) triples, "
          f"{v['single_paper']:,} of them on a single paper. Real where it fires",
          "and NOT a coverage measure, per the control above.", "",
          "| drug | drug | gene | variant | papers |", "|---|---|---|---|--:|"]
    for row in v["top"][:15]:
        L.append(f"| {row['drug_a']} | {row['drug_b']} | {row['target']} | "
                 f"`{row['variant']}` | {row['papers']} |")
    L += ["", "## What neither layer can say", "",
          "* **`cotreat` is not efficacy.** It records that two drugs were given",
          "  together, not that the combination worked, and carries no direction",
          "  or outcome.",
          "* **A gene is not a target.** Docetaxel + estramustine ties to KLK3,",
          "  which is PSA: a response biomarker in prostate cancer, not the",
          "  drug's target. The join says the paper connects both drugs to that",
          "  gene and nothing more.",
          "* **Gene-symbol collisions pass straight through.** This reads",
          "  PubTator's own gene assignment, so the `atlas_ambiguity` blocklist,",
          "  which guards the alias-resolution path in `atlas_graph.resolve`,",
          "  does not apply. `dabrafenib + trametinib -> MAP2K7` is almost",
          "  certainly *MEK* resolving to MKK7 rather than MAP2K1/2. Flagged,",
          "  not fixed.",
          "* **Attention, not importance**, as everywhere else in this atlas.",
          f"* **Full table** is `analysis/atlas-combination-gaps.tsv.gz`; both",
          "  layers are in it, distinguished by the first column.", "",
          "## The transferable part", "",
          "The finding here is not the combination table. It is that **an",
          "absence measured by a join is a statement about the join until a",
          "positive control says otherwise** -- and that running one cost a few",
          "minutes and overturned the headline. This repository has made the",
          "opposite mistake before: a measured 24% performance deficit that",
          "turned out to be a configuration artifact, and a replication",
          "'collapse' that was a censoring artifact of the author's own window.", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
