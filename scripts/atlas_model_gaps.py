#!/usr/bin/env python3
"""Atlas: which ferroptosis genes does the engine not model? (#ATLAS-GAPS)

WHY
---
`ferroptosis-core` has grown roughly thirty off-by-default layers, each added
because the author read a paper. That is a reading list, not a priority order,
and nothing has ever asked the literature which mechanisms it actually
emphasises.

The census can. 13,346 cancer articles carry the MeSH `Ferroptosis` descriptor,
and PubTator annotates their genes, so the field's own attention is countable
and comparable against what the engine has a handle for.

THREE THINGS THAT MUST BE HANDLED OR THE ANSWER IS NONSENSE
------------------------------------------------------------
ORTHOLOGS. GPX4 appears as separate human, mouse and rat identifiers. Ranked
raw, GPX4 occupies three of the top eleven positions and two of them look
unmodelled. Species ambiguity is the benign class from
`analysis/atlas-ambiguity.md` and a prevalence question wants them merged.

METHOD GENES. GAPDH and beta-actin rank 11th and 17th. They are Western-blot
loading controls, and caspase-3 is measured in ferroptosis papers precisely to
show the death is NOT apoptotic. Reporting these as modelling gaps would be
absurd, so they are excluded by name with the reason recorded.

WHAT "MODELLED" MEANS. The 20 module CLAIMS in `atlas_module_support.py` are a
hand-written list for a different analysis, not a statement of engine coverage.
Comparing against it reports NRF2 and NCOA4 as gaps when `params.rs` carries
`nrf2_gsh_rate` and `ferritinophagy_release`. Coverage is read from the
parameter and module surface instead.

Usage:
    python scripts/atlas_model_gaps.py
"""

import collections
import glob
import gzip
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_graph import load_index  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-model-gaps.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-model-gaps.json"
CORE = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src"

# Genes that appear because of how experiments are REPORTED, not what they are
# about. Each is excluded with its reason rather than silently dropped.
METHOD_GENES = {
    "GAPDH": "Western-blot loading control",
    "BETA-ACTIN": "Western-blot loading control",
    "CASPASE-3": "apoptosis marker, measured to show the death is NOT apoptotic",
    "CD8": "immune-cell marker rather than a ferroptosis mechanism",
    "KI-67": "proliferation marker",
}

# Engine coverage: parameter or module name -> the gene concept it implements.
# Read from params.rs field names and the module list, because the CLAIMS list
# is not a coverage statement.
COVERAGE = {
    "GPX4": "gpx4_rate", "SLC7A11": "erastin_xc_inhib", "ACSL4": "acsl4_status_boost",
    "AIFM2": "fsp1_rate", "FSP1": "fsp1_rate", "DHODH": "dhodh_rate",
    "GCH1": "gch1_rate", "NCOA4": "ferritinophagy_release", "PROM2": "prom2_iron_efflux",
    "FTH1": "ferritinophagy_release", "NRF2": "nrf2_gsh_rate", "NFE2L2": "nrf2_gsh_rate",
    "ALOX15": "alox_propagation_boost", "SCD": "scd_mufa_rate", "MBOAT2": "mboat_mufa_boost",
    "POR": "por_h2o2_rate", "DHCR7": "dhc7_radical_trap", "VKORC1L1": "vitk_radical_trap",
    "IFN-GAMMA": "ifngamma module", "IFNG": "ifngamma module",
    "SLC3A2": "erastin_xc_inhib", "CDH1": "contact module", "YAP1": "contact module",
    "ATP7A": "copper module", "ATP7B": "copper module", "CHMP5": "escrt_repair_rate",
    "CHMP6": "escrt_repair_rate", "HDAC1": "hdac_inhibitor", "CD274": "pd1_brake",
    "PDCD1": "pd1_brake", "FAR1": "ether_pufa_fraction", "AGPS": "ether_pufa_fraction",
}
ALIAS = {"NRF2": "NFE2L2", "HO-1": "HMOX1", "TRANSFERRIN RECEPTOR": "TFRC",
         "P53": "TP53", "COX-2": "PTGS2", "CYCLOOXYGENASE-2": "PTGS2",
         "FSP1": "AIFM2", "HIF-1ALPHA": "HIF1A"}


def main() -> int:
    root = atlas_root()
    idx = load_index(root)

    print("collecting ferroptosis-indexed articles ...", flush=True)
    ferro = set()
    for f in sorted(glob.glob(str(root / "records" / "*.jsonl.gz"))):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                if "Ferroptosis" in (r.get("mesh") or []):
                    ferro.add(r["pmid"])
    print(f"  {len(ferro):,} articles", flush=True)

    raw = collections.Counter()
    with gzip.open(root / "entities" / "gene.tsv.gz", "rt", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[0] in ferro:
                raw[p[2]] += 1

    # collapse orthologs by canonical name
    merged = collections.Counter()
    for gid, c in raw.items():
        merged[(idx["canon"].get(gid, gid) or gid).upper()] += c

    norm = lambda n: ALIAS.get(n, n)  # noqa: E731
    covered = {norm(k) for k in COVERAGE}
    rows, excluded = [], []
    for name, c in merged.most_common(60):
        n = norm(name)
        if n in METHOD_GENES or name in METHOD_GENES:
            excluded.append({"gene": name, "papers": c,
                             "reason": METHOD_GENES.get(name) or METHOD_GENES.get(n)})
            continue
        rows.append({"gene": name, "normalised": n, "papers": c,
                     "modelled": n in covered,
                     "handle": COVERAGE.get(name) or COVERAGE.get(n)})
    gaps = [r for r in rows if not r["modelled"]][:12]

    L = [
        "# Which ferroptosis genes does the engine not model? (#ATLAS-GAPS)", "",
        "Generated by `scripts/atlas_model_gaps.py`.", "",
        "`ferroptosis-core` has around thirty off-by-default layers, each added because",
        "the author read a paper. That is a reading list, not a priority order, and",
        "nothing had asked the literature which mechanisms it actually emphasises.",
        f"{len(ferro):,} cancer articles carry the MeSH `Ferroptosis` descriptor and",
        "PubTator annotates their genes, so the field's attention is countable.", "",
        "## The gaps, by how much the literature studies them", "",
        "| gene | articles | engine handle |", "|---|---|---|",
    ]
    for r in gaps:
        L.append(f"| **{r['gene']}** | {r['papers']:,} | none |")

    L += [
        "", "## What the engine already covers", "",
        "| gene | articles | handle |", "|---|---|---|",
    ]
    for r in [x for x in rows if x["modelled"]][:14]:
        L.append(f"| {r['gene']} | {r['papers']:,} | `{r['handle']}` |")

    L += [
        "", "## Three corrections without which this list is nonsense", "",
        "**Orthologs.** GPX4 appears as separate human, mouse and rat identifiers.",
        "Ranked raw it occupies three of the top eleven positions and two of them",
        "read as unmodelled. Species ambiguity is the benign class measured in",
        "`analysis/atlas-ambiguity.md`, and a prevalence question wants them merged.", "",
        "**Method genes.** These rank high because of how experiments are reported,",
        "not what they are about, and are excluded by name:", "",
        "| gene | articles | why excluded |", "|---|---|---|",
    ] + [
        f"| {e['gene']} | {e['papers']:,} | {e['reason']} |" for e in excluded
    ] + [
        "",
        "Caspase-3 is the instructive one: ferroptosis papers measure it precisely to",
        "show the death is NOT apoptotic, so its rank is evidence of the field's",
        "controls rather than of a mechanism worth modelling.", "",
        "**What 'modelled' means.** The 20 module claims in",
        "`atlas_module_support.py` are a hand-written list for a different analysis,",
        "not a statement of engine coverage. Comparing against it reports NRF2 and",
        "NCOA4 as gaps when `params.rs` carries `nrf2_gsh_rate` and",
        "`ferritinophagy_release`. Coverage is read from the parameter and module",
        "surface instead.", "",
        "## How to use this", "",
        "* It ranks ATTENTION, not importance. A heavily studied gene may be heavily",
        "  studied because it is easy to assay.",
        "* A gap is not a defect. The engine is a ferroptosis model, not a cancer",
        "  model, and several of these are upstream regulators whose effect the engine",
        "  already absorbs into an existing parameter.",
        "* The repo's own layer-freeze policy (CONTRIBUTING.md) requires a named",
        "  calibration target before a new axis lands. Nothing here bypasses that;",
        "  this only says which axes the literature would support calibrating.", "",
        "## What happened when the top four were checked (#616)", "",
        "None of them became a layer, because no calibration target exists. The route",
        "that partially anchored ACSL4 (cBioPortal TCGA within-cohort z-scores, #462)",
        "turns out to recover the normal distribution for every gene tested -- HMOX1",
        "14.9%, TP53 14.1%, TFRC 15.2%, KEAP1 13.8% below z = -1, against the 15.87%",
        "a normal gives for anything -- so it cannot distinguish a gene whose",
        "low-expression tail matters from one with no ferroptosis role. Two of the",
        "four are also already absorbed by an existing parameter (TFRC by",
        "`ferritinophagy_release`, KEAP1 by `nrf2_gsh_rate`, which it regulates).", "",
        "The one signal that survived is TP53's DEEP tail (4.3% below z = -2 against",
        "an expected 2.28%, with much the widest inter-cancer spread), consistent with",
        "real deletion rather than standardisation. See",
        "`analysis/calibration/calibration-feasibility.md`, and the",
        "'Layers proposed and NOT built' rows in `CALIBRATION_STATUS.md`.", "",
        "So read this table as a map of where the literature's attention and the",
        "available data do NOT currently overlap, rather than as a backlog.", "",
        "## Limits", "",
        "* Gene mention is not mechanism. An article annotated for a gene may mention",
        "  it once in a discussion.",
        "* PubTator annotates abstracts, so a gene central to a paper's Results but",
        "  absent from its abstract is undercounted.",
        "* The coverage map is hand-written from parameter names and will drift as the",
        "  engine changes; it is the same class of hand-curated artifact this project",
        "  has twice found to be wrong, so treat a single row with suspicion.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "ferroptosis_articles": len(ferro), "rows": rows,
        "excluded_method_genes": excluded, "gaps": gaps,
    }, indent=2) + "\n")
    print("\ntop unmodelled: " + ", ".join(f"{r['gene']}({r['papers']:,})" for r in gaps[:6]))
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
