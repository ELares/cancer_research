#!/usr/bin/env python3
"""Atlas: corpus support for each simulation module's central claim (#ATLAS).

WHY
---
`ferroptosis-core` carries ~30 off-by-default "realism layers". Each was added
on the strength of one or two papers the author read, cited by PMID in the
module docs. That is the weakest form of grounding available: it cannot
distinguish a mechanism the field has replicated from one asserted once.

The atlas supplies the denominator. For each module, this asks the 7,951,325
typed relations over 1,603,105 cancer articles:

  * how many DISTINCT articles assert a relation between the module's two
    entities, under which predicates; and
  * whether the module's own cited PMID is among them.

That second check matters both ways. If the cited paper is present, the module's
grounding is confirmed to be in the machine-readable record. If it is absent
while other papers assert the same edge, the claim is corroborated but the
citation is not machine-visible. If nothing is found at all, the claim rests on
exactly the one reading the author did.

WHAT THIS CANNOT DO
-------------------
PubTator3's entity types are Gene, Chemical, Disease, Species, CellLine and
Variant. There is NO biological-process type, so `ferroptosis`, `apoptosis` and
`autophagy` are not first-class entities and "X induces ferroptosis" is not
expressible. Every claim below is therefore reduced to a gene-gene or
gene-chemical pair, which is a genuine loss of specificity: the graph can say
IFN-gamma and SLC7A11 are negatively correlated across 8 papers, not that this
sensitises cells to ferroptosis.

Counts are also not truth. Roughly half of all relations are the weakest
predicate, `associate`, closer to co-mention than to knowledge, and PubTator's
extractor scores ~79.6 F1 on BioRED. A count means the field discusses the pair.

Usage:
    python scripts/atlas_module_support.py
"""

import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_graph import load_index, support  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-module-support.md"

# (module, entity A, entity B, PMID cited in the module docs, one-line claim)
# Pairs are the closest gene/chemical proxy for each module's mechanism, since
# the process itself is not an addressable entity.
CLAIMS = [
    ("ifngamma", "IFNG", "SLC7A11", "31043744",
     "IFN-gamma represses System Xc- (SLC7A11), starving cystine"),
    ("ifngamma", "IFNG", "SLC3A2", "31043744",
     "IFN-gamma represses the System Xc- heavy chain SLC3A2"),
    ("acsl4", "ACSL4", "GPX4", "27842070",
     "ACSL4 sets the oxidisable-PUFA substrate GPX4 must defend"),
    ("fsp1", "AIFM2", "GPX4", "31634899",
     "FSP1/AIFM2 is the GPX4-independent parallel defence"),
    ("dhodh", "DHODH", "GPX4", "33864050",
     "DHODH repairs lipid peroxides in parallel to GPX4"),
    ("gch1", "GCH1", "GPX4", "31919077",
     "GCH1/BH4 is a GPX4-independent radical-trapping defence"),
    ("prom2", "PROM2", "FTH1", "31761539",
     "PROM2 exports ferritin-bound iron, draining the labile pool"),
    ("vitk", "VKORC1L1", "GPX4", "37467745",
     "VKORC1L1 reduces vitamin K to a GPX4-independent radical trap"),
    ("copper", "ATP7A", "GPX4", "34390123",
     "copper ionophores degrade ATP7A and deplete GSH/GPX4"),
    ("alox", "ALOX15", "ACSL4", "27506793",
     "ALOX15 enzymatically peroxidises the ACSL4-supplied PUFA"),
    ("por", "POR", "CYB5R1", "33860083",
     "POR and CYB5R1 generate the H2O2 the Fenton reaction needs"),
    ("dhc7", "DHCR7", "GPX4", "38297130",
     "7-DHC (consumed by DHCR7) is an endogenous radical trap"),
    ("repair", "CHMP5", "CHMP6", "31761326",
     "ESCRT-III membrane repair blocks death execution"),
    ("contact", "CDH1", "YAP1", "31341276",
     "E-cadherin -> NF2/YAP suppresses ACSL4/TFRC in dense cells"),
    ("mboat", "MBOAT2", "GPX4", "37267948",
     "MBOAT1/2 remodel phospholipids toward MUFA-PE, GPX4-independently"),
    ("ether_lipid", "FAR1", "AGPS", "32939090",
     "FAR1/AGPS make the ether-PUFA pool that promotes ferroptosis"),
    ("dc_ferroptosis", "CD274", "SLC7A11", "39423128",
     "PD-L1 loss on DCs downregulates System Xc-, killing the DCs"),
    ("system_xc", "SLC7A11", "GPX4", "22632970",
     "System Xc- supplies the cystine that sustains GSH for GPX4"),
    ("erastin", "erastin", "SLC7A11", "22632970",
     "erastin inhibits System Xc-"),
    ("hdac_persister", "HDAC1", "AIFM2", "41481741",
     "HDACs and FSP1 together suppress persister-cell ferroptosis"),
]


def load_comentions(root: Path) -> dict:
    """Sentence-level co-mention counts from full text, if that layer has been built.

    The PubTator relation graph is abstract-level and its edge recall is low, so a
    zero there is uninformative. Full-text co-mention is the recall complement:
    no predicate, no direction, but it answers "does the literature discuss this
    pair at all", which is the question a zero in the relation column raises.
    """
    f = root / "comention" / "pairs.tsv.gz"
    if not f.exists():
        return {}
    out = {}
    with gzip.open(f, "rt", encoding="utf-8") as fh:
        for line in fh:
            a, b, c = line.rstrip("\n").split("\t")
            out[(a, b)] = int(c)
    return out


def main() -> None:
    root = atlas_root()
    idx = load_index(root)
    comention = load_comentions(root)
    if comention:
        print(f"co-mention layer: {len(comention):,} pairs", flush=True)
    rows = []
    for module, a, b, pmid, claim in CLAIMS:
        r = support(idx, a, b)
        if r is None:
            rows.append(dict(module=module, a=a, b=b, pmid=pmid, claim=claim,
                             resolved=False, total=0, preds={}, cited_present=None))
            continue
        pos = r["predicates"].get("positive_correlate", 0)
        neg = r["predicates"].get("negative_correlate", 0)
        rows.append(dict(module=module, a=a, b=b, pmid=pmid, claim=claim,
                         resolved=True, total=r["total"], preds=r["predicates"],
                         a_name=r["a_name"], b_name=r["b_name"],
                         cited_present=(pmid in r["pmids"]) if r["pmids"] else None,
                         pmids=r["pmids"][:8],
                         pos=pos, neg=neg,
                         contested=bool(pos and neg),
                         balance=(min(pos, neg) / max(pos, neg)) if (pos and neg) else None,
                         comention=comention.get(tuple(sorted((r["a"], r["b"]))))))

    found = [r for r in rows if r["resolved"] and r["total"] > 0]
    none_ = [r for r in rows if r["resolved"] and r["total"] == 0]
    unres = [r for r in rows if not r["resolved"]]

    L = [
        "# Corpus support for each simulation module's central claim (#ATLAS)", "",
        "Generated by `scripts/atlas_module_support.py` against the atlas relation",
        "graph: 7,951,325 typed relations over 1,603,105 cancer-article PMIDs.", "",
        "## Why", "",
        "Each `ferroptosis-core` realism layer was added on the strength of one or two",
        "papers, cited by PMID in the module docs. This asks how many DISTINCT cancer",
        "articles assert the same entity relation, and whether the cited paper is among",
        "them.", "",
        "## What this cannot do -- read before using the numbers", "",
        "PubTator3 has entity types Gene, Chemical, Disease, Species, CellLine and",
        "Variant. There is **no biological-process type**, so `ferroptosis` is not an",
        "addressable entity and \"X induces ferroptosis\" cannot be expressed. Every claim",
        "is reduced to a gene-gene or gene-chemical pair, which loses the mechanism: the",
        "graph can say IFN-gamma and SLC7A11 are negatively correlated across 8 papers,",
        "not that this sensitises cells to ferroptosis.", "",
        "Counts are not truth. About half of all relations carry the weakest predicate,",
        "`associate`, which is nearer co-mention than knowledge, and the extractor scores",
        "~79.6 F1 on BioRED. A high count means the field discusses the pair.", "",
        f"## Result: {len(found)}/{len(rows)} claims have corpus support", "",
        "> **The denominator is hand-made.** These "
        f"{len(rows)} are author-written claims with",
        "> author-chosen proxy entity pairs, covering 19 of roughly 30 library modules",
        "> (`CLAIMS` in this script). The fraction is a statement about that curated list,",
        "> not a survey of the library, and a different choice of proxy pairs would give a",
        "> different fraction. Several corroborated rows also rest on a single extracted",
        "> assertion, so read the per-row counts rather than the headline.", "",
        "| module | pair | relation articles | full-text co-mentions | predicates | cited PMID? | contested |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in sorted(rows, key=lambda r: -r["total"]):
        pair = f"`{r['a']}` - `{r['b']}`"
        if not r["resolved"]:
            L.append(f"| {r['module']} | {pair} | _unresolved_ | - | - | - | - |")
            continue
        preds = ", ".join(f"{k} {v}" for k, v in
                          sorted(r["preds"].items(), key=lambda kv: -kv[1])[:3]) or "-"
        cited = {True: "yes", False: "no", None: "-"}[r["cited_present"]]
        con = (f"**yes** (+{r['pos']}/-{r['neg']}, bal {r['balance']:.2f})"
               if r.get("contested") else "no")
        cm = r.get("comention")
        cm_s = "-" if cm is None else f"**{cm:,}**"
        L.append(f"| {r['module']} | {pair} | {r['total']:,} | {cm_s} | {preds} | "
                 f"{cited} | {con} |")

    contested = [r for r in rows if r.get("contested")]
    if contested:
        L += ["", "## Claims that sit on a CONTESTED edge", "",
              "These pairs are asserted in BOTH directions by the literature",
              "(`positive_correlate` and `negative_correlate`), yet each module cites a",
              "single paper and its docs do not mention the disagreement. See",
              "`analysis/atlas-contradictions.md` for the general catalogue.", "",
              "| module | pair | + | - | balance | cited PMID |",
              "|---|---|---|---|---|---|"]
        for r in sorted(contested, key=lambda r: -(r["balance"] or 0)):
            L.append(f"| `{r['module']}` | {r['a']} - {r['b']} | {r['pos']} | {r['neg']} | "
                     f"{r['balance']:.2f} | {r['pmid']} |")
        L += ["",
              "A high balance means the field is genuinely split; a low one means the module",
              "is on the majority side of a mostly-settled question. Neither says the module",
              "is wrong -- it says the module docs should state which side they took.", ""]

    L += ["", "## Reading", "",
          f"* **{len(found)} of {len(rows)}** module claims are corroborated by at least one",
          "  other cancer article in the graph, so they are not single-paper assertions.",
          ]
    if none_:
        rescued = [r for r in none_ if r.get("comention")]
        L += [f"* **{len(none_)}** resolved to real entities but have NO asserted relation "
              "in the abstract-level graph:",
              ""] + [f"  * `{r['module']}`: {r['a']} - {r['b']} — {r['claim']}"
                     + (f"  _(but **{r['comention']:,}** full-text co-mentions)_"
                        if r.get("comention") else "")
                     for r in none_] + [""]
        if rescued:
            L += [f"  **{len(rescued)} of those {len(none_)} ARE discussed in full text.** A zero in",
                  "  the relation column is therefore not evidence against the mechanism.",
                  "  It has two distinct causes, and they were once conflated here:", "",
                  "  1. **Abstract-level extraction sparsity.** PubTator's relations come from",
                  "     abstracts, so a mechanism established in a Results section may never",
                  "     be asserted in an extractable sentence.",
                  "  2. **Entity-resolution collisions.** The relation WAS extracted and then",
                  "     filed under the wrong gene. `fsp1` (AIFM2-GPX4) is the measured case:",
                  "     it read zero not because nobody asserts GPX4-FSP1, but because every",
                  "     such relation in the census was filed under ATL1, a spastic-paraplegia",
                  "     gene sharing the `FSP1` alias. See `analysis/atlas-disambiguation.md`;",
                  "     the corrections are applied at index build, so this column now counts",
                  "     them. An earlier version of this report attributed that zero entirely",
                  "     to cause 1, which was wrong.", ""]
    if unres:
        L += [f"* **{len(unres)}** could not be resolved to an entity identifier at all:",
              ""] + [f"  * `{r['module']}`: {r['a']} - {r['b']}" for r in unres] + [""]

    L += ["## Detail", ""]
    for r in sorted(rows, key=lambda r: -r["total"]):
        if not r["resolved"] or r["total"] == 0:
            continue
        L += [f"### `{r['module']}` — {r['a_name']} / {r['b_name']}", "",
              f"{r['claim']}", "",
              f"* {r['total']:,} asserting relations",
              "* predicates: " + ", ".join(f"`{k}` {v}" for k, v in
                                           sorted(r["preds"].items(), key=lambda kv: -kv[1])),
              f"* module cites PMID {r['pmid']}"
              + (" — **present in the graph**" if r["cited_present"] else
                 " — not among the sampled supporting PMIDs"),
              f"* example PMIDs: {', '.join(r['pmids'])}", ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"{len(found)}/{len(rows)} claims corroborated; {len(none_)} with no relation; "
          f"{len(unres)} unresolved")
    for r in sorted(found, key=lambda r: -r["total"])[:10]:
        c = {True: "cited-present", False: "cited-absent", None: "-"}[r["cited_present"]]
        print(f"  {r['module']:<16}{r['a']:>10} - {r['b']:<10}{r['total']:>7,}  {c}")


if __name__ == "__main__":
    main()
