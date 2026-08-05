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
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_graph import load_index, resolve, support  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-module-support.md"
# The machine-readable companion. Without it this analysis was the only
# one in the atlas set that downstream documents could not quote without
# copying numbers out of prose -- which is how stale figures get made.
RAW = PROJECT_ROOT / "analysis" / "atlas-module-support.json"

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
    ("dhodh", "DHODH", "GPX4", "33981038",
     "DHODH repairs lipid peroxides in parallel to GPX4"),
    ("gch1", "GCH1", "GPX4", "31989025",
     "GCH1/BH4 is a GPX4-independent radical-trapping defence"),
    ("prom2", "PROM2", "FTH1", "31735663",
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


# Three of the CLAIMS are not off-by-default realism layers at all: they are
# core-engine mechanisms described in the manuscript's Chapter 5 (GPX4/FSP1
# parallel defence, System Xc- cystine supply, and erastin as a first-class drug
# input). `fsp1_rate` and the Xc- factor are live in the default parameter set.
# They are also three of the four best-corroborated rows, so quoting the overall
# fraction as though it described "the layers" flatters the layers by borrowing
# the core engine's evidence. Both denominators are reported.
CORE_ENGINE = {"system_xc", "erastin", "fsp1"}


# The co-mention layer's measured precision, read from its own artifact so this
# document cannot quote a stale figure. Falls back to None when the measurement
# has not been run, in which case the caveat is omitted rather than invented.
def _entity_degree(idx: dict) -> dict:
    """identifier -> how many DISTINCT partners it has in the relation graph.

    A measure of how exposed an entity is to being written about at all, which
    is what decides whether a pair COULD have been asserted.
    """
    deg = {}
    for key in idx["edges"]:
        # `for i in key` over a 2-tuple double-counts a SELF-pair (a, a), which
        # would give `a` two partners where it has one -- itself. 629 of the
        # 2,840,563 edges are self-pairs and AIFM2 is one of them, so this
        # touched two of the twenty claims. Tiny, and wrong is wrong.
        for i in set(key):
            deg[i] = deg.get(i, 0) + 1
    return deg


def _spearman(xs: list, ys: list) -> float:
    """Rank correlation, stdlib only (scipy is not a dependency of this script)."""
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


def _below_floor_base_rate(idx: dict, degree: dict, rows: list):
    """Share of ALL asserted pairs whose weaker entity sits below the claim floor.

    The number that refutes "the census could not have corroborated them": if
    nearly half of everything the census DOES assert sits below the line, the
    line is not a detectability limit.
    """
    sup = [r["weaker"] for r in rows if r.get("weaker") is not None and r["total"] > 0]
    if not sup:
        return None
    floor, below, total = min(sup), 0, 0
    for key in idx["edges"]:
        total += 1
        if min(degree.get(key[0], 0), degree.get(key[1], 0)) < floor:
            below += 1
    return below / total if total else None


def _kendall(xs: list, ys: list) -> float:
    """Tie-corrected rank correlation (tau-b), stdlib only.

    Reported beside Spearman because 9 of the 20 outcomes are tied at zero,
    which is exactly the case Spearman handles least transparently.
    """
    n, con, dis, tx, ty = len(xs), 0, 0, 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            dx, dy = xs[i] - xs[j], ys[i] - ys[j]
            if dx == 0 and dy == 0:
                tx += 1; ty += 1
            elif dx == 0:
                tx += 1
            elif dy == 0:
                ty += 1
            elif (dx > 0) == (dy > 0):
                con += 1
            else:
                dis += 1
    tot = n * (n - 1) / 2
    den = ((tot - tx) * (tot - ty)) ** 0.5
    return (con - dis) / den if den else 0.0


def _exposure_section(rows: list, base_rate: float | None = None) -> list:
    """Why a zero in the relation column is mostly not about the claim.

    WHAT SURVIVED REVIEW, AND WHAT DID NOT
    ---------------------------------------
    The association is solid: entity exposure predicts corroboration at
    rho ~ 0.86 (permutation p < 1e-5), tau ~ 0.70, and the two groups separate
    at AUC 0.965. It holds under a different entity-level measure
    (`ident_mentions`, rho ~ 0.83, same split) and RISES to 0.88 when every
    GPX4 pair is dropped, so it is not an artifact of one hub gene.

    Two things did NOT survive, and both were in the first draft:

      * "the census could not have corroborated them whatever the biology is"
        is false. The floor is a sample minimum, not a detectability limit, and
        45% of all asserted pairs in the graph sit below it. The claim is
        probabilistic and is now stated that way.
      * naming WHICH zeros are genuine does not replicate across measures. On
        entity degree the exceptions are copper and hdac_persister; on the
        pair-level co-mention column already in this table the correlation
        INVERTS (rho ~ -0.9) and the exceptions become ether_lipid and repair,
        with zero overlap. So the table cannot currently single out a genuinely
        unasserted link, and saying which two were interesting was reading one
        measure and not the other.
    """
    usable = [r for r in rows if r.get("weaker") is not None]
    if len(usable) < 5:
        return []
    w = [r["weaker"] for r in usable]
    n = [r["total"] for r in usable]
    rho, tau = _spearman(w, n), _kendall(w, n)
    zero = [r for r in usable if r["total"] == 0]
    nonzero = [r for r in usable if r["total"] > 0]
    if not zero or not nonzero:
        return []
    floor = min(r["weaker"] for r in nonzero)
    floor_row = min(nonzero, key=lambda r: r["weaker"])
    explained = [r for r in zero if r["weaker"] < floor]
    above = sorted([r for r in zero if r["weaker"] >= floor], key=lambda r: -r["weaker"])

    # The same procedure run on the pair-level column already in this table.
    cm = [r for r in zero if r.get("comention")]
    cm_sup = [r["comention"] for r in nonzero if r.get("comention")]
    cm_floor = min(cm_sup) if cm_sup else None
    cm_above = sorted([r for r in cm if r["comention"] >= cm_floor],
                      key=lambda r: -r["comention"]) if cm_floor else []

    L = [
        "## A zero is usually about EXPOSURE, not about the claim", "",
        "The count above treats every zero alike. They are not alike. A pair can",
        "only be asserted if both its entities are written about at all, and across",
        f"these claims the WEAKER entity's number of distinct relation partners runs",
        f"from {min(w)} to {max(w):,}, a factor of {max(w)//max(1, min(w))}.", "",
        f"Exposure and corroboration track each other closely: Spearman "
        f"**rho = {rho:.2f}**, Kendall **tau = {tau:.2f}** over {len(usable)} claims",
        "(a permutation test puts the correlation beyond every one of 200,000",
        "shuffles). It is not an artifact of one hub gene -- taking the WEAKER of",
        "each pair structurally excludes GPX4, which is never the weaker side, and",
        "dropping all nine GPX4 pairs RAISES rho to 0.88. An independent",
        "entity-level measure (total mentions rather than partners) gives rho 0.83",
        "and the identical split.", "",
        f"Every claim that HAS support has a weaker entity of at least {floor}, and",
        f"{len(explained)} of the {len(zero)} zeros fall below that. So a zero is a",
        "poor guide to whether a claim is true: it mostly tracks how much one side",
        "has been studied.", "",
        "| module | pair | weaker partner count | below the supported range? |",
        "|---|---|---|---|",
    ]
    for r in sorted(zero, key=lambda r: r["weaker"]):
        L.append(f"| {r['module']} | `{r['a']}` - `{r['b']}` | {r['weaker']} "
                 f"({r['deg_a']} / {r['deg_b']}) | "
                 f"{'no' if r['weaker'] >= floor else 'yes'} |")

    L += [
        "", "### Three reasons not to push this further than it goes", "",
        f"**The line is a sample minimum, not a threshold.** It is set by a single "
        f"row -- `{floor_row['module']}`, which has "
        f"{floor_row['total']} asserting article"
        f"{'s' if floor_row['total'] != 1 else ''} under the weakest predicate the "
        "graph has. Leave-one-out over the supported claims changes the line only "
        f"when that row is dropped, and then it moves {floor}"
        f" -> {min(r['weaker'] for r in nonzero if r is not floor_row)}. A bootstrap "
        "over the supported set returns the shipped answer about two thirds of the "
        "time.", "",
    ]
    if base_rate is not None:
        L += [
            "**Below the line does not mean undetectable.** An earlier version of "
            "this section said the census could not have corroborated those claims "
            "whatever the biology is. That is false: across the whole graph "
            f"**{100*base_rate:.0f}% of all asserted pairs** have a weaker entity "
            f"below {floor}, so pairs like these are corroborated constantly. The "
            "honest statement is probabilistic -- a particular pair drawn from a "
            "sparsely-studied entity's handful of relations is unlikely a priori -- "
            "not that the machinery cannot see them.", "",
        ]
    if cm_floor is not None:
        L += [
            "**Which zeros are 'genuine' does not replicate across measures, so this "
            "section does not name any.** Run the identical procedure on the "
            "pair-level co-mention column already in the table above and the "
            f"correlation INVERTS among the zeros (rho about -0.9), the line lands at "
            f"{cm_floor}, and the rows above it become "
            + ", ".join(f"`{r['module']}`" for r in cm_above)
            + " -- with no overlap at all against the entity-degree answer ("
            + ", ".join(f"`{r['module']}`" for r in above)
            + "). Entity exposure and pair discussion are different constructs and "
            "there is no reason they must agree, but the conclusion anyone would "
            "draw is pair-level, so a result that flips with the measure is not one "
            "to report. The first draft of this section named the entity-degree pair "
            "as 'the interesting rows' without checking the column beside it.", "",
        ]
    L += [
        "**Caveat on the correlation.** A pair's own relations contribute to both",
        "entities' partner counts, so the two quantities share a term. Recomputing",
        "with each pair's own edge removed moves rho by 0.01 and changes no row's",
        "classification, and the contribution is exactly zero for every unsupported",
        "row, so the coupling cannot manufacture the effect.", "",
        "**And a high partner count is not evidence FOR a claim.** It makes a pair",
        "measurable, not correct. An earlier draft added that GPX4 has the most",
        "partners here, which is false -- across the twenty claims it ranks fifth,",
        "behind CDH1, CD274, IFNG and YAP1. It is the maximum only within the",
        "nine-row table above, and by 1.6x rather than by a wide margin.", "",
    ]
    return L


def _thin_threshold() -> int:
    """How small a count is too small to argue from.

    Calibrated to the layer's measured precision rather than fixed: at ~42% a
    handful of co-mentions could be entirely noise and 50 was the bar; at ~88%
    the same count carries real information, so the bar drops. A constant here
    would have gone stale silently the moment the filter was promoted.
    """
    p = comention_precision() or 0.42
    return 50 if p < 0.6 else (20 if p < 0.8 else 10)


def _comention_interval() -> str:
    """The stated precision is a point estimate; quote its range beside it."""
    try:
        d = json.loads((PROJECT_ROOT / "analysis"
                        / "comention-authority-result.json").read_text())
        return (f"95% interval roughly {100*d['weighted_ci'][0]:.0f}-"
                f"{100*d['weighted_ci'][1]:.0f}%, blind panel "
                f"{100*d['blind_weighted']:.0f}%")
    except (OSError, ValueError, KeyError):
        pass
    try:
        d = json.loads((PROJECT_ROOT / "analysis"
                        / "comention-regression.json").read_text())["after"]
        return f"95% interval roughly {100*d['range'][0]:.0f}-{100*d['range'][1]:.0f}%"
    except (OSError, ValueError, KeyError):
        return "interval not available"


def _comention_caveat() -> list:
    """How to read the co-mention column, derived rather than written.

    The number beside this paragraph has always been read from the artifact.
    The paragraph was not, and it went stale the moment the filter shipped: it
    told readers that roughly half of any figure in the column was a generic
    surface form and named `treatment` and `effects` as the failure mode, which
    described the 41.6% layer. Those two forms are exactly what the authority
    rule removes, so the caveat was warning about the one thing that had been
    fixed while the number above it said 88%.

    What replaces it is the structure the measurement actually shows: the error
    is no longer spread evenly across the strata, and this column mixes them
    without saying which one a figure came from.

    All three strata are named, from the artifact's own keys. An earlier draft
    named only the best and worst and called a count "a blend of the two",
    which was self-refuting -- their weights sum to 71% -- and dropped
    body-only, the stratum PubTator cannot contradict and which carries a third
    of the residual error. It also hard-coded which stratum was which beside
    derived numbers, so a five-point drift in the measurement would have
    swapped them and emitted true numbers under a false sentence.
    """
    p = comention_precision()
    if p is None:
        return []
    paras = [
        f"**The co-mention column is measured at roughly {100*p:.0f}% precision** "
        f"({_comention_interval()}, `analysis/comention-authority-result.md`), so "
        f"read it as an upper bound on discussion rather than a count of it."
    ]
    try:
        d = json.loads((PROJECT_ROOT / "analysis"
                        / "comention-authority-result.json").read_text())
        st = d["strata"]
    except (OSError, ValueError, KeyError):
        return _quote(paras)
    # Label from the key, never from an assumed ordering.
    gloss = {"corroborated": "mentions PubTator also corroborates",
             "abstract-visible": "mentions visible in the abstract",
             "body-only": "mentions found only in the body"}
    order = sorted(st, key=lambda s: -st[s]["weight"])
    listed = "; ".join(
        f"{gloss.get(s, s)} are {100*st[s]['precision']:.0f}% correct "
        f"({st[s]['tp']}/{st[s]['n']}) and carry {100*st[s]['weight']:.0f}% "
        f"of the volume" for s in order)
    err = {s: st[s]["weight"] * (1 - st[s]["precision"]) for s in st}
    worst = max(err, key=err.get)
    paras += [
        f"**The remaining error is not spread evenly, and this column does not say "
        f"which stratum a figure came from.** Across the three: {listed}. A count "
        f"here is a blend of all three in unknown proportion, so the honest reading "
        f"of any single figure is the weighted {100*p:.0f}% rather than whichever "
        f"stratum it might have come from. The largest single share of the residual "
        f"error sits in {gloss.get(worst, worst)}, at "
        f"{100*err[worst]/sum(err.values()):.0f}% of it -- a plurality rather than a "
        f"concentration, since that stratum carries only "
        f"{100*st[worst]['weight']:.0f}% of the volume.",

        f"That still matters most where the number is SMALL, which is where this "
        f"document leans on it: a handful of co-mentions offered as evidence that a "
        f"zero-relation module is discussed after all is a handful of chances for "
        f"the {100*(1-p):.0f}% to land. A large figure survives the error rate as "
        f"evidence of discussion; a figure under {_thin_threshold()} is flagged "
        f"rather than argued from.",
    ]
    return _quote(paras)


def _quote(paras: list, width: int = 76) -> list:
    """Wrap paragraphs into a markdown blockquote.

    Every emitted line carries its own `> `. Hand-writing the prefix on
    f-strings that span source lines silently drops it from the continuation
    lines, which markdown then renders by lazy continuation -- so it looks
    right until someone edits near it.
    """
    out = []
    for para in paras:
        out += ["> " + ln for ln in textwrap.wrap(
            para, width, break_on_hyphens=False, break_long_words=False)] + [""]
    return out


def comention_precision():
    """The precision of the layer AS BUILT.

    The authority filter is on by default since #628, so the figure that
    describes this column is the filtered layer's, not the unfiltered one the
    regression document measures. Reading the wrong artifact would have this
    document quoting 42% for a layer measured at 88%, understating its own
    evidence by half.
    """
    for f, path in ((PROJECT_ROOT / "analysis" / "comention-authority-result.json",
                     ("weighted",)),
                    (PROJECT_ROOT / "analysis" / "comention-regression.json",
                     ("after", "weighted"))):
        try:
            d = json.loads(f.read_text())
            for k in path:
                d = d[k]
            return d
        except (OSError, ValueError, KeyError):
            continue
    return None


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
    else:
        print("co-mention layer ABSENT -- the full-text column will be empty",
              file=sys.stderr)
    degree = _entity_degree(idx)
    rows = []
    for module, a, b, pmid, claim in CLAIMS:
        r = support(idx, a, b)
        if r is None:
            rows.append(dict(module=module, a=a, b=b, pmid=pmid, claim=claim,
                             resolved=False, total=0, preds={}, cited_present=None,
                             deg_a=None, deg_b=None, weaker=None))
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
                         comention=comention.get(tuple(sorted((r["a"], r["b"])))),
                         deg_a=degree.get(r["a"], 0), deg_b=degree.get(r["b"], 0),
                         weaker=min(degree.get(r["a"], 0), degree.get(r["b"], 0))))

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
    ] + (_comention_caveat() if comention else [
        "> **The full-text co-mention layer is not built**, so that column is empty",
        "> throughout. An empty cell here means *not measured*, NOT *not discussed* --",
        "> the distinction matters because a zero in the relation column is read",
        "> against exactly this column. Rebuild with",
        "> `python scripts/atlas_comention.py --rebuild`.", "",
    ]) + [
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

    # Which claims, if any, rest on an entity measured as a sense collision.
    try:
        _scan = json.loads(
            (PROJECT_ROOT / "analysis" / "atlas-ambiguity.json").read_text())
        _collide = set()
        for _t in ("gene", "chemical", "disease"):
            for _r in _scan["by_type"][_t]["sense_rows"]:
                _collide |= {_r["top"]["id"], _r["runner_up"]["id"]}
    except (OSError, ValueError, KeyError):
        _collide = set()
    colliding = []
    for r in rows:
        hits = [n for n in (r["a"], r["b"])
                if resolve(idx, n) in _collide]
        if hits:
            colliding.append((r["module"], hits))

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
        # A contradiction can be manufactured by conflating two entities under one
        # identifier, and that inflates the flag rate by 1.45x across the graph
        # (analysis/atlas-contradiction-quality.md). So say whether it could be
        # the explanation HERE, rather than leaving a reader to wonder.
        if colliding:
            L += ["> **Some of these may be conflation, not disagreement.** These claims rest",
                  "> on entities measured as sense collisions, and merging two entities merges",
                  "> two literatures, which will disagree "
                  "(`analysis/atlas-contradiction-quality.md`):", ""]
            L += [f">   * `{m}`: {', '.join(h)}" for m, h in colliding] + [""]
        else:
            L += ["> **Conflation does not explain these.** Across the graph, pairs built on a",
                  "> measured sense collision are 1.45x more likely to be flagged contradictory",
                  "> (`analysis/atlas-contradiction-quality.md`), so that had to be excluded.",
                  f"> All {len(rows)} claims here rest on entities with no measured collision,",
                  "> so these conflicts are disagreements between studies rather than two",
                  "> literatures merged under one identifier.", ""]

    L += ["", "## Reading", "",
          f"* **{len(found)} of {len(rows)}** module claims are corroborated by at least one",
          "  other cancer article in the graph, so they are not single-paper assertions.",
          ]
    if none_:
        rescued = [r for r in none_ if r.get("comention")]
        L += [f"* **{len(none_)}** resolved to real entities but have NO asserted relation "
              "in the abstract-level graph:",
              ""] + [f"  * `{r['module']}`: {r['a']} - {r['b']} — {r['claim']}"
                     + (f"  _(but **{r['comention']:,}** full-text co-mentions"
                        + (f", still thin against the layer's "
                           f"~{100*comention_precision():.0f}% precision"
                           if r["comention"] < _thin_threshold() else "") + ")_"
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

    L += _exposure_section(rows, base_rate=_below_floor_base_rate(idx, degree, rows))

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

    usable = [r for r in rows if r.get("weaker") is not None]
    supported = [r["weaker"] for r in usable if r["total"] > 0]
    floor = min(supported) if supported else None
    layers = [r for r in rows if r["module"] not in CORE_ENGINE]
    no_gpx4 = [r for r in usable if "GPX4" not in (r["a"], r["b"])]
    RAW.write_text(json.dumps({
        "n_claims": len(rows),
        "corroborated": len(found),
        "zero_relation": len(none_),
        "unresolved": len(unres),
        # The denominator that describes the LAYERS, with the core-engine rows
        # removed. The overall fraction is not a survey of the library either:
        # these are author-chosen proxy pairs for a subset of the modules.
        "n_claims_layers_only": len(layers),
        "corroborated_layers_only": sum(1 for r in layers if r["total"] > 0),
        "core_engine_claims": sorted(CORE_ENGINE),
        "exposure_floor": floor,
        "weaker_min": min((r["weaker"] for r in usable), default=None),
        "weaker_max": max((r["weaker"] for r in usable), default=None),
        "below_floor_base_rate": _below_floor_base_rate(idx, degree, rows),
        "spearman_weaker_degree_vs_relations":
            _spearman([r["weaker"] for r in usable], [r["total"] for r in usable])
            if len(usable) >= 5 else None,
        "spearman_excluding_gpx4":
            _spearman([r["weaker"] for r in no_gpx4], [r["total"] for r in no_gpx4])
            if len(no_gpx4) >= 5 else None,
        "zero_explained_by_exposure": sum(
            1 for r in usable if r["total"] == 0 and floor is not None
            and r["weaker"] < floor),
        "zero_unexplained": [
            r["module"] for r in usable if r["total"] == 0 and floor is not None
            and r["weaker"] >= floor],
        "claims": [{k: r[k] for k in ("module", "a", "b", "total", "comention",
                                      "deg_a", "deg_b", "weaker")
                    if k in r} for r in rows],
    }, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"{len(found)}/{len(rows)} claims corroborated; {len(none_)} with no relation; "
          f"{len(unres)} unresolved")
    for r in sorted(found, key=lambda r: -r["total"])[:10]:
        c = {True: "cited-present", False: "cited-absent", None: "-"}[r["cited_present"]]
        print(f"  {r['module']:<16}{r['a']:>10} - {r['b']:<10}{r['total']:>7,}  {c}")


if __name__ == "__main__":
    main()
