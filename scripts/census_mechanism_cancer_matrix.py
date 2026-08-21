#!/usr/bin/env python3
"""The mechanism-by-cancer-site matrix, and what a "gap" means at census scale.

The manuscript reported a 19x22 mechanism-cancer matrix over a retrieved corpus
and counted 94 empty cells (22.5%) as gaps -- then measured, in its own
sensitivity analysis, that the count falls to 29-38 under a coarser taxonomy.
A quantity that moves by a factor of three with a labelling choice is a
property of the labelling.

AT CENSUS SCALE THE ZERO COUNT COLLAPSES: 6 of 288 cells. That is the first
result and it settles the original question -- the empty cells were
overwhelmingly a property of retrieval, not of the field. It also RETIRES ZERO
AS A MEASURE. With 4.4 million articles almost every mechanism-site pair has
been written about at least once, so counting zeros stops discriminating.

WHAT REPLACES IT is the comparison a count cannot make: how a cell compares to
what its own marginals predict. A mechanism studied in 30,000 articles and a
site carrying 300,000 records will meet often by arithmetic alone; a small
mechanism and a small site will meet rarely without anything being neglected.
The ratio of observed to expected separates those.

THE NULL IS KNOWN FALSE AND THAT IS THE POINT, so no p-value is computed here.
Mechanisms and anatomical sites are genuinely dependent -- CAR-T concentrates in
leukaemia because that is what it treats -- so "significantly different from
independence" is true of nearly every cell and would be a test of a hypothesis
nobody holds. The ratio is reported as an effect size, and the ranking is what
carries.

AND A DEPLETED CELL IS NOT A GAP. HIFU is depleted in leukaemia because
leukaemia has no focal target for an ablative beam; that is the modality
working as designed. This analysis therefore reports depletion and refuses to
name any cell an opportunity -- separating "nobody has tried this" from "this
cannot work here" needs knowledge of the mechanism, not of the counts.
"""
import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
SITE_MAP = REPO / "analysis/site-descriptor-map.tsv"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-mechanism-cancer-matrix.md"
OUT_JSON = REPO / "analysis/census-mechanism-cancer-matrix.json"

# The manuscript's corpus-scale figures, for the comparison this exists to make.
CORPUS_CELLS = 418
CORPUS_ZEROS = 94
# Below this many observed articles a ratio is a statement about a handful.
MIN_EXPECTED = 20.0
# A row's spread is only a property of the mechanism if the row has enough
# cells to have a shape. Microbiome clears the expectation floor in 2 of 18
# sites, and a standard deviation over two numbers is not a spread.
MIN_CELLS_FOR_SPREAD = 8


def load_sites() -> dict:
    out = defaultdict(set)
    for ln in SITE_MAP.read_text(encoding="utf-8").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split("\t")
        if len(p) >= 3:
            out[p[0]].add(p[2].strip().lower())
    return dict(out)


def scan(stride: int = 1) -> dict:
    import yaml

    sites = load_sites()
    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    mech = {k.lower(): {x.lower() for x in v["descriptors"]}
            for k, v in mp.items() if v["descriptors"]}

    cell = defaultdict(Counter)
    mech_n = Counter()
    site_n = Counter()
    universe = 0
    n = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if not ms:
                    continue
                hit_m = [k for k, d in mech.items() if ms & d]
                if not hit_m:
                    continue
                hit_s = [s for s, d in sites.items() if ms & d]
                if not hit_s:
                    continue
                # THE UNIVERSE IS ARTICLES CARRYING BOTH, and it has to be, or
                # the expectation is computed against a population most of the
                # matrix could never have entered. Marginals are counted over
                # this same universe for the same reason.
                universe += 1
                for k in hit_m:
                    mech_n[k] += 1
                for s in hit_s:
                    site_n[s] += 1
                for k in hit_m:
                    for s in hit_s:
                        cell[k][s] += 1
    return {
        "census": n,
        "universe": universe,
        "mechanism_totals": dict(mech_n),
        "site_totals": dict(site_n),
        "cells": {k: dict(v) for k, v in cell.items()},
        "min_expected": MIN_EXPECTED,
    }


def assemble(d: dict) -> dict:
    mt, st = d["mechanism_totals"], d["site_totals"]
    U = d["universe"]
    rows = []
    zeros = []
    for m in sorted(mt):
        for s in sorted(st):
            obs = d["cells"].get(m, {}).get(s, 0)
            exp = mt[m] * st[s] / U if U else 0.0
            row = {"mechanism": m, "site": s, "observed": obs,
                   "expected": round(exp, 1),
                   "ratio": round(obs / exp, 2) if exp else None,
                   # A ratio off a small expectation is a statement about a
                   # handful of articles, not about a field.
                   "interpretable": exp >= d["min_expected"]}
            rows.append(row)
            if obs == 0:
                zeros.append(row)
    interp = [r for r in rows if r["interpretable"] and r["ratio"] is not None]

    # PER-MECHANISM SPREAD. Arguably the better summary than either tail: it
    # asks how site-specific a mechanism's literature is, without needing a
    # cell to be extreme.
    by_m = {}
    for r in interp:
        by_m.setdefault(r["mechanism"], []).append(r["ratio"])
    spread = []
    for m, rat in sorted(by_m.items()):
        if len(rat) < MIN_CELLS_FOR_SPREAD:
            continue
        mean = sum(rat) / len(rat)
        sd = (sum((x - mean) ** 2 for x in rat) / len(rat)) ** 0.5
        spread.append({"mechanism": m, "cells": len(rat),
                       "articles": mt[m], "spread": round(sd, 2),
                       "min_ratio": round(min(rat), 2),
                       "max_ratio": round(max(rat), 2)})
    spread.sort(key=lambda r: -r["spread"])
    out = dict(d)
    out["rows"] = rows
    out["n_cells"] = len(rows)
    out["n_zero"] = len(zeros)
    out["zero_cells"] = [f"{r['mechanism']}/{r['site']}" for r in zeros]
    out["n_interpretable"] = len(interp)
    out["most_depleted"] = sorted(interp, key=lambda r: r["ratio"])[:12]
    out["most_enriched"] = sorted(interp, key=lambda r: -r["ratio"])[:12]
    out["spread"] = spread
    out["min_cells_for_spread"] = MIN_CELLS_FOR_SPREAD
    # THE CONTROL: are the tails just small mechanisms with noisy ratios? If
    # spread tracked size the ranking would be an artifact of the denominator
    # rather than a property of the field, so it is measured rather than
    # argued.
    import math

    if len(spread) >= 4:
        xs = [math.log10(r["articles"]) for r in spread]
        ys = [r["spread"] for r in spread]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
        den = (sum((a - mx) ** 2 for a in xs)
               * sum((b - my) ** 2 for b in ys)) ** 0.5
        out["size_vs_spread_r"] = round(num / den, 2) if den else None
    else:
        out["size_vs_spread_r"] = None
    # A mechanism in BOTH tails is showing conservation, not two findings:
    # concentrating in a few sites forces depletion elsewhere.
    dep = {r["mechanism"] for r in out["most_depleted"]}
    enr = {r["mechanism"] for r in out["most_enriched"]}
    out["in_both_tails"] = sorted(dep & enr)
    # The comparison this analysis exists to make.
    out["corpus_zero_share"] = round(100 * CORPUS_ZEROS / CORPUS_CELLS, 1)
    out["census_zero_share"] = round(100 * len(zeros) / len(rows), 1)
    return out


def render(d: dict) -> str:
    L = ["# The mechanism-by-site matrix at census scale\n"]
    L.append(
        f"Generated by `scripts/census_mechanism_cancer_matrix.py`. The universe "
        f"is the {d['universe']:,} census articles carrying BOTH a mechanism "
        f"descriptor and an anatomical-site descriptor; both marginals are "
        f"counted over that same universe, because an expectation computed "
        f"against a population most of the matrix could never enter is not an "
        f"expectation.\n"
    )
    L.append("## Zero cells: the original question, answered and retired\n")
    L.append(
        f"The manuscript counted **{CORPUS_ZEROS} of {CORPUS_CELLS} cells "
        f"({d['corpus_zero_share']}%)** empty over a retrieved corpus, and its "
        f"own sensitivity analysis found the count falls to 29-38 under a "
        f"coarser taxonomy. At census scale it is **{d['n_zero']} of "
        f"{d['n_cells']} ({d['census_zero_share']}%)**"
        + (": " + ", ".join(f"`{z}`" for z in d["zero_cells"]) if d["zero_cells"]
           else ".")
        + " So the empty cells were overwhelmingly a property of the retrieval, "
          "which is what the manuscript suspected and could not show.\n"
    )
    L.append(
        "That also RETIRES the measure. With 4.4 million articles nearly every "
        "mechanism-site pair has been written about at least once, so counting "
        "zeros no longer discriminates between a neglected combination and a "
        "well-studied one. A gap measure that returns 2% on the full "
        "literature is not detecting gaps.\n"
    )
    L.append("## What replaces it\n")
    L.append(
        f"How a cell compares to what its own marginals predict. A mechanism in "
        f"tens of thousands of articles and a site carrying hundreds of "
        f"thousands of records meet often by arithmetic alone. Ratios are shown "
        f"only where the expectation reaches {d['min_expected']:.0f} articles "
        f"({d['n_interpretable']} of {d['n_cells']} cells); below that a ratio "
        f"describes a handful.\n"
    )
    L.append("### Most depleted\n")
    L.append("| mechanism | site | observed | expected | ratio |")
    L.append("|---|---|--:|--:|--:|")
    for r in d["most_depleted"]:
        L.append(f"| {r['mechanism']} | {r['site']} | {r['observed']:,} | "
                 f"{r['expected']:,.0f} | {r['ratio']} |")
    L.append("")
    L.append("### Most concentrated\n")
    L.append("| mechanism | site | observed | expected | ratio |")
    L.append("|---|---|--:|--:|--:|")
    for r in d["most_enriched"]:
        L.append(f"| {r['mechanism']} | {r['site']} | {r['observed']:,} | "
                 f"{r['expected']:,.0f} | {r['ratio']} |")
    L.append("")
    L.append("### How site-specific each mechanism is\n")
    L.append(
        f"Spread is the standard deviation of a mechanism's ratios across the "
        f"sites where its expectation is interpretable, shown for the "
        f"{len(d['spread'])} mechanisms clearing {d['min_cells_for_spread']} "
        f"such cells. It summarises a whole row without needing any single cell "
        f"to be extreme, and it separates a therapy aimed at particular organs "
        f"from a platform or a process that appears wherever cancer is "
        f"studied.\n"
    )
    L.append("| mechanism | articles | cells | spread | lowest | highest |")
    L.append("|---|--:|--:|--:|--:|--:|")
    for r in d["spread"]:
        L.append(f"| {r['mechanism']} | {r['articles']:,} | {r['cells']} | "
                 f"{r['spread']} | {r['min_ratio']} | {r['max_ratio']} |")
    L.append("")
    if d["size_vs_spread_r"] is not None:
        rr = d["size_vs_spread_r"]
        L.append(
            f"**The obvious objection, tested.** A small mechanism spread over "
            f"eighteen sites has small expectations per cell and noisier "
            f"ratios, so the tails could be a property of the denominator "
            f"rather than of the field. Correlation between a mechanism's size "
            f"and its spread is **{rr:+.2f}** over {len(d['spread'])} "
            f"mechanisms"
            + (", so size does not predict spread and the ranking is not a "
               "small-number artifact.\n" if abs(rr) < 0.4 else
               ", which is large enough that the ranking is partly a "
               "measurement of mechanism size and should not be read as a "
               "property of the field alone.\n")
        )
    if d["in_both_tails"]:
        L.append(
            f"**{', '.join(f'`{m}`' for m in d['in_both_tails'])} appear in "
            f"BOTH tails, and that is one fact rather than two.** A row's "
            f"ratios are constrained: concentrating a fixed literature into a "
            f"few sites forces depletion in the rest. So for these mechanisms "
            f"the depleted cells are the arithmetic shadow of the concentrated "
            f"ones, and reading the two tables as independent evidence would "
            f"double-count a single pattern.\n"
        )
    L.append("## No p-value, deliberately\n")
    L.append(
        "The null here is that a mechanism's literature is spread across sites "
        "independently of which site it is, and that null is KNOWN FALSE before "
        "any data is read -- CAR-T concentrates in leukaemia because that is "
        "what it treats. A test against it would return significance for nearly "
        "every cell and would be testing a hypothesis nobody holds. The ratio "
        "is reported as an effect size and the ranking is what carries.\n"
    )
    L.append("## A depleted cell is not a gap\n")
    L.append(
        "This analysis names no opportunities, and the restraint is not "
        "modesty. An ablative modality is depleted in disseminated disease "
        "because there is nothing to ablate -- the modality working as designed, "
        "not an oversight. Separating *nobody has tried this* from *this cannot "
        "work here* needs knowledge of the mechanism, and the counts do not "
        "carry it. What the table supplies is a ranked list of places where the "
        "literature is thinner than arithmetic predicts, which is where such a "
        "judgement would start.\n"
    )
    L.append("## Two limits inherited from the labels\n")
    L.append(
        f"The matrix has {len(d['mechanism_totals'])} mechanisms, not the "
        f"manuscript's 19, because a mechanism MeSH cannot express cannot "
        f"appear in it at all. TTFields and bioelectric modulation are ABSENT "
        f"rather than empty, and their absence is invisible in a zero count -- "
        f"which is worth stating precisely because this page reports the zero "
        f"count falling.\n"
    )
    L.append(
        "Sites come from NLM's C04 tree and mechanisms from descriptor sets of "
        "varying breadth, so a cell inherits both. A broad mechanism descriptor "
        "spreads across sites more evenly than the therapy it names, which "
        "flattens its row toward a ratio of 1 and makes it look less "
        "concentrated than it is.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only
                 else scan(a.stride))
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  universe {d['universe']:,}  zeros {d['n_zero']}/{d['n_cells']}  "
          f"interpretable {d['n_interpretable']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
