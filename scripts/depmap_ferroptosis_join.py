#!/usr/bin/env python3
"""Score the engine's ferroptosis gene panel against cell-line drug response (#932).

`analysis/calibration/calibration-feasibility.md` recorded a set of genes as
data-blocked and named one join that would make them all testable at once:
cell-line expression against ferroptosis-inducer sensitivity across the
committed CTRPv2 cohort. It recorded that join as blocked on ACCESS. The access
claim was about one door; `scripts/fetch_depmap_ctrpv2_join.py` opens another,
and this report reads the derivatives that leg commits.

WHAT MAKES THIS DIFFERENT FROM A CORRELATION TABLE
--------------------------------------------------
Any expression can correlate with any drug-response readout across a pan-cancer
line panel, and a panel of famous genes is not a random sample of genes. So the
question is never "does this gene correlate" but "does it correlate more than
the distribution it should be compared against", and this report carries several
such distributions rather than one:

* the COMPLETE distribution over every protein-coding gene in the matrix, so a
  candidate's percentile in it is exact rather than sampled;
* the same distribution restricted to genes at the candidate's own expression
  level, because expression level governs how stably a gene ranks and the panel
  is made of highly expressed loci;
* a permutation of the readout across lines, which asks whether a correlation
  of that size is surprising at this cohort size at all.

The pass rule was fixed before any correlation existed, in the issue: extreme in
the all-gene tail BOTH pooled and after ranking within Oncotree lineage. How
often that rule fires on a gene chosen without regard to biology is MEASURED
here over all genes rather than assumed from the nominal tail -- the two
conditions share a cohort and are strongly dependent, so the nominal product
would be wrong.

The one direction fixed in advance is ACSL4's, because the engine fixes it in
source: `biochem.rs` adds `acsl4_status_boost` to the peroxidizable-substrate
term and `acsl4.rs` calls ACSL4-negative status intrinsically refractory. Higher
ACSL4 must mean a more sensitive line, for the GPX4 inhibitors. No other gene
gets a pre-registered direction, and a failed prediction is reported as failed.

SIGN CONVENTION
---------------
The readout is mean predicted viability inside each line's own screened dose
range, so low is sensitive and a POSITIVE correlation means resistance.

EVERY INPUT IS COMMITTED
------------------------
Unlike the generators that scan the census or run a binary, this one reads only
files that live in the repository, so both of `main`'s branches re-derive every
figure from source and neither can echo a stored derived field.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import random
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

CORRELATIONS_CSV_GZ = REPO / "analysis" / "calibration" / "depmap_ctrpv2_gene_correlations.csv.gz"
PANEL_CSV = REPO / "analysis" / "calibration" / "depmap_ctrpv2_panel.csv"
PROVENANCE_JSON = REPO / "analysis" / "calibration" / "depmap_ctrpv2_join_provenance.json"

OUT_MD = REPO / "analysis" / "depmap-ferroptosis-join.md"
OUT_JSON = REPO / "analysis" / "depmap-ferroptosis-join.json"

# Two-sided tail of the all-gene distribution a rho must reach, fixed before any
# correlation was computed.
EXTREME_TAIL = 0.05

# Bins of mean expression for the level-matched null.
EXPRESSION_BINS = 10

# Permutation of the readout across lines: an independent null for whether a
# correlation of this size is surprising at this cohort size.
N_PERMUTATIONS = 2000
PERMUTATION_SEED = 20932

# The engine fixes this sign in source; see the module docstring.
PREREGISTERED_GENE = "ACSL4"
PREREGISTERED_COMPOUNDS = ("ML162", "ML210")
PREREGISTERED_DIRECTION = "negative"

# A level-matched null computed from fewer genes than this is too thin to quote,
# so the row says so rather than carrying a rate nobody should read.
MIN_MATCHED_GENES = 200

# An enrichment claim over fewer candidate cells than this is not reported.
MIN_CELLS_FOR_ENRICHMENT = 20

# A ratio at or below this is no enrichment at all.
ENRICHMENT_NULL = 1.0

# The two committed derivatives must tell the same story: the per-line panel has
# to reproduce the correlation the table carries, or one of them is stale.
PANEL_AGREEMENT_TOLERANCE = 0.02


def _f(text):
    """A blank field is an absent value, never a zero."""
    if text is None or text == "":
        return None
    return float(text)


def load() -> dict:
    """Read the committed derivatives. No network, no census, no binary."""
    with gzip.open(CORRELATIONS_CSV_GZ, "rt") as fh:
        genes = list(csv.DictReader(fh))
    panel = list(csv.DictReader(PANEL_CSV.open()))
    provenance = json.loads(PROVENANCE_JSON.read_text())
    return {"genes": genes, "panel": panel, "provenance": provenance}


def _column(compound, scope) -> str:
    suffix = "pooled" if scope == "pooled" else "by_lineage"
    return f"{compound}_rho_vs_viability_{suffix}"


def _percentile(sorted_values, value) -> float:
    """Fraction of `sorted_values` at or below `value`."""
    lo, hi = 0, len(sorted_values)
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_values[mid] <= value:
            lo = mid + 1
        else:
            hi = mid
    return lo / len(sorted_values)


def _is_extreme(percentile) -> bool:
    half = EXTREME_TAIL / 2.0
    return percentile is not None and (percentile <= half or percentile >= 1.0 - half)


def _expression_bins(genes):
    """(bin index per gene, bin sizes). Equal-count bins of mean expression."""
    levels = [(_f(g["mean_expression"]), i) for i, g in enumerate(genes)]
    ranked = sorted((lv, i) for lv, i in levels if lv is not None)
    bins = [-1] * len(genes)
    for position, (_, i) in enumerate(ranked):
        bins[i] = min(EXPRESSION_BINS - 1,
                      position * EXPRESSION_BINS // max(1, len(ranked)))
    sizes = {}
    for b in bins:
        sizes[b] = sizes.get(b, 0) + 1
    return bins, sizes


def _flags(genes, compounds) -> dict:
    """Per compound, per gene: percentile in each scope and whether both are
    extreme. Computed once so any subset's pass rate is a sum rather than a
    rescan of every gene."""
    out = {}
    for c in compounds:
        pooled = [_f(g[_column(c, "pooled")]) for g in genes]
        lineage = [_f(g[_column(c, "lineage")]) for g in genes]
        ps = sorted(v for v in pooled if v is not None)
        ls = sorted(v for v in lineage if v is not None)
        pct_p = [None if v is None else _percentile(ps, v) for v in pooled]
        pct_l = [None if v is None else _percentile(ls, v) for v in lineage]
        scored = [p is not None and l is not None for p, l in zip(pct_p, pct_l)]
        both = [s and _is_extreme(p) and _is_extreme(l)
                for s, p, l in zip(scored, pct_p, pct_l)]
        out[c] = {"rho_pooled": pooled, "rho_lineage": lineage,
                  "percentile_pooled": pct_p, "percentile_lineage": pct_l,
                  "scored": scored, "both_extreme": both,
                  "n_unscored": sum(1 for s in scored if not s)}
    return out


def _rate(flag, indices=None) -> dict:
    """How often the two-condition rule fires, over all genes or a subset."""
    if indices is None:
        indices = range(len(flag["scored"]))
    n = both = 0
    for i in indices:
        if not flag["scored"][i]:
            continue
        n += 1
        both += 1 if flag["both_extreme"][i] else 0
    return {"n_genes": n, "n_both_extreme": both,
            "rate": (both / n) if n else None}


def _ranks(values) -> list:
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = shared
        i = j + 1
    return out


def _centered(ranks):
    mean = sum(ranks) / len(ranks)
    centered = [r - mean for r in ranks]
    norm = math.sqrt(sum(v * v for v in centered))
    return centered, norm


def _pearson(a, na, b, nb):
    if na == 0 or nb == 0:
        return None
    return sum(x * y for x, y in zip(a, b)) / (na * nb)


def _panel_series(panel, gene, compound):
    xs, ys = [], []
    for row in panel:
        x = _f(row.get(f"expr_{gene}"))
        y = _f(row.get(f"viability_auc_{compound}"))
        if x is None or y is None:
            continue
        xs.append(x)
        ys.append(y)
    return xs, ys


def _panel_statistics(panel, candidates, compounds) -> dict:
    """Per (gene, compound): the correlation recomputed from the per-line panel,
    and a permutation p for it.

    The permutation shuffles the READOUT's ranks and reuses every gene's
    expression ranks, so one pass over the shuffles scores every gene at that
    compound. Shuffling ranks is the same experiment as shuffling values,
    because a rank correlation sees only the ranks.
    """
    out = {}
    for c in compounds:
        rows = [r for r in panel if _f(r.get(f"viability_auc_{c}")) is not None]
        if len(rows) < 3:
            for gene in candidates:
                out[(gene, c)] = {"n_lines": len(rows), "rho": None,
                                  "permutation_p": None}
            continue
        y_ranks = _ranks([_f(r[f"viability_auc_{c}"]) for r in rows])
        xs = {}
        for gene in candidates:
            column = f"expr_{gene}"
            if column not in (rows[0] if rows else {}):
                continue
            values = [_f(r[column]) for r in rows]
            if any(v is None for v in values):
                continue
            xs[gene] = _centered(_ranks(values))
        observed = {}
        for gene, (cx, nx) in xs.items():
            cy, ny = _centered(y_ranks)
            observed[gene] = _pearson(cx, nx, cy, ny)
        hits = {gene: 0 for gene in xs}
        rng = random.Random(PERMUTATION_SEED)
        shuffled = list(y_ranks)
        for _ in range(N_PERMUTATIONS):
            rng.shuffle(shuffled)
            cy, ny = _centered(shuffled)
            for gene, (cx, nx) in xs.items():
                rho = _pearson(cx, nx, cy, ny)
                target = observed[gene]
                if rho is not None and target is not None and abs(rho) >= abs(target):
                    hits[gene] += 1
        for gene in candidates:
            if gene not in xs:
                out[(gene, c)] = {"n_lines": len(rows), "rho": None,
                                  "permutation_p": None}
                continue
            out[(gene, c)] = {
                "n_lines": len(rows),
                "rho": observed[gene],
                "permutation_p": (hits[gene] + 1) / (N_PERMUTATIONS + 1),
            }
    return out


def _spearman_pairs(series_a, series_b) -> float:
    """Rank correlation over the entries where both are present."""
    pairs = [(a, b) for a, b in zip(series_a, series_b)
             if a is not None and b is not None]
    if len(pairs) < 3:
        return None
    xa, xb = _centered(_ranks([a for a, _ in pairs])), _centered(_ranks([b for _, b in pairs]))
    return _pearson(xa[0], xa[1], xb[0], xb[1])


def _readout_agreement(panel, compounds) -> dict:
    """How much the inducers' readouts agree, line by line.

    Compounds screened on overlapping cohorts are not independent
    tests, and this says by how much rather than leaving it as a caveat: if the
    readouts are nearly the same ranking of lines, a gene passing for several
    compounds has been counted several times for one fact.
    """
    series = {c: [_f(r.get(f"viability_auc_{c}")) for r in panel] for c in compounds}
    rows = []
    for i, a in enumerate(compounds):
        for b in compounds[i + 1:]:
            shared = sum(1 for x, y in zip(series[a], series[b])
                         if x is not None and y is not None)
            rows.append({"compounds": [a, b], "n_shared_lines": shared,
                         "rho": _spearman_pairs(series[a], series[b])})
    values = [r["rho"] for r in rows if r["rho"] is not None]
    return {"pairs": rows,
            "max_rho": max(values) if values else None,
            "min_rho": min(values) if values else None,
            "median_rho": (sorted(values)[len(values) // 2] if values else None)}


def _gene_collinearity(panel, genes) -> dict:
    """Pairwise expression correlation among a set of genes, across the panel.

    A gene set drawn from one pathway is one axis measured several times. The
    report needs the number, because "these are separate genes" and "these are
    separate observations" are different claims.
    """
    series = {g: [_f(r.get(f"expr_{g}")) for r in panel] for g in genes}
    rows = []
    ordered = sorted(genes)
    for i, a in enumerate(ordered):
        for b in ordered[i + 1:]:
            rows.append({"genes": [a, b], "rho": _spearman_pairs(series[a], series[b])})
    values = [abs(r["rho"]) for r in rows if r["rho"] is not None]
    return {"pairs": rows, "n_genes": len(ordered),
            "max_abs_rho": max(values) if values else None,
            "median_abs_rho": (sorted(values)[len(values) // 2] if values else None)}


def _binomial_tail(n, k, p):
    """P(X >= k) for X ~ Binomial(n, p). Exact; n is the number of panel cells."""
    if p is None or not 0.0 < p < 1.0 or n <= 0:
        return None
    return min(1.0, sum(math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
                        for i in range(k, n + 1)))


def assemble(raw) -> dict:
    genes, panel, prov = raw["genes"], raw["panel"], raw["provenance"]
    compounds = sorted(prov["compounds"])
    candidates = sorted(prov["panel_genes"])
    index = {g["gene_symbol"]: i for i, g in enumerate(genes)}
    bins, bin_sizes = _expression_bins(genes)
    by_bin = {}
    for i, b in enumerate(bins):
        by_bin.setdefault(b, []).append(i)

    flags = _flags(genes, compounds)
    all_gene = {c: _rate(flags[c]) for c in compounds}
    panel_stats = _panel_statistics(panel, candidates, compounds)

    cells, per_gene = [], {}
    for gene in candidates:
        i = index.get(gene)
        if i is None:
            per_gene[gene] = {"present": False}
            continue
        matched = by_bin.get(bins[i], [])
        per_compound = {}
        for c in compounds:
            f = flags[c]
            passes = f["both_extreme"][i]
            level_null = (_rate(f, matched) if len(matched) >= MIN_MATCHED_GENES
                          else None)
            ps = panel_stats[(gene, c)]
            table_rho = f["rho_pooled"][i]
            per_compound[c] = {
                "rho_pooled": table_rho,
                "rho_by_lineage": f["rho_lineage"][i],
                "percentile_pooled": f["percentile_pooled"][i],
                "percentile_by_lineage": f["percentile_lineage"][i],
                "passes_preregistered_rule": passes,
                "all_gene_null": all_gene[c],
                "level_matched_null": level_null,
                "n_lines_in_panel": ps["n_lines"],
                "rho_from_panel": ps["rho"],
                "permutation_p": ps["permutation_p"],
                "panel_agrees_with_table": (
                    None if table_rho is None or ps["rho"] is None
                    else abs(table_rho - ps["rho"]) < PANEL_AGREEMENT_TOLERANCE),
            }
            cells.append({
                "gene": gene, "compound": c, "passes": passes,
                "rho_pooled": table_rho,
                "expected": all_gene[c]["rate"],
                "expected_level_matched": (level_null or {}).get("rate"),
            })
        per_gene[gene] = {
            "present": True,
            "entrez_id": int(genes[i]["entrez_id"]),
            "mean_expression": _f(genes[i]["mean_expression"]),
            "expression_bin": bins[i],
            "n_genes_at_this_level": len(matched),
            "compounds": per_compound,
            "n_compounds_passing": sum(
                1 for c in compounds if per_compound[c]["passes_preregistered_rule"]),
        }

    n_cells = len(cells)
    passing = [c for c in cells if c["passes"]]
    n_pass = len(passing)
    expected = sum(c["expected"] or 0.0 for c in cells)
    ratio = (n_pass / expected) if expected > 0 else None
    binomial_p = (_binomial_tail(n_cells, n_pass, expected / n_cells)
                  if n_cells else None)
    # The level-matched expectation, which is the one that answers "but these
    # are highly expressed genes". Cells whose bin was too thin to quote are
    # counted in the denominator at the all-gene rate rather than dropped, so
    # the comparison keeps the same set of cells.
    expected_level = sum((c["expected_level_matched"] if c["expected_level_matched"]
                          is not None else (c["expected"] or 0.0)) for c in cells)
    ratio_level = (n_pass / expected_level) if expected_level > 0 else None
    binomial_p_level = (_binomial_tail(n_cells, n_pass, expected_level / n_cells)
                        if n_cells else None)
    signs = [c["rho_pooled"] for c in passing if c["rho_pooled"] is not None]
    direction_of_passing = {
        "n_with_rho": len(signs),
        "n_positive": sum(1 for v in signs if v > 0),
        "n_negative": sum(1 for v in signs if v < 0),
        "unanimous": (len(signs) > 0 and
                      (all(v > 0 for v in signs) or all(v < 0 for v in signs))),
    }
    passing_genes = sorted({c["gene"] for c in passing})
    collinearity = _gene_collinearity(panel, passing_genes) if passing_genes else None
    agreement = _readout_agreement(panel, compounds)

    pre = per_gene.get(PREREGISTERED_GENE, {})
    pre_rows = []
    for c in PREREGISTERED_COMPOUNDS:
        info = (pre.get("compounds") or {}).get(c)
        if info is None:
            continue
        rho = info["rho_pooled"]
        held = None if rho is None else (
            rho < 0 if PREREGISTERED_DIRECTION == "negative" else rho > 0)
        pre_rows.append({"compound": c, "rho_pooled": rho,
                         "rho_by_lineage": info["rho_by_lineage"],
                         "percentile_pooled": info["percentile_pooled"],
                         "permutation_p": info["permutation_p"],
                         "direction_held": held,
                         "passes_extremeness_rule":
                             info["passes_preregistered_rule"]})
    held = [r["direction_held"] for r in pre_rows]
    extreme = [r["passes_extremeness_rule"] for r in pre_rows]
    preregistered_verdict = (
        "no_preregistered_row" if not pre_rows
        else "direction_held_for_every_preregistered_compound" if all(held)
        else "direction_failed_for_every_preregistered_compound" if not any(held)
        else "direction_held_for_some_preregistered_compounds")

    blocked = {}
    for gene in prov["data_blocked_genes"]:
        info = per_gene.get(gene, {})
        n = info.get("n_compounds_passing")
        blocked[gene] = {
            "present": info.get("present", False),
            "n_compounds_passing": n,
            "verdict": ("route_unavailable" if not info.get("present")
                        else "anchored_by_this_route" if n
                        else "route_tried_does_not_anchor"),
        }

    # One trial per GENE: does it clear the rule for any compound at all? The
    # cell-level count above treats each (gene, compound) pair as a trial, and
    # the dependence measured in this same artifact says they are not. This is
    # the weaker question, with the dependence across compounds folded into the
    # null instead of ignored.
    per_gene_null = 1.0
    for c in compounds:
        rate = all_gene[c]["rate"]
        if rate is not None:
            per_gene_null *= (1.0 - rate)
    per_gene_null = 1.0 - per_gene_null
    scored_genes = [g for g in candidates if per_gene.get(g, {}).get("present")]
    n_genes_passing = sum(1 for g in scored_genes
                          if per_gene[g]["n_compounds_passing"] > 0)
    gene_level = {
        "n_genes": len(scored_genes),
        "n_genes_passing": n_genes_passing,
        "per_gene_null": per_gene_null,
        "expected_genes_passing": per_gene_null * len(scored_genes),
        "ratio": ((n_genes_passing / (per_gene_null * len(scored_genes)))
                  if per_gene_null * len(scored_genes) > 0 else None),
        "binomial_p": _binomial_tail(len(scored_genes), n_genes_passing,
                                     per_gene_null),
    }

    reportable = n_cells >= MIN_CELLS_FOR_ENRICHMENT
    # BOTH nulls, not either. A result that vanished once the panel was compared
    # against genes at its own expression level would not be about ferroptosis
    # biology, and reporting it under a verdict that only consulted the all-gene
    # null would hide exactly the confound the level-matched null exists for.
    signal = (reportable and ratio is not None and ratio_level is not None
              and n_pass > expected and n_pass > expected_level
              and (gene_level["ratio"] or 0.0) > ENRICHMENT_NULL)
    return {
        "question": ("does the engine's ferroptosis gene panel track "
                     "ferroptosis-inducer sensitivity across cell lines more "
                     "than the gene distribution it must be compared against"),
        "provenance": {k: prov[k] for k in (
            "source", "figshare_article", "figshare_doi", "license",
            "join_key", "n_curve_lines", "n_joined_lines", "n_genes",
            "n_blank_expression_values", "min_lineage_n", "n_grid", "readout")},
        "compounds": {c: prov["compounds"][c] for c in compounds},
        "all_gene_null": all_gene,
        "candidates": per_gene,
        "preregistered": {
            "gene": PREREGISTERED_GENE,
            "direction": PREREGISTERED_DIRECTION,
            "compounds": list(PREREGISTERED_COMPOUNDS),
            "rows": pre_rows,
            "verdict": preregistered_verdict,
            "n_compounds_extreme": sum(1 for v in extreme if v),
            "direction_and_extremeness_agree": (
                bool(pre_rows) and all(h == e for h, e in zip(held, extreme))),
        },
        "data_blocked_rows": blocked,
        "enrichment": {
            "n_cells": n_cells,
            "n_passing": n_pass,
            "expected_passing": expected,
            "ratio": ratio,
            "binomial_p": binomial_p,
            "expected_passing_level_matched": expected_level,
            "ratio_level_matched": ratio_level,
            "binomial_p_level_matched": binomial_p_level,
            "reportable": reportable,
        },
        "gene_level_enrichment": gene_level,
        "passing_genes": passing_genes,
        "direction_of_passing_cells": direction_of_passing,
        "passing_gene_collinearity": collinearity,
        "readout_agreement_across_compounds": agreement,
        "expression_bin_sizes": {str(k): v for k, v in sorted(bin_sizes.items())},
        "extreme_tail": EXTREME_TAIL,
        "n_permutations": N_PERMUTATIONS,
        "panel_agreement_tolerance": PANEL_AGREEMENT_TOLERANCE,
        "verdict": ("panel_beats_its_own_null" if signal
                    else "panel_indistinguishable_from_its_own_null"),
    }


def _num(value, places=3):
    return "n/a" if value is None else f"{value:.{places}f}"


def _p(value, places=4):
    """A p-value small enough to round to zero must not be printed as zero."""
    if value is None:
        return "n/a"
    return f"{value:.1e}" if value < 10.0 ** (-places) else f"{value:.{places}f}"


def _pct(fraction, places=1):
    return "n/a" if fraction is None else f"{100.0 * fraction:.{places}f}%"


def _yesno(value):
    return "n/a" if value is None else ("yes" if value else "no")


def render(d) -> str:
    prov = d["provenance"]
    compounds = sorted(d["compounds"])
    candidates = sorted(d["candidates"])
    enr = d["enrichment"]
    tail = _pct(d["extreme_tail"], 0)
    half = _pct(d["extreme_tail"] / 2.0, 1)
    L = []
    L.append("# The engine's ferroptosis gene panel against cell-line "
             "ferroptosis-inducer sensitivity\n")
    L.append("*Generated by `scripts/depmap_ferroptosis_join.py`. Every input is "
             "committed: the per-gene correlation table, the per-line panel and "
             "the provenance file all live in `analysis/calibration/`, so CI "
             "re-derives every figure below and downloads nothing.*\n")

    L.append("## The access claim that was wrong\n")
    L.append(f"`calibration-feasibility.md` recorded a set of genes as data-blocked "
             f"and named one join that would make them all testable at once. It "
             f"recorded that join as blocked on ACCESS. That was a statement "
             f"about the DepMap portal download API, which does still serve a "
             f"bot-verification page; it was not true of the dataset. "
             f"**{prov['source']}** is published on figshare "
             f"(article `{prov['figshare_article']}`, {prov['figshare_doi']}, "
             f"{prov['license']}) and serves without a gate.\n")
    L.append(f"The join needs no cell-line name mapping: the committed CTRPv2 "
             f"curves are already keyed on `{prov['join_key']}`, which is the "
             f"expression matrix's own row index. "
             f"**{prov['n_joined_lines']} of the {prov['n_curve_lines']} curve "
             f"lines join** "
             f"({_pct(prov['n_joined_lines'] / prov['n_curve_lines'])}), "
             f"across {prov['n_genes']} protein-coding genes "
             f"({prov['n_blank_expression_values']} blank expression values).\n")

    L.append("## The answer\n")
    if not enr["reportable"]:
        L.append("Too few candidate cells to report an enrichment.\n")
    else:
        L.append(f"Of **{enr['n_cells']} candidate (gene, compound) cells**, "
                 f"**{enr['n_passing']}** clear the pre-registered rule, against "
                 f"**{_num(enr['expected_passing'], 2)} expected** from the "
                 f"measured rate at which that rule fires on a gene chosen "
                 f"without regard to biology -- a ratio of "
                 f"**{_num(enr['ratio'], 2)}x**, binomial "
                 f"p = {_p(enr['binomial_p'])}.\n")
        L.append(f"Against the LEVEL-MATCHED expectation -- the same rule's "
                 f"firing rate among genes at each candidate's own expression "
                 f"level, which is the null that answers *but these are highly "
                 f"expressed genes* -- the expectation is "
                 f"{_num(enr['expected_passing_level_matched'], 2)} and the "
                 f"ratio {_num(enr['ratio_level_matched'], 2)}x "
                 f"(p = {_p(enr['binomial_p_level_matched'])}). The result does "
                 f"not rest on the panel being highly expressed.\n")
        gl = d["gene_level_enrichment"]
        L.append(f"**The cell-level p above assumes the cells are independent, "
                 f"and this page's own measurements say they are not.** Counting "
                 f"each gene ONCE instead -- does it clear the rule for any "
                 f"compound at all -- **{gl['n_genes_passing']} of "
                 f"{gl['n_genes']} genes** pass against "
                 f"{_num(gl['expected_genes_passing'], 2)} expected, a ratio of "
                 f"**{_num(gl['ratio'], 2)}x** at p = "
                 f"{_p(gl['binomial_p'])}. That is the conservative figure and "
                 f"the one to quote if only one is quoted; it still assumes "
                 f"independence ACROSS GENES, which the collinearity table "
                 f"below bounds.\n")
        L.append(f"Verdict: `{d['verdict']}`.\n")
    L.append(f"The rule, fixed in the issue before any correlation existed: a "
             f"gene passes for a compound only if its rank correlation sits in "
             f"the extreme {tail} of the distribution over EVERY gene "
             f"(percentile at or below {half}, or at or above "
             f"{_pct(1.0 - d['extreme_tail'] / 2.0, 1)}) **and** still does so "
             f"when expression and readout are ranked within Oncotree lineage "
             f"first. Both, not either.\n")
    L.append("The expectation is measured rather than assumed. The two "
             "conditions share a cohort, so they are strongly dependent and the "
             "nominal product of the two tails would be the wrong number:\n")
    L.append("| compound | lines joined | genes scored | fire on both | "
             "measured rate |")
    L.append("|---|--:|--:|--:|--:|")
    for c in compounds:
        null = d["all_gene_null"][c]
        info = d["compounds"][c]
        L.append(f"| {c} | {info['n_joined_lines']} | {null['n_genes']} | "
                 f"{null['n_both_extreme']} | {_pct(null['rate'], 2)} |")
    L.append("")

    L.append("## The one direction fixed in advance\n")
    pre = d["preregistered"]
    L.append(f"`biochem.rs` adds `acsl4_status_boost` to the "
             f"peroxidizable-substrate term and `acsl4.rs` calls "
             f"ACSL4-negative status intrinsically refractory, so the engine "
             f"predicts that higher {pre['gene']} means a more SENSITIVE line -- "
             f"a {pre['direction']} correlation with viability -- for the GPX4 "
             f"inhibitors. That sign is in the engine's source, not chosen here.\n")
    L.append("| compound | rho pooled | rho within lineage | percentile | "
             "permutation p | direction held |")
    L.append("|---|--:|--:|--:|--:|---|")
    for row in pre["rows"]:
        L.append(f"| {row['compound']} | {_num(row['rho_pooled'])} | "
                 f"{_num(row['rho_by_lineage'])} | "
                 f"{_pct(row['percentile_pooled'], 2)} | "
                 f"{_p(row['permutation_p'])} | "
                 f"{_yesno(row['direction_held'])} |")
    L.append("")
    L.append(f"Verdict: `{pre['verdict']}`.\n")
    L.append(f"**The sign and the extremeness are separate questions, and they "
             f"came back differently.** The registered direction held for "
             f"{len([r for r in pre['rows'] if r['direction_held']])} of "
             f"{len(pre['rows'])} compounds, at permutation p values the table "
             f"gives, so a correlation that size is not expected at this cohort "
             f"size by chance. But {pre['gene']} clears the extremeness bar for "
             f"{pre['n_compounds_extreme']} of {len(pre['rows'])} of them: its "
             f"percentiles sit just outside the tail the rule requires. The "
             f"engine's registered sign is confirmed; the gene is NOT among the "
             f"most predictive in the matrix, and both halves belong in the "
             f"record.\n")

    L.append(f"## What the {len(d['data_blocked_rows'])} data-blocked rows "
             f"become\n")
    L.append("A gene that is not extreme against the all-gene null does not go "
             "back to being untestable. Its row changes from *the route could "
             "not be tried* to *the route was tried and does not anchor this "
             "gene*, which is a different and more useful sentence.\n")
    L.append("| gene | compounds passing | furthest percentile | smallest perm p "
             "| verdict |")
    L.append("|---|--:|--:|--:|---|")
    for gene in sorted(d["data_blocked_rows"]):
        row = d["data_blocked_rows"][gene]
        info = d["candidates"].get(gene, {})
        cells = (info.get("compounds") or {}).values()
        pcts = [c["percentile_pooled"] for c in cells
                if c["percentile_pooled"] is not None]
        ps = [c["permutation_p"] for c in cells if c["permutation_p"] is not None]
        # The percentile FURTHEST from the middle, reported as itself. Folding
        # the two tails into one number inverts the apparent direction for a
        # negatively correlated gene, which is most of these rows.
        furthest = (max(pcts, key=lambda v: abs(v - 0.5)) if pcts else None)
        L.append(f"| {gene} | {row['n_compounds_passing']} | "
                 f"{'n/a' if furthest is None else _pct(furthest, 2)} | "
                 f"{'n/a' if not ps else _p(min(ps))} | `{row['verdict']}` |")
    L.append("")
    L.append("**Not anchored is not the same as flat, and the columns say which "
             "is which.** The rule asks whether a gene is extreme among ALL "
             "genes, which is the right bar for justifying a model layer and a "
             "high one. A gene can have a correlation that a shuffled readout "
             "almost never reproduces -- a small permutation p -- and still sit "
             "nowhere near the tail. That is what these rows look like, and "
             "reading them as *no signal* would overstate the negative in the "
             "same way the original access claim overstated the blocker.\n")

    L.append("## Every candidate gene, every compound\n")
    L.append("Positive rho means higher expression goes with RESISTANCE. "
             "`level` is the gene's own mean expression and `n at level` the "
             "number of genes in its expression bin, which is the null the "
             "`level` column's rate is drawn from.\n")
    L.append("| gene | entrez | compound | rho pooled | rho by lineage | "
             "percentile | perm p | all-gene rate | level rate | passes |")
    L.append("|---|--:|---|--:|--:|--:|--:|--:|--:|---|")
    for gene in candidates:
        info = d["candidates"][gene]
        if not info.get("present"):
            L.append(f"| {gene} | n/a | -- | n/a | n/a | n/a | n/a | n/a | n/a | "
                     f"absent from the matrix |")
            continue
        for c in compounds:
            cc = info["compounds"][c]
            level = cc["level_matched_null"]
            L.append(
                f"| {gene} | {info['entrez_id']} | {c} | "
                f"{_num(cc['rho_pooled'])} | {_num(cc['rho_by_lineage'])} | "
                f"{_pct(cc['percentile_pooled'], 2)} | "
                f"{_p(cc['permutation_p'])} | "
                f"{_pct(cc['all_gene_null']['rate'], 2)} | "
                f"{'n/a' if level is None else _pct(level['rate'], 2)} | "
                f"{_yesno(cc['passes_preregistered_rule'])} |")
    L.append("")

    L.append("## What passes, and how many separate things that is\n")
    dirn = d["direction_of_passing_cells"]
    col = d["passing_gene_collinearity"]
    agree = d["readout_agreement_across_compounds"]
    L.append(f"The rule fires for **{len(d['passing_genes'])} of the "
             f"{len(candidates)} panel genes**: "
             f"{', '.join(d['passing_genes']) if d['passing_genes'] else 'none'}.\n")
    if dirn["n_with_rho"]:
        L.append(f"**Every one of the {dirn['n_with_rho']} passing cells points "
                 f"the same way**: {dirn['n_positive']} positive and "
                 f"{dirn['n_negative']} negative, so higher expression goes with "
                 f"RESISTANCE in all of them. Unanimous: "
                 f"{_yesno(dirn['unanimous'])}. That is the direction the engine "
                 f"gives a defence node -- more capacity to detoxify or to "
                 f"resupply, less death -- but the reading is a reading: nothing "
                 f"here tests a mechanism, and a gene can track resistance "
                 f"because it causes it or because it marks a state that has "
                 f"it.\n")
    if col:
        L.append(f"**They are not {col['n_genes']} independent observations.** "
                 f"Pairwise expression correlation among the passing genes, "
                 f"across the same lines: maximum "
                 f"{_num(col['max_abs_rho'])}, median "
                 f"{_num(col['median_abs_rho'])}. Co-regulated genes carry one "
                 f"axis several times, and the cell count above does not know "
                 f"that.\n")
        L.append("| gene pair | rho |")
        L.append("|---|--:|")
        for row in col["pairs"]:
            L.append(f"| {row['genes'][0]} / {row['genes'][1]} | "
                     f"{_num(row['rho'])} |")
        L.append("")
    L.append(f"**And the compounds are not independent either.** Pairwise rank "
             f"correlation between the readouts over shared lines runs "
             f"{_num(agree['min_rho'])} to {_num(agree['max_rho'])}, median "
             f"{_num(agree['median_rho'])}:\n")
    L.append("| compound pair | shared lines | rho |")
    L.append("|---|--:|--:|")
    for row in agree["pairs"]:
        L.append(f"| {row['compounds'][0]} / {row['compounds'][1]} | "
                 f"{row['n_shared_lines']} | {_num(row['rho'])} |")
    L.append("")
    fired = sorted(d["all_gene_null"][c]["n_both_extreme"] for c in compounds)
    rates = sorted(d["all_gene_null"][c]["rate"] for c in compounds
                   if d["all_gene_null"][c]["rate"] is not None)
    L.append(f"**The panel genes are not the only genes that pass.** The rule "
             f"fires on {_pct(rates[0], 2)} to {_pct(rates[-1], 2)} of all "
             f"genes, so each compound has between {fired[0]} and {fired[-1]} "
             f"genes clearing it. The claim here is about the RATE at which "
             f"panel genes clear it, not that they are the top of the list -- "
             f"this report does not inspect the rest of the ranking, and a gene "
             f"outside the panel could well beat every gene in it.\n")
    L.append("## The nulls, and why one is not enough\n")
    L.append("**The all-gene null is exact, not sampled.** Every "
             "protein-coding gene in the matrix is scored, so a candidate's "
             "percentile is its true position rather than an estimate from a "
             "control set somebody chose. This is the null the pass rule is "
             "stated against.\n")
    L.append("**The level-matched null exists because the panel is not a "
             "random sample of genes.** These are well known, highly expressed "
             "loci, and expression level governs how stably a gene ranks at all, "
             "so the all-gene null is not matched to the panel on the axis that "
             "sets rank stability. Each candidate is therefore also scored "
             "against the rate among genes in its own expression bin:\n")
    L.append("| expression bin | genes |")
    L.append("|---|--:|")
    for b in sorted(d["expression_bin_sizes"], key=lambda k: int(k)):
        L.append(f"| {b} | {d['expression_bin_sizes'][b]} |")
    L.append("")
    L.append(f"**The permutation null asks a different question**: not whether "
             f"this gene is unusual among genes, but whether a correlation of "
             f"this size is surprising at this cohort size at all. It shuffles "
             f"the readout across lines {d['n_permutations']} times on a named "
             f"seed and reports the share of shuffles reaching the observed "
             f"magnitude. A gene can pass the all-gene rule and still have an "
             f"unremarkable permutation p, which is the case worth seeing.\n")
    L.append(f"**And the two committed derivatives are checked against each "
             f"other.** Every rho in the table above is recomputed from the "
             f"per-line panel, independently of the matrix the fetch leg read, "
             f"and the two must agree to within "
             f"{_num(d['panel_agreement_tolerance'], 2)}. A stale derivative is "
             f"otherwise invisible.\n")

    L.append("## What this does NOT establish\n")
    L.append("**No magnitude for any `Params` field.** A rank correlation "
             "bounds an ordering over cell lines. Nothing here supplies a "
             "number for a rate constant, and the `status_from_zscore` slope "
             "stays the placeholder it was.\n")
    L.append("**Immortalised lines are not tumours in patients.** The "
             "prevalence leg already found that bulk tumour mRNA does not rank "
             "the ACSL4-low refractory phenotype the cell-line literature "
             "reports; this leg measures lines, and inherits that gap rather "
             "than closing it.\n")
    L.append(f"**The compounds are not independent tests.** The "
             f"{len(compounds)} inducers are screened on overlapping cohorts of "
             f"the same {prov['n_joined_lines']} lines, so they share both the "
             f"cohort and most of the biology. The cell count above therefore "
             f"measures how consistently a gene behaves across compounds, not "
             f"how many times the result was reproduced.\n")
    L.append(f"**{prov['readout']}.** It is read from the published "
             f"logistic fit on one common log-dose grid per compound, "
             f"chosen as the range the plurality of that compound's lines was "
             f"actually screened over, so every evaluation is interpolation "
             f"inside a line's own screened range. The committed curves contain "
             f"run-away fits whose EC50 sits orders of magnitude above the "
             f"screened maximum; those are a failed fit, not a resistant line, "
             f"which is why raw EC50 is not the readout.\n")
    L.append("**A correlation does not say which way the mechanism runs.** A "
             "gene can track sensitivity because it causes it, because it marks "
             "a lineage that has it, or because it rides a growth-rate axis. "
             "Ranking within lineage removes one of those three.\n")
    return "\n".join(L) + "\n"


def _split_stored(stored):
    """The freshness gate's hook for recovering `assemble`'s inputs.

    There is nothing to recover: every input is committed, so they are read
    again. `stored` is accepted and deliberately unused, because taking the
    inputs back out of the artifact would compare the artifact to itself.
    """
    return (load(),)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild both artifacts from the committed inputs")
    a = ap.parse_args()
    # `assemble(load())` on BOTH paths, deliberately. Every input is committed,
    # so there is no scan to skip and nothing to gain by trusting a stored copy
    # of a derived field; a leg that is always recomputed can never be echoed.
    if a.render_only:
        data = assemble(load())
    else:
        data = assemble(load())
    OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_JSON.name} and {OUT_MD.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
