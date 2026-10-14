#!/usr/bin/env python3
"""Join DepMap cell-line expression onto the committed CTRPv2 ferroptosis curves (#932).

WHY
---
`analysis/calibration/calibration-feasibility.md` recorded four genes as
data-blocked and named exactly one join that would make all four testable at
once: cell-line expression against ferroptosis-inducer sensitivity across the
committed CTRPv2 cohort. It recorded that join as blocked on ACCESS, because the
DepMap portal download API began serving a bot-verification page where the
documented CSV had been. That observation remains true of the portal API and is
false of the dataset: DepMap publishes every public release on figshare under
CC BY 4.0, served without a gate, and the committed curves are already keyed on
the matrix's own row index, so no cell-line name mapping is involved.

WHAT IT DOES (provenance tool, NOT run in CI)
---------------------------------------------
1. Resolve the figshare article's file list, download the expression matrix and
   the model metadata into a raw cache directory that is NOT committed, and
   md5-verify both against figshare's own supplied hashes.
2. Read sensitivity INSIDE each line's recorded dose range, on one common dose
   grid per compound, selected by `ctrp_dose_support.select_supported_cohort`.
   The committed curves contain run-away fits whose EC50 sits orders of
   magnitude above the screened maximum; those are a failed fit rather than a
   resistant line, so raw EC50 is not the readout.
3. Join on ModelID and compute, for EVERY protein-coding gene in the matrix, the
   rank correlation between its expression and that readout -- pooled, and again
   ranked within Oncotree lineage so the lineage confound is measured.
4. Write three small committed derivatives plus provenance.

WHY THE STATISTIC AND NOT THE MATRIX
------------------------------------
The null this analysis needs is the complete distribution over all genes, so a
candidate gene's percentile in it is exact rather than sampled. The matrix is
hundreds of megabytes; the per-gene correlation is the same information for that
purpose and keeps the existing offline contract, under which CI re-derives every
reported figure from committed files and downloads nothing.

SIGN CONVENTION
---------------
The readout is mean predicted VIABILITY across the common dose grid, so a low
value is a sensitive line. A POSITIVE correlation therefore means higher
expression goes with RESISTANCE. Every stored column says so in its name.

GENE IDENTITY
-------------
Genes are resolved by Entrez ID, not by symbol, and the declared symbol must
agree with the matrix's own header or the run fails. This repo has already been
bitten by a symbol-level mapping error in a third-party annotation, and a
candidate gene silently resolving to a different locus would be invisible in
every figure downstream.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

# ONE Python source for the curve formula: the fetch leg that wrote the committed
# curves also owns how to evaluate them.
from fetch_calibration_data import predicted_viability  # noqa: E402
from ctrp_dose_support import select_supported_cohort  # noqa: E402

OUT_DIR = REPO / "analysis" / "calibration"
CURVES_CSV = OUT_DIR / "ctrpv2_ferroptosis_curves.csv"
CORRELATIONS_CSV_GZ = OUT_DIR / "depmap_ctrpv2_gene_correlations.csv.gz"
PANEL_CSV = OUT_DIR / "depmap_ctrpv2_panel.csv"
PROVENANCE_JSON = OUT_DIR / "depmap_ctrpv2_join_provenance.json"

# DepMap 24Q4 Public on figshare. The portal download API is gated; this is not.
FIGSHARE_ARTICLE = 27993248
FIGSHARE_API = "https://api.figshare.com/v2/articles/{article}"
EXPRESSION_FILE = "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL_FILE = "Model.csv"

# Number of log-spaced doses in the common grid. The grid is the modal recorded
# range for the compound, so the policy is "the range the plurality of lines was
# actually screened over" rather than a number chosen to retain a cohort.
N_GRID = 64

# The seven genes the cBioPortal prevalence leg measured, carried forward with
# their Entrez IDs so this leg asks about the same loci, plus the engine's own
# lipid-peroxidation and iron nodes.
PANEL_GENES = {
    # measured by the prevalence leg (#462, #616): the four data-blocked rows
    # are HMOX1 / TP53 / TFRC / KEAP1.
    "ACSL4": 2182,
    "GPX4": 2879,
    "SLC7A11": 23657,
    "HMOX1": 3162,
    "TP53": 7157,
    "TFRC": 7037,
    "KEAP1": 9817,
    # engine nodes not previously measured against any data.
    "AIFM2": 84883,   # FSP1
    "LPCAT3": 10162,
    "ALOX15": 246,
    "ALOX15B": 247,
    "NCOA4": 8031,
    "FTH1": 2495,
    "SCD": 6319,      # SCD1
    "ACSL3": 2181,
    "GCLC": 2729,
    "GSS": 2937,
    "SLC3A2": 6520,
    "NFE2L2": 4780,
}

# The four rows `calibration-feasibility.md` records as data-blocked.
DATA_BLOCKED_GENES = ("HMOX1", "TP53", "TFRC", "KEAP1")

# Lineages smaller than this cannot support a within-lineage rank, so they are
# dropped from the stratified statistic and counted rather than pooled in.
MIN_LINEAGE_N = 8

HEADER_RE = re.compile(r"^(?P<symbol>.+?) \((?P<entrez>\d+)\)$")


def _get_json(url: str):
    req = urllib.request.Request(url, headers={"accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def resolve_figshare_files(article: int = FIGSHARE_ARTICLE) -> dict:
    """Return {filename: {id, md5, size, download_url}} for the files we need."""
    doc = _get_json(FIGSHARE_API.format(article=article))
    wanted = {}
    for f in doc.get("files", []):
        if f["name"] in (EXPRESSION_FILE, MODEL_FILE):
            wanted[f["name"]] = {
                "id": f["id"],
                "md5": f.get("supplied_md5") or f.get("computed_md5"),
                "size": f["size"],
                "download_url": f["download_url"],
            }
    missing = {EXPRESSION_FILE, MODEL_FILE} - set(wanted)
    if missing:
        raise SystemExit(f"figshare article {article} is missing {sorted(missing)}")
    return {
        "article": article,
        "doi": doc.get("doi"),
        "title": doc.get("title"),
        "license": (doc.get("license") or {}).get("name"),
        "published_date": doc.get("published_date"),
        "files": wanted,
    }


def download_verified(url: str, md5: str, dest: Path) -> str:
    """Download to `dest` unless it is already present with the right md5."""
    if dest.exists():
        got = hashlib.md5(dest.read_bytes()).hexdigest()
        if got == md5:
            print(f"cached, md5 verified: {dest.name}")
            return got
        print(f"cached copy of {dest.name} has md5 {got}, expected {md5}; refetching")
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {dest.name} ...")
    with urllib.request.urlopen(url, timeout=3600) as r, dest.open("wb") as out:
        while chunk := r.read(1 << 20):
            out.write(chunk)
    got = hashlib.md5(dest.read_bytes()).hexdigest()
    if md5 and got != md5:
        raise SystemExit(f"md5 mismatch for {dest.name}: figshare {md5} != download {got}")
    print(f"md5 verified: {got}")
    return got


def modal_range(rows) -> tuple:
    """The (minimum, maximum) recorded dose the plurality of this compound's
    lines was screened over. Ties break on the wider range, then on the lower
    minimum, so the choice does not depend on row order."""
    counts = {}
    for r in rows:
        key = (float(r["MinimumDose"]), float(r["MaximumDose"]))
        counts[key] = counts.get(key, 0) + 1
    return max(counts.items(), key=lambda kv: (kv[1], kv[0][1] - kv[0][0], -kv[0][0]))[0]


def dose_grid(low: float, high: float, n: int = N_GRID) -> list:
    """`n` log-spaced doses from `low` to `high`, with the ENDPOINTS EXACT.

    Computing the last point as `10 ** (log10(low) + (log10(high) - log10(low)))`
    lands a few ulps ABOVE `high`, which the cohort selector correctly reads as a
    dose outside every curve's recorded range -- so a whole compound is excluded
    by a rounding error rather than by any property of the data. The interior
    points are log-spaced; the two the selector compares against are the
    recorded bounds themselves.
    """
    lo, hi = math.log10(low), math.log10(high)
    grid = [10 ** (lo + (hi - lo) * i / (n - 1)) for i in range(n)]
    grid[0], grid[-1] = low, high
    return grid


def readout(row, grid) -> float:
    """Mean predicted viability over `grid`. Low is sensitive."""
    lower = float(row["LowerAsymptote"])
    upper = float(row["UpperAsymptote"])
    ec50 = float(row["EC50"])
    slope = float(row["Slope"])
    return sum(predicted_viability(d, lower, upper, ec50, slope) for d in grid) / len(grid)


def compound_cohorts(curves_path: Path = CURVES_CSV) -> dict:
    """Per compound: {grid, support, readout: {ModelID: viability_auc}}."""
    rows = list(csv.DictReader(curves_path.open()))
    by_compound = {}
    for r in rows:
        by_compound.setdefault(r["CompoundName"], []).append(r)
    out = {}
    for compound, crows in sorted(by_compound.items()):
        low, high = modal_range(crows)
        grid = dose_grid(low, high)
        retained, support = select_supported_cohort(crows, grid)
        out[compound] = {
            "grid_low_um": low,
            "grid_high_um": high,
            "n_grid": len(grid),
            "n_curves": len(crows),
            "n_retained": support["n_retained"],
            "n_excluded": support["n_excluded"],
            "readout": {r["ModelID"]: readout(r, grid) for r in retained},
        }
    return out


def parse_expression(path: Path, keep_models: set):
    """Stream the matrix. Returns (model_ids, symbols, entrez, matrix, n_blank)."""
    symbols, entrez = [], []
    models, rows, blanks = [], [], 0
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split(",")
        for col in header[1:]:
            m = HEADER_RE.match(col)
            if not m:
                raise SystemExit(f"unparsed expression column header: {col!r}")
            symbols.append(m.group("symbol"))
            entrez.append(int(m.group("entrez")))
        width = len(symbols)
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if parts[0] not in keep_models:
                continue
            values = parts[1:]
            if len(values) != width:
                raise SystemExit(f"{parts[0]}: {len(values)} values for {width} genes")
            blanks += sum(1 for v in values if v == "")
            models.append(parts[0])
            rows.append(np.array([np.nan if v == "" else v for v in values],
                                 dtype=np.float32))
    return models, symbols, np.array(entrez), np.vstack(rows), blanks


def _rank_columns(matrix: np.ndarray) -> np.ndarray:
    """Average-tie ranks down each column, as float64."""
    return stats.rankdata(matrix, axis=0).astype(np.float64)


def _pearson_against(ranked: np.ndarray, target_ranks: np.ndarray) -> np.ndarray:
    """Pearson r between every column of `ranked` and `target_ranks`.

    A column with zero variance has no correlation rather than a correlation of
    zero, so it returns NaN and is counted, not reported as a null result.
    """
    x = ranked - ranked.mean(axis=0)
    y = target_ranks - target_ranks.mean()
    sx = np.sqrt((x * x).sum(axis=0))
    sy = math.sqrt(float((y * y).sum()))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (x * y[:, None]).sum(axis=0) / (sx * sy)
    r[sx == 0] = np.nan
    return r


def stratified_ranks(matrix: np.ndarray, target: np.ndarray, strata: list):
    """Rank expression and target WITHIN each stratum, then stack the strata.

    Each stratum's ranks are CENTRED ON THAT STRATUM'S OWN MEAN before
    stacking, which is the whole point of stratifying: within-stratum ranks of
    differently sized strata have different means, so stacking them raw leaves a
    between-stratum component in the correlation -- exactly the lineage effect
    this is meant to remove. Centred, the pooled Pearson over the stack is the
    within-stratum correlation, weighted by stratum variance.

    Returns (centred_matrix, centred_target, kept_strata, n_dropped_lines). A
    stratum below `MIN_LINEAGE_N` cannot support a rank and is dropped.
    """
    keep_rows, ranked_cols, ranked_target, kept = [], [], [], {}
    for name in sorted(set(strata)):
        idx = [i for i, s in enumerate(strata) if s == name]
        if len(idx) < MIN_LINEAGE_N:
            continue
        kept[name] = len(idx)
        keep_rows.extend(idx)
        block = _rank_columns(matrix[idx, :])
        ranked_cols.append(block - block.mean(axis=0))
        y = stats.rankdata(target[idx]).astype(np.float64)
        ranked_target.append(y - y.mean())
    if not ranked_cols:
        return None, None, kept, len(strata)
    return (np.vstack(ranked_cols), np.concatenate(ranked_target), kept,
            len(strata) - len(keep_rows))


def correlate(cohorts, models, matrix, lineage_of):
    """Per compound: pooled and lineage-stratified rank correlation per gene."""
    index = {m: i for i, m in enumerate(models)}
    results, meta = {}, {}
    for compound, info in cohorts.items():
        ids = [m for m in models if m in info["readout"]]
        if len(ids) < MIN_LINEAGE_N:
            raise SystemExit(f"{compound}: only {len(ids)} joined lines")
        rows = [index[m] for m in ids]
        sub = matrix[rows, :]
        target = np.array([info["readout"][m] for m in ids], dtype=np.float64)
        pooled = _pearson_against(_rank_columns(sub), stats.rankdata(target))
        strata = [lineage_of.get(m, "") for m in ids]
        rk, rt, kept, dropped = stratified_ranks(sub, target, strata)
        strat = _pearson_against(rk, rt) if rk is not None else np.full(sub.shape[1], np.nan)
        results[compound] = {"pooled": pooled, "stratified": strat}
        meta[compound] = {
            "n_joined_lines": len(ids),
            "n_constant_genes": int(np.isnan(pooled).sum()),
            "lineages_kept": kept,
            "n_lines_dropped_small_lineage": int(dropped),
            **{k: v for k, v in info.items() if k != "readout"},
        }
    return results, meta


def write_correlations(symbols, entrez, results, matrix):
    """Per gene: its correlation with each compound's readout, and its own
    expression level.

    The level is carried because a candidate gene panel is not a random sample
    of genes -- these are well known, highly expressed loci, and expression
    level sets how stably a gene ranks at all. A null drawn from ALL genes is
    therefore not matched to the panel on the axis that governs rank stability,
    and the report needs both nulls to say which one the result survives.
    """
    compounds = sorted(results)
    header = ["gene_symbol", "entrez_id", "mean_expression", "sd_expression"]
    for c in compounds:
        header += [f"{c}_rho_vs_viability_pooled", f"{c}_rho_vs_viability_by_lineage"]
    means = np.nanmean(matrix, axis=0)
    sds = np.nanstd(matrix, axis=0)
    with gzip.open(CORRELATIONS_CSV_GZ, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for i, (sym, ez) in enumerate(zip(symbols, entrez)):
            row = [sym, int(ez), f"{means[i]:.6f}", f"{sds[i]:.6f}"]
            for c in compounds:
                for key in ("pooled", "stratified"):
                    v = results[c][key][i]
                    row.append("" if not np.isfinite(v) else f"{v:.6f}")
            w.writerow(row)
    return compounds


def resolve_panel(symbols, entrez) -> dict:
    """{symbol: column index}, requiring the declared Entrez ID AND symbol to
    agree with the matrix's own header."""
    by_entrez = {int(e): i for i, e in enumerate(entrez)}
    out, problems = {}, []
    for sym, ez in PANEL_GENES.items():
        i = by_entrez.get(ez)
        if i is None:
            problems.append(f"{sym}: Entrez {ez} absent from the matrix")
            continue
        if symbols[i] != sym:
            problems.append(f"Entrez {ez} is {symbols[i]!r} in the matrix, declared {sym!r}")
            continue
        out[sym] = i
    if problems:
        raise SystemExit("panel gene identity failed:\n  " + "\n  ".join(problems))
    return out


def write_panel(models, matrix, panel, cohorts, lineage_of):
    genes = sorted(panel)
    compounds = sorted(cohorts)
    header = (["ModelID", "OncotreeLineage"]
              + [f"expr_{g}" for g in genes]
              + [f"viability_auc_{c}" for c in compounds])
    with PANEL_CSV.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for i, m in enumerate(models):
            row = [m, lineage_of.get(m, "")]
            row += [f"{matrix[i, panel[g]]:.6f}" for g in genes]
            row += [(f"{cohorts[c]['readout'][m]:.6f}"
                     if m in cohorts[c]["readout"] else "") for c in compounds]
            w.writerow(row)
    return header


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", type=Path, default=Path("/tmp"),
                    help="cache dir for the ~507MB download (NOT committed)")
    a = ap.parse_args()
    raw = a.raw_dir

    prov = resolve_figshare_files()
    print(f"{prov['title']} | {prov['doi']} | {prov['license']}")
    md5s = {}
    for name in (MODEL_FILE, EXPRESSION_FILE):
        f = prov["files"][name]
        md5s[name] = download_verified(f["download_url"], f["md5"], raw / name)

    lineage_of = {}
    with (raw / MODEL_FILE).open() as fh:
        for r in csv.DictReader(fh):
            lineage_of[r["ModelID"]] = r.get("OncotreeLineage") or ""

    cohorts = compound_cohorts()
    wanted = set().union(*(set(c["readout"]) for c in cohorts.values()))
    print(f"{len(wanted)} cell lines across {len(cohorts)} compounds")

    models, symbols, entrez, matrix, blanks = parse_expression(raw / EXPRESSION_FILE, wanted)
    print(f"joined {len(models)} lines x {len(symbols)} genes ({blanks} blank values)")

    panel = resolve_panel(symbols, entrez)
    results, meta = correlate(cohorts, models, matrix, lineage_of)
    compounds = write_correlations(symbols, entrez, results, matrix)
    write_panel(models, matrix, panel, cohorts, lineage_of)

    PROVENANCE_JSON.write_text(json.dumps({
        "source": prov["title"],
        "figshare_article": prov["article"],
        "figshare_doi": prov["doi"],
        "license": prov["license"],
        "published_date": prov["published_date"],
        "files": {n: {"md5": md5s[n], "size_bytes": prov["files"][n]["size"]}
                  for n in (MODEL_FILE, EXPRESSION_FILE)},
        "join_key": "ModelID",
        "n_curve_lines": len(wanted),
        "n_joined_lines": len(models),
        "n_genes": len(symbols),
        "n_blank_expression_values": blanks,
        "min_lineage_n": MIN_LINEAGE_N,
        "n_grid": N_GRID,
        "panel_genes": {g: PANEL_GENES[g] for g in sorted(PANEL_GENES)},
        "data_blocked_genes": list(DATA_BLOCKED_GENES),
        "compounds": {c: meta[c] for c in compounds},
        "readout": "mean predicted viability over the common log-dose grid; low is sensitive",
    }, indent=2, sort_keys=True) + "\n")
    print(f"wrote {CORRELATIONS_CSV_GZ.name}, {PANEL_CSV.name}, {PROVENANCE_JSON.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
