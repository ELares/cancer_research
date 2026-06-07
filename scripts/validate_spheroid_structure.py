#!/usr/bin/env python3
"""Validate the spheroid module's zone geometry vs Browning 2021 size-resolved data (#333).

The `ferroptosis-core` `spheroid` module assigns radial phenotype zones using
FIXED fractional thresholds taken from the Browning 2021 LIMITING (large-spheroid)
structure: the glycolytic rim begins at ~0.90 of the radius and the persister/
necrotic-like core sits inside ~0.73 of the radius (`glycolytic_frac=0.73`,
`oxphos_frac=0.39` as VOLUME fractions; 0.90 = 0.73^(1/3), 0.73 = 0.39^(1/3)).
Those thresholds are size-INDEPENDENT in the model. Real spheroids are not: small
spheroids have a thick proliferating rim and NO necrotic core; the necrotic core
appears only above a critical size and then grows toward the limiting fraction.

This script quantifies that finite-size error against the experimental confocal
structure data from Browning et al. (eLife 2021, PMID 34842141; data repo
github.com/ap-browning/Spheroids, `Data/AllConfocalData.csv`), which reports, per
spheroid, the outer radius R (µm) and the two internal boundaries as fractions of
R: phi = proliferating-rim inner boundary, eta = necrotic-core outer boundary. It
bins by R, compares the measured phi/eta to the model's fixed thresholds, and
derives the SIZE RANGE in which the fixed-threshold model is trustworthy.

Offline contract: the raw CSV is NOT committed; CI never downloads. The committed
derived summary (`analysis/calibration/spheroid-structure-validation.json`) +
report are the reproducible artifacts; re-run this script to regenerate them.
"""

import argparse
import csv
import json
import statistics
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = REPO_ROOT / "analysis" / "calibration" / "spheroid-structure-validation.json"
OUT_MD = REPO_ROOT / "analysis" / "calibration" / "spheroid-structure-validation.md"

DATA_URL = "https://raw.githubusercontent.com/ap-browning/Spheroids/master/Data/AllConfocalData.csv"

# Model fixed thresholds (radius fractions) from spheroid::SpheroidConfig::literature():
# rim begins at glycolytic_frac^(1/3); core boundary at oxphos_frac^(1/3).
MODEL_RIM_BOUNDARY = 0.73 ** (1.0 / 3.0)   # ~0.9004  (phi analog: proliferating-rim inner edge)
MODEL_CORE_BOUNDARY = 0.39 ** (1.0 / 3.0)  # ~0.7306  (eta analog: necrotic/persister-core outer edge)

# Size bins (outer radius µm). The last open bin catches the few largest.
SIZE_BINS = ((0, 200), (200, 300), (300, 400), (400, 10000))
# A bin "matches" the fixed-threshold model when both boundaries are within this
# radius-fraction tolerance AND a necrotic core actually exists in most spheroids.
MATCH_TOL = 0.15
CORE_PRESENT_FRAC = 0.5


def fetch_rows(url=DATA_URL, raw_csv=None):
    if raw_csv:
        text = Path(raw_csv).read_text()
    else:
        with urllib.request.urlopen(url, timeout=60) as r:
            text = r.read().decode("utf-8", "replace")
    rows = []
    for rec in csv.DictReader(text.splitlines()):
        try:
            rows.append({
                "cell_line": rec["CellLine"],
                "R": float(rec["R"]),
                "phi": float(rec["ϕ"]),   # ϕ
                "eta": float(rec["η"]),    # η
            })
        except (KeyError, ValueError):
            continue
    return rows


def _q(vals, q):
    s = sorted(vals)
    return s[min(len(s) - 1, int(q * len(s)))]


def summarize(rows, bins=SIZE_BINS):
    out = []
    for lo, hi in bins:
        sub = [r for r in rows if lo <= r["R"] < hi]
        if not sub:
            continue
        phi = [r["phi"] for r in sub]
        eta = [r["eta"] for r in sub]
        core_present = sum(1 for e in eta if e > 0.01) / len(sub)
        out.append({
            "r_lo_um": lo, "r_hi_um": (None if hi >= 10000 else hi),
            "n": len(sub),
            "phi_median": round(statistics.median(phi), 3),
            "phi_q1": round(_q(phi, 0.25), 3), "phi_q3": round(_q(phi, 0.75), 3),
            "eta_median": round(statistics.median(eta), 3),
            "eta_q1": round(_q(eta, 0.25), 3), "eta_q3": round(_q(eta, 0.75), 3),
            "frac_with_necrotic_core": round(core_present, 3),
            "model_phi_abs_err": round(abs(MODEL_RIM_BOUNDARY - statistics.median(phi)), 3),
            "model_eta_abs_err": round(abs(MODEL_CORE_BOUNDARY - statistics.median(eta)), 3),
        })
    return out


def valid_size_range(summary):
    """Smallest R (µm) above which the fixed-threshold model matches the data:
    both boundary errors <= MATCH_TOL and a necrotic core present in most spheroids."""
    matching = [
        b for b in summary
        if b["model_phi_abs_err"] <= MATCH_TOL
        and b["model_eta_abs_err"] <= MATCH_TOL
        and b["frac_with_necrotic_core"] >= CORE_PRESENT_FRAC
    ]
    if not matching:
        return None
    return min(b["r_lo_um"] for b in matching)


def run(args):
    rows = fetch_rows(raw_csv=args.raw_csv)
    summary = summarize(rows)
    valid_lo = valid_size_range(summary)
    result = {
        "source": "Browning et al. 2021 eLife (PMID 34842141), Data/AllConfocalData.csv",
        "data_repo": "github.com/ap-browning/Spheroids",
        "n_spheroids": len(rows),
        "cell_lines": sorted({r["cell_line"] for r in rows}),
        "r_um_min": round(min(r["R"] for r in rows), 1),
        "r_um_max": round(max(r["R"] for r in rows), 1),
        "model_rim_boundary": round(MODEL_RIM_BOUNDARY, 4),
        "model_core_boundary": round(MODEL_CORE_BOUNDARY, 4),
        "match_tolerance": MATCH_TOL,
        "size_bins": summary,
        "valid_radius_um_min": valid_lo,
        "valid_diameter_um_min": (None if valid_lo is None else 2 * valid_lo),
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_report(result)
    print(f"n={len(rows)} spheroids, R {result['r_um_min']}-{result['r_um_max']} µm")
    for b in summary:
        hi = b["r_hi_um"] if b["r_hi_um"] else "+"
        print(f"  R {b['r_lo_um']}-{hi}: n={b['n']} phi={b['phi_median']} eta={b['eta_median']} "
              f"core_present={b['frac_with_necrotic_core']}")
    print(f"valid size range: R >= {valid_lo} µm (diameter >= {result['valid_diameter_um_min']} µm)")
    print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)} + {OUT_MD.relative_to(REPO_ROOT)}")
    return result


def write_report(r):
    lines = [
        "# Spheroid zone-geometry validation vs Browning 2021 (#333)",
        "",
        "Generated by `scripts/validate_spheroid_structure.py` (offline; the raw data is",
        "not committed, the derived summary is). Validates the `spheroid` module's FIXED",
        "radial zone thresholds against size-resolved confocal structure data.",
        "",
        "## Source",
        "",
        f"Browning et al. 2021 eLife (PMID 34842141), `Data/AllConfocalData.csv` from",
        f"`github.com/ap-browning/Spheroids`: **{r['n_spheroids']} spheroids**, cell lines",
        f"{', '.join(r['cell_lines'])}, outer radius {r['r_um_min']}-{r['r_um_max']} µm. Each",
        "spheroid reports the outer radius R and two internal boundaries as fractions of R:",
        "`phi` (proliferating-rim inner edge) and `eta` (necrotic-core outer edge).",
        "",
        "## The model's fixed thresholds (size-independent)",
        "",
        f"- rim / proliferating-boundary analog: **{r['model_rim_boundary']}** of radius",
        f"  (`glycolytic_frac=0.73`, 0.73^(1/3)).",
        f"- core / necrotic-boundary analog: **{r['model_core_boundary']}** of radius",
        f"  (`oxphos_frac=0.39`, 0.39^(1/3)).",
        "",
        "## Measured boundaries vs size (the finite-size effect)",
        "",
        "| R bin (µm) | n | phi median (model 0.90) | eta median (model 0.73) | frac w/ necrotic core |",
        "|---|---:|---:|---:|---:|",
    ]
    for b in r["size_bins"]:
        hi = b["r_hi_um"] if b["r_hi_um"] else "+"
        lines.append(
            f"| {b['r_lo_um']}-{hi} | {b['n']} | {b['phi_median']} | {b['eta_median']} | "
            f"{b['frac_with_necrotic_core']} |"
        )
    vd = r["valid_diameter_um_min"]
    lines += [
        "",
        "## Result: valid size range",
        "",
        f"The fixed-threshold model matches the data (both boundary errors <= "
        f"{r['match_tolerance']} of radius AND a necrotic core present in most spheroids)",
        f"only for **R >= {r['valid_radius_um_min']} µm (diameter >= {vd} µm)**.",
        "",
        "Below that size the model is wrong in a specific, quantified way:",
        "- Small spheroids (R < 300 µm) have essentially **no necrotic core** (eta ~ 0,",
        "  necrotic core present in a small minority), but the model places a persister/",
        "  necrotic-like core at 0.73 of the radius regardless. The model therefore",
        "  **over-predicts the resistant core** for small spheroids.",
        "- The proliferating rim is thick in small spheroids (low phi: most of the volume",
        "  proliferates) and thins toward the model's 0.90 only as the spheroid grows.",
        "- The necrotic core emerges near R ~ 300 µm (diameter ~ 600 µm), consistent with",
        "  an O2 diffusion length of ~150-200 µm.",
        "",
        "## Scope of what this constrains",
        "",
        "- This validates the zone **geometry** (where the boundaries sit vs size). It does",
        "  NOT calibrate the per-zone biochem gradient **strengths** (MUFA/GSH/iron), which",
        "  need spatially-resolved metabolomics, not structural boundaries, so those remain",
        "  uncalibrated placeholders (unchanged).",
        "- It validates STRUCTURE; size-resolved depth-dependent KILL data (the other half",
        "  of #333) is not in this dataset and is the remaining piece.",
        "- The model's fixed thresholds are the large-spheroid LIMITING structure, so the",
        "  practical guidance is: trust the spheroid layer's zones only for spheroids at or",
        f"  above ~{vd} µm diameter; treat smaller-spheroid runs as over-coring.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-csv", type=Path, default=None, help="use a local AllConfocalData.csv instead of fetching")
    args = ap.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
