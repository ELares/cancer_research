#!/usr/bin/env python3
"""Fetch the CTRPv2 ferroptosis-inducer dose-response calibration target (#330).

WHY
---
The ferroptosis-core kill switch has been honest-but-uncalibrated scaffolding
(see simulations/calibration/CALIBRATION_STATUS.md). #330 calibrates its kill
rates against PUBLIC drug-response data with held-out validation. GDSC2 does not
screen the canonical ferroptosis inducers (only Cisplatin/Elesclomol/Sorafenib),
so the right source is CTRPv2 (Cancer Therapeutics Response Portal v2;
Seashore-Ludlow 2015 Cancer Discov, Rees 2016 Nat Chem Biol), which screened
erastin, ML162, ML210 and related probes across ~860 cell lines.

The original NCI CTD2 portal URL is dead (301 -> studycatalog.cancer.gov). This
script pulls the reprocessed CTRPv2 that DepMap redistributes ("Harmonized CTD^2"
release) via the DepMap download API, which is login-free.

WHAT IT DOES (provenance tool, NOT run in CI)
---------------------------------------------
1. GET the DepMap download catalog, resolve the (time-limited, signed) URL for
   `CTRPResponseCurves.csv`, download it to a raw cache dir (NOT committed).
2. Verify the download MD5 against the catalog's `md5_hash`.
3. Filter to the ferroptosis inducers, write the small derived target
   `analysis/calibration/ctrpv2_ferroptosis_curves.csv` (committed) plus a
   per-compound summary `ctrpv2_ferroptosis_summary.json`.

The offline/reproducible CI contract is preserved: CI never downloads anything;
downstream code and tests read the committed derived CSV. This script is the
documented, re-runnable way to regenerate that derivative from the public source.

The committed CSV stores the per-cell-line 4-parameter logistic fit
(EC50, lower/upper asymptote, Hill slope) over the screened dose range. Use
`predicted_viability` to reconstruct viability(dose); note the fitted slope is
NEGATIVE and the curve is DECREASING, so the exponent is `-slope` (see the
function + its sanity assert).
"""

import argparse
import csv
import hashlib
import json
import statistics
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "analysis" / "calibration"
CURVES_CSV = OUT_DIR / "ctrpv2_ferroptosis_curves.csv"
SUMMARY_JSON = OUT_DIR / "ctrpv2_ferroptosis_summary.json"

DEPMAP_CATALOG = "https://depmap.org/portal/api/download/files"
TARGET_FILENAME = "CTRPResponseCurves.csv"
RELEASE_SUBSTR = "CTD"  # matches the "Harmonized CTD^2 ..." release

# Ferroptosis inducers screened in CTRPv2. GPX4 inhibitors (ML162, ML210) are the
# direct analog of the model's RSL3/GPX4i kill switch; erastin is the system-xc-
# (cystine-import) inhibitor; CIL55/CIL56 are additional ferroptosis probes.
FERROPTOSIS_COMPOUNDS = ("ML162", "ML210", "ERASTIN", "CIL56", "CIL55")

# Columns kept in the committed derived CSV.
KEEP_COLS = (
    "CompoundName", "ModelID", "EC50", "LowerAsymptote", "UpperAsymptote",
    "Slope", "MinimumDose", "MaximumDose", "DoseUnit",
)


def predicted_viability(dose: float, lower: float, upper: float, ec50: float, slope: float) -> float:
    """CTRPv2 4-parameter logistic viability at `dose`.

    The fitted `slope` is NEGATIVE for a cytotoxic curve, so the decreasing form
    uses exponent `-slope`: viability is `upper` (~1) at low dose and `lower`
    (the residual-viability / kill-ceiling complement) at saturating dose.

    Evaluated in log space and clamped so non-responder fits (EC50 far above the
    screened range, steep slope) cannot overflow float range.
    """
    import math

    log_term = (-slope) * math.log(dose / ec50)  # dose, ec50 > 0
    if log_term > 700:  # term -> +inf, curve at its low-dose... high-dose asymptote
        return lower
    if log_term < -700:
        return upper
    return lower + (upper - lower) / (1.0 + math.exp(log_term))


def auc_fraction(lower: float, upper: float, ec50: float, slope: float, dmin: float, dmax: float, n: int = 256) -> float:
    """Mean predicted viability across the screened log-dose range (area under the
    viability curve, normalized to [0,1]). Lower means more potent."""
    import math

    lo, hi = math.log10(dmin), math.log10(dmax)
    xs = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    vals = [predicted_viability(10 ** x, lower, upper, ec50, slope) for x in xs]
    return sum(vals) / len(vals)


def _quantile(sorted_vals, q):
    if not sorted_vals:
        return None
    idx = min(len(sorted_vals) - 1, int(q * len(sorted_vals)))
    return sorted_vals[idx]


def resolve_signed_url(catalog_text: str) -> "tuple[str, str, str]":
    """Return (url, md5_hash, release) for the target file from the catalog CSV."""
    reader = csv.DictReader(catalog_text.splitlines())
    for row in reader:
        if (row.get("filename") or "").strip() == TARGET_FILENAME and RELEASE_SUBSTR in (row.get("release") or ""):
            return row.get("url"), (row.get("md5_hash") or "").strip(), (row.get("release") or "").strip()
    raise SystemExit(f"could not find {TARGET_FILENAME} (release ~{RELEASE_SUBSTR}) in the DepMap catalog")


def filter_curves(raw_csv_path: Path, compounds=FERROPTOSIS_COMPOUNDS) -> "list[dict]":
    """Read the full curves CSV, keep only ferroptosis-inducer rows with finite fits."""
    wanted = {c.upper() for c in compounds}
    out = []
    with open(raw_csv_path, newline="") as f:
        for row in csv.DictReader(f):
            name = (row.get("CompoundName") or "").strip().upper()
            if name not in wanted:
                continue
            try:
                rec = {k: row[k] for k in KEEP_COLS}
                for num in ("EC50", "LowerAsymptote", "UpperAsymptote", "Slope", "MinimumDose", "MaximumDose"):
                    float(rec[num])
            except (KeyError, ValueError, TypeError):
                continue
            rec["CompoundName"] = name
            out.append(rec)
    out.sort(key=lambda r: (r["CompoundName"], r["ModelID"]))
    return out


def summarize(rows: "list[dict]") -> dict:
    """Per-compound empirical summary: n, EC50 quartiles, median residual viability, median AUC."""
    by = {}
    for r in rows:
        by.setdefault(r["CompoundName"], []).append(r)
    summary = {}
    for cmpd, recs in sorted(by.items()):
        ec50 = sorted(float(r["EC50"]) for r in recs)
        lowers = sorted(float(r["LowerAsymptote"]) for r in recs)
        aucs = sorted(
            auc_fraction(
                float(r["LowerAsymptote"]), float(r["UpperAsymptote"]), float(r["EC50"]),
                float(r["Slope"]), float(r["MinimumDose"]), float(r["MaximumDose"]),
            )
            for r in recs
        )
        summary[cmpd] = {
            "n_cell_lines": len(recs),
            "ec50_um_median": round(statistics.median(ec50), 4),
            "ec50_um_q1": round(_quantile(ec50, 0.25), 4),
            "ec50_um_q3": round(_quantile(ec50, 0.75), 4),
            "residual_viability_median": round(statistics.median(lowers), 4),
            "kill_ceiling_median": round(1.0 - statistics.median(lowers), 4),
            "auc_viability_median": round(statistics.median(aucs), 4),
            "dose_um_min": min(float(r["MinimumDose"]) for r in recs),
            "dose_um_max": max(float(r["MaximumDose"]) for r in recs),
        }
    return summary


def write_outputs(rows: "list[dict]", release: str, md5_hash: str) -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(CURVES_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(KEEP_COLS))
        w.writeheader()
        w.writerows(rows)
    summary = summarize(rows)
    meta = {
        "source": "DepMap-redistributed CTRPv2 (Cancer Therapeutics Response Portal v2)",
        "depmap_release": release,
        "depmap_file": TARGET_FILENAME,
        "depmap_file_md5": md5_hash,
        "primary_refs": ["Seashore-Ludlow 2015 Cancer Discov PMID 26181016", "Rees 2016 Nat Chem Biol PMID 26656090"],
        "n_curves": len(rows),
        "compounds": summary,
    }
    SUMMARY_JSON.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", type=Path, default=Path("/tmp"),
                    help="cache dir for the raw 45MB download (NOT committed)")
    ap.add_argument("--raw-csv", type=Path, default=None,
                    help="use an already-downloaded CTRPResponseCurves.csv instead of fetching")
    args = ap.parse_args()

    if args.raw_csv:
        raw = args.raw_csv
        release, md5_hash = "(local file)", ""
    else:
        print(f"fetching catalog {DEPMAP_CATALOG}")
        with urllib.request.urlopen(DEPMAP_CATALOG, timeout=120) as r:
            catalog = r.read().decode("utf-8", "replace")
        url, md5_hash, release = resolve_signed_url(catalog)
        print(f"release: {release}  md5: {md5_hash}")
        args.raw_dir.mkdir(parents=True, exist_ok=True)
        raw = args.raw_dir / TARGET_FILENAME
        print(f"downloading -> {raw}")
        urllib.request.urlretrieve(url, raw)
        got = hashlib.md5(raw.read_bytes()).hexdigest()
        if md5_hash and got != md5_hash:
            raise SystemExit(f"MD5 mismatch: catalog {md5_hash} != download {got}")
        print(f"md5 verified: {got}")

    rows = filter_curves(raw)
    meta = write_outputs(rows, release, md5_hash)
    print(f"wrote {CURVES_CSV.relative_to(REPO_ROOT)} ({len(rows)} curves)")
    print(f"wrote {SUMMARY_JSON.relative_to(REPO_ROOT)}")
    for cmpd, s in meta["compounds"].items():
        print(f"  {cmpd:<8} n={s['n_cell_lines']:<4} medEC50={s['ec50_um_median']} uM  "
              f"kill_ceiling={s['kill_ceiling_median']}  AUC={s['auc_viability_median']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
