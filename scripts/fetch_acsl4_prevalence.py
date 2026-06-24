#!/usr/bin/env python3
"""Fetch per-cancer-type ACSL4 (+ GPX4, SLC7A11) mRNA distributions from cBioPortal
TCGA PanCancer Atlas, and derive the ACSL4-low/negative prevalence that calibrates
the ACSL4-status biomarker layer (#444, ferroptosis-core acsl4 module).

The #444 layer maps a per-tumor ACSL4 expression *status* (1.0 = wild-type baseline,
< 1 = ferroptosis-refractory via a collapsed PUFA substrate) to a PUFA-incorporation
boost, but the per-cancer-type prevalence of ACSL4-low/negative tumors was flagged
DATA-GATED. This script closes that gap with the only login-free public source: the
cBioPortal REST API over the 32 TCGA PanCancer Atlas studies.

Two complementary readouts (both committed, with honest caveats in the report):

  1. WITHIN-cohort low-ACSL4 prevalence (z-score tails): the fraction of tumors per
     cancer type with ACSL4 mRNA z < -1 ("low") and z < -2 ("very low", a
     refractory-leaning proxy). The z-scores are computed WITHIN each study, so this
     is the within-cohort stratification prior (about a normal lower tail by
     construction; the cross-type spread is small here and that is the point).

  2. CROSS-cohort relative level (raw RSEM medians): the per-cancer-type median raw
     ACSL4 expression, used ONLY to RANK which cancer TYPES sit constitutively low.
     Cross-study RSEM medians carry batch effects, so this is read qualitatively (a
     ranking), and validated against the Doll-2017 literature claim that HCC (lihc)
     and AML (laml) are constitutively ACSL4-low.

OFFLINE CONTRACT (mirrors scripts/fetch_calibration_data.py, #330): this script
hits the network and is run LOCALLY; CI never runs it. The committed derived
artifacts (analysis/calibration/acsl4_prevalence_tcga.csv + .json) are the
reproducible result that the test suite and the calibration doc read.

Usage:
    python3 scripts/fetch_acsl4_prevalence.py            # fetch + write artifacts
    python3 scripts/fetch_acsl4_prevalence.py --check     # print summary, write nothing
"""
import argparse
import json
import statistics as st
import sys
import urllib.request
from pathlib import Path

CBIOPORTAL = "https://www.cbioportal.org/api"
REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "analysis" / "calibration"

# Hugo symbol -> Entrez Gene ID for the ferroptosis lipid/defense axis genes.
GENES = {"ACSL4": 2182, "GPX4": 2879, "SLC7A11": 23657}

# z-score thresholds for the within-cohort low-expression tails.
LOW_Z = -1.0
VERYLOW_Z = -2.0


def _get(path: str):
    req = urllib.request.Request(
        CBIOPORTAL + path, headers={"accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def _post(path: str, body: dict):
    data = json.dumps(body).encode()
    req = urllib.request.Request(
        CBIOPORTAL + path,
        data=data,
        headers={"accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def pan_can_atlas_studies() -> list[dict]:
    studies = _get("/studies?pageSize=2000")
    pan = [s for s in studies if "pan_can_atlas" in s.get("studyId", "")]
    return sorted(pan, key=lambda s: s["studyId"])


def _profile_ids(study_id: str) -> tuple[str | None, str | None]:
    """Return (zscore_profile_id, raw_rsem_profile_id) for a study, or None each."""
    profs = _get(f"/studies/{study_id}/molecular-profiles")
    z = raw = None
    for p in profs:
        if p.get("molecularAlterationType") != "MRNA_EXPRESSION":
            continue
        pid = p["molecularProfileId"]
        if p.get("datatype") == "Z-SCORE" and pid.endswith("_all_sample_Zscores"):
            z = pid
        elif p.get("datatype") == "CONTINUOUS" and pid.endswith("_rna_seq_v2_mrna"):
            raw = pid
    return z, raw


def _fetch_values(profile_id: str, study_id: str, entrez: int) -> list[float]:
    body = {"entrezGeneIds": [entrez], "sampleListId": f"{study_id}_all"}
    rows = _post(f"/molecular-profiles/{profile_id}/molecular-data/fetch", body)
    return [r["value"] for r in rows if r.get("value") is not None]


def fetch_all() -> list[dict]:
    out = []
    studies = pan_can_atlas_studies()
    print(f"Fetching {len(studies)} TCGA PanCancer Atlas studies from cBioPortal...")
    for s in studies:
        sid = s["studyId"]
        ctype = sid.replace("_tcga_pan_can_atlas_2018", "")
        z_pid, raw_pid = _profile_ids(sid)
        row: dict = {
            "cancer_type": ctype,
            "study_id": sid,
            "cancer_type_id": s.get("cancerTypeId"),
        }
        for gene, entrez in GENES.items():
            zvals = _fetch_values(z_pid, sid, entrez) if z_pid else []
            rawvals = _fetch_values(raw_pid, sid, entrez) if raw_pid else []
            n = len(zvals)
            row[f"{gene}_n"] = n
            row[f"{gene}_median_z"] = round(st.median(zvals), 4) if zvals else None
            row[f"{gene}_frac_low"] = (
                round(sum(1 for v in zvals if v < LOW_Z) / n, 4) if n else None
            )
            row[f"{gene}_frac_verylow"] = (
                round(sum(1 for v in zvals if v < VERYLOW_Z) / n, 4) if n else None
            )
            row[f"{gene}_median_rsem"] = (
                round(st.median(rawvals), 2) if rawvals else None
            )
        out.append(row)
        a = row.get("ACSL4_frac_low")
        print(
            f"  {ctype:8s} n={row.get('ACSL4_n'):4d}  ACSL4 frac_low={a}  "
            f"median_rsem={row.get('ACSL4_median_rsem')}"
        )
    return out


def _validation(rows: list[dict]) -> dict:
    """Validate the Doll-2017 direction: HCC (lihc) and AML (laml) should rank among
    the lower cancer types by raw ACSL4 RSEM. Reported qualitatively (batch-effect
    caveat)."""
    ranked = sorted(
        [r for r in rows if r.get("ACSL4_median_rsem") is not None],
        key=lambda r: r["ACSL4_median_rsem"],
    )
    order = [r["cancer_type"] for r in ranked]
    n = len(order)

    def pct_rank(ct):
        return round(order.index(ct) / (n - 1), 3) if ct in order and n > 1 else None

    return {
        "n_cancer_types_ranked": n,
        "lowest_5_acsl4_rsem": order[:5],
        "highest_5_acsl4_rsem": order[-5:],
        "lihc_hcc_rank_percentile": pct_rank("lihc"),
        "laml_aml_rank_percentile": pct_rank("laml"),
        "note": (
            "Lower percentile = lower raw ACSL4 RSEM (more refractory-leaning). "
            "Cross-study RSEM medians carry batch effects, so this ranking is "
            "qualitative; the Doll-2017 claim is that HCC/AML are constitutively "
            "ACSL4-low."
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="Fetch + print, write nothing")
    args = ap.parse_args()

    rows = fetch_all()
    validation = _validation(rows)
    print("\n=== Cross-type ACSL4 raw-RSEM ranking (Doll-2017 validation) ===")
    print("  lowest 5:", validation["lowest_5_acsl4_rsem"])
    print("  highest 5:", validation["highest_5_acsl4_rsem"])
    print(f"  lihc (HCC) percentile: {validation['lihc_hcc_rank_percentile']}")
    print(f"  laml (AML) percentile: {validation['laml_aml_rank_percentile']}")

    if args.check:
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # CSV
    cols = ["cancer_type", "study_id", "cancer_type_id"]
    for g in GENES:
        cols += [
            f"{g}_n",
            f"{g}_median_z",
            f"{g}_frac_low",
            f"{g}_frac_verylow",
            f"{g}_median_rsem",
        ]
    csv_path = OUT_DIR / "acsl4_prevalence_tcga.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(c, "")) for c in cols) + "\n")
    # JSON summary
    acsl4_low = [r["ACSL4_frac_low"] for r in rows if r.get("ACSL4_frac_low") is not None]
    summary = {
        "source": "cBioPortal REST API, TCGA PanCancer Atlas (_tcga_pan_can_atlas_2018)",
        "endpoint": CBIOPORTAL,
        "genes": GENES,
        "low_z_threshold": LOW_Z,
        "verylow_z_threshold": VERYLOW_Z,
        "n_cancer_types": len(rows),
        "acsl4_frac_low_median_across_types": round(st.median(acsl4_low), 4)
        if acsl4_low
        else None,
        "acsl4_frac_low_min_max": [min(acsl4_low), max(acsl4_low)] if acsl4_low else None,
        "validation": validation,
    }
    json_path = OUT_DIR / "acsl4_prevalence_tcga.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nWrote {csv_path} ({len(rows)} cancer types) and {json_path}")


if __name__ == "__main__":
    main()
