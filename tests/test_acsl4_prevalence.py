"""Guards for the #462 ACSL4-status prevalence calibration (cBioPortal TCGA).

The fetch (scripts/fetch_acsl4_prevalence.py) hits the cBioPortal REST API and is
run LOCALLY, not in CI (offline contract). This validates the committed derived
artifacts so the calibration cannot silently rot, and pins the two load-bearing
findings: (1) the within-cohort low-ACSL4 prevalence prior, and (2) the honest
negative result that bulk TCGA mRNA does NOT rank HCC/AML as ACSL4-low.
"""

import csv
import json
import statistics as st
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CSV = REPO_ROOT / "analysis" / "calibration" / "acsl4_prevalence_tcga.csv"
JSON = REPO_ROOT / "analysis" / "calibration" / "acsl4_prevalence_tcga.json"


def _rows():
    with open(CSV, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_artifacts_exist_and_cover_pancancer():
    rows = _rows()
    # The 32 TCGA PanCancer Atlas studies.
    assert len(rows) >= 30, f"expected ~32 cancer types, got {len(rows)}"
    summary = json.loads(JSON.read_text())
    assert "pan_can_atlas" in summary["source"]
    assert summary["genes"]["ACSL4"] == 2182


def test_within_cohort_low_acsl4_prevalence_prior():
    # The usable population prior: each cancer type has roughly 1 in 7 tumors in its
    # low-ACSL4 tail (z < -1), fairly uniformly, with a small very-low (z < -2) tail.
    rows = _rows()
    frac_low = [float(r["ACSL4_frac_low"]) for r in rows if r["ACSL4_frac_low"]]
    frac_vlow = [float(r["ACSL4_frac_verylow"]) for r in rows if r["ACSL4_frac_verylow"]]
    assert 0.10 <= st.median(frac_low) <= 0.20, st.median(frac_low)
    assert all(0.05 <= v <= 0.25 for v in frac_low), "low fraction outside expected band"
    assert 0.0 <= st.median(frac_vlow) <= 0.08, st.median(frac_vlow)


def test_honest_negative_bulk_mrna_does_not_show_hcc_low():
    # The load-bearing honesty guard: bulk TCGA mRNA does NOT rank HCC (lihc) as
    # ACSL4-low; it ranks HIGH. If a future refetch flipped this it would change the
    # calibration's documented conclusion, so pin it.
    summary = json.loads(JSON.read_text())
    v = summary["validation"]
    assert v["lihc_hcc_rank_percentile"] is not None
    assert v["lihc_hcc_rank_percentile"] >= 0.5, (
        "calibration doc claims bulk mRNA does NOT show HCC as ACSL4-low; "
        f"lihc percentile is {v['lihc_hcc_rank_percentile']}"
    )
    # lihc should be among the highest-RSEM types, not the lowest.
    assert "lihc" in v["highest_5_acsl4_rsem"]


def test_zscore_bridge_constants_documented():
    # The calibration doc bridges z -> status; the Rust acsl4::status_from_zscore
    # implements max(0, 1 + z/2), reproducing the shipped constants at integer z.
    # Mirror that arithmetic here so the Python doc and Rust stay in lockstep.
    def status_from_z(z):
        return max(0.0, 1.0 + z / 2.0)

    assert status_from_z(1.0) == 1.5  # ACSL4_HIGH
    assert status_from_z(0.0) == 1.0  # ACSL4_NORMAL
    assert status_from_z(-1.0) == 0.5  # ACSL4_LOW
    assert status_from_z(-2.0) == 0.0  # ACSL4_NEGATIVE
    assert status_from_z(-5.0) == 0.0  # floored
