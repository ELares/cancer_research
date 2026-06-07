"""Tests for the CTRPv2 ferroptosis calibration-target tooling (#330).

Covers the dose-response math (overflow-safe 4-parameter logistic, AUC), the
derive/summary logic on synthetic input, the catalog-URL resolver, and structural
self-consistency of the committed artifacts (the derived CSV + summary JSON that
downstream calibration reads, since CI never re-downloads).
"""

import csv
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import fetch_calibration_data as fcd  # noqa: E402

CURVES_CSV = REPO_ROOT / "analysis" / "calibration" / "ctrpv2_ferroptosis_curves.csv"
SUMMARY_JSON = REPO_ROOT / "analysis" / "calibration" / "ctrpv2_ferroptosis_summary.json"


# --------------------------------------------------------------------------
# Dose-response math
# --------------------------------------------------------------------------


def test_viability_is_decreasing_with_dose():
    lo, hi, ec50, slope = 0.05, 1.0, 0.5, -6.0
    v_low = fcd.predicted_viability(1e-3, lo, hi, ec50, slope)
    v_mid = fcd.predicted_viability(ec50, lo, hi, ec50, slope)
    v_high = fcd.predicted_viability(33.0, lo, hi, ec50, slope)
    assert v_low > v_mid > v_high  # cytotoxic curve falls with dose
    assert v_low == pytest.approx(hi, abs=0.02)
    assert v_high == pytest.approx(lo, abs=0.02)
    assert v_mid == pytest.approx((lo + hi) / 2, abs=1e-6)  # midpoint at EC50


def test_viability_overflow_safe_for_nonresponders():
    # Huge EC50 + steep slope (a real non-responder fit) must not overflow.
    v = fcd.predicted_viability(10.0, 0.5, 1.01, 4.75e8, -0.45)
    assert 0.5 <= v <= 1.01
    # extreme steep slope at high relative dose -> clamps to an asymptote, no crash
    v2 = fcd.predicted_viability(1e6, 0.02, 1.0, 1e-3, -50.0)
    assert v2 == pytest.approx(0.02, abs=1e-9)


def test_auc_bounds_and_potency_ordering():
    # More potent (lower EC50) => lower AUC over the same dose range.
    potent = fcd.auc_fraction(0.02, 1.0, 0.1, -6.0, 1e-3, 33.0)
    weak = fcd.auc_fraction(0.02, 1.0, 10.0, -6.0, 1e-3, 33.0)
    assert 0.0 <= potent <= 1.0 and 0.0 <= weak <= 1.0
    assert potent < weak


# --------------------------------------------------------------------------
# Derive / summarize / resolve (synthetic)
# --------------------------------------------------------------------------


def _write_raw(tmp_path, rows):
    p = tmp_path / "raw.csv"
    cols = ["ModelID", "SampleID", "CompoundName", "CompoundID", "EC50",
            "LowerAsymptote", "UpperAsymptote", "Slope", "MinimumDose", "MaximumDose", "DoseUnit"]
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return p


def test_filter_curves_keeps_only_ferroptosis_and_finite(tmp_path):
    rows = [
        dict(ModelID="ACH-1", SampleID="s", CompoundName="ML210", CompoundID="d", EC50="0.5",
             LowerAsymptote="0.1", UpperAsymptote="1.0", Slope="-6", MinimumDose="0.001", MaximumDose="33", DoseUnit="uM"),
        dict(ModelID="ACH-2", SampleID="s", CompoundName="paclitaxel", CompoundID="d", EC50="0.5",
             LowerAsymptote="0.1", UpperAsymptote="1.0", Slope="-6", MinimumDose="0.001", MaximumDose="33", DoseUnit="uM"),
        dict(ModelID="ACH-3", SampleID="s", CompoundName="ML162", CompoundID="d", EC50="not_a_number",
             LowerAsymptote="0.1", UpperAsymptote="1.0", Slope="-6", MinimumDose="0.001", MaximumDose="33", DoseUnit="uM"),
        dict(ModelID="ACH-4", SampleID="s", CompoundName="erastin", CompoundID="d", EC50="4.0",
             LowerAsymptote="0.2", UpperAsymptote="1.0", Slope="-3", MinimumDose="0.002", MaximumDose="66", DoseUnit="uM"),
    ]
    out = fcd.filter_curves(_write_raw(tmp_path, rows))
    names = {r["CompoundName"] for r in out}
    assert names == {"ML210", "ERASTIN"}  # paclitaxel dropped, bad-EC50 ML162 dropped, names upper-cased
    assert all(set(fcd.KEEP_COLS).issubset(r.keys()) for r in out)


def test_summarize_matches_hand_computation():
    rows = [
        dict(CompoundName="ML210", ModelID="a", EC50="0.4", LowerAsymptote="0.1", UpperAsymptote="1.0", Slope="-6", MinimumDose="0.001", MaximumDose="33", DoseUnit="uM"),
        dict(CompoundName="ML210", ModelID="b", EC50="0.6", LowerAsymptote="0.2", UpperAsymptote="1.0", Slope="-6", MinimumDose="0.001", MaximumDose="33", DoseUnit="uM"),
        dict(CompoundName="ML210", ModelID="c", EC50="0.5", LowerAsymptote="0.15", UpperAsymptote="1.0", Slope="-6", MinimumDose="0.001", MaximumDose="33", DoseUnit="uM"),
    ]
    s = fcd.summarize(rows)["ML210"]
    assert s["n_cell_lines"] == 3
    assert s["ec50_um_median"] == pytest.approx(0.5)
    assert s["residual_viability_median"] == pytest.approx(0.15)
    assert s["kill_ceiling_median"] == pytest.approx(0.85)


def test_resolve_signed_url_picks_target_row():
    catalog = (
        "release,release_date,filename,url,md5_hash\n"
        "Other Release,2024,SomethingElse.csv,http://x/other,aaa\n"
        "Harmonized CTD^2 25Q2,2025,CTRPResponseCurves.csv,http://x/signed,bbb\n"
    )
    url, md5, release = fcd.resolve_signed_url(catalog)
    assert url == "http://x/signed" and md5 == "bbb" and "CTD" in release


def test_resolve_signed_url_raises_when_missing():
    with pytest.raises(SystemExit):
        fcd.resolve_signed_url("release,filename,url,md5_hash\nX,Nope.csv,u,h\n")


# --------------------------------------------------------------------------
# Committed-artifact self-consistency
# --------------------------------------------------------------------------


def test_committed_curves_schema_and_compounds():
    assert CURVES_CSV.exists()
    with open(CURVES_CSV, newline="") as f:
        reader = csv.DictReader(f)
        assert list(reader.fieldnames) == list(fcd.KEEP_COLS)
        names = set()
        n = 0
        for row in reader:
            names.add(row["CompoundName"])
            float(row["EC50"]); float(row["LowerAsymptote"]); float(row["Slope"])  # finite
            n += 1
    assert names.issubset(set(fcd.FERROPTOSIS_COMPOUNDS))
    assert n > 2000  # ~3021 curves


def test_committed_summary_consistent_with_csv():
    summary = json.loads(SUMMARY_JSON.read_text())
    assert summary["depmap_file"] == "CTRPResponseCurves.csv"
    assert summary["depmap_file_md5"]  # provenance recorded
    # n_curves matches the CSV row count
    with open(CURVES_CSV, newline="") as f:
        csv_rows = list(csv.DictReader(f))
    assert summary["n_curves"] == len(csv_rows)
    # per-compound n_cell_lines sums to the total
    assert sum(c["n_cell_lines"] for c in summary["compounds"].values()) == len(csv_rows)
    # the GPX4 inhibitors are present with plausible sub-micromolar median EC50
    for gpx4i in ("ML162", "ML210"):
        assert summary["compounds"][gpx4i]["ec50_um_median"] < 2.0
        assert 0.0 < summary["compounds"][gpx4i]["kill_ceiling_median"] <= 1.0
