"""Assay-support guards for fitted CTRPv2 calibration targets."""

import csv
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from ctrp_dose_support import select_supported_cohort  # noqa: E402


def curve(model, minimum="0.1", maximum="30", **overrides):
    return {"ModelID": model, "SampleID": f"sample-{model}",
            "CompoundName": "ERASTIN", "MinimumDose": minimum,
            "MaximumDose": maximum, "DoseUnit": "uM", **overrides}


def test_common_cohort_preserves_boundaries_objects_and_denominator():
    rows = [curve("low-only", maximum="1"), curve("both"),
            curve("high-only", minimum="1"), curve("also-both")]
    retained, meta = select_supported_cohort(iter(rows), (30, 0.1, 1))
    assert len(retained) == 2
    assert retained[0] is rows[1] and retained[1] is rows[3]
    assert meta["n_original"] == 4 and meta["n_retained"] == 2
    assert meta["n_excluded"] == 2
    # Support varies by dose in the source; all targets use the same two rows.
    assert meta["per_dose_support"] == [
        {"dose_um": 30.0, "n_supported": 3},
        {"dose_um": 0.1, "n_supported": 3},
        {"dose_um": 1.0, "n_supported": 4},
    ]
    assert meta["excluded_rows"][0]["ModelID"] == "low-only"
    assert meta["excluded_rows"][0]["reasons"] == ["requested_dose_above_maximum"]
    assert meta["excluded_rows"][1]["reasons"] == ["requested_dose_below_minimum"]
    assert meta == json.loads(json.dumps(meta))
    assert meta == select_supported_cohort(rows, (30, 0.1, 1))[1]


def test_individually_supported_doses_do_not_make_a_common_cohort():
    rows = [curve("low", maximum="1"), curve("high", minimum="10")]
    with pytest.raises(ValueError, match="no curve supports the entire"):
        select_supported_cohort(rows, (0.1, 30))


def test_zero_support_dose_reports_its_count():
    with pytest.raises(ValueError, match="100 uM: 0"):
        select_supported_cohort([curve("a")], (0.1, 100))


@pytest.mark.parametrize("doses", [(), (0,), (-1,), (float("nan"),),
                                  (float("inf"),), ("bad",), (True,)])
def test_invalid_grid_rejected(doses):
    with pytest.raises(ValueError, match="dose"):
        select_supported_cohort([curve("a")], doses)


@pytest.mark.parametrize("field,value", [
    ("MinimumDose", "nan"), ("MaximumDose", "inf"),
    ("MinimumDose", "-1"), ("MaximumDose", "0"),
    ("MinimumDose", None), ("MaximumDose", "bad"),
    ("MinimumDose", True), ("MaximumDose", "0.01"),
    ("DoseUnit", "nM"), ("DoseUnit", None),
])
def test_invalid_metadata_rejected_even_if_another_row_is_supported(field, value):
    # An invalid observation must not silently disappear from the denominator.
    rows = [curve("good"), curve("bad", **{field: value})]
    with pytest.raises(ValueError, match=f"row 1:.*{field}|row 1: MinimumDose exceeds"):
        select_supported_cohort(rows, (0.1, 30))


def test_missing_range_rejected():
    row = curve("missing")
    del row["MaximumDose"]
    with pytest.raises(ValueError, match="MaximumDose"):
        select_supported_cohort([row], (0.1, 30))


def test_empty_input_and_mixed_compounds_rejected():
    with pytest.raises(ValueError, match="input curve cohort is empty"):
        select_supported_cohort([], (1,))
    with pytest.raises(ValueError, match="one compound at a time"):
        select_supported_cohort([curve("a"), curve("b", CompoundName="ML162")], (1,))


def test_committed_erastin_cannot_support_100_um():
    with (REPO / "analysis/calibration/ctrpv2_ferroptosis_curves.csv").open(newline="") as f:
        rows = [row for row in csv.DictReader(f) if row["CompoundName"] == "ERASTIN"]
    assert len(rows) == 795
    assert max(float(row["MaximumDose"]) for row in rows) == 66
    with pytest.raises(ValueError, match="100 uM: 0"):
        select_supported_cohort(rows, (0.1, 0.3, 1, 3, 10, 30, 100))
    retained, meta = select_supported_cohort(rows, (0.1, 0.3, 1, 3, 10, 30))
    assert meta["n_retained"] == 785
    assert meta["n_excluded"] == 10
    assert all(float(row["MinimumDose"]) <= 0.1 and float(row["MaximumDose"]) >= 30
               for row in retained)
