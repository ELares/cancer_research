"""
Unit tests for `simulations/calibration/calibrate.py`'s 3D extractor.

Locks the row-filter invariants of `_find_3d_condition` and the
named-metric dispatch of `extract_tme_3d_json`. Fixtures build minimal
`summary.json`-shaped dicts in-memory (no on-disk dependency) so the
tests pin the contract: if sim-tme-3d's `ConditionResult` schema drifts
and the extractor silently returns `None`, these tests fail loudly
instead.

Covers issue #222 item 1.

Run: pytest tests/test_calibrate.py -v
"""

import json
import sys
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
CALIBRATION_DIR = REPO_ROOT / "simulations" / "calibration"
sys.path.insert(0, str(CALIBRATION_DIR))

import calibrate  # noqa: E402


# ============================================================
# Fixture: minimal summary.json-shaped dict
# ============================================================

def _condition(
    *,
    treatment,
    o2_condition="uniform",
    o2_lambda_um=None,
    immune_mode="off",
    stromal_mode=None,
    ph_mode=None,
    normoxic_kill_rate=0.0,
    hypoxic_kill_rate=0.0,
    immune_kills=None,
    stromal_adjacent_kill_rate=None,
):
    """Build one `ConditionResult`-shaped dict.

    Mirrors the fields of `simulations/sim-tme-3d/src/main.rs`'s
    `ConditionResult` struct that the extractor actually reads. Other
    fields (total_tumor, peak_damp, etc.) are omitted because the
    extractor doesn't consult them — this keeps the fixture minimal and
    makes the test failures point to schema fields that actually matter.
    """
    return {
        "treatment": treatment,
        "o2_condition": o2_condition,
        "o2_lambda_um": o2_lambda_um,
        "immune_mode": immune_mode,
        "stromal_mode": stromal_mode,
        "ph_mode": ph_mode,
        "normoxic_kill_rate": normoxic_kill_rate,
        "hypoxic_kill_rate": hypoxic_kill_rate,
        "immune_kills": immune_kills,
        "stromal_adjacent_kill_rate": stromal_adjacent_kill_rate,
    }


@pytest.fixture
def summary_dict():
    """A summary.json-shaped dict with rows for all three Q1/Q2/Q3 metrics."""
    return {
        "conditions": [
            # Row 0: RSL3, λ=120, immune_off — for rsl3_o2_collapse_ratio
            _condition(
                treatment="RSL3",
                o2_condition="gradient",
                o2_lambda_um=120.0,
                immune_mode="off",
                normoxic_kill_rate=0.80,
                hypoxic_kill_rate=0.20,
            ),
            # Row 1: SDT, λ=120, immune_on — numerator of immune_sdt_rsl3_ratio
            _condition(
                treatment="SDT",
                o2_condition="gradient",
                o2_lambda_um=120.0,
                immune_mode="immune_on",
                immune_kills=400,
            ),
            # Row 2: RSL3, λ=120, immune_on, no stromal — denominator of immune ratio
            # AND no-stromal baseline for stromal_shielding_ratio
            _condition(
                treatment="RSL3",
                o2_condition="gradient",
                o2_lambda_um=120.0,
                immune_mode="immune_on",
                immune_kills=100,
                stromal_adjacent_kill_rate=0.50,
            ),
            # Row 3: RSL3, λ=120, immune_on, stromal_on — numerator of stromal ratio
            _condition(
                treatment="RSL3",
                o2_condition="gradient",
                o2_lambda_um=120.0,
                immune_mode="immune_on",
                stromal_mode="stromal_on",
                stromal_adjacent_kill_rate=0.10,
            ),
        ],
        "note": "test fixture",
    }


def _write_summary(tmp_path, summary):
    """Write the fixture dict to disk and return a target dict pointing at it."""
    out_dir = tmp_path / "tme-3d"
    out_dir.mkdir()
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary))
    return summary_path


def _target(metric, summary_path):
    """Build a target dict in the shape `extract_tme_3d_json` expects."""
    return {
        "binary": "sim-tme-3d",
        "output_file": str(summary_path),
        "extraction": {"metric": metric},
    }


@pytest.fixture(autouse=True)
def _patch_resolve(monkeypatch):
    """Make `_resolve_output_path` honor an absolute path in output_file."""
    def fake_resolve(target):
        return Path(target["output_file"])
    monkeypatch.setattr(calibrate, "_resolve_output_path", fake_resolve)


# ============================================================
# Metric 1: rsl3_o2_collapse_ratio (single-row, two-field derivation)
# Pins: stromal_mode=None / ph_mode=None filter, immune_off filter,
# hypoxic/normoxic derivation.
# ============================================================

class TestRsl3O2CollapseRatio:
    def test_returns_hypoxic_over_normoxic(self, tmp_path, summary_dict):
        path = _write_summary(tmp_path, summary_dict)
        result = calibrate.extract_tme_3d_json(_target("rsl3_o2_collapse_ratio", path))
        assert result == pytest.approx(0.20 / 0.80)

    def test_skips_immune_on_rows(self, tmp_path, summary_dict):
        # Remove the immune_off RSL3 row — only immune_on RSL3 left.
        # Extractor must return None (not silently pick the immune_on row).
        summary_dict["conditions"] = [
            c for c in summary_dict["conditions"]
            if not (c["treatment"] == "RSL3" and c["immune_mode"] == "off")
        ]
        path = _write_summary(tmp_path, summary_dict)
        result = calibrate.extract_tme_3d_json(_target("rsl3_o2_collapse_ratio", path))
        assert result is None

    def test_skips_stromal_on_rows(self, tmp_path, summary_dict):
        # Add a stromal_on RSL3 immune_off row and remove the no-stromal one.
        # `_find_3d_condition` defaults stromal_mode=None — must not match.
        summary_dict["conditions"] = [
            c for c in summary_dict["conditions"]
            if not (c["treatment"] == "RSL3" and c["immune_mode"] == "off")
        ]
        summary_dict["conditions"].append(_condition(
            treatment="RSL3",
            o2_condition="gradient",
            o2_lambda_um=120.0,
            immune_mode="off",
            stromal_mode="stromal_on",
            normoxic_kill_rate=0.80,
            hypoxic_kill_rate=0.20,
        ))
        path = _write_summary(tmp_path, summary_dict)
        result = calibrate.extract_tme_3d_json(_target("rsl3_o2_collapse_ratio", path))
        assert result is None


# ============================================================
# Metric 2: immune_sdt_rsl3_ratio (two-row, single-field per row)
# Pins: matched-context filter (both rows at immune_on, λ=120,
# stromal=None, ph=None), divide-by-zero guard, missing-row handling.
# ============================================================

class TestImmuneSdtRsl3Ratio:
    def test_returns_sdt_over_rsl3_immune_kills(self, tmp_path, summary_dict):
        path = _write_summary(tmp_path, summary_dict)
        result = calibrate.extract_tme_3d_json(_target("immune_sdt_rsl3_ratio", path))
        assert result == pytest.approx(400.0 / 100.0)

    def test_skips_when_rsl3_immune_kills_zero(self, tmp_path, summary_dict):
        # Zero RSL3 immune kills must not raise ZeroDivisionError.
        for c in summary_dict["conditions"]:
            if c["treatment"] == "RSL3" and c["immune_mode"] == "immune_on":
                c["immune_kills"] = 0
        path = _write_summary(tmp_path, summary_dict)
        result = calibrate.extract_tme_3d_json(_target("immune_sdt_rsl3_ratio", path))
        assert result is None


# ============================================================
# Metric 3: stromal_shielding_ratio (two-row, stromal_mode filter)
# Pins: stromal_mode=None vs "stromal_on" partitioning — the row-filter
# invariant most likely to silently break if the schema changes.
# ============================================================

class TestStromalShieldingRatio:
    def test_returns_stromal_on_over_no_stromal(self, tmp_path, summary_dict):
        path = _write_summary(tmp_path, summary_dict)
        result = calibrate.extract_tme_3d_json(_target("stromal_shielding_ratio", path))
        assert result == pytest.approx(0.10 / 0.50)

    def test_skips_when_no_stromal_baseline_missing(self, tmp_path, summary_dict):
        # Drop the no-stromal RSL3 immune_on row. Without a denominator
        # the metric must SKIP, not silently fall back to stromal_on.
        summary_dict["conditions"] = [
            c for c in summary_dict["conditions"]
            if not (
                c["treatment"] == "RSL3"
                and c["immune_mode"] == "immune_on"
                and c["stromal_mode"] is None
            )
        ]
        path = _write_summary(tmp_path, summary_dict)
        result = calibrate.extract_tme_3d_json(_target("stromal_shielding_ratio", path))
        assert result is None


# ============================================================
# Issue #222 item 3: duplicate-row warning
# ============================================================

class TestDuplicateRowWarning:
    def test_warns_on_duplicate_match(self, tmp_path, summary_dict):
        # Inject a duplicate of the RSL3 immune_off row. The extractor
        # must still return a value (first-match-wins) but warn loudly.
        dup = dict(summary_dict["conditions"][0])
        summary_dict["conditions"].append(dup)
        path = _write_summary(tmp_path, summary_dict)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = calibrate.extract_tme_3d_json(
                _target("rsl3_o2_collapse_ratio", path)
            )
        assert result == pytest.approx(0.20 / 0.80)
        assert any("matched" in str(w.message) for w in caught), (
            f"expected duplicate-row warning, got: {[str(w.message) for w in caught]}"
        )
