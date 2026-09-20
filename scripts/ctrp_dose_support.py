"""Select a fixed CTRPv2 cohort supported across a complete dose grid.

The committed CTRPv2 rows contain fitted dose-response curves, not raw assay
replicates. Evaluating a curve inside its recorded dose range is interpolation
of that fit; outside it is extrapolation. Calibration targets use only the
former. Selecting rows once for the entire grid keeps the cell-line denominator
fixed across doses, including when individual lines have different ranges.
"""

import math


def _positive_finite(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a positive finite number") from exc
    if isinstance(value, bool) or not math.isfinite(number) or number <= 0:
        raise ValueError(f"{label} must be a positive finite number")
    return number


def select_supported_cohort(rows, doses):
    """Return ``(retained_rows, provenance)`` for one compound's dose grid.

    Retained rows are the original objects, in their original order. Both
    recorded range endpoints are inclusive. Every row must have valid positive
    finite ``MinimumDose``/``MaximumDose`` values and exactly ``DoseUnit='uM'``;
    malformed metadata raises rather than silently removing observations.

    Provenance reports support at each requested dose BEFORE selecting the
    common cohort, plus the rows excluded from that cohort and their reasons.
    A valid row can support individual doses without supporting the whole grid.
    Empty grids, input collections, or resulting cohorts raise ``ValueError``.
    This selects assay support only; it does not establish independent cell-line
    holdouts or make fitted curve values into raw experimental measurements.
    """
    grid = [_positive_finite(d, f"requested dose {i}") for i, d in enumerate(doses)]
    if not grid:
        raise ValueError("requested dose grid is empty")
    curve_rows = list(rows)
    if not curve_rows:
        raise ValueError("input curve cohort is empty")

    compounds = {r["CompoundName"] for r in curve_rows if "CompoundName" in r}
    if len(compounds) > 1:
        raise ValueError("select one compound at a time; mixed CompoundName values")

    ranges = []
    for i, row in enumerate(curve_rows):
        if row.get("DoseUnit") != "uM":
            raise ValueError(f"row {i}: DoseUnit must be 'uM'")
        low = _positive_finite(row.get("MinimumDose"), f"row {i}: MinimumDose")
        high = _positive_finite(row.get("MaximumDose"), f"row {i}: MaximumDose")
        if low > high:
            raise ValueError(f"row {i}: MinimumDose exceeds MaximumDose")
        ranges.append((low, high))

    per_dose_support = [
        {"dose_um": dose, "n_supported": sum(low <= dose <= high for low, high in ranges)}
        for dose in grid
    ]
    grid_low, grid_high = min(grid), max(grid)
    retained, excluded = [], []
    for i, (row, (low, high)) in enumerate(zip(curve_rows, ranges)):
        reasons = []
        if grid_low < low:
            reasons.append("requested_dose_below_minimum")
        if grid_high > high:
            reasons.append("requested_dose_above_maximum")
        if not reasons:
            retained.append(row)
        else:
            exclusion = {
                "row_index": i,
                "minimum_dose_um": low,
                "maximum_dose_um": high,
                "reasons": reasons,
            }
            for key in ("ModelID", "SampleID"):
                if key in row:
                    exclusion[key] = row[key]
            excluded.append(exclusion)

    if not retained:
        support = ", ".join(f"{p['dose_um']:g} uM: {p['n_supported']}"
                            for p in per_dose_support)
        raise ValueError("no curve supports the entire requested dose grid; "
                         f"per-dose support: {support}")

    metadata = {
        "policy": "entire_grid_within_recorded_dose_range",
        "n_original": len(curve_rows),
        "n_retained": len(retained),
        "n_excluded": len(excluded),
        "requested_doses_um": grid,
        "per_dose_support": per_dose_support,
        "excluded_rows": excluded,
    }
    if compounds:
        metadata["compound"] = next(iter(compounds))
    return retained, metadata
