#!/usr/bin/env python3
"""Pure data-loading + aggregation helpers for the corpus dashboard (#354).

Separated from the Streamlit app (`scripts/dashboard.py`) so the logic is
importable and unit-tested in CI WITHOUT Streamlit (a UI-only, non-pinned
dependency). Everything here is stdlib + reads the committed corpus index, so it
runs offline.
"""

import json
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
INDEX = REPO_ROOT / "corpus" / "INDEX.jsonl"

# Multi-valued list fields a record can be filtered/aggregated on.
LIST_FIELDS = ("mechanisms", "cancer_types", "biology_processes", "tissue_categories", "pathway_targets")


def load_index(path=INDEX):
    """Load corpus/INDEX.jsonl into a list of record dicts (skips blank/bad lines)."""
    records = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


def _as_list(rec, field):
    v = rec.get(field)
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


def value_counts(records, field):
    """Count occurrences of each value of `field` across records. For list fields
    each element is counted once per record; for scalar fields the value is counted."""
    c = Counter()
    for r in records:
        if field in LIST_FIELDS:
            for v in _as_list(r, field):
                c[v] += 1
        else:
            v = r.get(field)
            if v not in (None, ""):
                c[v] += 1
    return dict(c.most_common())


def year_histogram(records):
    """{year: count} over records with a usable integer year, sorted by year."""
    c = Counter()
    for r in records:
        y = r.get("year")
        try:
            c[int(y)] += 1
        except (TypeError, ValueError):
            continue
    return dict(sorted(c.items()))


def mechanism_cancer_matrix(records, top_mech=None, top_cancer=None):
    """Co-occurrence counts {(mechanism, cancer): n}. Optionally restrict to the
    top-N mechanisms / cancers by frequency (for a readable heatmap)."""
    mech_keep = set(list(value_counts(records, "mechanisms"))[:top_mech]) if top_mech else None
    canc_keep = set(list(value_counts(records, "cancer_types"))[:top_cancer]) if top_cancer else None
    matrix = Counter()
    for r in records:
        mechs = [m for m in _as_list(r, "mechanisms") if mech_keep is None or m in mech_keep]
        cancers = [c for c in _as_list(r, "cancer_types") if canc_keep is None or c in canc_keep]
        for m in mechs:
            for c in cancers:
                matrix[(m, c)] += 1
    return dict(matrix)


def filter_records(records, mechanisms=None, cancer_types=None, evidence_levels=None, year_range=None):
    """Return records matching ALL provided filters (AND across filter types, OR
    within a list filter). `year_range` is an inclusive (lo, hi) tuple or None."""
    def ok(r):
        if mechanisms and not (set(_as_list(r, "mechanisms")) & set(mechanisms)):
            return False
        if cancer_types and not (set(_as_list(r, "cancer_types")) & set(cancer_types)):
            return False
        if evidence_levels and r.get("evidence_level") not in evidence_levels:
            return False
        if year_range:
            try:
                y = int(r.get("year"))
            except (TypeError, ValueError):
                return False
            if not (year_range[0] <= y <= year_range[1]):
                return False
        return True

    return [r for r in records if ok(r)]


def summary_stats(records):
    """Headline counts for the dashboard header."""
    years = [int(r["year"]) for r in records if str(r.get("year", "")).strip().isdigit()]
    return {
        "n_records": len(records),
        "n_mechanisms": len(value_counts(records, "mechanisms")),
        "n_cancer_types": len(value_counts(records, "cancer_types")),
        "n_evidence_tagged": sum(1 for r in records if r.get("evidence_level")),
        "year_min": min(years) if years else None,
        "year_max": max(years) if years else None,
    }
