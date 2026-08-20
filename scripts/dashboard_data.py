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


# --- census layer (#RETIRE-FROZEN) ----------------------------------------
#
# The corpus layer above reads 4,830 RECORDS. The census is 5,187,265 and is
# gitignored, so record-level browsing of it is impossible in a browser and
# impossible for anyone who has not run the multi-hour ingest.
#
# What IS shippable is the committed AGGREGATES -- the analysis JSON the census
# scripts write, ~62 KB in total. That is the whole census's shape at a size a
# Pyodide page can load, and it is the layer a reader actually wants: nobody
# browses four million records one at a time.
#
# FAIL-SOFT, not fail-open. A missing artifact returns None and the caller says
# so, because rendering an empty census panel that looks populated is worse
# than rendering a notice. This is the same reason the simulation tab degrades
# to committed intervals rather than showing a blank sweep.

ANALYSIS = REPO_ROOT / "analysis"
CENSUS_ARTIFACTS = {
    "profile": "census-mechanism-profile.json",
    "growth": "census-mechanism-growth.json",
    "design": "census-evidence-design.json",
    "sites": "census-mechanism-sites.json",
    "chains": "census-diagnostic-chains.json",
}


def load_census(names=None, base=None):
    """Load the committed census aggregates. Returns {name: dict-or-None}."""
    base = Path(base) if base is not None else ANALYSIS
    want = names or list(CENSUS_ARTIFACTS)
    out = {}
    for key in want:
        path = base / CENSUS_ARTIFACTS[key]
        try:
            out[key] = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            out[key] = None
    return out


def census_headline(design):
    """The two denominators, both of them.

    Quoting the trial share against the whole census understates it and against
    classifiable records alone overstates it, so neither is returned alone --
    the caller gets both or nothing. Returns None if the artifact is missing.
    """
    if not design:
        return None
    total = design.get("census")
    classifiable = design.get("classifiable")
    trials = (design.get("classes") or {}).get("trial")
    if not (total and classifiable and trials is not None):
        return None
    return {
        "census": total,
        "classifiable": classifiable,
        "undetermined": total - classifiable,
        "trials": trials,
        "share_of_census": round(100 * trials / total, 2),
        "share_of_classifiable": round(100 * trials / classifiable, 2),
    }


def census_mechanism_rows(profile):
    """Per-mechanism rows for display, sorted by trial share.

    SORTED BY TRIAL SHARE, NOT BY VOLUME, and the reason is a finding rather
    than a preference: descriptor breadth varies enormously across mechanisms,
    so a volume ordering is substantially an ordering of how broad each
    descriptor is. Trial share is a ratio within a mechanism and does not have
    that problem.
    """
    if not profile:
        return []
    rows = []
    for r in profile.get("rows", []):
        rows.append({
            "mechanism": r["mechanism"],
            "census articles": r["census"],
            "clinical trials": r["trials"],
            "trial share %": r["trial_share"],
            "growth 2015-2025": r.get("growth"),
            "top site": (r["top_sites"][0]["site"] if r.get("top_sites") else None),
            "top partner": (r["top_partners"][0]["mechanism"]
                            if r.get("top_partners") else None),
        })
    rows.sort(key=lambda r: -(r["trial share %"] or 0))
    return rows
