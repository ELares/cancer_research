#!/usr/bin/env python3
"""
Rebuild INDEX.jsonl from all corpus articles.

Usage:
    python build_index.py
"""

import json
import sys
from pathlib import Path

from tqdm import tqdm

from article_io import load_frontmatter
from config import INDEX_FILE, PMID_DIR
from provenance import append_provenance_record


# Fields to include in the index (keep it lean for fast reads)
INDEX_FIELDS = [
    "pmid", "doi", "pmcid", "title", "journal", "year", "month",
    "is_oa", "oa_status", "cited_by_count",
    "mechanisms", "biology_processes", "pathway_targets", "radioligand_targets",
    "cancer_types", "cancer_subtypes", "tissue_categories", "evidence_level", "resistant_states", "combination_evidence",
    "icite_rcr", "icite_percentile", "icite_is_clinical",
    "date_added",
]


def sanitize_index_value(value):
    """Normalize strings so each JSONL record stays on one physical line."""
    if isinstance(value, str):
        return " ".join(value.split())
    return value


def main():
    files = sorted(PMID_DIR.glob("*.md"))
    if not files:
        print("No articles found in corpus.")
        return

    print(f"Building index from {len(files)} articles...")

    entries = []
    for filepath in tqdm(files, desc="  Reading"):
        fm = load_frontmatter(filepath)
        if not fm:
            continue

        entry = {}
        for field in INDEX_FIELDS:
            if field in fm:
                entry[field] = sanitize_index_value(fm[field])

        # Add author count (not full list — too large for index)
        authors = fm.get("authors", [])
        entry["author_count"] = len(authors)
        if authors:
            entry["first_author"] = sanitize_index_value(authors[0])

        entries.append(entry)

    # Sort by year (newest first), then by citations
    entries.sort(key=lambda e: (-(e.get("year") or 0), -(e.get("cited_by_count") or 0)))

    # Write INDEX.jsonl
    INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(INDEX_FILE, "w", encoding="utf-8") as f:
        for entry in entries:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print(f"\nIndex written to {INDEX_FILE}")
    print(f"  Total articles: {len(entries)}")

    # Stats
    years = [e["year"] for e in entries if e.get("year")]
    if years:
        print(f"  Year range: {min(years)} - {max(years)}")

    oa_count = sum(1 for e in entries if e.get("is_oa"))
    print(f"  Open access: {oa_count}/{len(entries)}")

    tagged = sum(1 for e in entries if e.get("mechanisms"))
    print(f"  With mechanism tags: {tagged}/{len(entries)}")

    with_rcr = sum(1 for e in entries if e.get("icite_rcr"))
    print(f"  With iCite RCR: {with_rcr}/{len(entries)}")

    append_provenance_record(
        "build_index.py",
        {
            "index_entry_count": len(entries),
            "index_year_min": min(years) if years else None,
            "index_year_max": max(years) if years else None,
            "index_open_access_count": oa_count,
            "index_mechanism_tagged_count": tagged,
            "index_rcr_count": with_rcr,
        },
    )
    print("  Provenance appended to analysis/provenance.jsonl")


if __name__ == "__main__":
    main()
