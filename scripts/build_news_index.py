#!/usr/bin/env python3
"""Build claim-centric NEWS_INDEX.jsonl from processed articles.

Usage:
    python build_news_index.py
    python build_news_index.py --approved-only
"""

import argparse
import json
import sys
from pathlib import Path

from tqdm import tqdm

# Allow imports from the scripts directory
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import (
    SOURCE_TIER_DEFINITIONS,
    NEWS_RATE,
    NEWS_DIR,
    CLAIM_FACTUAL_MARKERS,
    CLAIM_TYPE_MARKERS,
    resilient_get,
    PUBMED_ESEARCH,
    NCBI_API_KEY,
    NCBI_RATE,
    INDEX_FILE,
)
from article_io import load_article, save_article


# ---------------------------------------------------------------------------
# Index building
# ---------------------------------------------------------------------------

def build_index(approved_only: bool = False) -> int:
    """Build NEWS_INDEX.jsonl from all processed news articles.

    Each JSONL entry represents one claim, denormalised with its parent
    article metadata so consumers can filter/sort without joining files.

    Args:
        approved_only: If True, skip articles whose ``review_status``
            is not ``"approved"``.

    Returns:
        Number of claim entries written.
    """
    source_dir = NEWS_DIR / "by-source"
    if not source_dir.exists():
        print(f"No news articles found at {source_dir}")
        return 0

    article_paths = sorted(source_dir.glob("**/*.md"))
    if not article_paths:
        print("No .md files found under news/by-source/")
        return 0

    entries: list[dict] = []

    for article_path in tqdm(article_paths, desc="  Reading articles"):
        fm, _ = load_article(article_path)
        if not fm:
            continue

        if approved_only and fm.get("review_status") != "approved":
            continue

        claims = fm.get("claims", [])
        if not claims:
            continue

        # Shared article-level fields
        article_meta = {
            "article_url": fm.get("url", ""),
            "source_domain": fm.get("source_domain", ""),
            "tier": fm.get("tier"),
            "article_title": fm.get("title", ""),
            "author": fm.get("author", ""),
            "date_published": fm.get("date_published"),
            "review_status": fm.get("review_status", "pending"),
            "credibility_score": fm.get("credibility_score"),
        }

        for claim in claims:
            entry = {
                "claim_id": claim.get("id", ""),
                "claim_text": claim.get("text", ""),
                "claim_type": claim.get("type", ""),
                "claim_category": claim.get("category", ""),
                "verification_status": claim.get("verification_status", "unverified"),
                "linked_pmids": claim.get("linked_pmids", []),
            }
            entry.update(article_meta)
            entries.append(entry)

    # Sort: dated entries newest-first, undated entries at the end
    dated = [e for e in entries if e.get("date_published")]
    undated = [e for e in entries if not e.get("date_published")]
    dated.sort(key=lambda e: e["date_published"], reverse=True)
    entries = dated + undated

    # Write index
    index_path = NEWS_DIR / "NEWS_INDEX.jsonl"
    NEWS_DIR.mkdir(parents=True, exist_ok=True)
    with open(index_path, "w", encoding="utf-8") as f:
        for entry in entries:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print(f"\nIndex written to {index_path}")
    print(f"  Total claims: {len(entries)}")
    print(f"  Articles processed: {len(article_paths)}")

    # Summary stats
    verified = sum(1 for e in entries if e.get("verification_status") == "verified")
    factual = sum(1 for e in entries if e.get("claim_category") == "FACTUAL")
    print(f"  Verified claims: {verified}/{factual} factual")

    tiers = {}
    for e in entries:
        t = e.get("tier", 0)
        tiers[t] = tiers.get(t, 0) + 1
    for t in sorted(tiers):
        print(f"  Tier {t}: {tiers[t]} claims")

    return len(entries)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build claim-centric NEWS_INDEX.jsonl."
    )
    parser.add_argument("--approved-only", action="store_true",
                        help="Include only articles with review_status='approved'")
    args = parser.parse_args()

    build_index(approved_only=args.approved_only)


if __name__ == "__main__":
    main()
