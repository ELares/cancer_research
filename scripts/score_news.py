#!/usr/bin/env python3
"""Compute credibility score for news articles.

Usage:
    python score_news.py article.md
    python score_news.py --all
"""

import argparse
import sys
from datetime import date, datetime
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
# Scoring
# ---------------------------------------------------------------------------

def _months_since(date_str: str | None) -> float | None:
    """Return the number of months between *date_str* (YYYY-MM-DD) and today.

    Returns None if *date_str* is missing or unparseable.
    """
    if not date_str:
        return None
    try:
        pub_date = datetime.strptime(str(date_str)[:10], "%Y-%m-%d").date()
    except (ValueError, TypeError):
        return None
    delta = date.today() - pub_date
    return delta.days / 30.44  # average days per month


def compute_score(fm: dict) -> float:
    """Compute a 0-100 credibility score from article frontmatter.

    Weights:
        40 % -- verified-claim ratio (among FACTUAL claims)
        30 % -- author credentialing
        20 % -- recency
        10 % -- cross-citation (deferred; always 0.0 in v1)

    The final score is multiplied by a tier weight:
        tier 1 -> 1.0, tier 2 -> 0.8, tier 3 -> 0.6, other -> 0.3
    """
    tier = fm.get("tier", 0)
    tier_weight: float = {1: 1.0, 2: 0.8, 3: 0.6}.get(tier, 0.3)

    # --- Verified-claim ratio ---
    claims = fm.get("claims", [])
    factual_claims = [c for c in claims if c.get("category") == "FACTUAL"]
    if factual_claims:
        # Both "verified" and "self-referencing" count as verified for scoring.
        # Per criteria doc: self-referencing sources ARE the authority.
        verified_count = sum(
            1 for c in factual_claims
            if c.get("verification_status") in ("verified", "self-referencing")
        )
        verified_ratio: float = verified_count / len(factual_claims)
    else:
        verified_ratio = 1.0  # no factual claims -> not penalised

    # --- Author score ---
    if fm.get("author_credentialed"):
        author_score: float = 1.0
    elif fm.get("author"):
        author_score = 0.7
    else:
        author_score = 0.3

    # --- Recency ---
    months = _months_since(fm.get("date_published"))
    if months is None:
        recency: float = 0.5  # unknown date -- middle ground
    elif months < 6:
        recency = 1.0
    elif months < 12:
        recency = 0.8
    elif months < 36:
        recency = 0.5
    else:
        recency = 0.2

    # --- Cross-citation (deferred for v1) ---
    cross_citation: float = 0.0

    score = tier_weight * (
        40 * verified_ratio
        + 30 * author_score
        + 20 * recency
        + 10 * cross_citation
    )
    return round(score, 1)


# ---------------------------------------------------------------------------
# Article-level scoring
# ---------------------------------------------------------------------------

def score_article(article_path: Path) -> float | None:
    """Load an article, compute its credibility score, and save back.

    Returns:
        The computed score, or None if the article could not be loaded.
    """
    fm, body = load_article(article_path)
    if not fm:
        print(f"  skipping (no frontmatter): {article_path.name}")
        return None

    score = compute_score(fm)
    fm["credibility_score"] = score
    save_article(article_path, fm, body)
    return score


# ---------------------------------------------------------------------------
# Batch helpers
# ---------------------------------------------------------------------------

def find_all_articles() -> list[Path]:
    """Return all news article paths under news/by-source/."""
    source_dir = NEWS_DIR / "by-source"
    if not source_dir.exists():
        return []
    return sorted(source_dir.glob("**/*.md"))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Score news article credibility.")
    parser.add_argument("article", nargs="?", help="Path to a single article")
    parser.add_argument("--all", action="store_true",
                        help="Score all news articles")
    args = parser.parse_args()

    if not args.article and not args.all:
        parser.error("Provide an article path or --all")

    if args.article:
        path = Path(args.article).resolve()
        score = score_article(path)
        if score is not None:
            print(f"{path.name}: credibility_score = {score}")
        return

    # --all mode
    articles = find_all_articles()
    print(f"Scoring {len(articles)} articles...")

    scores: list[float] = []
    for article_path in tqdm(articles, desc="  Scoring"):
        s = score_article(article_path)
        if s is not None:
            scores.append(s)

    if scores:
        avg = sum(scores) / len(scores)
        print(f"\nDone. Scored {len(scores)} articles.")
        print(f"  Mean score: {avg:.1f}")
        print(f"  Range: {min(scores):.1f} - {max(scores):.1f}")
    else:
        print("\nNo articles to score.")


if __name__ == "__main__":
    main()
