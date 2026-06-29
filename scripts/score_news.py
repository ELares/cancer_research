#!/usr/bin/env python3
"""Compute credibility score for news articles.

Usage:
    python score_news.py article.md
    python score_news.py --all
"""

import argparse
import json
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

def _parse_date(date_str: str | None) -> "date | None":
    """Parse a YYYY-MM-DD string to a date, or None if missing/unparseable."""
    if not date_str:
        return None
    try:
        return datetime.strptime(str(date_str)[:10], "%Y-%m-%d").date()
    except (ValueError, TypeError):
        return None


def _months_since(date_str: str | None, as_of: "date | None" = None) -> float | None:
    """Months between *date_str* (YYYY-MM-DD) and *as_of* (default today).

    Pass a fixed *as_of* (e.g. the article's recorded ``scored_at``) to make the
    recency term DETERMINISTIC instead of wall-clock-dependent: with the default
    today, an article silently crosses the 6/12/36-month recency buckets as time
    passes, which would mutate its committed credibility_score on a later
    ``--all`` rerun (#587).

    Returns None if *date_str* is missing or unparseable.
    """
    pub_date = _parse_date(date_str)
    if pub_date is None:
        return None
    ref = as_of or date.today()
    delta = ref - pub_date
    return delta.days / 30.44  # average days per month


_CORPUS_PMIDS: set[str] | None = None


def _corpus_pmids() -> set[str]:
    """Lazy-load the set of PMIDs in the local corpus (cached). Used by the
    cross-citation term to reward news claims anchored to a paper we actually
    hold (#532)."""
    global _CORPUS_PMIDS
    if _CORPUS_PMIDS is None:
        _CORPUS_PMIDS = set()
        if INDEX_FILE.exists():
            for line in INDEX_FILE.read_text(encoding="utf-8").splitlines():
                try:
                    _CORPUS_PMIDS.add(str(json.loads(line)["pmid"]))
                except (json.JSONDecodeError, KeyError):
                    continue
    return _CORPUS_PMIDS


def compute_score(fm: dict, as_of: "date | None" = None) -> float:
    """Compute a 0-100 credibility score from article frontmatter.

    *as_of* is the reference date for the recency term; pass a fixed date (or let
    it resolve from the article's recorded ``scored_at``) to keep the score
    deterministic rather than wall-clock-dependent (#587).

    Weights:
        40 % -- verified-claim ratio (among FACTUAL claims)
        30 % -- author credentialing
        20 % -- recency
        10 % -- cross-citation: count of DISTINCT corpus PMIDs the article's claims
                cite, full credit at >= 3 (#532; was a 0.0 stub; re-anchored off
                the verified-ratio-redundant anchored-fraction in #571)

    The final score is multiplied by a tier weight:
        tier 1 -> 1.0, tier 2 -> 0.8, tier 3 -> 0.6, other -> 0.3
    """
    tier = fm.get("tier", 0)
    tier_weight: float = {1: 1.0, 2: 0.8, 3: 0.6}.get(tier, 0.3)

    # --- Verified-claim ratio ---
    # `or []` (not a default arg): a present-but-null `claims:` field returns
    # None, which would crash both this comprehension and the cross-citation
    # one below; one hand-edited / extraction-skipped article must not kill the
    # whole `--all` batch.
    claims = fm.get("claims") or []
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

    # --- Recency (frozen, not wall-clock: #587) ---
    # Reference date precedence: an explicit `as_of`, else the article's recorded
    # `scored_at`, else today (first scoring). This is what makes a later `--all`
    # rerun reproduce the same recency bucket instead of drifting.
    ref = as_of or _parse_date(fm.get("scored_at")) or date.today()
    months = _months_since(fm.get("date_published"), ref)
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

    # --- Cross-citation: BREADTH of corpus grounding (#532, re-anchored #571) ---
    # The earlier version (fraction of linked claims anchored to the corpus) was
    # largely redundant with the verified-ratio term: verify_news_claims sets a
    # corpus-verified claim's linked_pmids TO corpus PMIDs, so those claims always
    # scored 1.0. Measure instead the number of DISTINCT corpus papers the article
    # cites across all its claims — independent of per-claim verification: an
    # article grounded in several different corpus papers scores higher than one
    # that re-cites the same paper. Full credit at >= 3 distinct corpus papers.
    # Caveat (honesty): verify_news_claims links each corpus-verified claim to up
    # to 5 of its title-matched corpus_hits, so a single strongly-matched claim can
    # already supply >= 3 distinct PMIDs and saturate this term. It therefore
    # rewards corpus-citation breadth but is not a pure "distinct claims" measure;
    # it is a 10% heuristic, and the de-duplication is what removes the old
    # verified_ratio overlap (#571).
    corpus = _corpus_pmids()
    distinct_corpus = {
        str(p) for c in claims for p in (c.get("linked_pmids") or [])
        if str(p) in corpus
    }
    cross_citation: float = min(1.0, len(distinct_corpus) / 3.0)

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

    # Freeze the recency reference date at first scoring (#587): once `scored_at`
    # is recorded, `compute_score` reads it instead of the wall clock, so a later
    # rerun reproduces the same credibility_score. The one-time migration stamps
    # today; that preserves an already-scored article's score precisely when its
    # recency bucket has not shifted since it was last scored (verified
    # empirically at migration — all 38 scores byte-identical).
    fm.setdefault("scored_at", date.today().isoformat())
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
