#!/usr/bin/env python3
"""Extract claims from fetched news articles.

Usage:
    python extract_claims.py article.md      # Process one article
    python extract_claims.py --all           # Process all articles without claims
"""

import argparse
import re
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
# Body text cleaning (pre-sentence-splitting)
# ---------------------------------------------------------------------------

# Matches ScienceDaily release-path fragments: "260420014746.htm" etc.
_HTM_PATH_RE = re.compile(r"\d{8,}\.htm\b")
# Matches bare year or month lines that are URL path residue
_BARE_YEAR_RE = re.compile(r"^\d{4}$")
_BARE_MONTH_RE = re.compile(r"^\d{1,2}$")


def _clean_body_text(text: str) -> str:
    """Remove URL path fragments and other extraction artifacts from body text.

    ScienceDaily pages embed related-story URLs whose path components
    (e.g., "2026", "04", "260420014746.htm") survive HTML extraction as
    separate lines.  This function strips them before sentence splitting
    so they don't contaminate claim text.
    """
    lines = text.split("\n")
    cleaned = []
    for line in lines:
        stripped = line.strip()
        # Remove lines that are just URL path fragments
        if _HTM_PATH_RE.search(stripped):
            continue
        if _BARE_YEAR_RE.match(stripped):
            continue
        if _BARE_MONTH_RE.match(stripped):
            continue
        cleaned.append(line)
    return "\n".join(cleaned)


# ---------------------------------------------------------------------------
# Sentence splitting
# ---------------------------------------------------------------------------

# Abbreviations that should NOT trigger a sentence break
_ABBREV_PATTERN = re.compile(
    r"\b(?:Dr|Mr|Mrs|Ms|Prof|Sr|Jr|vs|etc|al|approx|dept|est|govt|incl)\."
)


def split_sentences(text: str) -> list[str]:
    """Split text into sentences on '.', '!', '?' boundaries.

    Handles common abbreviations (Dr., et al., vs.) and decimal numbers
    (e.g. 3.7%) so they are not treated as sentence boundaries.
    """
    # Replace abbreviation dots with a placeholder
    protected = _ABBREV_PATTERN.sub(lambda m: m.group(0).replace(".", "\x00"), text)

    # Protect decimal numbers  (e.g. "3.7%", "0.05", "1.2-fold")
    protected = re.sub(r"(\d)\.(\d)", lambda m: m.group(1) + "\x00" + m.group(2), protected)

    # Split on sentence-ending punctuation followed by whitespace + uppercase
    # or end-of-string
    parts = re.split(r'(?<=[.!?])\s+(?=[A-Z"])', protected)

    sentences: list[str] = []
    for part in parts:
        restored = part.replace("\x00", ".").strip()
        if restored:
            sentences.append(restored)

    return sentences


# ---------------------------------------------------------------------------
# Claim detection helpers
# ---------------------------------------------------------------------------

# Pre-compile factual marker regexes for speed
_COMPILED_FACTUAL = [(pat, re.compile(pat, re.IGNORECASE)) for pat in CLAIM_FACTUAL_MARKERS]


def detect_factual_markers(sentence: str) -> list[str]:
    """Return the names (regex patterns) of all factual markers that match."""
    matched: list[str] = []
    for pattern_str, compiled in _COMPILED_FACTUAL:
        if compiled.search(sentence):
            matched.append(pattern_str)
    return matched


def classify_claim_type(sentence: str) -> str:
    """Classify a sentence's claim type using CLAIM_TYPE_MARKERS.

    Checks in priority order: event > result > mechanism > opinion >
    speculation.  Returns the type of the first match.  Defaults to
    "result" when factual markers are present but no type-specific
    keyword matches.
    """
    lower = sentence.lower()
    for claim_type in ("event", "result", "mechanism", "opinion", "speculation"):
        keywords = CLAIM_TYPE_MARKERS.get(claim_type, [])
        for kw in keywords:
            if kw in lower:
                return claim_type
    return "result"  # default when factual markers triggered extraction


def classify_claim_category(sentence: str, has_factual_markers: bool) -> str:
    """Classify a sentence into FACTUAL / SPECULATIVE / INTERPRETIVE.

    * FACTUAL -- factual markers detected.
    * SPECULATIVE -- speculation keywords and no factual markers.
    * INTERPRETIVE -- everything else.
    """
    if has_factual_markers:
        return "FACTUAL"

    lower = sentence.lower()
    for kw in CLAIM_TYPE_MARKERS.get("speculation", []):
        if kw in lower:
            return "SPECULATIVE"

    return "INTERPRETIVE"


# ---------------------------------------------------------------------------
# Claim extraction
# ---------------------------------------------------------------------------

def _make_claim_id(domain: str, date: str, slug: str, index: int) -> str:
    """Build a deterministic claim ID like ``statnews.com-2024-06-01-fda-approves-drug-001``."""
    date_part = date or "undated"
    return f"{domain}-{date_part}-{slug}-{index:03d}"


def _slug_from_path(article_path: Path) -> str:
    """Derive a short slug from the article filename."""
    stem = article_path.stem
    # Strip date prefix if present (YYYY-MM-DD-)
    stem = re.sub(r"^\d{4}-\d{2}-\d{2}-", "", stem)
    # Strip version suffix
    stem = re.sub(r"-v\d+$", "", stem)
    return stem[:40]


def extract_claims(article_path: Path) -> list[dict]:
    """Extract and tag claims from a single news article.

    Updates the article's YAML frontmatter ``claims`` field in-place and
    writes the file back to disk.

    Returns:
        List of claim dicts.
    """
    fm, body = load_article(article_path)
    if not fm:
        print(f"  skipping (no frontmatter): {article_path.name}")
        return []

    domain = fm.get("source_domain", "unknown")
    date = fm.get("date_published", "") or ""
    slug = _slug_from_path(article_path)

    # Clean body text: remove URL path fragments, bare year/month lines,
    # and other artifacts from HTML extraction before sentence splitting.
    body = _clean_body_text(body)

    sentences = split_sentences(body)
    claims: list[dict] = []
    claim_index = 1

    # Import high-specificity triggers for non-factual claims.
    from config import CLAIM_OPINION_TRIGGERS, CLAIM_SPECULATION_TRIGGERS

    for sentence in sentences:
        factual_markers = detect_factual_markers(sentence)
        has_factual = bool(factual_markers)

        # Path 1: Factual claims — triggered by factual markers (numbers,
        # endpoints, trial phases, etc.).
        if has_factual:
            claim_type = classify_claim_type(sentence)
            category = classify_claim_category(sentence, has_factual)
            claim = {
                "id": _make_claim_id(domain, date, slug, claim_index),
                "text": sentence,
                "type": claim_type,
                "category": category,
                "verification_status": "unverified",
                "verification_source": None,
                "linked_pmids": [],
            }
            claims.append(claim)
            claim_index += 1
            continue

        # Path 2: Non-factual claims — triggered by high-specificity
        # multi-word phrases only.  Single-word markers like "could" or
        # "said" are too broad and would match almost every sentence.
        sentence_lower = sentence.lower()

        is_opinion = any(t in sentence_lower for t in CLAIM_OPINION_TRIGGERS)
        is_speculation = any(t in sentence_lower for t in CLAIM_SPECULATION_TRIGGERS)

        if is_opinion or is_speculation:
            if is_speculation:
                category = "SPECULATIVE"
                claim_type = "speculation"
            else:
                category = "INTERPRETIVE"
                claim_type = "opinion"

            claim = {
                "id": _make_claim_id(domain, date, slug, claim_index),
                "text": sentence,
                "type": claim_type,
                "category": category,
                "verification_status": None,  # non-factual: not verified
                "verification_source": None,
                "linked_pmids": [],
            }
            claims.append(claim)
            claim_index += 1

    fm["claims"] = claims
    save_article(article_path, fm, body)

    return claims


# ---------------------------------------------------------------------------
# Batch helpers
# ---------------------------------------------------------------------------

def _needs_extraction(article_path: Path) -> bool:
    """True if the article has no claims extracted yet."""
    fm, _ = load_article(article_path)
    if not fm:
        return False
    return not fm.get("claims")


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
    parser = argparse.ArgumentParser(description="Extract claims from news articles.")
    parser.add_argument("article", nargs="?", help="Path to a single article")
    parser.add_argument("--all", action="store_true",
                        help="Process all articles that have no claims yet")
    args = parser.parse_args()

    if not args.article and not args.all:
        parser.error("Provide an article path or --all")

    if args.article:
        path = Path(args.article).resolve()
        claims = extract_claims(path)
        print(f"Extracted {len(claims)} claims from {path.name}")
        for c in claims:
            print(f"  [{c['category']}:{c['type']}] {c['text'][:80]}...")
        return

    # --all mode
    articles = find_all_articles()
    pending = [a for a in articles if _needs_extraction(a)]
    print(f"Articles needing claim extraction: {len(pending)}/{len(articles)}")

    total_claims = 0
    for article_path in tqdm(pending, desc="  Extracting"):
        claims = extract_claims(article_path)
        total_claims += len(claims)

    print(f"\nDone. Extracted {total_claims} claims from {len(pending)} articles.")


if __name__ == "__main__":
    main()
