#!/usr/bin/env python3
"""Verify factual claims against local corpus and PubMed.

Usage:
    python verify_news_claims.py article.md
    python verify_news_claims.py --all
"""

import argparse
import json
import re
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

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
# Corpus index
# ---------------------------------------------------------------------------

_corpus_cache: list[dict] | None = None


def load_corpus_index() -> list[dict]:
    """Load INDEX_FILE (corpus/INDEX.jsonl) and return a list of dicts.

    The result is cached after the first call so repeated verifications
    within the same process don't re-read the file.
    """
    global _corpus_cache
    if _corpus_cache is not None:
        return _corpus_cache

    if not INDEX_FILE.exists():
        print(f"  warning: corpus index not found at {INDEX_FILE}")
        _corpus_cache = []
        return _corpus_cache

    entries: list[dict] = []
    with open(INDEX_FILE, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                entries.append(json.loads(line))
    _corpus_cache = entries
    return _corpus_cache


# ---------------------------------------------------------------------------
# Local corpus search
# ---------------------------------------------------------------------------

def search_corpus(keywords: list[str], index: list[dict]) -> list[dict]:
    """Search the corpus index for entries whose title contains any keyword.

    Case-insensitive substring matching.  Returns a deduplicated list of
    matching entries (each containing at least ``pmid`` and ``title``).
    """
    if not keywords or not index:
        return []

    lower_keywords = [kw.lower() for kw in keywords if len(kw) >= 3]
    if not lower_keywords:
        return []

    seen_pmids: set[str] = set()
    matches: list[dict] = []

    for entry in index:
        title = (entry.get("title") or "").lower()
        if not title:
            continue
        for kw in lower_keywords:
            if kw in title:
                pmid = entry.get("pmid", "")
                if pmid and pmid not in seen_pmids:
                    seen_pmids.add(pmid)
                    matches.append(entry)
                break  # one keyword match is enough per entry

    return matches


# ---------------------------------------------------------------------------
# PubMed search
# ---------------------------------------------------------------------------

def search_pubmed(query: str, max_results: int = 5) -> list[dict]:
    """Search PubMed via ESearch and return a list of {pmid, title} dicts.

    Uses resilient_get with NCBI_RATE for polite access.
    """
    params: dict = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        resp = resilient_get(PUBMED_ESEARCH, params=params, rate_limiter=NCBI_RATE)
        resp.raise_for_status()
    except Exception as exc:
        print(f"  PubMed search failed: {exc}")
        return []

    data = resp.json()
    pmids = data.get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return []

    # Fetch titles via efetch
    return _fetch_titles(pmids)


def _fetch_titles(pmids: list[str]) -> list[dict]:
    """Fetch article titles from PubMed for a list of PMIDs."""
    from config import PUBMED_EFETCH  # avoid circular at module level

    params: dict = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "rettype": "abstract",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        resp = resilient_get(PUBMED_EFETCH, params=params, rate_limiter=NCBI_RATE)
        resp.raise_for_status()
    except Exception as exc:
        print(f"  PubMed efetch failed: {exc}")
        return [{"pmid": p, "title": ""} for p in pmids]

    results: list[dict] = []
    try:
        root = ET.fromstring(resp.text)
        for article_el in root.findall(".//PubmedArticle"):
            pmid_el = article_el.find(".//PMID")
            title_el = article_el.find(".//ArticleTitle")
            results.append({
                "pmid": pmid_el.text if pmid_el is not None else "",
                "title": title_el.text if title_el is not None else "",
            })
    except ET.ParseError:
        results = [{"pmid": p, "title": ""} for p in pmids]

    return results


# ---------------------------------------------------------------------------
# Search-term extraction
# ---------------------------------------------------------------------------

# Words too common to be useful as search terms
_STOP_WORDS = {
    "the", "a", "an", "and", "or", "but", "in", "on", "at", "to", "for",
    "of", "with", "by", "from", "is", "are", "was", "were", "be", "been",
    "has", "have", "had", "that", "this", "it", "its", "as", "not", "can",
    "will", "new", "also", "may", "more", "than", "they", "their", "which",
    "who", "said", "would", "could", "about", "into", "over", "after",
    "before", "between", "such", "most", "only", "other", "some", "all",
}


def extract_search_terms(claim_text: str) -> list[str]:
    """Extract focused search terms from a claim sentence.

    Looks for:
    - Capitalised multi-word terms (proper nouns, drug names)
    - Numbers paired with clinical keywords (e.g. "Phase 3")
    """
    terms: list[str] = []

    # Multi-word capitalised phrases (e.g. "Keytruda", "FDA", "Phase III")
    caps = re.findall(r"\b[A-Z][a-z]*(?:\s+[A-Z][a-z]*)+\b", claim_text)
    terms.extend(caps)

    # Individual capitalised words that are likely proper nouns / drug names
    single_caps = re.findall(r"\b[A-Z][a-z]{2,}\b", claim_text)
    for w in single_caps:
        if w.lower() not in _STOP_WORDS:
            terms.append(w)

    # All-caps acronyms (FDA, OS, PFS, etc.)
    acronyms = re.findall(r"\b[A-Z]{2,6}\b", claim_text)
    for a in acronyms:
        if a.lower() not in _STOP_WORDS and len(a) >= 2:
            terms.append(a)

    # Deduplicate while preserving order
    seen: set[str] = set()
    unique: list[str] = []
    for t in terms:
        key = t.lower()
        if key not in seen:
            seen.add(key)
            unique.append(t)

    return unique


# ---------------------------------------------------------------------------
# Claim verification
# ---------------------------------------------------------------------------

def verify_claim(claim: dict, corpus_index: list[dict]) -> dict:
    """Verify a single claim against local corpus, then PubMed.

    Only FACTUAL claims are verified.  Other categories are returned
    unchanged.

    Updates ``verification_status``, ``verification_source``, and
    ``linked_pmids`` on the claim dict (mutated in place and returned).
    """
    if claim.get("category") != "FACTUAL":
        return claim

    terms = extract_search_terms(claim.get("text", ""))
    if not terms:
        return claim

    # --- Local corpus search ---
    corpus_hits = search_corpus(terms, corpus_index)
    if corpus_hits:
        claim["verification_status"] = "verified"
        claim["verification_source"] = "corpus"
        claim["linked_pmids"] = [h["pmid"] for h in corpus_hits[:5]]
        return claim

    # --- PubMed fallback ---
    query = " ".join(terms[:5])  # keep query short
    pubmed_hits = search_pubmed(query, max_results=5)
    if pubmed_hits:
        claim["verification_status"] = "verified"
        claim["verification_source"] = "pubmed"
        claim["linked_pmids"] = [h["pmid"] for h in pubmed_hits if h.get("pmid")]
        return claim

    # No match
    claim["verification_status"] = "unverified"
    return claim


# ---------------------------------------------------------------------------
# Article-level verification
# ---------------------------------------------------------------------------

def verify_article(article_path: Path) -> int:
    """Load an article, verify each claim, save back.

    Returns:
        Number of claims whose status changed.
    """
    fm, body = load_article(article_path)
    if not fm:
        print(f"  skipping (no frontmatter): {article_path.name}")
        return 0

    claims = fm.get("claims", [])
    if not claims:
        return 0

    corpus_index = load_corpus_index()
    changed = 0

    for claim in claims:
        old_status = claim.get("verification_status")
        verify_claim(claim, corpus_index)
        if claim.get("verification_status") != old_status:
            changed += 1

    fm["claims"] = claims
    save_article(article_path, fm, body)
    return changed


# ---------------------------------------------------------------------------
# Batch helpers
# ---------------------------------------------------------------------------

def _has_unverified_factual(article_path: Path) -> bool:
    """True if the article has any FACTUAL claims still unverified."""
    fm, _ = load_article(article_path)
    if not fm:
        return False
    for claim in fm.get("claims", []):
        if claim.get("category") == "FACTUAL" and claim.get("verification_status") == "unverified":
            return True
    return False


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
    parser = argparse.ArgumentParser(description="Verify news article claims.")
    parser.add_argument("article", nargs="?", help="Path to a single article")
    parser.add_argument("--all", action="store_true",
                        help="Verify all articles with unverified factual claims")
    args = parser.parse_args()

    if not args.article and not args.all:
        parser.error("Provide an article path or --all")

    if args.article:
        path = Path(args.article).resolve()
        changed = verify_article(path)
        print(f"Verified {path.name}: {changed} claim(s) updated")
        return

    # --all mode
    articles = find_all_articles()
    pending = [a for a in articles if _has_unverified_factual(a)]
    print(f"Articles with unverified factual claims: {len(pending)}/{len(articles)}")

    total_changed = 0
    for article_path in tqdm(pending, desc="  Verifying"):
        changed = verify_article(article_path)
        total_changed += changed

    print(f"\nDone. Updated {total_changed} claims across {len(pending)} articles.")


if __name__ == "__main__":
    main()
