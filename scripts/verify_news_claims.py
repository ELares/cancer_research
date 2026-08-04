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
        hit_count = sum(1 for kw in lower_keywords if kw in title)
        if hit_count >= 2:
            pmid = entry.get("pmid", "")
            if pmid and pmid not in seen_pmids:
                seen_pmids.add(pmid)
                matches.append(entry)

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

# Capitalised words that are never proper nouns, so a sentence-initial one
# carries no search specificity even if the sentence-start rule is relaxed.
_NON_SPECIFIC = {
    "one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
    "ten", "eleven", "twelve", "first", "second", "third", "half", "both",
    "many", "several", "these", "those", "there", "here", "when", "while",
    "however", "although", "because", "since", "results", "result", "patients",
    "researchers", "scientists", "study", "studies", "data", "among", "during",
}

# The claim and a candidate paper must share at least this many content words
# before the paper counts as supporting it. Without a check of this kind, ANY
# non-empty search result was accepted as verification.
MIN_TITLE_OVERLAP = 2

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

    # Individual capitalised words that are likely proper nouns / drug names.
    # A word is NOT a proper noun merely because it opens a sentence: this rule
    # used to accept "Seven" from "Seven of these 26 patients had inoperable
    # tumors", producing the single-term query `Seven`, which matches ~836,000
    # PubMed records and returned the five most recently indexed.
    for m in re.finditer(r"\b[A-Z][a-z]{2,}\b", claim_text):
        before = claim_text[:m.start()].rstrip()
        sentence_initial = (not before) or before[-1] in ".!?;:"
        w = m.group(0)
        if sentence_initial or w.lower() in _STOP_WORDS or w.lower() in _NON_SPECIFIC:
            continue
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


def supports_claim(claim_text: str, title: str) -> bool:
    """Does a candidate paper share enough subject matter to count as support?

    Deliberately crude, and deliberately applied: a supporting citation shares
    vocabulary with the claim it supports. One that shares nothing is not
    evidence of anything, whatever the search engine returned.
    """
    def words(s: str) -> set:
        # Compared on a 5-character prefix so morphological variants match:
        # without it "tumors"/"tumor" and "immune"/"immunity" count as
        # different words and a genuinely supporting paper is rejected.
        return {w[:5] for w in re.findall(r"[a-z]{4,}", (s or "").lower())
                if w not in _STOP_WORDS and w not in _NON_SPECIFIC}
    return len(words(claim_text) & words(title)) >= MIN_TITLE_OVERLAP


# ---------------------------------------------------------------------------
# Claim verification
# ---------------------------------------------------------------------------

# Tier 1 institutional domains whose claims cite their own data.
# Claims from these sources are "self-referencing" — the institution IS
# the authority, so verification is legitimate but not independent.
_SELF_REFERENCING_DOMAINS = frozenset([
    "who.int", "cancer.gov", "fda.gov", "gco.iarc.fr",
    "clinicaltrials.gov", "nih.gov",
])


def _mark_unverified(claim: dict) -> dict:
    """Mark a claim unverified and DROP the evidence a previous run attached.

    Setting the status without clearing the fields leaves failed evidence in
    place: ``news_verification_audit.py`` reads every claim that carries
    ``linked_pmids`` regardless of status, and a consumer reading the index has
    no reason to expect identifiers hanging off a claim the pipeline has just
    declined to verify. The first re-run left 17 claims reading
    ``unverified`` while still advertising ``verification_source: pubmed`` and
    five PMIDs each.
    """
    claim["verification_status"] = "unverified"
    claim["verification_source"] = None
    claim["linked_pmids"] = []
    return claim


def verify_claim(
    claim: dict,
    corpus_index: list[dict],
    source_domain: str = "",
) -> dict:
    """Verify a single claim against local corpus, then PubMed.

    Only FACTUAL claims are verified.  Other categories are returned
    unchanged.

    Verification statuses (per criteria doc Step 3):
      - verified:          primary source found (PMID/DOI linked)
      - unverified:        no primary source found
      - self-referencing:  the news source IS the primary authority
                           (e.g., WHO citing IARC data)
      - disputed:          reserved for future use when primary source
                           contradicts the claim

    Updates ``verification_status``, ``verification_source``, and
    ``linked_pmids`` on the claim dict (mutated in place and returned).
    """
    if claim.get("category") != "FACTUAL":
        return claim

    # --- Self-referencing check ---
    # Tier 1 institutional sources that cite their own data are
    # "self-referencing" — legitimate but not independently verified.
    # Checked BEFORE term extraction: whether the publisher is its own authority
    # is a property of the publisher, not of what can be pulled out of the
    # sentence. Ordered the other way round, three WHO claims that yield no
    # search terms ("Approximately 38% of cancers can currently be prevented…")
    # fall through and are labelled unverified instead.
    domain_base = source_domain.removeprefix("www.")
    if domain_base in _SELF_REFERENCING_DOMAINS:
        claim["verification_status"] = "self-referencing"
        claim["verification_source"] = f"institutional authority ({domain_base})"
        return claim

    terms = extract_search_terms(claim.get("text", ""))
    if not terms:
        # No distinguishing term, so there is nothing to search for and this
        # pipeline cannot verify the claim. Returning it untouched would leave
        # whatever a previous run wrote — and these are precisely the sentences
        # the sentence-initial-capital bug reached ("Seven of these 26 patients
        # had inoperable tumors" → query `Seven`), so a stale verdict here is
        # the one least worth preserving. 13 claims kept a bogus "verified"
        # this way through the first re-run.
        return _mark_unverified(claim)

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
    # A search that RETURNS something is not a search that FOUND something.
    # Accepting any non-empty result marked 44 claims "verified" against papers
    # sharing no subject matter at all -- an electric-fields brain-cancer claim
    # against freshwater fish biodiversity and speech-language pathology --
    # because a degenerate query returns the most recently indexed records.
    # See analysis/news-verification-audit.md.
    supporting = [h for h in pubmed_hits
                  if h.get("pmid") and supports_claim(claim.get("text", ""), h.get("title", ""))]
    if supporting:
        claim["verification_status"] = "verified"
        claim["verification_source"] = "pubmed"
        claim["linked_pmids"] = [h["pmid"] for h in supporting]
        return claim

    # No match
    return _mark_unverified(claim)


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
    source_domain = fm.get("source_domain", "")
    changed = 0

    for claim in claims:
        old_status = claim.get("verification_status")
        verify_claim(claim, corpus_index, source_domain=source_domain)
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
                        help="Verify articles that still have unverified factual claims")
    parser.add_argument("--force", action="store_true",
                        help="Re-verify EVERY article, including ones whose factual "
                             "claims are all already verified (use after a linker fix)")
    args = parser.parse_args()

    if not args.article and not (args.all or args.force):
        parser.error("Provide an article path, --all, or --force")

    if args.article:
        path = Path(args.article).resolve()
        changed = verify_article(path)
        print(f"Verified {path.name}: {changed} claim(s) updated")
        return

    articles = find_all_articles()
    pending = [a for a in articles if _has_unverified_factual(a)]

    if args.force:
        targets = articles
        print(f"Re-verifying all {len(targets)} articles (--force)")
    else:
        # --all skips any article whose factual claims are ALL already verified,
        # which is exactly the wrong set to skip after fixing the linker: those
        # articles hold the verdicts the broken linker was most confident about.
        # State the size of the blind spot rather than leaving it silent.
        targets = pending
        skipped_verified = 0
        for a in articles:
            if a in pending:
                continue
            fm, _ = load_article(a)
            skipped_verified += sum(
                1 for c in fm.get("claims", [])
                if c.get("category") == "FACTUAL"
                and c.get("verification_status") == "verified"
            )
        print(f"Articles with unverified factual claims: {len(pending)}/{len(articles)}")
        if skipped_verified:
            print(f"  note: {len(articles) - len(pending)} article(s) skipped, holding "
                  f"{skipped_verified} already-`verified` claim(s) that will NOT be "
                  f"re-checked. Use --force to revisit them.")

    total_changed = 0
    for article_path in tqdm(targets, desc="  Verifying"):
        changed = verify_article(article_path)
        total_changed += changed

    print(f"\nDone. Updated {total_changed} claims across {len(targets)} articles.")


if __name__ == "__main__":
    main()
