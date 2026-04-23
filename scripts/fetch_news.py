#!/usr/bin/env python3
"""Fetch and store news articles for the authentication pipeline.

Usage:
    python fetch_news.py --url URL           # Fetch single article
    python fetch_news.py --rss FEED_URL      # Fetch from RSS feed
    python fetch_news.py --rss FEED_URL --limit 10
"""

import argparse
import hashlib
import re
import sys
from pathlib import Path
from urllib.parse import urlparse

import yaml
from bs4 import BeautifulSoup
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
# Source classification
# ---------------------------------------------------------------------------

def classify_source(url: str) -> tuple[int, str, str]:
    """Match URL against SOURCE_TIER_DEFINITIONS path_prefixes.

    Matches the LONGEST prefix to avoid ambiguity when a domain appears in
    multiple tiers with different path qualifiers.

    Returns:
        (tier_num, tier_name, source_domain) -- (0, "Excluded", domain) if
        no prefix matches.
    """
    parsed = urlparse(url)
    domain = parsed.netloc.removeprefix("www.")
    # Build a normalised URL fragment: "domain/path" without scheme or www.
    url_path = f"{domain}{parsed.path}".rstrip("/")

    best_tier: int = 0
    best_name: str = "Excluded"
    best_len: int = 0

    for tier_key, tier_def in SOURCE_TIER_DEFINITIONS.items():
        tier_num = int(tier_key[-1])  # "tier1" -> 1
        tier_name = tier_def["name"]
        for prefix in tier_def["path_prefixes"]:
            prefix_norm = prefix.removeprefix("www.").rstrip("/")
            if url_path.startswith(prefix_norm) and len(prefix_norm) > best_len:
                best_tier = tier_num
                best_name = tier_name
                best_len = len(prefix_norm)

    return best_tier, best_name, domain


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

def fetch_url(url: str) -> tuple[str, str, bool]:
    """Fetch a news URL using the shared rate limiter.

    Returns:
        (html, final_url, is_paywall) -- is_paywall is True when the
        response indicates a paywall (HTTP 402 or body < 500 chars).
    """
    resp = resilient_get(url, rate_limiter=NEWS_RATE)
    if resp.status_code >= 400 and resp.status_code != 402:
        raise RuntimeError(f"HTTP {resp.status_code} fetching {url}")
    final_url = resp.url
    is_paywall = False

    if resp.status_code == 402:
        is_paywall = True
        return resp.text, final_url, is_paywall

    html_text = resp.text
    if len(html_text) < 500:
        is_paywall = True

    return html_text, final_url, is_paywall


# ---------------------------------------------------------------------------
# Text extraction
# ---------------------------------------------------------------------------

def extract_text(html: str) -> tuple[str, str, str, str]:
    """Extract article metadata and body from raw HTML.

    Returns:
        (title, author, date, body_text)
    """
    soup = BeautifulSoup(html, "html.parser")

    # --- Title ---
    title = ""
    h1 = soup.find("h1")
    if h1:
        title = h1.get_text(strip=True)
    if not title:
        title_tag = soup.find("title")
        if title_tag:
            title = title_tag.get_text(strip=True)

    # --- Author ---
    author = ""
    meta_author = soup.find("meta", attrs={"name": "author"})
    if meta_author and meta_author.get("content"):
        author = meta_author["content"].strip()
    if not author:
        # Look for common byline patterns
        for cls in ("byline", "author", "author-name"):
            el = soup.find(class_=re.compile(cls, re.IGNORECASE))
            if el:
                author = el.get_text(strip=True)
                break

    # --- Date ---
    date_str = ""
    time_el = soup.find("time")
    if time_el:
        date_str = time_el.get("datetime", "") or time_el.get_text(strip=True)
    if not date_str:
        meta_date = soup.find("meta", attrs={"property": "article:published_time"})
        if meta_date and meta_date.get("content"):
            date_str = meta_date["content"].strip()
    # Normalise to YYYY-MM-DD if possible
    date_match = re.search(r"(\d{4}-\d{2}-\d{2})", date_str)
    if date_match:
        date_str = date_match.group(1)

    # --- Body ---
    body_text = ""
    article_el = soup.find("article")
    if article_el:
        body_text = article_el.get_text(separator="\n", strip=True)
    if not body_text:
        main_el = soup.find("main")
        if main_el:
            body_text = main_el.get_text(separator="\n", strip=True)
    if not body_text:
        # Fallback: largest text block among <div>s
        divs = soup.find_all("div")
        if divs:
            body_text = max(
                (d.get_text(separator="\n", strip=True) for d in divs),
                key=len,
                default="",
            )

    return title, author, date_str, body_text


# ---------------------------------------------------------------------------
# Hashing & slugifying
# ---------------------------------------------------------------------------

# Phrases that indicate page chrome / navigation, not article content.
_BOILERPLATE_PATTERNS = [
    "Explore More", "RELATED STORIES", "RELATED TOPICS",
    "READ MORE", "ADVERTISEMENT", "Subscribe to",
    "Share this story", "Share on Facebook", "Share on Twitter",
    "Sign up for", "Newsletter", "Terms of Use",
    "Privacy Policy", "Cookie Policy", "ScienceDaily.",
    "About ScienceDaily", "Free Newsletters",
    "Materials provided by", "Note: Content may be edited",
    "MOST POPULAR THIS WEEK", "Strange & Offbeat",
]


def _strip_boilerplate(text: str) -> str:
    """Remove lines containing known page-chrome phrases."""
    lines = text.split("\n")
    cleaned = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        # Skip lines that are entirely boilerplate
        if any(bp.lower() in stripped.lower() for bp in _BOILERPLATE_PATTERNS):
            continue
        # Skip very short lines (navigation crumbs, labels)
        if len(stripped) < 20 and not any(c.isdigit() for c in stripped):
            continue
        cleaned.append(line)
    return "\n".join(cleaned)


def compute_content_hash(text: str, url: str = "") -> str:
    """SHA-256 of whitespace-normalised text, prefixed with ``sha256:``.

    When *text* is empty or whitespace-only the *url* is included in the
    hash so that distinct empty-body articles do not collide.
    """
    normalised = " ".join(text.split())
    if not normalised:
        normalised = url + normalised
    digest = hashlib.sha256(normalised.encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


def slugify(title: str) -> str:
    """Lowercase, replace non-alphanum with hyphens, truncate to 60 chars."""
    slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
    return slug[:60]


# ---------------------------------------------------------------------------
# Storage
# ---------------------------------------------------------------------------

def _find_existing(canonical_url: str, domain: str) -> Path | None:
    """Return the first article whose canonical_url matches, or None."""
    source_dir = NEWS_DIR / "by-source" / domain
    if not source_dir.exists():
        return None
    for md_path in source_dir.glob("*.md"):
        fm, _ = load_article(md_path)
        if fm.get("canonical_url") == canonical_url:
            return md_path
    return None


def store_article(
    url: str,
    canonical_url: str,
    tier: int,
    tier_name: str,
    domain: str,
    title: str,
    author: str,
    date: str,
    text: str,
    paywall: bool,
    content_hash: str,
) -> Path | None:
    """Write a fetched article to ``news/by-source/{domain}/``.

    Handles deduplication via canonical_url + content_hash.  When the
    content hash changes, a version-2 file is created and both files
    cross-reference each other.

    Returns:
        The path written, or None if the article was skipped.
    """
    existing = _find_existing(canonical_url, domain)

    if existing is not None:
        fm, body = load_article(existing)
        if fm.get("content_hash") == content_hash:
            print(f"  already fetched: {existing.name}")
            return None
        # Content changed -- create a v2 file, link both
        slug = slugify(title)
        date_prefix = date or "undated"
        v2_name = f"{date_prefix}-{slug}-v2.md"
        dest_dir = NEWS_DIR / "by-source" / domain
        dest_dir.mkdir(parents=True, exist_ok=True)
        v2_path = dest_dir / v2_name

        # Update old file
        fm["superseded_by"] = v2_name
        save_article(existing, fm, body)

        frontmatter = _build_frontmatter(
            url, canonical_url, tier, tier_name, domain,
            title, author, date, paywall, content_hash,
        )
        frontmatter["supersedes"] = existing.name
        save_article(v2_path, frontmatter, text)
        print(f"  updated version: {v2_path}")
        return v2_path

    # New article
    slug = slugify(title)
    date_prefix = date or "undated"
    filename = f"{date_prefix}-{slug}.md"
    dest_dir = NEWS_DIR / "by-source" / domain
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / filename

    frontmatter = _build_frontmatter(
        url, canonical_url, tier, tier_name, domain,
        title, author, date, paywall, content_hash,
    )
    save_article(dest, frontmatter, text)
    print(f"  stored: {dest}")
    return dest


def _build_frontmatter(
    url: str,
    canonical_url: str,
    tier: int,
    tier_name: str,
    domain: str,
    title: str,
    author: str,
    date: str,
    paywall: bool,
    content_hash: str,
) -> dict:
    """Construct the YAML frontmatter dict for a news article."""
    fm: dict = {
        "url": url,
        "canonical_url": canonical_url,
        "source_domain": domain,
        "tier": tier,
        "tier_name": tier_name,
        "title": title,
        "date_published": date or None,
        "content_hash": content_hash,
        "paywall": paywall,
        "claims": [],
        "review_status": "pending",
        "credibility_score": None,
    }
    if author:
        fm["author"] = author
    return fm


# ---------------------------------------------------------------------------
# RSS helpers
# ---------------------------------------------------------------------------

def fetch_rss_entries(feed_url: str, limit: int | None = None) -> list[dict]:
    """Parse an RSS/Atom feed and return a list of entry dicts.

    Each dict has at minimum ``link`` (str).  Optional: ``title``, ``published``.
    """
    import feedparser  # lazy import -- not every invocation uses RSS

    feed = feedparser.parse(feed_url)
    entries: list[dict] = []
    for entry in feed.entries:
        link = entry.get("link")
        if not link:
            continue
        # Use published_parsed (time struct) for reliable date extraction
        pub_date = ""
        parsed_time = entry.get("published_parsed")
        if parsed_time:
            pub_date = f"{parsed_time.tm_year:04d}-{parsed_time.tm_mon:02d}-{parsed_time.tm_mday:02d}"
        elif entry.get("published"):
            pub_date = entry["published"]

        entries.append({
            "link": link,
            "title": entry.get("title", ""),
            "published": pub_date,
        })
        if limit and len(entries) >= limit:
            break
    return entries


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def process_url(url: str, rss_date: str = "") -> Path | None:
    """End-to-end: fetch, extract, classify, store a single URL.

    *rss_date* is an optional fallback date from the RSS feed entry,
    used when the HTML extraction cannot find a date.
    """
    tier, tier_name, domain = classify_source(url)
    if tier == 0:
        print(f"  excluded (no tier match): {url}")
        return None

    print(f"  tier {tier} ({tier_name}) | {domain}")

    html_text, final_url, is_paywall = fetch_url(url)
    title, author, date_str, body_text = extract_text(html_text)

    # Use RSS date as fallback if HTML extraction didn't find one
    if not date_str and rss_date:
        date_match = re.search(r"(\d{4}-\d{2}-\d{2})", rss_date)
        if date_match:
            date_str = date_match.group(1)
        else:
            # Try parsing common RSS date formats (e.g., "Wed, 16 Apr 2026 ...")
            for fmt in ("%a, %d %b %Y %H:%M:%S %z", "%a, %d %b %Y %H:%M:%S %Z",
                        "%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%dT%H:%M:%SZ"):
                try:
                    from datetime import datetime
                    dt = datetime.strptime(rss_date.strip(), fmt)
                    date_str = dt.strftime("%Y-%m-%d")
                    break
                except ValueError:
                    continue

    # Strip common boilerplate from body text
    body_text = _strip_boilerplate(body_text)

    content_hash = compute_content_hash(body_text, url=url)

    return store_article(
        url=url,
        canonical_url=final_url,
        tier=tier,
        tier_name=tier_name,
        domain=domain,
        title=title,
        author=author,
        date=date_str,
        text=body_text,
        paywall=is_paywall,
        content_hash=content_hash,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch and store news articles.")
    parser.add_argument("--url", help="Fetch a single article by URL")
    parser.add_argument("--rss", help="Fetch articles from an RSS feed URL")
    parser.add_argument("--limit", type=int, default=None,
                        help="Max entries to fetch from RSS feed")
    args = parser.parse_args()

    if not args.url and not args.rss:
        parser.error("Provide --url or --rss")

    if args.url:
        print(f"Fetching: {args.url}")
        process_url(args.url)

    if args.rss:
        entries = fetch_rss_entries(args.rss, limit=args.limit)
        print(f"RSS feed: {len(entries)} entries")
        for entry in tqdm(entries, desc="  Fetching"):
            process_url(entry["link"], rss_date=entry.get("published", ""))


if __name__ == "__main__":
    main()
