#!/usr/bin/env python3
"""
Fetch additional cancer research articles from Semantic Scholar.

Supplements the PubMed corpus with:
1. Citation-graph discovery: find highly-cited papers our PubMed queries missed
2. TLDR summaries for existing corpus articles
3. Broader search coverage (S2 indexes more sources than PubMed alone)

Usage:
    python fetch_semantic_scholar.py --mode search --max 500
    python fetch_semantic_scholar.py --mode enrich-tldr
    python fetch_semantic_scholar.py --mode citation-discovery --seed-pmids 29260225 37165196
"""

import argparse
import json
import re
import sys
import time
from datetime import date
from pathlib import Path

import requests
import yaml
from tqdm import tqdm

from config import (
    CORPUS_DIR, DOI_LOOKUP, NCBI_API_KEY, NCBI_RATE,
    OPENALEX_API_KEY, OPENALEX_EMAIL, OPENALEX_RATE, OPENALEX_WORKS,
    PMC_BIOC, PMC_BIOC_RATE, PMID_DIR, PROJECT_ROOT, PUBMED_EFETCH,
    SEMANTIC_SCHOLAR_API_KEY, resilient_get,
)
from article_io import load_article, load_frontmatter, save_article

S2_API = "https://api.semanticscholar.org/graph/v1"
S2_FIELDS = "title,abstract,year,citationCount,isOpenAccess,openAccessPdf,tldr,externalIds,publicationTypes"
S2_HEADERS = {"x-api-key": SEMANTIC_SCHOLAR_API_KEY} if SEMANTIC_SCHOLAR_API_KEY else {}

# Rate limit: 10 req/s with key
S2_INTERVAL = 0.12 if SEMANTIC_SCHOLAR_API_KEY else 1.1


def s2_get(url: str, params: dict = None) -> dict | None:
    """GET from Semantic Scholar API with rate limiting and retry."""
    time.sleep(S2_INTERVAL)
    try:
        resp = requests.get(url, params=params, headers=S2_HEADERS, timeout=30)
        if resp.status_code == 429:
            print("  Rate limited, waiting 5s...", file=sys.stderr)
            time.sleep(5)
            resp = requests.get(url, params=params, headers=S2_HEADERS, timeout=30)
        if resp.status_code != 200:
            return None
        return resp.json()
    except (requests.RequestException, json.JSONDecodeError):
        return None


# ============================================================
# Mode 1: Search — find papers S2 has that PubMed may not
# ============================================================

SEARCH_QUERIES = [
    "tumor treating fields cancer immunogenic cell death",
    "sonodynamic therapy cancer clinical trial",
    "irreversible electroporation immunotherapy synergy",
    "bioelectric membrane potential cancer therapy",
    "cancer convergent therapy multi-modal combination",
    "physical modality immunotherapy cancer synergy",
    "nanoparticle sonodynamic therapy cancer",
    "CRISPR CAR-T engineering cancer",
    "mRNA cancer vaccine neoantigen clinical",
    "oncolytic virus checkpoint inhibitor combination",
    "PARP inhibitor myeloma synthetic lethality",
    "frequency therapy pulsed electromagnetic cancer",
    "HIFU immunogenic cell death tumor ablation",
    "antibody drug conjugate bispecific combination",
    "microbiome fecal transplant immunotherapy cancer",
    "electrochemical therapy tumor ablation clinical",
    "epigenetic therapy immune checkpoint cancer",
    "metabolic reprogramming Warburg immunotherapy",
]


def search_s2(query: str, limit: int = 100) -> list[dict]:
    """Search Semantic Scholar and return papers with PMIDs."""
    papers = []
    offset = 0

    while offset < limit:
        batch_size = min(100, limit - offset)
        data = s2_get(
            f"{S2_API}/paper/search",
            params={"query": query, "limit": batch_size, "offset": offset, "fields": S2_FIELDS},
        )
        if not data or not data.get("data"):
            break

        for paper in data["data"]:
            ext = paper.get("externalIds") or {}
            pmid = ext.get("PubMed")
            doi = ext.get("DOI", "")
            if pmid:
                papers.append({
                    "pmid": str(pmid),
                    "doi": doi,
                    "s2_id": paper.get("paperId", ""),
                    "title": paper.get("title", ""),
                    "abstract": paper.get("abstract", ""),
                    "year": paper.get("year"),
                    "citation_count": paper.get("citationCount", 0),
                    "is_oa": paper.get("isOpenAccess", False),
                    "tldr": (paper.get("tldr") or {}).get("text", ""),
                    "pub_types": paper.get("publicationTypes") or [],
                })

        offset += batch_size
        if offset >= data.get("total", 0):
            break

    return papers


# ============================================================
# Mode 2: Enrich existing corpus with TLDR summaries
# ============================================================

def enrich_tldr(batch_size: int = 500) -> None:
    """Add Semantic Scholar TLDR summaries to existing corpus articles."""
    files = sorted(PMID_DIR.glob("*.md"))
    print(f"Checking {len(files)} articles for TLDR enrichment...")

    # Find articles missing TLDR
    needs_tldr = []
    for fp in files:
        fm = load_frontmatter(fp)
        if fm and not fm.get("s2_tldr"):
            needs_tldr.append((fp, fm.get("pmid", fp.stem)))

    print(f"  {len(needs_tldr)} articles need TLDR enrichment")
    if not needs_tldr:
        return

    # Batch lookup by PMID
    updated = 0
    for i in tqdm(range(0, len(needs_tldr), batch_size), desc="  S2 batch lookup"):
        batch = needs_tldr[i:i + batch_size]
        pmid_list = [f"PMID:{pmid}" for _, pmid in batch]

        data = s2_get(
            f"{S2_API}/paper/batch",
            params={"fields": "tldr,citationCount,externalIds"},
        )
        # Batch endpoint uses POST
        time.sleep(S2_INTERVAL)
        try:
            resp = requests.post(
                f"{S2_API}/paper/batch",
                json={"ids": pmid_list},
                params={"fields": "tldr,citationCount,externalIds"},
                headers=S2_HEADERS,
                timeout=60,
            )
            if resp.status_code != 200:
                print(f"  Batch failed: HTTP {resp.status_code}", file=sys.stderr)
                continue
            results = resp.json()
        except (requests.RequestException, json.JSONDecodeError) as e:
            print(f"  Batch failed: {e}", file=sys.stderr)
            continue

        for (fp, pmid), result in zip(batch, results):
            if result is None:
                continue
            tldr = (result.get("tldr") or {}).get("text", "")
            s2_citations = result.get("citationCount", 0)

            if tldr:
                fm, body = load_article(fp)
                fm["s2_tldr"] = tldr
                fm["s2_citation_count"] = s2_citations
                save_article(fp, fm, body)
                updated += 1

    print(f"\nDone. Added TLDR to {updated} articles.")


# ============================================================
# Mode 3: Citation discovery — find papers citing our top articles
# ============================================================

def citation_discovery(seed_pmids: list[str], max_per_seed: int = 100) -> list[dict]:
    """Find highly-cited papers that cite our seed articles."""
    discovered = []
    existing = {p.stem for p in PMID_DIR.glob("*.md")}

    for pmid in tqdm(seed_pmids, desc="  Discovering citations"):
        data = s2_get(
            f"{S2_API}/paper/PMID:{pmid}/citations",
            params={"fields": S2_FIELDS, "limit": max_per_seed},
        )
        if not data or not data.get("data"):
            continue

        for entry in data["data"]:
            citing = entry.get("citingPaper", {})
            ext = citing.get("externalIds") or {}
            citing_pmid = ext.get("PubMed")
            if not citing_pmid or str(citing_pmid) in existing:
                continue

            cites = citing.get("citationCount", 0)
            if cites < 10:  # Only care about impactful citing papers
                continue

            discovered.append({
                "pmid": str(citing_pmid),
                "doi": ext.get("DOI", ""),
                "s2_id": citing.get("paperId", ""),
                "title": citing.get("title", ""),
                "abstract": citing.get("abstract", ""),
                "year": citing.get("year"),
                "citation_count": cites,
                "is_oa": citing.get("isOpenAccess", False),
                "tldr": (citing.get("tldr") or {}).get("text", ""),
                "pub_types": citing.get("publicationTypes") or [],
                "discovered_via": f"cites PMID:{pmid}",
            })
            existing.add(str(citing_pmid))

    # Deduplicate and sort by citation count
    seen = set()
    unique = []
    for p in sorted(discovered, key=lambda x: -(x.get("citation_count", 0))):
        if p["pmid"] not in seen:
            seen.add(p["pmid"])
            unique.append(p)

    return unique


def fetch_and_save_new_articles(papers: list[dict]) -> int:
    """Fetch PubMed metadata for discovered papers and save to corpus."""
    from fetch_articles import fetch_best_fulltext, fetch_pubmed_metadata, enrich_with_openalex, save_article as save_new, update_doi_lookup

    existing = {p.stem for p in PMID_DIR.glob("*.md")}
    new_pmids = [p["pmid"] for p in papers if p["pmid"] not in existing]

    if not new_pmids:
        print("No new articles to add.")
        return 0

    print(f"Fetching PubMed metadata for {len(new_pmids)} new articles...")
    articles = fetch_pubmed_metadata(new_pmids)
    print(f"  Got metadata for {len(articles)} articles")

    # Merge S2 data (TLDR, S2 citation count)
    s2_map = {p["pmid"]: p for p in papers}
    for art in articles:
        s2 = s2_map.get(art["pmid"])
        if s2:
            art["s2_tldr"] = s2.get("tldr", "")
            art["s2_citation_count"] = s2.get("citation_count", 0)
            art["discovered_via"] = s2.get("discovered_via", "semantic_scholar_search")

    print("Enriching with OpenAlex...")
    enrich_with_openalex(articles)

    # Download full text for OA articles with PMC IDs
    oa_candidates = [a for a in articles if a.get("pmcid") or a.get("oa_url")]
    if oa_candidates:
        print(f"Downloading full text ({len(oa_candidates)} OA/PMC candidates)...")
        for art in tqdm(oa_candidates, desc="  Full text"):
            art["_full_text"] = fetch_best_fulltext(art)

    print("Saving articles...")
    saved = 0
    for art in articles:
        full_text = art.pop("_full_text", None)
        save_new(art, full_text)
        saved += 1

    update_doi_lookup(articles)
    return saved


def main():
    parser = argparse.ArgumentParser(description="Fetch articles from Semantic Scholar")
    parser.add_argument("--mode", choices=["search", "enrich-tldr", "citation-discovery"], required=True)
    parser.add_argument("--max", type=int, default=200, help="Max articles per query (search mode)")
    parser.add_argument("--seed-pmids", nargs="+", help="Seed PMIDs for citation discovery")
    parser.add_argument("--seed-top-n", type=int, default=50, help="Use top N corpus articles as seeds")
    args = parser.parse_args()

    if not SEMANTIC_SCHOLAR_API_KEY:
        print("ERROR: SEMANTIC_SCHOLAR_API_KEY not set in .env", file=sys.stderr)
        sys.exit(1)

    print(f"Semantic Scholar API key: ...{SEMANTIC_SCHOLAR_API_KEY[-6:]}")
    print(f"Mode: {args.mode}\n")

    if args.mode == "search":
        all_papers = []
        for query in SEARCH_QUERIES:
            print(f"Searching: {query}")
            papers = search_s2(query, limit=args.max)
            print(f"  Found {len(papers)} papers with PMIDs")
            all_papers.extend(papers)

        # Deduplicate
        seen = set()
        unique = []
        existing = {p.stem for p in PMID_DIR.glob("*.md")}
        for p in all_papers:
            if p["pmid"] not in seen and p["pmid"] not in existing:
                seen.add(p["pmid"])
                unique.append(p)

        print(f"\nTotal unique new papers: {len(unique)}")
        if unique:
            saved = fetch_and_save_new_articles(unique)
            print(f"Saved {saved} new articles to corpus")

    elif args.mode == "enrich-tldr":
        enrich_tldr()

    elif args.mode == "citation-discovery":
        if args.seed_pmids:
            seeds = args.seed_pmids
        else:
            # Use top N articles by citation count as seeds
            index = [json.loads(l) for l in (CORPUS_DIR / "INDEX.jsonl").read_text().splitlines() if l.strip()]
            index.sort(key=lambda x: -(x.get("cited_by_count") or 0))
            seeds = [e["pmid"] for e in index[:args.seed_top_n]]
            print(f"Using top {len(seeds)} articles as seeds")

        papers = citation_discovery(seeds, max_per_seed=200)
        print(f"\nDiscovered {len(papers)} new highly-cited papers")

        if papers:
            # Show top discoveries
            print("\nTop 20 discoveries:")
            for p in papers[:20]:
                print(f"  PMID {p['pmid']} ({p['citation_count']} cites) — {p['title'][:80]}")
                if p.get("tldr"):
                    print(f"    TLDR: {p['tldr'][:120]}")

            saved = fetch_and_save_new_articles(papers)
            print(f"\nSaved {saved} new articles to corpus")

    print("\nNext steps:")
    print("  python enrich_metadata.py    # PubTator + iCite for new articles")
    print("  python tag_articles.py       # Re-tag all articles")
    print("  python build_index.py        # Rebuild INDEX.jsonl")
    print("  python analyze_corpus.py     # Regenerate analysis files")


if __name__ == "__main__":
    main()
