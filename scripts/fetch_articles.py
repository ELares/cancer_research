#!/usr/bin/env python3
"""
Fetch cancer research articles from PubMed + PMC + OpenAlex.

Usage:
    python fetch_articles.py "cancer immunotherapy" --max 500
    python fetch_articles.py '"Tumor Treating Fields"[MeSH]' --max 1000
    python fetch_articles.py --query-file queries.txt
"""

import argparse
import json
import re
import sys
import time
from datetime import date
from pathlib import Path
from xml.etree import ElementTree as ET

import requests
import yaml
from tqdm import tqdm

from config import (
    CORPUS_DIR, DOI_LOOKUP, INDEX_FILE, NCBI_API_KEY, NCBI_RATE,
    OPENALEX_API_KEY, OPENALEX_EMAIL, OPENALEX_RATE, OPENALEX_WORKS,
    PMC_BIOC, PMC_BIOC_RATE, PMID_DIR, PUBMED_EFETCH, PUBMED_ESEARCH,
    resilient_get,
)


def pubmed_search(query: str, max_results: int = 500) -> list[str]:
    """Search PubMed and return list of PMIDs."""
    pmids = []
    retstart = 0
    batch_size = min(max_results, 10000)

    while retstart < max_results:
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": batch_size,
            "retstart": retstart,
            "retmode": "json",
            "sort": "relevance",
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        resp = resilient_get(PUBMED_ESEARCH, params=params, timeout=30, rate_limiter=NCBI_RATE)
        resp.raise_for_status()
        data = resp.json()

        result = data.get("esearchresult", {})
        batch_ids = result.get("idlist", [])
        total = int(result.get("count", 0))

        if not batch_ids:
            break

        pmids.extend(batch_ids)
        retstart += len(batch_ids)

        if retstart >= total:
            break

    return pmids[:max_results]


def fetch_pubmed_metadata(pmids: list[str]) -> list[dict]:
    """Fetch metadata + abstracts from PubMed for a batch of PMIDs."""
    articles = []
    batch_size = 200

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]

        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "rettype": "xml",
            "retmode": "xml",
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        resp = resilient_get(PUBMED_EFETCH, params=params, timeout=60, rate_limiter=NCBI_RATE)
        resp.raise_for_status()

        root = ET.fromstring(resp.content)
        for article_elem in root.findall(".//PubmedArticle"):
            article = parse_pubmed_xml(article_elem)
            if article:
                articles.append(article)

    return articles


def parse_pubmed_xml(elem) -> dict | None:
    """Parse a single PubmedArticle XML element into a dict."""
    medline = elem.find("MedlineCitation")
    if medline is None:
        return None

    pmid_elem = medline.find("PMID")
    if pmid_elem is None:
        return None

    pmid = pmid_elem.text
    article = medline.find("Article")
    if article is None:
        return None

    # Title
    title_elem = article.find("ArticleTitle")
    title = _text(title_elem) or ""

    # Abstract
    abstract_parts = []
    abstract_elem = article.find("Abstract")
    if abstract_elem is not None:
        for text_elem in abstract_elem.findall("AbstractText"):
            label = text_elem.get("Label", "")
            text = _text(text_elem) or ""
            if label and text:
                abstract_parts.append(f"**{label}**: {text}")
            elif text:
                abstract_parts.append(text)
    abstract = "\n\n".join(abstract_parts)

    # Authors
    authors = []
    for author in article.findall(".//Author"):
        last = _text(author.find("LastName")) or ""
        fore = _text(author.find("ForeName")) or ""
        if last:
            authors.append(f"{last} {fore}".strip())

    # Journal
    journal_elem = article.find("Journal")
    journal_title = ""
    year, month, volume, issue = None, None, "", ""
    if journal_elem is not None:
        journal_title = _text(journal_elem.find("Title")) or ""
        ji = journal_elem.find("JournalIssue")
        if ji is not None:
            volume = _text(ji.find("Volume")) or ""
            issue = _text(ji.find("Issue")) or ""
            pub_date = ji.find("PubDate")
            if pub_date is not None:
                y = _text(pub_date.find("Year"))
                m = _text(pub_date.find("Month"))
                if y:
                    year = int(y)
                if m:
                    month = _month_to_int(m)

    # Pages
    pages = _text(article.find("Pagination/MedlinePgn")) or ""

    # DOI
    doi = ""
    for id_elem in article.findall("ELocationID"):
        if id_elem.get("EIdType") == "doi":
            doi = id_elem.text or ""
            break
    if not doi:
        article_ids = elem.find("PubmedData/ArticleIdList")
        if article_ids is not None:
            for aid in article_ids.findall("ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = aid.text or ""
                    break

    # PMC ID
    pmcid = ""
    article_ids = elem.find("PubmedData/ArticleIdList")
    if article_ids is not None:
        for aid in article_ids.findall("ArticleId"):
            if aid.get("IdType") == "pmc":
                pmcid = aid.text or ""
                break

    # MeSH terms
    mesh_terms = []
    for mesh in medline.findall(".//MeshHeading/DescriptorName"):
        term = mesh.text
        if term:
            mesh_terms.append(term)

    # Publication types
    pub_types = []
    for pt in article.findall(".//PublicationType"):
        if pt.text:
            pub_types.append(pt.text)

    return {
        "pmid": pmid,
        "doi": doi,
        "pmcid": pmcid,
        "title": title,
        "authors": authors,
        "journal": journal_title,
        "year": year,
        "month": month,
        "volume": volume,
        "issue": issue,
        "pages": pages,
        "abstract": abstract,
        "mesh_terms": mesh_terms,
        "pub_types": pub_types,
    }


def enrich_with_openalex(articles: list[dict]) -> None:
    """Add OA status and PDF URLs from OpenAlex (in-place)."""
    # Build DOI → article index
    doi_map = {}
    for art in articles:
        if art.get("doi"):
            doi_map[art["doi"].lower()] = art

    if not doi_map:
        return

    # Query OpenAlex in batches using DOI filter
    dois = list(doi_map.keys())
    batch_size = 50

    for i in range(0, len(dois), batch_size):
        batch = dois[i:i + batch_size]
        doi_filter = "|".join(f"https://doi.org/{d}" for d in batch)

        params = {
            "filter": f"doi:{doi_filter}",
            "per_page": 50,
            "select": "doi,open_access,cited_by_count,primary_topic,biblio",
        }
        if OPENALEX_EMAIL:
            params["mailto"] = OPENALEX_EMAIL
        if OPENALEX_API_KEY:
            params["api_key"] = OPENALEX_API_KEY

        try:
            resp = resilient_get(OPENALEX_WORKS, params=params, timeout=30, rate_limiter=OPENALEX_RATE)
            resp.raise_for_status()
            data = resp.json()

            for work in data.get("results", []):
                work_doi = (work.get("doi") or "").replace("https://doi.org/", "").lower()
                art = doi_map.get(work_doi)
                if not art:
                    continue

                oa = work.get("open_access", {})
                art["is_oa"] = oa.get("is_oa", False)
                art["oa_status"] = oa.get("oa_status", "closed")
                art["oa_url"] = oa.get("oa_url") or ""
                art["cited_by_count"] = work.get("cited_by_count", 0)

                topic = work.get("primary_topic")
                if topic:
                    art["openalex_topic"] = topic.get("display_name", "")

        except requests.RequestException as e:
            print(f"  Warning: OpenAlex batch failed: {e}", file=sys.stderr)


def fetch_pmc_fulltext(pmcid: str) -> str | None:
    """Fetch full text from PMC BioC API. Returns plain text or None."""
    if not pmcid:
        return None

    url = f"{PMC_BIOC}/{pmcid}/unicode"

    try:
        resp = resilient_get(url, timeout=60, rate_limiter=PMC_BIOC_RATE)
        if resp.status_code != 200:
            return None

        data = resp.json()
        texts = []
        for doc in data if isinstance(data, list) else [data]:
            for document in doc.get("documents", [doc]):
                for passage in document.get("passages", []):
                    section = passage.get("infons", {}).get("section_type", "")
                    text = passage.get("text", "")
                    if text and section not in ("REF", "SUPPL", "TABLE", "FIG"):
                        texts.append(text)

        return "\n\n".join(texts) if texts else None

    except (requests.RequestException, json.JSONDecodeError, KeyError):
        return None


def save_article(article: dict, full_text: str | None) -> Path:
    """Save article as markdown with YAML frontmatter."""
    pmid = article["pmid"]
    filepath = PMID_DIR / f"{pmid}.md"

    # Build frontmatter
    frontmatter = {
        "pmid": pmid,
        "doi": article.get("doi", ""),
        "pmcid": article.get("pmcid", ""),
        "title": article.get("title", ""),
        "authors": article.get("authors", []),
        "journal": article.get("journal", ""),
        "year": article.get("year"),
        "month": article.get("month"),
        "volume": article.get("volume", ""),
        "issue": article.get("issue", ""),
        "pages": article.get("pages", ""),
        "is_oa": article.get("is_oa", False),
        "oa_status": article.get("oa_status", "unknown"),
        "oa_url": article.get("oa_url", ""),
        "cited_by_count": article.get("cited_by_count", 0),
        "mesh_terms": article.get("mesh_terms", []),
        "pub_types": article.get("pub_types", []),
        "mechanisms": [],
        "cancer_types": [],
        "evidence_level": "",
        "genes": [],
        "drugs": [],
        "date_added": str(date.today()),
    }

    if article.get("openalex_topic"):
        frontmatter["openalex_topic"] = article["openalex_topic"]

    # Build markdown body
    body_parts = []
    body_parts.append(f"# {article.get('title', 'Untitled')}\n")

    if article.get("abstract"):
        body_parts.append("## Abstract\n")
        body_parts.append(article["abstract"])
        body_parts.append("")

    if full_text:
        body_parts.append("## Full Text\n")
        body_parts.append(full_text)
        body_parts.append("")
    else:
        access = "open access" if article.get("is_oa") else "paywalled"
        body_parts.append(f"## Full Text\n\nFull text not downloaded ({access}).\n")

    # Write file
    yaml_str = yaml.dump(frontmatter, default_flow_style=False, allow_unicode=True, sort_keys=False)
    content = f"---\n{yaml_str}---\n\n" + "\n".join(body_parts)
    filepath.write_text(content, encoding="utf-8")

    return filepath


def update_doi_lookup(articles: list[dict]) -> None:
    """Append DOI→PMID mappings to DOI_LOOKUP.jsonl."""
    # Load existing mappings
    existing = set()
    if DOI_LOOKUP.exists():
        for line in DOI_LOOKUP.read_text().splitlines():
            if line.strip():
                try:
                    entry = json.loads(line)
                    existing.add(entry.get("pmid", ""))
                except json.JSONDecodeError:
                    pass

    with open(DOI_LOOKUP, "a", encoding="utf-8") as f:
        for art in articles:
            if art.get("doi") and art["pmid"] not in existing:
                entry = {"doi": art["doi"], "pmid": art["pmid"]}
                f.write(json.dumps(entry) + "\n")


def _text(elem) -> str | None:
    """Extract text content from an XML element, including mixed content."""
    if elem is None:
        return None
    return "".join(elem.itertext()).strip() or None


def _month_to_int(month_str: str) -> int | None:
    """Convert month name or number string to int."""
    if not month_str:
        return None
    if month_str.isdigit():
        return int(month_str)
    months = {
        "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
        "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
    }
    return months.get(month_str[:3].lower())


def main():
    parser = argparse.ArgumentParser(description="Fetch cancer research articles from PubMed")
    parser.add_argument("query", nargs="?", help="PubMed search query (MeSH or keyword)")
    parser.add_argument("--query-file", help="File with one query per line")
    parser.add_argument("--max", type=int, default=500, help="Max articles per query (default: 500)")
    parser.add_argument("--skip-fulltext", action="store_true", help="Skip PMC full-text download")
    parser.add_argument("--skip-openalex", action="store_true", help="Skip OpenAlex OA enrichment")
    args = parser.parse_args()

    # Collect queries
    queries = []
    if args.query:
        queries.append(args.query)
    if args.query_file:
        qf = Path(args.query_file)
        if qf.exists():
            queries.extend(line.strip() for line in qf.read_text().splitlines() if line.strip() and not line.startswith("#"))

    if not queries:
        parser.error("Provide a query or --query-file")

    # Ensure directories exist
    PMID_DIR.mkdir(parents=True, exist_ok=True)
    DOI_LOOKUP.parent.mkdir(parents=True, exist_ok=True)

    # Check which PMIDs we already have
    existing_pmids = {p.stem for p in PMID_DIR.glob("*.md")}
    print(f"Existing articles in corpus: {len(existing_pmids)}")

    total_new = 0

    for query in queries:
        print(f"\n{'='*60}")
        print(f"Query: {query}")
        print(f"{'='*60}")

        # Step 1: Search PubMed
        print(f"Searching PubMed (max {args.max})...")
        pmids = pubmed_search(query, max_results=args.max)
        print(f"  Found {len(pmids)} PMIDs")

        # Filter out already-fetched
        new_pmids = [p for p in pmids if p not in existing_pmids]
        print(f"  New (not in corpus): {len(new_pmids)}")

        if not new_pmids:
            print("  All articles already fetched. Skipping.")
            continue

        # Step 2: Fetch metadata + abstracts
        print(f"Fetching metadata from PubMed...")
        articles = fetch_pubmed_metadata(new_pmids)
        print(f"  Got metadata for {len(articles)} articles")

        # Step 3: Enrich with OpenAlex
        if not args.skip_openalex:
            print("Enriching with OpenAlex (OA status, citations)...")
            enrich_with_openalex(articles)
            oa_count = sum(1 for a in articles if a.get("is_oa"))
            print(f"  Open access: {oa_count}/{len(articles)}")

        # Step 4: Download full text for OA articles with PMC IDs
        if not args.skip_fulltext:
            oa_with_pmc = [a for a in articles if a.get("pmcid")]
            print(f"Downloading full text from PMC ({len(oa_with_pmc)} articles with PMCID)...")
            for art in tqdm(oa_with_pmc, desc="  PMC full text"):
                full_text = fetch_pmc_fulltext(art["pmcid"])
                art["_full_text"] = full_text

        # Step 5: Save articles
        print("Saving articles...")
        saved = 0
        for art in tqdm(articles, desc="  Saving"):
            full_text = art.pop("_full_text", None)
            save_article(art, full_text)
            existing_pmids.add(art["pmid"])
            saved += 1

        # Step 6: Update DOI lookup
        update_doi_lookup(articles)

        print(f"  Saved {saved} new articles")
        total_new += saved

    print(f"\n{'='*60}")
    print(f"Done. Total new articles: {total_new}")
    print(f"Total in corpus: {len(existing_pmids)}")
    print(f"\nNext steps:")
    print(f"  python enrich_metadata.py   # Add PubTator + iCite data")
    print(f"  python tag_articles.py      # Auto-tag mechanisms + cancer types")
    print(f"  python build_index.py       # Rebuild INDEX.jsonl")


if __name__ == "__main__":
    main()
