#!/usr/bin/env python3
"""
Enrich existing corpus articles with PubTator3 annotations and NIH iCite metrics.

Usage:
    python enrich_metadata.py              # Enrich all un-enriched articles
    python enrich_metadata.py --force      # Re-enrich all articles
    python enrich_metadata.py --pmids 123 456  # Enrich specific PMIDs
"""

import argparse
import json
import sys
from pathlib import Path

import requests
from tqdm import tqdm

from article_io import load_article, save_article
from config import (
    ICITE_API, ICITE_RATE, PMID_DIR, PUBTATOR_API, PUBTATOR_RATE,
    resilient_get,
)


def fetch_pubtator_batch(pmids: list[str]) -> dict:
    """Fetch PubTator3 annotations for a batch of PMIDs. Returns {pmid: {genes, diseases, chemicals, mutations}}."""
    results = {}
    batch_size = 100

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]

        try:
            resp = resilient_get(
                PUBTATOR_API,
                params={"pmids": ",".join(batch)},
                timeout=30,
                rate_limiter=PUBTATOR_RATE,
            )
            if resp.status_code != 200:
                continue

            data = resp.json()
            entries = data.get("PubTator3", data) if isinstance(data, dict) else data

            for entry in entries:
                pmid = str(entry.get("pmid", entry.get("id", ""))).split("|")[0]
                if not pmid:
                    continue

                annotations = {"genes": set(), "diseases": set(), "chemicals": set(), "mutations": set(), "species": set()}

                for passage in entry.get("passages", []):
                    for ann in passage.get("annotations", []):
                        ann_type = ann.get("infons", {}).get("type", "").lower()
                        text = ann.get("text", "")
                        if not text:
                            continue
                        if ann_type == "gene":
                            annotations["genes"].add(text)
                        elif ann_type == "disease":
                            annotations["diseases"].add(text)
                        elif ann_type == "chemical":
                            annotations["chemicals"].add(text)
                        elif ann_type in ("mutation", "variant", "snp"):
                            annotations["mutations"].add(text)
                        elif ann_type == "species":
                            annotations["species"].add(text)

                results[pmid] = {k: sorted(v) for k, v in annotations.items()}

        except (requests.RequestException, json.JSONDecodeError) as e:
            print(f"  Warning: PubTator batch failed: {e}", file=sys.stderr)

    return results


def fetch_icite_batch(pmids: list[str]) -> dict:
    """Fetch iCite citation metrics for a batch of PMIDs. Returns {pmid: metrics_dict}."""
    results = {}
    batch_size = 200

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]

        try:
            resp = resilient_get(
                ICITE_API,
                params={"pmids": ",".join(batch)},
                timeout=30,
                rate_limiter=ICITE_RATE,
            )
            if resp.status_code != 200:
                continue

            data = resp.json()
            for pub in data.get("data", []):
                pmid = str(pub.get("pmid", ""))
                if not pmid:
                    continue

                results[pmid] = {
                    "citation_count": pub.get("citation_count", 0),
                    "relative_citation_ratio": pub.get("relative_citation_ratio"),
                    "nih_percentile": pub.get("nih_percentile"),
                    "expected_citations_per_year": pub.get("expected_citations_per_year"),
                    "is_clinical": pub.get("is_clinical", False),
                    "apt": pub.get("apt"),  # Approximate Potential to Translate
                }

        except (requests.RequestException, json.JSONDecodeError) as e:
            print(f"  Warning: iCite batch failed: {e}", file=sys.stderr)

    return results


def main():
    parser = argparse.ArgumentParser(description="Enrich corpus with PubTator + iCite data")
    parser.add_argument("--force", action="store_true", help="Re-enrich all articles")
    parser.add_argument("--pmids", nargs="+", help="Only enrich specific PMIDs")
    parser.add_argument("--skip-pubtator", action="store_true", help="Skip PubTator annotations")
    parser.add_argument("--skip-icite", action="store_true", help="Skip iCite metrics")
    args = parser.parse_args()

    # Find articles to enrich
    if args.pmids:
        files = [PMID_DIR / f"{pmid}.md" for pmid in args.pmids]
        files = [f for f in files if f.exists()]
    else:
        files = sorted(PMID_DIR.glob("*.md"))

    if not files:
        print("No articles found to enrich.")
        return

    # Filter to un-enriched unless --force
    if not args.force:
        unenriched = []
        for f in files:
            fm, _ = load_article(f)
            needs_pubtator = not args.skip_pubtator and not fm.get("pubtator_enriched")
            needs_icite = not args.skip_icite and not fm.get("icite_enriched")
            if needs_pubtator or needs_icite:
                unenriched.append(f)
        files = unenriched

    print(f"Articles to enrich: {len(files)}")
    if not files:
        print("All articles already enriched. Use --force to re-enrich.")
        return

    pmids = [f.stem for f in files]

    # Step 1: PubTator
    pubtator_data = {}
    if not args.skip_pubtator:
        print("Fetching PubTator3 annotations...")
        pubtator_data = fetch_pubtator_batch(pmids)
        print(f"  Got annotations for {len(pubtator_data)} articles")

    # Step 2: iCite
    icite_data = {}
    if not args.skip_icite:
        print("Fetching NIH iCite metrics...")
        icite_data = fetch_icite_batch(pmids)
        print(f"  Got metrics for {len(icite_data)} articles")

    # Step 3: Update articles
    print("Updating articles...")
    updated = 0
    for filepath in tqdm(files, desc="  Enriching"):
        pmid = filepath.stem
        fm, body = load_article(filepath)
        if not fm:
            continue

        changed = False

        # PubTator
        pt = pubtator_data.get(pmid)
        if pt:
            fm["genes"] = pt.get("genes", [])
            fm["drugs"] = pt.get("chemicals", [])
            fm["diseases_annotated"] = pt.get("diseases", [])
            fm["mutations"] = pt.get("mutations", [])
            fm["pubtator_enriched"] = True
            changed = True
        elif not args.skip_pubtator:
            fm["pubtator_enriched"] = True
            changed = True

        # iCite
        ic = icite_data.get(pmid)
        if ic:
            fm["icite_citation_count"] = ic.get("citation_count", 0)
            fm["icite_rcr"] = ic.get("relative_citation_ratio")
            fm["icite_percentile"] = ic.get("nih_percentile")
            fm["icite_is_clinical"] = ic.get("is_clinical", False)
            fm["icite_apt"] = ic.get("apt")
            fm["icite_enriched"] = True
            changed = True
        elif not args.skip_icite:
            fm["icite_enriched"] = True
            changed = True

        if changed:
            save_article(filepath, fm, body)
            updated += 1

    print(f"\nDone. Updated {updated} articles.")
    print(f"\nNext steps:")
    print(f"  python tag_articles.py   # Auto-tag mechanisms + cancer types")
    print(f"  python build_index.py    # Rebuild INDEX.jsonl")


if __name__ == "__main__":
    main()
