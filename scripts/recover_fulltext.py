#!/usr/bin/env python3
"""Audit and recover missing article full text in the local corpus.

This script scans `corpus/by-pmid/*.md` for the standard placeholder:
`Full text not downloaded (open access|paywalled).`

It can:
- write a CSV audit of incomplete records
- retry open-access recovery from PMC BioC or the stored OA URL

Usage:
    python scripts/recover_fulltext.py --report analysis/fulltext_audit.csv
    python scripts/recover_fulltext.py --retry-oa --limit 50
"""

from __future__ import annotations

import argparse
import csv
import io
import re
from pathlib import Path

import fitz
from lxml import html

from article_io import load_article, save_article
from config import PMID_DIR, resilient_get
from fetch_articles import fetch_pmc_fulltext


PLACEHOLDER_RE = re.compile(r"Full text not downloaded \((open access|paywalled)\)\.")


def find_placeholder(body: str) -> str | None:
    match = PLACEHOLDER_RE.search(body)
    return match.group(1) if match else None


def extract_text_from_pdf_bytes(pdf_bytes: bytes) -> str | None:
    try:
        with fitz.open(stream=pdf_bytes, filetype="pdf") as doc:
            chunks = []
            for page in doc:
                text = page.get_text("text").strip()
                if text:
                    chunks.append(text)
        joined = "\n\n".join(chunks).strip()
        return joined if len(joined) >= 1000 else None
    except Exception:
        return None


def extract_text_from_html(html_bytes: bytes) -> str | None:
    try:
        doc = html.fromstring(html_bytes)
    except Exception:
        return None

    for xpath in [
        "//main//p",
        "//article//p",
        "//div[contains(@class,'article')]//p",
        "//body//p",
    ]:
        paragraphs = [p.text_content().strip() for p in doc.xpath(xpath)]
        paragraphs = [p for p in paragraphs if len(p) > 40]
        if paragraphs:
            joined = "\n\n".join(paragraphs).strip()
            return joined if len(joined) >= 1000 else None
    return None


def fetch_from_oa_url(url: str) -> tuple[str | None, str]:
    if not url:
        return None, "no_oa_url"

    try:
        resp = resilient_get(
            url,
            timeout=60,
            retries=1,
        )
    except Exception as exc:
        return None, f"request_failed:{type(exc).__name__}"

    if resp.status_code != 200:
        return None, f"http_{resp.status_code}"

    content_type = (resp.headers.get("content-type") or "").lower()
    final_url = resp.url.lower()

    if "application/pdf" in content_type or final_url.endswith(".pdf") or "pdf" in final_url:
        text = extract_text_from_pdf_bytes(resp.content)
        return text, "oa_url_pdf" if text else "oa_url_pdf_unreadable"

    text = extract_text_from_html(resp.content)
    return text, "oa_url_html" if text else "oa_url_html_unreadable"


def replace_fulltext_placeholder(body: str, full_text: str) -> str:
    replacement = f"## Full Text\n\n{full_text.strip()}\n"
    return re.sub(
        r"## Full Text\s+Full text not downloaded \((?:open access|paywalled)\)\.\s*",
        replacement,
        body,
        count=1,
        flags=re.DOTALL,
    )


def audit_rows() -> list[dict]:
    rows = []
    for path in sorted(PMID_DIR.glob("*.md")):
        fm, body = load_article(path)
        missing_kind = find_placeholder(body)
        if not missing_kind:
            continue
        rows.append({
            "pmid": fm.get("pmid", path.stem),
            "title": fm.get("title", ""),
            "year": fm.get("year", ""),
            "journal": fm.get("journal", ""),
            "missing_kind": missing_kind,
            "is_oa": fm.get("is_oa", False),
            "oa_status": fm.get("oa_status", ""),
            "pmcid": fm.get("pmcid", ""),
            "oa_url": fm.get("oa_url", ""),
            "path": str(path),
        })
    return rows


def write_csv(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "pmid", "title", "year", "journal", "missing_kind",
                "is_oa", "oa_status", "pmcid", "oa_url", "path",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def retry_open_access(limit: int | None = None) -> list[dict]:
    rows = audit_rows()
    open_access_rows = [r for r in rows if r["missing_kind"] == "open access"]
    open_access_rows.sort(
        key=lambda r: (
            0 if r["pmcid"] else 1,
            0 if r["oa_url"] else 1,
            str(r["pmid"]),
        )
    )
    if limit is not None:
        open_access_rows = open_access_rows[:limit]

    results = []
    for row in open_access_rows:
        path = Path(row["path"])
        fm, body = load_article(path)

        full_text = None
        source = ""

        pmcid = (fm.get("pmcid") or "").strip()
        if pmcid:
            full_text = fetch_pmc_fulltext(pmcid)
            source = "pmc_bioc" if full_text else "pmc_bioc_failed"

        if not full_text:
            full_text, source = fetch_from_oa_url((fm.get("oa_url") or "").strip())

        result = {
            "pmid": row["pmid"],
            "title": row["title"],
            "status": "recovered" if full_text else "still_missing",
            "source": source,
            "path": row["path"],
        }

        if full_text:
            fm["fulltext_source"] = source
            body = replace_fulltext_placeholder(body, full_text)
            save_article(path, fm, body)

        results.append(result)

    return results


def write_retry_report(results: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["pmid", "title", "status", "source", "path"],
        )
        writer.writeheader()
        writer.writerows(results)


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit/recover missing full text in the corpus.")
    parser.add_argument("--report", default="analysis/fulltext_audit.csv", help="CSV path for audit output")
    parser.add_argument("--retry-oa", action="store_true", help="Retry recovery for open-access missing full text")
    parser.add_argument("--retry-report", default="analysis/fulltext_retry_results.csv", help="CSV path for retry output")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of OA records retried")
    args = parser.parse_args()

    rows = audit_rows()
    write_csv(rows, Path(args.report))
    print(f"Wrote audit: {args.report} ({len(rows)} incomplete articles)")

    if args.retry_oa:
        results = retry_open_access(limit=args.limit)
        write_retry_report(results, Path(args.retry_report))
        recovered = sum(1 for r in results if r["status"] == "recovered")
        print(f"Retried OA recovery: {len(results)} articles, recovered {recovered}")
        print(f"Wrote retry results: {args.retry_report}")


if __name__ == "__main__":
    main()
