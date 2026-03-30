#!/usr/bin/env python3
"""Search local textbook PDFs for grounding snippets.

Usage:
    python scripts/search_books.py "oxidative phosphorylation"
    python scripts/search_books.py "apoptosis" --book Biology2e-WEB.pdf --limit 5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import fitz


BOOKS_DIR = Path(__file__).resolve().parent.parent / "books"


def snippet_around(text: str, query: str, radius: int = 220) -> str:
    lower = text.lower()
    idx = lower.find(query.lower())
    if idx == -1:
        return ""
    start = max(0, idx - radius)
    end = min(len(text), idx + len(query) + radius)
    snippet = text[start:end].replace("\n", " ")
    return " ".join(snippet.split())


def search_pdf(pdf_path: Path, query: str, limit: int) -> list[dict]:
    matches = []
    with fitz.open(pdf_path) as doc:
        for page_idx, page in enumerate(doc):
            text = page.get_text("text")
            if query.lower() not in text.lower():
                continue
            matches.append({
                "book": pdf_path.name,
                "page": page_idx + 1,
                "snippet": snippet_around(text, query),
            })
            if len(matches) >= limit:
                break
    return matches


def main() -> None:
    parser = argparse.ArgumentParser(description="Search local textbook PDFs for a query.")
    parser.add_argument("query", help="Case-insensitive text query")
    parser.add_argument("--book", help="Optional single PDF filename in books/")
    parser.add_argument("--limit", type=int, default=10, help="Max hits per book")
    args = parser.parse_args()

    if args.book:
        pdfs = [BOOKS_DIR / args.book]
    else:
        pdfs = sorted(BOOKS_DIR.glob("*.pdf"))

    total = 0
    for pdf in pdfs:
        if not pdf.exists():
            continue
        matches = search_pdf(pdf, args.query, args.limit)
        if not matches:
            continue
        print(f"\n== {pdf.name} ==")
        for match in matches:
            print(f"p.{match['page']}: {match['snippet']}")
            total += 1

    if total == 0:
        print("No matches found.")


if __name__ == "__main__":
    main()
