#!/usr/bin/env python3
"""
Verify all PMID references in the article draft against the corpus.

Checks:
- PMID exists in corpus
- Title matches (or flags mismatch)
- Journal matches
- Year matches
- First author matches

Outputs a report and optionally generates corrected reference lines.

Usage:
    python verify_references.py
"""

import re
import json
from pathlib import Path

import yaml

from config import PROJECT_ROOT

ARTICLE = PROJECT_ROOT / "article" / "drafts" / "v1.md"
PMID_DIR = PROJECT_ROOT / "corpus" / "by-pmid"


def load_frontmatter(filepath: Path) -> dict | None:
    with open(filepath, "r", encoding="utf-8") as f:
        first_line = f.readline()
        if first_line.strip() != "---":
            return None
        yaml_lines = []
        for line in f:
            if line.strip() == "---":
                break
            yaml_lines.append(line)
    if not yaml_lines:
        return None
    return yaml.safe_load("".join(yaml_lines))


def extract_references(article_text: str) -> list[dict]:
    """Extract reference list entries from the article."""
    refs = []
    # Match lines like: 1. PMID: 29260225 -- Stupp R et al. Title. *Journal* (Year). N citations.
    pattern = re.compile(
        r'^\d+\.\s+PMID:\s*(\d+)\s*--\s*(.*?)\.\s*\*(.+?)\*\s*\((\d{4})\)\.\s*([\d,]+)\s*citations\.',
        re.MULTILINE
    )
    for m in pattern.finditer(article_text):
        refs.append({
            "pmid": m.group(1),
            "text": m.group(2).strip(),
            "journal_in_article": m.group(3).strip(),
            "year_in_article": int(m.group(4)),
            "citations_in_article": m.group(5).replace(",", ""),
            "full_line": m.group(0),
        })
    return refs


def normalize_journal(j: str) -> str:
    """Normalize journal name for comparison."""
    j = j.lower().strip()
    j = re.sub(r'[.,;:()]', '', j)
    j = re.sub(r'\s+', ' ', j)
    # Common abbreviation mappings
    mappings = {
        "nat rev drug discov": "nature reviews drug discovery",
        "nat rev cancer": "nature reviews cancer",
        "nat rev clin oncol": "nature reviews clinical oncology",
        "nat med": "nature medicine",
        "nat biotechnol": "nature biotechnology",
        "nat commun": "nature communications",
        "nat chem biol": "nature chemical biology",
        "lancet oncol": "the lancet oncology",
        "j clin oncol": "journal of clinical oncology",
        "j hematol oncol": "journal of hematology & oncology",
        "ca cancer j clin": "ca a cancer journal for clinicians",
        "curr oncol": "current oncology toronto ont",
        "ann surg": "annals of surgery",
        "ann surg oncol": "annals of surgical oncology",
        "ann biomed eng": "annals of biomedical engineering",
        "ann oncol": "annals of oncology",
        "eur urol": "european urology",
        "eur urol oncol": "european urology oncology",
        "eur urol focus": "european urology focus",
        "eur j med chem": "european journal of medicinal chemistry",
        "adv mater": "advanced materials",
        "adv healthc mater": "advanced healthcare materials",
        "adv exp med biol": "advances in experimental medicine and biology",
        "adv drug deliv rev": "advanced drug delivery reviews",
        "mol cancer": "molecular cancer",
        "mol ther": "molecular therapy",
        "int j mol sci": "international journal of molecular sciences",
        "int j hyperthermia": "international journal of hyperthermia",
        "j immunother cancer": "journal for immunotherapy of cancer",
        "signal transduct target ther": "signal transduction and targeted therapy",
        "cell death dis": "cell death & disease",
        "cell death discov": "cell death discovery",
        "cell mol immunol": "cellular & molecular immunology",
        "cancer med": "cancer medicine",
        "cancer immunol res": "cancer immunology research",
        "cns oncol": "cns oncology",
        "hum gene ther": "human gene therapy",
        "gut microbes": "gut microbes",
        "trends cancer": "trends in cancer",
        "trends biochem sci": "trends in biochemical sciences",
        "curr res transl med": "current research in translational medicine",
        "blood cancer j": "blood cancer journal",
        "ultrasound med biol": "ultrasound in medicine & biology",
    }
    for abbr, full in mappings.items():
        if j == abbr:
            return full
    return j


def extract_first_author(text: str) -> str | None:
    """Extract first author from reference text like 'Stupp R et al. Title...'"""
    m = re.match(r'^([A-Z][a-z]+(?:\s+[A-Z]+)?)\s+et\s+al', text)
    if m:
        return m.group(1).strip()
    # Single author: "Levin M. Title"
    m = re.match(r'^([A-Z][a-z]+\s+[A-Z]+)\.\s', text)
    if m:
        return m.group(1).strip()
    return None


def format_author(authors: list[str]) -> str:
    """Format first author from corpus author list."""
    if not authors:
        return ""
    first = authors[0]
    parts = first.split()
    if len(parts) >= 2:
        last = parts[0]
        initials = "".join(p[0] for p in parts[1:] if p)
        return f"{last} {initials}"
    return first


def main():
    article_text = ARTICLE.read_text(encoding="utf-8")
    refs = extract_references(article_text)
    print(f"Found {len(refs)} references in article\n")

    errors = []
    corrected_lines = []

    for ref in refs:
        pmid = ref["pmid"]
        fp = PMID_DIR / f"{pmid}.md"

        if not fp.exists():
            errors.append(f"MISSING: PMID {pmid} not in corpus")
            corrected_lines.append(ref["full_line"])
            continue

        fm = load_frontmatter(fp)
        if not fm:
            errors.append(f"PARSE ERROR: PMID {pmid} frontmatter unreadable")
            corrected_lines.append(ref["full_line"])
            continue

        issues = []

        # Check journal
        corpus_journal = (fm.get("journal") or "").lower().strip()
        article_journal = normalize_journal(ref["journal_in_article"])
        corpus_journal_norm = normalize_journal(corpus_journal)

        # Fuzzy match: check if one contains the other
        journal_match = (
            article_journal == corpus_journal_norm or
            article_journal in corpus_journal_norm or
            corpus_journal_norm in article_journal or
            corpus_journal.replace("the ", "") in article_journal.replace("the ", "")
        )
        if not journal_match:
            issues.append(f"JOURNAL: article='{ref['journal_in_article']}' corpus='{fm.get('journal')}'")

        # Check year
        corpus_year = fm.get("year")
        if corpus_year and corpus_year != ref["year_in_article"]:
            issues.append(f"YEAR: article={ref['year_in_article']} corpus={corpus_year}")

        # Check first author
        article_author = extract_first_author(ref["text"])
        corpus_authors = fm.get("authors", [])
        corpus_first = format_author(corpus_authors) if corpus_authors else ""

        if article_author and corpus_first:
            # Compare last names
            article_last = article_author.split()[0].lower()
            corpus_last = corpus_first.split()[0].lower()
            if article_last != corpus_last:
                issues.append(f"AUTHOR: article='{article_author}' corpus='{corpus_first}' (full: {corpus_authors[0]})")

        # Check title (just first 50 chars)
        corpus_title = (fm.get("title") or "").strip().rstrip(".")
        # Extract title from ref text (after author, before journal)
        ref_title_match = re.search(r'(?:et al\.|[A-Z]\.) (.+)$', ref["text"])
        if ref_title_match:
            ref_title = ref_title_match.group(1).strip().rstrip(".")
            if ref_title[:40].lower() != corpus_title[:40].lower():
                issues.append(f"TITLE: article='{ref_title[:60]}...' corpus='{corpus_title[:60]}...'")

        if issues:
            errors.append(f"PMID {pmid}: " + " | ".join(issues))

        # Generate corrected line
        correct_author = corpus_first
        correct_journal = fm.get("journal", ref["journal_in_article"])
        correct_year = corpus_year or ref["year_in_article"]
        correct_title = corpus_title
        cited = fm.get("cited_by_count", ref["citations_in_article"])

        # Abbreviate journal for reference
        corrected = f"PMID: {pmid} -- {correct_author} et al. {correct_title}. *{correct_journal}* ({correct_year}). {cited} citations."
        corrected_lines.append(corrected)

    # Print report
    print(f"{'='*60}")
    print(f"VERIFICATION REPORT")
    print(f"{'='*60}")
    print(f"Total references: {len(refs)}")
    print(f"Errors found: {len(errors)}")
    print()

    if errors:
        print("ERRORS:")
        for e in errors:
            print(f"  - {e}")
    else:
        print("All references verified clean!")

    # Write corrected references
    outfile = PROJECT_ROOT / "article" / "references" / "corrected_references.txt"
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with open(outfile, "w", encoding="utf-8") as f:
        for i, line in enumerate(corrected_lines, 1):
            f.write(f"{i}. {line}\n")
    print(f"\nCorrected reference list written to: {outfile}")


if __name__ == "__main__":
    main()
