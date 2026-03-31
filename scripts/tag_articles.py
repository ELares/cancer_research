#!/usr/bin/env python3
"""
Auto-tag articles by mechanism, cancer type, and evidence level.
Reads title + abstract + MeSH terms, matches against keyword dictionaries.

Usage:
    python tag_articles.py          # Tag all articles, rebuild tag index files
    python tag_articles.py --dry-run  # Show what would be tagged without writing
"""

import argparse
import re
import sys
from pathlib import Path

from tqdm import tqdm

from article_io import load_article, save_article
from config import (
    BIOLOGY_PROCESS_KEYWORDS, CANCER_TYPE_KEYWORDS, EVIDENCE_LEVEL_KEYWORDS, MECHANISM_KEYWORDS,
    PMID_DIR, TAGS_DIR,
    RESISTANT_STATE_RULES,
)
from evidence_utils import is_protocol_like, is_review_like, normalize_text

GENERIC_CANCER_TERMS = (
    "cancer", "neoplasm", "carcinoma", "tumor", "tumour",
    "oncology", "malign", "glioblastoma", "melanoma", "sarcoma",
    "leukemia", "lymphoma", "myeloma",
)

EVIDENCE_PUBTYPE_MARKERS = {
    "phase3-clinical": ("clinical trial, phase iii", "clinical trial, phase 3", "clinical trial, phase iv", "clinical trial, phase 4"),
    "phase2-clinical": ("clinical trial, phase ii", "clinical trial, phase 2"),
    "phase1-clinical": ("clinical trial, phase i", "clinical trial, phase 1", "clinical trial, phase i/ii"),
}


def has_cancer_context(text: str) -> bool:
    """Require generic cancer context for broad mechanism tags."""
    return any(term in text for term in GENERIC_CANCER_TERMS)


def get_searchable_text(fm: dict, body: str) -> str:
    """Combine title, abstract, MeSH terms, and body into searchable text."""
    parts = [
        fm.get("title", ""),
        " ".join(fm.get("mesh_terms", [])),
        " ".join(fm.get("diseases_annotated", [])),
        " ".join(fm.get("genes", [])),
        " ".join(fm.get("drugs", [])),
    ]

    # Extract abstract from body (between ## Abstract and next ##)
    abstract_match = re.search(r"## Abstract\n\n?(.*?)(?=\n## |\Z)", body, re.DOTALL)
    if abstract_match:
        parts.append(abstract_match.group(1))

    return normalize_text(" ".join(parts))


def match_keywords(text: str, keyword_dict: dict) -> list[str]:
    """Match text against keyword dictionary using word-boundary matching.

    Short keywords (<=4 chars) use strict word-boundary regex to avoid
    false positives like 'all' matching 'overall'.
    """
    matched = []
    for tag, keywords in keyword_dict.items():
        for kw in keywords:
            kw_lower = kw.lower()
            if len(kw_lower) <= 4:
                # Strict word-boundary match for short keywords
                if re.search(r'\b' + re.escape(kw_lower) + r'\b', text):
                    matched.append(tag)
                    break
            else:
                if kw_lower in text:
                    matched.append(tag)
                    break
    return sorted(matched)


def text_matches_keyword(text: str, keyword: str) -> bool:
    """Return True when a keyword is present in text with stricter handling for short terms."""
    kw_lower = keyword.lower()
    if len(kw_lower) <= 4:
        return bool(re.search(r'\b' + re.escape(kw_lower) + r'\b', text))
    return kw_lower in text


def match_resistant_states(text: str) -> list[str]:
    """Match resistant states using composite rules rather than single-keyword OR logic."""
    matched = []
    for state, rule in RESISTANT_STATE_RULES.items():
        all_groups = rule.get("all_of", [])
        if all(any(text_matches_keyword(text, kw) for kw in group) for group in all_groups):
            matched.append(state)
    return sorted(matched)


def match_mechanisms(text: str) -> list[str]:
    """Match mechanisms with a coarse cancer-context gate to reduce off-target tags."""
    if not has_cancer_context(text):
        return []
    return match_keywords(text, MECHANISM_KEYWORDS)


def match_evidence_level(fm: dict, text: str) -> str:
    """Match text against evidence level keywords, return best match.

    Uses word-boundary matching for short keywords to avoid false positives.
    Priority order: phase3 > phase2 > phase1 > invivo > invitro > theoretical.
    """
    if is_review_like(fm) or is_protocol_like(fm):
        return ""

    pub_types = [normalize_text(p) for p in fm.get("pub_types", [])]
    for level in ["phase3-clinical", "phase2-clinical", "phase1-clinical"]:
        if any(marker in pub_types for marker in EVIDENCE_PUBTYPE_MARKERS[level]):
            return level

    for level in ["phase3-clinical", "phase2-clinical", "phase1-clinical",
                   "preclinical-invivo", "preclinical-invitro", "theoretical"]:
        for kw in EVIDENCE_LEVEL_KEYWORDS[level]:
            kw_lower = kw.lower()
            if len(kw_lower) <= 4:
                if re.search(r'\b' + re.escape(kw_lower) + r'\b', text):
                    return level
            else:
                if kw_lower in text:
                    return level
    return ""


def write_tag_files(tag_type: str, tag_pmids: dict[str, list[str]]) -> None:
    """Write tag index files. Each file contains one PMID per line."""
    tag_dir = TAGS_DIR / tag_type
    tag_dir.mkdir(parents=True, exist_ok=True)

    # Clear existing tag files in this directory
    for f in tag_dir.glob("*.txt"):
        f.unlink()

    for tag, pmids in sorted(tag_pmids.items()):
        if pmids:
            filepath = tag_dir / f"{tag}.txt"
            filepath.write_text("\n".join(sorted(pmids)) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Auto-tag articles and build tag indexes")
    parser.add_argument("--dry-run", action="store_true", help="Show tags without writing")
    args = parser.parse_args()

    files = sorted(PMID_DIR.glob("*.md"))
    if not files:
        print("No articles found in corpus.")
        return

    print(f"Tagging {len(files)} articles...")

    # Accumulators for tag index files
    mechanism_pmids: dict[str, list[str]] = {k: [] for k in MECHANISM_KEYWORDS}
    biology_process_pmids: dict[str, list[str]] = {k: [] for k in BIOLOGY_PROCESS_KEYWORDS}
    cancer_pmids: dict[str, list[str]] = {k: [] for k in CANCER_TYPE_KEYWORDS}
    evidence_pmids: dict[str, list[str]] = {k: [] for k in EVIDENCE_LEVEL_KEYWORDS}
    resistant_state_pmids: dict[str, list[str]] = {k: [] for k in RESISTANT_STATE_RULES}
    journal_pmids: dict[str, list[str]] = {}

    stats = {"mechanisms": 0, "biology_processes": 0, "cancer_types": 0, "evidence": 0, "resistant_states": 0}

    for filepath in tqdm(files, desc="  Tagging"):
        fm, body = load_article(filepath)
        if not fm:
            continue

        pmid = fm.get("pmid", filepath.stem)
        text = get_searchable_text(fm, body)

        # Match
        mechanisms = match_mechanisms(text)
        biology_processes = match_keywords(text, BIOLOGY_PROCESS_KEYWORDS)
        cancer_types = match_keywords(text, CANCER_TYPE_KEYWORDS)
        evidence = match_evidence_level(fm, text)
        resistant_states = match_resistant_states(text)

        # Update frontmatter
        fm["mechanisms"] = mechanisms
        fm["biology_processes"] = biology_processes
        fm["cancer_types"] = cancer_types
        fm["evidence_level"] = evidence
        fm["resistant_states"] = resistant_states

        if not args.dry_run:
            save_article(filepath, fm, body)

        # Accumulate for tag files
        for m in mechanisms:
            mechanism_pmids[m].append(pmid)
        for b in biology_processes:
            biology_process_pmids[b].append(pmid)
        for c in cancer_types:
            cancer_pmids[c].append(pmid)
        if evidence:
            evidence_pmids[evidence].append(pmid)
        for r in resistant_states:
            resistant_state_pmids[r].append(pmid)

        # Journal tag
        journal = fm.get("journal", "")
        if journal:
            journal_key = re.sub(r"[^a-z0-9]+", "-", journal.lower()).strip("-")
            if journal_key not in journal_pmids:
                journal_pmids[journal_key] = []
            journal_pmids[journal_key].append(pmid)

        # Stats
        if mechanisms:
            stats["mechanisms"] += 1
        if biology_processes:
            stats["biology_processes"] += 1
        if cancer_types:
            stats["cancer_types"] += 1
        if evidence:
            stats["evidence"] += 1
        if resistant_states:
            stats["resistant_states"] += 1

    # Write tag index files
    if not args.dry_run:
        print("\nWriting tag index files...")
        write_tag_files("by-mechanism", mechanism_pmids)
        write_tag_files("by-biology-process", biology_process_pmids)
        write_tag_files("by-cancer-type", cancer_pmids)
        write_tag_files("by-evidence-level", evidence_pmids)
        write_tag_files("by-resistant-state", resistant_state_pmids)
        write_tag_files("by-journal", journal_pmids)

    # Print summary
    print(f"\nTagging complete:")
    print(f"  Articles with mechanism tags: {stats['mechanisms']}/{len(files)}")
    print(f"  Articles with biology-process tags: {stats['biology_processes']}/{len(files)}")
    print(f"  Articles with cancer type tags: {stats['cancer_types']}/{len(files)}")
    print(f"  Articles with evidence level: {stats['evidence']}/{len(files)}")
    print(f"  Articles with resistant-state tags: {stats['resistant_states']}/{len(files)}")

    print(f"\nMechanism distribution:")
    for tag, pmids in sorted(mechanism_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nBiology-process distribution:")
    for tag, pmids in sorted(biology_process_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nCancer type distribution:")
    for tag, pmids in sorted(cancer_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nEvidence level distribution:")
    for tag, pmids in sorted(evidence_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nResistant-state distribution:")
    for tag, pmids in sorted(resistant_state_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nJournals represented: {len([j for j, p in journal_pmids.items() if p])}")

    if not args.dry_run:
        print(f"\nNext step:")
        print(f"  python build_index.py    # Rebuild INDEX.jsonl")


if __name__ == "__main__":
    main()
