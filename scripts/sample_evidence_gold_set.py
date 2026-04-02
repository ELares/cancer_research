#!/usr/bin/env python3
"""
Build a deterministic gold-label sampling sheet for evidence-tier evaluation.

The sample is stratified by:
- target mechanism
- whether the current heuristic already assigned an evidence level

Usage:
    python scripts/sample_evidence_gold_set.py
"""

import csv
import random
import re
from pathlib import Path

from article_io import load_article
from config import PMID_DIR, PROJECT_ROOT
from evidence_utils import is_protocol_like, is_review_like, normalize_text

TARGET_MECHANISMS = [
    "immunotherapy",
    "mRNA-vaccine",
    "electrochemical-therapy",
    "ttfields",
    "synthetic-lethality",
]

EVIDENCE_LEVELS = [
    "phase3-clinical",
    "phase2-clinical",
    "phase1-clinical",
    "clinical-other",
    "preclinical-invivo",
    "preclinical-invitro",
    "theoretical",
    "none-applicable",
]

SAMPLES_PER_BUCKET = 10
RANDOM_SEED = 33
OUTPUT_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-set-v1.csv"


def extract_abstract(body: str) -> str:
    match = re.search(r"## Abstract\n\n?(.*?)(?=\n## |\Z)", body, re.DOTALL)
    if not match:
        return ""
    return normalize_text(match.group(1))


def is_primary_study_like(frontmatter: dict) -> bool:
    return not is_review_like(frontmatter) and not is_protocol_like(frontmatter)


def build_row(frontmatter: dict, body: str, sample_mechanism: str, sample_bucket: str) -> dict:
    pub_types = frontmatter.get("pub_types", [])
    abstract = extract_abstract(body)
    abstract_excerpt = abstract[:500]
    if len(abstract) > 500:
        abstract_excerpt += "..."

    return {
        "pmid": frontmatter.get("pmid", ""),
        "sample_mechanism": sample_mechanism,
        "sample_bucket": sample_bucket,
        "predicted_evidence_level": frontmatter.get("evidence_level", ""),
        "title": frontmatter.get("title", ""),
        "year": frontmatter.get("year", ""),
        "journal": frontmatter.get("journal", ""),
        "pub_types": " | ".join(pub_types),
        "mechanisms": " | ".join(frontmatter.get("mechanisms", [])),
        "abstract_excerpt": abstract_excerpt,
        "gold_evidence_level": "",
        "gold_label_status": "pending",
        "gold_notes": "",
    }


def main() -> None:
    rng = random.Random(RANDOM_SEED)
    rows = []

    articles = []
    for filepath in sorted(PMID_DIR.glob("*.md")):
        frontmatter, body = load_article(filepath)
        if not frontmatter:
            continue
        if not is_primary_study_like(frontmatter):
            continue
        articles.append((frontmatter, body))

    for mechanism in TARGET_MECHANISMS:
        tagged_pool = []
        untagged_pool = []

        for frontmatter, body in articles:
            if mechanism not in frontmatter.get("mechanisms", []):
                continue
            if frontmatter.get("evidence_level"):
                tagged_pool.append((frontmatter, body))
            else:
                untagged_pool.append((frontmatter, body))

        for bucket_name, pool in (("predicted-tagged", tagged_pool), ("predicted-untagged", untagged_pool)):
            if len(pool) < SAMPLES_PER_BUCKET:
                raise SystemExit(
                    f"Not enough articles for {mechanism} / {bucket_name}: "
                    f"need {SAMPLES_PER_BUCKET}, found {len(pool)}"
                )
            sampled = rng.sample(pool, SAMPLES_PER_BUCKET)
            for frontmatter, body in sorted(sampled, key=lambda item: int(item[0].get("pmid", 0))):
                rows.append(build_row(frontmatter, body, mechanism, bucket_name))

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pmid",
                "sample_mechanism",
                "sample_bucket",
                "predicted_evidence_level",
                "title",
                "year",
                "journal",
                "pub_types",
                "mechanisms",
                "abstract_excerpt",
                "gold_evidence_level",
                "gold_label_status",
                "gold_notes",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} sampled rows to {OUTPUT_FILE}")
    print(f"  Target mechanisms: {', '.join(TARGET_MECHANISMS)}")
    print(f"  Rows per mechanism: {SAMPLES_PER_BUCKET * 2}")


if __name__ == "__main__":
    main()
