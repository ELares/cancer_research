#!/usr/bin/env python3
"""
Build a deterministic gold-label sampling sheet for evidence-tier evaluation.

The sample is stratified by:
- target mechanism
- whether the current heuristic already assigned an evidence level

Usage:
    python scripts/sample_evidence_gold_set.py               # v1 (100 articles, 5 mechanisms)
    python scripts/sample_evidence_gold_set.py --version v2   # v2 (250 articles, 10 mechanisms)
"""

import argparse
import csv
import random
import re

from article_io import load_article
from config import PMID_DIR, PROJECT_ROOT
from evidence_utils import is_protocol_like, is_review_like, normalize_text

TARGET_MECHANISMS_V1 = [
    "immunotherapy",
    "mRNA-vaccine",
    "electrochemical-therapy",
    "ttfields",
    "synthetic-lethality",
]

TARGET_MECHANISMS_V2 = [
    "immunotherapy",
    "nanoparticle",
    "car-t",
    "oncolytic-virus",
    "synthetic-lethality",
    "crispr",
    "antibody-drug-conjugate",
    "ttfields",
    "sonodynamic",
    "electrochemical-therapy",
]

SAMPLES_PER_BUCKET_V1 = 10
SAMPLES_PER_BUCKET_V2 = 12  # 12 tagged + 13 untagged = 25 per mechanism
SAMPLES_UNTAGGED_V2 = 13

RANDOM_SEED = 33
RANDOM_SEED_V2 = 44
OUTPUT_FILE_V1 = PROJECT_ROOT / "analysis" / "evidence-gold-set-v1.csv"
OUTPUT_FILE_V2 = PROJECT_ROOT / "analysis" / "evidence-gold-set-v2.csv"


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


FIELDNAMES = [
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
]


def load_all_articles():
    """Load all primary-study-like articles from the corpus."""
    articles = []
    for filepath in sorted(PMID_DIR.glob("*.md")):
        frontmatter, body = load_article(filepath)
        if not frontmatter:
            continue
        if not is_primary_study_like(frontmatter):
            continue
        articles.append((frontmatter, body))
    return articles


def sample_v1(articles):
    """Original v1 sampling: 5 mechanisms × 10 tagged + 10 untagged = 100."""
    rng = random.Random(RANDOM_SEED)
    rows = []
    sampled_pmids = set()

    for mechanism in TARGET_MECHANISMS_V1:
        tagged_pool = [a for a in articles if mechanism in a[0].get("mechanisms", []) and a[0].get("evidence_level")]
        untagged_pool = [a for a in articles if mechanism in a[0].get("mechanisms", []) and not a[0].get("evidence_level")]

        for bucket_name, pool, n in [("predicted-tagged", tagged_pool, SAMPLES_PER_BUCKET_V1),
                                     ("predicted-untagged", untagged_pool, SAMPLES_PER_BUCKET_V1)]:
            if len(pool) < n:
                raise SystemExit(f"Not enough articles for {mechanism} / {bucket_name}: need {n}, found {len(pool)}")
            sampled = rng.sample(pool, n)
            for fm, body in sorted(sampled, key=lambda item: int(item[0].get("pmid", 0))):
                pmid = fm.get("pmid", "")
                if pmid in sampled_pmids:
                    raise SystemExit(f"Duplicate PMID: {pmid}")
                sampled_pmids.add(pmid)
                rows.append(build_row(fm, body, mechanism, bucket_name))

    return rows


def sample_v2(articles):
    """Expanded v2 sampling: 10 mechanisms × (12 tagged + 13 untagged) = 250.

    Includes all v1 articles (with their existing gold labels) plus 150 new articles.
    """
    rng = random.Random(RANDOM_SEED_V2)
    rows = []
    sampled_pmids = set()

    # Load v1 gold labels to carry forward
    v1_labels = {}
    label_file = PROJECT_ROOT / "analysis" / "evidence-gold-labels-v1.csv"
    if label_file.exists():
        with open(label_file, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                v1_labels[row["pmid"]] = row

    # Include all v1 articles first (preserve their labels)
    v1_set_file = PROJECT_ROOT / "analysis" / "evidence-gold-set-v1.csv"
    if v1_set_file.exists():
        with open(v1_set_file, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                pmid = row["pmid"]
                # Carry forward v1 gold labels
                if pmid in v1_labels:
                    row["gold_evidence_level"] = v1_labels[pmid].get("gold_evidence_level", "")
                    row["gold_label_status"] = v1_labels[pmid].get("gold_label_status", "manual-v1")
                    row["gold_notes"] = v1_labels[pmid].get("gold_notes", "")
                rows.append(row)
                sampled_pmids.add(pmid)

    # Sample new articles from the expanded mechanism pool
    for mechanism in TARGET_MECHANISMS_V2:
        tagged_pool = [a for a in articles
                       if mechanism in a[0].get("mechanisms", [])
                       and a[0].get("evidence_level")
                       and a[0].get("pmid", "") not in sampled_pmids]
        untagged_pool = [a for a in articles
                         if mechanism in a[0].get("mechanisms", [])
                         and not a[0].get("evidence_level")
                         and a[0].get("pmid", "") not in sampled_pmids]

        # How many more do we need for this mechanism?
        existing_for_mech = sum(1 for r in rows if r.get("sample_mechanism") == mechanism)
        target_tagged = SAMPLES_PER_BUCKET_V2
        target_untagged = SAMPLES_UNTAGGED_V2
        need_tagged = max(0, target_tagged - sum(1 for r in rows
                                                  if r.get("sample_mechanism") == mechanism
                                                  and r.get("sample_bucket") == "predicted-tagged"))
        need_untagged = max(0, target_untagged - sum(1 for r in rows
                                                      if r.get("sample_mechanism") == mechanism
                                                      and r.get("sample_bucket") == "predicted-untagged"))

        for bucket_name, pool, n in [("predicted-tagged", tagged_pool, need_tagged),
                                     ("predicted-untagged", untagged_pool, need_untagged)]:
            if n <= 0:
                continue
            available = min(len(pool), n)
            if available < n:
                print(f"  Warning: {mechanism}/{bucket_name}: need {n} more, only {available} available")
            sampled = rng.sample(pool, available)
            for fm, body in sorted(sampled, key=lambda item: int(item[0].get("pmid", 0))):
                pmid = fm.get("pmid", "")
                if pmid in sampled_pmids:
                    continue
                sampled_pmids.add(pmid)
                rows.append(build_row(fm, body, mechanism, bucket_name))

    return rows


def write_csv(rows, output_file):
    """Write gold-set rows to CSV."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Sample evidence gold set.")
    parser.add_argument("--version", choices=["v1", "v2"], default="v1",
                        help="v1: 100 articles (5 mechanisms), v2: 250 articles (10 mechanisms)")
    args = parser.parse_args()

    print("Loading articles...")
    articles = load_all_articles()
    print(f"  Loaded {len(articles)} primary-study-like articles")

    if args.version == "v1":
        rows = sample_v1(articles)
        output = OUTPUT_FILE_V1
        mechs = TARGET_MECHANISMS_V1
    else:
        rows = sample_v2(articles)
        output = OUTPUT_FILE_V2
        mechs = TARGET_MECHANISMS_V2

    write_csv(rows, output)
    print(f"Wrote {len(rows)} sampled rows to {output}")
    print(f"  Target mechanisms: {', '.join(mechs)}")
    v1_count = sum(1 for r in rows if r.get("gold_label_status", "").startswith("manual"))
    pending = sum(1 for r in rows if r.get("gold_label_status", "") == "pending")
    print(f"  Already labeled (from v1): {v1_count}")
    print(f"  Pending labeling: {pending}")


if __name__ == "__main__":
    main()
