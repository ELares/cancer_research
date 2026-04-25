#!/usr/bin/env python3
"""Compute inter-rater agreement (Cohen's kappa) for evidence gold-set labeling.

Reads primary and second-rater label CSVs, computes agreement metrics,
and outputs a markdown report. See analysis/evidence-adjudication-protocol.md
for the adjudication process.

Usage:
    python scripts/compute_inter_rater_agreement.py \
        --primary analysis/evidence-gold-labels-v2.csv \
        --second  analysis/evidence-gold-labels-v2-second-rater.csv
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

LABELS = [
    "phase3-clinical",
    "phase2-clinical",
    "phase1-clinical",
    "clinical-other",
    "preclinical-invivo",
    "preclinical-invitro",
    "theoretical",
    "none-applicable",
]

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_FILE = PROJECT_ROOT / "analysis" / "inter-rater-agreement.md"


def cohens_kappa(primary, second, labels):
    """Compute Cohen's kappa (unweighted) for two rater label lists."""
    n = len(primary)
    if n == 0:
        return 0.0

    # Observed agreement
    agree = sum(1 for a, b in zip(primary, second) if a == b)
    p_o = agree / n

    # Expected agreement by chance
    primary_counts = Counter(primary)
    second_counts = Counter(second)
    p_e = sum(primary_counts.get(l, 0) * second_counts.get(l, 0) for l in labels) / (n * n)

    if p_e >= 1.0:
        return 1.0
    return (p_o - p_e) / (1.0 - p_e)


def cohens_kappa_weighted(primary, second, labels):
    """Compute linearly weighted Cohen's kappa (ordinal-aware)."""
    n = len(primary)
    if n == 0:
        return 0.0

    n_labels = len(labels)
    label_idx = {l: i for i, l in enumerate(labels)}

    # Build confusion matrix
    matrix = [[0] * n_labels for _ in range(n_labels)]
    for a, b in zip(primary, second):
        i = label_idx.get(a, n_labels - 1)
        j = label_idx.get(b, n_labels - 1)
        matrix[i][j] += 1

    # Weight matrix (linear): w_ij = |i - j| / (n_labels - 1)
    w_o = 0.0  # observed disagreement
    w_e = 0.0  # expected disagreement
    row_sums = [sum(row) for row in matrix]
    col_sums = [sum(matrix[i][j] for i in range(n_labels)) for j in range(n_labels)]

    for i in range(n_labels):
        for j in range(n_labels):
            weight = abs(i - j) / max(n_labels - 1, 1)
            w_o += weight * matrix[i][j] / n
            w_e += weight * row_sums[i] * col_sums[j] / (n * n)

    if w_e == 0:
        return 1.0
    return 1.0 - w_o / w_e


def main():
    parser = argparse.ArgumentParser(description="Compute inter-rater agreement.")
    parser.add_argument("--primary", required=True, help="Primary rater labels CSV")
    parser.add_argument("--second", required=True, help="Second rater labels CSV")
    parser.add_argument("--output", default=str(OUTPUT_FILE), help="Output markdown file")
    args = parser.parse_args()

    # Load labels
    primary_labels = {}
    with open(args.primary, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            primary_labels[row["pmid"]] = row.get("gold_evidence_level", "")

    second_labels = {}
    with open(args.second, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            second_labels[row["pmid"]] = row.get("gold_evidence_level", "")

    # Find overlap
    overlap_pmids = sorted(set(primary_labels) & set(second_labels))
    if not overlap_pmids:
        print("ERROR: No overlapping PMIDs between primary and second rater.")
        raise SystemExit(1)

    primary = [primary_labels[p] for p in overlap_pmids]
    second = [second_labels[p] for p in overlap_pmids]

    # Filter to labeled (non-empty) pairs only
    pairs = [(a, b) for a, b in zip(primary, second) if a and b]
    if not pairs:
        print("ERROR: No labeled overlapping articles.")
        raise SystemExit(1)

    p_labels, s_labels = zip(*pairs)

    # Compute metrics
    kappa = cohens_kappa(p_labels, s_labels, LABELS)
    kappa_w = cohens_kappa_weighted(p_labels, s_labels, LABELS)
    agree_count = sum(1 for a, b in pairs if a == b)
    n = len(pairs)

    # Per-tier agreement
    tier_counts = Counter()
    tier_agree = Counter()
    for a, b in pairs:
        tier_counts[a] += 1
        if a == b:
            tier_agree[a] += 1

    # Build report
    lines = ["# Inter-Rater Agreement Report\n"]
    lines.append(f"Overlap articles: {n}\n")
    lines.append(f"## Overall Metrics\n")
    lines.append(f"- Raw agreement: {agree_count}/{n} ({100 * agree_count / n:.1f}%)")
    lines.append(f"- Cohen's kappa (unweighted): {kappa:.3f}")
    lines.append(f"- Cohen's kappa (linear weighted): {kappa_w:.3f}\n")
    lines.append(f"## Per-Tier Agreement\n")
    lines.append("| Tier | N (primary) | Agreed | Agreement % |")
    lines.append("|------|-------------|--------|-------------|")
    for tier in LABELS:
        tc = tier_counts.get(tier, 0)
        ta = tier_agree.get(tier, 0)
        pct = f"{100 * ta / tc:.0f}%" if tc >= 5 else "N<5"
        lines.append(f"| {tier} | {tc} | {ta} | {pct} |")

    lines.append(f"\n## Confusion Matrix (Primary × Second Rater)\n")
    label_idx = {l: i for i, l in enumerate(LABELS)}
    matrix = [[0] * len(LABELS) for _ in range(len(LABELS))]
    for a, b in pairs:
        i = label_idx.get(a, len(LABELS) - 1)
        j = label_idx.get(b, len(LABELS) - 1)
        matrix[i][j] += 1

    header = "| Primary \\ Second | " + " | ".join(l[:10] for l in LABELS) + " |"
    sep = "|---|" + "|".join(["---"] * len(LABELS)) + "|"
    lines.append(header)
    lines.append(sep)
    for i, label in enumerate(LABELS):
        row = f"| {label[:10]} | " + " | ".join(str(matrix[i][j]) for j in range(len(LABELS))) + " |"
        lines.append(row)

    output_path = Path(args.output)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Report written to {output_path}")
    print(f"  Cohen's kappa: {kappa:.3f} (weighted: {kappa_w:.3f})")
    print(f"  Agreement: {agree_count}/{n} ({100 * agree_count / n:.1f}%)")


if __name__ == "__main__":
    main()
