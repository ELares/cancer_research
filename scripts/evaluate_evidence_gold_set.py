#!/usr/bin/env python3
"""
Evaluate heuristic evidence tags against a manually labeled gold set.

Usage:
    python scripts/evaluate_evidence_gold_set.py
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

from config import PROJECT_ROOT

INPUT_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-set-v1.csv"
LABEL_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-labels-v1.csv"
OUTPUT_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-eval.md"

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


def f1(tp: int, fp: int, fn: int) -> float:
    if tp == 0:
        return 0.0
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    return 2 * precision * recall / (precision + recall)


def safe_ratio(num: int, den: int) -> float:
    return num / den if den else 0.0


def main() -> None:
    with open(INPUT_FILE, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    with open(LABEL_FILE, newline="", encoding="utf-8") as handle:
        labels = {row["pmid"]: row for row in csv.DictReader(handle)}

    labeled = []
    unlabeled_pmids = []
    for row in rows:
        label = labels.get(row["pmid"])
        if not label:
            unlabeled_pmids.append(row["pmid"])
            continue
        merged = dict(row)
        merged["gold_evidence_level"] = label["gold_evidence_level"]
        merged["gold_label_status"] = label.get("gold_label_status", "")
        merged["gold_notes"] = label.get("gold_notes", "")
        labeled.append(merged)

    if not labeled:
        raise SystemExit(f"No labeled rows found in {LABEL_FILE}")

    exact = sum(1 for row in labeled if row["predicted_evidence_level"] == row["gold_evidence_level"])
    gold_positive = sum(1 for row in labeled if row["gold_evidence_level"] != "none-applicable")
    pred_positive = sum(1 for row in labeled if row["predicted_evidence_level"])
    tp_binary = sum(
        1
        for row in labeled
        if row["gold_evidence_level"] != "none-applicable" and row["predicted_evidence_level"]
    )
    fp_binary = sum(
        1
        for row in labeled
        if row["gold_evidence_level"] == "none-applicable" and row["predicted_evidence_level"]
    )
    fn_binary = sum(
        1
        for row in labeled
        if row["gold_evidence_level"] != "none-applicable" and not row["predicted_evidence_level"]
    )

    per_label = {}
    for label in LABELS:
        tp = sum(
            1
            for row in labeled
            if row["gold_evidence_level"] == label and row["predicted_evidence_level"] == label
        )
        fp = sum(
            1
            for row in labeled
            if row["gold_evidence_level"] != label and row["predicted_evidence_level"] == label
        )
        fn = sum(
            1
            for row in labeled
            if row["gold_evidence_level"] == label and row["predicted_evidence_level"] != label
        )
        per_label[label] = {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": safe_ratio(tp, tp + fp),
            "recall": safe_ratio(tp, tp + fn),
            "f1": f1(tp, fp, fn),
        }

    per_mechanism = defaultdict(list)
    for row in labeled:
        per_mechanism[row["sample_mechanism"]].append(row)

    confusion = Counter()
    confusion_examples = defaultdict(list)
    for row in labeled:
        predicted = row["predicted_evidence_level"] or "unclassified"
        gold = row["gold_evidence_level"]
        if gold == predicted:
            continue
        key = (gold, predicted)
        confusion[key] += 1
        if len(confusion_examples[key]) < 3:
            confusion_examples[key].append(f"{row['pmid']} ({row['sample_mechanism']})")

    lines = ["# Evidence Gold-Set Evaluation\n"]
    lines.append(
        f"Manual labels are present for **{len(labeled)}** sampled records from "
        f"`analysis/evidence-gold-set-v1.csv`, with label assignments stored in "
        f"`analysis/evidence-gold-labels-v1.csv`.\n"
    )
    lines.append(
        "- Sampling design: 10 `predicted-tagged` and 10 `predicted-untagged` rows for each of "
        "`immunotherapy`, `mRNA-vaccine`, `electrochemical-therapy`, `ttfields`, and `synthetic-lethality`."
    )
    if unlabeled_pmids:
        lines.append(f"- Unlabeled sampled rows still pending: **{len(unlabeled_pmids)}**")
    else:
        lines.append("- All sampled rows currently have manual labels.")
    lines.append("")

    lines.append("## Overall Metrics\n")
    lines.append(f"- Exact-label accuracy: **{exact}/{len(labeled)} ({exact/len(labeled):.1%})**")
    lines.append(
        f"- Binary evidence-detection precision: **{safe_ratio(tp_binary, tp_binary + fp_binary):.1%}** "
        f"({tp_binary}/{tp_binary + fp_binary or 1})"
    )
    lines.append(
        f"- Binary evidence-detection recall: **{safe_ratio(tp_binary, tp_binary + fn_binary):.1%}** "
        f"({tp_binary}/{tp_binary + fn_binary or 1})"
    )
    lines.append(
        f"- Binary evidence-detection F1: **{f1(tp_binary, fp_binary, fn_binary):.3f}**"
    )
    lines.append(
        f"- Gold positive rows: **{gold_positive}**; predicted positive rows: **{pred_positive}**\n"
    )

    lines.append("## Per-Label Metrics\n")
    lines.append("| Label | TP | FP | FN | Precision | Recall | F1 |")
    lines.append("|---|---|---|---|---|---|---|")
    for label in LABELS:
        metrics = per_label[label]
        lines.append(
            f"| **{label}** | {metrics['tp']} | {metrics['fp']} | {metrics['fn']} | "
            f"{metrics['precision']:.1%} | {metrics['recall']:.1%} | {metrics['f1']:.3f} |"
        )

    lines.append("\n## Per-Mechanism Exact Accuracy\n")
    lines.append("| Mechanism | Labeled rows | Exact accuracy | Predicted positive | Gold positive |")
    lines.append("|---|---|---|---|---|")
    for mechanism in sorted(per_mechanism):
        subset = per_mechanism[mechanism]
        subset_exact = sum(1 for row in subset if row["predicted_evidence_level"] == row["gold_evidence_level"])
        subset_pred_positive = sum(1 for row in subset if row["predicted_evidence_level"])
        subset_gold_positive = sum(1 for row in subset if row["gold_evidence_level"] != "none-applicable")
        lines.append(
            f"| **{mechanism}** | {len(subset)} | {subset_exact/len(subset):.1%} | "
            f"{subset_pred_positive} | {subset_gold_positive} |"
        )

    lines.append("\n## Most Common Confusions\n")
    for (gold_label, predicted_label), count in confusion.most_common(12):
        examples = ", ".join(confusion_examples[(gold_label, predicted_label)])
        lines.append(f"- **{gold_label} -> {predicted_label}**: {count}  ")
        lines.append(f"  Example PMIDs: {examples}")

    lines.append("\n## Interpretation\n")
    lines.append(
        "- The current evidence tagger behaves like a conservative detector: it rarely assigns evidence to rows manually labeled `none-applicable`, but it misses a large share of valid evidence-bearing records."
    )
    lines.append(
        "- The largest blind spots in this sample are unlabeled `theoretical`, `clinical-other`, and `preclinical-invitro` studies. That lines up with earlier qualitative audit notes."
    )
    lines.append(
        "- The gold set supports using coverage-aware manuscript language. The current heuristic is much more reliable for `if tagged, usually real` than for `if untagged, probably absent`."
    )

    OUTPUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote evaluation report to {OUTPUT_FILE}")
    print(f"  Labeled rows: {len(labeled)}")
    print(f"  Exact accuracy: {exact/len(labeled):.1%}")


if __name__ == "__main__":
    main()
