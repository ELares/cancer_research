#!/usr/bin/env python3
"""Re-measure the evidence tagger on the gold set with the MeSH fallback (#346).

Re-tags the 100 manually-labeled gold PMIDs twice — with the MeSH fallback OFF
(the committed baseline) and ON (FERRO_MESH_EXPANSION semantics) — by toggling
``tag_articles.EVIDENCE_USE_MESH_FALLBACK`` in-process and running the SAME
``match_evidence_level`` over each article's corpus record. Emits:

  - analysis/evidence-gold-set-v1-retag.csv : per-PMID predicted_baseline /
    predicted_mesh / gold, with a provenance header (does NOT touch the frozen
    predicted_evidence_level column in evidence-gold-set-v1.csv).
  - analysis/evidence-gold-mesh-eval.md     : BEFORE/AFTER binary
    precision+recall, per-level recall lift, and the empty-MeSH-floor split.

Offline and deterministic (no corpus mutation, no network). The keyword tagger
remains the reproducible baseline; this only reports what the off-by-default
MeSH layer would change.

Usage:
    python scripts/retag_gold_set.py

`compute_metrics()` is imported by tests/test_mesh_gold_regression.py as the CI
precision-floor guard, so the measurement and the guard never diverge.
"""

import csv
import sys
from pathlib import Path

# Allow `python scripts/retag_gold_set.py` and `from retag_gold_set import ...`.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import tag_articles
from article_io import load_article
from config import PROJECT_ROOT

INPUT_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-set-v1.csv"
LABEL_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-labels-v1.csv"
PMID_DIR = PROJECT_ROOT / "corpus" / "by-pmid"
RETAG_CSV = PROJECT_ROOT / "analysis" / "evidence-gold-set-v1-retag.csv"
EVAL_MD = PROJECT_ROOT / "analysis" / "evidence-gold-mesh-eval.md"

# Recall lift is concentrated in these "soft" levels (the pub_type-driven
# clinical phases are already at 100% recall and are untouched by the fallback).
SOFT_LEVELS = ["clinical-other", "preclinical-invivo", "preclinical-invitro", "theoretical"]


def _load_gold():
    """Return [(pmid, sample_mechanism, gold_level)] for the manually-labeled rows."""
    labels = {}
    with open(LABEL_FILE, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            labels[row["pmid"]] = row["gold_evidence_level"]
    rows = []
    with open(INPUT_FILE, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pmid = row["pmid"]
            if pmid in labels:
                rows.append((pmid, row.get("sample_mechanism", ""), labels[pmid]))
    return rows


def _predict(pmid: str, use_mesh: bool) -> str:
    """Run the tagger's evidence matcher on one corpus record with the MeSH
    fallback toggled. Returns "" when the record is absent (counted as a miss)."""
    path = PMID_DIR / f"{pmid}.md"
    if not path.exists():
        return ""
    fm, body = load_article(path)
    text = tag_articles.get_searchable_text(fm, body)
    prev = tag_articles.EVIDENCE_USE_MESH_FALLBACK
    tag_articles.EVIDENCE_USE_MESH_FALLBACK = use_mesh
    try:
        return tag_articles.match_evidence_level(fm, text)
    finally:
        tag_articles.EVIDENCE_USE_MESH_FALLBACK = prev


def _has_mesh(pmid: str) -> bool:
    path = PMID_DIR / f"{pmid}.md"
    if not path.exists():
        return False
    fm, _ = load_article(path)
    return bool(fm.get("mesh_terms"))


def _binary(predictions, gold):
    """Binary evidence-detection TP/FP/FN/recall/precision. A gold level other
    than none-applicable is a positive; any non-empty prediction is a positive."""
    tp = fp = fn = 0
    for pred, g in zip(predictions, gold):
        gold_pos = g not in ("", "none-applicable")
        pred_pos = bool(pred)
        if gold_pos and pred_pos:
            tp += 1
        elif (not gold_pos) and pred_pos:
            fp += 1
        elif gold_pos and (not pred_pos):
            fn += 1
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "recall": recall, "precision": precision}


def compute_metrics():
    """Re-tag the gold set both ways and return a dict with baseline + mesh binary
    metrics, the per-level recall lift, and the empty-MeSH-floor split. Pure
    (no file writes), so the CI regression test can call it directly."""
    rows = _load_gold()
    pmids = [r[0] for r in rows]
    gold = [r[2] for r in rows]
    base = [_predict(p, False) for p in pmids]
    mesh = [_predict(p, True) for p in pmids]

    # Empty-MeSH-floor split among the binary false-negatives of the baseline.
    fn_pmids = [p for p, b, g in zip(pmids, base, gold)
                if g not in ("", "none-applicable") and not b]
    fn_empty_mesh = sum(1 for p in fn_pmids if not _has_mesh(p))

    per_level = {}
    for level in SOFT_LEVELS:
        gold_n = sum(1 for g in gold if g == level)
        base_hit = sum(1 for pr, g in zip(base, gold) if g == level and pr == level)
        mesh_hit = sum(1 for pr, g in zip(mesh, gold) if g == level and pr == level)
        per_level[level] = {"gold": gold_n, "baseline": base_hit, "mesh": mesh_hit}

    # Count gold PMIDs whose corpus record actually resolved, so the CI guard can
    # assert the offline-corpus dependency directly rather than relying on the
    # recall floor to incidentally catch only large corpus loss (a missing record
    # is silently a miss in `_predict`).
    records_found = sum(1 for p in pmids if (PMID_DIR / f"{p}.md").exists())

    return {
        "n": len(rows),
        "records_found": records_found,
        "baseline": _binary(base, gold),
        "mesh": _binary(mesh, gold),
        "per_level": per_level,
        "fn_total": len(fn_pmids),
        "fn_empty_mesh": fn_empty_mesh,
        "rows": list(zip(pmids, [r[1] for r in rows], base, mesh, gold)),
    }


def _pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def write_outputs(m):
    # Per-PMID retag CSV (frozen baseline column in the source CSV untouched).
    with open(RETAG_CSV, "w", newline="", encoding="utf-8") as f:
        f.write(
            "# Regenerate: python scripts/retag_gold_set.py  (#346 MeSH evidence fallback).\n"
            "# predicted_baseline = tagger with FERRO_MESH_EXPANSION off (the committed baseline);\n"
            "# predicted_mesh = with the MeSH fallback on. gold_evidence_level is the manual label.\n"
        )
        w = csv.writer(f)
        w.writerow(["pmid", "sample_mechanism", "predicted_baseline", "predicted_mesh", "gold_evidence_level"])
        for pmid, mech, base, mesh, gold in m["rows"]:
            w.writerow([pmid, mech, base, mesh, gold])

    b, x = m["baseline"], m["mesh"]
    lines = []
    lines.append("# Evidence-tagger gold-set re-measurement: MeSH fallback (#346)\n")
    lines.append(
        f"Generated by `scripts/retag_gold_set.py` over the {m['n']} manually-labeled "
        "gold rows (`analysis/evidence-gold-labels-v1.csv`). BEFORE = the keyword/pub_type "
        "tagger (the byte-identical production default); AFTER = the same tagger with the "
        "off-by-default MeSH-descriptor fallback enabled (`FERRO_MESH_EXPANSION=1`). The "
        "MeSH layer is NOT on in the production corpus run; these AFTER numbers report what "
        "enabling it would change.\n"
    )
    lines.append("## Binary evidence detection (gold-positive vs predicted-positive)\n")
    lines.append("| | TP | FP | FN | recall | precision |")
    lines.append("|---|---|---|---|---|---|")
    lines.append(f"| **baseline (MeSH off)** | {b['tp']} | {b['fp']} | {b['fn']} | {_pct(b['recall'])} | {_pct(b['precision'])} |")
    lines.append(f"| **MeSH fallback (on)** | {x['tp']} | {x['fp']} | {x['fn']} | {_pct(x['recall'])} | {_pct(x['precision'])} |")
    lines.append("")
    lines.append(
        f"Recall lifts from **{_pct(b['recall'])}** to **{_pct(x['recall'])}** "
        f"(+{x['tp'] - b['tp']} true positives) at **{_pct(x['precision'])}** precision "
        f"(baseline {_pct(b['precision'])}; {x['fp'] - b['fp']} net-new false positive(s)).\n"
    )
    lines.append("## Per-level recall lift (the pub_type-driven clinical phases are already 100% and untouched)\n")
    lines.append("| level | gold rows | baseline hits | MeSH-on hits |")
    lines.append("|---|---|---|---|")
    for level in SOFT_LEVELS:
        d = m["per_level"][level]
        lines.append(f"| {level} | {d['gold']} | {d['baseline']} | {d['mesh']} |")
    lines.append("")
    lines.append("## The unrecoverable MeSH floor\n")
    lines.append(
        f"Of the {m['fn_total']} baseline binary false-negatives, **{m['fn_empty_mesh']} carry empty "
        "`mesh_terms`** and so cannot be recovered by any MeSH layer (recent / not-yet-indexed "
        "records). That residue, plus the gold-positives whose MeSH is present but non-discriminative "
        "(e.g. in-vitro studies PubMed rarely tags with a Cell Line descriptor), is the empirical case "
        "for the deferred embedding/semantic-retrieval leg of #346. The keyword tagger remains the "
        "reproducible baseline; corpus-level non-detection stays provisional.\n"
    )
    EVAL_MD.write_text("\n".join(lines), encoding="utf-8")


def main():
    m = compute_metrics()
    write_outputs(m)
    b, x = m["baseline"], m["mesh"]
    print(
        f"gold rows={m['n']}  baseline recall={_pct(b['recall'])}/prec={_pct(b['precision'])}  "
        f"MeSH-on recall={_pct(x['recall'])}/prec={_pct(x['precision'])}  "
        f"(+{x['tp'] - b['tp']} TP, {x['fp'] - b['fp']} net-new FP)"
    )
    print(f"wrote {RETAG_CSV.relative_to(PROJECT_ROOT)} and {EVAL_MD.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
