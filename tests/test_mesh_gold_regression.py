#!/usr/bin/env python3
"""CI precision-floor regression guard for the MeSH evidence fallback (#346).

The biggest risk of the MeSH layer is a precision regression: a descriptor that
encodes clinical CONTEXT rather than evidence DESIGN (or a future broad-descriptor
edit) would erode the keyword tagger's 96% precision. This test re-measures the
fallback (flag ON) on the gold set via the SAME `compute_metrics` the report uses
and fails loudly if precision falls into the danger zone or the recall gain
evaporates — so a leaky descriptor change cannot land silently.

Thresholds are deliberately robust, not point-pinned: the gold set is small (87
positives / 13 negatives), so one FP is a ~1.5-2 point precision swing and
pinning to the measured 95.2% would false-fail on a benign 1-FP change.
precision >= 0.93 resolves a >=2-net-new-FP regression on this 13-negative set
(at the measured fp=3, +2 FP -> 0.922 fails); a single-FP leak is below the set's
resolution, so this guard catches a broad/leaky descriptor change, not a subtle
one-article leak. recall >= 0.62 ensures the layer still materially beats the
~55% keyword baseline. (A larger, independent gold set would sharpen both floors.)

Run: pytest tests/test_mesh_gold_regression.py -v
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

PRECISION_FLOOR = 0.93
RECALL_FLOOR = 0.62


def test_mesh_fallback_holds_precision_and_lifts_recall():
    from retag_gold_set import compute_metrics

    m = compute_metrics()
    base, mesh = m["baseline"], m["mesh"]

    # The harness must actually find the gold rows in the CSV...
    assert m["n"] >= 90, f"expected ~100 gold rows, got {m['n']}"
    # ...AND resolve their corpus records: a missing record is silently a miss in
    # _predict, so assert corpus presence directly rather than relying on the
    # recall floor to incidentally catch only large corpus loss.
    assert m["records_found"] >= 90, (
        f"only {m['records_found']}/{m['n']} gold PMIDs resolved to a corpus record — "
        f"the offline corpus is incomplete; the measurement would be silently understated."
    )
    assert base["tp"] > 0 and mesh["tp"] > 0, "no positives detected — measurement is broken"

    assert mesh["precision"] >= PRECISION_FLOOR, (
        f"MeSH-fallback precision {mesh['precision']:.3f} < floor {PRECISION_FLOOR}: a descriptor "
        f"edit likely promoted clinical-context / opinion articles. Baseline precision "
        f"{base['precision']:.3f}. Tighten the descriptor sets in tag_articles.EVIDENCE_MESH_MARKERS."
    )
    assert mesh["recall"] >= RECALL_FLOOR, (
        f"MeSH-fallback recall {mesh['recall']:.3f} < floor {RECALL_FLOOR}: the layer no longer "
        f"materially beats the keyword baseline ({base['recall']:.3f}). Did a descriptor set lose "
        f"its corpus coverage?"
    )
    # The layer must add recall over the baseline (its entire purpose).
    assert mesh["recall"] > base["recall"], (
        f"MeSH fallback did not improve recall: {mesh['recall']:.3f} <= baseline {base['recall']:.3f}"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
