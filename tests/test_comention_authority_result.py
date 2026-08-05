"""Guards for the filtered layer's measured precision (#628).

This is the one result in the thread that runs IN FAVOUR of a change I made,
which makes it the one most in need of a guard that fails when it is overstated.
"""

import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

RAW = REPO_ROOT / "analysis" / "comention-authority-result.json"
DOC = REPO_ROOT / "analysis" / "comention-authority-result.md"
JUDGED = REPO_ROOT / "analysis" / "comention"
STRATA = ["corroborated", "abstract-visible", "body-only"]


def _raw():
    return json.loads(RAW.read_text())


def test_every_stratum_judgement_is_committed_and_recomputable():
    d = _raw()
    for s in STRATA:
        f = JUDGED / f"{s}-authority-judgements.csv"
        rows = [r for r in csv.DictReader(f.open()) if r["verdict"] in ("TP", "FP")]
        tp = sum(1 for r in rows if r["verdict"] == "TP")
        assert d["strata"][s]["n"] == len(rows) and d["strata"][s]["tp"] == tp
        assert abs(d["strata"][s]["precision"] - tp / len(rows)) < 1e-12
        # Each verdict must be checkable by a reader.
        assert all("matched_span" in r and "sentence" in r for r in rows)


def test_the_weighted_total_is_the_weighted_total():
    """It must be the strata combined, not a separately asserted number."""
    d = _raw()
    got = sum(v["weight"] * v["precision"] for v in d["strata"].values())
    assert abs(got - d["weighted"]) < 1e-9
    # And the weights must be a partition of the sampled mentions.
    assert abs(sum(v["weight"] for v in d["strata"].values()) - 1.0) < 1e-6


def test_the_interval_is_reported_beside_the_point_estimate():
    """88% on 90 judged mentions is not 88% known to a point."""
    d = _raw()
    lo, hi = d["weighted_ci"]
    assert lo < d["weighted"] < hi
    assert hi - lo > 0.05, "the interval is implausibly tight for n=90"
    flat = " ".join(DOC.read_text().split())
    assert f"[{100*lo:.1f}, {100*hi:.1f}]" in flat


def test_the_confound_and_the_unblinding_are_both_declared():
    """The comparison carries a confound and the judging carries a bias, and
    this is the one result in the thread that favours a change I made."""
    flat = " ".join(DOC.read_text().split())
    assert "The comparison is not." in flat
    assert "unblinded and I knew which layer this was" in flat
    assert "runs the opposite way to the earlier measurements" in flat
    # Recall cost must be stated, not just precision.
    assert "What it cost" in flat and "half the layer's output is gone" in flat.lower()


def test_the_improvement_is_real_across_every_stratum():
    """If any stratum stops improving, the headline is no longer a summary."""
    d = _raw()
    for s in STRATA:
        assert d["strata"][s]["precision"] > d["unfiltered"][s], (
            f"{s} no longer improves under the filter")
    assert d["weighted"] > d["unfiltered"]["weighted"]


def test_the_unblinding_is_bounded_not_only_declared():
    """Declaring a bias is weaker than bounding it.

    Every judgement that could not be settled mechanically is enumerated and
    resolved both ways. If the adverse bound ever falls to the unfiltered
    layer's precision, the improvement is no longer robust to my own judgement
    and the headline must be re-read.
    """
    d = _raw()
    assert "weighted_adverse" in d and "weighted_favourable" in d
    assert d["weighted_adverse"] <= d["weighted"] <= d["weighted_favourable"]
    assert d["weighted_adverse"] > d["unfiltered"]["weighted"], (
        f"resolving every borderline call against the filter gives "
        f"{d['weighted_adverse']:.3f}, at or below the unfiltered "
        f"{d['unfiltered']['weighted']:.3f}; the improvement no longer survives "
        "a hostile reading of the judging")
    # The borderline set must be enumerated, not summarised.
    import comention_authority_result as c
    assert sum(len(v) for v in c.BORDERLINE.values()) >= 10, (
        "too few borderline calls declared for the bound to mean anything")
    flat = " ".join(DOC.read_text().split())
    assert "Resolving every borderline call against the filter" in flat
