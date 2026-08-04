"""Guards for the thesis-leg trajectory analysis (#ATLAS-TRAJ).

The counts in `atlas-thesis-position.md` are a snapshot, and a snapshot cannot
separate "small because it is new and rising" from "small because nobody went
there". The trajectory answers that by comparing each leg's SHARE of the
ferroptosis field between two pooled windows.

The load-bearing result is a NEGATIVE -- every leg is flat -- so these guard
against it silently becoming a positive.
"""

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from atlas_thesis_position import trajectories, wilson  # noqa: E402

RAW = REPO_ROOT / "analysis" / "atlas-thesis-position.json"
DOC = REPO_ROOT / "analysis" / "atlas-thesis-position.md"


def _raw():
    return json.loads(RAW.read_text())


def test_wilson_is_a_real_interval():
    """Bounds stay in [0,1], bracket the estimate, and tighten with n."""
    assert wilson(0, 0) == (0.0, 0.0)
    lo, hi = wilson(0, 176)
    # Zero observations: the lower bound is zero up to floating point, and the
    # upper bound is the useful part (what the share could be and still show 0).
    assert lo < 1e-12 and 0 < hi < 0.05, (lo, hi)
    for k, n in [(4, 176), (15, 4496), (219, 4496), (575, 4496)]:
        lo, hi = wilson(k, n)
        assert 0.0 <= lo <= k / n <= hi <= 1.0, (k, n, lo, hi)
    assert (wilson(50, 5000)[1] - wilson(50, 5000)[0]) < (
        wilson(5, 500)[1] - wilson(5, 500)[0])


def test_trajectory_is_recomputable_from_the_committed_counts():
    """The report's table must be derivable from by_year, not stored prose."""
    d = _raw()
    assert d["trajectory"], "trajectory missing from the committed JSON"
    recomputed = trajectories(d["by_year"])
    for leg, tr in d["trajectory"].items():
        assert recomputed[leg]["moved"] == tr["moved"], leg
        assert abs(recomputed[leg]["late"]["share"] - tr["late"]["share"]) < 1e-12


def test_every_thesis_leg_is_flat_except_the_vocabulary_artifact():
    """The finding. If a leg ever starts moving, this must fail and be re-read."""
    tr = _raw()["trajectory"]
    movers = {k for k, v in tr.items() if v["moved"]}
    assert movers == {"lipid peroxidation"}, (
        f"expected only the MeSH-vocabulary artifact to move, got {movers}")
    assert tr["lipid peroxidation"]["direction"] == "falls"


def test_the_sonodynamic_leg_stayed_small_rather_than_accelerating():
    """The correction this analysis makes to the earlier snapshot.

    A 32-article leg could have been small-and-accelerating. It was not: the
    share is statistically indistinguishable between the windows, and its late
    upper confidence bound stays under 1% of the field.
    """
    sdt = _raw()["trajectory"]["sonodynamic therapy"]
    assert not sdt["moved"], "the SDT leg is no longer flat; re-read the finding"
    assert sdt["late"]["ci"][1] < 0.01, sdt["late"]["ci"]


def test_the_doc_does_not_claim_a_growing_thesis():
    """Prose guard: the field grew 25x, which invites a claim the data refuses."""
    txt = DOC.read_text()
    assert "persistently unexplored" in txt
    assert "Nothing moved except lipid peroxidation" in txt
    # The MeSH date was asserted from memory once and was wrong (2018 vs 2020).
    assert "2020-01-01" in txt and "D000079403" in txt
