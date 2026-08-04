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
    # The late upper bound is 0.36%. A 0.01 ceiling would let the share nearly
    # triple and still pass, which is not a guard on "stayed small".
    assert sdt["late"]["ci"][1] < 0.005, sdt["late"]["ci"]
    assert sdt["late"]["share"] < 0.004, sdt["late"]["share"]
    # And the claim must stay bounded by what the test can see.
    assert sdt["min_detectable_ratio"] > 2, (
        "the detection floor dropped; the 'no detectable acceleration' framing "
        "may now understate what this data could have shown")


def test_the_report_states_what_the_test_could_not_have_detected():
    """A flat verdict from a test with no power is not evidence of no change.

    The disjoint-interval test cannot call a rise in the sonodynamic leg below
    about 4.7x, so "flat" there means "indistinguishable from a 4.6x
    acceleration". An earlier draft of the report turned that into "the field is
    measurably not moving toward the thesis", which the data cannot support.
    """
    txt = DOC.read_text()
    d = _raw()["trajectory"]
    for leg, tr in d.items():
        if tr.get("min_detectable_ratio"):
            assert f"a {tr['min_detectable_ratio']:.1f}x change" in txt, (
                f"{leg}'s detection floor is not reported alongside its verdict")
    assert "measurably not moving toward the thesis" not in txt
    assert "no detectable acceleration" in txt


def test_the_prose_tracks_the_data_rather_than_being_hardcoded():
    """The report's narrative must be conditional on the computed verdicts.

    A previous version hardcoded "Nothing moved except lipid peroxidation" and
    computed the mover list into a variable it never used, so a different window
    choice produced a document asserting that nothing moved three lines under a
    table saying otherwise. Rebuilt from the committed JSON under a forced
    verdict; the summary line must follow the data.
    """
    import atlas_thesis_position as m

    d = _raw()
    movers = sorted(k for k, v in d["trajectory"].items() if v["moved"])
    txt = DOC.read_text()
    assert f"**Moved: {', '.join(movers)}." in txt, (
        "the report's summary line does not name the computed movers")

    # And the generator must react: re-render with a window pair that moves a
    # second leg, and the narrative has to change with it.
    per_year = {int(y): dict(c) for y, c in d["by_year"].items()}
    original = (m.EARLY_WINDOW, m.LATE_WINDOW)
    try:
        m.EARLY_WINDOW, m.LATE_WINDOW = ("2021", "2022"), ("2024", "2025")
        alt = m.trajectories({str(y): c for y, c in per_year.items()})
    finally:
        m.EARLY_WINDOW, m.LATE_WINDOW = original
    alt_movers = {k for k, v in alt.items() if v["moved"]}
    assert alt_movers != set(movers), (
        "fixture assumption gone: the alternate window no longer differs, so this "
        "cannot show the prose is conditional")


def test_the_mesh_date_is_cited_not_recalled():
    """It was asserted from memory once and was wrong (2018 rather than 2020)."""
    txt = DOC.read_text()
    assert "2020-01-01" in txt and "D000079403" in txt
    # And the explanation must not be presented as settled -- the series keeps
    # falling three years past the introduction and then reverses, which
    # indexers dropping a general term cannot produce.
    assert "REVERSES" in txt or "reverses" in txt
    assert "open question" in txt
