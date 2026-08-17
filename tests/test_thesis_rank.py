"""Guards for the all-modality thesis ranking (#725).

THE CLAIM
---------
Ranking the MeSH `Ferroptosis` intersection against every modality descriptor in
the committed partition universe puts the thesis's own central mechanism --
sonodynamic therapy -- at rank 22 of 166, on 32 articles. The largest
intersection is `antineoplastic agents` at 954, which is 29.8x the sonodynamic
leg, and `immunotherapy` is 370.

So the unexamined precedent is chemotherapy, not radiotherapy.

WHY THIS NEEDS GUARDING
-----------------------
1. THE UNIVERSE MUST NOT BE HAND-PICKED. The point is that a hand-written list
   cannot surface a leg nobody thought of. If this script ever builds its own
   modality list, it reproduces the defect it exists to correct.

2. THE COMPARISON MUST BE DESCRIPTOR-TO-DESCRIPTOR. The withdrawn version of
   this issue claimed radiotherapy was ~10x the sonodynamic leg, from a wide
   text-stem count against a narrow descriptor count. That reproduced under no
   definition. Both sides must come from the same axis.

3. THE WHOLE RANKING MUST BE PUBLISHED. A single pair quoted from a ranking
   nobody can see is a conclusion chosen first and evidenced afterwards.

4. UNMEASURABLE MODALITIES MUST BE NAMED, NOT ZEROED. A modality with no usable
   descriptor reported as 0 is a different and false claim.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_thesis_rank.py"
MD = REPO_ROOT / "analysis" / "atlas-thesis-rank.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-thesis-rank.json"
PARTITIONS = REPO_ROOT / "analysis" / "modality-partitions.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def test_the_universe_comes_from_the_committed_partitions():
    """A self-built list reproduces the defect this exists to correct."""
    src = SCRIPT.read_text()
    assert "d = json.loads(PARTITIONS.read_text())" in src, (
        "the universe is no longer LOADED from the partition file; naming the "
        "path in a docstring while building the list inline is the same defect "
        "wearing the right words")
    assert "modality-partitions.json" in src, (
        "the modality universe is no longer the committed partition file, so "
        "this ranks against a list written by the same person making the claim")
    assert "def modality_universe" in src
    d = _doc()
    parts = json.loads(PARTITIONS.read_text())
    expected = {x.strip().lower() for spec in parts.values()
                for side in ("pharmacological", "physical")
                for x in spec.get(side, [])}
    assert d["universe_size"] == len(expected), (
        f"the audit used {d['universe_size']} descriptors, the partition file "
        f"holds {len(expected)}")


def test_the_comparison_is_descriptor_to_descriptor():
    """The withdrawn 10x claim came from mixing a text stem with a descriptor."""
    src = SCRIPT.read_text()
    # everything counted must come from the record's `mesh` field
    assert 'r.get("mesh")' in src, "the scan no longer reads the descriptor axis"
    assert "title" not in src.split("def scan(")[1].split("def render")[0], (
        "the scan reads titles; mixing a text count with a descriptor count is "
        "exactly the asymmetry that made the withdrawn 10x claim irreproducible")


def test_the_whole_ranking_is_published():
    """A pair quoted from an unpublished ranking is a chosen conclusion."""
    d, md = _doc(), MD.read_text()
    assert len(d["intersections"]) > 50, (
        f"only {len(d['intersections'])} intersections recorded; the ranking is "
        "too short to locate a leg within it")
    assert "| rank |" in md, "the ranking table is gone"
    # the leg table must report a rank OUT OF the full ranking, not in isolation
    assert f"rank of {len(d['intersections'])}" in md, (
        "the leg table does not state the size of the ranking its ranks are "
        "drawn from")


def test_the_thesis_legs_are_located_not_asserted():
    """Their rank must come from the computed ranking."""
    d, md = _doc(), MD.read_text()
    counts = dict(d["intersections"])
    rank = {k: i + 1 for i, (k, _v) in enumerate(d["intersections"])}
    sdt = "ultrasonic therapy"
    assert sdt in counts, "the sonodynamic descriptor is not in the universe"
    assert f"| sonodynamic therapy | {sdt} | {counts[sdt]:,} | {rank[sdt]} |" in md, (
        "the sonodynamic row is not the one the ranking supports")
    # the headline multiple must be derived from the same two numbers
    top_name, top_n = d["intersections"][0]
    assert f"{top_n/counts[sdt]:.1f}x" in md, (
        "the headline multiple is not top-count over sonodynamic count")


def test_the_top_intersection_is_not_the_one_the_thesis_names():
    """If it were, this analysis would have found nothing worth reporting."""
    d = _doc()
    top = d["intersections"][0][0]
    named = {"ultrasonic therapy", "photochemotherapy", "photosensitizing agents",
             "high-intensity focused ultrasound ablation"}
    assert top not in named, (
        f"the largest intersection ({top}) is one the thesis already names, so "
        "the finding as written -- that something unexamined is larger -- does "
        "not hold and the report needs rewriting")


def test_unmeasurable_modalities_are_named_not_zeroed():
    src, md = SCRIPT.read_text(), MD.read_text()
    assert "NOT_MEASURABLE" in src and "named rather than shown as zero" in md, (
        "modalities with no usable descriptor are no longer named, so their "
        "absence reads as a measured zero")
    assert "OVER-estimate" in md and "OVER-estimate" in src, (
        "the Ultrasonic Therapy breadth caveat is gone; without it every ratio "
        "against the sonodynamic leg is quoted as if exact")


def test_attention_is_not_reported_as_endorsement():
    md = MD.read_text()
    assert "ATTENTION, not endorsement" in md, (
        "the report no longer distinguishes a co-indexed literature from "
        "evidence that a combination works")
    assert "ATTENTION, not endorsement" in SCRIPT.read_text(), (
        "the renderer no longer emits that distinction")


def test_an_empty_scan_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'if d["ferroptosis_total"] == 0:' in src
    assert "is not a finding" in src and "raise SystemExit" in src
