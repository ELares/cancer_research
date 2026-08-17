"""Guards for the modality-class sensitivity measurement (#724).

THE CLAIM
---------
The manuscript's pharmacological:physical ratio is reported as 9.1:1 by its own
method and 17.6:1 on the census. Five partitions, each built by ONE principle
applied to both classes, give 1.32:1 to 3.93:1. The direction survives; the
order-of-magnitude does not, and the gap is a property of a PHYSICAL class that
contains three mechanism tags and excludes radiotherapy.

WHAT MAKES THIS EASY TO FAKE, AND THEREFORE WORTH GUARDING
-----------------------------------------------------------
1. ONE-SIDED WIDENING. Adding radiotherapy to PHYSICAL while leaving cytotoxic
   chemotherapy out of PHARMACOLOGICAL produces an inversion by construction.
   Every partition must move both sides, so each is required to name a
   substantial member on BOTH sides that the project's own tag-based classes
   lack.

2. THE INPUT CONTAINING THE RESULT. The partitions are committed as an input
   file. If that file also carried each panel's reported ratio, the measurement
   could quietly validate itself against a number it was handed. The loader
   rejects a partition file containing a ratio.

3. A SINGLE FLATTERING PARTITION. One definition proves nothing about class
   sensitivity -- the whole point is the spread, so several are required and
   they must not all agree by being near-identical.

4. REPORTING A REPLACEMENT NUMBER. The deliverable is that the ratio depends on
   the boundary, not a better point estimate. The report must keep saying the
   direction holds.
"""

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_modality_ratio.py"
PARTITIONS = REPO_ROOT / "analysis" / "modality-partitions.json"
MD = REPO_ROOT / "analysis" / "atlas-modality-ratio.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-modality-ratio.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _parts():
    return json.loads(PARTITIONS.read_text())


def test_the_input_does_not_carry_the_result():
    """A partition file holding its own ratio lets the measurement self-validate."""
    for name, spec in _parts().items():
        assert "ratio" not in spec and "reported_ratio" not in spec, (
            f"{name} carries a ratio in the input file")
    src = SCRIPT.read_text()
    assert '"ratio" in spec or "reported_ratio" in spec' in src, (
        "the loader no longer rejects a partition file containing a result")


def test_every_partition_moves_both_classes():
    """One-sided widening manufactures whatever inversion it finds.

    Checked by requiring each partition to contain, on BOTH sides, a modality
    the project's own tag-based classes lack: radiotherapy on the physical side,
    cytotoxic chemotherapy on the pharmacological side. A partition that adds
    radiotherapy and not chemotherapy is the exact error a reviewer killed.
    """
    for name, spec in _parts().items():
        phys = " | ".join(spec["physical"]).lower()
        pharm = " | ".join(spec["pharmacological"]).lower()
        assert "radio" in phys or "irradiat" in phys or "brachytherapy" in phys, (
            f"{name}: the physical class has no radiotherapy-family member, so "
            "it repeats the omission this analysis exists to correct")
        assert ("chemotherapy" in pharm or "antineoplastic" in pharm
                or "drug therapy" in pharm), (
            f"{name}: the pharmacological class has no cytotoxic-chemotherapy "
            "member while the physical class gained radiotherapy -- that is "
            "the one-sided widening that manufactures an inversion")


def test_several_partitions_and_they_genuinely_differ():
    """The spread is the deliverable; near-identical partitions do not give one."""
    d = _doc()["partitions"]
    assert len(d) >= 4, f"only {len(d)} partitions; the spread is not established"
    ratios = sorted(v["ratio"] for v in d.values())
    assert ratios[-1] / max(ratios[0], 1e-9) > 1.5, (
        f"all partitions land within {ratios[-1]/ratios[0]:.2f}x of each other, "
        "so this does not demonstrate sensitivity to the class boundary")
    # and they must not be the same descriptor lists wearing different names
    parts = _parts()
    sizes = {k: (len(v["pharmacological"]), len(v["physical"]))
             for k, v in parts.items()}
    assert len(set(sizes.values())) >= 3, (
        f"partition sizes are near-identical {sizes}; these are not "
        "independent definitions")


def test_the_arithmetic_is_internally_consistent():
    d = _doc()
    assert d["census"] > 4_000_000
    for k, c in d["partitions"].items():
        assert c["both"] <= min(c["pharm"], c["phys"]), (
            f"{k}: more articles in both classes than in one of them")
        assert c["pharm"] <= d["census"] and c["phys"] <= d["census"]
        expect = c["pharm"] / c["phys"] if c["phys"] else 0
        assert abs(expect - c["ratio"]) < 1e-6, f"{k}: ratio is not pharm/phys"


def test_the_report_keeps_the_direction_and_refuses_a_replacement_number():
    """The finding is that the magnitude depends on the boundary."""
    d, md = _doc(), MD.read_text()
    ratios = [v["ratio"] for v in d["partitions"].values()]
    assert all(r > 1.0 for r in ratios), (
        "a partition puts physical above pharmacological; the report claims the "
        "direction survives everywhere and would need rewriting")
    assert "direction survives and the magnitude does not" in md, (
        "the report no longer states that the sign is robust and the size is "
        "not, which is the whole finding")
    # Pinned in the SOURCE too. Fifth appearance today of a guard that reads a
    # committed artifact, which is static and therefore cannot fail when the
    # generator that produced it is edited.
    assert "direction survives and the magnitude does not" in SCRIPT.read_text(), (
        "the renderer no longer emits the direction-vs-magnitude finding")
    lo, hi = min(ratios), max(ratios)
    assert f"{lo:.2f}:1 and {hi:.2f}:1" in md, (
        "the report does not state its own measured range")


def test_the_understatement_caveat_points_the_right_way():
    """The known bias cuts against this finding, and saying so is the point."""
    md = MD.read_text()
    assert "further down" in md, (
        "the report no longer says the descriptor-only bias would move the "
        "ratio further DOWN. A caveat that happens to favour the author's "
        "conclusion has to be stated in the direction it actually runs")
    assert "12.4%" in md and "22.6%" in md, (
        "the measured qualifier-axis figures backing that caveat are gone")


def test_both_column_is_reported_not_resolved():
    """Combined-modality treatment is real; forcing one class is the arbitrary act."""
    d, md = _doc(), MD.read_text()
    assert all(v["both"] > 0 for v in d["partitions"].values()), (
        "no partition reports articles in both classes, which is implausible "
        "for combined-modality treatment and suggests exclusivity was forced")
    assert "both" in md.lower() and "rather than resolved" in md


def test_an_empty_match_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'v["phys"] == 0' in src, "the empty check no longer tests the count"
    assert "is not a finding" in src and "raise SystemExit" in src


def test_render_only_works_without_the_census():
    # Run against a COPY of the tree's artifacts. Invoking the generator here
    # rewrites the committed report, which makes the suite non-idempotent and --
    # measured the hard way -- lets a mutation sweep corrupt a committed file:
    # a mutated renderer ran through this test and its text landed on disk.
    import shutil, tempfile
    with tempfile.TemporaryDirectory() as td:
        shutil.copy2(MD, Path(td) / MD.name)
        try:
            res = subprocess.run([sys.executable, str(SCRIPT), "--render-only"],
                                 cwd=REPO_ROOT, capture_output=True, text=True)
        finally:
            shutil.copy2(Path(td) / MD.name, MD)
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
