"""Guards for the assembled census findings page (#CENSUS-FINDINGS).

WHY THIS FILE EXISTS. `census_findings.py` assembles the campaign's
conclusions from other analyses' committed JSON, so a wrong figure here
propagates as a headline. An audit found every headline number on it carried
as a LITERAL while every other number was read: setting a source artifact's
field to nonsense left the page printing the old value, rewriting 17.6:1 to
1.6:1 in the generator passed every guard, and 15 of 16 planted mutations
survived -- including deleting a whole section.

Deriving those numbers made them CORRECT. It did not make them CHECKABLE.
That is what this file is for.
"""

import importlib.util
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "census_findings.py"
MD = REPO_ROOT / "analysis" / "census-findings.md"
LAND = REPO_ROOT / "analysis" / "atlas-landscape.json"


def _mod(path=SCRIPT, name="cf"):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_volume_figures_are_recomputed_not_read_back():
    """Recomputed HERE from the landscape artifact, so a mutation in the
    generator's own derivation cannot agree with itself.
    """
    if not LAND.exists():
        pytest.skip("landscape artifact absent")
    m = _mod()
    al = _mod(REPO_ROOT / "scripts" / "atlas_landscape.py", "al")
    rows = json.loads(LAND.read_text())["rows"]
    cen = {r["mechanism"].lower(): r for r in rows}

    def tot(keys, field="mesh_census"):
        return sum((cen.get(k) or {}).get(field) or 0 for k in keys)

    v = m._volume(json.loads(LAND.read_text()))
    assert abs(v["ratio"] - tot(al.PHARMACOLOGICAL) / tot(al.PHYSICAL)) < 1e-9
    assert abs(v["ratio_mesh_frozen"]
               - tot(al.PHARMACOLOGICAL, "mesh_frozen")
               / tot(al.PHYSICAL, "mesh_frozen")) < 1e-9
    # CAPTURE is mesh_frozen over mesh_census. `keyword_frozen` is the
    # project's own tagger and gives a different answer, which is the whole
    # reason a literal could be right and unattributable at the same time.
    cap_ph = tot(al.PHARMACOLOGICAL, "mesh_frozen") / tot(al.PHARMACOLOGICAL)
    cap_py = tot(al.PHYSICAL, "mesh_frozen") / tot(al.PHYSICAL)
    assert abs(v["oversample"] - cap_py / cap_ph) < 1e-9
    alt = tot(al.PHARMACOLOGICAL, "keyword_frozen") / tot(al.PHARMACOLOGICAL)
    assert abs(v["oversample"] - (cap_py / alt)) > 0.5, (
        "the capture field no longer distinguishes mesh_frozen from "
        "keyword_frozen, so the derivation has stopped naming its source")
    # the ratio must be pharmacological OVER physical, not the reverse
    assert v["ratio"] > 1
    for k in ("precise", "precise_alt", "criterion_restored"):
        assert v[k] > 0, f"{k} was not computed"


def test_all_three_restrictions_are_reported_not_only_the_one_that_inverts():
    """An earlier version quoted only the restriction that falls below the
    manuscript's figure, out of the three the sibling page publishes. Picking
    the reading that supports the conclusion is the defect this campaign keeps
    finding, and it does not stop being one when the conclusion is
    self-critical.
    """
    if not LAND.exists():
        pytest.skip("landscape artifact absent")
    m = _mod()
    # the renderer hard-wraps, so compare against whitespace-normalised text
    md = " ".join(MD.read_text().split())
    v = m._volume(json.loads(LAND.read_text()))
    shown = [v["precise"], v["precise_alt"], v["criterion_restored"]]
    for x in shown:
        assert f"{x:.2f} : 1" in md, (
            f"restriction {x:.2f}:1 is computed and not rendered")
    below = [x for x in shown if x < m.MANUSCRIPT_RATIO]
    words = {1: "one", 2: "two", 3: "three"}
    assert f"{words.get(len(below), len(below))} readings of three" in md, (
        f"{len(below)} of {len(shown)} restrictions fall below the "
        f"manuscript's {m.MANUSCRIPT_RATIO}:1 and the page does not say so")
    assert "quoted only" in md and "the one that inverts" in md


def test_the_two_factors_are_not_welded_into_one_derivation():
    """3.3x is the frozen->census step; 1.9x is the net against the
    manuscript's KEYWORD figure. They differ because MeSH labelling alone
    moves the ratio the other way, so `because` was a false derivation.
    """
    if not LAND.exists():
        pytest.skip("landscape artifact absent")
    m = _mod()
    md = " ".join(MD.read_text().split())
    v = m._volume(json.loads(LAND.read_text()))
    net = v["ratio"] / m.MANUSCRIPT_RATIO
    label = v["ratio_mesh_frozen"] / m.MANUSCRIPT_RATIO
    if abs(net - v["oversample"]) > 0.05 * v["oversample"]:
        assert "NOT ON COMPARABLE LABELS" in md, (
            f"the net factor {net:.2f}x and the over-sampling factor "
            f"{v['oversample']:.2f}x differ, so one does not explain the "
            "other, and the page presents them as though it did")
        assert f"{v['ratio_mesh_frozen']:.1f} : 1" in md, (
            "the MeSH-labelled frozen ratio is what separates the two factors "
            "and is not shown")
        assert f"{label:.2f}x" in md


def test_the_committed_page_is_what_the_generator_produces():
    """No freshness gate existed, so every renderer edit was invisible."""
    src = SCRIPT.read_text()
    assert "def _volume(" in src
    # the page is assembled from many artifacts; re-render and compare
    import subprocess
    import sys
    import shutil
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        shutil.copy2(MD, Path(td) / MD.name)
        try:
            r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO_ROOT,
                               capture_output=True, text=True)
            fresh = MD.read_text()
        finally:
            shutil.copy2(Path(td) / MD.name, MD)
    assert r.returncode == 0, f"{r.stdout}\n{r.stderr}"
    assert fresh == MD.read_text(), (
        "analysis/census-findings.md is not what the generator produces -- "
        "re-run `python scripts/census_findings.py`")


def test_no_headline_figure_is_a_literal():
    """The failure this file exists for: a number typed into the renderer
    cannot go stale, cannot be wrong about its own provenance, and cannot be
    checked.
    """
    src = SCRIPT.read_text()
    for bad in ("**17.6 : 1**", "**9.1 : 1**", "3.3x --", "6 to 2,792",
                "**6.43 : 1**", "runs from 6 to"):
        assert bad not in src, (
            f"{bad!r} is typed into the generator; derive it from the "
            "artifact it came from")
