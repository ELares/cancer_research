"""Guards for Figure 30 and the manuscript section that reads it.

A diagram summarising a measurement is the easiest place for a stale number to
hide: this repository has shipped that defect in prose more than once, and a
picture is harder to notice it in. So nothing on Figure 30 is typed, and the
manuscript paragraph beside it is checked against the same JSON.

The figure is DELIBERATELY a moving one — as arms land the grey shrinks — which
is exactly why the numbers have to be derived rather than captioned.
"""
import ast
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON_ = REPO / "analysis/modality-coverage.json"
DIAGRAMS = REPO / "scripts/generate_conceptual_diagrams.py"
MANUSCRIPT = REPO / "article/drafts/v1.md"
FIGURES = REPO / "article/figures"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


@pytest.fixture(scope="module")
def manuscript():
    return " ".join(MANUSCRIPT.read_text().split())


def test_the_figure_exists_in_both_formats():
    for ext in ("pdf", "png"):
        f = FIGURES / f"fig30_modality_landscape.{ext}"
        assert f.exists(), f"{f} is missing"
        assert f.stat().st_size > 5_000, f"{f} is suspiciously small"


def test_the_figure_reads_the_artifact_and_types_no_number():
    """Every quantity must come from the JSON.

    Parsed rather than eyeballed: the generator's function body must contain
    no bare integer above 1,000 (the scale every count in this figure lives
    at), because such a literal is either a typed census number or a typed
    threshold, and both go stale.
    """
    tree = ast.parse(DIAGRAMS.read_text())
    fn = next((n for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef)
               and n.name == "fig30_modality_landscape"), None)
    assert fn, "fig30_modality_landscape is gone"
    src = ast.get_source_segment(DIAGRAMS.read_text(), fn)
    assert "modality-coverage.json" in src, (
        "the figure no longer reads the artifact, so its numbers are typed")
    big = [n.value for n in ast.walk(fn)
           if isinstance(n, ast.Constant) and isinstance(n.value, int)
           and not isinstance(n.value, bool) and n.value > 1_000]
    assert not big, (
        f"the figure body contains typed integers {big}; every count on it "
        "must be read from analysis/modality-coverage.json")
    # And the legend's tier names must be the tiers the artifact actually uses.
    for tier in ("treatment", "modifier", "absent"):
        assert f'"{tier}"' in src, f"the {tier} tier is no longer drawn"


def test_the_figure_is_registered_with_its_input():
    """An unregistered figure is invisible to the freshness machinery."""
    reg = (REPO / "FIGURES.yaml").read_text()
    assert "fig30_modality_landscape" in reg
    block = reg[reg.index("fig30_modality_landscape"):]
    block = block[:block.index("\n  - ") if "\n  - " in block else len(block)]
    assert "analysis/modality-coverage.json" in block, (
        "the register does not record the figure's input, so a stale input "
        "cannot be traced to it")
    assert "generate_conceptual_diagrams.py" in block


def test_the_manuscript_section_matches_the_measurement(d, manuscript):
    """The prose beside the figure is the repo's dominant defect class.

    Every number in Section 10.2 is checked against the artifact, so the
    section fails the day an arm lands and it is not re-derived.
    """
    absent = [r for r in d["rows"] if r["engine_tier"] == "absent"]
    treat = [r for r in d["rows"] if r["engine_tier"] == "treatment"]
    assert "What This Engine Cannot Be Asked" in manuscript, (
        "Section 10.2 is gone; the criticism it answers has not gone away")
    if absent:
        a_vol = sum(r["census"] for r in absent)
        a_tr = sum(r["trials"] for r in absent)
        pct = f"{a_vol / d['total_census'] * 100:.0f}%"
        assert f"**{len(absent)} have no engine representation at all**" in manuscript
        assert f"{a_vol:,} census articles" in manuscript
        assert f"{a_tr:,} registered trials" in manuscript
        assert f"({pct} of the taxonomy's volume)" in manuscript
        assert f"{len(absent)} of sixteen mechanisms remain absent" in manuscript
    else:
        # The count reached zero, which is where the section has to change
        # what it is ABOUT rather than change a number. Presence is not
        # applicability, and the harder count must now be the stated one.
        assert "That column is now empty" in manuscript, (
            "nothing is absent any more and the section still reports an "
            "absence count; it has to move to the applicability count")
        assert "Presence is not applicability" in manuscript
        assert "no mechanism remains absent, but fifteen of sixteen remain " \
               "inapplicable" in manuscript
        assert len(treat) == 1, (
            f"{len(treat)} mechanisms are applicable and the section still "
            "says one")
        # The opening figure must still be quoted, because it is what makes
        # the closing one mean anything.
        assert "thirteen had **no engine representation at all**" in manuscript
        assert "90,019 census articles" in manuscript


def test_the_section_keeps_the_refusal_and_the_limits(manuscript):
    """It would be easy to write this section as a victory lap. These are the
    sentences that stop it being one, and each names a limit that would change
    a claim if lifted."""
    for frag in (
        "Volume is not importance, and the table is not a ranking",
        "arms is not parity",
        "not wired into any binary",
        "chemotherapy and radiotherapy are not even in the taxonomy",
    ):
        assert frag in manuscript, f"the section no longer says: {frag}"


def test_the_named_arms_are_the_arms_that_exist(d, manuscript):
    """Each arm the section describes must be findable in the engine, and the
    section must not describe one that was removed."""
    core = REPO / "simulations/ferroptosis-core/src"
    assert "Radiation" in d["treatment_variants"]
    assert (core / "radiation.rs").exists()
    params = (core / "params.rs").read_text()
    for field, label in (("parp_alpha_boost", "Synthetic lethality"),
                         ("baseline_antigenicity", "Checkpoint blockade")):
        assert field in params, f"{field} is gone but the section describes it"
        assert label in manuscript
    assert "liposomal_nanocarrier" in (core / "drug_transport.rs").read_text()
    assert "Nanocarrier delivery" in manuscript
    # The depth numbers quoted in the section are the ones the run produced.
    dr = json.loads((REPO / "analysis/depth-reach-comparison.json").read_text())
    rad = next(r for r in dr["rows"] if r["treatment"] == "Radiation")
    pdt = next(r for r in dr["rows"] if r["treatment"] == "PDT")
    assert f"{rad['surface_kill_pct']:.1f}% at the surface" in manuscript
    assert f"{rad['deep_kill_pct']:.1f}% at {rad['deep_mm']:.1f} mm" in manuscript
    assert f"{pdt['surface_kill_pct']:.1f}% to {pdt['deep_kill_pct']:.1f}%" in manuscript
