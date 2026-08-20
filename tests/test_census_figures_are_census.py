"""The census figures must plot what their captions claim.

THE DEFECT THIS FILE EXISTS FOR was one I introduced: five figure CAPTIONS were
rewritten to describe census measurements while their generators still produced
corpus figures. That is worse than leaving both stale, because a reader checks
the caption and trusts the picture, so the disagreement is invisible from
either side alone.

The guards therefore tie three things together that can each drift alone: the
FIGURES.yaml entry, the generator function, and the committed asset.
"""
import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
FIGURES = REPO / "FIGURES.yaml"
GEN = REPO / "scripts/generate_census_figures.py"
FIG_DIR = REPO / "article/figures"
LATEX = REPO / "scripts/generate_latex.py"

# manuscript figure -> the census asset that must back it
CENSUS_FIGURES = {
    1: "fig1c_ratio_straddle",
    2: "fig2c_census_volume",
    10: "fig9c_design_composition",
    11: "fig14c_class_by_site",
    12: "fig15c_mechanism_pairs",
    13: "fig16c_trial_share",
}
# The corpus figures they replaced. If one comes back while a census caption
# stands above it, the caption/asset disagreement has returned.
SUPERSEDED = {
    "fig5_publication_trends", "fig2_mechanism_heatmap", "fig9_evidence_tiers",
    "fig14_tissue_mechanism_heatmap", "fig15_designed_combinations",
    "fig16_weighted_evidence",
}


@pytest.fixture(scope="module")
def entries():
    d = yaml.safe_load(FIGURES.read_text())
    figs = d["figures"] if isinstance(d, dict) and "figures" in d else d
    return {e.get("manuscript_figure"): e for e in figs}


def test_each_census_figure_is_registered_generated_and_present(entries):
    src = GEN.read_text()
    for num, name in CENSUS_FIGURES.items():
        e = entries.get(num)
        assert e is not None, f"manuscript figure {num} has no FIGURES.yaml entry"
        assert e["filename"] == name, (
            f"figure {num} is registered as {e['filename']}, not {name}")
        assert e["type"] == "census-derived", (
            f"figure {num} is typed {e['type']}; it reads committed census JSON")
        assert e["generator"]["script"] == "scripts/generate_census_figures.py"
        fn = e["generator"]["function"]
        assert f"def {fn}(" in src, f"{fn} is not defined in {GEN.name}"
        for ext in ("pdf", "png"):
            assert (FIG_DIR / f"{name}.{ext}").exists(), (
                f"{name}.{ext} is missing; run scripts/generate_census_figures.py")


def test_the_latex_caption_points_at_the_same_asset():
    """The caption and the \\includegraphics come from the same dict, so a
    mismatch here is what silently put a census caption over a corpus image."""
    src = LATEX.read_text()
    for num, name in CENSUS_FIGURES.items():
        m = re.search(rf"'{num}': \('([^']+)'", src)
        assert m, f"generate_latex.py has no entry for figure {num}"
        assert m.group(1) == name, (
            f"figure {num}'s LaTeX entry points at {m.group(1)} while "
            f"FIGURES.yaml says {name}")


def test_no_superseded_corpus_figure_is_still_registered(entries):
    registered = {e["filename"] for e in entries.values() if e}
    back = registered & SUPERSEDED
    assert not back, (
        f"{sorted(back)} are the corpus figures the census ones replaced; a "
        "census caption standing over one of them is the exact "
        "caption/asset disagreement these figures were built to fix")


def test_every_census_figure_reads_only_committed_json():
    """The point of a separate generator: it runs offline in seconds.

    If one of these starts reading the corpus or a census shard, the figure
    stops being reproducible by a reader who has neither.
    """
    src = GEN.read_text()
    for banned in ("corpus/atlas/records", "corpus/by-pmid", "jsonl.gz"):
        assert banned not in src, (
            f"{GEN.name} reads {banned}; these figures must read only committed "
            "analysis JSON so anyone can regenerate them")


def test_the_frozen_index_read_is_optional_and_fails_soft():
    """fig2c overlays the corpus's own volume, and the corpus is large.

    A contrast panel missing one of its two curves must say so rather than
    rendering as a complete single-series chart.
    """
    src = GEN.read_text()
    assert "_corpus_year_counts" in src
    assert "if not FROZEN_INDEX.exists():" in src
    assert "overlay is absent" in src, (
        "fig2c does not tell the reader when the comparison curve is missing")


def test_the_undetermined_class_is_drawn_not_dropped():
    """It is the largest study-design class. A chart of only the classified
    records would imply the census assigns a design to everything."""
    src = GEN.read_text()
    m = re.search(r'order = \[([^\]]+)\]', src)
    assert m and "undetermined" in m.group(1), (
        "fig9c no longer plots the undetermined class")


# --- the generated LaTeX must reference files that exist -------------------

def test_every_includegraphics_in_the_tex_resolves():
    """The check no other guard makes.

    FIGURES.yaml is checked against the figures directory, and the LaTeX figure
    map is checked against FIGURES.yaml, but nothing checked the GENERATED
    v1.tex against the filesystem. A stale .tex -- regenerated before a figure
    was renamed, or not regenerated at all -- produces a document that fails at
    build time with a missing-file error, and there is no LaTeX toolchain here
    to discover that.
    """
    import re

    tex_path = REPO / "article/drafts/v1.tex"
    tex = tex_path.read_text()
    paths = re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", tex)
    assert paths, "v1.tex contains no figures, which means it did not generate"
    missing = [p for p in paths
               if not (tex_path.parent / p).resolve().exists()]
    assert not missing, (
        f"v1.tex references figures that do not exist: {missing}. Regenerate "
        "with scripts/generate_latex.py, or the build fails on a missing file.")


def test_the_tex_is_not_stale_against_the_figure_map():
    """A .tex generated before a figure was repointed still compiles -- against
    the OLD image. That is the caption/asset disagreement at file level, and it
    is invisible to a build that succeeds."""
    import re

    tex = (REPO / "article/drafts/v1.tex").read_text()
    for num, name in CENSUS_FIGURES.items():
        assert f"{name}.pdf" in tex, (
            f"v1.tex does not include {name}.pdf, so it predates figure {num} "
            "being repointed; run scripts/generate_latex.py")
    for old in SUPERSEDED:
        assert f"{old}.pdf" not in tex, (
            f"v1.tex still includes the superseded {old}.pdf")
