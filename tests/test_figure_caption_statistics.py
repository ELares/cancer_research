"""Manuscript captions that quote a statistic must match what the figure draws.

`fig1_ferroptosis_comparison` computes a chi-squared test and PRINTS it; the
caption in `article/drafts/v1.tex` carried the printed value by hand. The
corpus was then halved -- 10,413 files to 4,830 -- and retagged, and both the
figure and the caption went on stating numbers from the superseded corpus for
five months. The figure's own denominators (`n=617`, `n=563`) match the parent
of the commit that introduced it exactly, so it was stale the day it landed.

The figure is now regenerated and the caption reads chi^2 = 38.8. This derives
that number from the corpus rather than trusting either, which is the only
arrangement in which the two cannot drift apart again: a hand-written number
beside a computed one is this repository's most repeated defect.
"""
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TEX = REPO / "article/drafts/v1.tex"
GEN = REPO / "scripts/generate_figures.py"


def _corpus_derived_stems():
    """Figure stems FIGURES.yaml declares corpus-derived from this generator."""
    import yaml

    figs = yaml.safe_load((REPO / "FIGURES.yaml").read_text())["figures"]
    return sorted(
        f["filename"] for f in figs
        if f.get("type") == "corpus-derived"
        and str(f.get("generator", {}).get("script", "")).endswith(
            "generate_figures.py"))


def _generator():
    sys.path.insert(0, str(REPO / "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import generate_figures
    return generate_figures


def _caption_statistic():
    """The chi-squared value and p the caption states, as floats."""
    text = TEX.read_text()
    m = re.search(
        r"\\caption\{Ferroptosis engagement \(\$\\chi\^2=([\d.]+)\$, "
        r"\$p=([\d.]+)\\times10\^\{(-?\d+)\}\$",
        text)
    assert m, (
        "the fig1 caption no longer states a chi-squared value in the form "
        "this test reads; if the statistic was dropped, drop this test with it")
    return float(m.group(1)), float(m.group(2)) * 10 ** int(m.group(3))


def _figure_annotation():
    """The chi-squared and p the committed fig1 actually draws."""
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            return None
    doc = pymupdf.open(REPO / "article/figures/fig1_ferroptosis_comparison.pdf")
    try:
        text = " ".join(" ".join(page.get_text().split()) for page in doc)
    finally:
        doc.close()
    m = re.search(r"\u03c7\u00b2\s*=\s*([\d.]+),\s*p\s*=\s*([\d.]+e?-?\d*)", text)
    if not m:
        return None
    return float(m.group(1)), m.group(2)


def test_the_fig1_caption_matches_the_corpus():
    g = _generator()
    from scipy import stats as sp

    data = g.classify_ferroptosis(g.load_corpus())
    modalities = ["SDT", "IRE", "HIFU", "TTFields", "Frequency"]
    sdt_ferro = data["SDT"]["ferroptosis"]
    sdt_total = data["SDT"]["total"]
    other_ferro = sum(data[m]["ferroptosis"] for m in modalities if m != "SDT")
    other_total = sum(data[m]["total"] for m in modalities if m != "SDT")
    chi2, p, _, _ = sp.chi2_contingency(
        [[sdt_ferro, sdt_total - sdt_ferro],
         [other_ferro, other_total - other_ferro]])

    stated_chi2, stated_p = _caption_statistic()
    assert round(chi2, 1) == stated_chi2, (
        f"the caption states chi^2 = {stated_chi2} and the corpus gives "
        f"{chi2:.1f}. Regenerate fig1 and update the caption together -- the "
        "previous mismatch survived a corpus halving because nothing checked")
    # TWO-SIDED, to two significant figures. The first version asserted
    # `abs(p - stated_p) < stated_p`, which holds for EVERY stated_p >= p --
    # a caption claiming `p = 1.0`, i.e. no significance at all, passed it.
    # Only an understated p could fail, which is the wrong half: a caption
    # that overstates its own p is the one that misleads a reader.
    assert f"{p:.1e}" == f"{stated_p:.1e}", (
        f"the caption states p = {stated_p:.1e} and the corpus gives "
        f"{p:.1e}. These must agree to the digits both display")
    # AND THE FIGURE MUST DISPLAY THE SAME NUMBERS. The caption sits beside a
    # figure that draws its own annotation, and a reviewer found the two
    # disagreeing in the digit both print (caption 4.8, figure 4.7) in the very
    # commit written to stop a caption drifting from its figure. Read out of
    # the committed PDF, so this compares the artifact a reader sees.
    drawn = _figure_annotation()
    assert drawn is not None, "fig1 no longer annotates its chi-squared test"
    drawn_chi2, drawn_p = drawn
    assert drawn_chi2 == stated_chi2 and drawn_p == f"{stated_p:.1e}", (
        f"the figure draws chi^2 = {drawn_chi2}, p = {drawn_p} and the caption "
        f"states chi^2 = {stated_chi2}, p = {stated_p:.1e}")


def test_the_caption_statistic_is_the_one_the_figure_draws():
    """The test above recomputes the contingency table. If the GENERATOR ever
    builds a different one, the caption would agree with this test and disagree
    with the figure -- two artifacts consistent with each other and both wrong,
    which is the shape of defect this file exists to catch."""
    src = GEN.read_text()
    body = src[src.index("def fig1_ferroptosis_comparison"):]
    body = body[:body.index("\ndef ")]
    flat = " ".join(body.split())
    # THE OTHERS ARM TOO. Pinning only the SDT arm left `other_ferro`,
    # `other_total` and the `modalities` list free -- and those decide the
    # table. Measured: rewriting the others arm to drop one modality makes the
    # figure draw chi^2 = 36.3 while the caption and this test both still say
    # 38.8, which is verbatim the "two artifacts agreeing while both diverge
    # from the figure" failure this test is here to exclude.
    for fragment in (
            'modalities = ["SDT", "IRE", "HIFU", "TTFields", "Frequency"]',
            "sdt_ferro = data[\"SDT\"][\"ferroptosis\"]",
            "sdt_total = data[\"SDT\"][\"total\"]",
            "other_ferro = sum(data[m][\"ferroptosis\"] for m in modalities if m != \"SDT\")",
            "other_total = sum(data[m][\"total\"] for m in modalities if m != \"SDT\")",
            "contingency = [[sdt_ferro, sdt_total - sdt_ferro], [other_ferro, other_total - other_ferro]]",
            "chi2, p_value, _, _ = stats.chi2_contingency(contingency)"):
        assert " ".join(fragment.split()) in flat, (
            "fig1 no longer computes its chi-squared the way the caption test "
            f"reproduces it; missing: {fragment}")


def test_the_regenerated_corpus_figures_are_what_the_generator_draws(tmp_path):
    """The corpus figures with TRACKED inputs, regenerated and compared.

    The set is DISCOVERED from `FIGURES.yaml` -- the entries whose `type` is
    `corpus-derived` and whose generator is this script -- and asserted to be
    exactly what this test compares. Hard-coding it meant `tracked[:0]` passed,
    and a new corpus-derived figure would have been added to the manuscript and
    silently skipped here.

    The generator writes 20 figure stems, 13 of which are committed. Eight of
    those thirteen read `simulations/output/`, which is gitignored -- nothing
    under it is tracked but `.gitkeep` -- so regenerating them would commit
    plots drawn from data nobody else has and CI cannot reproduce:
    fig8_simulation_by_treatment, fig10_invivo_comparison, fig11_mufa_sweep,
    fig17_damp_heatmap, fig24_hypoxia_killcurve, fig25_bliss_synergy,
    fig26_vulnerability_window, fig27_resistance_asymmetry. An earlier version
    of this docstring said "three", naming only the last three of the eight.
    Seven more stems are drawn on every run and have never been committed.
    Both sets are reported on issue #788 rather than swept in here.

    PNG CONTENT IS NOT CHECKED, only the PDFs -- the same hole the census gate
    documents, for the same reason (PNG bytes are not portable across
    platforms, and CI compares against figures authored on macOS). Measured:
    restoring a superseded PNG beside a fresh PDF passes. The PNGs are
    committed and hashed into MANIFEST.sha256, so a swap is caught there as a
    changed file, but not here as a stale figure.
    """
    import os
    import subprocess

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_figure_freshness import _drawing

    out = tmp_path / "figs"
    out.mkdir()
    res = subprocess.run([sys.executable, str(GEN)], cwd=REPO,
                         capture_output=True, text=True,
                         env={**os.environ, "MPLBACKEND": "Agg",
                              "FERRO_FIG_DIR": str(out)})
    assert res.returncode == 0, res.stderr[-800:]
    tracked = _corpus_derived_stems()
    assert tracked == sorted([
        "fig1_ferroptosis_comparison", "fig4_molecular_overlap",
        "fig6_sdt_chain_evidence", "fig12_pathway_targets",
        "fig13_gold_set_eval"]), (
        f"FIGURES.yaml now declares {tracked} as corpus-derived from this "
        "generator. Add the new figure to this comparison deliberately, or "
        "correct FIGURES.yaml -- do not let the set drift silently")
    compared = 0
    for stem in tracked:
        produced = out / f"{stem}.pdf"
        assert produced.exists(), f"the generator no longer draws {stem}"
        committed = REPO / "article/figures" / f"{stem}.pdf"
        a, b = _drawing(produced), _drawing(committed)
        if a is None:
            pytest.skip("no PDF reader available")
        assert a == b, (
            f"article/figures/{stem}.pdf is not what "
            "scripts/generate_figures.py draws. Re-run it.")
        compared += 1
    # A COUNT, because slicing the loop to `tracked[:0]` left this green. The
    # census gate learned the same lesson through six review rounds: a check
    # that never runs reports exactly what a check that passes reports.
    assert compared == len(tracked), (
        f"compared {compared} figures, not {len(tracked)}")
