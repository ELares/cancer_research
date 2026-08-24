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
    # Order of magnitude, not the mantissa: the caption rounds and this test
    # should not fail on a rounding convention.
    assert abs(p - stated_p) < stated_p, (
        f"the caption states p = {stated_p:.2e} and the corpus gives {p:.2e}")


def test_the_caption_statistic_is_the_one_the_figure_draws():
    """The test above recomputes the contingency table. If the GENERATOR ever
    builds a different one, the caption would agree with this test and disagree
    with the figure -- two artifacts consistent with each other and both wrong,
    which is the shape of defect this file exists to catch."""
    src = GEN.read_text()
    body = src[src.index("def fig1_ferroptosis_comparison"):]
    body = body[:body.index("\ndef ")]
    flat = " ".join(body.split())
    for fragment in (
            "sdt_ferro = data[\"SDT\"][\"ferroptosis\"]",
            "sdt_total = data[\"SDT\"][\"total\"]",
            "contingency = [[sdt_ferro, sdt_total - sdt_ferro], [other_ferro, other_total - other_ferro]]",
            "chi2, p_value, _, _ = stats.chi2_contingency(contingency)"):
        assert " ".join(fragment.split()) in flat, (
            "fig1 no longer computes its chi-squared the way the caption test "
            f"reproduces it; missing: {fragment}")


def test_the_regenerated_corpus_figures_are_what_the_generator_draws(tmp_path):
    """The five corpus figures with TRACKED inputs, regenerated and compared.

    Only these five. The generator also writes three figures from
    `simulations/output/`, which is gitignored -- regenerating those would
    commit plots drawn from data nobody else has and CI cannot reproduce -- and
    seven that are not committed at all. Both sets are measured and reported in
    issue #788 rather than swept into this comparison.
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
    tracked = ["fig1_ferroptosis_comparison", "fig4_molecular_overlap",
               "fig6_sdt_chain_evidence", "fig12_pathway_targets",
               "fig13_gold_set_eval"]
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
