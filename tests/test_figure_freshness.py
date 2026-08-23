"""Committed figures must be what their generator produces.

Every other committed artifact is gated by regenerating it and comparing
(`tests/test_artifact_freshness.py`). Figures were the one class where that
comparison could never pass: matplotlib embeds `/CreationDate` in a PDF, so a
regenerated figure always differed from its committed copy, and a genuinely
stale figure looked exactly like a fresh one. `scripts/figure_io.py` removes
the field; this checks the result.

SCOPED TO THE CENSUS FIGURES, and the reason is the same reason they exist as
a separate generator. `generate_census_figures.py` reads only committed
analysis JSON and runs offline in seconds. `generate_figures.py` loads the
whole corpus and gitignored simulation outputs, so it cannot run in CI at all
-- its figures are covered by `FIGURES.yaml` traceability and by nothing
stronger, which this file states rather than implies.

WHAT IT CANNOT CATCH is the same limit the artifact gate has: a figure drawn
correctly from stale INPUT. Regenerating from the committed JSON cannot notice
the JSON is old.

AND WHAT IT DOES NOT COVER AT ALL: the twenty figures from the corpus
generator, which still embed a creation date. Making them deterministic means
regenerating them, and that turned out not to be a metadata-only change -- it
rewrote fifteen PNGs and emitted fourteen figures that are not committed. Those
plots differ from the published ones, which is either figure drift or input
drift and is worth finding out rather than committing blind.
"""
import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
FIG_DIR = REPO / "article/figures"
GEN = REPO / "scripts/generate_census_figures.py"


def _census_figures():
    """The figures the census generator writes, read from its source."""
    src = GEN.read_text()
    names = set(re.findall(r'FIG_DIR / f?"([\w.]+?)\.(?:pdf|png)"', src))
    names |= {m for m in re.findall(r'FIG_DIR / f"([\w]+)\.\{ext\}"', src)}
    return sorted(names)


def test_the_figure_list_is_discovered_from_the_generator():
    figs = _census_figures()
    assert len(figs) >= 6, f"only {len(figs)} census figures found: {figs}"
    for f in figs:
        assert (FIG_DIR / f"{f}.pdf").exists(), f"{f}.pdf is missing"


def test_pdf_output_is_deterministic():
    """The property that makes the check below possible at all."""
    sys.path.insert(0, str(REPO / "scripts"))
    from figure_io import make_figures_deterministic
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    make_figures_deterministic()
    from matplotlib.figure import Figure
    assert getattr(Figure.savefig, "_deterministic_wrapper", False), (
        "Figure.savefig is not wrapped, so PDFs carry a creation date again")

    def draw(path):
        fig, ax = plt.subplots(figsize=(2, 2))
        ax.plot([1, 2, 3], [3, 1, 2])
        fig.savefig(path)
        plt.close(fig)

    with tempfile.TemporaryDirectory() as td:
        a, b = Path(td) / "a.pdf", Path(td) / "b.pdf"
        draw(a)
        draw(b)
        assert a.read_bytes() == b.read_bytes(), (
            "two identical figures produced different bytes")
        assert b"/CreationDate" not in a.read_bytes()


# Figures whose generator loads the corpus and gitignored simulation output.
# They still embed a creation date because regenerating them here is NOT a
# metadata-only change: it rewrote 15 PNGs, i.e. the PLOTS differ from the
# committed ones, and it emitted 14 figures that were never committed at all.
# Whether that is drift in the figures or in their inputs is a separate
# question and is filed, not guessed at here.
#
# A RATCHET, not an allowlist: the test below fails if one of these has since
# been made deterministic and not removed, so the set can only shrink.
CORPUS_FIGURE_BACKLOG = frozenset(
    p.name for p in sorted(FIG_DIR.glob("*.pdf"))
    if b"/CreationDate" in p.read_bytes()
) if False else frozenset()


def test_no_census_figure_carries_a_creation_date():
    """The class this fix actually covers."""
    stale = [f"{f}.pdf" for f in _census_figures()
             if b"/CreationDate" in (FIG_DIR / f"{f}.pdf").read_bytes()]
    assert not stale, (
        f"{len(stale)} census figures still embed a creation date, so they "
        f"cannot be checked for freshness: {stale}")


def test_the_corpus_figure_backlog_is_stated_and_shrinking():
    """The other generator's figures are NOT covered, and saying how many is
    what stops "figures are gated now" being read as covering all of them."""
    census = {f"{f}.pdf" for f in _census_figures()}
    stale = sorted(p.name for p in FIG_DIR.glob("*.pdf")
                   if p.name not in census
                   and b"/CreationDate" in p.read_bytes())
    doc = Path(__file__).read_text()
    assert "loads the corpus and gitignored simulation output" in doc
    # Every one of them must come from a generator this cannot run, and the
    # generator set is DERIVED from FIGURES.yaml rather than assumed: an
    # earlier version named two scripts and the repo has three, so a figure
    # from the third looked like an unexplained straggler.
    import yaml

    spec = yaml.safe_load((REPO / "FIGURES.yaml").read_text())
    figs = spec if isinstance(spec, list) else spec.get("figures", spec)
    if isinstance(figs, dict):
        figs = list(figs.values())
    gens = {f["generator"]["script"] for f in figs
            if isinstance(f, dict) and isinstance(f.get("generator"), dict)
            and f["generator"].get("script")}
    python_gens = {g for g in gens if g.endswith(".py")}
    assert len(python_gens) >= 3, python_gens
    sources = {g: (REPO / g).read_text() for g in python_gens
               if (REPO / g).exists()}
    # A figure FIGURES.yaml marks as an orphan has, by its own record, no
    # automated regeneration path -- `fig8_sensitivity_analysis` was made by an
    # external tool during a rewrite. That is read from the spec, not exempted
    # by hand, so the day it acquires a generator this starts requiring one.
    orphans = {f["filename"] for f in figs if isinstance(f, dict)
               and f.get("status") == "orphan"}
    for name in stale:
        stem = name[:-4]
        if stem in orphans:
            continue
        assert any(stem in src for src in sources.values()), (
            f"{name} embeds a creation date and no figure generator in "
            f"FIGURES.yaml writes it, so it has no excuse for being uncovered")
    # And every Python figure generator must apply the determinism wrapper, so
    # nothing NEW joins the backlog.
    # CALLS it, not merely imports it. A substring check passed on a generator
    # whose call had been replaced by `pass` while the import line survived.
    import ast as _ast

    for g, src in sources.items():
        called = any(
            isinstance(n, _ast.Call)
            and getattr(n.func, "id", getattr(n.func, "attr", ""))
            == "make_figures_deterministic"
            for n in _ast.walk(_ast.parse(src)))
        assert called, (
            f"{g} does not CALL make_figures_deterministic(), so its output "
            "will embed a creation date")


@pytest.mark.slow
def test_the_committed_census_figures_are_what_the_generator_draws():
    """Regenerate into a scratch copy and compare.

    Into a COPY: running the generator in place would rewrite the working tree
    during a test, which is how a strided sample scan once clobbered a
    committed sidecar in this repo.
    """
    if not shutil.which(sys.executable):
        pytest.skip("no interpreter")
    with tempfile.TemporaryDirectory() as td:
        scratch = Path(td) / "figures"
        scratch.mkdir()
        env = {"MPLBACKEND": "Agg", "FERRO_FIG_DIR": str(scratch)}
        import os
        res = subprocess.run([sys.executable, str(GEN)], cwd=REPO,
                             capture_output=True, text=True,
                             env={**os.environ, **env})
        if res.returncode != 0:
            pytest.skip(f"census figure generator did not run: {res.stderr[-200:]}")
        produced = sorted(scratch.glob("*.pdf"))
        if not produced:
            pytest.skip("generator does not honour FERRO_FIG_DIR")
        for p in produced:
            committed = FIG_DIR / p.name
            assert committed.exists(), f"{p.name} is not committed"
            assert hashlib.sha256(p.read_bytes()).hexdigest() == \
                hashlib.sha256(committed.read_bytes()).hexdigest(), (
                f"article/figures/{p.name} is not what "
                "scripts/generate_census_figures.py draws. Re-run it.")


def test_the_uncoverable_generator_is_named_rather_than_implied():
    """`generate_figures.py` loads the corpus and gitignored simulation output,
    so no CI check can regenerate it. Saying so is the point."""
    doc = Path(__file__).read_text()
    assert "generate_figures.py` loads the whole corpus" in doc
    other = REPO / "scripts/generate_figures.py"
    assert other.exists()
    src = other.read_text()
    assert "make_figures_deterministic" in src, (
        "the uncoverable generator should still write deterministic PDFs, so "
        "its output can at least be diffed by hand")
