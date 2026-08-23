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

AND WHAT IT DOES NOT COVER AT ALL: the twenty non-census PDFs that still
embed a creation date. TWO DIFFERENT TWENTIES meet here and an earlier draft
of this sentence merged them. Twenty non-census PDFs carry a creation date;
`generate_figures.py` writes twenty figure stems; the overlap is 13, and the
other 7 stale ones come from `generate_conceptual_diagrams.py`,
`rare_event_analysis.py` and `sim-original`.

Making them deterministic means regenerating them, and that turned out not to
be a metadata-only change: it rewrote 8 PDFs and 8 PNGs -- the PLOTS differ
from the committed ones -- and emitted 7 figures (14 files) never committed at
all. That measurement covers the 15 stems of `generate_figures.py`'s twenty
that could be drawn here, so it is a FLOOR: five need gitignored simulation
output, and of the stale-PDF twenty only 8 were exercised at all. Those
plots differ from the published ones, which is either figure drift or input
drift and is worth finding out rather than committing blind.
"""
import hashlib
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
FIG_DIR = REPO / "article/figures"
GEN = REPO / "scripts/generate_census_figures.py"


def _drawing(path):
    """Every surface the page draws from: content, images, alpha planes,
    form xobjects, and the graphics-state / pattern / shading / annotation
    resources those streams reference by NAME.

    Not the file bytes. Removing `/CreationDate` makes a PDF reproducible on
    one machine and NOT across machines -- CI proved it, failing on Linux
    against figures written on macOS. A raw sha256 compares the toolchain as
    much as the figure. (What the CI log establishes is that the difference
    lies OUTSIDE everything hashed here; it does not identify the cause, and
    this docstring does not claim to.)

    THREE SURFACES WERE MISSED IN TURN, each found only by review:

    - Content streams alone made the check vacuous for the one census figure
      that is a raster: `fig5c_mechanism_site_matrix` is drawn with `imshow`,
      so 201 mutated matrix values passed as fresh.
    - Hashing the image then missed its `/SMask`. The heatmap's alpha plane is
      a separate 2443x1848 xref, returned as `img[1]` and discarded. Blanking
      it changes 48% of the rendered page.
    - Hashing both still missed `/ExtGState`. Content streams reference alphas
      by NAME (`/A2 gs`), so changing `alpha=0.18` to `0.90` in the generator
      leaves the stream byte-identical and moves 48,937 pixels.

    Each was the same defect one level down. `/Pattern`, `/Shading` and
    `/Annots` are hashed for the same reason even though every census page has
    them empty today.
    """
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            return None
    doc = pymupdf.open(path)
    try:
        out = []
        for pg in doc:
            out.append(("content", hashlib.sha256(pg.read_contents()).hexdigest()))
            for img in pg.get_images(full=True):
                xref, smask = img[0], img[1]
                out.append(("image",
                            hashlib.sha256(doc.extract_image(xref)["image"]).hexdigest()))
                if smask:
                    # DECOMPRESSED. `xref_stream_raw` returns the zlib-encoded
                    # bytes, so a different zlib build produces a different
                    # hash for identical pixels -- CI failed on Linux for
                    # exactly that. Every other surface here is decompressed
                    # too, which is what makes the comparison portable.
                    out.append(("smask", hashlib.sha256(
                        doc.xref_stream(smask)).hexdigest()))
            # Sorted by NAME. `sorted(get_xobjects())` orders by the object
            # NUMBER, an allocation artifact: a regeneration that permutes
            # numbering with identical content would reorder the list and
            # report a false "not fresh".
            for entry in sorted(pg.get_xobjects(), key=lambda e: e[1]):
                out.append((f"xobject:{entry[1]}", hashlib.sha256(
                    doc.xref_stream(entry[0])).hexdigest()))
            # Resources the streams name rather than inline.
            for key in ("ExtGState", "Pattern", "Shading"):
                val = doc.xref_get_key(pg.xref, f"Resources/{key}")
                out.append((key, hashlib.sha256(
                    _resolve(doc, val).encode()).hexdigest()))
            # `Annots` is a PAGE key, not a resource category. Read as
            # `Resources/Annots` it returns the literal "null" on every page
            # forever, so the surface was inert -- and it was the one hashed
            # surface with no positive control, which is why that shipped.
            out.append(("Annots", hashlib.sha256(_resolve(
                doc, doc.xref_get_key(pg.xref, "Annots")).encode()).hexdigest()))
        return out
    finally:
        doc.close()


def _resolve(doc, val) -> str:
    """Flatten a resource entry, following one level of indirection.

    `xref_get_key` returns `('xref', '4 0 R')` for an indirect dict, which is
    an object NUMBER -- stable across an unchanged regeneration but useless as
    a comparison, since it says nothing about the contents. Following it is
    what makes an `/ExtGState` alpha change visible.
    """
    kind, raw = (val if isinstance(val, tuple) else ("string", str(val)))
    if kind == "array":
        # Bare object numbers say nothing about contents, which is this
        # function's own complaint about `xref`. Resolve each element.
        nums = re.findall(r"(\d+) 0 R", str(raw))
        if nums:
            return "|".join(_resolve(doc, ("xref", f"{n} 0 R")) for n in nums)
        return str(raw)
    if kind == "xref":
        try:
            num = int(str(raw).split()[0])
            return doc.xref_object(num, compressed=True)
        except Exception:
            return str(raw)
    return str(raw)


def _census_figures():
    """The figures the census generator writes, read from its source."""
    src = GEN.read_text()
    names = set(re.findall(r'FIG_DIR / f?"([\w.]+?)\.(?:pdf|png)"', src))
    names |= {m for m in re.findall(r'FIG_DIR / f"([\w]+)\.\{ext\}"', src)}
    return sorted(names)


def test_the_figure_list_is_discovered_from_the_generator():
    figs = _census_figures()
    # EXACT. A floor let two savefig calls be replaced by `pass`, shrinking
    # coverage 25% while staying green -- and the flagship check derives its
    # file list from this same discovery.
    assert len(figs) == 8, f"{len(figs)} census figures found: {figs}"
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
    doc = _module_docstring()
    assert "loads the whole corpus and gitignored simulation outputs" in doc
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
    # EXACT. A floor let a generator vanish from FIGURES.yaml -- shrinking
    # coverage -- while still passing.
    assert len(python_gens) == 4, sorted(python_gens)
    missing = sorted(g for g in python_gens if not (REPO / g).exists())
    assert not missing, (
        f"FIGURES.yaml names generators that do not exist: {missing}. The "
        "previous version dropped them silently, so they were never checked.")
    sources = {g: (REPO / g).read_text() for g in python_gens}
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
        # A substring match is satisfied by a COMMENT saying the opposite:
        # `fig7_monte_carlo_simulation` passed only because
        # `generate_figures.py` carries a note explaining that the Rust binary
        # writes it and this script does not. Require the stem to appear in a
        # savefig-ish context, and fall back to the substring only when
        # FIGURES.yaml attributes the figure to a non-Python generator.
        attributed = {f["filename"]: f["generator"].get("script", "")
                      for f in figs if isinstance(f, dict)
                      and isinstance(f.get("generator"), dict)}
        if not str(attributed.get(stem, "")).endswith(".py"):
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


def test_the_committed_census_figures_are_what_the_generator_draws():
    """Regenerate into a scratch copy and compare.

    Deliberately NOT marked `slow`. `pytest.ini` documents that marker as
    meaning "the default run uses a cheaper check of the same property", and
    there is no cheaper counterpart here -- this is the only figure-freshness
    check. The marker is inert today (nothing passes `-m`), but if anyone
    implements what the comment describes, the marker would silence the gate.

    Into a COPY: running the generator in place would rewrite the working tree
    during a test, which is how a strided sample scan once clobbered a
    committed sidecar in this repo.
    """
    import os

    # NO SKIPS. Both conditions below were `pytest.skip` and both hid a real
    # defect: the check reported green while never running, which is the exact
    # thing this file argues about elsewhere. A generator that cannot run, or
    # that ignores the output override and would therefore write into the
    # working tree, is a failure.
    with tempfile.TemporaryDirectory() as td:
        scratch = Path(td) / "figures"
        scratch.mkdir()
        env = {**os.environ, "MPLBACKEND": "Agg",
               "FERRO_FIG_DIR": str(scratch)}
        res = subprocess.run([sys.executable, str(GEN)], cwd=REPO,
                             capture_output=True, text=True, env=env)
        assert res.returncode == 0, (
            f"the census figure generator failed:\n{res.stderr[-800:]}")
        # PDFs only, and the reason is measured rather than assumed. PNGs are
        # byte-identical when regenerated on the SAME machine -- which is what
        # I checked before gating them -- and are NOT across machines: CI
        # failed on Linux against PNGs authored on macOS, because font
        # rasterisation differs. A byte comparison of a PNG in CI compares the
        # font stack, not the figure. That is the same over-generalisation this
        # branch already made once with PDF bytes.
        #
        # A PLOT change is still caught: every census PNG is drawn by the
        # same code path as its PDF sibling, verified with a real edit
        # (`alpha=0.18`->`0.90` moved both fig28.pdf and fig28.png, and the
        # PDF failed).
        #
        # BUT "nothing is lost" WOULD BE FALSE and an earlier version of this
        # comment said it. A raster-only rcParam moves every PNG while leaving
        # vector PDFs untouched: measured, `savefig.dpi` 300 -> 150 in this
        # generator moves 8 of 8 PNGs and only 1 of 8 PDFs -- fig5c, the sole
        # figure carrying a raster. So editing one line here can stale seven
        # committed PNGs with this suite green. That residual is real, it is
        # the price of a portable comparison, and the `article/figures/**` CI
        # path is what puts such a change in front of a reviewer.
        produced = sorted(scratch.glob("*.pdf"))
        assert produced, (
            "the generator wrote no PDF into FERRO_FIG_DIR, so it is either "
            "ignoring the override and writing into the working tree, or "
            f"drawing nothing. stdout:\n{res.stdout[-500:]}")
        # EXACT, not a floor of one. With `assert produced` alone, a figure
        # that was never regenerated was simply never compared -- and the
        # comment claiming this check "derives its file list from the same
        # discovery" was false.
        assert {p.stem for p in produced} == set(_census_figures()), (
            "the generator drew "
            f"{sorted(p.stem for p in produced)} but the source declares "
            f"{_census_figures()}")
        for p in produced:
            committed = FIG_DIR / p.name
            assert committed.exists(), f"{p.name} is not committed"
            a, b = _drawing(p), _drawing(committed)
            assert a is not None, "no PDF reader available"
            assert a == b, (
                f"article/figures/{p.name} is not what "
                "scripts/generate_census_figures.py draws. Re-run it.")


def _module_docstring() -> str:
    """This file's docstring, parsed.

    NOT `Path(__file__).read_text()`. Both prose guards used to read the whole
    file, which contains the assert statement doing the reading -- so each
    literal occurred exactly once, inside its own assertion, and the ENTIRE
    module docstring could be deleted with all twelve tests green. A guard
    satisfied by its own source text is the vacuous-assertion class this repo
    keeps a meta-test for.
    """
    import ast as _ast

    # WHITESPACE-NORMALISED: every claim here wraps across lines, and a raw
    # substring search misses it. That has hidden a mutation seven times in
    # this session alone.
    return " ".join(
        (_ast.get_docstring(_ast.parse(Path(__file__).read_text())) or "").split())


def test_the_uncoverable_generator_is_named_rather_than_implied():
    """`generate_figures.py` loads the corpus and gitignored simulation output,
    so no CI check can regenerate it. Saying so is the point."""
    doc = _module_docstring()
    assert "loads the whole corpus and gitignored simulation outputs" in doc
    other = REPO / "scripts/generate_figures.py"
    assert other.exists()
    src = other.read_text()
    assert "make_figures_deterministic" in src, (
        "the uncoverable generator should still write deterministic PDFs, so "
        "its output can at least be diffed by hand")


# ---------------------------------------------------------------------------
# POSITIVE CONTROLS for `_drawing()`.
#
# Its absence is why three holes shipped in a row. A reviewer gutted the
# comparison four ways -- deleting the content hash, the image hash, the
# xobject hash, and finally replacing the whole function with a constant --
# and the suite stayed green every time. A checker that cannot be shown to
# detect anything is not a checker.
#
# Each control changes exactly ONE surface and asserts it is seen.
# ---------------------------------------------------------------------------

def _one_page_pdf(tmp, draw, name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(REPO / "scripts"))
    from figure_io import make_figures_deterministic
    make_figures_deterministic()
    fig, ax = plt.subplots(figsize=(2, 2))
    draw(ax)
    out = tmp / name
    fig.savefig(out)
    plt.close(fig)
    return out


@pytest.mark.parametrize("label,a,b", [
    ("content (a line moves)",
     lambda ax: ax.plot([0, 1], [0, 1]),
     lambda ax: ax.plot([0, 1], [1, 0])),
    ("image (raster values move)",
     lambda ax: ax.imshow([[0.1, 0.2], [0.3, 0.4]]),
     lambda ax: ax.imshow([[0.9, 0.2], [0.3, 0.4]])),
    # BOTH ARMS NON-OPAQUE. The first version used alpha 0.5 vs 1.0, and
    # matplotlib emits NO SMask for a fully opaque image -- so the control
    # compared surface-list LENGTHS and tested presence, never content.
    # Replacing the smask hash with the xref NUMBER, or with the image's own
    # hash, or with a constant, all passed the full suite, and blanking
    # fig5c's alpha plane passed as fresh. That is the exact defect this
    # control exists to prevent, on the exact surface it was written for.
    ("smask (only the alpha plane moves)",
     lambda ax: ax.imshow([[[1, 0, 0, 0.50], [0, 1, 0, 0.50]],
                           [[0, 0, 1, 0.50], [1, 1, 0, 0.50]]]),
     lambda ax: ax.imshow([[[1, 0, 0, 0.75], [0, 1, 0, 0.75]],
                           [[0, 0, 1, 0.75], [1, 1, 0, 0.75]]])),
    ("extgstate (only a named alpha moves)",
     lambda ax: (ax.plot([0, 1], [0, 1]), ax.grid(alpha=0.18)),
     lambda ax: (ax.plot([0, 1], [0, 1]), ax.grid(alpha=0.90))),
    # ISOLATING. Changing the marker SHAPE also changes the content stream, so
    # it passed with xobject hashing deleted -- a control that does not isolate
    # its surface proves nothing about that surface. Marker SIZE leaves the
    # stream byte-identical (`/M0 Do`) and rewrites the glyph the name points
    # at, which is measured: content same=True, xobjects same=False.
    # `Annots` had NO control, which is why it shipped read from the wrong
    # key ("Resources/Annots", always null) and hashed a constant forever.
    # A url= annotation becomes a page-level /Annots entry; the two urls
    # differ only inside it.
    ("annots (an annotation target changes)",
     lambda ax: (ax.plot([0, 1], [0, 1]),
                 ax.annotate("x", (0.5, 0.5), url="https://example.com/a")),
     lambda ax: (ax.plot([0, 1], [0, 1]),
                 ax.annotate("x", (0.5, 0.5), url="https://example.com/b"))),
    ("xobject (marker geometry changes, stream identical)",
     lambda ax: ax.plot([0, 1], [0, 1], marker="o", markersize=6),
     lambda ax: ax.plot([0, 1], [0, 1], marker="o", markersize=12)),
])
def test_the_comparison_detects_a_change_on_each_surface(tmp_path, label, a, b):
    """If any of these stops failing, `_drawing()` has a hole on that surface."""
    pa = _one_page_pdf(tmp_path, a, "a.pdf")
    pb = _one_page_pdf(tmp_path, b, "b.pdf")
    da, db = _drawing(pa), _drawing(pb)
    assert da is not None, "no PDF reader available"
    assert da != db, (
        f"_drawing() cannot see a change to {label}, so any figure whose only "
        "difference is on that surface passes as fresh")


def test_the_comparison_is_stable_for_an_unchanged_figure(tmp_path):
    """The other half: it must not fire on an honest regeneration, or the
    controls above would be satisfied by a function that always differs."""
    draw = lambda ax: (ax.imshow([[0.1, 0.2], [0.3, 0.4]]), ax.grid(alpha=0.3))
    pa = _one_page_pdf(tmp_path, draw, "a.pdf")
    pb = _one_page_pdf(tmp_path, draw, "b.pdf")
    assert _drawing(pa) == _drawing(pb), (
        "_drawing() differs for two identical figures, so every comparison it "
        "makes is meaningless")


def test_resolve_follows_arrays_and_indirect_references(tmp_path):
    """`_resolve`'s array branch is not reached by any current figure.

    Every census page's `/Annots` resolves as a single indirect reference, so
    the array path is defence-in-depth -- and an unmeasured branch in a guard
    is the same defect as an unmeasured sentence in a report. Exercised
    directly here so it cannot rot: hashing bare object NUMBERS is exactly the
    vacuity `_resolve` exists to avoid, and two different annotations
    routinely land on the same xref number.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(REPO / "scripts"))
    from figure_io import make_figures_deterministic
    make_figures_deterministic()

    def draw(url, name):
        fig, ax = plt.subplots(figsize=(2, 2))
        ax.plot([0, 1], [0, 1])
        ax.annotate("x", (0.5, 0.5), url=url)
        out = tmp_path / name
        fig.savefig(out)
        plt.close(fig)
        return out

    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    a = draw("https://example.com/a", "a.pdf")
    b = draw("https://example.com/b", "b.pdf")
    da, db = pymupdf.open(a), pymupdf.open(b)
    try:
        ka = da.xref_get_key(da[0].xref, "Annots")
        kb = db.xref_get_key(db[0].xref, "Annots")
        assert ka[1] == kb[1], (
            "the two annotations no longer share an object number, so this "
            "control no longer demonstrates that numbers are not contents")
        assert _resolve(da, ka) != _resolve(db, kb), (
            "_resolve returns the same value for different annotations, so it "
            "is comparing object numbers rather than contents")
        # And the array form is followed rather than stringified.
        arr = ("array", f"[{ka[1]}]")
        assert _resolve(da, arr) == _resolve(da, ka), (
            "_resolve does not follow a one-element array to the same content "
            "as the bare reference")
    finally:
        da.close()
        db.close()


def _source_of(func_name: str) -> str:
    """The source of one test function, comments included."""
    import ast as _ast

    src = Path(__file__).read_text()
    tree = _ast.parse(src)
    node = next(n for n in tree.body
                if isinstance(n, _ast.FunctionDef) and n.name == func_name)
    lines = src.splitlines()
    return "\n".join(lines[node.lineno - 1:node.end_lineno])


def _regen(tmp, edit=None):
    """Run the census generator into a scratch dir, optionally with one edit."""
    src = GEN.read_text()
    out = tmp / ("mut" if edit else "base")
    out.mkdir()
    env = {**os.environ, "MPLBACKEND": "Agg", "FERRO_FIG_DIR": str(out)}
    try:
        if edit:
            GEN.write_text(edit(src))
        res = subprocess.run([sys.executable, str(GEN)], cwd=REPO,
                             capture_output=True, text=True, env=env)
    finally:
        GEN.write_text(src)
    assert GEN.read_text() == src, "the generator was not restored"
    assert res.returncode == 0, res.stderr[-500:]
    return out


def test_the_stated_png_residual_is_measured_not_asserted():
    """The docstring quotes how many PNGs a raster-only rcParam can stale.

    A hand-written number beside a measured one is the defect this whole file
    retracts, and two reviews disagreed about this particular figure (1 of 8
    PDFs versus 2 of 8) -- which is itself the argument for deriving it. The
    residual is recomputed here and the prose must match.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        base = _regen(tmp)
        mut = _regen(tmp, lambda s: s.replace('"savefig.dpi": 300',
                                              '"savefig.dpi": 150'))
        moved_png, moved_pdf, caught = [], [], []
        for p in sorted(base.glob("*")):
            q = mut / p.name
            if p.read_bytes() == q.read_bytes():
                continue
            if p.suffix == ".png":
                moved_png.append(p.stem)
            else:
                moved_pdf.append(p.stem)
                if _drawing(p) != _drawing(q):
                    caught.append(p.stem)
        residual = len(moved_png) - len(caught)

    # Restricted to the FLAGSHIP's own source, not the whole file. Searching
    # the file let this guard be satisfied by ITS OWN DOCSTRING, which mentions
    # "1 of 8 PDFs" while describing the disagreement -- the vacuous-assertion
    # class, inside the test written to prevent it.
    #
    # Whitespace-normalised and comment markers stripped, because the claim
    # sits in a wrapped block comment.
    raw = _source_of("test_the_committed_census_figures_are_what_the_generator_draws")
    doc = " ".join(raw.replace("#", " ").split())
    assert f"moves {len(moved_png)} of 8 PNGs" in doc, (
        f"{len(moved_png)} PNGs move; the comment says otherwise")
    assert f"{len(moved_pdf)} of 8 PDFs" in doc, (
        f"{len(moved_pdf)} PDFs move; the comment says otherwise")
    assert residual > 0, "no residual, so the caveat is unnecessary"
    words = {6: "six", 7: "seven", 8: "eight"}
    assert f"stale {words.get(residual, residual)} committed PNGs" in doc, (
        f"{residual} PNGs would go stale undetected; the comment says otherwise")


def test_the_produced_set_equality_can_actually_fail():
    """A control for the one fix that had none.

    The exact produced-set assertion closed "a figure that was never
    regenerated was never compared", and deleting it left the suite green --
    so the hole could return silently. This is the same reasoning that gave
    `/Annots` a control after it shipped inert.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        out = _regen(Path(td))
        produced = sorted(out.glob("*.pdf"))
        assert {p.stem for p in produced} == set(_census_figures())
        # Drop one and the comparison must no longer hold.
        produced[0].unlink()
        assert {p.stem for p in sorted(out.glob("*.pdf"))} != set(_census_figures()), (
            "removing a regenerated figure leaves the set comparison "
            "satisfied, so it cannot detect one that was never drawn")
    # AND the flagship must actually make that comparison. Testing the idea
    # in isolation left the flagship's own assertion deletable with the suite
    # green -- which is the hole this control exists for.
    body = _source_of("test_the_committed_census_figures_are_what_the_generator_draws")
    flat = " ".join(body.split())
    assert "== set(_census_figures())" in flat, (
        "the flagship no longer compares its produced set to the discovered "
        "one, so a figure that was never regenerated is never compared")
    assert "assert True or" not in flat
