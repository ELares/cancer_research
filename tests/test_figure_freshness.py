"""Committed census figures must be what their generator draws.

Matplotlib embeds `/CreationDate` in a PDF, so a regenerated figure always
differed from its committed copy and a stale figure looked exactly like a fresh
one. `scripts/figure_io.py` removes the field; this checks the result.

WHAT IS CHECKED. The eight census figures, by regenerating them into a scratch
directory and comparing DRAWING SURFACES rather than file bytes -- content
streams, images and their alpha planes, form xobjects, fonts, page geometry,
and the graphics-state / pattern / shading / annotation resources those streams
reference by name. Bytes are not comparable across platforms; surfaces are, and
CI on two operating systems is the evidence.

WHAT IS NOT CHECKED, exhaustively, because every round of review of this file
found a confident sentence outrunning the behaviour:

- **No PNG content, at all.** The eight census PNGs are checked for existence
  and nothing else. Replacing one with an unrelated image passes. PNG bytes are
  not portable, and no portable comparison is implemented.
- **The twenty non-census PDFs.** They still embed a creation date. Regenerating
  them is not a metadata-only change -- it rewrites plots and emits figures that
  were never committed -- so it is filed (#788), not done here.
- **Stale inputs.** Regenerating from committed JSON cannot notice the JSON is
  old.
- **A correct-looking figure drawn from wrong data**, unless the difference
  reaches one of the hashed surfaces.

Each hashed surface has a positive control that changes ONLY that surface and
asserts it is seen; the flagship comparison has one that invokes it. Those
controls are the reason to believe any of the above.
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
            # `/Font` is the only POPULATED resource that was skipped, while
            # three empty ones were hashed as defence-in-depth. Content
            # streams name glyphs (`/F1 ... Tj`) exactly as they name alphas
            # (`/A2 gs`), so swapping the Type3 outlines for `zero` and `one`
            # in a figure made of counts changed what it reads and passed.
            # `/MediaBox` and `/Rotate` are page geometry, cheap, same class.
            for key in ("ExtGState", "Pattern", "Shading"):
                val = doc.xref_get_key(pg.xref, f"Resources/{key}")
                out.append((key, hashlib.sha256(
                    _resolve(doc, val).encode()).hexdigest()))
            # `/Font` needs the SUBTREE, not the dict. Hashing the resource
            # object alone lists font names and xref numbers, so swapping the
            # Type3 outlines for `zero` and `one` -- which changes what a
            # figure made of counts READS -- left it unchanged.
            out.append(("Font", _subtree_digest(
                doc, doc.xref_get_key(pg.xref, "Resources/Font"))))
            # `Annots` is a PAGE key, not a resource category. Read as
            # `Resources/Annots` it returns the literal "null" on every page
            # forever, so the surface was inert -- and it was the one hashed
            # surface with no positive control, which is why that shipped.
            for key in ("Annots", "MediaBox", "Rotate"):
                out.append((key, hashlib.sha256(_resolve(
                    doc, doc.xref_get_key(pg.xref, key)).encode()).hexdigest()))
        return out
    finally:
        doc.close()


def _subtree_digest(doc, val, depth: int = 6) -> str:
    """Hash every object and stream reachable from `val`, to a bounded depth.

    One level of indirection is not enough for a font: the resource dict names
    xrefs, and the glyph outlines are streams two levels below it.
    """
    seen, parts = set(), []

    def walk(v, d):
        if d < 0:
            return
        kind, raw = (v if isinstance(v, tuple) else ("string", str(v)))
        if kind not in ("xref", "array", "dict"):
            parts.append(str(raw))
            return
        nums = re.findall(r"(\d+) 0 R", str(raw))
        if not nums:
            parts.append(str(raw))
            return
        for n in nums:
            num = int(n)
            if num in seen:
                parts.append(f"@{num}")
                continue
            seen.add(num)
            try:
                obj = doc.xref_object(num, compressed=True)
            except Exception:
                parts.append(f"?{num}")
                continue
            parts.append(obj)
            try:
                parts.append(hashlib.sha256(doc.xref_stream(num)).hexdigest())
            except Exception:
                pass
            walk(("dict", obj), d - 1)

    walk(val, depth)
    return hashlib.sha256("|".join(parts).encode()).hexdigest()


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


def _assert_matches_committed(produced: Path) -> None:
    """The flagship's comparison, extracted so a control can INVOKE it.

    Inlined, it was unguarded: deleting its assertion, or pointing `committed`
    at the produced file, left the suite green with a genuinely stale figure
    on disk. The control below calls THIS.
    """
    committed = FIG_DIR / produced.name
    assert committed.exists(), f"{produced.name} is not committed"
    a, b = _drawing(produced), _drawing(committed)
    assert a is not None, "no PDF reader available"
    assert a == b, (
        f"article/figures/{produced.name} is not what "
        "scripts/generate_census_figures.py draws. Re-run it.")


def _census_figures():
    """The figures the census generator writes, read from its source."""
    return sorted(_census_outputs()["pdf"])


def _census_outputs() -> dict:
    """Census figure stems, PER EXTENSION.

    Folding both into one name set meant a deleted `.png` savefig left the
    count unchanged, the existence check only ever looked for `{stem}.pdf`,
    and the flagship globbed `*.pdf` -- so `fig2c_census_volume.png` could
    stop being generated with the whole suite green. That is the same floor
    defect this file retracts, still open on the PNG side.
    """
    src = GEN.read_text()
    out = {"pdf": set(), "png": set()}
    for stem, ext in re.findall(r'FIG_DIR / f?"([\w.]+?)\.(pdf|png)"', src):
        out[ext].add(stem)
    # `f"{stem}.{ext}"` loops write both.
    for stem in re.findall(r'FIG_DIR / f"([\w]+)\.\{ext\}"', src):
        out["pdf"].add(stem)
        out["png"].add(stem)
    return {k: sorted(v) for k, v in out.items()}


def test_the_figure_list_is_discovered_from_the_generator():
    outputs = _census_outputs()
    assert len(outputs["png"]) == 8, (
        f"{len(outputs['png'])} census PNGs are written: {outputs['png']}. "
        "PNG stems are counted separately from PDF stems precisely so a "
        "deleted `.png` savefig cannot hide behind the PDF count.")
    for stem in outputs["png"]:
        assert (FIG_DIR / f"{stem}.png").exists(), f"{stem}.png is missing"
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
    # The scope limits that remain, checked against the docstring. An earlier
    # version asserted "generate_figures.py cannot run in CI at all", which is
    # FALSE -- it exits 0 without the simulation outputs -- so a guard was
    # requiring a false sentence to stay present. That claim is gone.
    assert "No PNG content, at all" in doc
    assert "The twenty non-census PDFs" in doc
    import yaml as _yaml

    _spec = _yaml.safe_load((REPO / "FIGURES.yaml").read_text())
    figs = _spec if isinstance(_spec, list) else _spec.get("figures", _spec)
    if isinstance(figs, dict):
        figs = list(figs.values())
    attributed = {f["filename"]: str(f["generator"].get("script", ""))
                  for f in figs if isinstance(f, dict)
                  and isinstance(f.get("generator"), dict)}
    python_gens = {g for g in attributed.values() if g.endswith(".py")}
    missing = sorted(g for g in python_gens if not (REPO / g).exists())
    assert not missing, f"FIGURES.yaml names generators that do not exist: {missing}"
    sources = {g: (REPO / g).read_text() for g in python_gens}
    assert len(python_gens) == 4, sorted(python_gens)
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
        # writes it and this script does not.
        #
        # What this DOES, precisely, because an earlier version of this comment
        # described a "savefig-ish context" check that was never written: a
        # figure FIGURES.yaml attributes to a non-Python generator is skipped
        # outright, so fig7 no longer depends on a comment matching. For the
        # rest the check is still a substring, and a comment mentioning a stem
        # would still satisfy it -- stated rather than implied.
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
        # every vector PDF untouched: measured, turning `lines.antialiased`,
        # `text.antialiased` and `patch.antialiased` off in this generator
        # moves 8 of 8 PNGs and 0 of 8 PDFs, so editing one line here can
        # stale eight committed PNGs with this suite green.
        #
        # The first version of this comment used `savefig.dpi` 300 -> 150,
        # which was the wrong example for the claim: it moves 8 PNGs and 1 PDF
        # (fig5c, the only figure with a raster), and that PDF failing turns
        # the suite RED -- so the very line offered as "green" was caught. The
        # residual is real; the example had to be one that demonstrates it.
        #
        # That residual is the price of a portable comparison, and the
        # `article/figures/**` CI path is what puts such a change in front of
        # a reviewer.
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
            _assert_matches_committed(p)


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
    """If any of these stops failing, `_drawing()` has a hole on that surface.

    SELF-VERIFYING. Each control now proves its own isolation rather than
    leaving it to a reviewer's printout: the named surface must be PRESENT in
    both arms, and must be the ONLY one that moved. Both properties were
    hand-checked before and both had already failed once -- the smask control
    used `alpha=1.0`, which emits no SMask at all, so it compared list lengths
    and tested presence rather than content; and the xobject control changed
    the marker SHAPE, which moves the content stream too. Data-only fixes to
    either are one edit from regressing, and this is what stops that.
    """
    surface = label.split()[0]
    pa = _one_page_pdf(tmp_path, a, "a.pdf")
    pb = _one_page_pdf(tmp_path, b, "b.pdf")
    da, db = _drawing(pa), _drawing(pb)
    assert da is not None, "no PDF reader available"

    def keys(d):
        return [k.split(":")[0].lower() for k, _ in d]

    assert surface in keys(da) and surface in keys(db), (
        f"the {label} control does not emit a {surface!r} surface in both "
        f"arms, so it tests PRESENCE rather than content. Arms carry "
        f"{keys(da)} and {keys(db)}")
    assert len(da) == len(db), (
        f"the {label} control changes the surface LIST, not its contents: "
        f"{keys(da)} vs {keys(db)}")
    moved = {k.split(":")[0].lower()
             for (k, x), (_, y) in zip(da, db) if x != y}
    assert moved == {surface}, (
        f"the {label} control moves {sorted(moved) or 'nothing'}; a control "
        f"that does not isolate {surface!r} proves nothing about it")


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
    # Run a COPY. Writing the tracked generator and restoring it is the hazard
    # the flagship's own docstring says this avoids -- and killing a run
    # mid-test left the edit on disk, which is exactly how a strided sample
    # scan once clobbered a committed sidecar here.
    run_from = tmp / ("gen_mut.py" if edit else "gen_base.py")
    run_from.write_text(edit(src) if edit else src)
    env = {**os.environ, "MPLBACKEND": "Agg", "FERRO_FIG_DIR": str(out),
           "PYTHONPATH": str(REPO / "scripts")}
    res = subprocess.run([sys.executable, str(run_from)], cwd=REPO,
                         capture_output=True, text=True, env=env)
    assert GEN.read_text() == src, "the tracked generator was modified"
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
        mut = _regen(tmp, lambda s: s.replace(
            "plt.rcParams.update({",
            'plt.rcParams.update({\n    "lines.antialiased": False,\n'
            '    "text.antialiased": False,\n    "patch.antialiased": False,'))
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
    # And it must compare the PRODUCED file to the COMMITTED one. Rewriting
    # that line to compare the committed figure to itself passed every other
    # check, because all six surface controls exercise `_drawing()` and none
    # exercised the call site.
    #
    # Checked by AST on the ARGUMENTS, not by matching source text: this
    # assertion's own literal lives in the same file, so a rename that touches
    # both keeps them in step and passes. Structure cannot be renamed into
    # agreement.
    import ast as _ast

    node = next(n for n in _ast.parse(Path(__file__).read_text()).body
                if isinstance(n, _ast.FunctionDef)
                and n.name == "_assert_matches_committed")
    args = [_ast.unparse(c.args[0]) for c in _ast.walk(node)
            if isinstance(c, _ast.Call)
            and getattr(c.func, "id", "") == "_drawing" and c.args]
    assert len(args) >= 2 and len(set(args)) >= 2, (
        f"the flagship calls _drawing on {args}; comparing a figure against "
        "itself always succeeds")


def test_the_flagship_comparison_can_actually_fail(tmp_path):
    """END-TO-END control for the flagship, which had none.

    Every other control exercises `_drawing()`; nothing exercised the check
    that uses it. A reviewer voided both of the flagship's assertions with the
    suite green -- comparing the committed figure to ITSELF, and relaxing the
    set equality to a subset -- and a source-text grep does not cover the
    self-compare.

    So the flagship's own comparison is run here against a perturbed COPY of a
    committed figure, on the surface this whole branch was blocked over.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    assert src.exists()
    copy = tmp_path / src.name
    copy.write_bytes(src.read_bytes())

    doc = pymupdf.open(copy)
    try:
        smasks = [img[1] for pg in doc for img in pg.get_images(full=True)
                  if img[1]]
        assert smasks, (
            "fig5c no longer carries an alpha plane, so this control no longer "
            "perturbs the surface it was written for")
        length = len(doc.xref_stream(smasks[0]))
        doc.update_stream(smasks[0], b"\xff" * length)
        # Saved to a NEW path: pymupdf refuses a non-incremental save over the
        # file it opened.
        out = tmp_path / "perturbed.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    # INVOKE the flagship's comparison, do not re-implement it. The previous
    # version asserted `_drawing(out) != _drawing(src)` inline, which made it a
    # seventh `_drawing` control rather than a control of the check that uses
    # it -- so deleting the flagship's assertion, or pointing `committed` at
    # the produced file, both stayed green with a stale figure on disk.
    staged = tmp_path / src.name
    staged.write_bytes(out.read_bytes())
    with pytest.raises(AssertionError):
        _assert_matches_committed(staged)
    # And it must PASS on the unperturbed figure, or it would "detect"
    # everything and mean nothing.
    staged.write_bytes(src.read_bytes())
    _assert_matches_committed(staged)


def test_swapping_glyph_outlines_is_detected(tmp_path):
    """`/Font` control, crafted rather than drawn.

    A font-family change is not isolating -- glyph metrics move the content
    stream too, and the isolation assertion above correctly refuses it. The
    reachable attack is the one a reviewer used: swap the Type3 outlines for
    `zero` and `one` inside a committed figure, so a page made of counts reads
    differently while every other surface is untouched. Hashing the `/Font`
    RESOURCE did not catch it; hashing the subtree does.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    swapped = False
    try:
        fonts = doc.xref_get_key(doc[0].xref, "Resources/Font")
        assert fonts[0] == "xref", f"no font resource to perturb: {fonts}"
        num = int(str(fonts[1]).split()[0])
        for key in doc.xref_get_keys(num):
            val = doc.xref_get_key(num, key)
            if val[0] != "xref":
                continue
            procs = doc.xref_get_key(int(str(val[1]).split()[0]), "CharProcs")
            if procs[0] != "xref":
                continue
            cnum = int(str(procs[1]).split()[0])
            keys = doc.xref_get_keys(cnum)
            if "zero" in keys and "one" in keys:
                z = int(str(doc.xref_get_key(cnum, "zero")[1]).split()[0])
                o = int(str(doc.xref_get_key(cnum, "one")[1]).split()[0])
                zs, os_ = doc.xref_stream(z), doc.xref_stream(o)
                doc.update_stream(z, os_)
                doc.update_stream(o, zs)
                swapped = True
                break
        assert swapped, (
            "fig2c no longer carries Type3 `zero`/`one` glyphs, so this "
            "control no longer perturbs the surface it was written for")
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    moved = {k.split(":")[0] for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved == {"Font"}, (
        f"swapping glyph outlines moved {sorted(moved) or 'nothing'}; it must "
        "move the Font surface and only that")


@pytest.mark.parametrize("key,value", [("MediaBox", "[0 0 900 700]"),
                                       ("Rotate", "90")])
def test_page_geometry_changes_are_detected(tmp_path, key, value):
    """Crafted, for the same reason as the font control: no generator edit
    moves page geometry without also moving the content stream, and a surface
    hashed without a control is how `/Annots` shipped inert."""
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    try:
        doc.xref_set_key(doc[0].xref, key, value)
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    moved = {k.split(":")[0] for (k, x), (_, y) in zip(a, b) if x != y}
    assert key in moved, (
        f"changing /{key} moved {sorted(moved) or 'nothing'}; the surface is "
        "hashed but not detected")
