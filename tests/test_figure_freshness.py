"""Committed census figures must be what their generator draws.

Matplotlib embeds `/CreationDate` in a PDF, so a regenerated figure always
differed from its committed copy and a stale figure looked exactly like a fresh
one. `scripts/figure_io.py` removes the field; this checks the result.

WHAT IS CHECKED. The eight census figures, by regenerating them into a scratch
directory and comparing DRAWING SURFACES rather than file bytes -- content
streams, images and their alpha planes (stream AND dict), form xobjects
(stream AND the `/BBox`, `/Matrix`, `/Group`, `/OC`, `/Subtype` that place and
clip them), page geometry resolved through inheritance (`/MediaBox`,
`/CropBox`, `/Rotate`, `/UserUnit`), the page transparency `/Group`, the
document's optional-content configuration, and the graphics-state / pattern /
shading / annotation resources those streams reference by name. Fonts are NOT hashed -- see below. Bytes are not comparable across platforms; surfaces are, and
CI on two operating systems is the evidence.

WHAT IS NOT CHECKED. Not an exhaustive list -- an earlier version claimed
to be one and a reviewer immediately found four surfaces missing from it,
which is the shape of defect this file exists to retract. These are the
ones known and deliberate:

- **No PNG content, at all.** The eight census PNGs are checked for existence
  and nothing else. Replacing one with an unrelated image passes. PNG bytes are
  not portable, and no portable comparison is implemented.
- **The twenty non-census PDFs.** They still embed a creation date. Regenerating
  them is not a metadata-only change -- it rewrites plots and emits figures that
  were never committed -- so it is filed (#788), not done here.
- **Stale inputs.** Regenerating from committed JSON cannot notice the JSON is
  old.
- **The whole `/Font` subtree**, not only glyph outlines: an earlier version
  of this bullet said "glyph outlines", and a reviewer showed that understates
  it. Swapping the Type3 outlines for `zero` and `one` changes what a figure
  made of counts READS, but so does `/Encoding /Differences` swapping the same
  two names, and `/FontMatrix` doubles every character on the page. Anything a
  font dictionary controls is uncovered. The subtree comparison that would
  catch it is not portable across operating systems (CI failed on all eight
  figures). A generator edit cannot reach this -- font changes move glyph
  metrics and so the content stream -- but a hand-edited PDF can.
- **A correct-looking figure drawn from wrong data**, unless the difference
  reaches one of the hashed surfaces.

Each hashed surface has a positive control that changes ONLY that surface and
asserts it is seen; the flagship comparison has one that invokes it. Those
controls are the reason to believe any of the above.

Where several keys of ONE object are read through the same call, at least one
key of the group carries that control -- `/BBox`, `/Matrix` and `/Group` for
form xobjects, `/Decode` for the alpha plane -- and their siblings (`/OC`,
`/Subtype`, `/Width`, `/Height`, `/ColorSpace`, `/BitsPerComponent`,
`/ImageMask`, `/Matte`) are read by the same mechanism and are NOT separately
controlled. Said explicitly because "every surface has a control" would
otherwise be the kind of unmeasured sentence this file exists to catch.
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
                    # The DICT as well as the stream: setting `/Decode [1 0]`
                    # inverts the whole alpha plane without touching a byte of
                    # the stream, moving 57.7% of fig5c's pixels -- a larger
                    # hole than the 48% blanking this branch already blocked on.
                    #
                    # NAMED KEYS, not `xref_object`. The whole object carries
                    # `/Length` and its `/Filter` parameters -- compression
                    # artifacts, and `/Length` is exactly what forced the
                    # decompressed hash below -- so a re-save of an unchanged
                    # figure moved this surface and reported a false stale.
                    # Measured: `doc.save(deflate=True)` on fig5c moves it
                    # while nothing drawn changes.
                    for k in ("Decode", "Width", "Height", "ColorSpace",
                              "BitsPerComponent", "ImageMask", "Matte"):
                        out.append((f"smask:{k}", hashlib.sha256(
                            _resolve(doc, doc.xref_get_key(smask, k)).encode()
                        ).hexdigest()))
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
                # The DICT as well, exactly as for `/SMask` one level up, and
                # missed for the same reason -- "hash the xobject" was read as
                # "hash its stream". `/Matrix` decides WHERE and at what scale
                # a form draws and `/BBox` clips it, so neither touches a byte
                # of the stream: on `fig16c_trial_share.pdf` -- which
                # `article/drafts/v1.tex` \includegraphics at \textwidth, and
                # pdfTeX copies both keys verbatim -- `/Matrix [6 0 0 6 40 40]`
                # moved 13.06% of the page in two independent renderers and
                # `/BBox [0 0 .01 .01]` deleted every scatter marker.
                #
                # NAMED KEYS, not the whole object: `xref_object` carries
                # `/Length` and `/Filter`, which are compression artifacts and
                # not portable across platforms.
                for k in ("BBox", "Matrix", "Group", "OC", "Subtype"):
                    out.append((f"xobject:{entry[1]}:{k}", hashlib.sha256(
                        _resolve(doc, doc.xref_get_key(entry[0], k)).encode()
                    ).hexdigest()))
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
            # `/Font` IS DELIBERATELY NOT HASHED, and this is the gate's
            # sharpest limit. Hashing the font subtree does catch a glyph-
            # outline swap -- the reachable attack, where exchanging the Type3
            # outlines for `zero` and `one` changes what a figure made of
            # counts READS -- but it is NOT PORTABLE: font subsetting differs
            # between macOS and Linux, and CI failed on every figure. Hashing
            # the resource dict instead is portable and catches nothing, since
            # it lists only names and xref numbers.
            #
            # So glyph outlines are outside the gate, stated rather than
            # implied, and `test_the_font_surface_is_documented_as_uncovered`
            # keeps that honest. A generator edit cannot reach it -- font
            # changes move glyph metrics and therefore the content stream, so
            # all eight figures are caught that way; the exposure is hand-
            # edited or post-processed PDFs, which is what review is for.
            # `Annots` is a PAGE key, not a resource category. Read as
            # `Resources/Annots` it returns the literal "null" on every page
            # forever, so the surface was inert -- and it was the one hashed
            # surface with no positive control, which is why that shipped.
            out.append(("Annots", hashlib.sha256(_resolve(
                doc, doc.xref_get_key(pg.xref, "Annots")).encode()).hexdigest()))
            # Page geometry read through the RESOLVED properties, not raw page
            # keys. `xref_get_key(pg, "Rotate")` returns null when `/Rotate` is
            # set on the inherited `/Pages` node, so the page rendered rotated
            # and the surface saw nothing. And `/CropBox` -- the box viewers
            # display and `\includegraphics` selects -- was not hashed at all:
            # a committed figure could be cropped to 18% of its area, green.
            # `/Group` is the page's transparency group. Setting
            # `<</S/Transparency/CS/DeviceGray>>` greyscales fig5c's heatmap --
            # 59.18% of its pixels -- with every other surface identical.
            out.append(("Group", hashlib.sha256(_resolve(
                doc, doc.xref_get_key(pg.xref, "Group")).encode()).hexdigest()))
            for key, val in (("Rotate", pg.rotation),
                             ("MediaBox", tuple(pg.mediabox)),
                             ("CropBox", tuple(pg.cropbox)),
                             ("UserUnit", doc.xref_get_key(pg.xref, "UserUnit"))):
                out.append((key, hashlib.sha256(
                    str(val).encode()).hexdigest()))
        # OPTIONAL CONTENT, once per document. An OCG listed in the catalog's
        # `/OCProperties ... /D <</OFF[...]>>` and named by `/OC` on an xobject
        # hides that xobject entirely: applied to fig5c's two heatmap images it
        # renders the figure BLANK -- 49.33% of the page, larger than the 48%
        # alpha blanking this branch already blocked on -- while every drawing
        # surface stays byte-identical, because nothing drawn has changed.
        #
        # Renderer-dependent: MuPDF and Acrobat honour OCGs, Quartz ignores
        # them, and pdfTeX drops `/OCProperties`, so the built article is
        # probably unaffected. The COMMITTED STANDALONE FIGURE -- the artifact
        # MANIFEST.sha256 hashes and an archive would mint a DOI for -- renders
        # wrong, and that is the thing this file is about.
        out.append(("OCProperties", hashlib.sha256(_resolve(
            doc, doc.xref_get_key(doc.pdf_catalog(), "OCProperties")
        ).encode()).hexdigest()))
        return out
    finally:
        doc.close()


def _resolve(doc, val, depth: int = 0) -> str:
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
        if nums and depth < 4:
            return "|".join(
                _resolve(doc, ("xref", f"{n} 0 R"), depth + 1) for n in nums)
        return str(raw)
    if kind == "xref":
        try:
            num = int(str(raw).split()[0])
            text = doc.xref_object(num, compressed=True)
        except Exception:
            return str(raw)
        # An INDIRECT ARRAY resolves one level short otherwise. `/Annots` on
        # every census page is `10 0 R -> []`, so the array branch above --
        # which exists precisely because bare object numbers say nothing about
        # contents -- was never reached on a real figure. Adding an annotation
        # was still caught (the array text changes); mutating an existing one
        # in place would not have been.
        if depth < 4 and text.strip().startswith("[") and re.search(
                r"\d+ 0 R", text):
            return _resolve(doc, ("array", text), depth + 1)
        return text
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
    # And it must actually CALL the comparison on every produced file.
    # Slicing that loop to `produced[:0]` left 19 tests green with a wrong
    # figure on disk -- the one fix in this file that had no control.
    import ast as _ast

    fn = next(n for n in _ast.parse(Path(__file__).read_text()).body
              if isinstance(n, _ast.FunctionDef)
              and n.name == "test_the_committed_census_figures_are_what_the_generator_draws")
    # STATEMENT-LEVEL in the loop body, not anywhere in its subtree. The
    # first version of this guard asked only that the call appear somewhere
    # under the `For` node, which the node shape survives: `if False:` around
    # the call, and `try: ... except AssertionError: pass` around the body,
    # both shipped a wrong fig2c with all 23 tests green. Requiring the call
    # to be a direct child of the loop makes reaching it unconditional and
    # its failure unswallowable.
    looped = [
        node for node in _ast.walk(fn)
        if isinstance(node, _ast.For)
        and _ast.unparse(node.iter) == "produced"
        and any(isinstance(st, _ast.Expr) and isinstance(st.value, _ast.Call)
                and getattr(st.value.func, "id", "") == "_assert_matches_committed"
                for st in node.body)
    ]
    assert looped, (
        "the flagship does not iterate `produced` calling "
        "_assert_matches_committed unconditionally at the top of its loop "
        "body; slicing, renaming, or guarding that call behind `if`/`try` "
        "leaves the gate blind while every other check stays green")
    # AND `produced` must be the REGENERATED set. Nothing pinned it, so
    # repointing that one line at `FIG_DIR` turned the whole flagship into a
    # comparison of the committed figures against themselves -- with the
    # argument guard below still satisfied, because `_assert_matches_committed`
    # is unchanged and still reads two different paths.
    src = [_ast.unparse(t.value) for t in _ast.walk(fn)
           if isinstance(t, _ast.Assign)
           and [x for x in t.targets if _ast.unparse(x) == "produced"]]
    assert src and all("scratch.glob" in x for x in src), (
        f"the flagship binds `produced` to {src}; it must come from the "
        "scratch directory the generator wrote, or it compares the committed "
        "figures to themselves")
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


def test_the_font_surface_is_documented_as_uncovered(tmp_path):
    """Glyph outlines are OUTSIDE the gate, and this pins that.

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

    # NOT detected, and that is the documented limit rather than a surprise.
    # The assertion runs the other way so the day a portable font comparison
    # exists, this fails and the docstring gets corrected.
    a, b = _drawing(out), _drawing(src)
    assert a == b, (
        "the gate now detects a glyph-outline swap. That is an improvement -- "
        "update the docstring, which states glyph outlines are outside it, and "
        "invert this assertion.")
    doc_text = " ".join(_module_docstring().split())
    assert "The whole `/Font` subtree" in doc_text, (
        "the docstring no longer states that the font subtree is uncovered")
    # The WHOLE subtree, stated as such. Saying "glyph outlines" understates
    # it -- `/Encoding` and `/FontMatrix` are neither outlines nor covered --
    # and the test below is the control for the wider claim.
    for name in ("/Encoding", "/FontMatrix"):
        assert name in doc_text, (
            f"the docstring names the font limit without {name}, which "
            "understates it to a reader deciding whether to trust this gate")


@pytest.mark.parametrize("what", ["encoding", "fontmatrix"])
def test_the_font_limit_is_wider_than_glyph_outlines(tmp_path, what):
    """The `/Font` bullet said "glyph outlines"; a reviewer showed it is the
    whole dictionary. Swapping `/Encoding /Differences` for `zero` and `one`
    makes every digit in a figure made of counts read wrong (1.27% of the
    page), and doubling `/FontMatrix` doubles all its text (1.81%). Neither is
    an outline, and neither is detected.

    Asserted NOT detected, like the outline control: this documents the limit,
    and fails the day a portable font comparison closes it.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    touched = False
    try:
        fonts = doc.xref_get_key(doc[0].xref, "Resources/Font")
        assert fonts[0] == "xref", f"no font resource to perturb: {fonts}"
        num = int(str(fonts[1]).split()[0])
        for key in doc.xref_get_keys(num):
            val = doc.xref_get_key(num, key)
            if val[0] != "xref":
                continue
            fnum = int(str(val[1]).split()[0])
            if doc.xref_get_key(fnum, "FontMatrix")[0] == "null":
                continue
            if what == "fontmatrix":
                before = str(doc.xref_get_key(fnum, "FontMatrix")[1])
                doc.xref_set_key(fnum, "FontMatrix",
                                 "[0.002 0 0 0.002 0 0]")
                touched = True
            else:
                enc = doc.xref_get_key(fnum, "Encoding")
                if enc[0] == "xref":
                    target = int(str(enc[1]).split()[0])
                    path, raw = "Differences", str(
                        doc.xref_get_key(target, "Differences")[1])
                else:
                    # A DIRECT `/Encoding` dict, which is what matplotlib
                    # writes. Setting the whole dict as the value of
                    # `Differences` lands an edit and detects nothing about
                    # the surface named here, so take the ARRAY out of it.
                    target, path = fnum, "Encoding/Differences"
                    m = re.search(r"/Differences\s*(\[.*\])\s*>>\s*$",
                                  str(enc[1]), re.S)
                    if not m:
                        continue
                    raw = m.group(1)
                if "/zero" not in raw or "/one" not in raw:
                    continue
                before = str(doc.xref_get_key(fnum, "Encoding")[1])
                doc.xref_set_key(target, path,
                                 raw.replace("/zero", "/ZZZ")
                                    .replace("/one", "/zero")
                                    .replace("/ZZZ", "/one"))
                touched = True
            break
        assert touched, (
            f"fig2c no longer carries a Type3 font to perturb for {what}, so "
            "this control no longer exercises the surface it documents")
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    # The mutation must have LANDED before "not detected" means anything: an
    # edit that silently failed to apply is indistinguishable from a surface
    # the gate covers, and reads as reassurance either way.
    chk = pymupdf.open(out)
    try:
        key = "FontMatrix" if what == "fontmatrix" else "Encoding"
        assert str(chk.xref_get_key(fnum, key)[1]) != str(before), (
            f"the /{key} edit did not survive the save, so 'not detected' "
            "would be reporting an unlanded mutation as a covered surface")
    finally:
        chk.close()
    assert _drawing(out) == _drawing(src), (
        f"the gate now detects a /{what} change. That is an improvement -- "
        "narrow the docstring's font bullet and invert this assertion.")


def _forms(doc, pg):
    return [e for e in pg.get_xobjects()
            if str(doc.xref_get_key(e[0], "Subtype")[1]) == "/Form"]


@pytest.mark.parametrize("fig,key,value", [
    # `article/drafts/v1.tex` includes fig16c at \textwidth, and pdfTeX copies
    # a form xobject's `/Matrix` and `/BBox` verbatim, so this reaches the
    # manuscript. 13.06% of the page under two independent renderers.
    ("fig16c_trial_share", "Matrix", "[6 0 0 6 40 40]"),
    # Collapsing the clip box deletes every marker the form draws.
    ("fig16c_trial_share", "BBox", "[0 0 0.01 0.01]"),
    ("fig28_census_capture", "Matrix", "[3 0 0 3 10 10]"),
    # A transparency group on the FORM, the same attack as the page-level one
    # a level down: 0.70% of fig16c.
    ("fig16c_trial_share", "Group", "<</S/Transparency/CS/DeviceGray/I true>>"),
])
def test_form_xobject_placement_changes_are_detected(tmp_path, fig, key, value):
    """The dict, not just the stream -- the `/SMask` defect one level down.

    Neither key touches a byte of the stream that WAS hashed, so both shipped
    green until a reviewer tried them.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / f"{fig}.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    try:
        forms = _forms(doc, doc[0])
        assert forms, f"{fig} no longer draws through a form xobject"
        doc.xref_set_key(forms[0][0], key, value)
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    moved = {k for (k, v), (k2, v2) in zip(a, b) if v != v2}
    assert a != b, (
        f"a form xobject's /{key} can be rewritten without the gate noticing; "
        f"on {fig} that moves the drawing and ships green")
    # ISOLATION: the named surface moved and nothing else did. A control that
    # moves two surfaces proves neither, which is how the `alpha=1.0` smask
    # control passed while testing presence rather than content.
    assert moved == {f"xobject:{forms[0][1]}:{key}"}, (
        f"expected only the /{key} surface to move, got {sorted(moved)}")


def test_a_page_transparency_group_is_detected(tmp_path):
    """`/Group` greyscales fig5c's heatmap -- 59.18% of its pixels -- while
    every drawing surface stays byte-identical, because nothing drawn changed.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    try:
        assert doc.xref_get_key(doc[0].xref, "Group")[0] == "null", (
            "fig5c now carries a /Group, so this control no longer adds one")
        doc.xref_set_key(doc[0].xref, "Group",
                         "<</S/Transparency/CS/DeviceGray/I true>>")
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    moved = {k for (k, v), (k2, v2) in zip(a, b) if v != v2}
    assert moved == {"Group"}, (
        f"expected only the /Group surface to move, got {sorted(moved)}")


def test_hiding_content_with_an_optional_content_group_is_detected(tmp_path):
    """An OCG switched `/OFF` in the catalog and named by `/OC` renders fig5c
    BLANK -- 49.33% of the page, larger than the alpha blanking this branch
    already blocked on -- with every drawing surface identical.

    Renderer-dependent (MuPDF and Acrobat honour it, Quartz does not) and
    dropped by pdfTeX, so the exposure is the COMMITTED STANDALONE FIGURE,
    which is the artifact MANIFEST.sha256 hashes.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    try:
        assert doc.xref_get_key(doc.pdf_catalog(), "OCProperties")[0] == "null"
        x = doc.get_new_xref()
        doc.update_object(x, "<</Type/OCG/Name(hidden)>>")
        doc.xref_set_key(
            doc.pdf_catalog(), "OCProperties",
            f"<</OCGs[{x} 0 R]/D<</OFF[{x} 0 R]/Order[]/BaseState/ON>>>>")
        imgs = doc[0].get_images(full=True)
        assert imgs, "fig5c no longer carries the raster this control hides"
        for img in imgs:
            doc.xref_set_key(img[0], "OC", f"{x} 0 R")
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    moved = {k for (k, v), (k2, v2) in zip(a, b) if v != v2}
    assert moved == {"OCProperties"}, (
        f"expected only the optional-content surface to move, got {sorted(moved)}")


@pytest.mark.parametrize("key,value", [
    ("MediaBox", "[0 0 900 700]"),
    ("Rotate", "90"),
    # The box viewers display and `\includegraphics` selects. Cropping a
    # committed figure to 18% of its area passed.
    ("CropBox", "[100 100 300 300]"),
    ("UserUnit", "5"),
])
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


def test_inverting_the_alpha_decode_is_detected(tmp_path):
    """Control for the smask DICT, which the stream hash cannot see.

    `/Decode [1 0]` inverts the whole alpha plane without touching a byte of
    the stream: 57.7% of fig5c's pixels move. Larger than the 48% blanking
    this branch already treated as blocking, and invisible until the dict was
    hashed alongside the stream.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    touched = 0
    try:
        for pg in doc:
            for img in pg.get_images(full=True):
                if img[1]:
                    doc.xref_set_key(img[1], "Decode", "[1 0]")
                    touched += 1
        assert touched, "fig5c carries no alpha plane to invert"
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved == {"smask:Decode"}, (
        f"inverting the alpha /Decode moved {sorted(moved) or 'nothing'}; it "
        "must move the smask dict surface and only that")


def test_geometry_inherited_from_the_pages_node_is_detected(tmp_path):
    """Control for reading geometry RESOLVED rather than as raw page keys.

    `/Rotate` on the inherited `/Pages` node renders the page rotated while
    `xref_get_key(page, "Rotate")` returns null, so the surface saw nothing.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    try:
        parent = int(str(doc.xref_get_key(doc[0].xref, "Parent")[1]).split()[0])
        doc.xref_set_key(parent, "Rotate", "90")
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    moved = {k.split(":")[0] for (k, x), (_, y) in zip(a, b) if x != y}
    assert "Rotate" in moved, (
        f"an inherited /Rotate moved {sorted(moved) or 'nothing'}; the page "
        "renders rotated, so the geometry surface must see it")


def test_an_annotation_mutated_in_place_is_detected(tmp_path):
    """`_resolve` following an INDIRECT array, which no census figure exercises.

    `/Annots` is `10 0 R -> []` on all eight pages, so the array branch -- the
    one that exists because a bare object number says nothing about contents --
    was never reached on a real figure, and only the existing direct-array
    unit test covered it. ADDING an annotation was caught either way, because
    the array text changes; moving one already there was not.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    try:
        doc[0].add_rect_annot((10, 10, 100, 100))
        doc.save(tmp_path / "a.pdf", deflate=True)
    finally:
        doc.close()

    doc = pymupdf.open(tmp_path / "a.pdf")
    try:
        annot = doc[0].first_annot
        assert annot is not None, "the annotation was not written"
        # INDIRECT, as every census page has it. pymupdf writes `/Annots` back
        # as a DIRECT array, which the array branch already handled -- so the
        # first version of this control passed with the recursion removed,
        # testing the path that was never broken. Restore the real shape.
        arr = doc.get_new_xref()
        doc.update_object(arr, f"[{annot.xref} 0 R]")
        doc.xref_set_key(doc[0].xref, "Annots", f"{arr} 0 R")
        assert doc.xref_get_key(doc[0].xref, "Annots")[0] == "xref"
        doc.save(tmp_path / "a.pdf", incremental=True, encryption=0)
        before = str(doc.xref_get_key(annot.xref, "Rect")[1])
        doc.xref_set_key(annot.xref, "Rect", "[200 200 400 400]")
        assert str(doc.xref_get_key(annot.xref, "Rect")[1]) != before
        doc.save(tmp_path / "b.pdf", deflate=True)
    finally:
        doc.close()

    chk = pymupdf.open(tmp_path / "b.pdf")
    try:
        assert chk.xref_get_key(chk[0].xref, "Annots")[0] == "xref", (
            "the save flattened /Annots to a direct array, so this control "
            "exercises the branch that already worked")
    finally:
        chk.close()
    a, b = _drawing(tmp_path / "b.pdf"), _drawing(tmp_path / "a.pdf")
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved == {"Annots"}, (
        "moving an annotation already on the page is invisible: /Annots is an "
        f"indirect array and its object NUMBER did not change. Moved: {sorted(moved)}")
