"""Committed census figures must be what their generator draws.

Matplotlib embeds `/CreationDate` in a PDF, so a regenerated figure always
differed from its committed copy and a stale figure looked exactly like a fresh
one. `scripts/figure_io.py` removes the field; this checks the result.

WHAT IS CHECKED. The eight census figures, by regenerating them into a scratch
directory and comparing DRAWING SURFACES rather than file bytes -- content
streams; images, through their DECODED PIXELS plus the named dict keys that
decoding does not carry (`/OC` above all, which hides the image entirely);
their alpha planes (stream AND named dict keys); form xobjects (stream AND the
`/BBox`, `/Matrix`, `/Group`, `/OC`, `/Subtype` that place and clip them); page
geometry resolved through inheritance (`/MediaBox`, `/CropBox`, `/Rotate`,
`/UserUnit`); the page transparency `/Group`; the document's optional-content
configuration; and the graphics-state / pattern / shading / annotation
resources those streams name -- followed to their CONTENTS, streams included,
not compared as names and object numbers, and per NAME rather than as a set.
Following stops at `_MAX_DEPTH` (24) links, which is stated because "followed
to their contents" otherwise reads as unbounded and a limit that is not
written down is a limit nobody checks. Fonts are NOT hashed -- see below. Bytes are not comparable across platforms; surfaces are, and
CI is the evidence, with one
qualification worth stating: `test-linux` runs on every push and `test-macos`
is `if: github.event_name == 'schedule'`, so the two never run on the same
event and no pull request is checked on macOS. The per-PR signal is Linux
against figures written on macOS -- which is the cross-platform comparison
that matters here -- and the macOS lane is a weekly confirmation, not a gate.

WHAT IS NOT CHECKED. Not an exhaustive list -- an earlier version claimed
to be one and a reviewer immediately found four surfaces missing from it,
which is the shape of defect this file exists to retract. These are the
ones known and deliberate:

- **No PNG content, at all.** The eight census PNGs are checked for existence
  and nothing else. Replacing one with an unrelated image passes. PNG bytes are
  not portable, and no portable comparison is implemented.
- **Most non-census PDFs.** 35 committed figures come from other
  generators. Five of them -- the corpus-derived ones, whose inputs are tracked
  -- are regenerated and gated by `tests/test_figure_caption_statistics.py`;
  the other 30 are not, and 22 of those 30 `FIGURES.yaml` marks
  `type: simulation`. SIX of them are different in kind:
  `fig30_modality_landscape`, `fig31_modality_panel`, `fig32_modality_tme` and
  `fig33_adoptive_barriers`, `fig34_depth_reach`, `fig35_calibration_verdicts`,
  `fig36_fractionation`, `fig37_chemotherapy`, `fig38_checkpoint`,
  `fig39_adoptive_escalation`, `fig40_oncolytic_bind`,
  `fig41_adc_loading` and `fig42_ablation_sleeve` read COMMITTED
  artifacts (`analysis/modality-tme.json`, `analysis/modality-coverage.json`,
  `analysis/depth-reach-comparison.json`, `analysis/modality-calibration.json`,
  `analysis/modality-panel.json`), so unlike the rest of the backlog they
  regenerate offline and every number on them is pinned by
  `tests/test_modality_landscape_figure.py`, `tests/test_modality_panel.py`
  and `tests/test_chapter6_figures.py` -- the last of which reads numbers,
  positions, drawn colours and row pairings back out of the rendered PDFs.
  THIS SENTENCE HAS BEEN FALSE TWICE. It first claimed all four figures were
  number-pinned while no test named fig32 or fig33; the retraction was written
  for those two and the SAME claim was re-asserted for fig34 and fig35 in the
  same edit, and a reviewer then reversed fig34's order, halved its stated
  depth, rotated fig35's verdicts and inflated its title counts with all
  sixteen guards green. Both are pinned now. What is pinned is cell identity,
  sign, visible refusal text, colour direction, row pairing and the drawn
  order -- not every property a figure has: a reviewer defeated the first version of those
  checks ten ways, including reversing every row so each value landed on the
  wrong arm, because they collected drawn numbers into a SET. What is pinned
  is cell identity, sign, visible refusal text and colour direction -- not
  every property a figure has. That
  sentence was FALSE when written: it claimed all four were pinned while no
  test named fig32 or fig33 at all, and multiplying every drawn cell in fig32
  by a thousand left the suite green. The PDFs are still not byte-gated here. Eight of the 22 are drawn by `generate_figures.py` from
  `simulations/output/`, which is gitignored, so CI cannot regenerate them at
  all until #788's tracking decision is made; they were regenerated in the #790
  pass and no longer embed a creation date, which is necessary for a freshness
  check and not sufficient for one. The other 14 are not blocked that way, and the first
  version of this bullet wrongly said they were:
  `fig29_rare_event_resolution` reads the TRACKED `analysis/rare-event-sweep.jsonl`
  through a deterministic generator, so regenerating it once -- which it
  needs anyway, being one of the 14 below -- would bring it under a
  gate; `fig31_modality_panel` reads the TRACKED `analysis/modality-panel.json` and
  is one of five `type: simulation` figures (fig30_modality_landscape is `type: conceptual`) whose every number is pinned by a
  test (this bullet said "the only" while the bullet above it said four,
  contradicting itself in one docstring); and `fig32_modality_tme`, `fig33_adoptive_barriers`, `fig34_depth_reach` and
  `fig35_calibration_verdicts` read TRACKED analysis JSON through the same
  offline generator, and `tests/test_chapter6_figures.py` reads their
  numbers back out of the rendered PDFs; and `fig7_monte_carlo_simulation`
  reads the TRACKED
  `simulations/simulation_results.json`, but has no committed generator at all
  -- it is a matplotlib figure whose producing code is not in this repository.
  (It is skipped below by the non-Python-generator branch, not by the orphan
  exemption: `FIGURES.yaml` gives its status as `manuscript` and names a `.rs`
  file as its generator. An earlier version of this sentence said orphan,
  contradicting the code comment a few lines under it.) (An earlier version of
  this bullet said fig7 "is drawn by a Rust binary". It is not: `sim-original`
  prints JSON and links no plotting crate, and the committed PDF's `/Producer`
  is Matplotlib. The same false claim sits in a comment in
  `scripts/generate_figures.py` and is corrected there too.) Across all seventeen, two still embed a
  creation date and so cannot be compared at all: `fig7_monte_carlo_simulation`
  and `fig8_sensitivity_analysis`, which are precisely the two whose generator
  is not in this repository. Every figure that HAS a committed generator is now
  comparable.

  That count was seven until the detector was fixed, and the true figure was
  nine. The check was `b"/CreationDate" in blob`, and cairo -- which `dot`
  renders through -- writes the Info dictionary into a compressed object
  stream, so the literal never appears while the date is plainly there. The
  two it could not see were `fig19_immune_coupling_flow` and
  `fig22_decision_flowchart`, the only two drawn by graphviz, and therefore the
  only two the savefig wrapper could never have cleaned. A byte scan looking
  for the artifacts a wrapper cannot reach was blind in exactly that place.
- **Stale inputs.** Regenerating from committed JSON cannot notice the JSON is
  old.
- **THE GENERATOR ITSELF IS TRUSTED.** Everything here compares what
  `scripts/generate_census_figures.py` writes against what is committed, so a
  generator that draws honestly and then copies THE EIGHT CENSUS FIGURES from
  `article/figures` over its own output passes every test with all eight
  wrong. This is not closeable by comparing contents: on the authoring machine
  an honest regeneration and a copy are the same bytes by design, which is the
  whole point of removing `/CreationDate`. What stands between that and a merge
  is the `scripts/**` CI path filter putting the diff in front of a reviewer.

  The EIGHT is load-bearing and the first version of this bullet got it wrong,
  offering `article/figures/*.pdf` as the example -- that glob is 30 files, so
  it trips the produced-set equality, which is the very sentence this bullet
  ends with. Two weaker versions ARE caught that way: copying every PDF, and
  writing past `FERRO_FIG_DIR`, which trips `assert produced`. No test count is
  quoted here on purpose: the first version said "57 green" and the commit that
  wrote it added two tests.
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
- **A form xobject's own resource categories, TODAY.** They are hashed and
  controlled, but no census figure's form carries a `/Resources` of its own, so
  those surfaces hash `"null"` on both arms of every real comparison. That is
  the state `/Pattern`, `/Shading` and `/Annots` were in when each of them
  turned out to be inert, and it is written down here for the same reason.

Each hashed surface has a positive control that changes ONLY that surface and
asserts it is seen; the flagship comparison has one that invokes it. Those
controls are the reason to believe any of the above.

Where several keys of ONE object are read through the same call, at least one
key of the group carries that control -- `/BBox`, `/Matrix` and `/Group` for
form xobjects, `/Decode` for the alpha plane, `/OC` for an image -- and their
siblings (a form's own `/OC` and `/Subtype`; an IMAGE's `/Decode`, `/Width`,
`/Height`, `/ColorSpace`, `/BitsPerComponent`, `/ImageMask`, `/Mask`,
`/Interpolate`; and the alpha plane's `/Width`, `/Height`, `/ColorSpace`,
`/BitsPerComponent`, `/ImageMask`, `/Matte`) are read
by the same mechanism and are NOT separately controlled -- meaning each is an
independently deletable element of a literal tuple, and deleting any one of
them leaves this file green. "Read by the same mechanism" is a description of
how they are hashed, NOT an argument that a control on one covers another; the
list is a disclosure, and it is worth exactly that. A form's `/OC` was missing
from it until a reviewer removed the key and nothing failed, which is the
defect the list exists to prevent, inside the list. Said explicitly because "every surface has a
control" would otherwise be the kind of unmeasured sentence this file exists to
catch -- and it WAS one: a reviewer found three surfaces with no control at all
(the form xobject STREAM, `/Pattern`, `/Shading`), because the isolation
assertion collapsed `xobject:M0` and `xobject:M0:BBox` to one token and the
`/BBox` movement alone satisfied it. Each of the three could be replaced by a
constant with the file green. The isolation assertion now compares FULL keys
and requires the surface's own hash to move, not just a key hanging off it.
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
    """Every surface the page draws from: content, images and their dicts,
    alpha planes, form xobjects, and the CONTENTS of the graphics-state /
    pattern / shading / annotation resources those streams reference by name.

    "By NAME" was the earlier wording and was the bug. `_resolve` followed one
    indirection and stopped, so `/Pattern` hashed as `<</H1 R>>` -- a name and
    an object number, saying exactly as little about the pattern as `4 0 R`
    said about the ExtGState dict, which is the complaint the function was
    written for. A `hatch=` or a `hatch.linewidth` rcParam moved 7-13% of a
    page and a gouraud mesh's DATA moved 59%, all green.

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

    Each was the same defect one level down, and the pattern repeated after
    that: hashing `/Pattern` and `/Shading` at all still compared only the
    NAMES inside them, and an annotation's `/AP` -- the stream a renderer
    actually draws -- sat one reference beyond where `_resolve` stopped.
    `/Pattern`, `/Shading` and `/Annots` are empty on every census page today,
    which is exactly why those two shipped: an empty resource looks identical
    whether it is followed or not.
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
                out.append((f"image:{xref_name(doc, pg, xref)}",
                            hashlib.sha256(doc.extract_image(xref)["image"]).hexdigest()))
                # THE IMAGE DICT, which was not hashed at all -- "images and
                # their alpha planes (stream AND dict)" described the alpha
                # plane and was read as covering both. Most of the dict's
                # semantics survive into `extract_image`'s decoded PNG, so
                # `/Decode`, `/ColorSpace` and `/Mask` were caught through the
                # pixels; `/OC` is not, and hides the whole image.
                #
                # An image is NOT in `get_xobjects()` -- that returns form
                # xobjects only -- which is why the named-key loop below it
                # never reached these.
                for k in ("OC", "Decode", "Mask", "ColorSpace", "ImageMask",
                          "Width", "Height", "BitsPerComponent", "Interpolate"):
                    out.append((f"image:{xref_name(doc, pg, xref)}#{k}",
                                hashlib.sha256(_resolve(
                                    doc, doc.xref_get_key(xref, k)).encode()).hexdigest()))
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
                        out.append((f"smask#{k}", hashlib.sha256(
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
                    out.append((f"xobject:{entry[1]}#{k}", hashlib.sha256(
                        _resolve(doc, doc.xref_get_key(entry[0], k)).encode()
                    ).hexdigest()))
                # AND THE FORM'S OWN RESOURCES. Only the PAGE's were read, so
                # an alpha inside a form's own `/ExtGState` moved every scatter
                # marker on fig16c with the form's stream byte-identical and no
                # surface moving. Nested form STREAMS were already covered --
                # `get_xobjects()` recurses -- which made the gap exactly the
                # non-XObject categories hanging off a form.
                #
                # OUTSIDE the key loop. It was indented one level too far, so
                # every form's resources were walked five times per form and
                # emitted five identical surfaces -- equal on both arms, so
                # harmless to the verdict and five times the work.
                out.extend(_resource_surfaces(
                    doc, _materialize(doc, doc.xref_get_key(entry[0], "Resources")),
                    f"xobject:{entry[1]}:"))
            # Resources the streams name rather than inline.
            # `/Font` is the only POPULATED resource that was skipped, while
            # three empty ones were hashed as defence-in-depth. Content
            # streams name glyphs (`/F1 ... Tj`) exactly as they name alphas
            # (`/A2 gs`), so swapping the Type3 outlines for `zero` and `one`
            # in a figure made of counts changed what it reads and passed.
            # `/MediaBox` and `/Rotate` are page geometry, cheap, same class.
            out.extend(_resource_surfaces(doc, _effective_resources(doc, pg), ""))
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
        # OPTIONAL CONTENT, once per document. This is HALF the mechanism --
        # the catalog's configuration; the other half is `/OC` on the object it
        # hides, which is hashed with the image and the form xobject above and
        # was NOT hashed when this surface was added. The first version of the
        # control here wrote both and asserted only this one moved, which is a
        # statement that the `/OC` writes moved nothing: deleting them left it
        # green while the rendered difference fell from 49.21% to 0.00%.
        #
        # An OCG listed in the catalog's `/OCProperties ... /D <</OFF[...]>>`
        # and named by `/OC` hides what it names entirely: applied to fig5c's
        # two heatmap images it
        # renders the figure BLANK -- 49.21% of the page, larger than the 48%
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


# The resource categories a drawing can hide in, MINUS `/Font`, which is
# excluded on purpose and portably (see the long note in `_drawing`), and minus
# `/XObject`, which `get_xobjects()` already walks including nested forms.
_RESOURCE_CATEGORIES = ("ExtGState", "Pattern", "Shading", "ColorSpace",
                        "Properties")

# How far the `/Parent` chain is walked looking for an inherited `/Resources`.
# Stated, and exercised by a control that hoists it 34 nodes up, because the
# previous bound was neither: lowering it from 32 to 2 changed nothing anyone
# could observe.
_MAX_PAGE_TREE = 128


def _indirect(val):
    """The object number behind an indirect reference, or None."""
    if isinstance(val, tuple) and val[0] == "xref":
        try:
            return int(str(val[1]).split()[0])
        except Exception:
            return None
    return None


def _as_text(val) -> str:
    """A key's value as PDF source text, whatever shape pymupdf returned."""
    if isinstance(val, tuple):
        return str(val[1])
    return str(val)


def _materialize(doc, val):
    """An object NUMBER for a value, whether it is indirect or a direct dict.

    A direct `<<...>>` is written into a scratch object of the open document
    -- never saved -- so that every read below goes through pymupdf's own
    parser. The first version of this walked the dictionary TEXT with a
    hand-written scanner, and a reviewer defeated it twice in the ways a
    hand-written PDF scanner is always defeated:

    - `_balanced` counted parentheses without honouring `\\)`, so a literal
      string `(a\\) /ExtGState 99 0 R b)` ended early and the scanner resumed
      INSIDE the string, reading a decoy `/ExtGState` as a real entry. With
      `dict(...)` letting the later key win, a committed fig9c rendering 14.16%
      differently passed all 53 tests.
    - A bare-name value (`/Cs1 /DeviceRGB`, the ordinary spelling for a
      colour space) parsed as an empty value plus a phantom entry named
      `DeviceRGB`, so swapping which name mapped to which space compared equal.

    Both are properties of the SCANNER, not of the PDF, and neither exists if
    the object is handed to the parser that MuPDF already ships.
    """
    num = _indirect(val)
    if num is not None:
        return num
    text = _as_text(val).strip()
    if not text.startswith("<<"):
        return None
    try:
        num = doc.get_new_xref()
        doc.update_object(num, text)
        return num
    except Exception:
        return None


def _effective_resources(doc, pg):
    """`/Resources` RESOLVED THROUGH INHERITANCE, as an object number.

    `/Resources` is one of the four inheritable page attributes, alongside
    `/MediaBox`, `/CropBox` and `/Rotate` -- and geometry was already moved to
    the resolved properties for exactly this reason, while resources were not.
    Hoisting `/Resources` onto the `/Pages` node leaves the page rendering
    identically and made every category below read `null` FOREVER: measured,
    53.01% of fig5c moved with no surface moving, which is simultaneously the
    inherited-`/Rotate` defect and the inert-`Resources/Annots` defect.

    A DIRECT `/Resources` dict counts as present. Looking only for an object
    number treated an inlined one -- legal, identical rendering, and what
    several optimisers emit -- as absent, which was a third door to the same
    62.42%.
    """
    node, depth = pg.xref, 0
    while node and depth < _MAX_PAGE_TREE:
        val = doc.xref_get_key(node, "Resources")
        if _as_text(val).strip() not in ("null", ""):
            return _materialize(doc, val)
        node, depth = _indirect(doc.xref_get_key(node, "Parent")), depth + 1
    return None


def _resource_surfaces(doc, res, prefix: str):
    """One surface per NAME in each resource category.

    PER NAME, because a positional join compares the multiset of reachable
    contents and not the mapping: repointing `/A2` at a graphics state ALREADY
    referenced by `/A3` produced a byte-identical string and moved 2.78% of a
    page. A name-keyed surface cannot be satisfied by a permutation.

    NO BARE CATEGORY SURFACE. One was emitted alongside these, hashing the
    list of names, and it could be replaced by a constant -- or deleted -- with
    every test green, because the isolation rule deliberately excluded it. Its
    stated job, catching a name added or removed, is what the surface LIST
    already does: adding or removing a name adds or removes an entry here, and
    the comparison is of the whole list. An uncontrolled surface justified by a
    property something else provides is `/Annots` shipping inert, which is this
    file's signature defect.
    """
    out = []
    for cat in _RESOURCE_CATEGORIES:
        num = (_materialize(doc, doc.xref_get_key(res, cat))
               if res is not None else None)
        names = []
        if num is not None:
            try:
                names = sorted(doc.xref_get_keys(num))
            except Exception:
                names = []
        for name in names:
            out.append((f"{prefix}{cat}#{name}", hashlib.sha256(
                _resolve(doc, doc.xref_get_key(num, name)).encode()).hexdigest()))
    return out


_NAME_LIST_SURFACES = {c.lower() for c in _RESOURCE_CATEGORIES}


def _head(key: str) -> str:
    """The surface a key belongs to: `xobject:M0#BBox` -> `xobject`."""
    return key.split(":")[0].split("#")[0].lower()


def xref_name(doc, pg, xref) -> str:
    """The NAME a page's resources give an xref, falling back to the number.

    Keying image surfaces on the xref number alone would make them allocation
    artifacts -- the thing `sorted(get_xobjects(), key=...)` exists to avoid --
    so prefer the resource name, which is what the content stream says.
    """
    try:
        for item in pg.get_images(full=True):
            if item[0] == xref and len(item) > 7 and item[7]:
                return str(item[7])
    except Exception:
        pass
    return str(xref)


# How far a reference chain is followed. Not unbounded: a malformed or
# adversarial file can chain arbitrarily, and every level multiplies work.
#
# The number is HEADROOM, and the first version of this comment justified it
# with "chains of 11 exist in the test suite and 7 was the measured blind
# spot" -- a hand-written figure that measurement does not support, which is
# the defect this file exists to retract, in the sentence explaining a
# constant. Measured instead, and pinned by
# `test_the_traversal_depth_committed_figures_need_is_measured`: the eight
# committed figures reach depth 0 exactly, because no resource on any of them
# is nested at all, and the crafted probe reaches 6 exactly -- both pinned by
# equality below, so this comment cannot go stale silently. Measured by
# bisecting the constant, the suite needs a bound above 5: 4 and 5 fail
# `deep-only`, 6 passes. Everything above 6 is deliberate slack for files this
# project does not yet produce, and it is slack rather than evidence.
_MAX_DEPTH = 24
_MAX_DEPTH_SEEN = 0


def _resolve(doc, val, depth: int = 0, stack=None) -> str:
    """Flatten a resource entry, SUBSTITUTING contents where the reference is.

    `xref_get_key` returns `('xref', '4 0 R')` for an indirect dict, which is
    an object NUMBER -- stable across an unchanged regeneration but useless as
    a comparison, since it says nothing about the contents. Following it is
    what makes an `/ExtGState` alpha change visible.

    IN PLACE, not appended. An earlier version stripped every `N 0 R` and then
    joined the resolved contents positionally, which erased WHICH NAME MAPPED
    TO WHAT: repointing `/A2` from an opaque graphics state to a faint one
    already referenced by `/A3` gave a byte-identical string and moved 2.78% of
    a page. That is "compared as names and object numbers" restated, which is
    the exact property this function exists to avoid.

    Cycles are broken along the CURRENT CHAIN only. A global visited set also
    suppressed legitimate second visits: an object first reached down a long
    path had its children skipped forever, so a change was invisible when a
    deep path happened to be walked before a shallow one -- measured both ways.
    """
    kind, raw = (val if isinstance(val, tuple) else ("string", str(val)))
    stack = () if stack is None else stack
    if kind in ("array", "dict"):
        # A DIRECT dict too. Only arrays and indirect references were followed,
        # so a `("dict", ...)` value fell through to `_denumber`, which erases
        # the very references inside it: two `/OCProperties` differing only in
        # WHICH group is switched off both flattened to the same text, hiding
        # 48.70% of a page. That is weaker than comparing names and numbers,
        # which is the thing this function's docstring says it does not do.
        return _substitute(doc, str(raw), depth, stack)
    if kind == "xref":
        try:
            num = int(str(raw).split()[0])
        except Exception:
            return str(raw)
        if num in stack:
            # Along the CURRENT chain only. Deleting this line changes no
            # test, and that is not an oversight to fix by inventing one:
            # `_MAX_DEPTH` already terminates a cycle, so the guard's only
            # effect is bounding the WORK a cyclic file costs, which no
            # verdict can distinguish. Recorded rather than left implying
            # coverage that does not exist.
            return "<cycle>"
        try:
            text = doc.xref_object(num, compressed=True)
        except Exception:
            return str(raw)
        body = _substitute(doc, text, depth, stack + (num,))
        # THE STREAM, when the object has one. A `/Pattern` IS a stream: a
        # tile's geometry lives there and nowhere in the dict, so changing a
        # `hatch=` rewrote the stream while every hashed surface, this one
        # included, stayed identical.
        try:
            body += "|" + hashlib.sha256(doc.xref_stream(num)).hexdigest()
        except Exception:
            pass
        return body
    return _denumber(str(raw))


def _substitute(doc, text: str, depth: int, stack) -> str:
    """Replace each `N 0 R` in `text` with what it resolves to."""
    global _MAX_DEPTH_SEEN
    _MAX_DEPTH_SEEN = max(_MAX_DEPTH_SEEN, depth)
    if depth >= _MAX_DEPTH:
        return _denumber(text) + "|<depth-limit>"

    def one(m):
        return "{" + _resolve(doc, ("xref", m.group(0)), depth + 1, stack) + "}"

    return _denumber(re.sub(r"\d+ 0 R", one, str(text)))


def _denumber(text: str) -> str:
    """Strip object NUMBERS from resolved text, keeping the reference.

    Following nested references means their numbers appear in the hash, and a
    number is an allocation artifact: a regeneration that permutes numbering
    with identical content would report a false stale, which is the same
    complaint `sorted(get_xobjects(), key=...)` exists for one surface up.
    """
    # Word-boundaried, because `(page 12 0 Rev A)` had its middle eaten. The
    # residue is deliberate and stated: a PDF STRING containing the literal
    # text `/Length 5` or `12 0 R` is still rewritten. No drawing-relevant
    # string on a census page contains either, and the alternative is a full
    # PDF tokenizer for a hash input.
    text = re.sub(r"(?<![\w.])\d+ 0 R\b", "R", str(text))
    # AND the compression artifacts. `/Length` and the filter parameters
    # describe the ENCODING, not the drawing: the stream is hashed
    # decompressed a few lines up precisely so a different zlib build compares
    # equal. Leaving them in put the artifact back by another door, and it
    # was not theoretical -- with them in, a control that rewrote a pattern's
    # tile passed even with stream hashing removed, because the two tiles were
    # different LENGTHS. That is a control passing on a byte count.
    text = re.sub(r"/Length\s+(?:\d+|R)", "", text)
    text = re.sub(r"/(?:Filter|DecodeParms)\s*(?:/\w+|<<[^<>]*>>|\[[^\]]*\])", "", text)
    return text


_COMPARISONS = 0
_COMPARED = []


def _assert_matches_committed(produced: Path) -> None:
    """The flagship's comparison, extracted so a control can INVOKE it.

    Inlined, it was unguarded: deleting its assertion, or pointing `committed`
    at the produced file, left the suite green with a genuinely stale figure
    on disk. The control below calls THIS.

    Counts its calls, because a structural guard cannot express REACHABILITY
    and kept losing that race: the AST check was tightened to require this
    call as a direct child of the flagship's loop, and `if False:` AROUND the
    loop, `try/except AssertionError` around it, a decoy loop beside a sliced
    real one, and a bare `return` at the top of the test all left the whole
    suite green with a stale figure on disk. A count cannot be satisfied by a
    loop that is lexically present and never entered.
    """
    global _COMPARISONS
    _COMPARISONS += 1
    # WHICH file, not just how many. A count says comparisons ran; it cannot
    # say what they compared, and `_assert_matches_committed(produced[0])` --
    # one token -- left seven of eight figures untouched with the entire
    # repository suite green. The AST pin requires the call to be a direct
    # child of the loop and never looks at its ARGUMENT.
    _COMPARED.append(produced.resolve())
    # WHERE the file came from, checked at runtime. Pinning the source text of
    # `produced = sorted(scratch.glob(...))` left `scratch` itself unpinned, so
    # inserting `scratch = FIG_DIR` on the line above -- one line, in one test
    # -- made all eight comparisons compare a committed figure to ITSELF. That
    # satisfied every guard at once, the call counter included, because eight
    # real comparisons genuinely ran. No pin over the source can express this;
    # the property is about the path at the moment of the call.
    assert produced.resolve().parent != FIG_DIR.resolve(), (
        f"the flagship is comparing {produced.name} against itself: it was "
        "handed the committed figure, not the regenerated one")
    committed = FIG_DIR / produced.name
    assert committed.exists(), f"{produced.name} is not committed"
    a, b = _drawing(produced), _drawing(committed)
    assert a is not None, "no PDF reader available"
    assert a == b, (
        f"article/figures/{produced.name} is not what "
        "scripts/generate_census_figures.py draws. Re-run it.")


def _assert_compared_the_census_set(paths) -> None:
    """Which files were compared, and from where. Extracted so it has a control.

    Inlined, both halves were individually deletable with the suite green --
    the defect this file names as its own signature, applied to the fix the
    previous commit shipped for it.
    """
    assert {p.stem for p in paths} == set(_census_figures()), (
        f"the flagship compared {sorted(p.stem for p in paths)}, not the "
        f"census set {_census_figures()}")
    parents = {p.parent for p in paths}
    assert len(parents) == 1 and FIG_DIR.resolve() not in parents, (
        f"the compared files came from {sorted(map(str, parents))}; they must "
        "all come from one directory that is not article/figures")


def _run_generator(scratch: Path, env: dict):
    """Run the census generator and snapshot what it wrote, in one step.

    ONE STEP because the snapshot's whole job is to notice the scratch
    directory being repopulated afterwards, and any line between the two is
    outside the window it covers -- a single line copying the committed
    figures in, placed above the snapshot, passed with all eight genuinely
    stale. There is no gap IN THE CALLER to insert one into -- the two lines
    below are still a gap, and a copy loop between them still works. That is
    the boundary below, not a claim to have closed it.

    WHAT THIS DOES NOT CLAIM: an edit to THIS function, or to any of the
    guards it feeds, defeats it. Every attack in this file's history is an
    edit to the test suite, and no test can be robust against arbitrary edits
    to itself. What these guards buy is that the shortcut has to look like
    what it is in a diff -- `scratch = FIG_DIR`, a copy loop, a sliced
    range -- rather than passing silently. That is the boundary, stated
    rather than implied by an ever-tightening series of assertions.
    """
    res = subprocess.run([sys.executable, str(GEN)], cwd=REPO,
                         capture_output=True, text=True, env=env)
    assert res.returncode == 0, (
        f"the census figure generator failed:\n{res.stderr[-800:]}")
    return res, _snapshot(scratch)


def _snapshot(d: Path) -> dict:
    """Content hashes of a directory's PDFs, to detect it being repopulated."""
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(d).glob("*.pdf"))}


def _has_date(data) -> bool:
    """Does this PDF carry a creation or modification date?

    NOT `b"/CreationDate" in data`. That was the check here, at three sites,
    and it is blind to exactly the files that most needed it: cairo -- which
    is what `dot` renders through -- writes the Info dictionary into a
    COMPRESSED OBJECT STREAM, so the literal string never appears in the file
    while the date is plainly there to a parser that reads the Info
    dictionary.

    THAT IS THE ONE OF TWO PLACES A PDF DATE LIVES. An XMP packet can carry
    `xmp:CreateDate` with no Info entry, and this returns False for it, as the
    byte scan also did. Checked: none of the 36 committed figure PDFs carries
    any XMP packet, so the blind spot is not live -- but it is a blind spot,
    not the "any parser" the sentence used to claim. An encrypted PDF reports
    no metadata and also returns False silently; a non-PDF raises loudly. Measured on
    `fig19_immune_coupling_flow.pdf` at b8299eb8: the byte scan reports clean,
    `pymupdf` reports `D:20260422163455-07'00`. The backlog count this file
    states was 7 by the byte scan and 9 in fact, and the two it missed were
    the two whose date the savefig wrapper could never have removed.

    `test_the_byte_scan_this_replaced_is_blind_to_a_real_date` pins that case
    against the committed blob, so this is a measured difference and not a
    precaution.
    """
    import pymupdf

    doc = pymupdf.open(stream=data, filetype="pdf") if isinstance(data, bytes) \
        else pymupdf.open(data)
    try:
        md = doc.metadata or {}
        return bool((md.get("creationDate") or "").strip()
                    or (md.get("modDate") or "").strip())
    finally:
        doc.close()


def _census_figures():
    """The figures the census generator writes, read from its source."""
    return sorted(_census_outputs()["pdf"])


def _census_outputs(src: str = None) -> dict:
    """Census figure stems, PER EXTENSION.

    Folding both into one name set meant a deleted `.png` savefig left the
    count unchanged, the existence check only ever looked for `{stem}.pdf`,
    and the flagship globbed `*.pdf` -- so `fig2c_census_volume.png` could
    stop being generated with the whole suite green. That is the same floor
    defect this file retracts, still open on the PNG side.
    """
    src = GEN.read_text() if src is None else src
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
        assert not _has_date(a)


def test_no_census_figure_carries_a_creation_date():
    """The class this fix actually covers."""
    stale = [f"{f}.pdf" for f in _census_figures()
             if _has_date(FIG_DIR / f"{f}.pdf")]
    assert not stale, (
        f"{len(stale)} census figures still embed a creation date, so they "
        f"cannot be checked for freshness: {stale}")


def _corpus_derived_stems():
    """The figures the sibling gate covers, from FIGURES.yaml -- not a literal.

    This was `gated = 5` next to a derived total, in the test whose own comment
    argues a count pinned by its spelling cannot notice being wrong.
    """
    import yaml

    figs = yaml.safe_load((REPO / "FIGURES.yaml").read_text())["figures"]
    return sorted(
        f["filename"] for f in figs
        if f.get("type") == "corpus-derived"
        and str(f.get("generator", {}).get("script", "")).endswith(
            "generate_figures.py"))


def _spell(n: int) -> str:
    """Small numbers as the docstring writes them.

    The map used to be the only accepted spelling, so every figure added to
    the repository meant extending it and rewriting five sentences -- which is
    maintenance on a formatting convention, not on a claim. `_states` accepts
    either form now; this stays because the prose reads better in words where
    a word exists.
    """
    words = {0: "no", 1: "one", 2: "two", 3: "three", 5: "five", 7: "seven",
             8: "eight", 9: "nine", 10: "ten", 11: "eleven", 15: "fifteen",
             17: "seventeen", 18: "eighteen", 19: "nineteen", 22: "twenty-two",
             23: "twenty-three", 24: "twenty-four", 25: "twenty-five",
             26: "twenty-six", 27: "twenty-seven", 29: "twenty-nine",
             30: "thirty", 31: "thirty-one", 32: "thirty-two"}
    return words.get(n, str(n))


def _states(doc: str, n: int, template: str) -> bool:
    """Does the docstring state `n` in `template`, as a word or as a digit?

    The COUNT is the claim; whether it is spelled is not. Requiring one
    spelling made a correct docstring fail and sent the author to a word map
    instead of to the number.
    """
    forms = {str(n), _spell(n), _spell(n).capitalize()}
    return any(template.format(n=f) in doc for f in forms)


def test_the_corpus_figure_backlog_is_stated_and_shrinking():
    """The other generator's figures are NOT covered, and saying how many is
    what stops "figures are gated now" being read as covering all of them."""
    census = {f"{f}.pdf" for f in _census_figures()}
    # FROM GIT, not the working tree. Both counts here are statements about
    # what is COMMITTED, and running the generator -- which the README tells
    # you to do -- rewrites three simulation figures in place and drops the
    # stale count from 15 to 12 without committing anything.
    import subprocess

    listed = subprocess.run(
        ["git", "ls-files", "article/figures/*.pdf"],
        cwd=REPO, capture_output=True, text=True)
    assert listed.returncode == 0, listed.stderr
    paths = listed.stdout.split()
    assert paths, "git lists no committed figures"
    names = {Path(x).name for x in paths}
    stale = []
    for rel in sorted(paths):
        if Path(rel).name in census:
            continue
        blob = subprocess.run(["git", "show", f"HEAD:{rel}"],
                              cwd=REPO, capture_output=True)
        assert blob.returncode == 0, rel
        if _has_date(blob.stdout):
            stale.append(Path(rel).name)
    doc = _module_docstring()
    # The scope limits that remain, checked against the docstring. An earlier
    # version asserted "generate_figures.py cannot run in CI at all", which is
    # FALSE -- it exits 0 without the simulation outputs -- so a guard was
    # requiring a false sentence to stay present. That claim is gone.
    assert "No PNG content, at all" in doc
    # DERIVED, not a literal. The old assertion pinned the phrase "The twenty
    # non-census PDFs", so the sentence stayed green while this PR gated five
    # of them and removed the creation date from those five. A count stated in
    # prose and pinned by its own spelling cannot notice being wrong.
    assert "Most non-census PDFs" in doc
    # TRACKED, not on disk. `generate_figures.py` writes seven stems that are
    # never committed, so after the regeneration step the README documents,
    # this counted 29 and reported "29 non-census figures are committed" --
    # a false message from a guard about counting honestly.
    total = len([n for n in names if n not in census])
    gated = len(_corpus_derived_stems())
    # DERIVED from the count rather than pinned as a literal. The old form
    # hard-coded "Twenty-two" on BOTH sides of the `and`, so adding a figure
    # produced a failure naming the right number in the message and the wrong
    # one in the assertion, and the fix looked like editing two constants
    # instead of one.
    assert _states(doc, total, "{n} committed figures"), (
        f"{total} non-census figures are committed; the docstring does not "
        f"say so")
    assert _states(doc, total - gated, "the other {n} are not"), (
        f"{total - gated} non-census figures are ungated and the docstring "
        "says otherwise")
    import yaml as _yaml

    spec = _yaml.safe_load((REPO / "FIGURES.yaml").read_text())["figures"]
    sim = {f["filename"] for f in spec if f.get("type") == "simulation"}
    ungated = {Path(n).stem for n in names if n not in census} - set(
        _corpus_derived_stems())
    n_sim = len(ungated & sim)
    # WHICH of the ten are blocked by gitignored inputs, derived from
    # FIGURES.yaml rather than asserted. The first version of this bullet said
    # none of the ten is checkable, which is false: two do not read
    # `simulations/output/` at all.
    import subprocess as _sp

    tracked_files = set(_sp.run(["git", "ls-files"], cwd=REPO,
                                capture_output=True, text=True).stdout.split())
    blocked, free = [], []
    for f in spec:
        if f["filename"] not in (ungated & sim):
            continue
        ins = f.get("inputs") or []
        (blocked if any(not (i in tracked_files) for i in ins) else free).append(
            f["filename"])
    # SPELLED INTO THE PROSE, not compared against literals here. The first
    # version asserted `len(blocked) == 8 and len(free) == 2` -- two constants
    # in the test compared against FIGURES.yaml -- so the docstring could say
    # three and seven, or that none of the ten is checkable, and stay green.
    # That is the defect this whole paragraph exists to retract, in its guard.
    assert any(_states(doc, len(blocked), f"{{n}} of the {form}")
               for form in {str(len(blocked) + len(free)),
                            _spell(len(blocked) + len(free))}), (
        f"{len(blocked)} of the simulation figures read an untracked input and "
        f"the docstring says otherwise. blocked={sorted(blocked)}")
    assert _states(doc, len(free), "The other {n} are not blocked"), (
        f"{len(free)} of them do not, and the docstring says otherwise. "
        f"free={sorted(free)}")
    for name in free:
        assert name.split("_")[0] in doc, (
            f"{name} is not blocked by a gitignored input and the docstring "
            "does not name it among those that are not")
    # THE DENOMINATORS DIFFER, and the guard should say so rather than let a
    # reader assume they match. `stale` below is counted over all 22
    # non-census figures; this sentence is about the seventeen. They agree
    # today (both dated figures are inside the seventeen), and a dated
    # corpus-derived figure would make the guard demand a sentence that is
    # wrong about its own scope.
    # The denominator is DERIVED too. It was pinned as the literal
    # "seventeen" beside a derived numerator, so adding one ungated figure
    # made a true sentence fail for the wrong reason.
    assert any(_states(doc, n_sim, f"{{n}} of those {form}")
               for form in {str(total - gated), _spell(total - gated)}), (
        f"FIGURES.yaml marks {n_sim} of the ungated figures as simulation and "
        "the docstring says otherwise -- the first version said eight, which "
        "was the number this branch regenerated rather than the number that "
        "exist, and nothing measured it")
    assert f"{_spell(len(stale))} still embed a creation date" in doc, (
        f"{len(stale)} still embed a creation date and the docstring says "
        "otherwise")
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
        # `generate_figures.py` carries a note about it. That note used to say
        # the Rust binary writes the figure; it does not -- the binary prints
        # the DATA as JSON and the PDF's producer is Matplotlib -- and the note
        # says so now.
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
        res, before = _run_generator(scratch, env)
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
        # The scratch directory must not have been REPOPULATED between the
        # generator writing it and the comparisons reading it. Copying the
        # eight committed figures over the regenerated ones passes the produced
        # set equality, the path check and the call count, because every file
        # is still in the scratch directory and every comparison still runs --
        # it just compares each figure to itself by another route.
        assert _snapshot(scratch) == before, (
            "the regenerated figures were overwritten between the generator "
            "run and the comparison")


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
        return [_head(k) for k, _ in d]

    def primary(k):
        """The CONTENT-BEARING hash for this surface, not a weaker sibling.

        `xobject:M0` is the form's stream and `xobject:M0#BBox` is one key of
        its dict; collapsing both to `xobject` let the STREAM hash be replaced
        by a constant with this control green, because the `/BBox` movement
        alone satisfied the isolation assertion. So for those surfaces the
        content-bearing key is the one WITHOUT a `#`.

        For a resource category it is the other way round, and a depth rule
        got that backwards: `ExtGState` hashes the list of NAMES and
        `ExtGState#A2` hashes what the name resolves to. The name list cannot
        see an alpha change, so requiring it would accept the weaker hash.
        """
        if _head(k) != surface:
            return False
        return ("#" in k) if surface in _NAME_LIST_SURFACES else ("#" not in k)

    assert surface in keys(da) and surface in keys(db), (
        f"the {label} control does not emit a {surface!r} surface in both "
        f"arms, so it tests PRESENCE rather than content. Arms carry "
        f"{keys(da)} and {keys(db)}")
    assert len(da) == len(db), (
        f"the {label} control changes the surface LIST, not its contents: "
        f"{keys(da)} vs {keys(db)}")
    moved = {k for (k, x), (_, y) in zip(da, db) if x != y}
    assert {_head(k) for k in moved} == {surface}, (
        f"the {label} control moves {sorted(moved) or 'nothing'}; a control "
        f"that does not isolate {surface!r} proves nothing about it")
    assert any(primary(k) for k in moved), (
        f"the {label} control moves only {sorted(moved)}, which are named keys "
        f"hanging off {surface!r} rather than {surface!r} itself -- so that "
        "hash could be a constant with this control still green")


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
        # And the array form is followed rather than stringified. Contains
        # rather than equals: references are now substituted IN PLACE, so the
        # array's resolution carries the element's contents inside the
        # surrounding brackets rather than being identical to it. Equality
        # would be testing the punctuation.
        arr = ("array", f"[{ka[1]}]")
        assert _resolve(da, ka) in _resolve(da, arr), (
            "_resolve does not follow a one-element array to the content the "
            "bare reference resolves to")
        assert str(ka[1]) not in _resolve(da, arr), (
            "_resolve left an object number in the array's resolution, which "
            "is the allocation artifact it exists to strip")
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
    # STATEMENT LEVEL OF THE FUNCTION, not anywhere in its subtree: `walk`
    # yields the loop whatever encloses it, so `if False:` and a `def` that is
    # never called both satisfied a predicate written over the whole subtree.
    body = list(fn.body)
    for st in fn.body:
        if isinstance(st, (_ast.With, _ast.AsyncWith)):
            body.extend(st.body)
    looped = [
        node for node in body
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
    # And nothing may return before it. `if True: return` above the loop left
    # every test green with a stale figure, because the loop was still there
    # to be found.
    for st in body:
        if st is looped[0]:
            break
        assert not isinstance(st, _ast.Return), (
            "the flagship returns before its comparison loop")
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
        # `/Decode [1 0]` INVERTS the alpha plane without rewriting the
        # stream. The previous perturbation wrote `b"\xff" * length` through
        # `update_stream`, and the plane carries `/DecodeParms
        # <</Predictor 10>>`: that round trip is lossy (4,514,664 bytes in,
        # 4,512,816 out), so writing the stream's OWN bytes also moved the
        # surface and also rendered differently -- the control raised either
        # way and demonstrated a re-encoding artifact rather than a blanked
        # figure.
        before_dec = str(doc.xref_get_key(smasks[0], "Decode")[1])
        doc.xref_set_key(smasks[0], "Decode", "[1 0]")
        assert str(doc.xref_get_key(smasks[0], "Decode")[1]) != before_dec, (
            "the /Decode edit did not land, so nothing was perturbed")
        # Saved to a NEW path: pymupdf refuses a non-incremental save over the
        # file it opened.
        out = tmp_path / "perturbed.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    # THE PERTURBATION MUST RENDER DIFFERENTLY, checked without `_drawing`.
    # Replacing the blanking with a write of the stream's OWN bytes still made
    # this test raise: the alpha plane carries `/DecodeParms <</Predictor 10>>`
    # and pymupdf's `update_stream` round trip is lossy through it (4,514,664
    # bytes in, 4,512,816 out), so the surface moved whatever was written. The
    # control was passing on a re-encoding artifact, and its docstring claims a
    # perturbed COPY on the surface this branch was blocked over. Pixels are
    # the independent check: they come from the renderer, not from the hash
    # under test.
    def _pixels(path):
        d = pymupdf.open(path)
        try:
            return d[0].get_pixmap(dpi=72).samples
        finally:
            d.close()

    before, after = _pixels(src), _pixels(out)
    assert len(before) == len(after) and before != after, (
        "the perturbed copy renders identically to the committed figure, so "
        "this control demonstrates nothing about the comparison it invokes")

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
        assert zs != os_, (
            "the two glyph streams are identical, so swapping them is a no-op "
            "and this control documents a limit it never exercised")
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    # AND THE SWAP MUST HAVE SURVIVED THE SAVE. Its sibling below re-reads the
    # saved file for exactly this reason -- an edit that silently failed to
    # apply is indistinguishable from a surface the gate covers, and reads as
    # reassurance either way -- and this control never got that treatment:
    # replacing the swap with two writes of the ORIGINAL streams left it green.
    chk = pymupdf.open(out)
    try:
        assert chk.xref_stream(z) == os_ and chk.xref_stream(o) == zs, (
            "the glyph swap did not survive the save, so 'not detected' is "
            "reporting an unlanded mutation as a documented limit")
    finally:
        chk.close()

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
    assert len(a) == len(b), (
        "the arms carry different surface LISTS, so the isolation set below\n         is computed over a zip() that truncates silently")
    moved = {k for (k, v), (k2, v2) in zip(a, b) if v != v2}
    assert a != b, (
        f"a form xobject's /{key} can be rewritten without the gate noticing; "
        f"on {fig} that moves the drawing and ships green")
    # ISOLATION: the named surface moved and nothing else did. A control that
    # moves two surfaces proves neither, which is how the `alpha=1.0` smask
    # control passed while testing presence rather than content.
    assert moved == {f"xobject:{forms[0][1]}#{key}"}, (
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
    assert len(a) == len(b), (
        "the arms carry different surface LISTS, so the isolation set below\n         is computed over a zip() that truncates silently")
    moved = {k for (k, v), (k2, v2) in zip(a, b) if v != v2}
    assert moved == {"Group"}, (
        f"expected only the /Group surface to move, got {sorted(moved)}")


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
    assert len(a) == len(b), (
        "the arms carry different surface LISTS, so the isolation set below\n         is computed over a zip() that truncates silently")
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
    assert len(a) == len(b), (
        "the arms carry different surface LISTS, so the isolation set below\n         is computed over a zip() that truncates silently")
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved == {"smask#Decode"}, (
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
    assert len(a) == len(b), (
        "the arms carry different surface LISTS, so the isolation set below\n         is computed over a zip() that truncates silently")
    moved = {k.split(":")[0] for (k, x), (_, y) in zip(a, b) if x != y}
    assert "Rotate" in moved, (
        f"an inherited /Rotate moved {sorted(moved) or 'nothing'}; the page "
        "renders rotated, so the geometry surface must see it")


def _annotated(doc, tmp_path, name):
    """A copy of fig2c carrying one annotation, with `/Annots` INDIRECT.

    pymupdf writes `/Annots` back as a DIRECT array, which the array branch
    already handled -- so the first version of this control passed with the
    recursion removed, testing the path that was never broken. Every census
    page has it indirect (`10 0 R -> []`), and that is the shape to restore.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    copy = tmp_path / f"{name}_c.pdf"
    copy.write_bytes((FIG_DIR / "fig2c_census_volume.pdf").read_bytes())
    doc_ = pymupdf.open(copy)
    try:
        doc_[0].add_rect_annot((50, 50, 300, 300))
        out = tmp_path / f"{name}.pdf"
        doc_.save(out, deflate=True)
    finally:
        doc_.close()
    doc_ = pymupdf.open(out)
    try:
        annot = doc_[0].first_annot
        assert annot is not None, "the annotation was not written"
        arr = doc_.get_new_xref()
        doc_.update_object(arr, f"[{annot.xref} 0 R]")
        doc_.xref_set_key(doc_[0].xref, "Annots", f"{arr} 0 R")
        assert doc_.xref_get_key(doc_[0].xref, "Annots")[0] == "xref"
        doc_.save(out, incremental=True, encryption=0)
    finally:
        doc_.close()
    return out


@pytest.mark.parametrize("what", ["rect", "appearance"])
def test_an_annotation_mutated_in_place_is_detected(tmp_path, what):
    """`_resolve` following an INDIRECT reference to something that is not an
    array -- which is where every census page's `/Annots` actually leads.

    `/Annots` is `10 0 R -> []` on all eight pages, so the array branch -- the
    one that exists because a bare object number says nothing about contents --
    was never reached on a real figure. ADDING an annotation was caught either
    way, because the array text changes; moving one already there was not.

    `appearance` is the sharper case and was missed by the first version of
    this control: `/AP /N` is the stream a renderer actually DRAWS for an
    annotation, and the dict text stops one reference short of it. Rewriting
    it moved 17.29% of the page with every hashed surface identical, so
    "annotations are checked" covered the box and not the drawing.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    base = _annotated(None, tmp_path, "a")
    doc = pymupdf.open(base)
    try:
        annot = doc[0].first_annot
        if what == "rect":
            before = str(doc.xref_get_key(annot.xref, "Rect")[1])
            doc.xref_set_key(annot.xref, "Rect", "[200 200 400 400]")
            assert str(doc.xref_get_key(annot.xref, "Rect")[1]) != before
        else:
            ap = doc.xref_get_key(annot.xref, "AP/N")
            assert ap[0] == "xref", f"the annotation has no appearance stream: {ap}"
            num = int(str(ap[1]).split()[0])
            body = doc.xref_stream(num)
            assert body, "the appearance stream is empty"
            doc.update_stream(num, body + b"\nq 0 0 1 rg 0 0 50 50 re f Q\n")
        out = tmp_path / "b.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    chk = pymupdf.open(out)
    try:
        assert chk.xref_get_key(chk[0].xref, "Annots")[0] == "xref", (
            "the save flattened /Annots to a direct array, so this control "
            "exercises the branch that already worked")
    finally:
        chk.close()

    a, b = _drawing(out), _drawing(base)
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved == {"Annots"}, (
        f"an annotation's {what} can be rewritten with every hashed surface "
        f"identical. Moved: {sorted(moved) or 'nothing'}")


def _ocg(doc, name):
    x = doc.get_new_xref()
    doc.update_object(x, f"<</Type/OCG/Name({name})>>")
    return x


def test_hiding_content_with_an_optional_content_group_is_detected(tmp_path):
    """An OCG switched `/OFF` in the catalog and named by `/OC` renders fig5c
    BLANK -- 49.21% of the page, larger than the alpha blanking this branch
    already blocked on -- with every drawing surface identical.

    Renderer-dependent (MuPDF and Acrobat honour it, Quartz does not) and
    dropped by pdfTeX, so the exposure is the COMMITTED STANDALONE FIGURE,
    which is the artifact MANIFEST.sha256 hashes.

    ONE STEP: the catalog gains a configuration it did not have. Only the
    document-level surface may move.
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
        x = _ocg(doc, "hidden")
        doc.xref_set_key(
            doc.pdf_catalog(), "OCProperties",
            f"<</OCGs[{x} 0 R]/D<</OFF[{x} 0 R]/Order[]/BaseState/ON>>>>")
        out = tmp_path / "m.pdf"
        doc.save(out, deflate=True)
    finally:
        doc.close()

    a, b = _drawing(out), _drawing(src)
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, v), (k2, v2) in zip(a, b) if v != v2}
    assert moved == {"OCProperties"}, (
        f"expected only the optional-content surface to move, got {sorted(moved)}")


def test_repointing_an_images_oc_with_the_catalog_fixed_is_detected(tmp_path):
    """THE OTHER HALF, and the half that was missing.

    The control above wrote `/OC` onto both heatmap images and asserted that
    only `OCProperties` moved -- which is a statement that the `/OC` writes
    moved NOTHING. Deleting them left it green while the rendered difference
    fell from 49.21% to 0.00%: an OCG nothing references hides nothing. The
    verdict came entirely from the catalog string.

    So: identical catalog in both arms, carrying two OCGs of which one is
    `/OFF`, and the images' `/OC` pointing at the visible one in the first arm
    and the hidden one in the second. Nothing else differs, and 49.21% of the
    page does.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    arms = []
    for i in (0, 1):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            # Allocated in the same order in both arms, so the two catalogs
            # are byte-identical and only the reference moves.
            off, on = _ocg(doc, "off"), _ocg(doc, "on")
            doc.xref_set_key(
                doc.pdf_catalog(), "OCProperties",
                f"<</OCGs[{off} 0 R {on} 0 R]/D<</OFF[{off} 0 R]"
                "/Order[]/BaseState/ON>>>>")
            imgs = doc[0].get_images(full=True)
            assert imgs, "fig5c no longer carries the raster this control hides"
            for img in imgs:
                doc.xref_set_key(img[0], "OC", f"{(on, off)[i]} 0 R")
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a, b = _drawing(arms[1]), _drawing(arms[0])
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, v), (k2, v2) in zip(a, b) if v != v2}
    assert moved and all(k.endswith("#OC") for k in moved), (
        "an image's /OC can be repointed from a visible optional-content "
        "group to a hidden one with the catalog identical; that hides the "
        f"whole raster. Moved: {sorted(moved) or 'nothing'}")


def _with_resource(doc, category, name, obj, stream=None):
    """Register a new object under `Resources/<category>/<name>`."""
    xref = doc.get_new_xref()
    doc.update_object(xref, obj)
    if stream is not None:
        doc.update_stream(xref, stream)
    # ON THE CATEGORY OBJECT ITSELF. `/Resources` and each category under it
    # are INDIRECT on a matplotlib page, and pymupdf refuses a path through an
    # indirect ("path to 'P0' has indirects"); asking it to replace the whole
    # `Resources/<category>` dict instead succeeds and writes the literal
    # string "fitz: replace me!", so the two arms would differ only in a
    # placeholder and the control would pass while testing nothing.
    res = doc.xref_get_key(doc[0].xref, "Resources")
    assert res[0] == "xref", f"/Resources is not indirect: {res}"
    rnum = int(str(res[1]).split()[0])
    cat = doc.xref_get_key(rnum, category)
    if cat[0] != "xref":
        # `/ColorSpace` and `/Properties` are absent from a matplotlib page
        # entirely, which is the same reason `/Pattern` and `/Shading` shipped
        # uncontrolled: an ABSENT category looks identical whether it is
        # followed or not, so the control has to create it.
        holder = doc.get_new_xref()
        doc.update_object(holder, "<<>>")
        doc.xref_set_key(rnum, category, f"{holder} 0 R")
        cat = doc.xref_get_key(rnum, category)
    assert cat[0] == "xref", f"/{category} is not an indirect dict: {cat}"
    num = int(str(cat[1]).split()[0])
    doc.xref_set_key(num, name, f"{xref} 0 R")
    got = doc.xref_object(num, compressed=True)
    assert name in got and "replace me" not in got, (
        f"the resource was not registered: {got[:80]}")
    return xref


@pytest.mark.parametrize("category", ["Pattern", "Shading", "ColorSpace",
                                     "Properties"])
def test_a_named_resources_contents_are_followed(tmp_path, category):
    """`/Pattern` and `/Shading` were hashed by NAME, not contents.

    `_resolve` returned `<</P0 12 0 R>>`, which says exactly as little about P0
    as `4 0 R` said about the ExtGState dict -- the complaint the function was
    written for, one level short. A pattern IS a stream: its tile geometry
    lives there and nowhere in the dict, so a `hatch=` change rewrote the
    stream while every hashed surface stayed identical.

    Crafted, because no census figure carries a non-empty `/Pattern` today --
    every one is `<<>>`. That is precisely why it shipped: an empty resource
    hashed as defence-in-depth looks identical whether it is followed or not,
    which is how `/Annots` shipped inert too.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    tile = b"0.5 w 0 0 m 8 8 l S\n"
    obj = ("<</Type/Pattern/PatternType 1/PaintType 1/TilingType 1"
           "/BBox[0 0 8 8]/XStep 8/YStep 8/Resources<<>>>>")
    arms = []
    for i, body in enumerate((tile, tile.replace(b"0.5 w", b"3 w"))):
        copy = tmp_path / f"a{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            _with_resource(doc, category, "P0", obj, body)
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a, b = _drawing(arms[1]), _drawing(arms[0])
    assert len(a) == len(b), (
        "the arms carry different surface LISTS, so the isolation set below\n         is computed over a zip() that truncates silently")
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved and all(_head(k) == category.lower() for k in moved) \
            and any("#" in k for k in moved), (
        f"a /{category} entry's contents can be rewritten without the gate "
        f"noticing -- it compares the NAME and the object number. Moved: {sorted(moved)}")


def test_the_flagship_actually_performs_its_comparisons():
    """THE MERGE BLOCKER: reachability, measured rather than parsed.

    Six edits left the whole suite green with a genuinely stale committed
    figure -- `if False:` around the loop, `try/except AssertionError` around
    it, `while False:`, a decoy loop beside a sliced real one, the loop moved
    into a nested function, and `if True: return` above it -- and a seventh,
    a bare `return` as the flagship's first statement, needed no edit to the
    loop at all. Every one satisfies any predicate written over the source,
    because the loop is still lexically there. The guard's own comment claimed
    those exact attacks were closed.

    So this RUNS the flagship and counts what it compared. A count cannot be
    satisfied by a loop that is never entered, by a swallowed AssertionError,
    or by an early return.

    Disabling the gate now takes edits to two tests rather than one. That is a
    higher bar, not a proof: a `pytest.skip()` CALL as the flagship's first
    statement silences this test too and exits 0, and this file does not
    pretend to cover that.

    Measured, and the CONDITION matters -- two reviewers measured opposite
    verdicts for `@pytest.mark.skip` because one staged a stale figure and one
    did not, and the first version of this paragraph stated the verdict without
    the condition, which made it false as written:

    - WITH A STALE FIGURE ON DISK, a marker does not hide it. `@pytest.mark.skip`
      leaves this test RED (1 failed) because it calls the plain function, and
      `@pytest.mark.slow` with `-m "not slow"` likewise. That is the case that
      matters, and it is why the marker route is not the hole.
    - ON A CLEAN TREE the marked run still exits 0, reporting `1 skipped` or
      `1 deselected`. A marker cannot be read as evidence the gate ran.

    An in-body `pytest.skip()` call silences both tests in either case. So the
    hole is that one line, and the disclosure is about that line rather than
    about markers.
    """
    global _COMPARISONS
    before, mark = _COMPARISONS, len(_COMPARED)
    test_the_committed_census_figures_are_what_the_generator_draws()
    done = _COMPARISONS - before
    assert done == len(_census_figures()), (
        f"the flagship compared {done} figures, not {len(_census_figures())}; "
        "its loop is present in the source but not reaching "
        "_assert_matches_committed at runtime")
    # AND WHICH ONES. Eight calls all naming the same file satisfies a count
    # exactly as well as eight calls naming eight files, and
    # `_assert_matches_committed(produced[0])` did precisely that -- one token,
    # seven figures never compared, the whole repository suite green.
    _assert_compared_the_census_set(_COMPARED[mark:])


def _page_resources_xref(doc):
    res = doc.xref_get_key(doc[0].xref, "Resources")
    assert res[0] == "xref", f"/Resources is not indirect: {res}"
    return int(str(res[1]).split()[0])


def test_resources_inherited_from_the_pages_node_are_followed(tmp_path):
    """`/Resources` is INHERITABLE, and was read as a raw page key.

    Hoisting it onto the `/Pages` node renders the page identically and makes
    every category below it read `null` forever: measured, an alpha change with
    the dict hoisted moved 53% of fig5c with no surface moving, while the same
    change with `/Resources` left on the page is caught. That is simultaneously
    the inherited-`/Rotate` defect this file fixed one column over and the
    `Resources/Annots` defect it names as its signature failure.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    arms = []
    for i in (0, 1):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            rnum = _page_resources_xref(doc)
            parent = doc.xref_get_key(doc[0].xref, "Parent")
            assert parent[0] == "xref", "the page has no /Pages node to inherit from"
            doc.xref_set_key(int(str(parent[1]).split()[0]), "Resources",
                             f"{rnum} 0 R")
            doc.xref_set_key(doc[0].xref, "Resources", "null")
            assert doc.xref_get_key(doc[0].xref, "Resources")[0] == "null"
            if i:
                gs = int(str(doc.xref_get_key(rnum, "ExtGState")[1]).split()[0])
                for k in doc.xref_get_keys(gs):
                    doc.xref_set_key(gs, f"{k}/ca", "0.05")
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a, b = _drawing(arms[1]), _drawing(arms[0])
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved and all(_head(k) == "extgstate" for k in moved), (
        "an alpha inside an INHERITED /Resources moved nothing; the resource "
        f"categories are being read as raw page keys. Moved: {sorted(moved)}")


def test_a_form_xobjects_own_resources_are_followed(tmp_path):
    """A form carries its own `/Resources`, and only the PAGE's were read.

    An alpha inside a form's own graphics-state dict moved every scatter
    marker on fig16c -- 0.53% of the page -- with the form's stream
    byte-identical and no surface moving. Nested form STREAMS were already
    covered, which made the gap exactly the non-XObject categories on a form.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig16c_trial_share.pdf"
    arms = []
    for i, alpha in enumerate(("1", "0.05")):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            forms = _forms(doc, doc[0])
            assert forms, "fig16c no longer draws through a form xobject"
            state = doc.get_new_xref()
            doc.update_object(state, f"<</Type/ExtGState/ca {alpha}>>")
            fres = doc.xref_get_key(forms[0][0], "Resources")
            if fres[0] == "xref":
                rnum = int(str(fres[1]).split()[0])
            else:
                rnum = doc.get_new_xref()
                doc.update_object(rnum, "<<>>")
                doc.xref_set_key(forms[0][0], "Resources", f"{rnum} 0 R")
            holder = doc.get_new_xref()
            doc.update_object(holder, f"<</Q1 {state} 0 R>>")
            doc.xref_set_key(rnum, "ExtGState", f"{holder} 0 R")
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a, b = _drawing(arms[1]), _drawing(arms[0])
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved and all(k.startswith("xobject:") and "ExtGState" in k
                         for k in moved), (
        "a graphics state inside a form's OWN resources moved nothing. "
        f"Moved: {sorted(moved)}")


def test_marked_content_optional_content_is_detected(tmp_path):
    """THE THIRD OPTIONAL-CONTENT DOOR, and the comment said there were two.

    `/OC /MC0 BDC ... EMC` inside a content stream resolves `/MC0` through
    `/Resources/Properties`, which was not read anywhere. With the content
    stream and the catalog byte-identical in both arms and only that mapping
    repointed from a visible group to an `/OFF` one, 62.26% of fig5c changes --
    larger than the 49.21% this branch already blocked on, by the same
    mechanism, on the same figure.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    arms = []
    for i in (0, 1):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            off, on = _ocg(doc, "off"), _ocg(doc, "on")
            doc.xref_set_key(
                doc.pdf_catalog(), "OCProperties",
                f"<</OCGs[{off} 0 R {on} 0 R]/D<</OFF[{off} 0 R]"
                "/Order[]/BaseState/ON>>>>")
            props = doc.get_new_xref()
            doc.update_object(props, f"<</MC0 {(on, off)[i]} 0 R>>")
            doc.xref_set_key(_page_resources_xref(doc), "Properties",
                             f"{props} 0 R")
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a, b = _drawing(arms[1]), _drawing(arms[0])
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved and all(_head(k) == "properties" for k in moved), (
        "the marked-content optional-content mapping can be repointed from a "
        f"visible group to a hidden one and nothing moves. Moved: {sorted(moved)}")


def test_a_resource_name_repointed_to_a_sibling_is_detected(tmp_path):
    """WHICH NAME MAPS TO WHAT, not the multiset of reachable contents.

    `_resolve` used to strip every `N 0 R` and join the resolved contents
    positionally. Repointing `/A2` at a graphics state ALREADY referenced by
    another name therefore produced a byte-identical string -- 2.78% of a page,
    nothing moved -- which is "compared as names and object numbers" restated,
    the exact property the function exists to avoid.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    arms = []
    for i in (0, 1):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            gs = int(str(doc.xref_get_key(
                _page_resources_xref(doc), "ExtGState")[1]).split()[0])
            names = sorted(doc.xref_get_keys(gs))
            assert len(names) >= 2, f"fig2c carries too few graphics states: {names}"
            # Both arms hold the SAME two objects; only the mapping differs.
            objs = []
            for alpha in ("1", "0.05"):
                x = doc.get_new_xref()
                doc.update_object(x, f"<</Type/ExtGState/CA 1/ca {alpha}>>")
                objs.append(x)
            doc.xref_set_key(gs, names[0], f"{objs[0]} 0 R")
            doc.xref_set_key(gs, names[1], f"{objs[i]} 0 R")
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a, b = _drawing(arms[1]), _drawing(arms[0])
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved and all("#" in k and _head(k) == "extgstate" for k in moved), (
        "a resource NAME can be repointed to an object already reachable "
        "under another name with nothing moving, so the comparison is of "
        f"contents-as-a-set rather than of the mapping. Moved: {sorted(moved)}")


def test_denumber_strips_allocation_and_encoding_artifacts():
    """`_denumber` had no control at all, and it is load-bearing.

    Replacing it with the identity function left the whole suite green, and so
    did removing its `/Length` strip -- which matters because with `/Length`
    left in, a control that rewrites a stream passes on a BYTE COUNT even when
    the stream itself is not hashed. Two things that must both hold: object
    numbers are allocation artifacts, and `/Length` and the filter parameters
    describe the encoding rather than the drawing, since streams are hashed
    decompressed precisely so a different zlib build compares equal.
    """
    assert _denumber("<</X 12 0 R/Length 5>>") == _denumber("<</X 99 0 R/Length 7>>"), (
        "_denumber does not erase object numbers and stream lengths, so an "
        "unchanged figure can report stale after a renumbering re-save")
    assert _denumber("<</Filter/FlateDecode/ca 0.18>>") == _denumber("<</ca 0.18>>"), (
        "_denumber does not strip the filter, so re-encoding an identical "
        "stream reports a false stale")
    assert _denumber("<</ca 0.18>>") != _denumber("<</ca 0.90>>"), (
        "_denumber erases the values it exists to preserve")


def test_comparing_a_committed_figure_against_itself_is_refused(tmp_path):
    """CONTROL for the runtime path check, which had none.

    Removing that assertion left the whole file green, because nothing invoked
    `_assert_matches_committed` with a path inside `article/figures` -- and the
    attack it exists to stop is exactly that: one line, `scratch = FIG_DIR`,
    which every source-level pin still accepts because eight real comparisons
    do run. A guard with no control is how three surfaces shipped inert here.
    """
    committed = FIG_DIR / f"{_census_figures()[0]}.pdf"
    assert committed.exists()
    with pytest.raises(AssertionError, match="against itself"):
        _assert_matches_committed(committed)
    # And it must still ACCEPT an honest scratch copy, or the check would be
    # satisfied by refusing everything.
    copy = tmp_path / committed.name
    copy.write_bytes(committed.read_bytes())
    _assert_matches_committed(copy)


def _chain(doc, alpha, links):
    """`links` indirect objects ending at a graphics state with this alpha."""
    target = doc.get_new_xref()
    doc.update_object(target, f"<</Type/ExtGState/CA 1/ca {alpha}>>")
    link = target
    for _ in range(links):
        nxt = doc.get_new_xref()
        doc.update_object(nxt, f"<</Next {link} 0 R>>")
        link = nxt
    return link, target


@pytest.mark.parametrize("case", ["deep-only", "deep-then-shallow"])
def test_resolution_does_not_depend_on_the_path_that_reaches_an_object(
        tmp_path, case):
    """Two properties of the traversal, each with a way to fail.

    `deep-only` puts the change at the end of a long chain, so a traversal
    limit set too low stops seeing it -- the limit is real, and until this
    existed, lowering it from 24 to 3 changed nothing anyone could observe.

    `deep-then-shallow` reaches the SAME object twice inside ONE resolution,
    first down a chain long enough to exhaust the limit and then directly, and
    requires the same verdict. It pins traversal-order independence.

    WHAT THIS DOES NOT PROVE, stated because the obvious reading is wrong:
    cycles are broken per CHAIN rather than with a global visited set, and I
    could NOT build a case where swapping in a global set changes the verdict.
    Making the surfaces per-NAME removed the reachable version of that bug --
    each name is now its own resolution, so nothing carries between them. The
    per-chain rule is therefore a design choice this control does not exercise,
    and calling it "verified" would be the unmeasured sentence this file
    exists to catch. What IS exercised: lowering the limit from 24 to 3 fails
    `deep-only`, which was invisible before this test existed.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"

    def build(tag, alpha):
        copy = tmp_path / f"{tag}_c.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            links = 6 if case == "deep-only" else _MAX_DEPTH + 2
            head, target = _chain(doc, alpha, links)
            holder = doc.get_new_xref()
            if case == "deep-then-shallow":
                # BOTH references under ONE name, so they are resolved in one
                # traversal. Per-name surfaces mean two names are two separate
                # resolutions, which would hide the property this pins.
                inner = doc.get_new_xref()
                doc.update_object(inner, f"<</X {head} 0 R/Y {target} 0 R>>")
                doc.update_object(holder, f"<</A {inner} 0 R>>")
            else:
                doc.update_object(holder, f"<</A {head} 0 R>>")
            doc.xref_set_key(_page_resources_xref(doc), "Properties",
                             f"{holder} 0 R")
            out = tmp_path / f"{tag}.pdf"
            doc.save(out, deflate=True)
            return out
        finally:
            doc.close()

    a, b = _drawing(build(f"{case}_a", "1")), _drawing(build(f"{case}_b", "0.05"))
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved, (
        f"[{case}] a change is invisible because of the PATH that reaches it, "
        "not because of what it is")



def _ocg_pair(doc):
    return _ocg(doc, "a"), _ocg(doc, "b")


@pytest.mark.parametrize("shape", ["ocproperties", "ocmd", "extgstate"])
def test_a_direct_dict_value_is_followed(tmp_path, shape):
    """A value returned as a DIRECT dict, which was not followed at all.

    `_resolve` special-cased arrays and indirect references; a `("dict", ...)`
    value fell through to `_denumber`, which erases the references INSIDE it.
    Two `/OCProperties` differing only in which group is switched off both
    flattened to the same text -- so the comparison was weaker than comparing
    names and object numbers, which is the property the docstring claims.

    All three shapes are ordinary PDF that renders identically, and each moved
    48.70-62.42% of fig5c with no surface moving.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    arms = []
    for i in (0, 1):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            a, b = _ocg_pair(doc)
            if shape == "extgstate":
                state = doc.get_new_xref()
                doc.update_object(state, "<</Type/ExtGState/CA 1/ca %s>>"
                                  % ("1", "0.05")[i])
                doc.xref_set_key(_page_resources_xref(doc), "ExtGState",
                                 f"<</A1 {state} 0 R>>")
            else:
                off = b if shape == "ocmd" else (a, b)[i]
                doc.xref_set_key(
                    doc.pdf_catalog(), "OCProperties",
                    f"<</OCGs[{a} 0 R {b} 0 R]/D<</OFF[{off} 0 R]"
                    "/Order[]/BaseState/ON>>>>")
                for img in doc[0].get_images(full=True):
                    doc.xref_set_key(img[0], "OC",
                                     f"<</Type/OCMD/OCGs {(a, b)[i]} 0 R>>"
                                     if shape == "ocmd" else f"{b} 0 R")
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a_, b_ = _drawing(arms[1]), _drawing(arms[0])
    assert len(a_) == len(b_), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a_, b_) if x != y}
    assert moved, (
        f"[{shape}] a direct dict value can be repointed at different contents "
        "with nothing moving; direct dicts are not being followed")


@pytest.mark.parametrize("where", ["direct", "deep-parent"])
def test_resources_are_found_however_they_are_attached(tmp_path, where):
    """Two more doors to "every category reads null forever".

    `direct` inlines `/Resources` on the page -- legal, identical rendering,
    and what several optimisers emit -- which the first version of
    `_effective_resources` treated as ABSENT because it looked for an object
    number. `deep-parent` hoists it 34 nodes up a page tree, past the
    unstated 32-link bound, which was equally invisible: lowering that bound
    to 2 changed nothing anyone could observe.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    arms = []
    for i in (0, 1):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            rnum = _page_resources_xref(doc)
            state = doc.get_new_xref()
            doc.update_object(state,
                              "<</Type/ExtGState/CA 1/ca %s>>" % ("1", "0.05")[i])
            if where == "direct":
                text = doc.xref_object(rnum, compressed=True)
                inner = re.sub(r"/ExtGState\s+\d+ 0 R",
                               f"/ExtGState<</A1 {state} 0 R>>", text)
                assert inner != text, "fig5c's /Resources no longer names /ExtGState"
                doc.xref_set_key(doc[0].xref, "Resources", inner)
                assert doc.xref_get_key(doc[0].xref, "Resources")[0] == "dict", (
                    "the direct /Resources was not written as a dict, so this "
                    "control exercises the indirect path that already worked")
            else:
                gs = int(str(doc.xref_get_key(rnum, "ExtGState")[1]).split()[0])
                for k in doc.xref_get_keys(gs):
                    doc.xref_set_key(gs, f"{k}/ca", ("1", "0.05")[i])
                node = int(str(doc.xref_get_key(doc[0].xref, "Parent")[1]).split()[0])
                for _ in range(34):
                    up = doc.get_new_xref()
                    doc.update_object(
                        up, f"<</Type/Pages/Kids[{node} 0 R]/Count 1>>")
                    doc.xref_set_key(node, "Parent", f"{up} 0 R")
                    node = up
                doc.xref_set_key(node, "Resources", f"{rnum} 0 R")
                doc.xref_set_key(doc[0].xref, "Resources", "null")
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a, b = _drawing(arms[1]), _drawing(arms[0])
    assert len(a) == len(b), "the arms carry different surface lists"
    moved = {k for (k, x), (_, y) in zip(a, b) if x != y}
    assert moved and all(_head(k) == "extgstate" for k in moved), (
        f"[{where}] an alpha change is invisible because of HOW /Resources is "
        f"attached rather than what it says. Moved: {sorted(moved)}")


def test_image_surfaces_survive_a_renumbering(tmp_path):
    """CONTROL for `xref_name`, which had none.

    Replacing it with `str(xref)` left everything green, while its docstring
    claims image surfaces are keyed by the name the content stream uses so
    they are not allocation artifacts. A renumbering re-save then renames every
    image surface and the whole comparison reports stale -- loud, but a false
    stale is what this file's portability argument rests on avoiding.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig5c_mechanism_site_matrix.pdf"
    plain = tmp_path / "plain.pdf"
    copy = tmp_path / "c.pdf"
    copy.write_bytes(src.read_bytes())
    doc = pymupdf.open(copy)
    try:
        doc.save(plain, deflate=True)
    finally:
        doc.close()

    # Force different object numbers for the same drawing. The renumbering
    # comes from `garbage=4` compacting the table, NOT from the filler objects
    # -- measured, dropping them leaves the numbers just as different -- and
    # they are kept only to widen the shift. Said accurately because the first
    # version of this comment credited the fillers, and a control whose stated
    # mechanism is wrong is one edit from becoming a control that does nothing.
    shifted = tmp_path / "shifted.pdf"
    doc = pymupdf.open(copy)
    try:
        for n in range(20):
            x = doc.get_new_xref()
            doc.update_object(x, f"<</Type/Filler/N {n}>>")
        doc.save(shifted, deflate=True, garbage=4)
    finally:
        doc.close()

    def image_xrefs(path):
        d = pymupdf.open(path)
        try:
            return sorted(i[0] for i in d[0].get_images(full=True))
        finally:
            d.close()

    # FAIL, do not skip. This was a `pytest.skip`, so dropping `garbage=4` --
    # the control's own stated mechanism -- turned it into `56 passed,
    # 1 skipped`, exit 0, with `xref_name` then replaceable by `str(xref)` and
    # the suite green. The flagship's own comment two hundred lines up argues
    # exactly this: "NO SKIPS ... the check reported green while never
    # running." This was the last instance of it in the file.
    assert image_xrefs(plain) != image_xrefs(shifted), (
        "the objects were not renumbered, so this control cannot demonstrate "
        "that image surfaces are keyed on the resource name rather than the "
        "xref -- and it would have reported green while proving nothing")
    assert _drawing(plain) == _drawing(shifted), (
        "renumbering the objects of an identical drawing changed the surface "
        "list, so image surfaces are keyed on an allocation artifact")


def test_the_provenance_checks_can_actually_fail(tmp_path):
    """CONTROLS for the three guards the previous commit added.

    Each was individually deletable with the suite green: the stems assertion,
    the single-parent assertion, and `_COMPARED.append` recording the real
    path. That is this file's own signature defect -- "a guard with no control
    is how three surfaces shipped inert here" -- applied to the fix shipped
    FOR that defect, one commit later.
    """
    stems = _census_figures()
    good = [tmp_path / f"{s}.pdf" for s in stems]
    _assert_compared_the_census_set(good)          # the honest case passes

    # Eight calls naming ONE file satisfy a count exactly as well as eight
    # naming eight. This is the `_assert_matches_committed(produced[0])` edit.
    with pytest.raises(AssertionError, match="not the census set"):
        _assert_compared_the_census_set([good[0]] * len(stems))
    # A subset, which is what slicing the loop produces.
    with pytest.raises(AssertionError, match="not the census set"):
        _assert_compared_the_census_set(good[:-1])
    # The committed directory, which is what `scratch = FIG_DIR` produces.
    with pytest.raises(AssertionError, match="not article/figures"):
        _assert_compared_the_census_set(
            [FIG_DIR.resolve() / f"{s}.pdf" for s in stems])
    # And a mixture, which is what copying some figures in produces.
    with pytest.raises(AssertionError, match="one directory"):
        _assert_compared_the_census_set(
            good[:-1] + [FIG_DIR.resolve() / f"{stems[-1]}.pdf"])

    # `_COMPARED` must record the path actually handed to the comparison.
    copy = tmp_path / f"{stems[0]}.pdf"
    copy.write_bytes((FIG_DIR / copy.name).read_bytes())
    mark = len(_COMPARED)
    _assert_matches_committed(copy)
    assert _COMPARED[mark:] == [copy.resolve()], (
        "the comparison does not record the path it was given, so the "
        "provenance assertions above are checking a fabrication")


def test_the_scratch_snapshot_notices_a_repopulated_directory(tmp_path):
    """CONTROL for `_snapshot`, which had none.

    Replacing its body with `return {}` -- or with a name-keyed, content-blind
    map -- left the whole file green, while the assertion using it claims to
    notice the regenerated figures being overwritten. Both halves matter: it
    must see a content change under an unchanged NAME, which is exactly what
    copying a committed figure over a regenerated one looks like.
    """
    a, b = tmp_path / "one.pdf", tmp_path / "two.pdf"
    a.write_bytes(b"%PDF-1.4 first")
    b.write_bytes(b"%PDF-1.4 second")
    before = _snapshot(tmp_path)
    assert before == _snapshot(tmp_path), "_snapshot is not stable"
    a.write_bytes(b.read_bytes())          # same names, different contents
    assert _snapshot(tmp_path) != before, (
        "_snapshot cannot see a file being overwritten by another with the "
        "same name, which is the only thing it is there to notice")
    assert set(_snapshot(tmp_path)) == {"one.pdf", "two.pdf"}, (
        "_snapshot dropped a file")


@pytest.mark.parametrize("first,second", [
    ("/DeviceRGB", "/DeviceGray"),
    ("[/Indexed /DeviceRGB 1 <FF0000 00FF00>]", "[/Indexed /DeviceRGB 1 <0000FF FFFF00>]"),
])
def test_a_bare_name_resource_value_is_compared_by_mapping(tmp_path, first, second):
    """A resource value that is a NAME or an inline array, not a dictionary.

    `/Cs1 /DeviceRGB` is the ordinary spelling for a colour space, and the
    hand-written scanner this file used to walk resource dictionaries parsed it
    as an empty value plus a phantom entry called `DeviceRGB` -- so swapping
    WHICH name mapped to which space compared equal, the exact permutation the
    per-name surfaces exist to catch. Reading through pymupdf's parser instead
    fixes it, and this is what fails if anything goes back to reading text.
    """
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    src = FIG_DIR / "fig2c_census_volume.pdf"
    arms = []
    for i in (0, 1):
        copy = tmp_path / f"c{i}.pdf"
        copy.write_bytes(src.read_bytes())
        doc = pymupdf.open(copy)
        try:
            holder = doc.get_new_xref()
            a, b = (first, second) if i == 0 else (second, first)
            doc.update_object(holder, f"<</Cs1 {a}/Cs2 {b}>>")
            doc.xref_set_key(_page_resources_xref(doc), "ColorSpace",
                             f"{holder} 0 R")
            got = doc.xref_get_key(holder, "Cs1")
            assert "replace me" not in str(got), f"the value did not land: {got}"
            out = tmp_path / f"m{i}.pdf"
            doc.save(out, deflate=True)
            arms.append(out)
        finally:
            doc.close()

    a_, b_ = _drawing(arms[1]), _drawing(arms[0])
    assert len(a_) == len(b_), (
        f"the arms carry different surface lists: a name was invented or "
        f"dropped. {[k for k, _ in a_]} vs {[k for k, _ in b_]}")
    moved = {k for (k, x), (_, y) in zip(a_, b_) if x != y}
    assert moved == {"ColorSpace#Cs1", "ColorSpace#Cs2"}, (
        "swapping which NAME maps to which colour space moved "
        f"{sorted(moved) or 'nothing'}; the mapping is being compared as a set")


def _deep_probe(tmp: Path) -> Path:
    """A copy of fig2c whose `/Properties` hangs off a six-link chain."""
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf

    copy = tmp / "probe.pdf"
    copy.write_bytes((FIG_DIR / "fig2c_census_volume.pdf").read_bytes())
    doc = pymupdf.open(copy)
    try:
        head, _ = _chain(doc, "1", 6)
        holder = doc.get_new_xref()
        doc.update_object(holder, f"<</A {head} 0 R>>")
        doc.xref_set_key(_page_resources_xref(doc), "Properties", f"{holder} 0 R")
        out = tmp / "deep.pdf"
        doc.save(out, deflate=True)
        return out
    finally:
        doc.close()


def test_the_traversal_depth_committed_figures_need_is_measured():
    """`_MAX_DEPTH` is HEADROOM, and its justification used to be invented.

    The constant's comment claimed "chains of 11 exist in the test suite and 7
    was the measured blind spot". Neither number is real: instrumenting the
    traversal shows the eight committed figures reach depth 0 -- no resource on
    any of them is nested at all -- and the deepest crafted control reaches 6.

    So the number is derived here rather than asserted there, and this fails if
    a future figure ever gets close to the bound, which is the only condition
    under which the slack stops being slack.
    """
    global _MAX_DEPTH_SEEN
    _MAX_DEPTH_SEEN = 0
    for stem in _census_figures():
        assert _drawing(FIG_DIR / f"{stem}.pdf") is not None
    reached = _MAX_DEPTH_SEEN
    # The instrument first: a measurement of zero and a broken counter look
    # identical, and `_MAX_DEPTH_SEEN = 0` in place of the `max()` passed this
    # test unchanged. A crafted chain must move it.
    deep = _MAX_DEPTH_SEEN
    with tempfile.TemporaryDirectory() as td:
        probe = _deep_probe(Path(td))
        _MAX_DEPTH_SEEN = 0
        _drawing(probe)
        deep = _MAX_DEPTH_SEEN
    assert deep == 6, (
        f"the six-link probe registered depth {deep}, not 6; either the depth "
        "counter is not tracking -- in which case the measurement below is of "
        "the instrument rather than the file -- or the constant's comment, "
        "which states 6 exactly, is stale")
    assert reached == 0, (
        f"the committed figures now reach traversal depth {reached}, not 0; "
        "the constant's comment states 0 exactly, so update it together with "
        "the headroom check below")
    assert _MAX_DEPTH >= reached + 8, (
        f"_MAX_DEPTH ({_MAX_DEPTH}) leaves less than 8 levels of headroom "
        f"above the {reached} the committed figures actually need")


def test_the_png_and_pdf_output_sets_are_discovered_separately():
    """CONTROL for the split in `_census_outputs`, which had none.

    Folding both extensions into one set is green today, and it is the fix for
    a real defect: with one set, deleting a `.png` savefig leaves the count
    unchanged and the existence check -- which only ever looks for `{stem}.pdf`
    -- never notices, so a census PNG can stop being generated with the suite
    green. Two edits rather than one, and no control on either.
    """
    src = GEN.read_text()
    out = _census_outputs(src)
    assert set(out) == {"pdf", "png"}, f"the split is gone: {sorted(out)}"
    assert out["pdf"] and out["png"], "one of the extension sets is empty"
    # REMOVE a png savefig from a COPY of the source: the png set must shrink
    # and the pdf set must not. Folding the two into one makes both shrink or
    # neither, which is what leaves a deleted `.png` savefig invisible.
    stem = sorted(out["png"])[0]
    for spelling in (f'FIG_DIR / "{stem}.png"', f'FIG_DIR / f"{stem}.png"'):
        if spelling in src:
            doctored = src.replace(spelling, 'FIG_DIR / "_removed_.txt"')
            break
    else:
        raise AssertionError(
            f"{stem} is not written by a literal savefig, so this control "
            "cannot demonstrate the split")
    after = _census_outputs(doctored)
    assert stem not in after["png"], (
        "removing a png savefig did not shrink the png set")
    assert set(after["pdf"]) == set(out["pdf"]), (
        "removing a PNG savefig changed the PDF set, so the two extensions "
        "are discovered together and a deleted png savefig is invisible")


def test_the_byte_scan_this_replaced_is_blind_to_a_real_date():
    """The measured case that justifies `_has_date` opening a parser.

    Three sites here detected a creation date with `b"/CreationDate" in blob`.
    That is not a cheaper spelling of the same check -- it is a different and
    weaker one, and the difference is not hypothetical: it decided the backlog
    count this file states, which read seven while the truth was nine.

    Pinned against the committed blob at b8299eb8, the last commit before the
    dates were cleared, so this stays a measurement rather than a story. If the
    blob is ever unreachable the test fails rather than skipping: a positive
    control that quietly stops running is the failure mode it exists to
    prevent.
    """
    blob = subprocess.run(
        ["git", "show", "b8299eb8:article/figures/fig19_immune_coupling_flow.pdf"],
        cwd=REPO, capture_output=True)
    assert blob.returncode == 0, (
        "cannot read the pinned blob b8299eb8 -- CI needs fetch-depth: 0 for "
        f"this control to run: {blob.stderr.decode()[:200]}")
    data = blob.stdout
    assert data.startswith(b"%PDF"), "the pinned blob is not a PDF"
    assert b"CreationDate" not in data, (
        "the byte scan now FINDS a date in this blob, so it is no longer the "
        "example of the old check's blind spot and this control proves nothing")
    assert _has_date(data), (
        "the parser no longer reports a date for a blob that demonstrably "
        "carries one, so `_has_date` has stopped detecting the case the byte "
        "scan could not see")


def test_the_graphviz_path_behaves_when_dot_is_absent(tmp_path, monkeypatch):
    """Both halves of the graphviz fix are behavioural, and neither was pinned.

    Deleting the `shutil.which("dot")` guard -- restoring the crash this change
    exists to fix -- failed only `test_manifest_freshness.py`, i.e. the hash of
    the file, not its behaviour. Deleting the `_strip_pdf_date` call failed
    NOTHING. Measured with both mutants against the whole suite.

    Two claims, both cheap to hold:

      * with `dot` absent, `_render_graphviz` RETURNS FALSE rather than
        raising, and STRANDS NO `.gv` -- the pre-change code raised
        `FileNotFoundError` after fig18 and left `fig19_immune_coupling_flow.gv`
        untracked in `article/figures`;
      * `_strip_pdf_date` removes the date a real cairo file carries, checked
        against the committed blob rather than a synthetic one, so the case
        that motivated it is the case under test.
    """
    import importlib.util
    import sys

    # `scripts/` on the path first: the generator imports `figure_io` as a
    # sibling, so loading it by file location alone raises ModuleNotFoundError.
    sys.path.insert(0, str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        "_gcd", REPO / "scripts" / "generate_conceptual_diagrams.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_gcd"] = mod
    spec.loader.exec_module(mod)

    monkeypatch.setattr(mod.shutil, "which", lambda name: None)
    out_base = tmp_path / "figX_probe"
    rendered = mod._render_graphviz("digraph { a -> b }", out_base)

    assert rendered is False, (
        "_render_graphviz must report that it skipped when `dot` is absent, "
        "so main() can carry on and exit 0 instead of crashing part-way "
        "through the figure list")
    assert not list(tmp_path.iterdir()), (
        f"the graphviz path left {[p.name for p in tmp_path.iterdir()]} behind "
        "when `dot` was absent; the .gv is written before the check and must "
        "be cleaned up whether or not the render happens")

    # AND THE STRIP THROUGH ITS CALL SITE. Calling `_strip_pdf_date` directly
    # tests the helper and not the wiring: deleting the call from
    # `_render_graphviz` still passed. `dot` is faked here so the render path
    # runs end to end, and what it "renders" is the REAL cairo blob, so the
    # case under test is the case that motivated the fix.
    blob = subprocess.run(
        ["git", "show", "b8299eb8:article/figures/fig19_immune_coupling_flow.pdf"],
        cwd=REPO, capture_output=True).stdout
    assert blob, "could not read the committed cairo blob to test the strip"
    probe = tmp_path / "cairo.pdf"
    probe.write_bytes(blob)
    assert _has_date(probe), (
        "the fixture blob is supposed to be the dated cairo file; if it is not, "
        "this test is checking nothing")
    probe.unlink()

    def fake_dot(argv, **kw):
        out = Path(argv[argv.index("-o") + 1])
        out.write_bytes(blob if out.suffix == ".pdf" else b"\x89PNG\r\n")
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/dot")
    monkeypatch.setattr(mod.subprocess, "run", fake_dot)
    base2 = tmp_path / "figY_probe"
    assert mod._render_graphviz("digraph { a -> b }", base2) is True
    rendered = Path(str(base2) + ".pdf")
    assert rendered.exists(), "the faked render wrote no pdf"
    assert not _has_date(rendered), (
        "_render_graphviz returned a PDF that still carries a creation date, "
        "so the strip is not wired into the path that produces these figures")
