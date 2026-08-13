"""PDF text extraction must not drift silently under a PyMuPDF bump.

Two scripts turn PDFs into the corpus text every downstream number rests on --
`fetch_articles.py` and `recover_fulltext.py` -- and `search_books.py` reads
PDFs the same way to answer queries. All three go through one call,
`page.get_text("text")`. Until this file, no test imported any of them, so the
extraction path that produces the corpus had zero coverage.

The gap was found while bumping PyMuPDF 1.27.2.3 -> 1.28.0, which also moves the
bundled engine from MuPDF 1.27.2 to 1.29.0. The failure mode is silent: a shift
in whitespace, ligature handling or reading order does not raise, it just means
text fetched after the bump is not comparable with the files already committed.
Nothing red, a heterogeneous corpus.

WHAT WAS MEASURED, so this is a bound rather than a worry
----------------------------------------------------------
`fetch_articles.extract_text_from_pdf_bytes` was run on the fixture under BOTH
versions in isolated virtualenvs, and the output is byte-identical -- 1,522
characters including reading order and the page join, same SHA-256. So that
particular bump moved nothing on these shapes. These tests pin it forward.

WHY THE FIXTURE IS 1,522 CHARACTERS AND NOT 241
------------------------------------------------
Both extractors end with `return joined if len(joined) >= 1000 else None`. The
first version of this fixture was 241 characters, which sits inside the region
production code DISCARDS -- so the real functions could not be called on it, and
the tests instead asserted against a REIMPLEMENTATION of the call sites. Review
proved that vacuous: deleting the `if text:` skip from `fetch_articles.py` left
all of them green, because they were exercising the copy rather than the
original. The fixture is now long enough for the real functions to accept it,
and they are what is called below.

WHAT THE FIXTURE DOES AND DOES NOT COVER
-----------------------------------------
Covers: the real extractors end to end, multi-page joins, the empty-page skip,
line breaks within a page, and content-stream-versus-position ordering.

Does NOT cover: embedded or subset fonts, ligatures, CID encodings, columns,
tables or scanned images. The fixture is base-14 Helvetica and ASCII, because it
is authored byte-by-byte from stdlib to keep its provenance unambiguous. Real
papers carry all of the above, and that is where extraction risk actually
concentrates -- so a green run here does NOT license the claim that a bump is
safe for the corpus, only that it did not move these shapes.
"""

import hashlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

# A plain import, deliberately not `pytest.importorskip`. PyMuPDF is pinned in
# requirements-lock.txt and CI installs it, so a skip could never be legitimate,
# and this gate must go red rather than quietly skip if the import ever breaks.
#
# THE ALIAS THIS WARNED ABOUT HAS NOW BEEN DEPRECATED. An earlier version of this
# comment said PyMuPDF "has been steering users from `fitz` to `pymupdf`" and that
# all three scripts would break outright if the alias were dropped. 1.28.2 began
# emitting "The `fitz` API is deprecated and will be removed in future", so the
# scripts and this file were migrated to `import pymupdf` ahead of the removal.
# The migration was verified to change nothing: the real extractors return
# byte-identical text on the committed fixture before and after (1,522 chars,
# sha256 1449f164...), which is the only property that matters, since these
# functions are how PDFs become corpus text.
import pymupdf  # noqa: E402
from fetch_articles import extract_text_from_pdf_bytes as fetch_extract  # noqa: E402
from recover_fulltext import extract_text_from_pdf_bytes as recover_extract  # noqa: E402

FIXTURE = REPO_ROOT / "tests" / "fixtures" / "extraction-sample.pdf"

# Pinned so the fixture cannot be edited without the expectation being revisited.
# Without it, "regenerate the fixture until the test passes" silently becomes an
# option and the guard would assert only that two things agree.
FIXTURE_SHA256 = "498366c3111f83a0bde5a1859936ebcdd80e151c9f4a700ede52abbdce499345"

# SHA-256 of the extracted text. Identical under PyMuPDF 1.27.2.3 (MuPDF 1.27.2)
# and 1.28.0 (MuPDF 1.29.0), measured in isolated virtualenvs.
EXTRACTED_SHA256 = "1449f164601ede80"     # 16-char prefix; full value asserted via length + anchors
EXTRACTED_LEN = 1522

# The lines whose ORDER is the interesting property. `get_text("text")` follows
# content-stream order, not visual reading order: on page 2 the block that is
# LOWER on the page is emitted FIRST because it is drawn first. Confirmed causal
# -- flipping only the draw order in the generator, at identical coordinates,
# flips the extraction -- and `get_text("text", sort=True)` gives the opposite.
# If a future engine starts sorting by position, this is the assertion that says
# so, and the corpus would need re-extracting.
DRAWN_FIRST = "Second block, drawn first in the content stream."
HIGHER_ON_PAGE = "First block by position, drawn second."


def _text() -> str:
    return fetch_extract(FIXTURE.read_bytes())


def test_the_fixture_itself_has_not_been_edited():
    assert FIXTURE.exists(), "the extraction fixture is missing"
    got = hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    assert got == FIXTURE_SHA256, (
        "the fixture changed; the expectations below must be re-derived "
        "deliberately rather than by regenerating until green "
        "(see scripts/make_extraction_fixture.py)")


def test_the_production_extractor_output_is_unchanged():
    """The REAL function, not a copy of it -- see the module docstring."""
    got = _text()
    assert got is not None, (
        "the production extractor rejected the fixture; it is probably back "
        "under the 1,000-character floor")
    assert len(got) == EXTRACTED_LEN, f"extracted {len(got)} chars, expected {EXTRACTED_LEN}"
    assert hashlib.sha256(got.encode()).hexdigest().startswith(EXTRACTED_SHA256), (
        "extracted text changed; if this followed a PyMuPDF bump, the corpus "
        "fetched before and after it is no longer comparable")


def test_both_corpus_extractors_agree():
    """`fetch_articles` and `recover_fulltext` hold the same code separately.

    They are independent copies, so they can drift apart; the corpus is written
    by whichever ran, and nothing else would notice.
    """
    assert recover_extract(FIXTURE.read_bytes()) == _text()


def test_extraction_follows_content_stream_order_not_visual_order():
    got = _text()
    assert got.index(DRAWN_FIRST) < got.index(HIGHER_ON_PAGE), (
        "extraction is now position-sorted rather than content-stream ordered; "
        "every corpus record extracted before this change has a different "
        "reading order from every record extracted after it")


def test_empty_pages_are_skipped_not_joined_as_blanks():
    """Both extractors guard with `if text:`; the fixture has an empty page.

    Asserted on the REAL output: without the skip the join would emit a second
    blank-line separator, so a corpus record would carry structure that was
    never in the paper.
    """
    with pymupdf.open(stream=FIXTURE.read_bytes(), filetype="pdf") as doc:
        per_page = [p.get_text("text").strip() for p in doc]
    assert len(per_page) == 3, "the fixture should have three pages"
    assert per_page[1] == "", (
        "the empty page must be the MIDDLE one. As the last page it cannot "
        "detect the skip at all: without `if text:` the join appends a "
        "separator that the extractor's own closing .strip() removes, so the "
        "output is byte-identical either way -- measured, not supposed")
    assert _text().count("\n\n") == 1, (
        "two non-empty pages and one empty one should produce exactly one "
        "blank-line separator")


def test_the_path_open_shape_agrees_with_the_bytes_shape():
    """`search_books.py` opens by path where the extractors open from bytes."""
    with pymupdf.open(filename=str(FIXTURE)) as doc:
        by_path = "\n\n".join(
            t for t in (p.get_text("text").strip() for p in doc) if t).strip()
    assert by_path == _text()


def test_the_generator_reproduces_the_committed_fixture():
    """The fixture is authored, so its generator must still produce it.

    Deliberately NOT built with PyMuPDF: if the same library wrote and read the
    fixture, a change affecting both directions could cancel and the guard would
    pass through the drift it exists to catch.
    """
    import make_extraction_fixture as mk

    assert hashlib.sha256(mk.build()).hexdigest() == FIXTURE_SHA256, (
        "the generator no longer reproduces the committed fixture")


def test_every_pdf_reading_module_is_covered_here():
    """The inverse check: a NEW script that reads PDFs must not be invisible.

    Asserting that three named files use `get_text` would keep passing while a
    fourth extractor appeared beside them, unguarded. So the direction is
    reversed -- whatever imports PyMuPDF has to be accounted for. The match
    accepts BOTH `import fitz` and `import pymupdf`, which is what let the
    migration off the deprecated alias happen without blinding this check; a
    guard keyed to only the old spelling would have gone silently empty the
    moment the rename landed.
    """
    known = {"fetch_articles.py", "recover_fulltext.py", "search_books.py",
             "make_extraction_fixture.py"}
    importers = {p.name for p in (REPO_ROOT / "scripts").glob("*.py")
                 if "import fitz" in p.read_text() or "import pymupdf" in p.read_text()}
    assert importers <= known, (
        f"{sorted(importers - known)} read PDFs but are not covered by this "
        "gate; either add them here or note why they cannot drift the corpus")
    for name in importers & {"fetch_articles.py", "recover_fulltext.py",
                             "search_books.py"}:
        src = (REPO_ROOT / "scripts" / name).read_text()
        assert 'get_text("text")' in src, (
            f"{name} no longer uses get_text(\"text\"); this file pins a call it "
            "does not make, so the extraction path is unguarded again")
