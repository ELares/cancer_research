#!/usr/bin/env python3
"""Build the PDF text-extraction fixture (#PDF-FIXTURE).

WHY THIS EXISTS
---------------
Three scripts turn PDFs into the corpus text every downstream number is
computed from -- `fetch_articles.py` and `recover_fulltext.py` -- and `search_books.py` reads
PDFs the same way to answer queries. All three go through one call,
`page.get_text("text")`.
No test imported any of them. When PyMuPDF was bumped 1.27.2.3 -> 1.28.0
(bundling MuPDF 1.29.0) nothing in CI could have noticed if extraction had
changed, because nothing extracted anything.

That is a silent failure mode rather than a loud one. A shift in whitespace,
ligature handling or reading order does not raise; it just means text fetched
after the bump is not comparable with the 4,884 files already committed.

WHY THE FIXTURE IS AUTHORED HERE RATHER THAN DOWNLOADED
--------------------------------------------------------
A real paper would carry a licence and a provenance question, and this repo
tracks both (`PROVENANCE.yaml`). This PDF is written byte-by-byte from stdlib
with no library, so it is unambiguously the project's own, small enough to read,
and stable enough to pin.

Deliberately NOT generated with PyMuPDF: if the same library both wrote and read
the fixture, a change affecting both directions could cancel out and the guard
would pass through exactly the drift it exists to catch.

WHAT IT EXERCISES, and why each part is here
---------------------------------------------
  * THREE PAGES, two with text, because the extractors join pages with "\\n\\n" and a
    change in per-page stripping shows up only across a boundary.
  * SEVERAL LINES PER PAGE, because line breaks and inter-line spacing are the
    commonest thing to move between engine versions.
  * TEXT PLACED OUT OF READING ORDER -- page 2 draws its lower block FIRST in
    the content stream. An extractor that returns content-stream order gives a
    different answer from one that sorts by position, and that is precisely the
    "reading order" risk a MuPDF minor carries.
  * AN EMPTY PAGE IN THE MIDDLE, not at the end. Both extractors skip empty
    pages with `if text:`, and a TRAILING empty page cannot detect that skip:
    without it the join appends a separator that the extractor's own closing
    `.strip()` immediately removes, so the output is identical either way. The
    first version of this fixture put the empty page last and a mutation that
    deleted the skip from the real `fetch_articles.py` left every test green.
    Between two text pages the missing skip shows up as a doubled separator.
  * ENOUGH TEXT TO CLEAR 1,000 CHARACTERS. Both extractors end with
    `return joined if len(joined) >= 1000 else None`. A shorter fixture sits
    inside the region production code DISCARDS, so the real functions could not
    be called on it at all -- which is how the first version of this fixture
    ended up guarded by a reimplementation of the call sites instead of by the
    call sites themselves.

Usage:
    python scripts/make_extraction_fixture.py
"""

import pathlib
import sys

OUT = pathlib.Path(__file__).resolve().parent.parent / "tests" / "fixtures" / "extraction-sample.pdf"

# (x, y, text) in PDF user space, origin bottom-left. Page 2's blocks are listed
# bottom-first ON PURPOSE: see the module docstring.
PAGE_1 = [(72, 720, "Ferroptosis in cancer therapy"),
          (72, 700, "A fixture document for extraction drift."),
          (72, 670, "Lipid peroxidation drives an iron dependent death."),
          (72, 650, "GPX4 is the canonical defence against that peroxidation."),
          (72, 630, "System Xc minus imports cystine to sustain glutathione."),
          (72, 610, "Erastin inhibits that transporter and starves the cell."),
          (72, 590, "FSP1 and DHODH repair peroxides in parallel to GPX4."),
          (72, 570, "ACSL4 sets how much oxidisable substrate the membrane holds."),
          (72, 550, "Labile iron drives the Fenton chemistry that propagates."),
          (72, 530, "Ferritinophagy releases stored iron into that labile pool."),
          (72, 510, "PROM2 exports it again, draining the pool and resisting."),
          (72, 490, "The balance of those two decides whether a cell dies."),
          (72, 470, "None of these sentences is a claim about the literature."),
          (72, 450, "They exist to give the extractor enough text to work on.")]
PAGE_2 = [(72, 600, "Second block, drawn first in the content stream."),
          (72, 700, "First block by position, drawn second."),
          (72, 560, "Padding so the extracted document clears the one thousand"),
          (72, 540, "character floor that the production extractors apply before"),
          (72, 520, "they accept a PDF at all, since a shorter fixture would sit"),
          (72, 500, "inside the region that real code discards and could not be"),
          (72, 480, "passed through the real function under test."),
          (72, 450, "Trailing lines keep the page count and ordering intact."),
          (72, 430, "The reading order property above is unaffected by them."),
          (72, 410, "Each line is plain ASCII in the base fourteen Helvetica."),
          (72, 390, "No ligature, no embedded font, no column layout is present."),
          (72, 370, "Those shapes are named in the test as explicitly uncovered."),
          (72, 350, "A real paper carries all of them and this one does not."),
          (72, 330, "That limit is the honest boundary of what this guard proves.")]
PAGE_EMPTY = []      # deliberately empty: exercises the `if text:` skip


def _content(lines) -> bytes:
    """A content stream drawing each line in 12pt Helvetica."""
    if not lines:
        return b""
    out = ["BT", "/F1 12 Tf"]
    for x, y, text in lines:
        esc = text.replace("\\", r"\\").replace("(", r"\(").replace(")", r"\)")
        out.append(f"1 0 0 1 {x} {y} Tm ({esc}) Tj")
    out.append("ET")
    return ("\n".join(out) + "\n").encode("ascii")


def build() -> bytes:
    pages = [PAGE_1, PAGE_EMPTY, PAGE_2]
    n_pages = len(pages)
    # Object numbering: 1 catalog, 2 pages tree, 3 font,
    # then per page: page object and content stream.
    page_ids = [4 + 2 * i for i in range(n_pages)]
    objects = {
        1: b"<< /Type /Catalog /Pages 2 0 R >>",
        2: ("<< /Type /Pages /Count %d /Kids [%s] >>"
            % (n_pages, " ".join(f"{i} 0 R" for i in page_ids))).encode("ascii"),
        3: b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    }
    for i, lines in enumerate(pages):
        pid, cid = page_ids[i], page_ids[i] + 1
        objects[pid] = (
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] "
            "/Resources << /Font << /F1 3 0 R >> >> /Contents %d 0 R >>"
            % cid).encode("ascii")
        body = _content(lines)
        objects[cid] = (b"<< /Length %d >>\nstream\n" % len(body)) + body + b"endstream"

    out = bytearray(b"%PDF-1.4\n")
    offsets = {}
    for num in sorted(objects):
        offsets[num] = len(out)
        out += b"%d 0 obj\n" % num + objects[num] + b"\nendobj\n"
    xref_at = len(out)
    top = max(objects) + 1
    out += b"xref\n0 %d\n" % top
    out += b"0000000000 65535 f \n"
    for num in range(1, top):
        out += b"%010d 00000 n \n" % offsets[num]
    out += (b"trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n"
            % (top, xref_at))
    return bytes(out)


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    data = build()
    OUT.write_bytes(data)
    print(f"wrote {OUT} ({len(data):,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
