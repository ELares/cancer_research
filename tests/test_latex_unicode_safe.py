"""Guard: the generated LaTeX manuscript contains no pdflatex-unsafe Unicode.

`pdflatex` with `inputenc=utf8` maps a LIMITED set of non-ASCII code points to
glyphs (accented Latin, +/-, curly quotes, em/en dashes, ...). A code point
OUTSIDE that set — most notably any **Greek letter** (U+0370-03FF), e.g. a bare
`Σ` — is a hard `! LaTeX Error: Unicode character ...`. pdflatex recovers and
still emits a PDF, but it exits NON-ZERO, which broke the `release-pdf` workflow
(a literal `Σ` in the sensitivity paragraph; see scripts/generate_latex.py
`map_greek`). `scripts/generate_latex.py` now maps every Greek letter to a math
command, so the committed `article/drafts/v1.tex` must be free of bare Greek.

This catches a recurrence on the Python CI (PR-time) BEFORE it reaches the
post-merge PDF build: either a stale `v1.tex` or a generator regression that
stops mapping some Greek letter.
"""

import unicodedata
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
V1_TEX = REPO_ROOT / "article" / "drafts" / "v1.tex"

# Non-ASCII code points that `inputenc`/`fontenc` render WITHOUT error. Everything
# else in v1.tex must have been converted to a LaTeX command by generate_latex.py.
SAFE_NON_ASCII = {
    "±",  # ± plus-minus
    "—",  # — em dash
    "–",  # – en dash
    "“",  # " left double quote
    "”",  # " right double quote
    "‘",  # ' left single quote
    "’",  # ' right single quote
}
# Accented Latin (Latin-1 Supplement + Latin Extended-A) is handled by fontenc T1.
ACCENTED_LATIN_RANGES = ((0x00C0, 0x017F),)

GREEK_BLOCK = (0x0370, 0x03FF)


def _is_safe(ch: str) -> bool:
    cp = ord(ch)
    if cp < 0x80 or ch in SAFE_NON_ASCII:
        return True
    return any(lo <= cp <= hi for lo, hi in ACCENTED_LATIN_RANGES)


def test_v1_tex_is_committed_and_nonempty():
    assert V1_TEX.exists(), f"{V1_TEX} missing; run scripts/generate_latex.py"
    assert V1_TEX.stat().st_size > 100_000, "v1.tex unexpectedly small"


def test_v1_tex_has_no_bare_greek():
    """The exact class that broke release-pdf: a bare Greek code point."""
    text = V1_TEX.read_text(encoding="utf-8")
    greek = sorted(
        {ch for ch in text if GREEK_BLOCK[0] <= ord(ch) <= GREEK_BLOCK[1]}
    )
    assert not greek, (
        "Bare Greek letters in article/drafts/v1.tex break pdflatex "
        "(`! LaTeX Error: Unicode character ...`): "
        + ", ".join(f"U+{ord(ch):04X} {ch!r}" for ch in greek)
        + ". Re-run `python3 scripts/generate_latex.py` (its map_greek covers all "
        "Greek); if a NEW Greek letter is unmapped, add it to _GREEK_TO_LATEX."
    )


def test_v1_tex_has_no_other_pdflatex_unsafe_unicode():
    """Any non-ASCII char outside the inputenc/fontenc-safe set is a risk."""
    text = V1_TEX.read_text(encoding="utf-8")
    bad = {}
    for ch in text:
        if not _is_safe(ch):
            bad.setdefault(ch, 0)
            bad[ch] += 1
    if bad:
        detail = ", ".join(
            f"U+{ord(ch):04X} {ch!r} ({unicodedata.name(ch, '?')}) x{n}"
            for ch, n in sorted(bad.items(), key=lambda kv: -kv[1])
        )
        raise AssertionError(
            "pdflatex-unsafe non-ASCII in v1.tex (add a mapping in "
            f"scripts/generate_latex.py and regenerate): {detail}"
        )
