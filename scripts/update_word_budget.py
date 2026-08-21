#!/usr/bin/env python3
"""Rewrite book-outline.md's word-budget table from the manuscript.

REPLACES A SHELL HELPER THAT FALSIFIED THE TABLE. That helper kept the Total
correct by attributing every delta to Part II, so edits made in Chapters 7 and
10 and in the glossary were booked against Chapters 3-4 -- and the Total-only
guard stayed green throughout. book-outline.md's own note had already warned
about exactly this ("an earlier revision changed only the Total and left the
per-part rows stale") and nothing enforced it.

Every row is measured from its own section of the manuscript. Nothing is
derived from a delta, so there is no arithmetic in which a number can be booked
against the wrong row.

Counting and rounding are IMPORTED from the guard rather than reimplemented,
because a helper that counts words its own way will eventually disagree with
the test it exists to satisfy, and the disagreement would look like drift in
the manuscript.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUTLINE = REPO / "article/book-outline.md"
sys.path.insert(0, str(REPO / "tests"))
from test_book_outline_wordcount import (  # noqa: E402
    ROUND_TO, _part_word_counts, prose_words, MANUSCRIPT,
)


def main() -> int:
    measured = _part_word_counts()
    total = prose_words(MANUSCRIPT.read_text())

    # ROUNDING EACH ROW INDEPENDENTLY DOES NOT SUM. Every figure in the table
    # is rounded to the nearest hundred, and six independently-rounded rows can
    # miss the rounded total by up to three hundred -- which fails the
    # pre-existing invariant that the table adds up. A table whose rows do not
    # sum to its own Total is reporting two different documents, so the
    # residual is allocated by LARGEST REMAINDER: the rows whose exact values
    # sat closest to a rounding boundary absorb it, one step each. Every row
    # stays within one step of its measurement and the column adds up.
    step = 10 ** -ROUND_TO
    rounded = {k: round(v, ROUND_TO) for k, v in measured.items()}
    residual = round(total, ROUND_TO) - sum(rounded.values())
    if residual:
        direction = 1 if residual > 0 else -1
        # Distance to the boundary in the direction we must move.
        order = sorted(
            measured,
            key=lambda k: -direction * (measured[k] - rounded[k]))
        for k in order[:abs(residual) // step]:
            rounded[k] += direction * step
        assert sum(rounded.values()) == round(total, ROUND_TO), (
            "the residual allocation did not close the gap")
    measured = rounded
    s = OUTLINE.read_text()
    start = s.index("| Part | Chapters | Current words | Target words |")
    end = s.index("\n\n", start)
    table = s[start:end]
    lines = []
    changed = []
    for line in table.split("\n"):
        cells = [c.strip() for c in line.split("|")]
        if len(cells) < 5 or cells[1].startswith("---") or cells[1] == "Part":
            lines.append(line)
            continue
        label = cells[1].strip("* ")
        if label == "Total":
            new = round(total, ROUND_TO)
        elif label in measured:
            new = measured[label]        # already rounded and reconciled
        else:
            print(f"  ! no manuscript section for row {label!r}", file=sys.stderr)
            lines.append(line)
            continue
        old_cell = cells[3]
        bold = old_cell.startswith("**")
        new_cell = (f"**~{new:,}**" if bold else f"~{new:,}")
        if new_cell != old_cell:
            changed.append(f"{label}: {old_cell} -> {new_cell}")
        cells[3] = new_cell
        lines.append("| " + " | ".join(cells[1:-1]) + " |")
    OUTLINE.write_text(s[:start] + "\n".join(lines) + s[end:])
    import re

    # THE MANUSCRIPT'S OWN FRONT MATTER quotes the same total, and it is the
    # first factual statement a reader meets. It had drifted to ~46,500 against
    # a real 53,400 -- a 15% understatement -- because this helper updated the
    # outline and nothing updated the document the outline describes.
    #
    # AND EVERY OTHER SITE THAT QUOTES IT. The count appears in the manuscript's
    # front matter, the README and CLAUDE.md, and all three had drifted to
    # different wrong values -- 46,500 in one and 48,500 in two others against a
    # real 53,500. Updating one by hand is how three sites end up disagreeing;
    # a helper that owns all of them is why they cannot.
    want = f"~{round(total, ROUND_TO):,} words"
    for path in (MANUSCRIPT, REPO / "README.md", REPO / "CLAUDE.md"):
        if not path.exists():
            continue
        txt = path.read_text()
        new_txt, k = re.subn(r"~[\d,]+ words", want, txt)
        if k and new_txt != txt:
            path.write_text(new_txt)
            print(f"  {path.name} -> {want} ({k} site(s))")

    # The **Target:** line quotes the same total and goes stale on its own.
    s = OUTLINE.read_text()
    import re
    s2, n = re.subn(r"current: ~[\d,]+", f"current: ~{round(total, ROUND_TO):,}", s)
    if n:
        OUTLINE.write_text(s2)
    print(f"measured total {total:,} (rounds to {round(total, ROUND_TO):,})")
    for c in changed:
        print(f"  {c}")
    if not changed:
        print("  every row already matched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
