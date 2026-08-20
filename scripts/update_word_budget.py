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
            new = round(measured[label], ROUND_TO)
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
