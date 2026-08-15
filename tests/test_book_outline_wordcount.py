"""The book outline's word budget must add up to itself and to the manuscript.

WHY THIS EXISTS
---------------
`article/book-outline.md` carries a per-part word-budget table and a Total row.
A revision changed the Total from ~39,400 to ~48,400 and left all five per-part
rows untouched, so the table stopped summing to its own bottom line by 9,300
words -- and nothing failed, because no procedure in the repository produced
either number. The figure was described as "measured prose-only from v1.md"
while being reproducible by no committed measurement.

That is the shape this repo has been burned by more than once: a document that
computes some of its figures and hand-writes the sentence beside them reads as
though the whole line was measured. Two properties are pinned here.

1. THE TABLE SUMS TO ITSELF. Rows plus front matter equal the Total row. This
   is arithmetic on the document alone and cannot go stale.

2. THE TOTAL MATCHES THE MANUSCRIPT. Recomputed from v1.md at the same rounding
   the document uses, so editing the manuscript without updating the outline
   fails here rather than silently drifting.

WHAT "PROSE-ONLY" MEANS, exactly, since the phrase is doing the work: fenced
code blocks removed, then heading lines and markdown table rows dropped. A raw
`wc -w` over v1.md gives a different (larger) answer, so the definition has to
live somewhere executable rather than in an adjective.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTLINE = REPO_ROOT / "article" / "book-outline.md"
MANUSCRIPT = REPO_ROOT / "article" / "drafts" / "v1.md"

# The document rounds to the nearest hundred; the guard must round the same way
# or it would fail on a difference the document cannot express.
ROUND_TO = -2


def prose_words(markdown: str) -> int:
    """Words outside fenced code, heading lines and table rows."""
    body = re.sub(r"```.*?```", "", markdown, flags=re.S)
    kept = [
        line for line in body.split("\n")
        if not line.startswith("#") and not line.strip().startswith("|")
    ]
    return len(" ".join(kept).split())


def _budget_rows() -> "list[tuple[str, int]]":
    """(label, current-words) for every row of the word-budget table."""
    text = OUTLINE.read_text()
    start = text.index("| Part | Chapters | Current words | Target words |")
    end = text.index("\n\n", start)
    rows = []
    for line in text[start:end].split("\n")[2:]:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        found = re.search(r"([\d,]+)", cells[2])
        if found:
            rows.append((cells[0].replace("*", ""), int(found.group(1).replace(",", ""))))
    return rows


def test_the_budget_table_sums_to_its_own_total():
    """The failure that shipped: Total edited, rows left alone."""
    rows = _budget_rows()
    assert len(rows) >= 3, f"only parsed {len(rows)} budget rows; the table moved"
    totals = [(label, n) for label, n in rows if "Total" in label]
    assert len(totals) == 1, f"expected exactly one Total row, parsed {totals}"
    total = totals[0][1]
    parts = [(label, n) for label, n in rows if "Total" not in label]
    summed = sum(n for _, n in parts)
    assert summed == total, (
        f"the word-budget rows sum to {summed:,} but the Total row says "
        f"{total:,}, a {abs(total - summed):,}-word disagreement. Rows: "
        + ", ".join(f"{label}={n:,}" for label, n in parts))


def test_the_documented_total_matches_the_manuscript():
    """`measured prose-only from v1.md` has to be a measurement."""
    measured = prose_words(MANUSCRIPT.read_text())
    documented = [n for label, n in _budget_rows() if "Total" in label][0]
    assert round(measured, ROUND_TO) == documented, (
        f"book-outline.md documents {documented:,} words but v1.md measures "
        f"{measured:,} (rounds to {round(measured, ROUND_TO):,}). Update the "
        "Total row AND the per-part rows, then the **Target:** line.")


def test_the_headline_target_line_agrees_with_the_table():
    """The intro line and the table's Total target were 4,000 words apart."""
    text = OUTLINE.read_text()
    head = re.search(r"\*\*Target:\*\* ~([\d,]+) words \(current: ~([\d,]+)", text)
    assert head, "the **Target:** line no longer states a target and a current"
    target, current = (int(g.replace(",", "")) for g in head.groups())
    rows = _budget_rows()
    documented = [n for label, n in rows if "Total" in label][0]
    assert current == documented, (
        f"the **Target:** line says the manuscript is at {current:,} words "
        f"while the table's Total says {documented:,}")
    # the target column's own Total, parsed from the same row
    line = next(l for l in text.split("\n") if "**Total**" in l)
    target_cell = [c.strip() for c in line.strip().strip("|").split("|")][-1]
    table_target = int(re.search(r"([\d,]+)", target_cell).group(1).replace(",", ""))
    assert target == table_target, (
        f"the **Target:** line targets {target:,} words while the table's "
        f"Total target column says {table_target:,}")


def test_prose_only_is_not_the_same_as_a_raw_word_count():
    """Pins that the definition is doing real work.

    If `prose_words` ever degenerates to splitting the whole file, the two
    measurements above would still agree with each other while measuring
    something the document does not claim.
    """
    raw = MANUSCRIPT.read_text()
    assert prose_words(raw) < len(raw.split()), (
        "prose_words returns the raw word count, so the 'prose-only' "
        "qualifier is no longer excluding anything")
