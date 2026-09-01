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


# --- per-part rows, not only the Total ------------------------------------

PART_ROW_TO_HEADING = {
    "I: Why This Exists": "# Part I: Why This Exists",
    "II: What We Found": "# Part II: What We Found",
    "III: Simulations": "# Part III: What The Simulations Show",
    "IV: What's Next": "# Part IV: What Should Happen Next",
    "V: References/Tools": "# Part V: References and Tools",
    "front matter (title, abstract)": None,   # everything before Part I
}


def _part_word_counts() -> "dict[str, int]":
    lines = MANUSCRIPT.read_text().split("\n")
    starts = {}
    for i, line in enumerate(lines):
        for row, heading in PART_ROW_TO_HEADING.items():
            if heading is not None and line.strip() == heading:
                starts[row] = i
    assert len(starts) == len(PART_ROW_TO_HEADING) - 1, (
        f"manuscript part headings changed; found {sorted(starts)}")
    ordered = sorted(starts.items(), key=lambda kv: kv[1])
    out = {"front matter (title, abstract)": prose_words(
        "\n".join(lines[: ordered[0][1]]))}
    for k, (row, i) in enumerate(ordered):
        end = ordered[k + 1][1] if k + 1 < len(ordered) else len(lines)
        out[row] = prose_words("\n".join(lines[i:end]))
    return out


def test_every_per_part_row_matches_the_manuscript():
    """The Total alone is not enough, and this file already knew it.

    The outline's own note says the per-part figures are measured rather than
    planned "because an earlier revision changed only the Total and left the
    per-part rows stale". That is exactly what happened again: a helper kept
    the Total correct by attributing every delta to Part II, so edits made in
    Chapters 7 and 10 and in the glossary were booked against Chapters 3-4.
    The Total stayed green throughout.

    A row is allowed to be stale by less than the rounding step, since every
    figure in the table is rounded to the nearest hundred.
    """
    measured = _part_word_counts()
    documented = dict(_budget_rows())
    step = 10 ** -ROUND_TO
    # ONE FULL STEP, not half. Every figure is rounded to the nearest hundred
    # and six independently-rounded rows need not sum to the rounded total, so
    # the generator allocates the residual by largest remainder -- which can
    # move a row a full step away from its own measurement.
    #
    # That looks like a loosened guard and is not, because it does not stand
    # alone: the sum tests require the rows to add up EXACTLY. Booking one
    # part's edit against another has to break one or the other, since moving
    # words between rows changes two of them while the total stays put. The
    # defect this was written for -- a helper crediting every delta to Part II
    # -- fails on the first row it touches.
    bad = []
    for row, words in measured.items():
        if row not in documented:
            bad.append(f"{row}: no row in the budget table")
            continue
        if abs(documented[row] - words) > step:
            bad.append(f"{row}: documented {documented[row]:,} vs measured "
                       f"{words:,} (rounds to {round(words, ROUND_TO):,})")
    assert not bad, (
        "book-outline.md's per-part rows disagree with the manuscript:\n  "
        + "\n  ".join(bad)
        + "\nUpdate the rows the edit actually landed in, not whichever row "
          "makes the Total add up.")


def test_the_per_part_rows_sum_to_the_documented_total():
    """A table whose rows do not add up to its own Total is reporting two
    different documents. Checked against the DOCUMENTED figures rather than
    the measured ones, so it fails on a bookkeeping error even when every
    individual row happens to be within tolerance."""
    documented = dict(_budget_rows())
    total = documented["Total"]
    parts = sum(v for k, v in documented.items() if k != "Total")
    step = 10 ** -ROUND_TO
    assert abs(parts - total) <= step * 3, (
        f"the budget rows sum to {parts:,} against a documented Total of "
        f"{total:,}")


# --- the manuscript's claims about ITSELF ---------------------------------

def test_the_manuscript_states_its_own_length_correctly():
    """A document that misdescribes itself in its first line.

    The front matter claimed ~46,500 words while the manuscript measured
    53,400 -- a 15% understatement that had drifted quietly, because the
    word-budget guard checks `book-outline.md` and nothing checked the
    manuscript's claim about itself. It is the first factual statement a
    reader meets.
    """
    import re

    txt = MANUSCRIPT.read_text()
    m = re.search(r"~([\d,]+) words", txt)
    assert m, "the manuscript no longer states its own word count"
    claimed = int(m.group(1).replace(",", ""))
    measured = prose_words(txt)
    step = 10 ** -ROUND_TO
    assert abs(claimed - measured) <= step, (
        f"the manuscript says ~{claimed:,} words and measures {measured:,} "
        f"(rounds to {round(measured, ROUND_TO):,}). Run "
        "scripts/update_word_budget.py, which now updates this line too.")


def test_the_manuscript_counts_its_own_chapters_and_appendices():
    """The same class of self-description, checked the same way.

    Cheap to verify and easy to leave stale through a restructure -- and a
    reader has no way to notice, because the claim sits far from the thing it
    counts.
    """
    import re

    txt = MANUSCRIPT.read_text()
    chapters = len(re.findall(r"^## Chapter ", txt, re.M))
    appendices = len(re.findall(r"^## Appendix ", txt, re.M))
    m = re.search(r"(\d+) chapters \+ (\d+) appendices", txt)
    assert m, "the manuscript no longer states its chapter/appendix counts"
    assert int(m.group(1)) == chapters, (
        f"front matter says {m.group(1)} chapters; the manuscript has "
        f"{chapters}")
    assert int(m.group(2)) == appendices, (
        f"front matter says {m.group(2)} appendices; the manuscript has "
        f"{appendices}")


def test_no_reader_facing_file_quotes_a_stale_word_count():
    """Three files quoted the manuscript's length and all three disagreed.

    The front matter said ~46,500, the README and CLAUDE.md said ~48,500, and
    the manuscript measured 53,400. Each had been updated at a different time
    by someone fixing the one in front of them, which is exactly how three
    sites end up with three different wrong numbers.

    `scripts/update_word_budget.py` owns all of them now. This checks the
    property rather than trusting the helper.
    """
    import re

    measured = prose_words(MANUSCRIPT.read_text())
    step = 10 ** -ROUND_TO
    bad = []
    for name in ("article/drafts/v1.md", "README.md", "CLAUDE.md"):
        path = REPO_ROOT / name
        if not path.exists():
            continue
        for m in re.finditer(r"~([\d,]+) words", path.read_text()):
            claimed = int(m.group(1).replace(",", ""))
            if abs(claimed - measured) > step:
                bad.append(f"{name}: ~{claimed:,} against a measured "
                           f"{measured:,}")
    assert not bad, (
        "files quoting a stale manuscript length:\n  " + "\n  ".join(bad)
        + "\nRun scripts/update_word_budget.py, which updates every site.")


def test_the_outline_numbers_the_same_chapters_the_manuscript_does():
    """The outline is the manuscript's contract; it had silently gone off by one.

    Chapter 6 was inserted into `v1.md` and never into the outline, so from
    that point on every outline heading named the chapter after it, its
    per-part rows described the wrong chapter ranges, and its Total said
    eleven chapters against a manuscript with twelve. The word-budget
    arithmetic was right under the stale numbering, which is why the existing
    guard (rows must sum) stayed green through all of it -- a sum cannot see a
    label.
    """
    outline = (REPO_ROOT / "article/book-outline.md").read_text()
    manuscript = (REPO_ROOT / "article/drafts/v1.md").read_text()

    def chapters(text):
        return [(int(m.group(1)), m.group(2).strip())
                for m in re.finditer(r"^## Chapter (\d+): ([^(\n]+)", text, re.M)]

    out = chapters(outline)
    ms = chapters(manuscript)
    assert [n for n, _ in out] == list(range(1, len(out) + 1)), (
        f"the outline's chapter numbers are not consecutive: {[n for n, _ in out]}")
    assert len(out) == len(ms), (
        f"the outline describes {len(out)} chapters, the manuscript has "
        f"{len(ms)}: {[t for _, t in out]} vs {[t for _, t in ms]}")
    for (a, ta), (b, tb) in zip(out, ms):
        assert a == b and ta.lower() == tb.lower(), (
            f"outline Chapter {a} '{ta}' does not match manuscript "
            f"Chapter {b} '{tb}'")
    total = re.search(r"\| \*\*Total\*\* \| \*\*(\d+) \+ 3 apps\*\*", outline)
    assert total, "the outline's Total row no longer states a chapter count"
    assert int(total.group(1)) == len(ms), (
        f"the outline's Total says {total.group(1)} chapters; the manuscript "
        f"has {len(ms)}")


def test_each_part_target_equals_its_chapters_budgets():
    """The TARGET column was reconciled by nothing.

    Inserting Chapter 6 with a 5,000-word budget left Part III's target at
    16,500 while its chapters summed to 21,500, and the Total short by the same
    5,000 -- with the existing guard green throughout, because it reconciles
    the measured *Current words* column and never the planned one. A budget
    nothing adds up is not a budget.
    """
    text = (REPO_ROOT / "article/book-outline.md").read_text()
    budgets = {int(m.group(1)): int(m.group(2).replace(",", ""))
               for m in re.finditer(
                   r"^## Chapter (\d+): [^(]+\(~([\d,]+) words\)", text, re.M)}
    assert budgets, "no chapter budgets found in the outline"
    rows = re.findall(
        r"^\| (?!\*\*Total)([^|]+?) \| (\d+)-(\d+) \| ~[\d,]+ \| ([\d,]+) \|",
        text, re.M)
    assert rows, "no per-part rows with a chapter range found"
    covered = set()
    for name, lo, hi, target in rows:
        chapters = range(int(lo), int(hi) + 1)
        missing = [c for c in chapters if c not in budgets]
        assert not missing, f"{name.strip()} names chapters with no budget: {missing}"
        got = sum(budgets[c] for c in chapters)
        assert got == int(target.replace(",", "")), (
            f"{name.strip()} targets {target} but chapters {lo}-{hi} budget "
            f"{got:,}")
        covered.update(chapters)
    assert covered == set(budgets), (
        f"chapters in no part row: {sorted(set(budgets) - covered)}")


def test_the_issue_mapping_uses_the_live_chapter_numbers():
    """The mapping kept the pre-insertion numbering when the part table above
    it was updated, so it pointed at the wrong chapters and had no row at all
    for the last one."""
    text = (REPO_ROOT / "article/book-outline.md").read_text()
    n = len(re.findall(r"^## Chapter \d+: ", text, re.M))
    section = text.split("## Chapter → Issue Mapping", 1)
    assert len(section) == 2, "the issue-mapping section is gone"
    covered = set()
    for lo, hi in re.findall(r"^\| (\d+)(?:-(\d+))? \|", section[1], re.M):
        covered.update(range(int(lo), int(hi or lo) + 1))
    assert covered == set(range(1, n + 1)), (
        f"the issue mapping covers {sorted(covered)} against {n} chapters; "
        f"missing {sorted(set(range(1, n + 1)) - covered)}")
