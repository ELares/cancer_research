"""Every chapter and section the manuscript refers to must exist.

WHY THIS EXISTS, and it is not hypothetical. Inserting Chapter 6 (the
multi-modality engine) renumbered six chapters and twenty-seven subsections,
and every inline "Chapter 7" written before that insertion silently began
pointing at a different chapter. Nothing in this repository would have caught
it: the prose guards check what sentences SAY, and a cross-reference is
correct-looking text that resolves to the wrong place.

Two failure modes, and only the first is obvious:

* a reference to a chapter or section that does not exist -- loud once looked
  for, invisible otherwise;
* a reference that resolves to a real heading which is not the one meant.
  Nothing mechanical can detect that in general, so this file pins the
  handful of references whose TARGET is nameable, by requiring the heading
  they point at to contain the words the sentence promises.
"""
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "article/drafts/v1.md"


def _text() -> str:
    return MD.read_text()


def _chapters(s: str) -> dict:
    return {int(m.group(1)): m.group(2)
            for m in re.finditer(r"^## Chapter (\d+): (.+)$", s, re.M)}


def _sections(s: str) -> dict:
    return {m.group(1): m.group(2)
            for m in re.finditer(r"^### (\d+\.\d+) (.+)$", s, re.M)}


def _inline(pattern: str, s: str):
    """Matches that are NOT themselves headings."""
    for m in re.finditer(pattern, s):
        line_start = s.rfind("\n", 0, m.start()) + 1
        if s[line_start:line_start + 3] in ("## ", "###"):
            continue
        yield m


def test_every_chapter_reference_resolves():
    s = _text()
    chapters = _chapters(s)
    assert len(chapters) >= 10, chapters
    missing = {}
    for m in _inline(r"Chapters?\s+(\d+)", s):
        n = int(m.group(1))
        if n not in chapters:
            missing.setdefault(n, 0)
            missing[n] += 1
    assert not missing, (
        f"the manuscript refers to chapters that do not exist: {missing}. "
        f"Chapters present: {sorted(chapters)}")


def test_every_section_reference_resolves():
    s = _text()
    sections = _sections(s)
    assert len(sections) >= 20, len(sections)
    missing = {}
    for m in _inline(r"(?:Sections?|§)\s*(\d+\.\d+)", s):
        ref = m.group(1)
        if ref not in sections:
            missing.setdefault(ref, 0)
            missing[ref] += 1
    assert not missing, (
        f"the manuscript refers to sections that do not exist: {missing}")


def test_every_subsection_number_matches_its_own_chapter():
    """A renumbered chapter leaves its subsections behind.

    This is what actually happened: Chapter 7's subsections still read 6.1,
    6.2, 6.3 after the insertion, so every reference to "Section 6.2" pointed
    into a chapter that no longer contained it.
    """
    s = _text()
    current = None
    wrong = []
    for line in s.split("\n"):
        m = re.match(r"^## Chapter (\d+):", line)
        if m:
            current = int(m.group(1))
            continue
        m2 = re.match(r"^### (\d+)\.(\d+) ", line)
        if m2 and current is not None and int(m2.group(1)) != current:
            wrong.append((line[:50], current))
    assert not wrong, (
        f"{len(wrong)} subsections carry a chapter number that is not their "
        f"own: {wrong[:5]}")


# References whose TARGET is nameable. A resolving number is not enough --
# these check the heading it lands on is the one the sentence promises.
ANCHORED = [
    ("Chapter 5", "Ferroptosis Engine"),
    ("Chapter 6", "Multi-Modality Engine"),
]


@pytest.mark.parametrize("ref,expected", ANCHORED)
def test_the_anchored_references_point_where_they_claim(ref, expected):
    s = _text()
    n = int(ref.split()[-1])
    chapters = _chapters(s)
    assert n in chapters, f"{ref} does not exist"
    assert expected.lower() in chapters[n].lower(), (
        f"{ref} is '{chapters[n]}', not the {expected!r} the manuscript's "
        "cross-references assume. A renumbering has moved it and every "
        f"sentence saying '{ref}' now points somewhere else.")


def test_the_new_chapter_is_where_the_criticism_is_answered():
    """The chapter exists to answer one criticism, so it has to still contain
    the measurement that answers it rather than becoming a feature list."""
    s = _text()
    start = s.index("## Chapter 6: The Multi-Modality Engine")
    end = s.index("## Chapter 7:")
    body = s[start:end]
    for frag in ("90,019 census articles",
                 "That column is now empty",
                 "Breadth is not depth",
                 "used-in-any-reported-number status"):
        assert frag in body, f"Chapter 6 no longer says: {frag}"
    # And it must keep its refusals, which are what stop it being a feature
    # list. Each names a limit that would change a claim if lifted.
    for refusal in ("is not a ranking of therapies",
                    "are one function and a configuration struct",
                    "The taxonomy bounds everything above"):
        assert refusal in body, f"Chapter 6 dropped the refusal: {refusal}"
