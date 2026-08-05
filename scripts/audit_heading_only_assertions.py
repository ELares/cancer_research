#!/usr/bin/env python3
"""Find `assert "<literal>" in <doc>` where the literal lives ONLY in a heading.

ADVISORY, NOT A GATE, and the distinction is the point. Measured precision on
this repo is about half: of 13 hits, roughly six are genuine and the rest are
legitimate -- structural checks that a generated report HAS its sections, and
tests that run a string function over synthetic input which merely happens to
contain a heading-shaped literal. Shipping that as a blocking test would train
people to add exemptions, which is how a gate becomes decoration.

So it prints and exits 0. Run it when touching guards over generated prose:
    python scripts/audit_heading_only_assertions.py

WHY THIS SHAPE AND NOT ANOTHER. It is the one that actually recurred here: a
guard names a section's subject but is satisfied by the section's TITLE, which
the generator emits unconditionally. The body can be empty, or say the opposite,
and the assertion nods along. Found live in this repo -- `assert "did NOT
support" in md` matched only the heading of a section whose entire content sits
behind an `if`, so an empty section passed that assertion (the test as a whole
survived on its siblings, which is luck rather than design).

The sibling gate `tests/test_no_vacuous_assertions.py` covers a NARROWER class
that a machine can decide cleanly, and found zero. This one finds real defects
and cannot be made precise, which is why one blocks and the other advises.

Decidable: parse each test, collect `assert <str> in <name>` comparisons, then
search every committed analysis/*.md and article/*.md for that literal. If every
occurrence sits on a markdown heading line (or a bold-label line, which is the
same thing in this repo's generated prose), the assertion cannot distinguish a
populated section from an empty one.
"""
import ast, pathlib, re

DOCS = list(pathlib.Path("analysis").rglob("*.md")) + \
       list(pathlib.Path("article").rglob("*.md"))
TEXT = {p: p.read_text(errors="replace") for p in DOCS}

def is_heading(line: str) -> bool:
    s = line.strip()
    return s.startswith("#") or (s.startswith("**") and s.count("**") >= 2
                                 and len(s) < 200 and s.rstrip().endswith("**"))

def occurrences(lit: str):
    """(doc, n_total, n_heading_only) across the committed documents."""
    for p, t in TEXT.items():
        lines = [ln for ln in t.split("\n") if lit in ln]
        if lines:
            yield p, len(lines), sum(1 for ln in lines if is_heading(ln))

hits = []
for f in sorted(pathlib.Path("tests").glob("test_*.py")):
    try:
        tree = ast.parse(f.read_text())
    except SyntaxError:
        continue
    for fn in ast.walk(tree):
        if not isinstance(fn, ast.FunctionDef) or not fn.name.startswith("test"):
            continue
        for node in ast.walk(fn):
            if not isinstance(node, ast.Assert):
                continue
            for cmp_ in [n for n in ast.walk(node.test) if isinstance(n, ast.Compare)]:
                if not (len(cmp_.ops) == 1 and isinstance(cmp_.ops[0], ast.In)):
                    continue
                left = cmp_.left
                if not (isinstance(left, ast.Constant) and isinstance(left.value, str)):
                    continue
                lit = left.value
                if len(lit) < 6:
                    continue
                found = list(occurrences(lit))
                if not found:
                    continue
                if all(nh == n for _, n, nh in found):
                    where = ", ".join(f"{p.name}({n})" for p, n, _ in found)
                    hits.append((f.name, node.lineno, fn.name, lit, where))

for name, ln, fn, lit, where in hits:
    print(f"{name}:{ln}  {fn}")
    print(f"    assert {lit!r} in <doc>   -- heading-only in {where}")
print(f"\nTOTAL heading-only assertions: {len(hits)}")
