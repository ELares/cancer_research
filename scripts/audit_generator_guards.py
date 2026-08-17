#!/usr/bin/env python3
"""ADVISORY: guards that read a generated artifact but never their generator.

NON-BLOCKING BY DESIGN. The repo already learned this with
`audit_heading_only_assertions.py`: a gate that fires at roughly half precision
trains people to add exemptions rather than fix findings. This prints.

THE PATTERN. A script computes something and writes a committed report; the
guard reads the report. That guard cannot catch a change to the SCRIPT, because
the report is a committed file and does not move when the generator is edited.

Found six times in one day, every time by a mutation sweep rather than by the
suite: a hardcoded winner where a maximum used to be derived, a deleted
check-tag filter, an empty-result refusal replaced with `if False:`, a dropped
caveat, a ranking reverted to a dict that sort_keys reorders, and a central
finding removed from a renderer.

The fix in each case was one assertion against the generator's source alongside
the one against its output.

Usage:
    python scripts/audit_generator_guards.py
"""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"

SOURCE_MARKERS = ("SCRIPT.read_text()", "_src()", "SCRIPT).read_text()",
                  "script.read_text()", "SRC.read_text()")


def names_generator_and_artifact(tree):
    script = artifact = False
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and isinstance(n.value, str):
            v = n.value
            if v.endswith(".py") and "test" not in v:
                script = True
            if v.endswith((".md", ".json", ".tsv.gz")):
                artifact = True
    return script, artifact


def main():
    hits = []
    for p in sorted(TESTS.glob("test_*.py")):
        src = p.read_text(encoding="utf-8")
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        s, a = names_generator_and_artifact(tree)
        if s and a and not any(m in src for m in SOURCE_MARKERS):
            hits.append(p.name)
    print(f"{len(hits)} guard(s) read a generated artifact without pinning "
          f"their generator's source:\n")
    for h in hits:
        print(f"  {h}")
    print("\nADVISORY. Precision is roughly half: many of these are pure data "
          "checks where the artifact IS the subject and no generator assertion "
          "is owed. Read before acting.")


if __name__ == "__main__":
    main()
