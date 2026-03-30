#!/usr/bin/env python3
"""Move abstract-only corpus records into a separate directory.

Records are identified by the placeholder:
`Full text not downloaded (open access|paywalled).`
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

from config import ABSTRACT_PMID_DIR, PMID_DIR


PLACEHOLDER_RE = re.compile(r"Full text not downloaded \((open access|paywalled)\)\.")


def main() -> None:
    ABSTRACT_PMID_DIR.mkdir(parents=True, exist_ok=True)

    moved = 0
    for path in sorted(PMID_DIR.glob("*.md")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        if not PLACEHOLDER_RE.search(text):
            continue

        target = ABSTRACT_PMID_DIR / path.name
        if target.exists():
            target.unlink()
        shutil.move(str(path), str(target))
        moved += 1

    print(f"Moved {moved} abstract-only records to {ABSTRACT_PMID_DIR}")


if __name__ == "__main__":
    main()
