#!/usr/bin/env python3
"""Convert inline academic citations to Markdown footnote syntax.

Transforms:
  [PMID: XXXXX]           → [^pmidXXXXX]  + definition from reference list
  [PMID: A; PMID: B; ...] → [^groupN]     + combined definition
  [News: Author, "Title", URL, Date. Status] → [^newsN] + definition
  [Commentary: ...]        → [^commentaryN] + definition

Does NOT handle author-year (Porter et al., ...) or textbook (Biology2e Ch.7)
citations — those require manual prose editing.

Usage:
    python scripts/convert_citations_to_footnotes.py          # dry run (prints changes)
    python scripts/convert_citations_to_footnotes.py --apply  # writes to v1.md
"""

import argparse
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
V1_PATH = PROJECT_ROOT / "article" / "drafts" / "v1.md"


def parse_reference_list(md: str) -> dict[str, str]:
    """Build PMID → full citation text from the numbered reference list."""
    refs = {}
    for m in re.finditer(
        r"^\d+\.\s+PMID:\s+(\d+)\s+--\s+(.+?)$", md, re.MULTILINE
    ):
        pmid = m.group(1)
        # Clean up the reference text: remove citation count at end
        text = m.group(2).strip()
        text = re.sub(r"\.\s+\d+\s+citations?\.\s*$", ".", text)
        # Convert *Journal* to Journal (footnotes don't need markdown italics)
        text = text.replace("*", "")
        refs[pmid] = text
    return refs


def convert_single_pmid(md: str, refs: dict[str, str]) -> tuple[str, dict[str, str]]:
    """Convert [PMID: XXXXX] → [^pmidXXXXX] and collect definitions."""
    definitions = {}

    def repl(m):
        pmid = m.group(1)
        label = f"pmid{pmid}"
        if pmid in refs:
            definitions[label] = f"{refs[pmid]} PMID: {pmid}."
        else:
            definitions[label] = f"PMID: {pmid}."
        return f"[^{label}]"

    # Match single-PMID brackets: [PMID: 12345] or [PMID:12345]
    md = re.sub(r"\[PMID:\s*(\d+)\]", repl, md)
    return md, definitions


def convert_multi_pmid(md: str, refs: dict[str, str]) -> tuple[str, dict[str, str]]:
    """Convert [PMID: A; PMID: B; ...] → [^groupN] with combined definition."""
    definitions = {}
    counter = [0]

    def repl(m):
        bracket_content = m.group(1)
        pmids = re.findall(r"PMID:\s*(\d+)", bracket_content)
        if len(pmids) <= 1:
            return m.group(0)  # not a multi-PMID, skip

        counter[0] += 1
        label = f"refs_group{counter[0]}"
        parts = []
        for pmid in pmids:
            if pmid in refs:
                parts.append(f"{refs[pmid]} PMID: {pmid}")
            else:
                parts.append(f"PMID: {pmid}")
        definitions[label] = "; ".join(parts) + "."
        return f"[^{label}]"

    # Match brackets containing multiple PMIDs separated by semicolons
    md = re.sub(r"\[(PMID:\s*\d+(?:\s*;\s*PMID:\s*\d+)+)\]", repl, md)
    return md, definitions


def convert_news(md: str) -> tuple[str, dict[str, str]]:
    """Convert [News: ...] and [Commentary: ...] → [^newsN]."""
    definitions = {}
    counter = [0]

    def repl(m):
        counter[0] += 1
        tag = m.group(1).lower().replace(" ", "")  # "News" or "Commentary"
        content = m.group(2).strip()
        label = f"{tag}{counter[0]}"
        # Clean up the content for footnote
        # Remove structured format, keep readable
        definitions[label] = content.rstrip(".")  + "."
        return f"[^{label}]"

    # Match [News: ...] and [Commentary: ...]
    md = re.sub(
        r"\[(News|Commentary):\s*(.*?)\]",
        repl,
        md,
        flags=re.DOTALL,
    )
    return md, definitions


def insert_definitions_at_section_ends(md: str, all_defs: dict[str, str]) -> str:
    """Insert footnote definitions at the end of each ### section."""
    lines = md.split("\n")
    result = []
    pending_defs = {}

    # Track which definitions are referenced in the current section
    current_section_refs = set()

    for i, line in enumerate(lines):
        # Check if this line starts a new section (### or ## or #)
        is_section_boundary = (
            re.match(r"^#{1,4}\s", line)
            and i > 0  # don't trigger on the very first line
        )

        if is_section_boundary and current_section_refs:
            # Insert pending definitions before this section header
            result.append("")
            for label in sorted(current_section_refs):
                if label in all_defs:
                    result.append(f"[^{label}]: {all_defs[label]}")
                    del all_defs[label]
            result.append("")
            current_section_refs = set()

        # Track footnote references in this line
        for m in re.finditer(r"\[\^(\w+)\](?!:)", line):
            label = m.group(1)
            if label in all_defs:
                current_section_refs.add(label)

        result.append(line)

    # Insert any remaining definitions at the end
    if current_section_refs or all_defs:
        result.append("")
        remaining = current_section_refs | set(all_defs.keys())
        for label in sorted(remaining):
            if label in all_defs:
                result.append(f"[^{label}]: {all_defs[label]}")
        result.append("")

    return "\n".join(result)


def main():
    parser = argparse.ArgumentParser(description="Convert citations to footnotes.")
    parser.add_argument("--apply", action="store_true", help="Write changes to v1.md")
    args = parser.parse_args()

    md = V1_PATH.read_text(encoding="utf-8")

    # Parse reference list
    refs = parse_reference_list(md)
    print(f"Parsed {len(refs)} PMID references from reference list")

    # Convert multi-PMID brackets first (before singles, to avoid partial matches)
    md, multi_defs = convert_multi_pmid(md, refs)
    print(f"Converted {len(multi_defs)} multi-PMID brackets")

    # Convert single PMIDs
    md, single_defs = convert_single_pmid(md, refs)
    print(f"Converted {len(single_defs)} single PMID citations")

    # Convert News/Commentary
    md, news_defs = convert_news(md)
    print(f"Converted {len(news_defs)} News/Commentary citations")

    # Merge all definitions
    all_defs = {**multi_defs, **single_defs, **news_defs}
    print(f"Total footnote definitions: {len(all_defs)}")

    # Insert definitions at section ends
    md = insert_definitions_at_section_ends(md, all_defs)

    # Report remaining inline citations that need manual conversion
    remaining_pmid = len(re.findall(r"\[PMID:", md))
    remaining_news = len(re.findall(r"\[News:", md))
    remaining_commentary = len(re.findall(r"\[Commentary:", md))
    print(f"\nRemaining (need manual conversion):")
    print(f"  [PMID: ...]: {remaining_pmid}")
    print(f"  [News: ...]: {remaining_news}")
    print(f"  [Commentary: ...]: {remaining_commentary}")

    if args.apply:
        V1_PATH.write_text(md, encoding="utf-8")
        print(f"\nWritten to {V1_PATH}")
    else:
        print(f"\nDry run — use --apply to write changes")


if __name__ == "__main__":
    main()
