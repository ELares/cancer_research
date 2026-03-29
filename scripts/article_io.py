"""Shared article I/O utilities for reading/writing corpus markdown files."""

import re
from pathlib import Path

import yaml


def load_article(filepath: Path) -> tuple[dict, str]:
    """Load article, return (frontmatter_dict, body_text).

    Returns ({}, full_content) if frontmatter can't be parsed.
    """
    content = filepath.read_text(encoding="utf-8")
    match = re.match(r"^---\n(.*?\n)---\n\n?(.*)", content, re.DOTALL)
    if not match:
        return {}, content

    fm = yaml.safe_load(match.group(1)) or {}
    body = match.group(2)
    return fm, body


def save_article(filepath: Path, frontmatter: dict, body: str) -> None:
    """Save article with updated YAML frontmatter + markdown body."""
    yaml_str = yaml.dump(
        frontmatter,
        default_flow_style=False,
        allow_unicode=True,
        sort_keys=False,
    )
    content = f"---\n{yaml_str}---\n\n{body}"
    filepath.write_text(content, encoding="utf-8")


def load_frontmatter(filepath: Path) -> dict | None:
    """Load only the YAML frontmatter from an article file.

    More efficient than load_article when body content is not needed.
    """
    with open(filepath, "r", encoding="utf-8") as f:
        # Read until we find the closing ---
        first_line = f.readline()
        if first_line.strip() != "---":
            return None

        yaml_lines = []
        for line in f:
            if line.strip() == "---":
                break
            yaml_lines.append(line)

    if not yaml_lines:
        return None

    return yaml.safe_load("".join(yaml_lines))
