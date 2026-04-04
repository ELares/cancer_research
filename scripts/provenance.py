"""Helpers for lightweight pipeline provenance logging."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from config import INDEX_FILE, PMID_DIR, PROJECT_ROOT

ANALYSIS_DIR = PROJECT_ROOT / "analysis"
PROVENANCE_LOG = ANALYSIS_DIR / "provenance.jsonl"
QUERY_FILE = PROJECT_ROOT / "scripts" / "queries.txt"
REQUIREMENTS_FILE = PROJECT_ROOT / "requirements.txt"


def _parse_requirements_packages() -> list[str]:
    """Extract package names from requirements.txt, falling back to a static list."""
    if not REQUIREMENTS_FILE.exists():
        return ["PyYAML", "python-dotenv", "requests", "tqdm", "matplotlib", "numpy", "scipy"]
    packages = []
    for line in REQUIREMENTS_FILE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Strip version specifiers (>=, ==, ~=, etc.)
        name = line.split(">=")[0].split("==")[0].split("~=")[0].split("<=")[0].split("!=")[0].split("[")[0].strip()
        if name:
            packages.append(name)
    return packages


def _safe_package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for package in _parse_requirements_packages():
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _git_output(args: list[str]) -> str:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""
    return result.stdout.strip()


def _sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fingerprint_article_corpus() -> dict[str, str | int]:
    files = sorted(PMID_DIR.glob("*.md"))
    digest = hashlib.sha256()
    for path in files:
        stat = path.stat()
        digest.update(str(path.relative_to(PROJECT_ROOT)).encode("utf-8"))
        digest.update(str(stat.st_size).encode("utf-8"))
        digest.update(str(stat.st_mtime_ns).encode("utf-8"))
    return {
        "article_count": len(files),
        "article_metadata_fingerprint": digest.hexdigest(),
    }


def append_provenance_record(script_name: str, extra: dict | None = None) -> None:
    record: dict[str, object] = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "script": script_name,
        "git_commit": _git_output(["rev-parse", "HEAD"]),
        "git_dirty": bool(_git_output(["status", "--porcelain"])),
        "python_version": platform.python_version(),
        "package_versions": _safe_package_versions(),
        "query_file_hash": _sha256_file(QUERY_FILE),
        "index_hash": _sha256_file(INDEX_FILE),
    }
    record.update(_fingerprint_article_corpus())
    if extra:
        record.update(extra)

    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    with PROVENANCE_LOG.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")
