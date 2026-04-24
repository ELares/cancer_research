"""Generate a SHA256 manifest and optional filtered archive for DOI release.

Walks all git-tracked files, excludes non-redistributable content
(corpus full text, LFS book pointers), computes SHA256 checksums,
and optionally builds a filtered tar.gz archive suitable for Zenodo.

Usage:
    python scripts/generate_release_manifest.py               # manifest only
    python scripts/generate_release_manifest.py --build-archive  # manifest + archive
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Prefixes excluded from the release scope.
# corpus/by-pmid/  — 61 of 4,830 articles are non-OA; redistribution unsafe
# books/           — LFS pointers only; actual PDFs are 1.6 GB reference material
# news/            — fetched articles under fair use; not a redistributable license
EXCLUDE_PREFIXES = ("corpus/by-pmid/", "books/", "news/")

MANIFEST_FILE = PROJECT_ROOT / "MANIFEST.sha256"
ARCHIVE_FILE = PROJECT_ROOT / "release.tar.gz"
INDEX_FILE = PROJECT_ROOT / "corpus" / "INDEX.jsonl"


def get_tracked_files():
    """Return list of git-tracked file paths (relative to repo root)."""
    result = subprocess.run(
        ["git", "ls-files"],
        capture_output=True, text=True, cwd=PROJECT_ROOT,
    )
    if result.returncode != 0:
        print(f"Error running git ls-files: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return [f for f in result.stdout.strip().split("\n") if f]


def filter_files(files):
    """Exclude files outside the release scope."""
    return [f for f in files if not any(f.startswith(p) for p in EXCLUDE_PREFIXES)]


def compute_sha256(filepath):
    """Compute SHA256 hash of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def warn_non_oa():
    """Print non-OA article PMIDs as an informational warning."""
    if not INDEX_FILE.exists():
        print("  WARNING: INDEX.jsonl not found — cannot check OA status")
        return
    non_oa = []
    with open(INDEX_FILE) as f:
        for line in f:
            entry = json.loads(line)
            if not entry.get("is_oa", True):
                non_oa.append(entry["pmid"])
    print(f"\n  Non-OA articles excluded from release: {len(non_oa)}")
    if non_oa:
        print(f"  PMIDs: {', '.join(non_oa[:10])}{'...' if len(non_oa) > 10 else ''}")
    print(f"  These {len(non_oa)} full-text files in corpus/by-pmid/ are NOT in the archive.")
    print("  INDEX.jsonl (metadata) IS included for corpus reconstruction.\n")


def main():
    parser = argparse.ArgumentParser(description="Generate release manifest and optional archive.")
    parser.add_argument("--build-archive", action="store_true",
                        help="Also build a filtered release.tar.gz")
    args = parser.parse_args()

    print("Generating release manifest...")
    all_files = get_tracked_files()
    in_scope = filter_files(all_files)
    excluded_count = len(all_files) - len(in_scope)

    print(f"  Total tracked files: {len(all_files)}")
    print(f"  Excluded (corpus/by-pmid/ + books/ + news/): {excluded_count}")
    print(f"  In release scope: {len(in_scope)}")

    # Compute checksums and write manifest
    total_size = 0
    with open(MANIFEST_FILE, "w") as mf:
        for filepath in sorted(in_scope):
            full_path = PROJECT_ROOT / filepath
            if not full_path.is_file():
                continue
            sha = compute_sha256(full_path)
            total_size += full_path.stat().st_size
            mf.write(f"{sha}  {filepath}\n")

    size_mb = total_size / (1024 * 1024)
    print(f"  Manifest written: {MANIFEST_FILE.name} ({size_mb:.1f} MB across {len(in_scope)} files)")

    warn_non_oa()

    if args.build_archive:
        print(f"Building filtered archive: {ARCHIVE_FILE.name}...")
        with tarfile.open(ARCHIVE_FILE, "w:gz") as tar:
            for filepath in sorted(in_scope):
                full_path = PROJECT_ROOT / filepath
                if full_path.is_file():
                    tar.add(full_path, arcname=filepath)
        archive_mb = ARCHIVE_FILE.stat().st_size / (1024 * 1024)
        print(f"  Archive written: {ARCHIVE_FILE.name} ({archive_mb:.1f} MB)")
        print(f"\n  Verify with: tar tzf {ARCHIVE_FILE.name} | grep corpus/by-pmid/ | wc -l  # should be 0")

    print("\nDone.")


if __name__ == "__main__":
    main()
