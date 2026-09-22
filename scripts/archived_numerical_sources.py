"""Resolve exact historical numerical inputs without relabeling old archives.

The present checkout is preferred when its bytes match. Explicitly registered,
content-addressed source snapshots support older archives after code evolves.
This verifies original input identity; it does not claim that today's report
reader is the original numerical program or execute a historical simulation.
"""

import hashlib
import json
from pathlib import Path, PurePosixPath
import re

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOTS = Path("analysis/calibration/frozen-sources")


def resolve_source_bytes(relative: str, digest: str, *, root: Path = ROOT) -> bytes:
    """Return bytes only for the requested path and its exact archived SHA-256."""
    if (not isinstance(relative, str) or not relative
            or PurePosixPath(relative).is_absolute()
            or any(part in (".", "..") for part in relative.split("/"))
            or "\\" in relative or not isinstance(digest, str)
            or re.fullmatch(r"[0-9a-f]{64}", digest) is None):
        raise ValueError("invalid archived source identity")
    current = root / relative
    if current.resolve().is_relative_to(root.resolve()) and current.is_file():
        blob = current.read_bytes()
        if hashlib.sha256(blob).hexdigest() == digest:
            return blob
    manifest_path = root / SNAPSHOTS / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"no registered historical source: {relative} at {digest}")
    manifest = json.loads(manifest_path.read_text())
    if (not isinstance(manifest, dict) or type(manifest.get("schema_version")) is not int
            or manifest["schema_version"] != 1 or not isinstance(manifest.get("sources"), dict)):
        raise ValueError("invalid historical source registry schema")
    registered = manifest["sources"].get(relative)
    if not isinstance(registered, dict) or digest not in registered:
        raise ValueError(f"no registered historical source: {relative} at {digest}")
    record = registered[digest]
    if (not isinstance(record, dict) or not isinstance(record.get("recovered_from_commit"), str)
            or re.fullmatch(r"[0-9a-f]{40}", record["recovered_from_commit"]) is None):
        raise ValueError(f"invalid historical source registration: {relative} at {digest}")
    snapshot = root / SNAPSHOTS / f"{digest}.source"
    if (not snapshot.is_file() or snapshot.is_symlink()
            or not snapshot.resolve().is_relative_to(root.resolve())):
        raise ValueError(f"historical source snapshot unavailable: {relative} at {digest}")
    blob = snapshot.read_bytes()
    if hashlib.sha256(blob).hexdigest() != digest:
        raise ValueError(f"historical source snapshot hash mismatch: {relative} at {digest}")
    return blob


def verify_source_hashes(hashes: dict[str, str], *, root: Path = ROOT) -> None:
    """Verify every declared original source against actual available bytes."""
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError("nonempty archived source inventory required")
    for relative, digest in hashes.items():
        resolve_source_bytes(relative, digest, root=root)
