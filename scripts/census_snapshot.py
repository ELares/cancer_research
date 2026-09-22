"""Verify and fingerprint a complete indexed baseline before preparing review.

The manifest is an acquisition record, not evidence that its shards are still
readable or contain the declared counts. Each compressed shard is read once
into memory, hashed, and parsed from those same bytes. Only selected records
are retained; memory use is bounded by the largest compressed shard plus the
selected intersections. Update, unindexed and C04-only streams are not merged.
"""
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import stat
from typing import Callable
import zlib


def _canonical_sha256(value) -> str:
    encoded = json.dumps(value, sort_keys=True, ensure_ascii=False,
                         separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _signature(path: Path) -> tuple:
    value = path.stat()
    if not stat.S_ISREG(value.st_mode):
        raise ValueError(f"not a regular shard: {path.name}")
    return (value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)


def _nonfinite(value):
    raise ValueError(f"non-finite JSON number: {value}")


def validate_snapshot(snapshot: dict) -> None:
    """Check archived inventory/count consistency without rereading any source.

    A consistent digest is not authentication, and the manifest fingerprints
    cannot be recomputed without the original manifest. This is only an
    internal-consistency check on the archived snapshot description.
    """
    def digest(value) -> bool:
        return (isinstance(value, str) and len(value) == 64
                and all(character in "0123456789abcdef" for character in value))

    def count(value) -> bool:
        return type(value) is int and value >= 0

    try:
        if (type(snapshot["schema_version"]) is not int or snapshot["schema_version"] != 1
                or snapshot["stream"] != "records" or type(snapshot["stride"]) is not int
                or snapshot["stride"] != 1):
            raise ValueError("unsupported snapshot schema, stream or sampling stride")
        if not isinstance(snapshot["source_url"], str) or not snapshot["source_url"].strip():
            raise ValueError("missing snapshot source URL")
        for key in ("manifest_sha256", "baseline_manifest_sha256", "inventory_sha256"):
            if not digest(snapshot[key]):
                raise ValueError(f"invalid {key}")
        shards = snapshot["shards"]
        if not isinstance(shards, list) or not shards:
            raise ValueError("empty snapshot inventory")
        names = []
        for shard in shards:
            if set(shard) != {"name", "sha256", "bytes", "records"}:
                raise ValueError("invalid snapshot shard fields")
            name = shard["name"]
            if (not isinstance(name, str) or not name.endswith(".jsonl.gz")
                    or Path(name).name != name or "\\" in name):
                raise ValueError("invalid snapshot shard name")
            if not digest(shard["sha256"]) or not count(shard["bytes"]) or not shard["bytes"]:
                raise ValueError("invalid snapshot shard fingerprint or byte count")
            if not count(shard["records"]):
                raise ValueError("invalid snapshot shard record count")
            names.append(name)
        if names != sorted(set(names)):
            raise ValueError("snapshot shards must be unique and sorted")
        if snapshot["inventory_sha256"] != _canonical_sha256(shards):
            raise ValueError("snapshot inventory fingerprint does not match its entries")
        observed, declared = snapshot["observed_records"], snapshot["declared_records"]
        if (not count(observed) or not count(declared) or not observed
                or observed != declared or observed != sum(row["records"] for row in shards)):
            raise ValueError("snapshot record counts do not reconcile")
        basis = snapshot["cancer_basis"]
        if (not isinstance(basis, dict)
                or any(not isinstance(key, str) or not key.strip() or not count(value)
                       for key, value in basis.items())
                or sum(basis.values()) != observed):
            raise ValueError("snapshot cancer-basis counts do not reconcile")
    except (KeyError, TypeError, ValueError) as exc:
        raise SystemExit(f"Invalid census snapshot metadata: {exc}") from exc


def read_snapshot(root: Path, select: Callable[[dict], bool]) -> tuple[dict, list[dict]]:
    """Exhaust the indexed baseline, verify counts, and retain selected records.

    Fingerprints describe the bytes actually read, not a downloaded source's
    authenticity. The caller chooses the scientific intersection. This function
    neither adjudicates records nor establishes global PMID uniqueness.
    """
    root = Path(root)
    try:
        return _read_snapshot(root, select)
    except (OSError, UnicodeError, ValueError, EOFError, zlib.error) as exc:
        raise SystemExit(f"Invalid census snapshot at {root}: {exc}") from exc


def _read_snapshot(root: Path, select: Callable[[dict], bool]) -> tuple[dict, list[dict]]:
    manifest_path = root / "manifest.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes, parse_constant=_nonfinite)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("files"), dict):
        raise ValueError("manifest needs a files mapping")
    source = manifest.get("source")
    if not isinstance(source, str) or not source.strip():
        raise ValueError("manifest needs its source URL")

    baseline = {}
    expected = {}
    for name, entry in manifest["files"].items():
        if not isinstance(entry, dict):
            raise ValueError(f"invalid manifest entry: {name}")
        stream = entry.get("source", "baseline")
        if stream == "updatefiles":
            continue
        if stream != "baseline":
            raise ValueError(f"unknown manifest source for {name}: {stream}")
        count, shard_name = entry.get("cancer"), entry.get("records")
        if entry.get("parsed") is not True or type(count) is not int or count < 0:
            raise ValueError(f"baseline entry {name} needs parsed=true and a nonnegative count")
        if (not isinstance(shard_name, str) or not shard_name.endswith(".jsonl.gz")
                or Path(shard_name).name != shard_name or "\\" in shard_name):
            raise ValueError(f"unsafe or invalid shard name in {name}")
        if shard_name in expected:
            raise ValueError(f"duplicate manifest shard: {shard_name}")
        baseline[name] = entry
        expected[shard_name] = count
    if not expected:
        raise ValueError("no indexed baseline shards in the manifest")

    records = root / "records"

    def inventory() -> list[str]:
        return sorted(path.name for path in records.glob("*.jsonl.gz"))

    names = sorted(expected)
    if inventory() != names:
        raise ValueError("indexed shard inventory does not match the baseline manifest")
    signatures = {name: _signature(records / name) for name in names}
    if len({signature[:2] for signature in signatures.values()}) != len(names):
        raise ValueError("indexed shard paths alias the same file")

    retained = []
    shards = []
    basis_counts: Counter = Counter()
    for name in names:
        path = records / name
        if _signature(path) != signatures[name]:
            raise ValueError(f"shard changed during snapshot read: {name}")
        compressed = path.read_bytes()
        if not compressed.startswith(b"\x1f\x8b"):
            raise ValueError(f"invalid gzip header in {name}")
        count = 0
        with gzip.GzipFile(fileobj=io.BytesIO(compressed), mode="rb") as stream:
            with io.TextIOWrapper(stream, encoding="utf-8") as text:
                for line in text:
                    record = json.loads(line, parse_constant=_nonfinite)
                    if not isinstance(record, dict):
                        raise ValueError(f"non-object record in {name}")
                    basis = record.get("cancer_basis")
                    if basis is None:
                        basis = "unspecified"
                    if not isinstance(basis, str):
                        raise ValueError(f"non-text cancer_basis in {name}")
                    basis_counts[basis if basis.strip() else "unspecified"] += 1
                    count += 1
                    if select(record):
                        retained.append(record)
        if count != expected[name]:
            raise ValueError(f"record count mismatch in {name}: observed {count}, "
                             f"manifest declares {expected[name]}")
        shards.append({"name": name, "sha256": hashlib.sha256(compressed).hexdigest(),
                       "bytes": len(compressed), "records": count})

    observed = sum(shard["records"] for shard in shards)
    if not observed:
        raise ValueError("the indexed baseline contains no records")
    if inventory() != names or manifest_path.read_bytes() != manifest_bytes:
        raise ValueError("manifest or shard inventory changed during snapshot read")
    if any(_signature(records / name) != signatures[name] for name in names):
        raise ValueError("shard changed during snapshot read")
    snapshot = {
        "schema_version": 1, "stream": "records", "stride": 1,
        "source_url": source,
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "baseline_manifest_sha256": _canonical_sha256(baseline),
        "declared_records": sum(expected.values()), "observed_records": observed,
        "cancer_basis": dict(sorted(basis_counts.items())),
        "shards": shards, "inventory_sha256": _canonical_sha256(shards),
    }
    validate_snapshot(snapshot)
    return snapshot, retained
