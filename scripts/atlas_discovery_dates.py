"""Read discovery publication dates with a cache bound to local source metadata.

The four available census streams are merged using the earliest dated year for
each PMID. Missing streams are allowed, and absent or null years remain
undated. A present year must be a positive integer; zero is not an undated
sentinel. Dated records require a nonempty string PMID, without normalization.

Cache freshness uses relative paths and filesystem identity, size, modification
time and change time. This avoids rereading gigabytes of compressed records on
every hit. The fingerprint is metadata, not a content digest or proof that the
corpus is complete. Filesystem metadata guarantees and timestamp resolution
limit this check. Checks detect observed changes while loading or rebuilding;
they cannot prevent a producer changing files after the final check.

Publication-year support is not an observation-window end. No acquisition date,
follow-up horizon, or completeness claim is inferred from this year map.
"""

import gzip
import hashlib
import io
import json
import os
import pickle
from pathlib import Path
import stat
import tempfile
import zlib


YEAR_STREAMS = ("records", "records_c04only", "records_unindexed", "records_updates")
CACHE_VERSION = 1
FINGERPRINT_KIND = "filesystem-metadata-v1"


def _metadata(relative_path: str, info) -> dict:
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"date source is not a regular file: {relative_path}")
    return {"path": relative_path, "size": info.st_size,
            "dev": info.st_dev, "ino": info.st_ino,
            "mtime_ns": info.st_mtime_ns, "ctime_ns": info.st_ctime_ns}


def _inventory(root: Path) -> list[dict]:
    """Capture every available shard without reading its compressed contents."""
    entries = []
    for stream in YEAR_STREAMS:
        directory = root / stream
        try:
            info = directory.stat()
        except FileNotFoundError:
            continue
        if not stat.S_ISDIR(info.st_mode):
            raise ValueError(f"date-source stream is not a directory: {stream}")
        # Unlike glob, explicit enumeration does not turn permission errors
        # into an apparently absent source population.
        for path in sorted(directory.iterdir()):
            if path.name.startswith(".") or not path.name.endswith(".jsonl.gz"):
                continue
            relative = path.relative_to(root).as_posix()
            entries.append(_metadata(relative, path.stat()))
    return sorted(entries, key=lambda entry: entry["path"])


def _unchanged(root: Path, inventory: list[dict]) -> None:
    if _inventory(root) != inventory:
        raise ValueError("discovery date-source inventory changed during the read; retry after ingestion finishes")


def _provenance(inventory: list[dict]) -> dict:
    encoded = json.dumps(inventory, sort_keys=True, separators=(",", ":"),
                         allow_nan=False).encode("utf-8")
    counts = dict.fromkeys(YEAR_STREAMS, 0)
    for entry in inventory:
        counts[entry["path"].split("/", 1)[0]] += 1
    return {"inventory_sha256": hashlib.sha256(encoded).hexdigest(),
            "shards_per_stream": counts, "fingerprint_kind": FINGERPRINT_KIND}


def _dated_record(record, source: str):
    if not isinstance(record, dict):
        raise ValueError(f"{source}: date record must be an object")
    year = record.get("year")
    if year is None:
        return None
    if type(year) is not int or year <= 0:
        raise ValueError(f"{source}: a present publication year must be a positive integer or null")
    pmid = record.get("pmid")
    if not isinstance(pmid, str) or not pmid.strip():
        raise ValueError(f"{source}: dated record requires a nonempty string PMID")
    return pmid, year


def _scan_inventory(root: Path, inventory: list[dict]) -> dict:
    """Read exactly the captured inventory and check the versions opened."""
    years = {}
    for entry in inventory:
        path = root / entry["path"]
        with path.open("rb") as compressed:
            if _metadata(entry["path"], os.fstat(compressed.fileno())) != entry:
                raise ValueError(f"date source changed before opening: {entry['path']}")
            try:
                with gzip.GzipFile(fileobj=compressed, mode="rb") as zipped:
                    with io.TextIOWrapper(zipped, encoding="utf-8") as lines:
                        for line in lines:
                            dated = _dated_record(json.loads(line), entry["path"])
                            if dated is None:
                                continue
                            pmid, year = dated
                            if pmid not in years or year < years[pmid]:
                                years[pmid] = year
            except (EOFError, gzip.BadGzipFile, zlib.error) as error:
                raise ValueError(f"{entry['path']}: unreadable compressed date records: {error}") from error
            if _metadata(entry["path"], os.fstat(compressed.fileno())) != entry:
                raise ValueError(f"date source changed while reading: {entry['path']}")
    return years


def _scan_pmid_years(root: Path) -> dict:
    """Compatibility wrapper for an uncached, inventory-checked date scan."""
    root = Path(root)
    inventory = _inventory(root)
    years = _scan_inventory(root, inventory)
    _unchanged(root, inventory)
    return years


def _cached_years(payload, inventory):
    if (not isinstance(payload, dict)
            or type(payload.get("version")) is not int
            or payload["version"] != CACHE_VERSION
            or payload.get("streams") != list(YEAR_STREAMS)
            or payload.get("fingerprint_kind") != FINGERPRINT_KIND
            or payload.get("inventory") != inventory):
        return None
    years = payload.get("years")
    if not isinstance(years, dict):
        return None
    for pmid, year in years.items():
        if (not isinstance(pmid, str) or not pmid.strip()
                or type(year) is not int or year <= 0):
            return None
    return years


def _write_cache(cache: Path, payload: dict) -> None:
    """Replace a disposable cache only when its complete new payload is ready."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", dir=cache.parent,
                                         prefix=".pmid-years-", suffix=".tmp",
                                         delete=False) as handle:
            temporary = Path(handle.name)
            pickle.dump(payload, handle, protocol=5)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, cache)
    except (OSError, pickle.PickleError):
        # Cache persistence is optional, including on a read-only data root or
        # when only auxiliary streams exist and records/ is absent.
        pass
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass


def load_pmid_years(root: Path) -> tuple[dict, dict]:
    """Return publication years and the metadata provenance of the source set."""
    root = Path(root)
    inventory = _inventory(root)
    provenance = _provenance(inventory)
    if not inventory:
        # A leftover cache is not evidence that its absent sources are readable.
        _unchanged(root, inventory)
        return {}, provenance

    cache = root / "records" / ".pmid-years.pkl"
    try:
        with cache.open("rb") as handle:
            payload = pickle.load(handle)
    except Exception:
        # The cache is disposable. Source-reading and race errors remain
        # outside this catch and must not be replaced by stale cached dates.
        payload = None
    _unchanged(root, inventory)
    years = _cached_years(payload, inventory)
    if years is not None:
        _unchanged(root, inventory)
        print(f"  year map from cache ({len(years):,} PMIDs)", flush=True)
        return years, provenance

    years = _scan_inventory(root, inventory)
    _unchanged(root, inventory)
    _write_cache(cache, {"version": CACHE_VERSION, "streams": list(YEAR_STREAMS),
                         "fingerprint_kind": FINGERPRINT_KIND,
                         "inventory": inventory, "years": years})
    return years, provenance


def pmid_years(root: Path) -> dict:
    """Compatibility API used by discovery and the related temporal analyses."""
    return load_pmid_years(root)[0]
