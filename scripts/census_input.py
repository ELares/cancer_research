"""Read a declared census stream without mistaking unavailable input for zero matches.

Consumers must exhaust the iterator before writing reports: an existing set of
gzip shards can still contain no records. Counts are checked before any
consumer-specific filtering, so a readable census with no matching mechanisms
is valid. The default selects only indexed ``records/``. Callers may supply an
explicit directory to validate other streams separately; this reader never
merges them or changes a caller's population definition.
"""
import gzip
import json
from pathlib import Path
from typing import Iterator

from atlas_baseline import atlas_root


def _missing_input(records: Path) -> SystemExit:
    return SystemExit(
        f"No census records found at {records}; existing reports were not changed. "
        "Set FERRO_ATLAS_ROOT to the census data root, acquire the census with "
        "scripts/atlas_baseline.py, or use --render-only to regenerate from "
        "the committed counts."
    )


def census_shards(records: Path, stride: int = 1) -> list[Path]:
    """Select sorted shards; the stride samples files, not individual records."""
    if stride < 1:
        raise SystemExit("--stride must be a positive integer")
    shards = sorted(records.glob("*.jsonl.gz"))[::stride]
    if not shards:
        raise _missing_input(records)
    return shards


def iter_census_shards(shards: list[Path], records: Path) -> Iterator[dict]:
    """Read one selected inventory, including when its size is reported."""
    have_records = False
    for shard in shards:
        with gzip.open(shard, "rt", encoding="utf-8") as fh:
            for line in fh:
                record = json.loads(line)
                have_records = True
                yield record
    if not have_records:
        raise _missing_input(records)


def iter_census_records(records: Path | None = None, stride: int = 1) -> Iterator[dict]:
    """Stream parsed records from the selected census root, requiring input."""
    records = atlas_root() / "records" if records is None else records
    yield from iter_census_shards(census_shards(records, stride), records)
