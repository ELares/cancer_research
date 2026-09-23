"""Tiny date streams exercise cache freshness without opening the census."""

import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import atlas_discovery_dates as dates


def shard(root, stream="records", name="a", rows=None):
    path = root / stream / f"{name}.jsonl.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    content = "".join(json.dumps(row) + "\n" for row in
                      (rows if rows is not None else [{"pmid": "1", "year": 2025}]))
    previous = path.stat() if path.exists() else None
    path.write_bytes(gzip.compress(content.encode(), mtime=0))
    if previous is not None:
        # Make ordinary fixture changes visible even on coarse-clock filesystems.
        os.utime(path, ns=(previous.st_atime_ns, previous.st_mtime_ns + 2_000_000_000))
    return path


def cache_path(root):
    return root / "records" / ".pmid-years.pkl"


def fail(*args, **kwargs):
    pytest.fail("a warm cache read source contents")


@pytest.fixture
def populated(tmp_path):
    source = shard(tmp_path)
    years, provenance = dates.load_pmid_years(tmp_path)
    assert years == {"1": 2025}
    return tmp_path, source, provenance, cache_path(tmp_path).read_bytes()


def test_all_four_streams_merge_earliest_year_and_ignore_unselected_files(tmp_path):
    for stream, year in zip(dates.YEAR_STREAMS, [2025, 2020, 2018, 2010]):
        shard(tmp_path, stream, rows=[{"pmid": "shared", "year": year},
                                     {"pmid": stream, "year": year}])
    shard(tmp_path, "records", "z", [{"pmid": "shared", "year": 2027}])
    # Original glob.glob selection excluded hidden ingestion files.
    (tmp_path / "records" / ".partial.jsonl.gz").write_bytes(b"unfinished")
    (tmp_path / "records" / "notes.jsonl").write_bytes(b"unrelated")
    years, provenance = dates.load_pmid_years(tmp_path)
    assert years == {"shared": 2010, "records": 2025, "records_c04only": 2020,
                     "records_unindexed": 2018, "records_updates": 2010}
    assert provenance["shards_per_stream"] == dict(zip(dates.YEAR_STREAMS, [2, 1, 1, 1]))
    assert set(provenance) == {"inventory_sha256", "shards_per_stream", "fingerprint_kind"}
    assert provenance["fingerprint_kind"] == "filesystem-metadata-v1"
    assert dates.pmid_years(tmp_path) == years
    assert dates._scan_pmid_years(tmp_path) == years


def test_inventory_records_exact_file_metadata_in_sorted_order(tmp_path):
    paths = [shard(tmp_path, "records_updates", "z"), shard(tmp_path, "records", "a")]
    _, provenance = dates.load_pmid_years(tmp_path)
    payload = pickle.loads(cache_path(tmp_path).read_bytes())
    expected = []
    for path in sorted(paths):
        info = path.stat()
        expected.append({"path": path.relative_to(tmp_path).as_posix(),
                         "size": info.st_size, "dev": info.st_dev, "ino": info.st_ino,
                         "mtime_ns": info.st_mtime_ns, "ctime_ns": info.st_ctime_ns})
    assert payload["inventory"] == expected
    digest = hashlib.sha256(json.dumps(expected, sort_keys=True,
                                      separators=(",", ":")).encode()).hexdigest()
    assert provenance["inventory_sha256"] == digest


def test_warm_cache_never_opens_compressed_contents(populated, monkeypatch):
    root, source, provenance, before = populated
    original = Path.open

    def guarded(path, *args, **kwargs):
        if path.name.endswith(".jsonl.gz"):
            fail()
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded)
    monkeypatch.setattr(dates, "_scan_inventory", fail)
    assert dates.load_pmid_years(root) == ({"1": 2025}, provenance)
    assert cache_path(root).read_bytes() == before


@pytest.mark.parametrize("change", ["add_stream", "add_shard", "replace", "rewrite", "rename", "remove"])
def test_source_changes_invalidate_warm_cache(populated, change):
    root, source, provenance, before = populated
    expected = {"1": 2010}
    if change == "add_stream":
        shard(root, "records_updates", rows=[{"pmid": "1", "year": 2010}])
    elif change == "add_shard":
        shard(root, "records", "b", [{"pmid": "1", "year": 2010}])
    elif change == "replace":
        replacement = shard(root, "records", "replacement", [{"pmid": "1", "year": 2010}])
        info = source.stat()
        os.utime(replacement, ns=(info.st_atime_ns, info.st_mtime_ns))
        os.replace(replacement, source)
    elif change == "rewrite":
        info = source.stat()
        shard(root, rows=[{"pmid": "1", "year": 2010}])
        # Restoring mtime does not hide an in-place rewrite: ctime still changes.
        os.utime(source, ns=(info.st_atime_ns, info.st_mtime_ns))
        if dates._inventory(root) == pickle.loads(before)["inventory"]:
            pytest.skip("filesystem metadata cannot distinguish this rapid rewrite")
    elif change == "rename":
        source.rename(source.with_name("renamed.jsonl.gz"))
        expected = {"1": 2025}
    else:
        source.unlink()
        expected = {}
    years, after = dates.load_pmid_years(root)
    assert years == expected
    assert after["inventory_sha256"] != provenance["inventory_sha256"]
    if change != "remove":
        assert cache_path(root).read_bytes() != before


def test_empty_union_never_reads_stale_cache(populated, monkeypatch):
    root, source, _, before = populated
    source.unlink()
    monkeypatch.setattr(dates.pickle, "load", fail)
    years, provenance = dates.load_pmid_years(root)
    assert years == {}
    assert provenance["shards_per_stream"] == dict.fromkeys(dates.YEAR_STREAMS, 0)
    assert cache_path(root).read_bytes() == before
    assert dates.pmid_years(root / "missing") == {}


def test_removing_earliest_source_rebuilds_from_remaining_stream(populated):
    root, _, _, _ = populated
    update = shard(root, "records_updates", rows=[{"pmid": "1", "year": 2010}])
    assert dates.pmid_years(root) == {"1": 2010}
    update.unlink()
    assert dates.pmid_years(root) == {"1": 2025}


def test_missing_or_null_years_remain_undated_and_empty_map_can_be_cached(tmp_path, monkeypatch):
    shard(tmp_path, rows=[{}, {"pmid": "absent"}, {"year": None, "pmid": 4},
                         {"year": None, "pmid": ""}])
    first = dates.load_pmid_years(tmp_path)
    assert first[0] == {}
    monkeypatch.setattr(dates, "_scan_inventory", fail)
    assert dates.load_pmid_years(tmp_path) == first


def test_auxiliary_stream_alone_is_readable_without_cache_directory(tmp_path):
    shard(tmp_path, "records_updates", rows=[{"pmid": "1", "year": 2010}])
    assert dates.pmid_years(tmp_path) == {"1": 2010}
    assert not (tmp_path / "records").exists()


@pytest.mark.parametrize("mutation", ["corrupt", "legacy", "version", "bool_version", "streams",
                                       "kind", "inventory", "map", "year", "bool_year", "pmid"])
def test_invalid_cache_is_rebuilt_from_valid_sources(populated, mutation):
    root, _, _, before = populated
    payload = pickle.loads(before)
    payload["years"] = {"1": 1900}
    if mutation == "corrupt":
        cache_path(root).write_bytes(b"not a pickle")
    else:
        if mutation == "legacy":
            payload = {"streams": list(dates.YEAR_STREAMS), "years": {"1": 1900}}
        elif mutation == "version":
            payload["version"] += 1
        elif mutation == "bool_version":
            payload["version"] = True
        elif mutation == "streams":
            payload["streams"] = ["records"]
        elif mutation == "kind":
            payload["fingerprint_kind"] = "content"
        elif mutation == "inventory":
            payload["inventory"] = []
        elif mutation == "map":
            payload["years"] = []
        elif mutation == "year":
            payload["years"] = {"1": "1900"}
        elif mutation == "bool_year":
            payload["years"] = {"1": True}
        else:
            payload["years"] = {"": 1900}
        cache_path(root).write_bytes(pickle.dumps(payload))
    assert dates.pmid_years(root) == {"1": 2025}
    assert pickle.loads(cache_path(root).read_bytes())["years"] == {"1": 2025}


@pytest.mark.parametrize("row", [[], None, {"pmid": "1", "year": 0},
                                   {"pmid": "1", "year": -1}, {"pmid": "1", "year": True},
                                   {"pmid": "1", "year": 2020.0}, {"pmid": "1", "year": "2020"},
                                   {"year": 2020}, {"pmid": 1, "year": 2020},
                                   {"pmid": "", "year": 2020}, {"pmid": "  ", "year": 2020}])
def test_bad_raw_dates_fail_without_overwriting_previous_cache(populated, row):
    root, _, _, before = populated
    shard(root, rows=[row])
    with pytest.raises(ValueError):
        dates.load_pmid_years(root)
    assert cache_path(root).read_bytes() == before


@pytest.mark.parametrize("content", [b"not gzip", gzip.compress(b"{bad json}\n"),
                                     gzip.compress(b'{"pmid":"1","year":2010}\n')[:-4]])
def test_corrupt_or_truncated_shard_preserves_previous_cache(populated, content):
    root, source, _, before = populated
    source.write_bytes(content)
    with pytest.raises(ValueError):
        dates.load_pmid_years(root)
    assert cache_path(root).read_bytes() == before


@pytest.mark.parametrize("cache_state", ["absent", "previous", "matches_empty_source"])
def test_zero_byte_shard_is_rejected_even_when_other_stream_has_dates(tmp_path, cache_state):
    source = shard(tmp_path, rows=[{"pmid": "old", "year": 2000}])
    shard(tmp_path, "records_updates", rows=[{"pmid": "recent", "year": 2025}])
    assert dates.pmid_years(tmp_path) == {"old": 2000, "recent": 2025}
    cache = cache_path(tmp_path)
    source.write_bytes(b"")
    if cache_state == "absent":
        cache.unlink()
    elif cache_state == "matches_empty_source":
        # Earlier readers could cache a zero-byte gzip shard as an empty stream.
        payload = pickle.loads(cache.read_bytes())
        info = source.stat()
        entry = next(item for item in payload["inventory"]
                     if item["path"] == "records/a.jsonl.gz")
        entry.update(size=info.st_size, dev=info.st_dev, ino=info.st_ino,
                     mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns)
        payload["years"] = {"recent": 2025}
        cache.write_bytes(pickle.dumps(payload))
    before = cache.read_bytes() if cache.exists() else None
    with pytest.raises(ValueError, match="empty, not a valid gzip file: records/a.jsonl.gz"):
        dates.load_pmid_years(tmp_path)
    assert (cache.read_bytes() if cache.exists() else None) == before


def test_valid_empty_gzip_member_is_allowed_and_cached(tmp_path, monkeypatch):
    empty = shard(tmp_path, rows=[])
    assert empty.stat().st_size > 0
    shard(tmp_path, "records_updates", rows=[{"pmid": "recent", "year": 2025}])
    first = dates.load_pmid_years(tmp_path)
    assert first[0] == {"recent": 2025}
    assert first[1]["shards_per_stream"]["records"] == 1
    monkeypatch.setattr(dates, "_scan_inventory", fail)
    assert dates.load_pmid_years(tmp_path) == first


def test_non_directory_stream_and_non_regular_shard_are_not_absent(tmp_path):
    stream = tmp_path / "records"
    stream.write_text("not a directory")
    with pytest.raises(ValueError, match="not a directory"):
        dates.pmid_years(tmp_path)
    stream.unlink()
    (stream / "a.jsonl.gz").mkdir(parents=True)
    with pytest.raises(ValueError, match="not a regular file"):
        dates.pmid_years(tmp_path)


def test_unreadable_stream_is_not_absent(populated, monkeypatch):
    root, _, _, before = populated
    original = Path.iterdir

    def denied(path):
        if path == root / "records":
            raise PermissionError("unreadable directory")
        return original(path)

    monkeypatch.setattr(Path, "iterdir", denied)
    with pytest.raises(PermissionError):
        dates.pmid_years(root)
    assert cache_path(root).read_bytes() == before


def test_unreadable_shard_preserves_previous_cache(populated, monkeypatch):
    root, source, _, before = populated
    shard(root, rows=[{"pmid": "1", "year": 2010}])
    original = Path.open

    def denied(path, *args, **kwargs):
        if path == source:
            raise PermissionError("unreadable source")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", denied)
    with pytest.raises(PermissionError):
        dates.pmid_years(root)
    assert cache_path(root).read_bytes() == before


@pytest.mark.parametrize("stage", ["cache_load", "cache_validation"])
def test_source_change_during_cache_hit_is_rejected(populated, monkeypatch, stage):
    root, _, _, before = populated
    owner, name = ((dates.pickle, "load") if stage == "cache_load"
                   else (dates, "_cached_years"))
    original = getattr(owner, name)

    def raced(*args, **kwargs):
        result = original(*args, **kwargs)
        shard(root, "records_updates", rows=[{"pmid": "1", "year": 2010}])
        return result

    monkeypatch.setattr(owner, name, raced)
    with pytest.raises(ValueError, match="inventory changed"):
        dates.pmid_years(root)
    assert cache_path(root).read_bytes() == before


@pytest.mark.parametrize("stage", ["before_open", "during_read", "after_scan", "new_shard"])
def test_source_change_during_rebuild_is_rejected(populated, monkeypatch, stage):
    root, source, _, before = populated
    shard(root, rows=[{"pmid": "1", "year": 2018}])
    original_scan = dates._scan_inventory
    original_record = dates._dated_record

    def changed_record(*args):
        info = source.stat()
        os.utime(source, ns=(info.st_atime_ns, info.st_mtime_ns + 2_000_000_000))
        return original_record(*args)

    def changed_scan(*args):
        if stage == "before_open":
            shard(root, rows=[{"pmid": "1", "year": 2010}])
        if stage == "new_shard":
            # Must not be scanned opportunistically: it was not in inventory.
            (root / "records" / "new.jsonl.gz").write_bytes(b"unfinished")
        result = original_scan(*args)
        if stage == "after_scan":
            source.unlink()
        return result

    if stage == "during_read":
        monkeypatch.setattr(dates, "_dated_record", changed_record)
    else:
        monkeypatch.setattr(dates, "_scan_inventory", changed_scan)
    with pytest.raises(ValueError, match="changed"):
        dates.pmid_years(root)
    assert cache_path(root).read_bytes() == before


def test_cache_replacement_is_atomic_and_contains_complete_payload(populated, monkeypatch):
    root, _, _, before = populated
    shard(root, rows=[{"pmid": "1", "year": 2010}])
    original = os.replace
    calls = []

    def inspect(source, destination):
        assert Path(source).parent == cache_path(root).parent
        assert cache_path(root).read_bytes() == before
        assert pickle.loads(Path(source).read_bytes())["years"] == {"1": 2010}
        calls.append((source, destination))
        return original(source, destination)

    monkeypatch.setattr(dates.os, "replace", inspect)
    assert dates.pmid_years(root) == {"1": 2010}
    assert len(calls) == 1
    assert not list((root / "records").glob(".pmid-years-*.tmp"))


@pytest.mark.parametrize("stage", ["create", "serialize", "sync", "replace"])
def test_optional_cache_write_failure_preserves_old_cache_and_cleans_temp(populated, monkeypatch, stage):
    root, _, _, before = populated
    shard(root, rows=[{"pmid": "1", "year": 2010}])

    def denied(*args, **kwargs):
        raise PermissionError("cache is read-only")

    def interrupted(payload, handle, **kwargs):
        handle.write(b"partial pickle")
        raise pickle.PicklingError("serialization failed")

    owner, name, replacement = {
        "create": (dates.tempfile, "NamedTemporaryFile", denied),
        "serialize": (dates.pickle, "dump", interrupted),
        "sync": (dates.os, "fsync", denied),
        "replace": (dates.os, "replace", denied),
    }[stage]
    monkeypatch.setattr(owner, name, replacement)
    assert dates.pmid_years(root) == {"1": 2010}
    assert cache_path(root).read_bytes() == before
    assert not list((root / "records").glob(".pmid-years-*.tmp"))
