"""A review snapshot binds complete parsed input to its compressed source bytes."""
import builtins
import gzip
import hashlib
import io
import json
import os
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]


def write_snapshot(tmp_path, shards):
    root = tmp_path / "atlas"
    records = root / "records"
    records.mkdir(parents=True)
    entries = {}
    for index, (name, rows) in enumerate(shards.items()):
        payload = "".join(json.dumps(row) + "\n" for row in rows).encode()
        (records / name).write_bytes(gzip.compress(payload, mtime=0))
        entries[f"source-{index}.xml.gz"] = {
            "parsed": True, "cancer": len(rows), "records": name,
        }
    manifest = {"source": "https://example.org/pubmed/baseline/", "files": entries}
    (root / "manifest.json").write_text(json.dumps(manifest))
    return root


def _update_manifest(root, change):
    path = root / "manifest.json"
    manifest = json.loads(path.read_text())
    change(manifest)
    path.write_text(json.dumps(manifest))


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    separators=(",", ":"), allow_nan=False).encode()).hexdigest()


@pytest.fixture
def reader(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    import census_snapshot
    return census_snapshot


def test_snapshot_binds_all_parsed_shards_while_retaining_only_selected_records(reader, tmp_path):
    rows = [{"pmid": "1", "cancer_basis": "c04"},
            {"pmid": "2", "cancer_basis": "adjacent"}]
    root = write_snapshot(tmp_path, {"part-b.jsonl.gz": [rows[1]], "part-a.jsonl.gz": [rows[0]]})
    manifest_bytes = (root / "manifest.json").read_bytes()

    report, selected = reader.read_snapshot(root, lambda record: record["pmid"] == "2")

    assert selected == [rows[1]]
    assert report["schema_version"] == 1
    assert report["stream"] == "records"
    assert report["stride"] == 1
    assert report["declared_records"] == report["observed_records"] == 2
    assert report["cancer_basis"] == {"c04": 1, "adjacent": 1}
    assert report["manifest_sha256"] == hashlib.sha256(manifest_bytes).hexdigest()
    assert report["baseline_manifest_sha256"] == _digest(json.loads(manifest_bytes)["files"])
    expected = []
    for name in ("part-a.jsonl.gz", "part-b.jsonl.gz"):
        compressed = (root / "records" / name).read_bytes()
        expected.append({"name": name, "sha256": hashlib.sha256(compressed).hexdigest(),
                         "bytes": len(compressed), "records": 1})
    assert report["shards"] == expected
    assert report["inventory_sha256"] == _digest(expected)


def test_empty_shards_and_zero_selected_records_are_valid_in_a_nonempty_snapshot(reader, tmp_path):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [], "part-b.jsonl.gz": [{"pmid": "1"}]})
    report, selected = reader.read_snapshot(root, lambda record: False)
    assert selected == []
    assert report["declared_records"] == report["observed_records"] == 1
    assert [shard["records"] for shard in report["shards"]] == [0, 1]
    assert report["cancer_basis"] == {"unspecified": 1}


def test_update_entries_are_excluded_from_the_baseline_inventory_and_count(reader, tmp_path):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}]})
    original, _ = reader.read_snapshot(root, lambda record: True)
    _update_manifest(root, lambda m: m["files"].update({"update.xml.gz": {
        "parsed": True, "cancer": 999, "records": "update.jsonl.gz", "source": "updatefiles",
    }}))
    report, selected = reader.read_snapshot(root, lambda record: True)
    assert report["observed_records"] == report["declared_records"] == len(selected) == 1
    assert report["manifest_sha256"] != original["manifest_sha256"]
    assert report["baseline_manifest_sha256"] == original["baseline_manifest_sha256"]
    assert report["inventory_sha256"] == original["inventory_sha256"]


@pytest.mark.parametrize("fault", ["empty-baseline", "unknown-source", "unparsed", "negative-count",
                                  "boolean-count", "fractional-count", "traversal", "absolute-path",
                                  "duplicate-reference", "blank-source-url"])
def test_invalid_manifest_cannot_describe_a_review_snapshot(reader, tmp_path, fault):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}]})

    def change(manifest):
        entry = manifest["files"]["source-0.xml.gz"]
        if fault == "empty-baseline":
            manifest["files"] = {}
        elif fault == "blank-source-url":
            manifest["source"] = " "
        elif fault == "duplicate-reference":
            manifest["files"]["another.xml.gz"] = dict(entry)
        else:
            key, value = {
                "unknown-source": ("source", "other"), "unparsed": ("parsed", False),
                "negative-count": ("cancer", -1), "boolean-count": ("cancer", True),
                "fractional-count": ("cancer", 1.5),
                "traversal": ("records", "../records/part-a.jsonl.gz"),
                "absolute-path": ("records", str(root / "records/part-a.jsonl.gz")),
            }[fault]
            entry[key] = value

    _update_manifest(root, change)
    with pytest.raises(SystemExit):
        reader.read_snapshot(root, lambda record: True)


@pytest.mark.parametrize("fault", ["missing", "extra", "symlink-alias", "hardlink-alias"])
def test_inventory_must_contain_each_distinct_declared_file_exactly_once(reader, tmp_path, fault):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}]})
    original = root / "records/part-a.jsonl.gz"
    other = root / "records/part-b.jsonl.gz"
    if fault == "missing":
        original.unlink()
    elif fault == "extra":
        other.write_bytes(original.read_bytes())
    else:
        if fault == "symlink-alias":
            other.symlink_to(original)
        else:
            os.link(original, other)
        _update_manifest(root, lambda m: m["files"].update({"other.xml.gz": {
            "parsed": True, "cancer": 1, "records": other.name,
        }}))
    with pytest.raises(SystemExit):
        reader.read_snapshot(root, lambda record: True)


@pytest.mark.parametrize("declared", [0, 2])
def test_each_shard_count_must_match_its_manifest_even_with_no_selected_records(reader, tmp_path, declared):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}]})
    _update_manifest(root, lambda m: m["files"]["source-0.xml.gz"].update(cancer=declared))
    with pytest.raises(SystemExit):
        reader.read_snapshot(root, lambda record: False)


@pytest.mark.parametrize("payload", [b"not gzip", gzip.compress(b'{"pmid":"2"}\n{bad json\n'),
                                    gzip.compress(b"[]\n"), gzip.compress(b'{"value":NaN}\n')])
def test_late_unselected_input_must_still_parse_as_complete_record_objects(reader, tmp_path, payload):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}],
                                     "part-b.jsonl.gz": [{"pmid": "2"}]})
    (root / "records/part-b.jsonl.gz").write_bytes(payload)
    with pytest.raises(SystemExit):
        reader.read_snapshot(root, lambda record: False)


def test_an_entirely_empty_snapshot_is_unavailable(reader, tmp_path):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": []})
    with pytest.raises(SystemExit):
        reader.read_snapshot(root, lambda record: True)


def test_zero_byte_file_is_not_a_valid_empty_gzip_shard(reader, tmp_path):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}], "part-b.jsonl.gz": []})
    (root / "records/part-b.jsonl.gz").write_bytes(b"")
    with pytest.raises(SystemExit):
        reader.read_snapshot(root, lambda record: True)


@pytest.mark.parametrize("change", ["manifest", "inventory", "shard-metadata"])
def test_changes_during_record_selection_invalidate_the_snapshot(reader, tmp_path, change):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}]})

    def select(record):
        if change == "manifest":
            _update_manifest(root, lambda m: m.update(source="https://changed.example/"))
        elif change == "inventory":
            (root / "records/extra.jsonl.gz").write_bytes(gzip.compress(b""))
        else:
            path = root / "records/part-a.jsonl.gz"
            stat = path.stat()
            os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        return True

    with pytest.raises(SystemExit):
        reader.read_snapshot(root, select)


def test_compression_bytes_change_snapshot_identity_without_changing_records(reader, tmp_path):
    rows = [{"pmid": "1", "title": "Same parsed record"}]
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": rows})
    first, retained = reader.read_snapshot(root, lambda record: True)
    path = root / "records/part-a.jsonl.gz"
    path.write_bytes(gzip.compress(gzip.decompress(path.read_bytes()), mtime=123))
    second, reretained = reader.read_snapshot(root, lambda record: True)
    assert retained == reretained == rows
    assert first["manifest_sha256"] == second["manifest_sha256"]
    assert first["shards"][0]["sha256"] != second["shards"][0]["sha256"]
    assert first["inventory_sha256"] != second["inventory_sha256"]


def test_shard_is_opened_once_for_both_hashing_and_parsing(reader, tmp_path, monkeypatch):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1"}]})
    shard = root / "records/part-a.jsonl.gz"
    opened = []
    depth = 0

    def watch(original):
        def open_file(file, *args, **kwargs):
            nonlocal depth
            if depth == 0 and isinstance(file, (str, os.PathLike)) and Path(file) == shard:
                opened.append(Path(file))
                assert len(opened) == 1, "hashing and parsing must share one compressed read"
            depth += 1
            try:
                return original(file, *args, **kwargs)
            finally:
                depth -= 1
        return open_file

    # Older pathlib caches its opener; count public entry points without
    # double-counting their delegation to io.open on newer Python versions.
    monkeypatch.setattr(Path, "open", watch(Path.open))
    monkeypatch.setattr(builtins, "open", watch(builtins.open))
    monkeypatch.setattr(io, "open", watch(io.open))
    report, records = reader.read_snapshot(root, lambda record: True)
    assert opened == [shard]
    assert records == [{"pmid": "1"}]
    assert report["observed_records"] == 1
