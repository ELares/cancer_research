"""Run the design and access CLIs against isolated, real census shards.

Unavailable input must preserve published reports, while readable populations
with no trial, mechanism, or PMC matches remain valid measurements.
"""
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent
SCANNERS = ("census_evidence_design", "census_oa_bias")
OLD_JSON = b'{"existing": "published counts"}\n'
OLD_MD = b"Existing published interpretation.\n"
UNMATCHED = {"pmid": "42", "mesh": ["Unmapped test descriptor"]}
MATCHED = {
    "pmid": "43", "mesh": ["Nanoparticles"],
    "pub_types": ["Clinical Trial, Phase III"],
    "pmcid": "PMC43", "year": 2020,
}


def _write_shard(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record) + "\n")


def _preserved(module):
    assert module.OUT_JSON.read_bytes() == OLD_JSON
    assert module.OUT_MD.read_bytes() == OLD_MD


@pytest.fixture(params=SCANNERS)
def scanner(request, tmp_path, monkeypatch):
    name = request.param
    root = tmp_path / "external-atlas"
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(root))
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        name + "_design_access_inputs", REPO / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Check the public root setting before any CLI could read live bulk data.
    assert module.RECORDS == root / "records"
    for attr, content in (("OUT_JSON", OLD_JSON), ("OUT_MD", OLD_MD)):
        path = tmp_path / attr
        path.write_bytes(content)
        monkeypatch.setattr(module, attr, path)
    monkeypatch.setattr(sys, "argv", [f"{name}.py"])
    return name, module


@pytest.mark.parametrize(
    "state", ["missing", "empty-directory", "empty-shard", "empty-selection"])
def test_unavailable_input_preserves_both_reports(scanner, state, monkeypatch):
    name, module = scanner
    if state != "missing":
        module.RECORDS.mkdir(parents=True)
    if state in {"empty-shard", "empty-selection"}:
        _write_shard(module.RECORDS / "part-000.jsonl.gz", [])
    if state == "empty-selection":
        _write_shard(module.RECORDS / "part-001.jsonl.gz", [MATCHED])
        monkeypatch.setattr(sys, "argv", [f"{name}.py", "--stride", "2"])

    with pytest.raises(SystemExit) as exc:
        module.main()

    message = str(exc.value)
    assert "No census records" in message
    assert str(module.RECORDS) in message
    assert "FERRO_ATLAS_ROOT" in message
    assert "atlas_baseline.py" in message and "--render-only" in message
    _preserved(module)


@pytest.mark.parametrize("stride", [0, -1])
def test_invalid_stride_preserves_both_reports(scanner, stride, monkeypatch):
    name, module = scanner
    _write_shard(module.RECORDS / "part-000.jsonl.gz", [MATCHED])
    monkeypatch.setattr(sys, "argv", [f"{name}.py", "--stride", str(stride)])

    with pytest.raises(SystemExit, match="positive integer"):
        module.main()

    _preserved(module)


@pytest.mark.parametrize("population", ["unmatched", "no-trial", "no-pmc"])
def test_readable_sparse_populations_are_valid(scanner, population):
    name, module = scanner
    record = {
        "unmatched": UNMATCHED,
        "no-trial": {**MATCHED, "pub_types": ["Journal Article"]},
        "no-pmc": {**MATCHED, "pmcid": None},
    }[population]
    _write_shard(module.RECORDS / "part-000.jsonl.gz", [record])
    # These are distinct populations and must never join the indexed stream.
    for stream in ("records_unindexed", "records_updates", "records_c04only"):
        extra = module.RECORDS.parent / stream / "part.jsonl.gz"
        extra.parent.mkdir()
        extra.write_bytes(b"not gzip")

    assert module.main() == 0

    data = json.loads(module.OUT_JSON.read_text())
    assert data["census"] == 1
    if name == "census_evidence_design":
        assert sum(data["classes"].values()) == 1
        assert data["classes"]["trial"] == (1 if population == "no-pmc" else 0)
        assert data["classifiable"] == (1 if population == "no-pmc" else 0)
    else:
        assert data["with_pmcid"] == (1 if population == "no-trial" else 0)
        assert data["census_oa_rate"] == data["with_pmcid"]
        if population == "unmatched":
            assert data["mechanisms"] == []
        else:
            nanoparticle = next(row for row in data["mechanisms"]
                                if row["mechanism"] == "nanoparticle")
            assert nanoparticle["total"] == 1
            assert nanoparticle["with_fulltext"] == data["with_pmcid"]
            assert nanoparticle["without"] == 1 - data["with_pmcid"]
    assert module.OUT_MD.read_text().startswith("# ")
    assert module.OUT_MD.read_bytes() != OLD_MD


@pytest.mark.parametrize("corruption", ["later-json", "later-gzip"])
def test_bad_later_input_preserves_both_reports(scanner, corruption):
    _, module = scanner
    _write_shard(module.RECORDS / "part-000.jsonl.gz", [MATCHED])
    later = module.RECORDS / "part-001.jsonl.gz"
    if corruption == "later-json":
        with gzip.open(later, "wt", encoding="utf-8") as stream:
            stream.write("{invalid json\n")
        error = json.JSONDecodeError
    else:
        later.write_bytes(b"not gzip")
        error = gzip.BadGzipFile

    with pytest.raises(error):
        module.main()

    _preserved(module)


@pytest.mark.parametrize("render_only", [False, True], ids=["scan", "offline"])
def test_render_failure_preserves_both_reports(scanner, render_only, monkeypatch):
    name, module = scanner
    before_json = OLD_JSON
    if render_only:
        data = json.loads(
            (REPO / "analysis" / f"{name.replace('_', '-')}.json").read_text())
        # Noncanonical whitespace detects writes even when values stay equal.
        before_json = json.dumps(data, separators=(",", ":")).encode() + b"\n"
        module.OUT_JSON.write_bytes(before_json)
        monkeypatch.setattr(sys, "argv", [f"{name}.py", "--render-only"])
    else:
        _write_shard(module.RECORDS / "part-000.jsonl.gz", [MATCHED])

    def fail(_data):
        raise RuntimeError("intentional render failure")

    monkeypatch.setattr(module, "render", fail)
    with pytest.raises(RuntimeError, match="intentional render failure"):
        module.main()

    assert module.OUT_JSON.read_bytes() == before_json
    assert module.OUT_MD.read_bytes() == OLD_MD


def test_stride_selects_sorted_shards_not_records_or_a_prefix(
        scanner, tmp_path, monkeypatch):
    name, module = scanner
    # Unequal shard lengths and creation order distinguish all three choices.
    counts = {0: 1, 1: 7, 2: 2, 3: 9, 4: 4}
    for i in (3, 0, 4, 1, 2):
        _write_shard(module.RECORDS / f"part-{i:03}.jsonl.gz", [MATCHED] * counts[i])
    monkeypatch.setattr(sys, "argv", [f"{name}.py", "--stride", "2"])

    assert module.main() == 0

    selected = json.loads(module.OUT_JSON.read_text())
    assert selected["census"] == 7
    reference = tmp_path / "selected-records"
    _write_shard(reference / "selected.jsonl.gz", [MATCHED] * 7)
    monkeypatch.setattr(module, "RECORDS", reference)
    monkeypatch.setattr(sys, "argv", [f"{name}.py"])
    assert module.main() == 0
    assert json.loads(module.OUT_JSON.read_text()) == selected


def test_offline_render_preserves_committed_counts_without_raw_input(
        scanner, monkeypatch):
    name, module = scanner
    basename = name.replace("_", "-")
    original = (REPO / "analysis" / f"{basename}.json").read_bytes()
    module.OUT_JSON.write_bytes(original)
    monkeypatch.setattr(sys, "argv", [f"{name}.py", "--render-only"])

    def forbidden(*_args, **_kwargs):
        pytest.fail("Offline rendering must not scan raw census records")

    monkeypatch.setattr(module, "scan", forbidden)
    assert not module.RECORDS.exists()
    assert module.main() == 0
    assert json.loads(module.OUT_JSON.read_bytes()) == json.loads(original)
    assert module.OUT_MD.read_bytes() == (
        REPO / "analysis" / f"{basename}.md").read_bytes()
