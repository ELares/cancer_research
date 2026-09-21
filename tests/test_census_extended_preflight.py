"""Exercise census CLIs against unavailable, sparse, and malformed real inputs.

Each report is isolated from published artifacts. Only external PubMed responses
are stubbed; gzip parsing, classification, derived metrics and rendering run.
"""
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest

REPO = Path(__file__).resolve().parent.parent
SCANNERS = (
    "census_fulltext_ceiling",
    "census_mechanism_cancer_matrix",
    "census_external_check",
    "census_translation_lag",
    "census_synergy_metrics",
    "census_diagnostic_chains",
)
UNMATCHED = {"pmid": "42", "mesh": ["Unmapped test descriptor"]}
OLD_JSON = b'{"existing": "published counts"}\n'
OLD_MD = b"Existing published interpretation.\n"


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
    spec = importlib.util.spec_from_file_location(name + "_preflight", REPO / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Fail before invoking a scanner if it might read the real bulk store.
    assert module.RECORDS == root / "records"
    for attr, content in (("OUT_JSON", OLD_JSON), ("OUT_MD", OLD_MD)):
        path = tmp_path / attr
        path.write_bytes(content)
        monkeypatch.setattr(module, attr, path)
    if name == "census_external_check":
        monkeypatch.setattr(module, "pubmed_count", Mock(return_value=0))
        monkeypatch.setattr(module, "PAUSE", 0)
    if name == "census_diagnostic_chains":
        import tag_articles
        corpus = tmp_path / "frozen-corpus"
        corpus.mkdir()
        (corpus / "42.md").write_text(
            "---\npmid: '42'\ntitle: Unmatched fixture\nmesh_terms: []\n---\n\n## Abstract\nNo matching chain.\n")
        monkeypatch.setattr(tag_articles, "PMID_DIR", corpus)
    monkeypatch.setattr(sys, "argv", [f"{name}.py"])
    return name, module


@pytest.mark.parametrize("state", ["missing", "empty-directory", "empty-shard", "empty-selection"])
def test_unavailable_input_preserves_published_reports(scanner, state, monkeypatch):
    name, module = scanner
    if state != "missing":
        module.RECORDS.mkdir(parents=True)
    if state in {"empty-shard", "empty-selection"}:
        _write_shard(module.RECORDS / "part-000.jsonl.gz", [])
    if state == "empty-selection":
        _write_shard(module.RECORDS / "part-001.jsonl.gz", [UNMATCHED])
        monkeypatch.setattr(sys, "argv", [f"{name}.py", "--stride", "2"])
    with pytest.raises(SystemExit) as exc:
        module.main()
    message = str(exc.value)
    assert "No census records" in message
    assert str(module.RECORDS) in message
    assert "FERRO_ATLAS_ROOT" in message and "--render-only" in message
    _preserved(module)
    if name == "census_external_check":
        module.pubmed_count.assert_not_called()


@pytest.mark.parametrize("stride", [0, -1])
def test_invalid_stride_fails_before_writing_or_fetching(scanner, stride, monkeypatch):
    name, module = scanner
    _write_shard(module.RECORDS / "part-000.jsonl.gz", [UNMATCHED])
    monkeypatch.setattr(sys, "argv", [f"{name}.py", "--stride", str(stride)])
    with pytest.raises(SystemExit, match="positive integer"):
        module.main()
    _preserved(module)
    if name == "census_external_check":
        module.pubmed_count.assert_not_called()


def test_readable_input_with_zero_scientific_matches_is_valid(scanner):
    name, module = scanner
    _write_shard(module.RECORDS / "part-000.jsonl.gz", [UNMATCHED])
    # The text-recovered census is a separate stream and must not be included.
    extra = module.RECORDS.parent / "unindexed" / "part.jsonl.gz"
    extra.parent.mkdir()
    extra.write_bytes(b"not a gzip file")
    assert module.main() == 0
    data = json.loads(module.OUT_JSON.read_text())
    count_key = {
        "census_external_check": "census_records",
        "census_diagnostic_chains": "census_records",
        "census_synergy_metrics": "ferroptosis_articles",
    }.get(name, "census")
    assert data[count_key] == (0 if name == "census_synergy_metrics" else 1)
    assert module.OUT_MD.read_text().startswith("# ")
    assert module.OUT_MD.read_bytes() != OLD_MD


def test_offline_render_uses_committed_counts_without_census_or_network(scanner, monkeypatch):
    name, module = scanner
    original = (REPO / "analysis" / f"{name.replace('_', '-')}.json").read_bytes()
    module.OUT_JSON.write_bytes(original)
    monkeypatch.setattr(sys, "argv", [f"{name}.py", "--render-only"])
    assert not module.RECORDS.exists()
    assert module.main() == 0
    assert json.loads(module.OUT_JSON.read_bytes()) == json.loads(original)
    assert module.OUT_MD.read_text() == module.render(json.loads(module.OUT_JSON.read_text()))
    if name == "census_external_check":
        module.pubmed_count.assert_not_called()


@pytest.mark.parametrize("render_only", [False, True], ids=["scan", "offline"])
def test_render_errors_preserve_both_artifacts(scanner, monkeypatch, render_only):
    name, module = scanner
    before_json = OLD_JSON
    if render_only:
        data = json.loads((REPO / "analysis" / f"{name.replace('_', '-')}.json").read_text())
        before_json = json.dumps(data, separators=(",", ":")).encode() + b"\n"
        module.OUT_JSON.write_bytes(before_json)
        monkeypatch.setattr(sys, "argv", [f"{name}.py", "--render-only"])
    else:
        _write_shard(module.RECORDS / "part-000.jsonl.gz", [UNMATCHED])
    def fail(_data):
        raise RuntimeError("intentional render failure")
    monkeypatch.setattr(module, "render", fail)
    with pytest.raises(RuntimeError, match="intentional render failure"):
        module.main()
    assert module.OUT_JSON.read_bytes() == before_json
    assert module.OUT_MD.read_bytes() == OLD_MD


@pytest.mark.parametrize("corruption", ["later-json", "later-gzip"])
def test_bad_later_input_preserves_both_artifacts(scanner, corruption):
    name, module = scanner
    _write_shard(module.RECORDS / "part-000.jsonl.gz", [UNMATCHED])
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
    if name == "census_external_check":
        module.pubmed_count.assert_not_called()


def _scan(name, module, stride):
    entry = {"census_external_check": "census_counts", "census_diagnostic_chains": "scan_census"}.get(name, "scan")
    return getattr(module, entry)(stride)


def test_stride_selects_sorted_shards_not_records_or_a_prefix(scanner):
    name, module = scanner
    matched = {"pmid": "42", "mesh": ["Ferroptosis"], "cancer_basis": "C04"}
    # Creation order differs from filename order; unequal file lengths expose
    # record sampling, prefix sampling, and unsorted filesystem traversal.
    counts = {0: 1, 1: 7, 2: 2, 3: 9, 4: 4}
    for i in (3, 0, 4, 1, 2):
        _write_shard(module.RECORDS / f"part-{i:03}.jsonl.gz", [matched] * counts[i])
    data = _scan(name, module, 2)
    key = {"census_diagnostic_chains": "records", "census_synergy_metrics": "ferroptosis_articles"}.get(name, "census")
    assert data[key] == 7
    if name == "census_diagnostic_chains":
        assert data["shards"] == 3


@pytest.mark.parametrize("scanner", ["census_synergy_metrics"], indirect=True)
def test_synergy_control_keeps_its_own_stride(scanner):
    name, module = scanner
    matched = {"pmid": "42", "mesh": ["Ferroptosis"], "title": "Bliss independence"}
    for i in range(41):
        _write_shard(module.RECORDS / f"part-{i:03}.jsonl.gz", [matched])
    data = module.scan(2)
    assert data["ferroptosis_articles"] == 21
    assert data["control_records"] == 2
    assert data["control_stride"] == 40
    assert data["control_metric"] == {"Bliss independence": 2}


@pytest.mark.parametrize("scanner", ["census_synergy_metrics"], indirect=True)
@pytest.mark.parametrize("control_state", ["empty", "corrupt"])
def test_synergy_control_must_finish_before_either_report_is_written(
        scanner, control_state, monkeypatch):
    name, module = scanner
    if control_state == "empty":
        # The subject sees the second shard, but the control selects only the
        # empty first shard. Availability must be checked for both populations.
        _write_shard(module.RECORDS / "part-000.jsonl.gz", [])
        _write_shard(module.RECORDS / "part-001.jsonl.gz", [UNMATCHED])
        error = SystemExit
    else:
        for i in range(41):
            _write_shard(module.RECORDS / f"part-{i:03}.jsonl.gz", [UNMATCHED])
        with gzip.open(module.RECORDS / "part-040.jsonl.gz", "wt") as stream:
            stream.write("{invalid json\n")
        # The subject stride skips shard 40; the control includes it.
        monkeypatch.setattr(sys, "argv", [f"{name}.py", "--stride", "3"])
        error = json.JSONDecodeError
    with pytest.raises(error):
        module.main()
    _preserved(module)
