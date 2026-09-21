"""Fresh direction reports must use labels for exactly the scanned articles.

Tiny gzip fixtures exercise the public CLI, classification, candidate export,
adjudication validation and report writes for both direction measurements.
Published artifacts are used only as read-only historical render fixtures.
"""
import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO = Path(__file__).resolve().parent.parent
SCANNERS = ("census_hypoxia_direction", "census_thesis_direction")
CSV_FIELDS = [
    "cohort_sha256", "record_sha256", "pmid", "title", "abstract",
    "regex_label", "adjudicated", "reason",
]
OLD_JSON = b'{"existing": "published counts"}\n'
OLD_MD = b"Existing published interpretation.\n"
UNMATCHED = {"pmid": "99", "mesh": ["Unmapped fixture descriptor"]}


def _write_shard(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record) + "\n")


def _write_csv(path, rows, fields=CSV_FIELDS):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        return reader.fieldnames, list(reader)


def _snapshot(scanner):
    return scanner.json.read_bytes(), scanner.md.read_bytes()


def _reset_reports(scanner):
    scanner.json.write_bytes(OLD_JSON)
    scanner.md.write_bytes(OLD_MD)


def _run(scanner, monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", [
        scanner.name + ".py", "--output-dir", str(scanner.output), *map(str, args),
    ])
    return scanner.module.main()


@pytest.fixture(params=SCANNERS)
def scanner(request, tmp_path, monkeypatch):
    name = request.param
    root = tmp_path / "external-atlas"
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(root))
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        name + "_direction_input_test", REPO / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Check the public setting before any call could reach the real bulk store.
    assert module.RECORDS == root / "records"
    output = tmp_path / "reports"
    output.mkdir()
    basename = name.replace("_", "-")
    json_path, md_path = output / (basename + ".json"), output / (basename + ".md")
    monkeypatch.setattr(module, "OUT_JSON", json_path)
    monkeypatch.setattr(module, "OUT_MD", md_path)
    isolated_repo = tmp_path / "isolated-repository"
    analysis = isolated_repo / "analysis"
    analysis.mkdir(parents=True)
    legacy_name = name.removeprefix("census_").replace("_", "-") + "-adjudication.csv"
    legacy = analysis / legacy_name
    legacy.write_bytes((REPO / "analysis" / legacy_name).read_bytes())
    monkeypatch.setattr(module, "REPO", isolated_repo)
    if name == "census_hypoxia_direction":
        labels = ("protects", "sensitises")
        titles = ("Hypoxia protects tumour cells", "Hypoxia promotes ferroptosis")
        mesh = ["Ferroptosis", "Hypoxia"]
    else:
        labels = ("exploit", "obstacle")
        titles = ("Overcoming drug resistance through ferroptosis",
                  "Ferroptosis resistance limits treatment")
        mesh = ["Ferroptosis", "Drug Resistance, Neoplasm"]
    candidates = [
        {"pmid": str(i + 1), "year": 2025, "mesh": mesh,
         "title": title, "abstract": "An experimental observation."}
        for i, title in enumerate(titles)
    ]
    neither = {**candidates[0], "pmid": "3", "title": "A study of this pathway"}
    both = {**candidates[0], "pmid": "4", "title": ". ".join(titles)}
    result = SimpleNamespace(
        name=name, module=module, records=root / "records", output=output,
        json=json_path, md=md_path, legacy=legacy, labels=labels,
        candidates=candidates, neither=neither, both=both,
    )
    _reset_reports(result)
    return result


def _export_completed(scanner, monkeypatch, tmp_path):
    path = tmp_path / "current-labels.csv"
    assert _run(scanner, monkeypatch, "--export-candidates", path) == 0
    fields, rows = _read_csv(path)
    assert fields == CSV_FIELDS
    assert rows
    for row in rows:
        row["adjudicated"] = row["regex_label"]
        row["reason"] = "The title and abstract explicitly state this direction."
    _write_csv(path, rows)
    return path, rows


@pytest.mark.parametrize("state", ["missing", "empty-directory", "empty-shard", "empty-selection"])
def test_unavailable_census_preserves_reports(scanner, monkeypatch, state):
    if state != "missing":
        scanner.records.mkdir(parents=True)
    if state in {"empty-shard", "empty-selection"}:
        _write_shard(scanner.records / "part-000.jsonl.gz", [])
    args = []
    if state == "empty-selection":
        _write_shard(scanner.records / "part-001.jsonl.gz", [UNMATCHED])
        args = ["--stride", "2"]

    with pytest.raises(SystemExit) as exc:
        _run(scanner, monkeypatch, *args)

    message = str(exc.value)
    assert "No census records" in message
    assert str(scanner.records) in message
    assert "FERRO_ATLAS_ROOT" in message and "--render-only" in message
    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("stride", [0, -1])
def test_invalid_stride_preserves_reports(scanner, monkeypatch, capsys, stride):
    _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)

    with pytest.raises(SystemExit) as exc:
        _run(scanner, monkeypatch, "--stride", stride)

    assert "positive integer" in str(exc.value) + capsys.readouterr().err
    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("corruption", ["later-json", "later-gzip"])
def test_late_parser_failure_preserves_reports_and_export(scanner, monkeypatch, tmp_path, corruption):
    _write_shard(scanner.records / "part-000.jsonl.gz", scanner.candidates)
    later = scanner.records / "part-001.jsonl.gz"
    if corruption == "later-json":
        with gzip.open(later, "wt", encoding="utf-8") as stream:
            stream.write("{invalid json\n")
        error = json.JSONDecodeError
    else:
        later.write_bytes(b"not a gzip file")
        error = gzip.BadGzipFile
    export = tmp_path / "candidates.csv"
    export.write_bytes(b"Existing candidates.\n")

    with pytest.raises(error):
        _run(scanner, monkeypatch, "--export-candidates", export)

    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)
    assert export.read_bytes() == b"Existing candidates.\n"


@pytest.mark.parametrize("population", ["zero-match", "neither-only", "both-only"])
def test_no_candidates_is_valid_and_never_attaches_legacy_labels(
        scanner, monkeypatch, tmp_path, population):
    record = {"zero-match": UNMATCHED, "neither-only": scanner.neither,
              "both-only": scanner.both}[population]
    _write_shard(scanner.records / "part.jsonl.gz", [record])
    # A malformed unindexed stream must not affect the indexed census.
    unindexed = scanner.records.parent / "unindexed"
    unindexed.mkdir()
    (unindexed / "part.jsonl.gz").write_bytes(b"not gzip")
    export = tmp_path / "empty-candidates.csv"

    assert _run(scanner, monkeypatch, "--export-candidates", export) == 0

    data = json.loads(scanner.json.read_text())
    assert data["total"] == (0 if population == "zero-match" else 1)
    assert data["classified"] == 0
    assert data["adjudication"] == {}
    assert data["cohort"]
    fields, rows = _read_csv(export)
    assert fields == CSV_FIELDS and rows == []
    report = scanner.md.read_text()
    assert report.startswith("# ")
    assert "44.4%" not in report and "44%" not in report
    assert "None%" not in report
    assert "**Corrected, the exploit share" not in report


def test_fresh_candidates_preserve_raw_counts_without_attaching_historical_labels(
        scanner, monkeypatch, tmp_path):
    records = [*scanner.candidates, scanner.neither, scanner.both, UNMATCHED]
    _write_shard(scanner.records / "part.jsonl.gz", records)
    export = tmp_path / "candidates.csv"
    legacy_before = scanner.legacy.read_bytes()

    assert _run(scanner, monkeypatch, "--export-candidates", export) == 0

    data = json.loads(scanner.json.read_text())
    assert data["total"] == 4 and data["classified"] == 2
    assert data["counts"] == {scanner.labels[0]: 1, scanner.labels[1]: 1,
                              "both": 1, "neither": 1}
    assert data["adjudication"] == {}
    fields, rows = _read_csv(export)
    assert fields == CSV_FIELDS
    assert {row["pmid"] for row in rows} == {"1", "2"}
    assert {row["regex_label"] for row in rows} == set(scanner.labels)
    assert len({row["cohort_sha256"] for row in rows}) == 1
    for row in rows:
        assert len(row["cohort_sha256"]) == len(row["record_sha256"]) == 64
        assert row["title"] and row["abstract"]
        assert row["adjudicated"] == row["reason"] == ""
    assert scanner.legacy.read_bytes() == legacy_before


def test_complete_exported_csv_round_trips_into_a_fresh_measurement(scanner, monkeypatch, tmp_path):
    _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)
    path, rows = _export_completed(scanner, monkeypatch, tmp_path)
    # Assign every candidate the first direction, so labels visibly change the
    # result rather than merely reproducing the 50/50 keyword count.
    for row in rows:
        row["adjudicated"] = scanner.labels[0]
    _write_csv(path, rows)

    assert _run(scanner, monkeypatch, "--adjudication", path) == 0

    data = json.loads(scanner.json.read_text())
    assert data["classified"] == 2
    assert data["adjudication"]
    if scanner.name == "census_hypoxia_direction":
        assert data["adjudication"]["n"] == 2
        assert data["adj_protects_share"] == 100.0
    else:
        assert data["adjudication"]["rows"] == 2
        assert data["adjudication"]["exploit_share_of_decided"] == 100.0
    assert "100" in scanner.md.read_text()


@pytest.mark.parametrize("change", ["different-pmids", "changed-title", "changed-abstract", "stride"])
def test_equal_candidate_counts_do_not_authorize_labels_for_a_different_cohort(
        scanner, monkeypatch, tmp_path, change):
    shard = scanner.records / "part-000.jsonl.gz"
    _write_shard(shard, scanner.candidates)
    _write_shard(scanner.records / "part-001.jsonl.gz", [scanner.neither])
    path, _rows = _export_completed(scanner, monkeypatch, tmp_path)
    original_counts = json.loads(scanner.json.read_text())["counts"]
    changed = [dict(record) for record in scanner.candidates]
    args = []
    if change == "different-pmids":
        for record in changed:
            record["pmid"] = "9" + record["pmid"]
    elif change in {"changed-title", "changed-abstract"}:
        field = change.removeprefix("changed-")
        changed[0][field] += " Updated evidence."
    else:
        args = ["--stride", "2"]
    _write_shard(shard, changed)
    current = scanner.module.scan(2 if change == "stride" else 1)
    assert all(current["counts"][label] == original_counts[label] for label in scanner.labels)
    if change != "stride":
        assert current["counts"] == original_counts
    _reset_reports(scanner)

    with pytest.raises(SystemExit):
        _run(scanner, monkeypatch, "--adjudication", path, *args)

    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("invalid", [
    "blank-label", "unknown-label", "blank-reason", "duplicate-row",
    "missing-row", "extra-row", "changed-regex-label", "missing-hash-column",
])
def test_invalid_explicit_adjudication_preserves_reports(scanner, monkeypatch, tmp_path, invalid):
    _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)
    path, rows = _export_completed(scanner, monkeypatch, tmp_path)
    fields = CSV_FIELDS
    if invalid == "blank-label":
        rows[0]["adjudicated"] = ""
    elif invalid == "unknown-label":
        rows[0]["adjudicated"] = "certainly-valid"
    elif invalid == "blank-reason":
        rows[0]["reason"] = " \t "
    elif invalid == "duplicate-row":
        rows.append(dict(rows[0]))
    elif invalid == "missing-row":
        rows.pop()
    elif invalid == "extra-row":
        rows.append({**rows[0], "pmid": "unscanned", "record_sha256": "f" * 64})
    elif invalid == "changed-regex-label":
        rows[0]["regex_label"] = next(
            label for label in scanner.labels if label != rows[0]["regex_label"])
    elif invalid == "missing-hash-column":
        fields = [field for field in CSV_FIELDS if field != "cohort_sha256"]
        rows = [{key: value for key, value in row.items() if key in fields} for row in rows]
    _write_csv(path, rows, fields)
    _reset_reports(scanner)

    with pytest.raises(SystemExit):
        _run(scanner, monkeypatch, "--adjudication", path)

    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


def test_legacy_csv_is_rejected_when_explicitly_supplied(scanner, monkeypatch):
    _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)

    with pytest.raises(SystemExit):
        _run(scanner, monkeypatch, "--adjudication", scanner.legacy)

    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("legacy_state", ["modified", "missing"])
def test_offline_render_uses_stored_measurement_without_reading_sources_or_labels(
        scanner, monkeypatch, legacy_state):
    stored = json.loads((REPO / "analysis" / scanner.json.name).read_text())
    # Deliberately noncanonical formatting makes even a same-data rewrite fail.
    original = json.dumps(stored, separators=(", ", ": ")).encode() + b"\n\n"
    scanner.json.write_bytes(original)
    if legacy_state == "modified":
        scanner.legacy.write_text("this is no longer a CSV\n")
    else:
        scanner.legacy.unlink()

    def forbidden(*_args, **_kwargs):
        raise AssertionError("Offline rendering attempted to rebuild the measurement")

    monkeypatch.setattr(scanner.module, "scan", forbidden)
    monkeypatch.setattr(scanner.module, "assemble", forbidden)
    if hasattr(scanner.module, "load_adjudication"):
        monkeypatch.setattr(scanner.module, "load_adjudication", forbidden)
    assert not scanner.records.exists()

    assert _run(scanner, monkeypatch, "--render-only") == 0

    assert scanner.json.read_bytes() == original
    assert scanner.md.read_text() == scanner.module.render(stored)


@pytest.mark.parametrize("option", ["--export-candidates", "--adjudication"])
def test_offline_render_rejects_fresh_adjudication_options(scanner, monkeypatch, tmp_path, option):
    with pytest.raises(SystemExit):
        _run(scanner, monkeypatch, "--render-only", option, tmp_path / "labels.csv")

    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("render_only", [False, True], ids=["fresh", "offline"])
def test_render_failure_preserves_both_reports(scanner, monkeypatch, tmp_path, render_only):
    args = []
    export = tmp_path / "candidates.csv"
    export.write_bytes(b"Existing candidates.\n")
    if render_only:
        scanner.json.write_bytes((REPO / "analysis" / scanner.json.name).read_bytes())
        args = ["--render-only"]
    else:
        _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)
        args = ["--export-candidates", export]
    before = _snapshot(scanner)

    def fail(_data):
        raise RuntimeError("intentional render failure")

    monkeypatch.setattr(scanner.module, "render", fail)
    with pytest.raises(RuntimeError, match="intentional render failure"):
        _run(scanner, monkeypatch, *args)

    assert _snapshot(scanner) == before
    assert export.read_bytes() == b"Existing candidates.\n"


def test_output_directory_keeps_default_reports_untouched(scanner, monkeypatch, tmp_path):
    _write_shard(scanner.records / "part.jsonl.gz", [scanner.neither])
    defaults = tmp_path / "default-reports"
    defaults.mkdir()
    default_json, default_md = defaults / scanner.json.name, defaults / scanner.md.name
    default_json.write_bytes(OLD_JSON)
    default_md.write_bytes(OLD_MD)
    monkeypatch.setattr(scanner.module, "OUT_JSON", default_json)
    monkeypatch.setattr(scanner.module, "OUT_MD", default_md)

    assert _run(scanner, monkeypatch) == 0

    data = json.loads(scanner.json.read_text())
    assert data["total"] == data["counts"]["neither"] == 1
    assert data["adjudication"] == {}
    assert default_json.read_bytes() == OLD_JSON
    assert default_md.read_bytes() == OLD_MD
    assert "44.4%" not in scanner.md.read_text()


def test_complete_but_ambiguous_adjudication_has_no_directional_share(scanner, monkeypatch, tmp_path):
    _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)
    path, rows = _export_completed(scanner, monkeypatch, tmp_path)
    for row in rows:
        row["adjudicated"] = "ambiguous"
    _write_csv(path, rows)

    assert _run(scanner, monkeypatch, "--adjudication", path) == 0

    data = json.loads(scanner.json.read_text())
    assert data["adjudication"]
    if scanner.name == "census_hypoxia_direction":
        assert data["adjudication"]["directional"] == 0
        assert data["adj_protects_share"] is None
    else:
        assert data["adjudication"]["decided"] == 0
        assert data["adjudication"]["exploit_share_of_decided"] is None
    assert "None%" not in scanner.md.read_text()


def test_candidate_export_and_adjudication_are_mutually_exclusive(scanner, monkeypatch, tmp_path):
    with pytest.raises(SystemExit):
        _run(scanner, monkeypatch, "--export-candidates", tmp_path / "candidates.csv",
             "--adjudication", tmp_path / "labels.csv")

    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("link_kind", ["symlink", "hardlink"])
@pytest.mark.parametrize("render_only", [False, True], ids=["fresh", "offline"])
def test_report_outputs_cannot_alias_each_other(scanner, monkeypatch, link_kind, render_only):
    if render_only:
        stored = json.loads((REPO / "analysis" / scanner.json.name).read_text())
        scanner.json.write_bytes(json.dumps(stored, separators=(", ", ": ")).encode() + b"\n\n")
    else:
        _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)
    scanner.md.unlink()
    if link_kind == "symlink":
        scanner.md.symlink_to(scanner.json)
    else:
        scanner.md.hardlink_to(scanner.json)
    before = _snapshot(scanner)

    with pytest.raises(SystemExit):
        _run(scanner, monkeypatch, *(["--render-only"] if render_only else []))

    assert _snapshot(scanner) == before
    assert scanner.md.samefile(scanner.json)


@pytest.mark.parametrize("output_kind", ["json", "md"])
@pytest.mark.parametrize("link_kind", ["same-path", "symlink", "hardlink"])
def test_adjudication_input_cannot_alias_a_report_output(
        scanner, monkeypatch, tmp_path, output_kind, link_kind):
    _write_shard(scanner.records / "part.jsonl.gz", scanner.candidates)
    labels, _rows = _export_completed(scanner, monkeypatch, tmp_path)
    report = getattr(scanner, output_kind)
    report.unlink()
    if link_kind == "same-path":
        report.write_bytes(labels.read_bytes())
        labels = report
    elif link_kind == "symlink":
        report.symlink_to(labels)
    else:
        report.hardlink_to(labels)
    # The aliased report contains a valid completed CSV, so success without the
    # path guard would consume its labels and then overwrite that same input.
    before = _snapshot(scanner)
    labels_before = labels.read_bytes()

    with pytest.raises(SystemExit):
        _run(scanner, monkeypatch, "--adjudication", labels)

    assert _snapshot(scanner) == before
    assert labels.read_bytes() == labels_before
    assert report.samefile(labels)


def test_scan_stride_selects_sorted_shards_and_all_their_records(scanner):
    sizes = {0: 1, 1: 7, 2: 2, 3: 9, 4: 4}
    for i in (3, 0, 4, 1, 2):
        records = [{**scanner.candidates[0], "pmid": f"{i}-{j}"}
                   for j in range(sizes[i])]
        _write_shard(scanner.records / f"part-{i:03}.jsonl.gz", records)

    result = scanner.module.scan(2)

    assert result["total"] == 7
    assert result["counts"][scanner.labels[0]] == 7
    assert result["cohort"]
