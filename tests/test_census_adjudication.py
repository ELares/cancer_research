"""Identity and coverage checks for new, fully adjudicated direction cohorts."""
import csv
import importlib.util
import io
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "census_adjudication_test", REPO / "scripts/census_adjudication.py")
scanner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(scanner)


def _cohort(records, analysis="example"):
    cohort = scanner.AdjudicationCohort(analysis, ("exploit", "obstacle"))
    for record, label in records:
        cohort.add(record, label)
    return cohort.as_dict()


@pytest.fixture
def cohort():
    return _cohort([
        ({"pmid": "123", "title": 'A title with "quotes", commas and β',
          "abstract": "First line.\nSecond line.", "mesh": ["Ferroptosis"]}, "exploit"),
        ({"pmid": "456", "title": "Another title", "year": 2024}, "obstacle"),
    ])


def _completed(cohort):
    rows = list(csv.DictReader(io.StringIO(scanner.candidate_csv(cohort))))
    for row in rows:
        row.update(adjudicated=row["regex_label"], reason="Checked title and abstract.")
    return rows


def _save(path, rows, fields=scanner.CSV_FIELDS):
    with path.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_export_requires_explicit_decisions_and_preserves_full_text(cohort, tmp_path):
    exported = list(csv.DictReader(io.StringIO(scanner.candidate_csv(cohort))))
    assert len(exported) == 2
    assert all(row["adjudicated"] == row["reason"] == "" for row in exported)
    assert exported[0]["abstract"] == "First line.\nSecond line."
    path = _save(tmp_path / "complete.csv", _completed(cohort))
    assert scanner.read_adjudication(path, cohort, ("exploit", "obstacle", "ambiguous")) == _completed(cohort)


def test_cohort_digest_ignores_input_order_but_includes_unclassified_members():
    rows = [({"pmid": "123", "title": "A"}, "exploit"),
            ({"pmid": "456", "title": "B"}, "neither")]
    original = _cohort(rows)
    assert original["cohort_sha256"] == _cohort(list(reversed(rows)))["cohort_sha256"]
    assert original["cohort_sha256"] != _cohort(rows[:1])["cohort_sha256"]
    assert original["records"] == 2
    assert len(original["candidates"]) == 1


@pytest.mark.parametrize("change", ["title", "abstract", "mesh", "year", "pmid", "label", "analysis"])
def test_record_or_analysis_changes_invalidate_the_binding(change):
    record = {"pmid": "123", "title": "A", "abstract": "B", "mesh": ["C"], "year": 2024}
    original = _cohort([(record, "exploit")])
    altered = dict(record)
    label, analysis = "exploit", "example"
    if change == "label":
        label = "obstacle"
    elif change == "analysis":
        analysis = "another-analysis"
    else:
        altered[change] = {"mesh": ["D"], "year": 2025}.get(change, "changed")
    changed = _cohort([(altered, label)], analysis)
    assert changed["cohort_sha256"] != original["cohort_sha256"]
    if change not in ("label", "analysis"):
        assert changed["candidates"][0]["record_sha256"] != original["candidates"][0]["record_sha256"]


@pytest.mark.parametrize("second", [
    {"pmid": "123", "title": "A"},
    {"pmid": "123", "title": "Changed"},
    {"pmid": "00123", "title": "Changed"},
])
def test_duplicate_observations_are_rejected(second):
    with pytest.raises(SystemExit, match="Duplicate record"):
        _cohort([({"pmid": "123", "title": "A"}, "exploit"), (second, "neither")])


def test_identical_content_without_pmids_is_not_counted_twice():
    with pytest.raises(SystemExit, match="Duplicate record"):
        _cohort([({"title": "A"}, "exploit"), ({"title": "A"}, "exploit")])


@pytest.mark.parametrize("field", scanner.IDENTITY_FIELDS)
def test_editing_exported_identity_or_displayed_evidence_is_rejected(cohort, tmp_path, field):
    rows = _completed(cohort)
    rows[0][field] = "modified"
    path = _save(tmp_path / "changed.csv", rows)
    with pytest.raises(SystemExit, match="not a selected candidate|changes the exported"):
        scanner.read_adjudication(path, cohort, ("exploit", "obstacle", "ambiguous"))


def test_reordered_decisions_remain_valid(cohort, tmp_path):
    rows = list(reversed(_completed(cohort)))
    path = _save(tmp_path / "reordered.csv", rows)
    assert scanner.read_adjudication(path, cohort, ("exploit", "obstacle")) == rows


@pytest.mark.parametrize("contents", [
    "title,regex_label,adjudicated,reason\nA,exploit,exploit,Checked\n",
    ",".join(scanner.CSV_FIELDS[:-1] + ("title",)) + "\n",
    ",".join(scanner.CSV_FIELDS + ("extra",)) + "\n",
    ",".join(scanner.CSV_FIELDS) + '\n"unterminated',
])
def test_legacy_or_malformed_headers_and_csv_are_rejected(cohort, tmp_path, contents):
    path = tmp_path / "bad.csv"
    path.write_text(contents)
    with pytest.raises(SystemExit, match="Adjudication"):
        scanner.read_adjudication(path, cohort, ("exploit", "obstacle"))


def test_no_candidates_exports_a_header_but_cannot_produce_an_adjudicated_share(tmp_path):
    empty = _cohort([({"pmid": "123", "title": "A"}, "neither")])
    assert list(csv.DictReader(io.StringIO(scanner.candidate_csv(empty)))) == []
    path = _save(tmp_path / "empty.csv", [])
    with pytest.raises(SystemExit, match="no singly classified candidates"):
        scanner.read_adjudication(path, empty, ("exploit", "obstacle"))


@pytest.mark.parametrize("field", ["title", "abstract"])
@pytest.mark.parametrize("newline", ["\r", "\n", "\r\n"])
def test_exported_text_with_newlines_round_trips(tmp_path, field, newline):
    record = {"pmid": "123", "title": "Title", "abstract": "Abstract"}
    record[field] = f"First{newline}Second"
    current = _cohort([(record, "exploit")])
    rows = _completed(current)
    assert len(rows) == 1
    assert rows[0][field] == record[field]
    path = _save(tmp_path / "newlines.csv", rows)
    assert scanner.read_adjudication(path, current, ("exploit", "obstacle")) == rows


@pytest.mark.parametrize("field", ["title", "abstract", "reason"])
def test_large_text_import_restores_csv_field_limit(tmp_path, field):
    text = "Evidence β " * 15000
    record = {"pmid": "123", "title": "Title", "abstract": "Abstract"}
    if field != "reason":
        record[field] = text
    current = _cohort([(record, "exploit")])
    row = {**current["candidates"][0], "cohort_sha256": current["cohort_sha256"],
           "adjudicated": "exploit", "reason": "Reviewed complete evidence."}
    if field == "reason":
        row[field] = text
    else:
        assert text in scanner.candidate_csv(current)
    path = _save(tmp_path / "long.csv", [row])
    previous_limit = csv.field_size_limit()
    assert scanner.read_adjudication(path, current, ("exploit", "obstacle")) == [row]
    assert csv.field_size_limit() == previous_limit


@pytest.mark.parametrize("failure", ["header", "parsing", "validation"])
def test_failed_import_restores_csv_field_limit(cohort, tmp_path, failure):
    rows = _completed(cohort)
    rows[0]["reason"] = "Long reviewed reason. " * 10000
    path = _save(tmp_path / "invalid.csv", rows)
    if failure == "header":
        path.write_text(path.read_text().replace("cohort_sha256", "unknown_column", 1))
    elif failure == "parsing":
        with path.open("a") as target:
            target.write('"unterminated')
    else:
        rows[0]["adjudicated"] = "invalid"
        _save(path, rows)
    previous_limit = csv.field_size_limit()
    with pytest.raises(SystemExit, match="Adjudication"):
        scanner.read_adjudication(path, cohort, ("exploit", "obstacle"))
    assert csv.field_size_limit() == previous_limit


def test_import_preserves_an_existing_higher_csv_field_limit(cohort, tmp_path):
    rows = _completed(cohort)
    path = _save(tmp_path / "complete.csv", rows)
    previous_limit = csv.field_size_limit()
    higher_limit = max(previous_limit, path.stat().st_size) + 1000
    try:
        csv.field_size_limit(higher_limit)
        assert scanner.read_adjudication(path, cohort, ("exploit", "obstacle")) == rows
        assert csv.field_size_limit() == higher_limit
    finally:
        csv.field_size_limit(previous_limit)
