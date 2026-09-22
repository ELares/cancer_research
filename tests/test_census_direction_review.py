"""Local review packets preserve complete evidence without inventing decisions."""
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys

import pytest

from test_census_snapshot import write_snapshot


REPO = Path(__file__).resolve().parents[1]
PATHS = {
    "intersection-records.jsonl", "readiness.json", "readiness.md",
    "hypoxia/census-hypoxia-direction.json", "hypoxia/census-hypoxia-direction.md",
    "hypoxia/candidates.csv", "thesis/census-thesis-direction.json",
    "thesis/census-thesis-direction.md", "thesis/candidates.csv",
}
PRIVATE_TEXT = "PRIVATE-EVIDENCE-TEXT must stay in the local packet."


@pytest.fixture
def review(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    import census_direction_review
    return census_direction_review


@pytest.fixture
def root(tmp_path):
    mesh = ["Ferroptosis", "Hypoxia", "Drug Resistance, Neoplasm"]
    titles = ["Hypoxia protects cells", None,
              "Hypoxia protects cells; hypoxia promotes ferroptosis", None,
              "Overcoming drug resistance", "Ferroptosis resistance",
              "Overcoming drug resistance and ferroptosis resistance"]
    records = [{"pmid": str(index), "mesh": mesh, "title": title,
                "abstract": PRIVATE_TEXT, "cancer_basis": "c04"}
               for index, title in enumerate(titles, 1)]
    records[1]["abstract"] = "Hypoxia promotes ferroptosis"
    records[3]["abstract"] = None
    records[4]["abstract"] = None
    records.append({"pmid": "outside", "mesh": [], "title": "Out of scope"})
    return write_snapshot(tmp_path, {"part-a.jsonl.gz": records[:4],
                                     "part-b.jsonl.gz": records[4:]})


@pytest.fixture
def prepared(review, root):
    return review.prepare_bundle(root)


def _save_bundle(path, files):
    for name, content in files.items():
        destination = path / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    return path


def _rewrite_readiness(review, bundle, report):
    (bundle / "readiness.json").write_text(json.dumps(report, indent=1) + "\n")
    (bundle / "readiness.md").write_text(review.render(report))


def test_prepared_packet_covers_both_complete_intersections_without_public_article_text(prepared):
    report, files = prepared
    assert set(files) == PATHS
    assert all(type(content) is bytes for content in files.values())
    assert report["schema_version"] == 1
    assert report["status"] == "unreviewed"
    assert report["snapshot"]["observed_records"] == 8
    assert report["union_records"] == 7
    records = [json.loads(line) for line in files["intersection-records.jsonl"].splitlines()]
    assert {record["pmid"] for record in records} == {str(i) for i in range(1, 8)}
    assert PRIVATE_TEXT in files["intersection-records.jsonl"].decode()
    for name in ("readiness.json", "readiness.md"):
        assert PRIVATE_TEXT not in files[name].decode()
        assert "Hypoxia protects cells" not in files[name].decode()
    for name, labels in (("hypoxia", ("protects", "sensitises")),
                         ("thesis", ("exploit", "obstacle"))):
        summary = report["analyses"][name]
        assert summary["counts"] == {labels[0]: 1, labels[1]: 1, "both": 1, "neither": 4}
        assert summary["intersection_records"] == 7
        assert summary["candidate_records"] == 2
        assert summary["excluded_records"] == 5
        local = json.loads(files[f"{name}/census-{name}-direction.json"])
        assert summary["cohort_sha256"] == local["cohort"]["cohort_sha256"]
    assert set(report["payload_sha256"]) == PATHS - {"readiness.json", "readiness.md"}
    for name, digest in report["payload_sha256"].items():
        assert digest == hashlib.sha256(files[name]).hexdigest()
    for name, digest in report["provenance"]["source_sha256"].items():
        assert digest == hashlib.sha256((REPO / name).read_bytes()).hexdigest()


def test_missing_text_and_blank_worksheets_remain_distinct_from_adjudication(prepared):
    report, files = prepared
    assert report["analyses"]["hypoxia"]["missing_text"]["sensitises"] == {
        "title": 1, "abstract": 0, "both": 0,
    }
    assert report["analyses"]["thesis"]["missing_text"]["exploit"] == {
        "title": 0, "abstract": 1, "both": 0,
    }
    for name in ("hypoxia", "thesis"):
        assert report["analyses"][name]["missing_text"]["neither"]["both"] == 1
        rows = list(csv.DictReader(io.StringIO(files[f"{name}/candidates.csv"].decode())))
        assert len(rows) == 2
        assert all(row["adjudicated"] == row["reason"] == "" for row in rows)
        local = json.loads(files[f"{name}/census-{name}-direction.json"])
        assert local["adjudication"] == {}
        assert {row["record_sha256"] for row in rows} == {
            candidate["record_sha256"] for candidate in local["cohort"]["candidates"]
        }


def test_nonempty_census_with_no_eligible_intersection_produces_empty_review_worksheets(review, tmp_path):
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [{"pmid": "1", "mesh": []}]})
    report, files = review.prepare_bundle(root)
    assert report["union_records"] == 0
    assert files["intersection-records.jsonl"] == b""
    for name, summary in report["analyses"].items():
        assert summary["intersection_records"] == summary["candidate_records"] == 0
        assert sum(summary["counts"].values()) == 0
        assert list(csv.DictReader(io.StringIO(files[f"{name}/candidates.csv"].decode()))) == []
    assert review.verify_bundle(_save_bundle(tmp_path / "bundle", files)) == report


def test_packet_verifies_offline_after_the_original_census_is_unavailable(review, root, prepared, tmp_path):
    report, files = prepared
    bundle = _save_bundle(tmp_path / "bundle", files)
    shutil.rmtree(root)
    assert review.verify_bundle(bundle) == report
    wording = review.render(report)
    assert "does not reread the full census" in " ".join(wording.split())


@pytest.mark.parametrize("separator", ["\u0085", "\u2028", "\u2029"],
                         ids=["next-line", "line-separator", "paragraph-separator"])
def test_unicode_line_separators_in_article_text_are_not_jsonl_record_boundaries(
        review, tmp_path, separator):
    title = f"Hypoxia protects cells{separator}across conditions"
    abstract = f"First paragraph.{separator}Second paragraph."
    records = [
        {"pmid": "1", "mesh": ["Ferroptosis", "Hypoxia"],
         "title": title, "abstract": abstract},
        {"pmid": "2", "mesh": ["Ferroptosis", "Drug Resistance, Neoplasm"],
         "title": "Overcoming drug resistance", "abstract": abstract},
    ]
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": records})
    report, files = review.prepare_bundle(root)
    assert separator.encode("utf-8") in files["intersection-records.jsonl"]
    assert report["union_records"] == 2
    shutil.rmtree(root)
    assert review.verify_bundle(_save_bundle(tmp_path / "bundle", files)) == report
    rows = list(csv.DictReader(io.StringIO(files["hypoxia/candidates.csv"].decode())))
    assert len(rows) == 1
    assert rows[0]["title"] == title
    assert rows[0]["abstract"] == abstract


@pytest.mark.parametrize("name", ["intersection-records.jsonl", "hypoxia/candidates.csv",
                                  "thesis/census-thesis-direction.json",
                                  "hypoxia/census-hypoxia-direction.md", "readiness.md"])
def test_changed_payload_bytes_are_rejected(review, prepared, tmp_path, name):
    _, files = prepared
    bundle = _save_bundle(tmp_path / "bundle", files)
    path = bundle / name
    path.write_bytes(path.read_bytes() + b"changed")
    with pytest.raises(SystemExit):
        review.verify_bundle(bundle)


@pytest.mark.parametrize("change", ["source-hash", "summary", "snapshot-count",
                                   "snapshot-inventory", "snapshot-shard-hash", "decision", "evidence"])
def test_refreshed_hashes_cannot_hide_inconsistent_evidence_or_invented_decisions(
        review, prepared, tmp_path, change):
    report, files = prepared
    bundle = _save_bundle(tmp_path / "bundle", files)
    if change == "source-hash":
        name = next(iter(report["provenance"]["source_sha256"]))
        report["provenance"]["source_sha256"][name] = "0" * 64
    elif change == "summary":
        report["analyses"]["hypoxia"]["counts"]["protects"] += 1
    elif change == "snapshot-count":
        report["snapshot"]["observed_records"] += 1
    elif change == "snapshot-inventory":
        report["snapshot"]["inventory_sha256"] = "0" * 64
    elif change == "snapshot-shard-hash":
        report["snapshot"]["shards"][0]["sha256"] = "0" * 64
    else:
        name = "hypoxia/candidates.csv" if change == "decision" else "intersection-records.jsonl"
        if change == "decision":
            rows = list(csv.DictReader(io.StringIO(files[name].decode())))
            rows[0].update(adjudicated="protects", reason="Invented decision")
            output = io.StringIO(newline="")
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
            content = output.getvalue().encode()
        else:
            records = [json.loads(line) for line in files[name].splitlines()]
            records[0]["abstract"] = "Changed evidence without changing its lexical direction."
            content = "".join(json.dumps(record) + "\n" for record in records).encode()
        (bundle / name).write_bytes(content)
        report["payload_sha256"][name] = hashlib.sha256(content).hexdigest()
    _rewrite_readiness(review, bundle, report)
    with pytest.raises(SystemExit):
        review.verify_bundle(bundle)


@pytest.mark.parametrize("field,value", [("mesh", "Ferroptosis"), ("mesh", ["Ferroptosis", 7]),
                                        ("title", 7), ("abstract", [])])
def test_malformed_record_fields_are_rejected_before_scientific_filtering(review, tmp_path, field, value):
    record = {"pmid": "1", "mesh": [], "title": "Outside the selected intersection"}
    record[field] = value
    root = write_snapshot(tmp_path, {"part-a.jsonl.gz": [record]})
    with pytest.raises(SystemExit):
        review.prepare_bundle(root)


@pytest.mark.parametrize("fault", ["missing", "extra"])
def test_packet_requires_the_complete_fixed_file_inventory(review, prepared, tmp_path, fault):
    _, files = prepared
    bundle = _save_bundle(tmp_path / "bundle", files)
    if fault == "missing":
        (bundle / "thesis/candidates.csv").unlink()
    else:
        (bundle / "extra.txt").write_text("unlisted evidence")
    with pytest.raises(SystemExit):
        review.verify_bundle(bundle)


def test_publishing_does_not_replace_an_existing_packet(review, prepared, tmp_path):
    report, files = prepared
    bundle = tmp_path / "bundle"
    review.publish_bundle(bundle, files)
    assert review.verify_bundle(bundle) == report
    before = {name: (bundle / name).read_bytes() for name in PATHS}
    with pytest.raises(SystemExit):
        review.publish_bundle(bundle, files)
    assert {name: (bundle / name).read_bytes() for name in PATHS} == before


@pytest.mark.parametrize("via_symlink", [False, True])
def test_cli_rejects_article_text_packet_inside_repository_before_scanning(
        review, root, tmp_path, monkeypatch, via_symlink):
    repository = tmp_path / "repository"
    repository.mkdir()
    monkeypatch.setattr(review, "REPO", repository)
    parent = repository
    if via_symlink:
        parent = tmp_path / "alias"
        parent.symlink_to(repository, target_is_directory=True)
    bundle = parent / "review-packet"

    def forbidden(*args):
        pytest.fail("invalid packet destination must fail before scanning")

    monkeypatch.setattr(review, "prepare_bundle", forbidden)
    monkeypatch.setattr(sys, "argv", ["census_direction_review.py", "--atlas-root", str(root),
                                     "--bundle-dir", str(bundle)])
    with pytest.raises(SystemExit, match="outside the repository"):
        review.main()
    assert not bundle.exists()


@pytest.mark.parametrize("failure", ["late-json", "render"])
def test_cli_preparation_failures_preserve_public_outputs_and_leave_no_packet(
        review, root, tmp_path, monkeypatch, failure):
    out_json, out_md = tmp_path / "public.json", tmp_path / "public.md"
    out_json.write_bytes(b"existing JSON")
    out_md.write_bytes(b"existing Markdown")
    monkeypatch.setattr(review, "OUT_JSON", out_json)
    monkeypatch.setattr(review, "OUT_MD", out_md)
    if failure == "late-json":
        import gzip
        (root / "records/part-b.jsonl.gz").write_bytes(gzip.compress(b"{bad JSON\n"))
    else:
        def fail_render(report):
            raise ValueError("render failed")
        monkeypatch.setattr(review, "render", fail_render)
    bundle = tmp_path / "bundle"
    monkeypatch.setattr(sys, "argv", ["census_direction_review.py", "--atlas-root", str(root),
                                     "--bundle-dir", str(bundle)])
    with pytest.raises(SystemExit):
        review.main()
    assert out_json.read_bytes() == b"existing JSON"
    assert out_md.read_bytes() == b"existing Markdown"
    assert not bundle.exists()


def test_render_only_uses_the_public_report_without_reading_census_or_writing_json(
        review, prepared, tmp_path, monkeypatch):
    report, files = prepared
    out_json, out_md = tmp_path / "public.json", tmp_path / "public.md"
    original = json.dumps(report, separators=(",", ":")).encode()
    out_json.write_bytes(original)
    monkeypatch.setattr(review, "OUT_JSON", out_json)
    monkeypatch.setattr(review, "OUT_MD", out_md)

    def forbidden(*args):
        pytest.fail("offline rendering must not read the census")

    monkeypatch.setattr(review, "read_snapshot", forbidden)
    monkeypatch.setattr(sys, "argv", ["census_direction_review.py", "--render-only"])
    assert review.main() == 0
    assert out_json.read_bytes() == original
    assert out_md.read_bytes() == files["readiness.md"]


def test_public_rendering_does_not_depend_on_dictionary_insertion_order(review, prepared):
    report, _ = prepared

    def reverse_dicts(value):
        if isinstance(value, dict):
            return {key: reverse_dicts(item) for key, item in reversed(list(value.items()))}
        if isinstance(value, list):
            return [reverse_dicts(item) for item in value]
        return value

    assert review.render(reverse_dicts(report)) == review.render(report)
