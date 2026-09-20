"""A missing census must never replace published counts with empty reports.

Exercise each real CLI, scanner and gzip parser against isolated tiny inputs.
The valid zero-match cases distinguish unavailable input from a scientific
result of zero, especially where a scanner counts only ferroptosis articles.
"""
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent
SCANNERS = (
    "census_mechanism_growth",
    "census_unnamed_modalities",
    "census_protocol_precedent",
    "census_modality_comparison",
    "census_normal_tissue",
    "census_mechanism_sites",
)
UNMATCHED = {"pmid": "42", "mesh": ["Unmapped test descriptor"]}
EXISTING_JSON = b'{"existing": "published counts"}\n'
EXISTING_MD = b"Existing published interpretation.\n"


def _load(name, monkeypatch):
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        f"{name}_input_test", REPO / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _isolate_reports(module, tmp_path, monkeypatch):
    for attr, content in (("OUT_JSON", EXISTING_JSON), ("OUT_MD", EXISTING_MD)):
        path = tmp_path / attr
        path.write_bytes(content)
        monkeypatch.setattr(module, attr, path)


def _assert_preserved(module):
    assert module.OUT_JSON.read_bytes() == EXISTING_JSON
    assert module.OUT_MD.read_bytes() == EXISTING_MD


def _write_shard(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for record in records:
            fh.write(json.dumps(record) + "\n")


@pytest.fixture(params=SCANNERS)
def scanner_name(request):
    return request.param


@pytest.fixture
def scanner(scanner_name, tmp_path, monkeypatch):
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(tmp_path / "atlas"))
    module = _load(scanner_name, monkeypatch)
    monkeypatch.setattr(module, "RECORDS", tmp_path / "atlas" / "records")
    _isolate_reports(module, tmp_path, monkeypatch)
    monkeypatch.setattr(sys, "argv", [f"{scanner_name}.py"])
    return module


@pytest.mark.parametrize("input_state", ["missing", "empty-directory", "empty-shard"])
def test_unavailable_census_preserves_reports(scanner, input_state):
    records = scanner.RECORDS
    if input_state != "missing":
        records.mkdir(parents=True)
    if input_state == "empty-shard":
        _write_shard(records / "part.jsonl.gz", [])

    with pytest.raises(SystemExit) as exc:
        scanner.main()

    message = str(exc.value)
    assert "No census records" in message
    assert str(records) in message
    assert "atlas_baseline.py" in message
    assert "--render-only" in message
    _assert_preserved(scanner)


@pytest.mark.parametrize("stride", [0, -1])
def test_invalid_stride_preserves_reports(scanner, scanner_name, monkeypatch, stride):
    _write_shard(scanner.RECORDS / "part.jsonl.gz", [UNMATCHED])
    monkeypatch.setattr(sys, "argv", [f"{scanner_name}.py", "--stride", str(stride)])

    with pytest.raises(SystemExit, match="positive integer"):
        scanner.main()

    _assert_preserved(scanner)


def test_external_census_with_no_scientific_matches_is_valid(
        scanner_name, tmp_path, monkeypatch):
    """Use the public root setting without patching RECORDS or scan()."""
    root = tmp_path / "external-atlas"
    records = root / "records"
    _write_shard(records / "part.jsonl.gz", [UNMATCHED])
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(root))
    module = _load(scanner_name, monkeypatch)
    # Check before invoking the CLI so a broken path cannot read real bulk data.
    assert module.RECORDS == records
    _isolate_reports(module, tmp_path, monkeypatch)
    monkeypatch.setattr(sys, "argv", [f"{scanner_name}.py"])

    assert module.main() == 0

    result = json.loads(module.OUT_JSON.read_text())
    if scanner_name in {"census_protocol_precedent", "census_normal_tissue"}:
        assert result["ferroptosis_articles"] == 0
    else:
        assert result["census"] == 1
    if scanner_name == "census_protocol_precedent":
        assert result["arm1_precedent"] == result["arm2_precedent"] == 0
        assert result["thin_arm"] is None
        assert "precedent cannot be compared" in module.OUT_MD.read_text()
        assert "abundant precedent" not in module.OUT_MD.read_text()
    elif scanner_name == "census_normal_tissue":
        assert result["organ_toxicity_articles"] == 0
        assert result["adjudicated_n"] == 0
        assert "verdicts" not in result and "harm_share" not in result
        assert "no normal-tissue harm comparison" in module.OUT_MD.read_text()
    elif scanner_name == "census_modality_comparison":
        assert result["arms"] == {"mesh": {}, "text": {}}
        assert "cannot be compared" in module.OUT_MD.read_text()
        assert "radioligand therapy is further" not in module.OUT_MD.read_text()
    else:
        assert result["rows"] == []
        if scanner_name == "census_unnamed_modalities":
            assert "No trial shares" in module.OUT_MD.read_text()
            assert "hundreds of indexed trials" not in module.OUT_MD.read_text()
    assert module.OUT_MD.read_text().startswith("# ")


def test_render_only_works_without_raw_records(scanner, scanner_name, monkeypatch):
    committed = REPO / "analysis" / f"{scanner_name.replace('_', '-')}.json"
    scanner.OUT_JSON.write_bytes(committed.read_bytes())
    monkeypatch.setattr(sys, "argv", [f"{scanner_name}.py", "--render-only"])
    assert not scanner.RECORDS.exists()

    assert scanner.main() == 0

    result = json.loads(scanner.OUT_JSON.read_text())
    assert result == json.loads(committed.read_text())
    assert scanner.OUT_MD.read_text() == scanner.render(result)
    assert scanner.OUT_MD.read_bytes() != EXISTING_MD


@pytest.mark.parametrize("render_only", [False, True], ids=["scan", "render-only"])
def test_render_failure_preserves_both_reports(
        scanner, scanner_name, monkeypatch, render_only):
    committed = REPO / "analysis" / f"{scanner_name.replace('_', '-')}.json"
    raw = json.loads(committed.read_text())

    def fail_render(_data):
        raise RuntimeError("intentional rendering failure")

    monkeypatch.setattr(scanner, "render", fail_render)
    if render_only:
        # Distinct formatting detects a JSON write even if its counts agree.
        before_json = json.dumps(raw, separators=(",", ":")).encode() + b"\n"
        scanner.OUT_JSON.write_bytes(before_json)
        monkeypatch.setattr(sys, "argv", [f"{scanner_name}.py", "--render-only"])
    else:
        before_json = EXISTING_JSON
        monkeypatch.setattr(scanner, "scan", lambda stride: raw)

    with pytest.raises(RuntimeError, match="intentional rendering failure"):
        scanner.main()

    assert scanner.OUT_JSON.read_bytes() == before_json
    assert scanner.OUT_MD.read_bytes() == EXISTING_MD


def test_parser_failure_after_valid_record_preserves_reports(scanner):
    records = scanner.RECORDS
    records.mkdir(parents=True)
    with gzip.open(records / "part.jsonl.gz", "wt", encoding="utf-8") as fh:
        fh.write(json.dumps(UNMATCHED) + "\n")
        fh.write("{invalid json\n")

    with pytest.raises(json.JSONDecodeError):
        scanner.main()

    _assert_preserved(scanner)


@pytest.mark.parametrize("scanner_name", ["census_mechanism_growth"], indirect=True)
def test_growth_with_dated_input_and_no_mechanisms_has_no_ratio(scanner):
    """A measurable field trend cannot supply a missing mechanism denominator."""
    _write_shard(scanner.RECORDS / "part.jsonl.gz", [
        {**UNMATCHED, "pmid": str(year), "year": year}
        for year in range(2015, 2026)
    ])

    assert scanner.main() == 0

    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["census"] == 11
    assert result["field_start"] == result["field_end"] == 1
    assert result["field_growth"] == 1.0
    assert result["union_growth"] is None
    assert result["mechanisms_over_field"] is None
    assert result["rows"] == []
    assert "growth ratios are unavailable" in scanner.OUT_MD.read_text()


@pytest.mark.parametrize("scanner_name", ["census_mechanism_growth"], indirect=True)
@pytest.mark.parametrize("observations", [
    [(2025, True)],
    [(2015, False), (2025, True)],
    [(2015, True)],
    [(2020, True)],
], ids=["missing-start-year", "zero-mechanism-start", "zero-field-end", "neither-endpoint"])
def test_growth_unavailable_baselines_do_not_produce_comparative_claims(
        scanner, observations):
    _write_shard(scanner.RECORDS / "part.jsonl.gz", [
        {"pmid": str(i), "year": year,
         "mesh": ["Gene Editing"] if matched else UNMATCHED["mesh"]}
        for i, (year, matched) in enumerate(observations)
    ])

    assert scanner.main() == 0

    result = json.loads(scanner.OUT_JSON.read_text())
    markdown = scanner.OUT_MD.read_text()
    assert result["rows"], "Exercise matched input rather than the no-matches guard."
    assert result["mechanisms_over_field"] is None
    assert "growth comparison is unavailable" in markdown
    assert "None" not in markdown
    assert "The mechanisms this project tracks grew" not in markdown
    assert "field is close to flat" not in markdown
    assert "No mechanism is smaller" not in markdown


@pytest.mark.parametrize("scanner_name", ["census_mechanism_growth"], indirect=True)
@pytest.mark.parametrize("field_counts,mechanism_counts,ratio,direction", [
    ((100, 200), (40, 20), 0.25, "lower than"),
    ((100, 200), (30, 60), 1.0, "equal to"),
    ((100, 200), (30, 90), 1.5, "higher than"),
    ((200, 100), (100, 80), 1.6, "higher than"),
    ((100, 100), (30, 30), 1.0, "equal to"),
    ((100, 100), (30, 0), 0.0, "lower than"),
], ids=["mechanisms-decline-field-grows", "equal-growth", "faster-growth",
        "both-decline", "both-constant", "observed-zero-mechanisms-at-end"])
def test_growth_comparison_follows_observed_factors(
        scanner, field_counts, mechanism_counts, ratio, direction):
    records = []
    for year, field, mechanism in zip((2015, 2025), field_counts, mechanism_counts):
        records.extend(
            {"pmid": f"{year}-{i}", "year": year,
             "mesh": ["Gene Editing"] if i < mechanism else UNMATCHED["mesh"]}
            for i in range(field)
        )
    _write_shard(scanner.RECORDS / "part.jsonl.gz", records)

    assert scanner.main() == 0

    result = json.loads(scanner.OUT_JSON.read_text())
    markdown = scanner.OUT_MD.read_text()
    assert result["mechanisms_over_field"] == ratio
    assert f"**x{ratio}** the field growth factor, {direction}" in markdown
    assert "does not by itself establish" in markdown
    assert "growth comparison is unavailable" not in markdown
    assert "field is close to flat" not in markdown
    assert "The mechanisms this project tracks grew" not in markdown
    if not mechanism_counts[1]:
        assert result["union_growth"] == 0.0
        assert "**x0.0**" in markdown


@pytest.mark.parametrize("scanner_name", ["census_mechanism_growth"], indirect=True)
@pytest.mark.parametrize("assessable", [False, True])
def test_growth_absence_of_assessed_declines_does_not_include_thin_baselines(
        scanner, assessable):
    mechanisms = {"thin-declining": {"2015": 20, "2025": 10}}
    if assessable:
        mechanisms["stable"] = {"2015": 30, "2025": 30}
    result = scanner.assemble({
        "census": 200, "shards": 1, "start_year": 2015, "end_year": 2025,
        "field_by_year": {"2015": 100, "2025": 100},
        "union_by_year": {"2015": 50 if assessable else 20,
                          "2025": 40 if assessable else 10},
        "mechanism_by_year": mechanisms,
    })

    markdown = scanner.render(result)

    assert "No mechanism is smaller" not in markdown
    if assessable:
        assert "Among mechanisms meeting the 30-article starting baseline" in markdown
    else:
        assert "No mechanism meets the 30-article starting baseline" in markdown


@pytest.mark.parametrize("scanner_name", ["census_mechanism_growth"], indirect=True)
@pytest.mark.parametrize("field,mechanisms,expected", [
    ((100000, 499), (10000, 100), 2.0),
    ((133204, 147008), (6050, 16776), 2.51),
], ids=["small-positive-field-factor", "round-final-ratio-only"])
def test_growth_comparison_uses_raw_endpoints(scanner, field, mechanisms, expected):
    result = scanner.assemble({
        "census": sum(field), "shards": 1, "start_year": 2015, "end_year": 2025,
        "field_by_year": dict(zip(("2015", "2025"), field)),
        "union_by_year": dict(zip(("2015", "2025"), mechanisms)),
        "mechanism_by_year": {"ferroptosis": dict(zip(("2015", "2025"), mechanisms))},
    })

    assert result["mechanisms_over_field"] == expected
    markdown = scanner.render(result)
    assert "growth comparison is unavailable" not in markdown
    if field[1] < field[0] / 200:
        assert "x0.00499" in markdown


@pytest.mark.parametrize("scanner_name", ["census_mechanism_growth"], indirect=True)
def test_growth_detects_small_declines_hidden_by_display_rounding(scanner):
    result = scanner.assemble({
        "census": 120000, "shards": 1, "start_year": 2015, "end_year": 2025,
        "field_by_year": {"2015": 40000, "2020": 40000, "2025": 40000},
        "union_by_year": {"2015": 20029, "2020": 20019, "2025": 20017},
        "mechanism_by_year": {
            "slightly-declining": {"2015": 20000, "2020": 20000, "2025": 19999},
            "below-both-floors": {"2015": 29, "2020": 19, "2025": 18},
        },
    })

    assert result["shrinking"] == ["slightly-declining"]
    assert result["recent_shrinking"] == ["slightly-declining"]
    markdown = scanner.render(result)
    assert "1 mechanism(s) are SMALLER" in markdown
    assert "No mechanism is smaller" not in markdown


@pytest.mark.parametrize("scanner_name", ["census_normal_tissue"], indirect=True)
def test_normal_tissue_reassembly_does_not_reuse_another_samples_labels(scanner):
    committed = REPO / "analysis/census-normal-tissue.json"
    raw = json.loads(committed.read_text())
    raw["sample"] = [{"pmid": "unadjudicated-new-input"}]

    result = scanner.assemble(raw)

    assert result["organ_toxicity_articles"] == raw["organ_toxicity_articles"]
    assert result["adjudicated_n"] == 0
    for key in ("verdicts", "harm_share", "inhibitor_named_in_sample",
                "inhibitor_is_probe", "probe_share_of_inhibitor_mentions"):
        assert key not in result
    assert "No adjudications match" in scanner.render(result)


@pytest.mark.parametrize("scanner_name", ["census_protocol_precedent"], indirect=True)
@pytest.mark.parametrize("title,thin_role", [
    ("RSL3", "FSP1 / DHODH inhibitor (arm 2)"),
    ("iFSP1", "GPX4 inhibitor (arm 1)"),
    ("RSL3 and iFSP1", None),
])
def test_protocol_thin_arm_uses_counts_including_zero_and_ties(scanner, title, thin_role):
    _write_shard(scanner.RECORDS / "part.jsonl.gz", [
        {"pmid": "42", "mesh": ["Ferroptosis"], "title": title},
    ])

    assert scanner.main() == 0

    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["ferroptosis_articles"] == 1
    assert result["thin_arm"] == thin_role
    report = scanner.OUT_MD.read_text()
    assert "factor of None" not in report
    if thin_role is None:
        assert "neither arm has less precedent" in report


@pytest.mark.parametrize("scanner_name", ["census_mechanism_sites"], indirect=True)
def test_assigned_site_without_class_matches_has_no_comparison(scanner):
    _write_shard(scanner.RECORDS / "part.jsonl.gz", [
        {"pmid": "42", "mesh": ["Breast Neoplasms"]},
    ])

    assert scanner.main() == 0

    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["census"] == 1
    assert result["site_totals"] == {"breast": 1}
    assert result["physical_total"] == result["pharmacological_total"] == 0
    assert [row["site"] for row in result["rows"]] == ["breast"]
    assert "comparison of class enrichment is unavailable" in scanner.OUT_MD.read_text()


@pytest.mark.parametrize("scanner_name", ["census_mechanism_sites"], indirect=True)
def test_site_with_zero_physical_count_keeps_class_comparison(scanner):
    """Both classes are measurable even when a site's physical count is zero."""
    _write_shard(scanner.RECORDS / "part.jsonl.gz", [
        {"pmid": "41", "mesh": ["Breast Neoplasms",
                                   "High-Intensity Focused Ultrasound Ablation"]},
        {"pmid": "42", "mesh": ["Lung Neoplasms", "Immune Checkpoint Inhibitors"]},
    ])

    assert scanner.main() == 0

    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["census"] == 2
    assert result["site_totals"] == {"breast": 1, "lung": 1}
    assert result["physical_total"] == result["pharmacological_total"] == 1
    rows = {row["site"]: row for row in result["rows"]}
    assert rows["breast"]["physical_enrichment"] == 2.0
    assert rows["lung"]["physical_enrichment"] == 0.0
    assert rows["lung"]["pharmacological_enrichment"] == 2.0
    report = scanner.OUT_MD.read_text()
    assert "| breast | 1 | 2.00x | 0 | 0.00x |" in report
    assert "| lung | 0 | 0.00x | 1 | 2.00x |" in report
    assert "ratio of highest to lowest physical enrichment is undefined" in report


def test_helper_samples_sorted_shards_and_keeps_all_records(tmp_path, monkeypatch):
    module = _load("census_input", monkeypatch)
    # Deliberately create out of order, with multiple records in each shard.
    for name, pmids in (("part-c", ["c1", "c2"]), ("part-a", ["a1", "a2"]),
                        ("part-b", ["b1", "b2"])):
        _write_shard(tmp_path / f"{name}.jsonl.gz", [{"pmid": p} for p in pmids])
    (tmp_path / "ignored.jsonl").write_text("invalid uncompressed input\n")

    result = list(module.iter_census_records(tmp_path, stride=2))

    assert [r["pmid"] for r in result] == ["a1", "a2", "c1", "c2"]


def test_helper_accepts_empty_shards_before_real_records(tmp_path, monkeypatch):
    module = _load("census_input", monkeypatch)
    _write_shard(tmp_path / "part-a.jsonl.gz", [])
    _write_shard(tmp_path / "part-b.jsonl.gz", [UNMATCHED])

    assert list(module.iter_census_records(tmp_path)) == [UNMATCHED]


def test_helper_rejects_empty_sample_even_when_unsampled_records_exist(
        tmp_path, monkeypatch):
    module = _load("census_input", monkeypatch)
    _write_shard(tmp_path / "part-a.jsonl.gz", [])
    _write_shard(tmp_path / "part-b.jsonl.gz", [UNMATCHED])

    with pytest.raises(SystemExit, match="No census records"):
        list(module.iter_census_records(tmp_path, stride=2))


def test_helper_default_root_follows_environment(tmp_path, monkeypatch):
    module = _load("census_input", monkeypatch)
    root = tmp_path / "external-atlas"
    _write_shard(root / "records" / "part.jsonl.gz", [UNMATCHED])
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(root))

    assert list(module.iter_census_records()) == [UNMATCHED]
