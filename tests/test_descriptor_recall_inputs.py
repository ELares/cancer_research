"""Descriptor reports distinguish missing inputs from valid zero observations."""
from copy import deepcopy
import gzip
import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPO = Path(__file__).resolve().parent.parent
OLD_JSON = b'Previous published JSON.\n'
OLD_MD = b'Previous published Markdown.\n'


@pytest.fixture
def scanner(tmp_path, monkeypatch):
    atlas = tmp_path / "external-atlas"
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(atlas))
    spec = importlib.util.spec_from_file_location(
        "descriptor_input_test", REPO / "scripts" / "atlas_descriptor_recall.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.ATLAS == atlas
    module.OUT_JSON = tmp_path / "report.json"
    module.OUT_MD = tmp_path / "report.md"
    module.OUT_JSON.write_bytes(OLD_JSON)
    module.OUT_MD.write_bytes(OLD_MD)
    monkeypatch.setattr(sys, "argv", ["atlas_descriptor_recall.py"])
    return module


def _write_shard(scanner, records, name="part-000.jsonl.gz"):
    path = scanner.ATLAS / "records" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record) + "\n")
    return path


def _snapshot(scanner):
    return scanner.OUT_JSON.read_bytes(), scanner.OUT_MD.read_bytes()


def _counts(scanner, pdt, sdt, *, total=None):
    arms = {}
    for name, counts in zip(("PDT", "SDT"), (pdt, sdt)):
        text, descriptor, both = counts
        arms[name] = {
            "text": text, "descriptor": descriptor, "both": both,
            "descriptors": sorted(scanner.ARMS[name]["descriptors"]),
            "label": scanner.ARMS[name]["label"], "descriptor_not_text": [],
        }
    return {
        "subject": "ferroptosis",
        "subject_articles": (sum(t + d - b for t, d, b in (pdt, sdt))
                             if total is None else total),
        "arms": arms, "ratio_pair": ["PDT", "SDT"], "manuscript_ratio": 2.93,
    }


@pytest.mark.parametrize("state", ["missing", "directory", "empty-shard"])
def test_unavailable_input_preserves_existing_reports(scanner, state):
    if state == "directory":
        (scanner.ATLAS / "records").mkdir(parents=True)
    elif state == "empty-shard":
        _write_shard(scanner, [])
    with pytest.raises(SystemExit, match="No census records") as error:
        scanner.main()
    assert str(scanner.ATLAS / "records") in str(error.value)
    assert "FERRO_ATLAS_ROOT" in str(error.value)
    assert "--render-only" in str(error.value)
    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("corruption", ["json", "gzip", "truncated-gzip"])
def test_later_corruption_preserves_both_reports(scanner, corruption):
    _write_shard(scanner, [{"mesh": ["Ferroptosis"], "title": "PDT SDT"}])
    later = _write_shard(scanner, [{"mesh": []}], "part-001.jsonl.gz")
    if corruption == "json":
        with gzip.open(later, "wt") as stream:
            stream.write("{broken JSON\n")
    elif corruption == "gzip":
        later.write_bytes(b"not gzip")
    else:
        later.write_bytes(later.read_bytes()[:-8])
    with pytest.raises((ValueError, OSError, EOFError)):
        scanner.main()
    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("invalid", [
    [], None, "record", {"mesh": "Ferroptosis"},
    {"mesh": {"Ferroptosis": True}}, {"mesh": [None]},
    {"mesh": [], "title": []}, {"mesh": [], "abstract": 0},
])
def test_record_shapes_are_validated_before_subject_filtering(scanner, invalid):
    _write_shard(scanner, [{"mesh": ["Ferroptosis"]}, invalid])
    with pytest.raises(ValueError, match="Census"):
        scanner.main()
    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("record, subject_count", [
    ({"mesh": ["Unrelated"], "title": None, "abstract": None}, 0),
    ({"mesh": ["FERROPTOSIS"], "title": None, "abstract": None}, 1),
])
def test_readable_zero_matches_are_written(scanner, record, subject_count, capsys):
    _write_shard(scanner, [record])
    scanner.main()
    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["subject_articles"] == subject_count
    assert all(s["text"] == s["descriptor"] == s["both"] == 0
               for s in result["arms"].values())
    assert result["ratio_by_text"] is None
    assert not result["descriptor_route_significantly_inflates"]
    assert not result["symmetric_ratio_covers_manuscript"]
    markdown = scanner.OUT_MD.read_text()
    assert "| PDT | `Photochemotherapy` | 0 | 0 | 0 | **n/a** | n/a |" in markdown
    assert "a recall difference cannot be assessed" in markdown
    assert "recall n/a" in capsys.readouterr().out
    if subject_count == 0:
        assert "no article matched the subject descriptor" in markdown


def test_scan_counts_title_abstract_descriptor_and_overlap_independently(scanner):
    _write_shard(scanner, [
        {"mesh": ["Ferroptosis", "PHOTOCHEMOTHERAPY"], "title": "photo-dynamic"},
        {"mesh": ["Ferroptosis"], "abstract": "PDT and sonodynamic"},
        {"mesh": ["Ferroptosis", "Ultrasonic Therapy"], "title": "An unrelated word"},
        {"mesh": ["Ferroptosis", "Ultrasonic Therapy"], "abstract": "SDT"},
        {"mesh": ["Unrelated", "Photochemotherapy"], "title": "PDT"},
    ])
    scanner.main()
    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["subject_articles"] == 4
    assert [(s["text"], s["descriptor"], s["both"])
            for s in result["arms"].values()] == [(2, 1, 1), (2, 2, 1)]
    assert result["descriptor_inflation"] == 0.5
    assert result["recall_asymmetry"] == 1
    assert "observed recalls are equal at 50.0%" in scanner.OUT_MD.read_text()


@pytest.mark.parametrize("pdt,sdt,expected", [
    ((1, 1, 1), (1, 1, 1), "observed recalls are equal at 100.0%"),
    ((2, 1, 0), (2, 1, 0), "observed recalls are equal at 0.0%"),
])
def test_equal_recall_and_ratio_cannot_produce_inflation_claim(scanner, pdt, sdt, expected):
    result = scanner.assemble(_counts(scanner, pdt, sdt))
    assert result["descriptor_inflation"] == 1
    assert not result["descriptor_route_significantly_inflates"]
    markdown = scanner.render(result)
    assert expected in markdown
    assert "same observed ratio (multiplier **1.00x**)" in markdown
    assert "recalls differ" not in markdown
    assert "the descriptor route inflates this ratio" not in markdown


@pytest.mark.parametrize("pdt,sdt,text_ratio,descriptor_ratio,multiplier", [
    ((2, 0, 0), (1, 1, 1), 2, 0, 0),
    ((0, 1, 0), (1, 1, 1), 0, 1, None),
    ((1, 1, 1), (0, 0, 0), None, None, None),
    ((2, 0, 0), (1, 0, 0), 2, None, None),
])
def test_sparse_axes_preserve_zero_and_undefined_distinction(
    scanner, pdt, sdt, text_ratio, descriptor_ratio, multiplier,
):
    result = scanner.assemble(_counts(scanner, pdt, sdt))
    assert result["ratio_by_text"] == text_ratio
    assert result["ratio_by_descriptor"] == descriptor_ratio
    assert result["descriptor_inflation"] == multiplier
    assert result["descriptor_inflation_ci"] is None
    assert not result["descriptor_route_significantly_inflates"]
    markdown = scanner.render(result)
    assert "multiplier interval is unavailable" in markdown
    if multiplier == 0:
        assert "multiplier **0.00x**" in markdown
        assert "| PDT | `Photochemotherapy` | 2 | 0 | 0 | **0.0%** | n/a |" in markdown


@pytest.mark.parametrize("pdt,sdt,inflates,claim", [
    ((100, 100, 100), (100, 50, 50), True, "inflates"),
    ((100, 50, 50), (100, 100, 100), False, "reduces"),
])
def test_direction_and_significance_follow_counts(scanner, pdt, sdt, inflates, claim):
    result = scanner.assemble(_counts(scanner, pdt, sdt))
    assert result["descriptor_route_significantly_inflates"] is inflates
    lo, hi = result["descriptor_inflation_ci"]
    assert (lo > 1) if inflates else (hi < 1)
    markdown = scanner.render(result)
    assert f"the descriptor route {claim} this ratio" in markdown
    assert "under the stated approximation" in markdown


def test_nonsignificant_point_difference_does_not_become_a_resolved_effect(scanner):
    result = scanner.assemble(_counts(scanner, (2, 2, 2), (2, 1, 1)))
    assert result["descriptor_inflation"] == 2
    assert not result["descriptor_route_significantly_inflates"]
    assert result["descriptor_inflation_ci"][0] < 1 < result["descriptor_inflation_ci"][1]
    markdown = scanner.render(result)
    assert "higher point ratio" in markdown
    assert "difference between routes is not resolved" in markdown


@pytest.mark.parametrize("pdt_text,sdt_text,claim", [
    (1000, 100, "understates the text-matched ratio"),
    (100, 100, "overstates the text-matched ratio"),
    (293, 100, "cannot distinguish"),
])
def test_manuscript_comparison_can_point_in_either_direction(scanner, pdt_text, sdt_text, claim):
    raw = _counts(scanner, (pdt_text, pdt_text, pdt_text), (sdt_text, sdt_text, sdt_text))
    assert claim in scanner.render(raw)


def test_render_rebuilds_stale_derived_flags_and_does_not_mutate_input(scanner):
    raw = _counts(scanner, (100, 50, 50), (100, 100, 100))
    raw.update(descriptor_route_significantly_inflates=True,
               symmetric_ratio_covers_manuscript=True, descriptor_inflation=100,
               descriptor_inflation_ci=[10, 200], recall_asymmetry=400)
    original = deepcopy(raw)
    markdown = scanner.render(raw)
    assert "the descriptor route reduces this ratio" in markdown
    assert "overstates the text-matched ratio" in markdown
    assert "factor of **2.00**" in markdown
    assert raw == original


def test_interval_and_ground_truth_limitations_are_explicit(scanner):
    markdown = scanner.render(_counts(scanner, (100, 100, 100), (100, 50, 50)))
    assert "not accuracy against independently adjudicated ground truth" in markdown
    assert "covariance within each arm" in markdown
    assert "does not account for cross-arm overlap" in markdown
    assert "Neither route supplies upper or lower bounds" in markdown
    assert "independent Poisson variables" in markdown


def test_offline_render_preserves_historical_json_byte_for_byte(scanner, monkeypatch, capsys):
    original = (REPO / "analysis" / "atlas-descriptor-recall.json").read_bytes()
    scanner.OUT_JSON.write_bytes(original)
    raw = json.loads(original)
    assert scanner.assemble(raw) == raw
    monkeypatch.setattr(sys, "argv", ["atlas_descriptor_recall.py", "--render-only"])
    monkeypatch.setattr(scanner, "scan", lambda: pytest.fail("render-only tried to scan"))
    scanner.main()
    assert scanner.OUT_JSON.read_bytes() == original
    assert scanner.OUT_MD.read_text() == scanner.render(raw)
    assert f"wrote {scanner.OUT_JSON}" not in capsys.readouterr().out


@pytest.mark.parametrize("failure", ["render", "serialize", "stage"])
def test_preparation_failures_preserve_both_reports(scanner, monkeypatch, failure):
    _write_shard(scanner, [{"mesh": ["Ferroptosis"], "title": "PDT SDT"}])
    def fail(*args, **kwargs):
        raise ValueError("injected preparation failure")
    if failure == "render":
        monkeypatch.setattr(scanner, "render", fail)
    elif failure == "serialize":
        monkeypatch.setattr(scanner.json, "dumps", fail)
    else:
        create = scanner.tempfile.NamedTemporaryFile
        calls = 0
        def fail_second_stage(*args, **kwargs):
            nonlocal calls
            calls += 1
            return fail() if calls == 2 else create(*args, **kwargs)
        monkeypatch.setattr(scanner.tempfile, "NamedTemporaryFile", fail_second_stage)
    with pytest.raises(ValueError, match="injected preparation failure"):
        scanner.main()
    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)
    assert not list(scanner.OUT_JSON.parent.glob(".report.*"))


def test_nonfinite_example_metadata_cannot_publish_nonstandard_json(scanner):
    _write_shard(scanner, [{
        "mesh": ["Ferroptosis", "Ultrasonic Therapy"],
        "title": "No modality text", "pmid": float("nan"),
    }])
    with pytest.raises(ValueError, match="Out of range float"):
        scanner.main()
    assert _snapshot(scanner) == (OLD_JSON, OLD_MD)


@pytest.mark.parametrize("change", [
    {"subject_articles": -1}, {"subject_articles": True},
    {"ratio_pair": ["PDT", "PDT"]}, {"ratio_pair": ["PDT", "OTHER"]},
    {"ratio_pair": "PDT"}, {"manuscript_ratio": 0},
    {"manuscript_ratio": True}, {"manuscript_ratio": float("nan")},
    {"manuscript_ratio": float("inf")},
])
def test_invalid_stored_summary_preserves_markdown_and_json(scanner, monkeypatch, change):
    raw = _counts(scanner, (1, 1, 1), (1, 1, 1))
    raw.update(change)
    scanner.OUT_JSON.write_text(json.dumps(raw))
    before = _snapshot(scanner)
    monkeypatch.setattr(sys, "argv", ["atlas_descriptor_recall.py", "--render-only"])
    with pytest.raises(ValueError):
        scanner.main()
    assert _snapshot(scanner) == before


@pytest.mark.parametrize("counts", [
    {"text": -1}, {"text": 0.5}, {"text": True}, {"both": 2},
    {"text": 2, "descriptor": 2, "both": 1},
])
def test_inconsistent_raw_counts_are_rejected(scanner, counts):
    raw = _counts(scanner, (1, 1, 1), (1, 1, 1))
    raw["arms"]["PDT"].update(counts)
    with pytest.raises(ValueError):
        scanner.assemble(raw)
