"""Missing collision evidence must not become a negative module-claim result."""

import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import atlas_module_support as module_support


def ambiguity():
    return {"by_type": {kind: {"sense_rows": []}
                        for kind in ("gene", "chemical", "disease")}}


@pytest.fixture
def report(monkeypatch, tmp_path):
    """Exercise the producer with one contested edge and temporary artifacts."""
    (tmp_path / "analysis").mkdir()
    idx = {"edges": {("1", "2"): {"positive_correlate": 3,
                                    "negative_correlate": 3}}}
    monkeypatch.setattr(module_support, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(module_support, "OUT", tmp_path / "report.md")
    monkeypatch.setattr(module_support, "RAW", tmp_path / "report.json")
    monkeypatch.setattr(module_support, "atlas_root", lambda: tmp_path)
    monkeypatch.setattr(module_support, "load_index", lambda root: idx)
    monkeypatch.setattr(module_support, "load_comentions", lambda root: {})
    monkeypatch.setattr(module_support, "CLAIMS", [
        ("synthetic", "A", "B", "1", "Synthetic claim")])
    monkeypatch.setattr(module_support, "resolve", lambda idx, name: {
        "A": "1", "B": "2"}.get(name))
    monkeypatch.setattr(module_support, "support", lambda idx, a, b: {
        "a": "1", "b": "2", "a_name": "A", "b_name": "B", "total": 6,
        "predicates": idx["edges"][("1", "2")], "pmids": ["1"]})
    monkeypatch.setattr(module_support, "_relation_scale", lambda: "6")
    monkeypatch.setattr(module_support, "_exposure_section", lambda *a, **kw: [])
    return tmp_path / "analysis/atlas-ambiguity.json"


@pytest.mark.parametrize("problem", [
    "missing", "malformed", "null", "missing-group", "invalid-rows", "invalid-id",
])
def test_unavailable_collision_scan_is_not_reported_as_no_matches(report, problem):
    scan = ambiguity()
    if problem == "null":
        scan = None
    elif problem == "missing-group":
        del scan["by_type"]["chemical"]
    elif problem == "invalid-rows":
        scan["by_type"]["gene"]["sense_rows"] = None
    elif problem == "invalid-id":
        scan["by_type"]["gene"]["sense_rows"] = [
            {"top": {"id": "1"}, "runner_up": {"id": 2}}]
    if problem != "missing":
        report.write_text("{broken" if problem == "malformed" else json.dumps(scan))

    module_support.main()
    text = module_support.OUT.read_text()
    assert "collision check is unavailable" in text
    assert "not checked for measured collisions" in text
    assert "This check identified no measured collision" not in text
    assert "matched the collision set used here" not in text


@pytest.mark.parametrize("ids, has_collision", [
    (None, False), (("3", "4"), False), (("1", "3"), True),
])
def test_valid_scan_distinguishes_matches_from_no_matches(report, ids, has_collision):
    scan = ambiguity()
    if ids:
        scan["by_type"]["gene"]["sense_rows"] = [
            {"top": {"id": ids[0]}, "runner_up": {"id": ids[1]}}]
    report.write_text(json.dumps(scan))

    module_support.main()
    text = module_support.OUT.read_text()
    assert "collision check is unavailable" not in text
    if has_collision:
        assert "These claims need a conflation check" in text
        assert "`synthetic`: A" in text
        assert "This check identified no measured collision" not in text
    else:
        assert "This check identified no measured collision" in text
        assert "None of the 1 claims matched the collision set used here" in text
        assert "does not exclude unmeasured" in text


def test_fixed_association_is_identified_as_historical(report):
    report.write_text(json.dumps(ambiguity()))
    module_support.main()
    for text in (module_support.OUT.read_text(),
                 (ROOT / "analysis/atlas-module-support.md").read_text()):
        assert "historical raw directional diagnostic reported a 1.45x association" in text
