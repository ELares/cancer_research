"""The manuscript comparison must consume counts, including counterevidence."""
import copy
import gzip
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def report_with_recall_counts(data):
    """A manuscript comparison over the same synthetic descriptor cohort."""
    report = json.loads((ROOT / "analysis/manuscript-vs-census.json").read_text())
    report["ferroptosis_records"] = data["subject_articles"]
    table = report["modality_table"]
    for row in table["rows"]:
        if row["modality"] in data["arms"]:
            row["census_ferroptosis"] = data["arms"][row["modality"]]["descriptor"]
    pdt, sdt = (data["arms"][arm]["descriptor"] for arm in ("PDT", "SDT"))
    ratio = pdt / sdt if sdt else None
    measurable = min(pdt, sdt) >= table["min_for_a_ratio"]
    table.update(census_pdt_sdt_ratio=ratio, ratio_is_measurable=measurable,
                 census_exceeds_manuscript=bool(measurable and ratio > 2.93),
                 direction_holds=bool(measurable and ratio > 1))
    for variant in table["descriptor_variants"]:
        variant.update(pdt=pdt, sdt=sdt, ratio=ratio)
    table.update(ratio_range=[ratio, ratio],
                 direction_holds_under_every_variant=bool(ratio and ratio > 1),
                 understatement_holds_under_every_variant=bool(ratio and ratio > 2.93))
    pdt_both, sdt_both = (data["arms"][arm]["both"] for arm in ("PDT", "SDT"))
    pdt_pct = round(100 * pdt_both / pdt, 1) if pdt else 0
    sdt_pct = round(100 * sdt_both / sdt, 1) if sdt else 0
    gap = round(abs(pdt_pct - sdt_pct), 1)
    table["on_modality_and_tumour"].update(
        pdt_n=pdt_both, sdt_n=sdt_both, pdt_pct=pdt_pct, sdt_pct=sdt_pct,
        gap_points=gap, symmetric_within_5_points=gap <= 5,
        filtered_ratio=pdt_both / sdt_both if sdt_both else None)
    return report


@pytest.fixture
def consumer(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "manuscript_recall_consumer", ROOT / "scripts/manuscript_vs_census.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "PROJECT_ROOT", tmp_path)
    directory = tmp_path / "analysis"
    directory.mkdir()
    path = directory / "atlas-descriptor-recall.json"
    raw = json.loads((ROOT / "analysis/atlas-descriptor-recall.json").read_text())
    raw["subject_articles"] = 10000

    def write(counts, pair=("PDT", "SDT")):
        data = copy.deepcopy(raw)
        data["ratio_pair"] = list(pair)
        for arm, (text, descriptor, both) in counts.items():
            data["arms"][arm].update(text=text, descriptor=descriptor, both=both)
        # Poisoning these cached values must not affect consumption.
        data.update(ratio_by_text=999, ratio_by_text_ci=[998, 1000],
                    symmetric_ratio_covers_manuscript=False)
        for arm in data["arms"].values():
            arm.update(recall=999, precision=999)
        path.write_text(json.dumps(data))
        return data

    return module, path, write


@pytest.mark.parametrize("text_pdt,relation,phrase", [
    (600, "above", "interval lies above"),
    (300, "overlaps", "cannot distinguish"),
    (100, "below", "interval lies below"),
    (0, "unavailable", "interval is unavailable"),
])
@pytest.mark.parametrize("descriptor_ratio", [1.0, 4.75])
def test_interval_direction_controls_all_recall_consuming_prose(
        consumer, text_pdt, relation, phrase, descriptor_ratio):
    module, _path, write = consumer
    pdt_descriptor = int(descriptor_ratio * 100)
    data = write({"PDT": (text_pdt, pdt_descriptor,
                          min(text_pdt, pdt_descriptor)),
                  "SDT": (100, 100, 80)})
    checked = module._recall_check()
    assert checked["symmetric_ratio"] == text_pdt / 100
    assert checked["interval_relation"] == relation
    assert checked["symmetric_agrees"] == (relation == "above")
    report = report_with_recall_counts(data)
    rendered = module.render(report).lower()
    assert phrase in module._recall_caveat().lower()
    assert phrase in rendered
    headline = module._headline(report).lower()
    if relation != "above":
        assert phrase in headline
        assert "survives, understated by the manuscript" not in rendered
    else:
        assert "supports understatement under the count model" in headline
        assert "supports understatement for the text-count comparison" in rendered
    assert "omits covariance between arms" in rendered
    if relation in {"below", "unavailable"}:
        assert "cannot distinguish" not in rendered


def test_equal_recalls_keep_both_arms_and_independently_report_precision(consumer):
    module, _path, write = consumer
    write({"PDT": (100, 200, 100), "SDT": (100, 100, 100)})
    caveat = module._recall_caveat()
    assert "(PDT): 100.0% recall, 50.0% precision" in caveat
    assert "(SDT): 100.0% recall, 100.0% precision" in caveat
    assert "observed recalls are equal" in caveat
    assert "lopsided" not in caveat and "precision really is symmetric" not in caveat


def test_zero_and_undefined_recalls_are_distinct(consumer):
    module, _path, write = consumer
    write({"PDT": (100, 100, 0), "SDT": (100, 100, 60)})
    assert "0.0% recall" in module._recall_caveat()
    assert "finite fold gap is unavailable" in module._recall_caveat()
    write({"PDT": (0, 5, 0), "SDT": (0, 0, 0)})
    caveat = module._recall_caveat()
    assert "(PDT): recall unavailable, 0.0% precision" in caveat
    assert "(SDT): recall unavailable, precision unavailable" in caveat
    assert "symmetric interval is unavailable" in caveat


def test_incompatible_ratio_order_is_rejected(consumer):
    module, _path, write = consumer
    write({"PDT": (300, 300, 300), "SDT": (100, 100, 100)},
          pair=("SDT", "PDT"))
    with pytest.raises(SystemExit, match="required by Section"):
        module._recall_check()


@pytest.mark.parametrize("field,value", [("subject", "apoptosis"),
                                        ("manuscript_ratio", 1.0)])
def test_incompatible_subject_or_manuscript_ratio_is_rejected(consumer, field, value):
    module, path, write = consumer
    data = write({"PDT": (300, 300, 300), "SDT": (100, 100, 100)})
    data[field] = value
    path.write_text(json.dumps(data))
    with pytest.raises(SystemExit, match="Section"):
        module._recall_check()


def test_required_recall_evidence_fails_closed(consumer):
    module, path, write = consumer
    with pytest.raises(SystemExit, match="missing"):
        module._recall_check()
    path.write_text("not json")
    with pytest.raises(SystemExit, match="unreadable"):
        module._recall_check()
    write({"PDT": (100, 100, 101), "SDT": (100, 100, 60)})
    with pytest.raises((SystemExit, ValueError)):
        module._recall_check()


@pytest.mark.parametrize("mismatch", ["subject_articles", "PDT", "SDT"])
def test_recall_from_a_different_census_cannot_restore_understatement(consumer, mismatch):
    module, path, write = consumer
    data = write({"PDT": (600, 600, 600), "SDT": (100, 100, 100)})
    report = report_with_recall_counts(data)
    assert "supports understatement" in module._headline(report)
    if mismatch == "subject_articles":
        data[mismatch] += 1
    else:
        # The same subject total does not excuse different descriptor margins.
        data["arms"][mismatch]["descriptor"] += 1
    path.write_text(json.dumps(data))
    with pytest.raises(SystemExit, match="census counts"):
        module._headline(report)
    with pytest.raises(SystemExit, match="census counts"):
        module.render(report)


@pytest.mark.parametrize("evidence", ["missing", "subject_articles", "PDT", "SDT"])
def test_unusable_required_recall_preserves_existing_output_pair(
        consumer, monkeypatch, evidence):
    module, path, write = consumer
    records = path.parent / "records"
    records.mkdir()
    with gzip.open(records / "fixture.jsonl.gz", "wt") as stream:
        stream.write(json.dumps({
            "pmid": "1", "year": 2015,
            "mesh": ["Ferroptosis", "Photochemotherapy", "Ultrasonic Therapy"],
        }) + "\n")
    monkeypatch.setattr(module, "RECORDS", records)
    measured = module.scan()
    monkeypatch.setattr(module, "scan", lambda: measured)
    if evidence != "missing":
        data = write({"PDT": (1, 1, 1), "SDT": (1, 1, 1)})
        data["subject_articles"] = 2 if evidence == "subject_articles" else 1
        if evidence in {"PDT", "SDT"}:
            data["arms"][evidence].update(descriptor=0, both=0)
        path.write_text(json.dumps(data))
    outputs = [path.parent / "prior.json", path.parent / "prior.md"]
    monkeypatch.setattr(module, "OUT_JSON", outputs[0])
    monkeypatch.setattr(module, "OUT_MD", outputs[1])
    for output in outputs:
        output.write_bytes(b"published result sentinel\n")
    with pytest.raises(SystemExit, match="missing" if evidence == "missing"
                       else "census counts"):
        module.main()
    assert all(output.read_bytes() == b"published result sentinel\n"
               for output in outputs)
