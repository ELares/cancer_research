"""Guard the duplicate-corpus audit script (#535).

These tests keep the detector lightweight: verify the title normalizer and the
load-bearing report split between preprint/published twins and other same-title
collisions.
"""

import importlib.util
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "detect_corpus_duplicates", REPO / "scripts" / "detect_corpus_duplicates.py"
)
dcd = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dcd)


def test_norm_title_lowercases_and_strips_non_alnum():
    assert dcd.norm_title(" Ferroptosis: In Cancer?! ") == "ferroptosisincancer"
    assert dcd.norm_title("") == ""
    assert dcd.norm_title(None) == ""


def test_main_writes_preprint_twin_and_other_collision_sections(tmp_path, monkeypatch, capsys):
    index_path = tmp_path / "INDEX.jsonl"
    report_path = tmp_path / "corpus-duplicate-audit.md"

    records = [
        {
            "pmid": "1001",
            "title": "Long shared ferroptosis title in pancreatic cancer",
            "journal": "bioRxiv",
            "year": 2024,
            "oa_status": "green",
        },
        {
            "pmid": "1002",
            "title": "Long shared ferroptosis title in pancreatic cancer",
            "journal": "Cancer Discovery",
            "year": 2025,
            "oa_status": "gold",
        },
        {
            "pmid": "2001",
            "title": "Independent same title collision in glioblastoma models",
            "journal": "Journal A",
            "year": 2021,
            "oa_status": "bronze",
        },
        {
            "pmid": "2002",
            "title": "Independent same title collision in glioblastoma models",
            "journal": "Journal B",
            "year": 2022,
            "oa_status": "hybrid",
        },
        {
            "pmid": "3001",
            "title": "short title",
            "journal": "medRxiv",
            "year": 2023,
            "oa_status": "green",
        },
        {
            "pmid": "3002",
            "title": "short title",
            "journal": "Nature",
            "year": 2024,
            "oa_status": "gold",
        },
    ]
    index_path.write_text(
        "\n".join(json.dumps(record) for record in records) + "\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(dcd, "INDEX", index_path)
    monkeypatch.setattr(dcd, "REPORT", report_path)
    monkeypatch.setattr(dcd, "REPO_ROOT", tmp_path)

    dcd.main()

    out = capsys.readouterr().out
    report = report_path.read_text(encoding="utf-8")

    assert "1 preprint/published twin groups, 1 other same-title groups." in out
    assert "## Preprint/published twins" in report
    assert "## Other same-title collisions (reviewed individually, #567)" in report
    assert "| 1001 | bioRxiv | 2024 | green | preprint |" in report
    assert "| 1002 | Cancer Discovery | 2025 | gold | published |" in report
    assert "| 2001 | Journal A | 2021 | bronze | published |" in report
    assert "| 2002 | Journal B | 2022 | hybrid | published |" in report
    assert "short title" not in report
    # the synthetic "other" group has no hand-verified verdict -> flagged for review (#567)
    assert "not yet reviewed" in report


def test_verdict_for_known_and_unreviewed_groups():
    # known same-title group -> the hand-verified verdict (order-independent: frozenset)
    known = dcd.verdict_for([{"pmid": "35433483"}, {"pmid": "38487722"}])
    assert "two distinct corrigenda" in known and "Published Erratum" in known
    # int vs str PMIDs reconcile, and the letters group resolves too
    letters = dcd.verdict_for([{"pmid": 19997112}, {"pmid": 19997110}])
    assert "two independent letters" in letters
    # an unknown group is flagged not-yet-reviewed rather than presumed benign
    assert "not yet reviewed" in dcd.verdict_for([{"pmid": "111"}, {"pmid": "222"}])

