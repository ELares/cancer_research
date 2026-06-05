"""Guard the #349 living-review delta logic + the frozen/living separation.

The PubMed-touching parts are not exercised here (they need the network and run
in the scheduled Action); these tests cover the PURE delta / changelog / query
logic and assert the load-bearing invariant that the living review writes ONLY
under the living directories, never the frozen manuscript corpus.
"""

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "living_review_update", REPO / "scripts" / "living_review_update.py"
)
lr = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lr)


def test_load_queries_skips_comments_and_blanks():
    text = '# header\n\n"A"[All Fields] AND cancer\n  \n# section\n"B" AND neoplasms\n'
    assert lr.load_queries(text) == ['"A"[All Fields] AND cancer', '"B" AND neoplasms']


def test_windowed_query_appends_open_ended_date_lower_bound():
    q = lr.windowed_query('"sonodynamic therapy" AND cancer', "2025-01-01")
    assert '"sonodynamic therapy" AND cancer' in q
    assert '"2025/01/01"[Date - Publication]' in q
    assert '"3000"[Date - Publication]' in q  # open upper bound


def test_compute_delta_excludes_frozen_and_reports_per_query():
    found = {"qA": ["1", "2", "3"], "qB": ["3", "4"]}
    frozen = {"2"}  # PMID 2 already in the frozen corpus
    new_by_query, all_new = lr.compute_delta(found, frozen)
    assert new_by_query["qA"] == ["1", "3"]  # 2 dropped (frozen)
    assert new_by_query["qB"] == ["3", "4"]
    # all_new is the de-duplicated union (3 appears in both queries, counted once).
    assert all_new == ["1", "3", "4"]


def test_compute_delta_dedups_within_a_query():
    found = {"qA": ["1", "1", "2"]}
    new_by_query, all_new = lr.compute_delta(found, set())
    assert new_by_query["qA"] == ["1", "2"]
    assert all_new == ["1", "2"]


def test_is_landmark_flags_high_evidence_pubtypes_and_high_cites():
    assert lr.is_landmark({"pub_types": ["Journal Article", "Clinical Trial, Phase III"]})
    assert lr.is_landmark({"pub_types": ["Randomized Controlled Trial"]})
    assert lr.is_landmark({"pub_types": ["Journal Article"], "cited_by_count": 250})
    assert not lr.is_landmark({"pub_types": ["Journal Article"]})
    assert not lr.is_landmark({})  # missing fields must not raise


def test_format_changelog_structure_and_landmark_section():
    new_by_query = {"qA": ["10", "11"], "qB": []}
    records = [
        {"pmid": "10", "year": 2025, "title": "A phase III ferroptosis trial",
         "mechanisms": ["ferroptosis"], "pub_types": ["Clinical Trial, Phase III"]},
        {"pmid": "11", "year": 2024, "title": "A small study", "mechanisms": [],
         "pub_types": ["Journal Article"]},
    ]
    md = lr.format_changelog("2025-06-05", "2025-01-01", new_by_query, records, 4830)
    assert "Living-review update — 2025-06-05" in md
    assert "New records this window: 2." in md
    assert "New landmark detections (1)" in md  # only the phase-III record
    assert "PMID 10" in md
    assert "4,830 PMIDs" in md  # frozen count, comma-formatted
    assert "**The frozen manuscript corpus is unchanged**" in md
    assert "| qA |" in md and "| 2 |" in md  # per-query table row


def test_living_review_writes_only_under_living_dirs_not_frozen():
    """The load-bearing separation: the script's output constants point under the
    living dirs (corpus/living + analysis/living-review), never the frozen corpus
    files (corpus/by-pmid, corpus/abstracts, corpus/INDEX.jsonl)."""
    assert lr.LIVING_DIR == REPO / "corpus" / "living"
    assert lr.REPORT_DIR == REPO / "analysis" / "living-review"
    frozen = {lr.FROZEN_INDEX, *lr.FROZEN_PMID_DIRS}
    # No output path may be inside, equal to, or a parent of a frozen path.
    for out in (lr.LIVING_DIR, lr.REPORT_DIR):
        for fp in frozen:
            assert out != fp
            assert fp not in out.parents
            assert out not in fp.parents
