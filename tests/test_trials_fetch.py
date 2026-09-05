"""Guards for the ClinicalTrials.gov fetcher.

This source is public domain, so the licence question that dominates every
other fetcher in this project does not arise. What DOES carry over is the
durability contract -- write before index, never the reverse -- and the reason
these tests exist is that the contract is IMPORTED rather than reimplemented.
A guard here that passed on a local copy of the ordering would be worthless;
these check the imported machinery is what runs.
"""
from __future__ import annotations

import gzip
import json
import sqlite3
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import corpus_expand_fetch as fx  # noqa: E402
import corpus_identity_index as ix  # noqa: E402
import trials_fetch as tf  # noqa: E402


def _study(nct, *, results=False, pmids=(), phases=("PHASE2",)):
    s = {"protocolSection": {
        "identificationModule": {"nctId": nct, "briefTitle": f"T{nct}"},
        "statusModule": {"overallStatus": "COMPLETED"},
        "designModule": {"phases": list(phases), "studyType": "INTERVENTIONAL",
                         "enrollmentInfo": {"count": 10}},
        "conditionsModule": {"conditions": ["Neoplasms"]},
        "armsInterventionsModule": {"interventions": [{"name": "drug"}]},
        "outcomesModule": {"primaryOutcomes": [{"measure": "OS"}]},
        "descriptionModule": {"briefSummary": "s"},
        "referencesModule": {"references": [{"pmid": p} for p in pmids]},
    }}
    if results:
        s["resultsSection"] = {"x": 1}
    return s


def _drive(tmp_path, monkeypatch, pages, limit=None):
    monkeypatch.setattr(tf, "OUT_ROOT", tmp_path / "out")
    monkeypatch.setattr(ix, "DB", tmp_path / "id.sqlite")
    monkeypatch.setattr(tf, "SLEEP", 0)
    monkeypatch.setattr(tf.time, "sleep", lambda *_: None)
    seq = list(pages)
    monkeypatch.setattr(tf, "_page",
                        lambda token, cond: seq.pop(0) if seq else {"studies": []})
    totals = tf.run("cancer", limit=limit, verbose=False)
    stored = []
    for f in sorted((tmp_path / "out").rglob("*.jsonl.gz")):
        try:
            with gzip.open(f, "rt") as fh:
                stored += [json.loads(l)["nct_id"] for l in fh]
        except EOFError:
            pass
    return totals, stored


def test_a_trial_is_stored_and_indexed_under_its_nct_id(tmp_path, monkeypatch):
    totals, stored = _drive(tmp_path, monkeypatch,
                            [{"studies": [_study("NCT1"), _study("NCT2")],
                              "totalCount": 2}])
    assert stored == ["NCT1", "NCT2"]
    c = sqlite3.connect(tmp_path / "id.sqlite")
    held = {k for (k,) in c.execute("SELECT key FROM held WHERE kind='nct'")}
    c.close()
    assert held == {"NCT1", "NCT2"}, (
        f"stored {stored} but indexed {held}; the next run would re-download")


def test_a_trial_already_held_is_not_fetched_again(tmp_path, monkeypatch):
    """The user requirement that governs every fetcher here."""
    page = {"studies": [_study("NCT1"), _study("NCT2")], "totalCount": 2}
    _drive(tmp_path, monkeypatch, [page])
    totals, stored = _drive(tmp_path, monkeypatch, [page])
    assert totals["new"] == 0, "a held trial was downloaded a second time"
    assert stored == ["NCT1", "NCT2"], "the second run duplicated the shard"


def test_the_nct_namespace_cannot_collide_with_article_keys(tmp_path, monkeypatch):
    """A trial is not an article. Forcing NCT ids into the PMID or DOI kind
    would make an article look held because a trial number matched it."""
    c = sqlite3.connect(tmp_path / "id.sqlite")
    c.executescript(ix.SCHEMA)
    c.execute("INSERT INTO held VALUES ('pmid','12345','census',0)")
    c.commit()
    c.close()
    totals, stored = _drive(tmp_path, monkeypatch,
                            [{"studies": [_study("12345")], "totalCount": 1}])
    assert stored == ["12345"], (
        "a trial was skipped because an ARTICLE with the same digits is held")
    assert tf.KIND == "nct" and tf.KIND not in ("pmid", "pmcid", "doi", "epmc")


def test_a_study_without_an_nct_id_is_refused(tmp_path, monkeypatch):
    """Nothing to dedup on later, so storing it guarantees a re-download."""
    bad = {"protocolSection": {"identificationModule": {"briefTitle": "no id"}}}
    totals, stored = _drive(tmp_path, monkeypatch,
                            [{"studies": [bad, _study("NCT9")], "totalCount": 2}])
    assert stored == ["NCT9"]


def test_the_durability_contract_is_imported_not_copied():
    """The write-then-index ordering, the flushed shards and the transient
    failure type were each got wrong at least once and fixed under review. A
    second copy would be a second chance to get them wrong, and the two would
    drift apart silently."""
    assert tf.Shards is fx.Shards
    assert tf.TransientFetchError is fx.TransientFetchError
    assert tf._RateLimit is fx._RateLimit
    src = Path(tf.__file__).read_text()
    assert "class Shards" not in src, "the shard writer was reimplemented here"
    assert "def _get(" not in src, "the fetch/retry loop was reimplemented here"


def test_a_registry_search_refusal_is_not_read_as_the_end_of_the_data(monkeypatch):
    """`_get` returns None for a refusal, which for ONE item is an answer. For
    a paged search it would silently truncate the crawl at that page."""
    monkeypatch.setattr(tf, "_get", lambda *a, **k: None)
    try:
        tf._page(None, "cancer")
    except tf.TransientFetchError:
        return
    raise AssertionError("a refused search page returned quietly")


def test_linked_pmids_are_kept_because_they_are_the_join(tmp_path, monkeypatch):
    """A registry entry citing its own publications is what makes a trial
    findable from a paper, and vice versa."""
    totals, _ = _drive(tmp_path, monkeypatch,
                       [{"studies": [_study("NCT1", pmids=["111", "222"])],
                         "totalCount": 1}])
    assert totals["linked_pmids"] == 2
