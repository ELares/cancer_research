"""Guards for the duplicate audit.

The audit's whole job is to answer "is every downloaded article new?", and the
two ways it can lie are opposite: counting the crawl's own writes as duplicates
(circularity), or never finding a duplicate at all (a dead predicate). Both
produce clean-looking output, so both are tested.
"""
from __future__ import annotations

import gzip
import json
import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from corpus_expand_verify import held_by_others, iter_records, verify  # noqa: E402
from corpus_identity_index import SCHEMA  # noqa: E402


def _index(tmp_path: Path, rows) -> Path:
    p = tmp_path / "identity.sqlite"
    c = sqlite3.connect(p)
    c.executescript(SCHEMA)
    c.executemany(
        "INSERT INTO held(kind,key,source,has_fulltext) VALUES (?,?,?,?)", rows)
    c.commit()
    c.close()
    return p


def _shards(tmp_path: Path, records) -> Path:
    d = tmp_path / "shards"
    d.mkdir()
    with gzip.open(d / "expanded-00000.jsonl.gz", "wt") as fh:
        for r in records:
            fh.write(json.dumps(r) + "\n")
    return d


def test_the_crawls_own_rows_are_not_duplicates(tmp_path):
    """THE circularity trap.

    The crawl indexes each record in the pass that writes it, so by audit time
    every record is present under an `expand-*` source. Counting those makes a
    perfectly clean crawl report 100% duplicates -- observed live as 6,807 of
    7,472. Only pre-existing provenance may answer the question.
    """
    idx = _index(tmp_path, [("pmid", "40000001", "expand-MED", 1)])
    c = sqlite3.connect(f"file:{idx}?mode=ro", uri=True)
    assert held_by_others(c, "40000001", None, None) is False
    c.close()


def test_a_genuinely_preexisting_article_is_a_duplicate(tmp_path):
    """The other direction: real prior holdings must still be caught."""
    idx = _index(tmp_path, [("pmid", "40000001", "census-mesh", 0)])
    c = sqlite3.connect(f"file:{idx}?mode=ro", uri=True)
    assert held_by_others(c, "40000001", None, None) is True
    c.close()


def test_any_one_key_is_enough(tmp_path):
    """The same article arrives under different identifiers from different
    sources. A PMID miss says nothing about whether we hold it by DOI."""
    idx = _index(tmp_path, [("doi", "10.1234/abc", "census-text", 1)])
    c = sqlite3.connect(f"file:{idx}?mode=ro", uri=True)
    assert held_by_others(c, "99999999", None, "https://doi.org/10.1234/ABC") is True
    c.close()


def test_control_failing_fails_the_whole_audit(tmp_path):
    """A predicate that never returns True reports zero duplicates -- which
    looks like success. The control is what separates the two, so an index with
    no pre-existing rows to control against must NOT pass."""
    idx = _index(tmp_path, [("pmid", "40000001", "expand-MED", 1)])
    sd = _shards(tmp_path, [{"pmid": "40000002"}])
    r = verify(sd, idx)
    assert r["duplicates"] == 0
    assert r["control_ok"] is False, "no pre-existing rows means nothing was proved"


def test_clean_crawl_passes_with_a_live_control(tmp_path):
    idx = _index(tmp_path, [
        ("pmid", "30000001", "census-mesh", 0),
        ("pmid", "30000002", "census-mesh", 0),
        ("pmid", "40000001", "expand-MED", 1),
    ])
    sd = _shards(tmp_path, [{"pmid": "40000001", "doi": "10.9/new"}])
    r = verify(sd, idx)
    assert (r["duplicates"], r["control_ok"], r["written"]) == (0, True, 1)


def test_a_real_duplicate_is_reported(tmp_path):
    idx = _index(tmp_path, [("pmid", "30000001", "census-mesh", 0)])
    sd = _shards(tmp_path, [{"pmid": "30000001"}])
    r = verify(sd, idx)
    assert r["duplicates"] == 1 and r["control_ok"] is True


def test_records_with_no_identifier_never_count_as_held(tmp_path):
    """Nothing to match on must mean 'not held', never a silent match on NULL."""
    idx = _index(tmp_path, [("pmid", "30000001", "census-mesh", 0)])
    c = sqlite3.connect(f"file:{idx}?mode=ro", uri=True)
    assert held_by_others(c, None, None, None) is False
    c.close()


def test_a_shard_still_being_written_does_not_abort_the_audit(tmp_path):
    """The crawl runs for days; auditing mid-run must read what is there rather
    than raise on the truncated gzip at the end."""
    d = tmp_path / "shards"
    d.mkdir()
    with gzip.open(d / "expanded-00000.jsonl.gz", "wt") as fh:
        fh.write(json.dumps({"pmid": "1"}) + "\n")
    (d / "expanded-00001.jsonl.gz").write_bytes(b"\x1f\x8b\x08\x00truncated")
    assert [r["pmid"] for r in iter_records(d)] == ["1"]


# --- Order-independence for the crawl report -------------------------------
#
# `corpus_expand_report` is EXEMPT in tests/test_artifact_freshness.py, because
# its input is a corpus that lives outside the repository and CI has no copy to
# regenerate from. Exemption there removes it from LIVE, which also removes it
# from that file's order-independence gate -- and order-dependence is a defect
# it actually had: both its tables published whatever order the counting dict
# happened to have. The gate is reinstated here over SYNTHETIC data, so it needs
# no shards and no committed artifact.

def _shuffled(obj, rng):
    if isinstance(obj, dict):
        items = list(obj.items())
        rng.shuffle(items)
        return {k: _shuffled(v, rng) for k, v in items}
    if isinstance(obj, list):
        return [_shuffled(v, rng) for v in obj]
    return obj


def _report_fixture():
    """Shaped from the generator's real output, not invented -- a fixture that
    omits a key the renderer reads tests nothing and raises KeyError."""
    return {
        "generated": "2026-01-01 00:00 UTC",
        "crawl_running": False,
        "records": 60, "with_fulltext": 12, "fulltext_chars": 1234, "shards": 2,
        "redistributable": 40, "redistributable_with_text": 8,
        "by_source": {"MED": 30, "PPR": 20, "PMC": 10},
        "fulltext_by_source": {"MED": 8, "PMC": 4},
        "licences": {"cc by": 25, "cc0": 25, "unknown": 10},
        "years": {"2026": 40, "2025": 20},
        "slices": [{"src": "MED", "year": 2026, "seen": 40, "kept": 30,
                    "fulltext": 8, "done": 1, "slices": 1}],
    }


def test_the_crawl_report_order_is_established_by_the_renderer():
    """Shuffling every dict must not change a single character of the output."""
    import random
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import corpus_expand_report as m

    base = m.render(_report_fixture())
    rng = random.Random(20260904)
    for _ in range(20):
        assert m.render(_shuffled(_report_fixture(), rng)) == base


def test_equal_counts_still_have_one_published_order():
    """The hardest case to notice: two rows tie, so the count cannot rank them
    and the label has to. Without a tie-break they swap between runs."""
    import random
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import corpus_expand_report as m

    d = _report_fixture()
    d["by_source"] = {"MED": 10, "PMC": 10, "PPR": 10}   # a three-way tie
    d["licences"] = {"cc by": 5, "cc0": 5}               # and a two-way tie
    base = m.render(d)
    rng = random.Random(7)
    for _ in range(20):
        assert m.render(_shuffled(d, rng)) == base
