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
        # THREE slices, two of them TIED on every ranked quantity. One slice
        # made the shuffle test vacuous for this table: _shuffled preserves
        # list order, so a single-element list can hide any ordering defect,
        # and a tie is the case a count alone cannot resolve.
        "slices": [
            {"src": "MED", "year": 2026, "seen": 40, "kept": 30, "fulltext": 8,
             "done": 1, "slices": 1},
            {"src": "PMC", "year": 2026, "seen": 20, "kept": 10, "fulltext": 2,
             "done": 1, "slices": 1},
            {"src": "PPR", "year": 2026, "seen": 20, "kept": 10, "fulltext": 2,
             "done": 1, "slices": 1},
        ],
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


# --- The index must not claim articles whose bytes are absent ---------------
#
# The failure these pin actually happened to real data: the fetcher indexed a
# record while its bytes sat in an unflushed gzip buffer, two kills discarded
# the buffers, and 374 index rows were left describing articles that are not on
# disk. Such a row is worse than a missing one -- it says "held", so every
# later run skips the article forever, silently and permanently.

def _shard_file(d: Path, name: str, records, truncate: bool = False):
    d.mkdir(parents=True, exist_ok=True)
    raw = b"".join((json.dumps(r) + "\n").encode() for r in records)
    import io
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as g:
        g.write(raw)
    data = buf.getvalue()
    if truncate:
        data = data[:-8]  # drop the gzip trailer, as an unclosed writer does
    (d / name).write_bytes(data)


def test_an_index_row_with_no_stored_record_is_reported(tmp_path):
    from corpus_expand_verify import orphaned_index_rows
    idx = _index(tmp_path, [
        ("pmid", "111", "expand-MED", 0),   # stored below
        ("pmid", "222", "expand-MED", 0),   # NOT stored -> orphan
        ("pmid", "333", "census-mesh", 0),  # not ours; never an orphan
    ])
    sd = tmp_path / "shards"
    _shard_file(sd, "expanded-00000.jsonl.gz", [{"pmid": "111"}])
    orphans, truncated = orphaned_index_rows(idx, sd)
    assert orphans == [("pmid", "222")], orphans
    assert truncated == []


def test_a_truncated_shard_is_named_rather_than_skipped(tmp_path):
    """Silently reading the readable prefix is how three truncated shards sat
    unnoticed while the index vouched for their lost tails."""
    from corpus_expand_verify import stored_keys
    sd = tmp_path / "shards"
    _shard_file(sd, "expanded-00000.jsonl.gz",
                [{"pmid": str(i)} for i in range(1, 6)], truncate=True)
    keys, truncated = stored_keys(sd)
    assert ("pmid", "1") in keys, "recoverable records must still be read"
    assert len(truncated) == 1 and truncated[0][0] == "expanded-00000.jsonl.gz"
    assert truncated[0][1] >= 1, "the count of what WAS recovered is part of the report"


def test_verify_fails_when_the_index_overclaims(tmp_path):
    idx = _index(tmp_path, [
        ("pmid", "30000001", "census-mesh", 0),   # control fodder
        ("pmid", "222", "expand-MED", 0),         # orphan
    ])
    sd = tmp_path / "shards"
    _shard_file(sd, "expanded-00000.jsonl.gz", [{"pmid": "111"}])
    r = verify(sd, idx)
    assert r["orphaned_index_rows"] == 1
    assert r["duplicates"] == 0, "an overclaim is not a duplicate; they are different faults"


def test_the_slice_table_order_survives_a_reordered_list():
    """Shuffling a LIST is not what `_shuffled` does, so this reverses it
    explicitly -- the one transformation that exposes an order inherited from
    the query rather than established by the renderer."""
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import corpus_expand_report as m

    d = _report_fixture()
    base = m.render(d)
    rev = _report_fixture()
    rev["slices"] = list(reversed(rev["slices"]))
    assert m.render(rev) == base, (
        "the slice table's order comes from the stored list, so two tied "
        "sources swap places whenever the query plan changes")


# --- What may be REPUBLISHED is narrower than what may be READ --------------

def _R(lic):
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from corpus_expand_report import _redistributable
    return _redistributable(lic)


def test_a_noncommercial_licence_is_not_redistributable():
    """The defect this pins shipped a wrong number into a published column.

    `_redistributable` tested `startswith`, and `"cc by-nc".startswith("cc by")`
    is True, so every NC licence counted as redistributable -- 237 records
    inside a reported 2,102 on the live crawl, in the very column the page
    tells readers to use INSTEAD of the total.
    """
    for lic in ("cc by-nc", "cc by-nc-sa", "cc by-nc-nd", "CC BY-NC 4.0",
                "cc by-nc/4.0"):
        assert not _R(lic), f"{lic} counted as redistributable"


def test_no_derivatives_is_caught_in_every_spelling():
    """The old ND guard fired only when the string ended exactly in `-nd`, so a
    version suffix or a space walked straight past it."""
    for lic in ("cc by-nd", "cc by-nd/4.0", "CC BY-ND 4.0", "cc-by-nd",
                "cc by-nc-nd/4.0"):
        assert not _R(lic), f"{lic} counted as redistributable"


def test_genuinely_open_licences_still_qualify():
    """A filter that rejects everything would pass both tests above, so the
    permissive direction is checked too."""
    for lic in ("cc by", "cc-by", "CC BY 4.0", "cc0", "CC0 1.0", "cc by-sa",
                "public domain"):
        assert _R(lic), f"{lic} wrongly excluded"


def test_an_absent_or_unknown_licence_is_not_assumed_open():
    for lic in ("", None, "none", "unknown", "all rights reserved", "copyright"):
        assert not _R(lic), f"{lic!r} treated as redistributable"


def test_a_clause_after_a_version_number_is_still_seen():
    """The old tokeniser split on `" 4."` and DISCARDED the rest, so a
    restrictive clause written after the version vanished and the licence was
    published as redistributable."""
    for lic in ("cc by 4.0 nd", "cc by 3.0 nd", "cc by 2.0 nc",
                "CC BY-NonCommercial 4.0", "CC BY-NoDerivatives 4.0"):
        assert not _R(lic), f"{lic} counted as redistributable"


def test_a_licence_url_is_read_rather_than_truncated():
    """Splitting on `/` reduced every CC URL to `https:`, which classified as
    nothing at all -- safe, but it meant the field was simply not read."""
    assert _R("https://creativecommons.org/licenses/by/4.0/")
    assert _R("https://creativecommons.org/licenses/by-sa/4.0/")
    assert _R("https://creativecommons.org/publicdomain/zero/1.0/")
    assert not _R("https://creativecommons.org/licenses/by-nc/4.0/")
    assert not _R("https://creativecommons.org/licenses/by-nc-nd/4.0/")


def test_a_bare_attribution_word_is_not_a_licence():
    """`by` on its own is the English word, not CC BY."""
    assert not _R("by")
    assert not _R("published by elsevier")


# --- --repair must not be able to empty the index ---------------------------

def test_repair_refuses_when_the_shard_directory_is_empty(tmp_path, capsys):
    """THE HAZARD THE REPAIR TOOL INTRODUCED, which was worse than the defect
    it fixes.

    `stored_keys` globs a directory. A directory that is missing, unmounted or
    simply not the one the crawl wrote to yields zero records WITHOUT raising,
    and every stored row then looks orphaned. The first version had no check:
    pointed at an empty directory it deleted the entire `expand-*` index --
    40,814 rows on the live index -- printed a success line and exited 0. And
    the resulting re-downloads are invisible to the duplicate audit, which
    ignores `expand-*` provenance by design.
    """
    from corpus_expand_verify import repair
    idx = _index(tmp_path, [("pmid", str(i), "expand-MED", 0) for i in range(1, 51)])
    empty = tmp_path / "no-shards"
    empty.mkdir()

    assert repair(idx, empty) == 2
    # The EMPTY-read message specifically, not just any refusal. Asserting
    # "REFUSING" alone let the ceiling check satisfy this test, so deleting
    # the empty check outright left it green -- verified by mutation. Two
    # guards covering one case are only two guards if the test can tell them
    # apart, and the empty read is the case that needs its own diagnosis:
    # the ceiling is skipped entirely when the index has no crawl rows yet.
    out = capsys.readouterr().out
    assert "no records were read from the shard directory" in out, out

    c = sqlite3.connect(idx)
    assert c.execute("SELECT COUNT(*) FROM held").fetchone()[0] == 50, (
        "rows were deleted despite the refusal")
    c.close()


def test_repair_refuses_when_the_shard_directory_is_missing(tmp_path):
    from corpus_expand_verify import repair
    idx = _index(tmp_path, [("pmid", "1", "expand-MED", 0)])
    assert repair(idx, tmp_path / "nope") == 2


def test_repair_refuses_when_almost_everything_looks_orphaned(tmp_path, capsys):
    """A real repair touches a small fraction (374 of 40,814 in the incident
    this was written for). A wrong directory touches nearly all of them, so
    the SHAPE of the damage is the signal, not just an empty read."""
    from corpus_expand_verify import repair
    idx = _index(tmp_path, [("pmid", str(i), "expand-MED", 0) for i in range(50)])
    sd = tmp_path / "shards"
    _shard_file(sd, "expanded-00000.jsonl.gz", [{"pmid": "1"}])  # 1 of 50 stored
    assert repair(idx, sd) == 2
    assert "ceiling" in capsys.readouterr().out


def test_repair_proceeds_on_a_genuine_small_repair(tmp_path):
    """The refusals must not block the case the tool exists for."""
    from corpus_expand_verify import repair
    idx = _index(tmp_path, [("pmid", str(i), "expand-MED", 0) for i in range(1, 51)])
    sd = tmp_path / "shards"
    _shard_file(sd, "expanded-00000.jsonl.gz",
                [{"pmid": str(i)} for i in range(1, 49)])  # 2 orphans of 50 = 4%
    assert repair(idx, sd) == 0
    c = sqlite3.connect(idx)
    assert c.execute("SELECT COUNT(*) FROM held").fetchone()[0] == 48
    c.close()


def test_repair_dry_run_deletes_nothing(tmp_path):
    from corpus_expand_verify import repair
    idx = _index(tmp_path, [("pmid", str(i), "expand-MED", 0) for i in range(1, 51)])
    sd = tmp_path / "shards"
    _shard_file(sd, "expanded-00000.jsonl.gz", [{"pmid": str(i)} for i in range(1, 49)])
    assert repair(idx, sd, dry_run=True) == 0
    c = sqlite3.connect(idx)
    assert c.execute("SELECT COUNT(*) FROM held").fetchone()[0] == 50
    c.close()


def test_the_shard_path_follows_the_same_env_var_as_the_crawl(monkeypatch):
    """The one tool that can DELETE data was the only one that could not
    follow a relocation, so pointing the crawl elsewhere armed the wipe."""
    import importlib
    import corpus_expand_verify as v
    monkeypatch.setenv("FERRO_EXPAND_OUT", "/tmp/somewhere-else")
    importlib.reload(v)
    try:
        assert str(v.DEFAULT_SHARDS).startswith("/tmp/somewhere-else")
    finally:
        monkeypatch.delenv("FERRO_EXPAND_OUT", raising=False)
        importlib.reload(v)


def test_repair_refuses_an_empty_read_even_with_no_rows_to_compare_against():
    """The ceiling check is skipped when there are no crawl rows (`if total
    and ...`), so the empty-read refusal is the only thing standing between a
    wrong --shards path and a silent no-op that reports success."""
    import tempfile
    from corpus_expand_verify import repair
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        idx = _index(d, [("pmid", "1", "census-mesh", 0)])  # no expand-* rows
        empty = d / "no-shards"
        empty.mkdir()
        assert repair(idx, empty) == 2


# --- Guards for the mutations that survived round 3 -------------------------
#
# A third review ran nine mutations against these files and SEVEN survived, so
# the previous commit's claim that every new guard was mutation-tested was
# false: the guards that were tested were the ones chosen for testing. Each
# test below exists because a specific mutation stayed green.

def test_a_clause_glued_to_its_version_is_still_a_clause():
    """`cc by-nd4.0` reaches an exact-token test as the opaque word `nd4.0`
    and sails past a check for `nd`."""
    for lic in ("cc by-nd4.0", "CC BY-ND4.0", "cc by-nc4.0", "cc by sa4.0 nd"):
        assert not _R(lic), f"{lic} counted as redistributable"


def test_concatenated_clause_codes_are_decomposed():
    """`CC BY-NCND` is one token to a naive splitter and two clauses in fact."""
    for lic in ("CC BY-NCND", "cc bync", "cc by-ncnd"):
        assert not _R(lic), f"{lic} counted as redistributable"


def test_an_open_word_inside_a_closed_sentence_is_not_a_licence():
    """A positive marker is not enough on its own: these both matched a rule
    that only looked for `zero` or `public domain`."""
    assert not _R("Zero rights granted, all rights reserved")
    assert not _R("not in the public domain")
    assert not _R("no reuse permitted")


def test_the_version_filter_is_load_bearing():
    """Removing `_VERSION` was a surviving mutation. Without it a bare version
    token stays in the set, and `4.0` is not a clause -- the failure it causes
    is that a decomposed clause can no longer be told from a number."""
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import corpus_expand_report as m
    assert m._licence_tokens("cc by 4.0") == {"cc", "by"}, (
        "a bare version number is being kept as a licence clause")
    assert m._licence_tokens("cc0 1.0") == {"cc0"}


def test_the_noise_filter_is_load_bearing():
    """Removing `_NOISE` was a surviving mutation: URL scheme and path words
    would enter the clause set, where `license` and `https` are not clauses."""
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import corpus_expand_report as m
    assert m._licence_tokens("https://creativecommons.org/licenses/by/4.0/") == {
        "cc", "by"}, "URL scaffolding is leaking into the clause set"


def test_the_restrictive_check_runs_before_any_positive_marker():
    """Moving the cc0/zero shortcut above the restrictive check was a
    surviving mutation. A string carrying BOTH must resolve to refusal."""
    assert not _R("cc0 nc")
    assert not _R("publicdomain zero nd")
    assert not _R("cc by-sa nc")


def test_the_repair_ceiling_actually_bounds_what_can_be_deleted(tmp_path):
    """REPAIR_MAX_FRACTION's VALUE was untested: the two existing tests bound
    it only to (0.04, 0.98), so setting it to 0.90 -- permitting deletion of
    ninety per cent of the crawl index -- passed the entire suite."""
    import corpus_expand_verify as v
    assert v.REPAIR_MAX_FRACTION <= 0.25, (
        f"a ceiling of {v.REPAIR_MAX_FRACTION:.0%} permits deleting most of the "
        "index; the incident this tool was built for orphaned 0.9%")
    assert v.REPAIR_MAX_FRACTION > 0.02, (
        "a ceiling this tight refuses the repair the tool exists to perform")

    # And the boundary is exercised, not just asserted: 15% must refuse.
    idx = _index(tmp_path, [("pmid", str(i), "expand-MED", 0) for i in range(1, 101)])
    sd = tmp_path / "shards"
    _shard_file(sd, "expanded-00000.jsonl.gz",
                [{"pmid": str(i)} for i in range(1, 86)])  # 15 orphans of 100
    assert v.repair(idx, sd) == 2


def test_a_missing_shard_directory_is_refused_for_being_missing(tmp_path, capsys):
    """Asserting only `== 2` let the empty-read refusal satisfy this, so
    deleting the is_dir check left it green -- the same vacuity the previous
    commit claimed to have fixed one function up."""
    import corpus_expand_verify as v
    idx = _index(tmp_path, [("pmid", "1", "expand-MED", 0)])
    assert v.repair(idx, tmp_path / "nope") == 2
    assert "is not a directory" in capsys.readouterr().out


def test_force_cannot_delete_against_a_path_that_does_not_exist(tmp_path, capsys):
    """--force restored the round-2 full-wipe: all three refusals sat inside
    one `if not force`, so a typo'd --shards with --force emptied the index
    and exited 0. A path that does not exist cannot be a considered override."""
    import corpus_expand_verify as v
    idx = _index(tmp_path, [("pmid", str(i), "expand-MED", 0) for i in range(1, 51)])
    assert v.repair(idx, tmp_path / "typo", force=True) == 2
    assert "not overridable by --force" in capsys.readouterr().out
    c = sqlite3.connect(idx)
    assert c.execute("SELECT COUNT(*) FROM held").fetchone()[0] == 50
    c.close()


def test_live_shard_survives_a_missing_pgrep(monkeypatch, tmp_path):
    """It used to raise FileNotFoundError and take the whole audit down on an
    image without procps -- for a cosmetic label."""
    import corpus_expand_verify as v

    def boom(*a, **k):
        raise FileNotFoundError("pgrep")

    monkeypatch.setattr(v.subprocess, "run", boom)
    assert v.live_shard(tmp_path) is None
