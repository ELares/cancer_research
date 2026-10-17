"""Guards for the duplicate audit.

The audit's whole job is to answer "is every downloaded article new?", and the
two ways it can lie are opposite: counting the crawl's own writes as duplicates
(circularity), or never finding a duplicate at all (a dead predicate). Both
produce clean-looking output, so both are tested.
"""
from __future__ import annotations

import gzip
import json
import pathlib
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
    orphans, truncated, _ = orphaned_index_rows(idx, sd)
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


# --- Torn writes are not truncations (#960) --------------------------------
#
# The tolerance above was written for a shard the crawl is still appending to.
# That shard's deflate stream simply ends, which raises EOFError, and a bad
# header raises gzip.BadGzipFile, an OSError. Both were caught.
#
# A MOUNT THAT DISAPPEARS UNDER AN OPEN HANDLE PRODUCES NEITHER, because an
# SMB client holds written bytes of its own and loses them with the
# connection, leaving whatever was already in those blocks. Measured on the
# #913 recovery walk, three distinct failures come out of that one cause and
# only one of them is an OSError:
#
#   bytes that do not inflate          -> zlib.error        (bare Exception)
#   bytes that inflate but are not text-> UnicodeDecodeError (a ValueError)
#   a stream that just stops           -> EOFError
#
# So `except (EOFError, OSError)` let two of the three through: the audit
# written to survive shard damage died on the first real corruption it met,
# and took `--repair` with it, leaving in place the orphaned index rows that
# repair exists to drop. expanded-01009 raised UnicodeDecodeError after 830
# readable records.


def _damaged(d: Path, name: str, records, tail: bytes) -> Path:
    """A shard whose bytes INFLATE but whose decoded text is not UTF-8.

    This is the real incident's shape and it is the one no amount of gzip
    corruption reproduces: the compressed stream is perfectly valid, so the
    decompressor is happy and the DECODER is what fails. Modelled by writing
    the junk through the compressor rather than appending it to the file.
    """
    p = d / name
    import io
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for r in records:
            gz.write((json.dumps(r) + "\n").encode())
            gz.flush()          # per-record, exactly as the crawl writes
        gz.write(tail)
    p.write_bytes(buf.getvalue())
    return p


# The text decoder reads in blocks, so every record sharing a block with the
# junk is lost with it and earlier blocks decode cleanly. That is why the real
# shard recovered 830 of ~950 rather than 0 or all of them -- and why a
# fixture has to span several blocks to show recovery at all. A 3-record
# fixture puts everything in one block and recovers nothing, which looks like
# the fix failing when it is the fixture being too small.
_BLOCKS = [{"pmid": str(i), "pad": "x" * 200} for i in range(400)]


def _direct_read_error(p: Path) -> tuple[int, str]:
    """What a plain reader actually raises, so the guard's expectation is
    DERIVED rather than a hardcoded guess about CPython's gzip layering.

    The first version of this guard asserted UnicodeDecodeError for a fixture
    that raises BadGzipFile, which tested my prediction instead of the fix.
    """
    n = 0
    try:
        with gzip.open(p, "rt", encoding="utf-8") as fh:
            for _ in fh:
                n += 1
    except Exception as e:
        return n, type(e).__name__
    return n, ""


def test_a_torn_tail_is_recorded_and_its_readable_prefix_recovered(tmp_path):
    d = tmp_path / "shards"
    d.mkdir()
    # 0x89 is what the real incident left behind: inflates, is not UTF-8.
    p = _damaged(d, "expanded-00000.jsonl.gz", _BLOCKS,
                 b"\x89\x50\x4e\x47\xff\xfe")
    direct_n, expected = _direct_read_error(p)
    assert expected == "UnicodeDecodeError", (
        f"this fixture must exercise the decoder arm, not another; got {expected}")
    damage: list = []
    got = [r["pmid"] for r in iter_records(d, damage)]
    assert got == [str(i) for i in range(len(got))], (
        "the recovered records must be the UNBROKEN PREFIX -- a reader that "
        "resynchronised past the damage would silently skip records while "
        "reporting a count")
    assert 0 < len(got) < len(_BLOCKS), (
        f"partial recovery is the whole point: got {len(got)} of {len(_BLOCKS)}")
    assert len(got) == direct_n, (
        "the audit must recover exactly what a plain reader does, no more")
    assert len(damage) == 1 and damage[0][0] == "expanded-00000.jsonl.gz"
    assert damage[0][2] == expected, (
        "the recorded error type must be the one actually raised, so a torn "
        "write stays distinguishable from a truncation in the report")


def test_one_damaged_shard_does_not_hide_the_other_thousand(tmp_path):
    """The regression that matters most. Before the fix the exception escaped
    the generator, so every shard AFTER the damaged one went unread -- and
    `stored_keys` builds the set `--repair` deletes against, so a crash part
    way through would have made the remaining shards' records look orphaned."""
    d = tmp_path / "shards"
    d.mkdir()
    _damaged(d, "expanded-00000.jsonl.gz", _BLOCKS, b"\x89\xff\xfe")
    with gzip.open(d / "expanded-00001.jsonl.gz", "wt") as fh:
        fh.write(json.dumps({"pmid": "LATER"}) + "\n")
    damage: list = []
    got = [r["pmid"] for r in iter_records(d, damage)]
    assert "LATER" in got, (
        "the shard AFTER the damaged one must still be read; before the fix "
        "the exception escaped the generator and every later shard went "
        "unread, which would make their records look orphaned to --repair")
    assert [x[0] for x in damage] == ["expanded-00000.jsonl.gz"]


def test_all_three_storage_failures_are_recorded_and_stay_distinct(tmp_path):
    """One cause, three exception types, two of them outside OSError.

    Collapsing them would make the report say less than it knows: a
    truncated shard is a run still in progress or a killed process, while a
    torn one is storage that went away mid-write. Every expectation here is
    derived from a direct read of the same fixture.
    """
    d = tmp_path / "shards"
    d.mkdir()
    recs = [{"pmid": str(i), "pad": "x" * 40} for i in range(400)]

    # (a) inflates, not text
    a = _damaged(d, "expanded-00000.jsonl.gz", recs, b"\x89\xff\xfe")
    # (b) does not inflate: corrupt the compressed bytes in place
    b = d / "expanded-00001.jsonl.gz"
    with gzip.open(b, "wt", encoding="utf-8") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
            fh.flush()
    raw = bytearray(b.read_bytes())
    for i in range(int(len(raw) * 0.85), int(len(raw) * 0.85) + 8):
        raw[i] ^= 0xFF
    b.write_bytes(bytes(raw))
    # (c) stops early
    c = d / "expanded-00002.jsonl.gz"
    with gzip.open(c, "wt", encoding="utf-8") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
            fh.flush()
    full = c.read_bytes()
    c.write_bytes(full[:int(len(full) * 0.70)])

    expected = {p.name: _direct_read_error(p)[1] for p in (a, b, c)}
    assert len(set(expected.values())) == 3, (
        f"the three fixtures must raise three different types; got {expected}")
    assert not set(expected.values()) & {""}, "every fixture must actually fail"

    damage: list = []
    recovered = list(iter_records(d, damage))
    assert len(damage) == 3, (
        f"all three must be reported rather than raised; got {damage}")
    assert {name: err for name, _, err in damage} == expected, (
        "each shard's recorded type must be the one it actually raised")
    assert len(recovered) >= 3, (
        "the readable prefixes are the point of tolerating damage at all")
    # NON-VACUITY, and the first version of this block WAS vacuous: it
    # resolved each recorded name with `getattr(__builtins__, name)`, and
    # inside an imported module `__builtins__` is a DICT, so every lookup
    # returned None, fell through to `Exception`, and reported False -- even
    # for OSError itself. `oserror_like` was always 0 and the assertion could
    # not fail. Resolved against real types now.
    import builtins
    import zlib
    types = {"error": zlib.error}
    resolved = {n: types.get(n) or getattr(builtins, n) for n in expected.values()}
    assert all(isinstance(t, type) for t in resolved.values()), resolved
    # Self-check the resolver on a name whose answer is known, or a resolver
    # that silently returned the wrong type would make the count meaningless.
    assert issubclass(getattr(builtins, "OSError"), OSError)
    not_oserror = [n for n, t in resolved.items() if not issubclass(t, OSError)]
    assert len(not_oserror) >= 2, (
        "at least two of the three fixtures must raise something OUTSIDE "
        "OSError, or this guard would pass against the original handler and "
        f"prove nothing; resolved {resolved}")


def test_repair_does_not_read_every_shard_twice():
    """STRUCTURAL, because the cost is invisible to any assertion about output.

    `repair` printed `len(stored)` from a SECOND full `stored_keys` pass over
    a corpus that is 51 GB on a network mount -- doubling the wall time of
    the one tool that runs when storage is already known to be unreliable,
    and widening the window in which a drop can interrupt it. The set is
    returned by `orphaned_index_rows` now; this fails if the rescan returns.
    """
    import ast
    src = (Path(__file__).resolve().parent.parent
           / "scripts" / "corpus_expand_verify.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "repair")
    called = {n.func.id for n in ast.walk(fn)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "stored_keys" not in called, (
        "repair must reuse the stored set from orphaned_index_rows rather "
        "than rescanning every shard to print its size")
    assert "orphaned_index_rows" in called


# --- repair must not race a live crawl (#962) ------------------------------


def _repair_fixture(tmp_path):
    """Many stored records and ONE index row that genuinely has no bytes.

    The orphan has to be a small SHARE, not just a small count: with one
    orphan out of two rows the existing REPAIR_MAX_FRACTION ceiling fires
    first and the repair refuses for an unrelated reason, which would make
    every guard below pass without exercising anything. 1 of 21 is 4.8%.
    """
    stored = [{"pmid": str(i)} for i in range(1, 21)]
    sd = _shards(tmp_path, stored)
    rows = [("pmid", r["pmid"], "expand-MED", 0) for r in stored]
    rows.append(("pmid", "999", "expand-MED", 0))
    return sd, _index(tmp_path, rows)


def test_repair_refuses_while_a_crawl_is_running(tmp_path, monkeypatch, capsys):
    """The delete is gated, not the scan.

    `verify` already consulted `live_shard`; `repair` -- the only function that
    DELETES -- did not, so correctness depended on remembering to sequence the
    two by hand.
    """
    import corpus_expand_verify as cv
    sd, idx = _repair_fixture(tmp_path)
    monkeypatch.setattr(cv, "live_shard", lambda d: "expanded-00007.jsonl.gz")
    rc = cv.repair(idx, sd)
    out = capsys.readouterr().out
    assert rc == 2, "a running crawl must stop the repair"
    assert "REFUSING" in out and "expanded-00007.jsonl.gz" in out, (
        "the refusal must name the shard being written, so an operator can "
        f"check it against the process table; got {out!r}")
    c = sqlite3.connect(idx)
    try:
        assert c.execute("SELECT COUNT(*) FROM held").fetchone()[0] == 21, (
            "the refusal must delete NOTHING")
    finally:
        c.close()


def test_repair_proceeds_when_no_crawl_is_running(tmp_path, monkeypatch):
    """The control. Without this the guard above is satisfied by a repair that
    never works at all."""
    import corpus_expand_verify as cv
    sd, idx = _repair_fixture(tmp_path)
    monkeypatch.setattr(cv, "live_shard", lambda d: None)
    assert cv.repair(idx, sd) == 0
    c = sqlite3.connect(idx)
    try:
        rows = [r[0] for r in c.execute("SELECT key FROM held")]
    finally:
        c.close()
    assert "999" not in rows, "the row with no stored record must go"
    assert len(rows) == 20, f"every backed row must survive; got {len(rows)}"


def test_force_overrides_the_live_crawl_refusal(tmp_path, monkeypatch):
    """Overridable, like the other two refusals that describe a recoverable
    situation -- an operator who has just stopped a crawl and checked the
    process table can say so. The non-directory refusal stays absolute."""
    import corpus_expand_verify as cv
    sd, idx = _repair_fixture(tmp_path)
    monkeypatch.setattr(cv, "live_shard", lambda d: "expanded-00007.jsonl.gz")
    assert cv.repair(idx, sd, force=True) == 0


def test_the_delete_is_gated_and_the_scan_is_not(tmp_path, monkeypatch):
    """STRUCTURAL. `verify` must stay runnable against a live crawl -- it is
    the tool for watching one -- so the gate belongs in `repair` alone. A fix
    that put it in `iter_records` or `stored_keys` would break that and still
    pass every behavioural guard above.
    """
    import ast
    src = (Path(__file__).resolve().parent.parent
           / "scripts" / "corpus_expand_verify.py").read_text()
    tree = ast.parse(src)

    def calls(name):
        fn = next(n for n in ast.walk(tree)
                  if isinstance(n, ast.FunctionDef) and n.name == name)
        return {n.func.id for n in ast.walk(fn)
                if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}

    assert "live_shard" in calls("repair"), "the delete must consult it"
    for scanner in ("iter_records", "stored_keys"):
        assert "live_shard" not in calls(scanner), (
            f"{scanner} must stay usable while a crawl runs; gating the scan "
            "would break `verify`, whose job is to report on a live crawl")


def test_repair_records_what_it_deletes_before_deleting_it(tmp_path, monkeypatch):
    """An orphan row is the ONLY record of a lost article.

    The bytes are gone from the shards -- that is what makes the row an orphan
    -- so once the row goes too, nothing anywhere names the article. Deleting
    restores eligibility to re-fetch, which pays off only when the slice
    holding it is walked again, and a slice already marked `done` never is. So
    the dump is what makes the repair recoverable rather than merely tidy.
    """
    import corpus_expand_verify as cv
    sd, idx = _repair_fixture(tmp_path)
    monkeypatch.setattr(cv, "live_shard", lambda d: None)
    assert cv.repair(idx, sd) == 0
    dumps = sorted(idx.parent.glob("orphans-*.json"))
    assert len(dumps) == 1, f"exactly one dump expected; got {dumps}"
    rec = json.loads(dumps[0].read_text())
    assert rec["deleted"] == [{"kind": "pmid", "key": "999"}], (
        f"the dump must name the rows actually deleted; got {rec['deleted']}")
    assert rec["index"] == str(idx) and rec["shard_dir"] == str(sd), (
        "a dump that does not say which index and shard directory it describes "
        "cannot be acted on later")


def test_repair_deletes_nothing_when_it_cannot_record_the_rows(tmp_path, monkeypatch):
    """A PRECONDITION, not a side effect. Recording the loss after destroying
    the evidence is the shape this repo already rejected for provenance."""
    import corpus_expand_verify as cv
    sd, idx = _repair_fixture(tmp_path)
    monkeypatch.setattr(cv, "live_shard", lambda d: None)

    real = pathlib.Path.write_text

    def refuse(self, *a, **k):
        if self.name.startswith("orphans-"):
            raise OSError("no space left on device")
        return real(self, *a, **k)

    monkeypatch.setattr(pathlib.Path, "write_text", refuse)
    assert cv.repair(idx, sd) == 2, "an unwritable dump must stop the delete"
    c = sqlite3.connect(idx)
    try:
        assert c.execute("SELECT COUNT(*) FROM held").fetchone()[0] == 21, (
            "the orphan row must survive a repair that could not record it")
    finally:
        c.close()


def test_the_dump_precedes_the_delete_in_source_order():
    """STRUCTURAL. Both behavioural guards above pass if the dump is written
    AFTER the delete, because they only observe the end state -- and a dump
    written afterwards is exactly the ordering that loses the rows when the
    process dies between the two."""
    import ast
    src = (Path(__file__).resolve().parent.parent
           / "scripts" / "corpus_expand_verify.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "repair")
    dump_line = delete_line = None
    for node in ast.walk(fn):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr == "write_text" and dump_line is None:
                dump_line = node.lineno
            if node.func.attr == "executemany":
                delete_line = node.lineno
    assert dump_line and delete_line, (dump_line, delete_line)
    assert dump_line < delete_line, (
        f"the dump (line {dump_line}) must be written before the DELETE "
        f"(line {delete_line}), or a crash between them loses the rows")


def test_the_orphan_dump_is_not_committable():
    """The dump lands beside the index, and everything else beside the index is
    gitignored -- the ledgers, the raw shards, the manifests. The dump was not,
    so a `git add -A` would have committed per-run operational data into a
    public repository. The pattern is DERIVED from the writer's own format
    string, so moving or renaming the dump fails here rather than silently
    leaving it committable.
    """
    import re
    import subprocess
    src = (Path(__file__).resolve().parent.parent
           / "scripts" / "corpus_expand_verify.py").read_text()
    m = re.search(r'f"(orphans-\{[a-z_]+\})\.json"', src)
    assert m, "the dump filename is no longer an f-string this guard can read"
    sample = Path("corpus/atlas") / (m.group(1).split("{")[0] + "20260101-000000.json")
    repo = Path(__file__).resolve().parent.parent
    r = subprocess.run(["git", "check-ignore", "-q", str(sample)],
                       cwd=repo, capture_output=True)
    if r.returncode == 128:
        pytest.skip("not a git work tree")
    assert r.returncode == 0, (
        f"{sample} is not gitignored, so a repair's orphan dump can be "
        "committed by accident; add a pattern covering it to .gitignore")
