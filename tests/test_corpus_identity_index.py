"""Guards for the deduplication index.

The user requirement this serves is one sentence -- "every time you download an
article make sure it's new" -- and it is enforced entirely by identifier
NORMALISATION. A DOI that arrives as a URL from one source and a bare string
from another is one article; if the two normalise differently the index says
"new" and the fetcher downloads what it already has. Every test here is
ultimately about that.
"""
import gzip
import json
import sqlite3
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import corpus_identity_index as ix  # noqa: E402


def _mem():
    c = sqlite3.connect(":memory:")
    c.executescript(ix.SCHEMA)
    return c


def test_a_doi_is_the_same_doi_however_it_arrives():
    """The failure this prevents is silent and expensive: the same paper
    downloaded once per spelling of its DOI."""
    forms = [
        "10.1038/S41586-023-06983-9",
        "10.1038/s41586-023-06983-9",
        "https://doi.org/10.1038/s41586-023-06983-9",
        "http://dx.doi.org/10.1038/S41586-023-06983-9",
        "doi:10.1038/s41586-023-06983-9",
        "  10.1038/s41586-023-06983-9.  ",
    ]
    got = {ix.norm_doi(f) for f in forms}
    assert got == {"10.1038/s41586-023-06983-9"}, got


def test_a_non_doi_is_refused_rather_than_stored():
    """Storing junk under the DOI key would make a later lookup miss."""
    for bad in ("", None, "n/a", "10.", "not-a-doi", "PMC12345", "12345"):
        assert ix.norm_doi(bad) is None, bad


def test_a_bare_number_is_never_read_as_a_pmc_id():
    """A bare number is indistinguishable from a PMID, and guessing would
    collide two identifier namespaces -- an article would be 'held' because
    some OTHER article's PMID happened to match its PMC number."""
    assert ix.norm_pmcid("PMC12345") == "PMC12345"
    assert ix.norm_pmcid("pmc12345") == "PMC12345"
    assert ix.norm_pmcid(" PMC12345 ") == "PMC12345"
    assert ix.norm_pmcid("12345") is None
    assert ix.norm_pmcid("PMC") is None
    assert ix.norm_pmid("0") is None
    assert ix.norm_pmid("12a") is None
    assert ix.norm_pmid(" 998 ") == "998"


def test_any_one_key_is_enough_to_call_an_article_held():
    """Sources disagree about which identifiers they carry: a preprint has a
    DOI and no PMID, a full-text shard has a PMC id, a census record may have
    only a PMID. Indexing one key would let the same paper back in."""
    c = _mem()
    c.execute("INSERT INTO held VALUES ('doi','10.1/x','t',0)")
    c.execute("INSERT INTO held VALUES ('pmcid','PMC9','t',0)")
    c.execute("INSERT INTO held VALUES ('pmid','777','t',0)")
    assert ix.is_held(c, doi="https://doi.org/10.1/X")
    assert ix.is_held(c, pmcid="pmc9")
    assert ix.is_held(c, pmid=777)
    # A record sharing NONE of them is new.
    assert not ix.is_held(c, pmid="778", pmcid="PMC10", doi="10.1/y")
    # And an all-empty candidate must not read as held.
    assert not ix.is_held(c, None, None, None)


def test_knowing_an_article_is_not_knowing_we_have_its_text():
    """`held` and `has_fulltext` answer different questions. The census knows
    4.4M articles and holds text for a fraction; a fetcher after TEXT that
    read `held` as `has text` would skip everything."""
    c = _mem()
    c.execute("INSERT INTO held VALUES ('pmid','1','census',0)")
    assert ix.is_held(c, pmid="1")
    assert c.execute("SELECT has_fulltext FROM held WHERE key='1'").fetchone()[0] == 0


def test_learning_about_text_is_never_forgotten():
    """A later source that has only metadata must not clear the flag set by an
    earlier one that had the text."""
    c = _mem()
    rows = [("pmid", "5", "fulltext-store", 1)]
    ix._flush(c, list(rows))
    ix._flush(c, [("pmid", "5", "metadata-only", 0)])
    assert c.execute("SELECT has_fulltext FROM held WHERE key='5'").fetchone()[0] == 1


def test_a_rescan_is_skipped_only_while_the_file_has_not_moved(tmp_path):
    """Keyed on mtime and size: a regenerated shard must be re-read, and an
    untouched one must not be. Keying on the NAME alone would make a
    regenerated shard invisible."""
    c = _mem()
    p = tmp_path / "s.jsonl.gz"
    with gzip.open(p, "wt") as f:
        f.write(json.dumps({"pmid": "42", "pmcid": "PMC42", "doi": "10.1/a"}) + "\n")
    assert ix.scan_jsonl_gz(c, p, "t", fulltext=False, verbose=False) == 1
    assert ix.scan_jsonl_gz(c, p, "t", fulltext=False, verbose=False) == 0
    with gzip.open(p, "wt") as f:
        f.write(json.dumps({"pmid": "43"}) + "\n")
        f.write(json.dumps({"pmid": "44"}) + "\n")
    assert ix.scan_jsonl_gz(c, p, "t", fulltext=False, verbose=False) == 2


def test_a_corrupt_line_does_not_abort_the_shard(tmp_path):
    """One bad line in a million-line shard must cost one record, not the shard.

    The docstring here used to add "and must not mark the shard scanned, which
    would hide the loss". That was FALSE -- the scan `continue`s past the bad
    line, exits the loop normally and records the shard unconditionally -- and
    the test never checked it, so a sentence describing behaviour the code did
    not have sat directly above assertions that could not detect the
    difference.

    Not marking it would also be the wrong fix: the same corrupt file would be
    rescanned on every run forever. The loss is COUNTED instead, so it is
    visible without being re-suffered.
    """
    c = _mem()
    p = tmp_path / "s.jsonl.gz"
    with gzip.open(p, "wt") as f:
        f.write(json.dumps({"pmid": "1"}) + "\n")
        f.write("{ this is not json\n")
        f.write(json.dumps({"pmid": "2"}) + "\n")
    assert ix.scan_jsonl_gz(c, p, "t", fulltext=False, verbose=False) == 2
    assert ix.is_held(c, pmid="1") and ix.is_held(c, pmid="2")

    row = c.execute("SELECT records, bad FROM scanned WHERE path=?",
                    (str(p),)).fetchone()
    assert row is not None, "the shard was not marked scanned, so it will be rescanned forever"
    assert row == (2, 1), (
        f"expected 2 records and 1 unparseable line recorded, got {row}; an "
        "uncounted bad line is a loss nothing can find later")

    # And a clean shard must record ZERO, or the column would be satisfied by
    # any constant and could not distinguish a loss from no loss.
    q = tmp_path / "clean.jsonl.gz"
    with gzip.open(q, "wt") as f:
        f.write(json.dumps({"pmid": "3"}) + "\n")
    ix.scan_jsonl_gz(c, q, "t", fulltext=False, verbose=False)
    assert c.execute("SELECT bad FROM scanned WHERE path=?",
                     (str(q),)).fetchone() == (0,)


def test_the_live_index_actually_covers_the_census():
    """The index is worthless if it is empty or stale, and 'it ran' is not the
    same as 'it holds the corpus'."""
    if not ix.DB.exists():
        pytest.skip("identity index not built in this checkout")
    c = ix.connect()
    s = ix.summary(c)
    assert s["pmid"] > 4_000_000, f"only {s['pmid']:,} PMIDs indexed; the census alone has ~4.4M"
    assert s["doi"] > 1_000_000
    assert s["pmcid"] > 500_000
    assert s["files_scanned"] > 2_000
    c.close()
