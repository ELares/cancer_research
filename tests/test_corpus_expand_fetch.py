"""Guards for the expansion fetcher.

Two requirements drive every test here, and both were stated as one sentence
each: never download an article we already have, and do not take paywalled
content. The first is enforced by asking the identity index before fetching
and updating it in the same pass; the second by asking a source that refuses,
rather than by this code deciding what is free.
"""
import gzip
import json
import sqlite3
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import corpus_expand_fetch as fx  # noqa: E402
import corpus_identity_index as ix  # noqa: E402


def test_text_is_asked_for_on_an_open_licence_not_only_the_oa_flag():
    """The flag and the licence are different facts, and the first version
    conflated them.

    `isOpenAccess` means "the full text is in Europe PMC's OA subset". It is
    `N` for about 118,000 cancer records that are nonetheless CC-BY, and every
    one of them was being skipped without being asked for.
    """
    assert fx._text_permitted({"isOpenAccess": "Y"})
    assert fx._text_permitted({"isOpenAccess": "N", "license": "cc by"})
    assert fx._text_permitted({"isOpenAccess": "N", "license": "CC BY-NC-ND"})
    assert fx._text_permitted({"license": "cc0"})
    # And a record with neither is not asked for at all.
    assert not fx._text_permitted({"isOpenAccess": "N", "license": None})
    assert not fx._text_permitted({"isOpenAccess": "N", "license": ""})
    assert not fx._text_permitted({})


def test_availability_and_permission_are_asked_separately():
    """Two different questions, and collapsing them cost a day of 404s.

    Preprints are overwhelmingly CC-BY, so a permission-only gate says yes to
    all 108,000 of them while Europe PMC holds text for barely a sixth. The
    rest answer 404 -- correctly, but only after a request and a sleep each.
    """
    here_and_allowed = {"inEPMC": "Y", "license": "cc by"}
    allowed_not_here = {"inEPMC": "N", "isOpenAccess": "N", "license": "cc by"}
    assert fx._text_available(here_and_allowed) and fx._text_permitted(here_and_allowed)
    assert fx._text_permitted(allowed_not_here), "a CC-BY record is permitted"
    assert not fx._text_available(allowed_not_here), (
        "a record Europe PMC does not hold is still being asked for, which is "
        "a guaranteed 404 per record")
    src = Path(fx.__file__).read_text()
    assert "_text_available(rec) and _text_permitted(rec)" in src, (
        "the fetch is gated on only one of the two questions")


def test_a_restrictive_licence_is_not_read_as_open():
    """The gate must not be a substring free-for-all: 'no licence' or a
    publisher-specific term is not CC."""
    for lic in ("all rights reserved", "copyright acs", "subscription",
                "unknown", "none"):
        assert not fx._text_permitted({"isOpenAccess": "N", "license": lic}), lic


def test_the_source_is_what_refuses_paywalled_text(monkeypatch):
    """This code never decides what is free. It asks, and a 404 is the answer.

    That is the whole reason no publisher page is ever fetched: the endpoint
    that serves open access is the only one asked, and it declines the rest.
    """
    calls = []

    def fake_get(url, tries=4, timeout=120):
        calls.append(url)
        return None  # what _get returns for a 404

    monkeypatch.setattr(fx, "_get", fake_get)
    assert fx.fetch_fulltext({"pmcid": "PMC1", "source": "MED"}) is None
    assert calls and calls[0].endswith("/PMC1/fullTextXML")
    assert "ebi.ac.uk" in calls[0], "text is fetched from somewhere other than Europe PMC"


def test_only_europe_pmc_is_ever_fetched_from():
    """A publisher-page fetch would be scraping, which MISSION.md rules out.

    Checked against the SOURCE rather than trusted: no other host may appear
    as a fetch target anywhere in the module.
    """
    src = Path(fx.__file__).read_text()
    hosts = set()
    for line in src.splitlines():
        if "http" not in line or line.strip().startswith(("#", "*", '"')):
            continue
        for tok in line.split('"'):
            if tok.startswith("http"):
                hosts.add(tok.split("/")[2])
    assert hosts <= {"www.ebi.ac.uk"}, f"fetches from unexpected hosts: {hosts}"


def test_a_record_with_no_identifier_is_refused():
    """Writing a record we could never recognise again would guarantee we
    re-download it, which is the one thing this fetcher must not do."""
    src = Path(fx.__file__).read_text()
    assert 'if not (pmid or pmcid or doi or rec.get("id")):' in src
    assert "continue  # nothing to dedup on later; refuse it" in src


def test_every_written_record_is_indexed_in_the_same_pass():
    """A crash between writing and indexing would re-fetch the record. The
    index write must sit inside the same loop as the shard write."""
    src = Path(fx.__file__).read_text()
    body = src[src.index("shards.write(out)"):src.index("if limit_new and totals")]
    assert "INSERT INTO held" in body, (
        "a record is written to a shard without being added to the index in "
        "the same pass, so an interrupted run will re-download it")
    assert "MAX(held.has_fulltext" in body


def test_shards_rotate_into_complete_files(tmp_path):
    """An interrupted run must leave readable shards, not one truncated file."""
    fx.SHARD_RECORDS = 3
    s = fx.Shards(tmp_path, "t")
    for i in range(7):
        s.write({"i": i})
    s.close()
    files = sorted(tmp_path.joinpath("t").glob("*.jsonl.gz"))
    assert len(files) >= 2
    seen = []
    for f in files:
        with gzip.open(f, "rt") as fh:
            seen += [json.loads(l)["i"] for l in fh]
    assert sorted(seen) == list(range(7)), "records were lost across a rotation"


def test_a_resumed_run_continues_rather_than_restarting(tmp_path, monkeypatch):
    """Shard numbering must resume past what is on disk, or a restart
    overwrites the previous run's output."""
    d = tmp_path / "t"
    d.mkdir()
    (d / "t-00000.jsonl.gz").write_bytes(b"")
    (d / "t-00003.jsonl.gz").write_bytes(b"")
    assert fx.Shards(tmp_path, "t").idx == 4


def test_a_finished_slice_is_not_walked_again():
    src = Path(fx.__file__).read_text()
    assert "if row and row[1]:" in src and "continue" in src
    assert "UPDATE slice SET done=1" in src


def test_a_persistent_http_failure_is_loud():
    """A silent give-up is indistinguishable from a slice with nothing new,
    which is the one confusion this ledger cannot afford."""
    src = Path(fx.__file__).read_text()
    assert "raise RuntimeError" in src, (
        "_get swallows a persistent failure, so a network outage would be "
        "recorded as 'no new records'")
