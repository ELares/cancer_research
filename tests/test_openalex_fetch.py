"""Guards for the OpenAlex fetcher.

The reason this source is here at all is a measurement that was wrong the first
time: `default.search:cancer` reads FULL TEXT and reported 4.9M PMID-less
works, and sampling them found hop powdery mildew and a paper titled "PEOPLE".
`title_and_abstract.search` is what makes the claim defensible, so the filter
is the first thing tested -- a silent revert would refill the corpus with noise
and nothing else would notice.
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
import openalex_fetch as oaf  # noqa: E402
from corpus_expand_fetch import TransientFetchError  # noqa: E402


def _work(oid, *, doi=None, pmid=None, abstract=None, wtype="article"):
    ids = {"openalex": f"https://openalex.org/{oid}"}
    if doi:
        ids["doi"] = doi
    if pmid:
        ids["pmid"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
    inv = None
    if abstract:
        inv = {}
        for i, word in enumerate(abstract.split()):
            inv.setdefault(word, []).append(i)
    return {"id": ids["openalex"], "ids": ids, "doi": doi, "title": f"T{oid}",
            "publication_year": 2024, "type": wtype, "language": "en",
            "primary_location": {"source": {"display_name": "J"}},
            "open_access": {"is_oa": True, "oa_url": "u"},
            "cited_by_count": 3, "topics": [{"display_name": "Oncology"}],
            "abstract_inverted_index": inv}


def _drive(tmp_path, monkeypatch, pages, limit=None):
    monkeypatch.setattr(oaf, "OUT_ROOT", tmp_path / "out")
    monkeypatch.setattr(ix, "DB", tmp_path / "id.sqlite")
    monkeypatch.setattr(oaf, "SLEEP", 0)
    monkeypatch.setattr(oaf.time, "sleep", lambda *_: None)
    seq = list(pages)
    monkeypatch.setattr(oaf, "_page",
                        lambda cursor, extra="": seq.pop(0) if seq else
                        {"results": [], "meta": {"count": 0, "next_cursor": None}})
    totals = oaf.run(limit=limit, verbose=False)
    stored = []
    for f in sorted((tmp_path / "out").rglob("*.jsonl.gz")):
        try:
            with gzip.open(f, "rt") as fh:
                stored += [json.loads(l)["openalex_id"] for l in fh]
        except EOFError:
            pass
    return totals, stored


def _page_of(works, nxt=None):
    return {"results": works, "meta": {"count": len(works), "next_cursor": nxt}}


def test_the_subject_filter_reads_title_and_abstract_not_full_text():
    """THE MEASUREMENT THAT DECIDED THIS SOURCE EXISTS.

    `default.search` matches anything mentioning cancer once -- sampling it
    returned hop powdery mildew and AIDS activism. Reverting to it would refill
    the corpus with noise while every other test stayed green.
    """
    assert oaf.SUBJECT.startswith("title_and_abstract.search"), oaf.SUBJECT
    assert "default.search" not in oaf.SUBJECT


def test_an_abstract_reconstructs_exactly_from_the_inverted_index():
    text = "Tamoxifen shows good response in advanced breast cancer"
    inv = {}
    for i, w in enumerate(text.split()):
        inv.setdefault(w, []).append(i)
    assert oaf.abstract_text(inv) == text
    assert oaf.abstract_text(None) is None
    assert oaf.abstract_text({}) is None


def test_a_repeated_word_lands_in_every_position():
    """The inverted index maps one word to MANY positions; taking only the
    first would silently shorten every abstract containing a repeat."""
    text = "cancer of the breast and cancer of the lung"
    inv = {}
    for i, w in enumerate(text.split()):
        inv.setdefault(w, []).append(i)
    assert oaf.abstract_text(inv) == text


def test_a_work_already_held_by_doi_is_not_stored_again(tmp_path, monkeypatch):
    """Cross-source dedup: Europe PMC and OpenAlex overlap heavily, and the
    same paper must not be stored once per source."""
    c = sqlite3.connect(tmp_path / "id.sqlite")
    c.executescript(ix.SCHEMA)
    c.execute("INSERT INTO held VALUES ('doi','10.1/known','census-mesh',0)")
    c.commit()
    c.close()
    totals, stored = _drive(tmp_path, monkeypatch, [_page_of([
        _work("W1", doi="https://doi.org/10.1/known"),
        _work("W2", doi="https://doi.org/10.1/new")])])
    assert stored == ["W2"], stored
    assert totals["held_already"] == 1


def test_a_work_with_no_doi_or_pmid_is_still_dedupable(tmp_path, monkeypatch):
    """Conference abstracts and dissertations -- the reason this source is
    worth having -- frequently carry neither. Storing one that cannot be
    recognised again guarantees re-downloading it forever."""
    totals, stored = _drive(tmp_path, monkeypatch, [_page_of([_work("W9")])])
    assert stored == ["W9"]
    c = sqlite3.connect(tmp_path / "id.sqlite")
    held = c.execute("SELECT kind, key FROM held").fetchall()
    c.close()
    assert (oaf.KIND, "W9") in held, held

    # ...and a second run must skip it.
    totals2, stored2 = _drive(tmp_path, monkeypatch, [_page_of([_work("W9")])])
    assert totals2["new"] == 0 and totals2["held_already"] == 1


def test_the_openalex_namespace_does_not_collide_with_article_keys():
    assert oaf.KIND not in ("pmid", "pmcid", "doi", "epmc", "nct")


def test_works_reachable_through_europe_pmc_are_skipped_by_default():
    """A work WITH a PMID is reachable through Europe PMC, which can also give
    us its text. Fetching it here spends requests to arrive at a worse copy."""
    import inspect
    assert "has_pmid:false" in inspect.getsource(oaf.run)
    sig = inspect.signature(oaf.run)
    assert sig.parameters["only_missing_pmid"].default is True


def test_a_paged_search_refusal_is_not_read_as_the_end_of_the_data(monkeypatch):
    monkeypatch.setattr(oaf, "_get", lambda *a, **k: None)
    try:
        oaf._page("*")
    except TransientFetchError:
        return
    raise AssertionError("a refused page returned quietly, truncating the crawl")


def test_the_durability_contract_is_imported_not_copied():
    assert oaf.Shards is fx.Shards
    assert oaf.TransientFetchError is fx.TransientFetchError
    src = Path(oaf.__file__).read_text()
    assert "class Shards" not in src and "def _get(" not in src

