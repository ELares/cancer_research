"""Guards for the retrieval index.

The requirement this serves is that a paper the corpus HOLDS can be FOUND.
That was measured to be false: 10 of 11 known-relevant papers on bacterial
cancer immunotherapy were already in the corpus, and a naive query over the
corpus returned none of them. Every test here defends one step of the path
from "the bytes are on disk" to "the paper is in the top ten".

Four of these tests exist because the first version of the indexer failed
them, and each failure was silent -- it produced an index, reported a count,
and lost or hid data without erroring.
"""
import gzip
import json
import sqlite3
import sys
import zlib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import corpus_retrieval_index as ri  # noqa: E402
import corpus_retrieval_search as rs  # noqa: E402


def _db(tmp_path) -> sqlite3.Connection:
    return ri.connect(tmp_path / "t.sqlite")


def _shard(path: Path, records: list[dict]) -> None:
    with gzip.open(path, "wt") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")


def _untrailered_shard(path: Path, records: list[dict]) -> None:
    """A gzip stream flushed with Z_SYNC_FLUSH and never closed.

    This is not a synthetic corruption -- it is exactly what the crawler's
    durability fix leaves on disk, because bytes are flushed so the index
    cannot claim a record before it is readable, and the trailer is only
    written on close(). 649 real shards are in this state.
    """
    co = zlib.compressobj(6, zlib.DEFLATED, 31)
    buf = b""
    for r in records:
        buf += co.compress((json.dumps(r) + "\n").encode())
    buf += co.flush(zlib.Z_SYNC_FLUSH)
    path.write_bytes(buf)


def test_markup_does_not_break_a_phrase_that_is_plainly_there():
    """THE defect that made gene-name queries fail.

    Publishers italicise gene symbols, so the highest-value tokens in cancer
    literature are the ones most often wrapped in markup. FTS5's unicode61
    tokeniser treats `<`, `>` and `/` as separators, so `<i>KRAS</i> -mutant`
    becomes the stream `i · kras · i · mutant` and the stray `i` sits between
    the two words a phrase query needs adjacent.

    Measured over the real corpus: 5.9% of titles and 53.7% of abstracts
    carry markup.
    """
    raw = "Efficacy in <i>KRAS</i> -mutant Pancreatic Cancer"
    c = sqlite3.connect(":memory:")
    c.execute("CREATE VIRTUAL TABLE d USING fts5(t, tokenize='porter unicode61 remove_diacritics 2')")
    c.execute("INSERT INTO d(t) VALUES (?)", (raw,))
    assert c.execute("SELECT COUNT(*) FROM d WHERE d MATCH '\"KRAS mutant\"'").fetchone()[0] == 0, (
        "if this passes, markup no longer breaks adjacency and this guard is vacuous")
    c.execute("DELETE FROM d")
    c.execute("INSERT INTO d(t) VALUES (?)", (ri.clean_text(raw),))
    assert c.execute("SELECT COUNT(*) FROM d WHERE d MATCH '\"KRAS mutant\"'").fetchone()[0] == 1
    c.close()


def test_entities_are_unescaped_before_tags_are_stripped():
    """The same markup arrives as `<i>` from one source and `&lt;i&gt;` from
    another, and some records are doubly escaped. Stripping tags first would
    leave the escaped form untouched and index the literal text `lt i gt`."""
    assert ri.clean_text("&lt;i&gt;In silico&lt;/i&gt; design") == "In silico design"
    assert ri.clean_text("&amp;lt;i&amp;gt;BRAF&amp;lt;/i&amp;gt; V600E") == "BRAF V600E"
    assert ri.clean_text("<sub>2</sub>O chemistry") == "2 O chemistry"


def test_a_shard_with_no_gzip_trailer_keeps_every_readable_record(tmp_path):
    """The first version LOST these records.

    A shard holding fewer records than one BATCH was parsed in full, every key
    added to the dedup set, and then discarded when EOFError skipped the flush
    at the end of the loop. The real shard that exposed it holds 392 records.
    """
    recs = [{"pmid": str(9000 + i), "title": f"paper {i}", "abstract": "x"} for i in range(392)]
    p = tmp_path / "trunc.jsonl.gz"
    _untrailered_shard(p, recs)
    with pytest.raises(EOFError):
        with gzip.open(p, "rt") as f:
            for _ in f:
                pass

    c = _db(tmp_path)
    added, bad = ri.index_shard(c, p, set(), verbose=False)
    assert added == 392, f"kept {added} of 392 readable records from an untrailered shard"
    assert bad == 0
    assert c.execute("SELECT COUNT(*) FROM doc").fetchone()[0] == 392
    c.close()


def test_an_untrailered_shard_is_not_rescanned_forever(tmp_path):
    """The first version refused to mark these shards scanned, reasoning that
    the rest should be picked up later. `scanned` is keyed on mtime AND size,
    so an extended shard is re-read anyway -- refusing to mark it only
    guaranteed a permanent re-read of a 49 MB file over the network, on every
    run, for every one of the 649 shards in this state."""
    recs = [{"pmid": str(7000 + i), "title": f"t{i}"} for i in range(10)]
    p = tmp_path / "live.jsonl.gz"
    _untrailered_shard(p, recs)
    c = _db(tmp_path)

    first, _ = ri.index_shard(c, p, set(), verbose=False)
    assert first == 10

    # THE CLAIM, asserted directly. An earlier version of this test asserted
    # that a second pass added no records, which a mutation removing the
    # `scanned` write SURVIVED: dedup already suppresses the records, so
    # "nothing was added" is true whether the shard was skipped or read in
    # full. The property is that the shard is MARKED, so it is never opened
    # again -- and only the row can witness that.
    row = c.execute("SELECT records, bad FROM scanned WHERE path=?",
                    (str(p),)).fetchone()
    assert row == (10, 0), (
        f"untrailered shard left unmarked (row={row}); it will be re-read on "
        "every run, and the real ones are 49 MB each over the network")

    # And the marking must not be permanent: a shard the crawler EXTENDS has a
    # new size, so it is re-read and its tail is picked up. Checked by reading
    # with an EMPTY seen set, so the re-read is visible in the return value
    # rather than being masked by dedup.
    _untrailered_shard(p, recs + [{"pmid": "7999", "title": "appended"}])
    grew, _ = ri.index_shard(c, p, set(), verbose=False)
    assert grew == 11, (
        f"extended shard yielded {grew}, expected all 11 re-read; a shard "
        "frozen at its first size loses every record appended after it")
    c.close()


def test_connect_honours_a_redirected_db(tmp_path):
    """`def connect(db=DB)` bound the default at definition time, so a test
    that reassigned the module global wrote 18,000 rows into the REAL index
    while reporting that it had used a scratch file."""
    ri_db = tmp_path / "redirected.sqlite"
    c = ri.connect(ri_db)
    c.execute("INSERT INTO doc(title) VALUES ('x')")
    c.commit()
    c.close()
    assert ri_db.exists(), "connect() ignored the path it was given"
    assert sqlite3.connect(str(ri_db)).execute(
        "SELECT COUNT(*) FROM doc").fetchone()[0] == 1


def test_the_fulltext_summary_can_exceed_zero(tmp_path):
    """`has_text` is stored as an INTEGER and the summary queried `='1'`.
    SQLite does not coerce INTEGER to TEXT, so the count was 0 over a corpus
    holding 1.5 million full texts. A summary that cannot rise above zero
    reports nothing, and would have been read as "no full text was indexed"."""
    p = tmp_path / "s.jsonl.gz"
    _shard(p, [{"pmid": "1", "title": "a", "has_text": True},
               {"pmid": "2", "title": "b", "has_text": False}])
    c = _db(tmp_path)
    ri.index_shard(c, p, set(), verbose=False)
    s = ri.summary(c)
    assert s["with_fulltext"] == 1, f"expected 1 full text, summary said {s['with_fulltext']}"
    assert s["docs"] == 2
    c.close()


def test_one_paper_from_two_sources_takes_one_result_slot(tmp_path):
    """A duplicate does not merely waste a row. It occupies one of the ten
    slots a human reads, so it DISPLACES a different paper from the answer."""
    a = tmp_path / "a.jsonl.gz"
    b = tmp_path / "b.jsonl.gz"
    _shard(a, [{"pmid": "555", "doi": "10.1/x", "title": "Engineered bacteria in tumours"}])
    # Same paper, arriving from another source with the DOI spelled as a URL
    # and no PMID at all -- the shape an OpenAlex record actually has.
    _shard(b, [{"doi": "https://doi.org/10.1/X", "title": "Engineered bacteria in tumours"}])
    c = _db(tmp_path)
    seen: set[str] = set()
    ri.index_shard(c, a, seen, verbose=False)
    ri.index_shard(c, b, seen, verbose=False)
    # The DOI-only record keys on the DOI; the first keyed on its PMID, so
    # these are two different keys and BOTH are indexed. This is the
    # documented limit of key-order dedup, asserted so it cannot change
    # silently -- not a claim that it is ideal.
    n = c.execute("SELECT COUNT(*) FROM doc").fetchone()[0]
    assert n == 2, (
        f"expected the known key-order limitation (2 rows), got {n}; if dedup "
        "improved to catch this, update the test and the index docstring")
    c.close()


def test_a_record_with_no_text_is_not_indexed(tmp_path):
    """`fulltext` shards carry {pmcid, pmid, text, licence} and no title. A row
    with neither title nor abstract can only ever add noise to a ranking."""
    assert ri.to_row({"pmid": "1", "text": "body only"}) is None
    assert ri.to_row({"pmid": "1", "title": "", "abstract": "   "}) is None
    assert ri.to_row({"pmid": "1", "abstract": "real abstract"}) is not None


def test_missing_citation_data_is_not_stored_as_never_cited(tmp_path):
    """Europe PMC records here carry no citation count. Storing 0 would make
    2.5 million papers look uncited next to 18,938 OpenAlex ones, and any
    later popularity ranking would read that as a measurement."""
    i = ri.COLUMNS.index("cited_by")   # by NAME: the position moved once
    assert ri.to_row({"pmid": "1", "title": "t"})[i] == -1
    assert ri.to_row({"pmid": "1", "title": "t", "cited_by": 0})[i] == 0
    assert ri.to_row({"pmid": "1", "title": "t", "cited_by": 42})[i] == 42


def test_better_matches_come_first(tmp_path):
    """bm25() returns a NEGATIVE score where more negative is better, so the
    ordering must be ASCENDING. Descending would return the WORST matches in
    the pool, which looks plausible until you read the results."""
    p = tmp_path / "s.jsonl.gz"
    _shard(p, [
        {"pmid": "1", "title": "Engineered bacteria for cancer immunotherapy",
         "abstract": "engineered bacteria"},
        {"pmid": "2", "title": "Unrelated survey of hospital logistics",
         "abstract": "a passing mention of bacteria somewhere in the text"},
    ])
    c = _db(tmp_path)
    ri.index_shard(c, p, set(), verbose=False)
    rows = rs.search(c, rs.phrase_query(["bacteria"]), rs.W_TITLE_HEAVY, limit=5)
    assert rows[0]["pmid"] == "1", (
        f"ranking put pmid={rows[0]['pmid']} first; the title match must win")
    assert rows[0]["score"] <= rows[-1]["score"], "scores are not ascending"
    c.close()


def test_a_phrase_is_narrower_than_the_same_words_loose(tmp_path):
    """The whole premise. If a phrase query did not cut the pool, there would
    be nothing to build: the measured comparison was 18,766 records against 8
    on the same subject."""
    p = tmp_path / "s.jsonl.gz"
    _shard(p, [
        {"pmid": "1", "title": "Engineered bacteria tag tumours"},
        {"pmid": "2", "title": "Bacteria in soil", "abstract": "engineered crops"},
        {"pmid": "3", "title": "Tumour staging", "abstract": "bacteria unrelated"},
    ])
    c = _db(tmp_path)
    ri.index_shard(c, p, set(), verbose=False)
    loose = rs.count(c, rs.bag_of_words("engineered bacteria"))
    phrase = rs.count(c, rs.phrase_query(["engineered bacteria"]))
    assert phrase < loose, f"phrase pool {phrase} did not narrow the loose pool {loose}"
    assert phrase == 1
    c.close()


def test_fts_syntax_in_a_question_does_not_crash_the_search(tmp_path):
    """A user question containing `AND`, a quote or a `*` is a question, not
    FTS5 syntax. An unescaped one is a syntax error, so the search would fail
    rather than return results."""
    # Against the REAL schema. An earlier version built a two-column synthetic
    # table, which stopped working the moment the queries became
    # column-restricted -- the fixture, not the code, was wrong, and a fixture
    # that does not match the schema cannot test the query it is given.
    p = tmp_path / "s.jsonl.gz"
    _shard(p, [{"pmid": "1", "title": "bacteria tumor", "abstract": "x"}])
    c = _db(tmp_path)
    ri.index_shard(c, p, set(), verbose=False)
    for q in ['bacteria AND "tumor"', "what about NOT bacteria?", "bacteria*",
              "(engineered) bacteria", "p53 / KRAS", "{title}: injection"]:
        m = rs.bag_of_words(q)
        # Must not raise. The count itself is not the subject; surviving the
        # user's punctuation is.
        rs.count(c, m)
    c.close()


def test_the_evaluation_report_order_is_established_by_the_renderer():
    """The compensating guard for `corpus_retrieval_eval`'s freshness exemption.

    That generator reads a ~5 GB gitignored index built from shards on external
    storage, so `test_artifact_freshness.py` cannot regenerate it and exempts
    it -- and exemption removes it from that file's LIVE list, which is what
    its order gate iterates. So the property that gate would have checked is
    checked HERE, where it needs neither the index nor a committed artifact:
    `render()` is handed a synthetic result dict twice, with the mapping keys
    in opposite orders, and must produce byte-identical output.

    The defect this prevents has happened twice in this repo: `sort_keys=True`
    reordered a stored dict and `--render-only` silently regenerated a
    DIFFERENT document, in one case dropping the finding the page existed to
    report.
    """
    import corpus_retrieval_eval as ev

    def spec(order):
        qs = []
        for qid in order:
            qs.append({
                "id": qid, "question": f"q {qid}", "concepts": ["x"],
                "note": "", "gold": ["1"],
                "a": {"arm": "a", "match": "x", "pool": 10, "ranks": {"1": 3},
                      "found": 1, "gold_total": 1, "recall_at_k": 1.0, "rr": 0.33,
                      "top3": [{"pmid": "1", "title": "t"}]},
                "b": {"arm": "b", "match": "x", "pool": 2, "ranks": {"1": 1},
                      "found": 1, "gold_total": 1, "recall_at_k": 1.0, "rr": 1.0,
                      "top3": [{"pmid": "1", "title": "t"}]},
            })
        return {
            "docs_indexed": 100, "queries": len(qs), "gold_total": len(qs),
            "depth": 50, "at_k": 10,
            "aggregate": {
                "a": {"mrr": 0.33, "recall_at_k": 1.0, "gold_found": 1, "median_pool": 10},
                "b": {"mrr": 1.0, "recall_at_k": 1.0, "gold_found": 1, "median_pool": 2},
            },
            "per_query": qs, "provenance": "p", "known_bias": "b",
        }

    first = ev.render(spec(["alpha", "beta"]))
    # Same queries, reversed. The renderer must impose the order it was given,
    # not inherit whatever order the mapping happened to carry.
    second = ev.render(spec(["beta", "alpha"]))
    assert first != second, (
        "render() produced identical output for two different query orders, so "
        "it is not reading the order it was given and this guard is vacuous")
    assert ev.render(spec(["alpha", "beta"])) == first, "render() is not deterministic"
    assert first.index("alpha") < first.index("beta")
    assert second.index("beta") < second.index("alpha")


def test_the_evaluation_states_its_verdict_from_the_numbers():
    """A report whose headline is a fixed string cannot deliver a negative
    result. This repo shipped exactly that in `manuscript_vs_census.py`: the
    conclusion was a static sentence saying the claim survived, and forcing the
    comparison to fail produced a document whose body contradicted its own
    opening line. So the verdict sentence must FLIP when arm B loses."""
    import corpus_retrieval_eval as ev

    def spec(a_mrr, b_mrr):
        return {
            "docs_indexed": 100, "queries": 1, "gold_total": 1, "depth": 50, "at_k": 10,
            "aggregate": {
                "a": {"mrr": a_mrr, "recall_at_k": 0.5, "gold_found": 1, "median_pool": 10},
                "b": {"mrr": b_mrr, "recall_at_k": 0.5, "gold_found": 1, "median_pool": 2},
            },
            "per_query": [], "provenance": "p", "known_bias": "b",
        }

    won = ev.render(spec(0.2, 0.9))
    lost = ev.render(spec(0.9, 0.2))
    assert "does not beat" not in won
    assert "does not beat" in lost, (
        "the report cannot report its own negative result; the verdict is not "
        "derived from the measured MRRs")


def test_the_row_shape_matches_the_declared_columns():
    """`to_row` returns a positional tuple and `_flush` inserts it against a
    hand-written column list. Inserting `mesh` at index 2 moved `cited_by` from
    8 to 9 and a test asserting index 8 kept passing against the wrong field.
    The declared names are the contract now, and this pins them to the tuple."""
    row = ri.to_row({"pmid": "1", "title": "t", "abstract": "a",
                     "mesh": ["Neoplasms"], "cited_by": 7})
    assert len(row) == len(ri.COLUMNS), (
        f"to_row returns {len(row)} values against {len(ri.COLUMNS)} declared "
        f"columns {ri.COLUMNS}")
    named = dict(zip(ri.COLUMNS, row))
    assert named["title"] == "t"
    assert named["mesh"] == "Neoplasms"
    assert named["cited_by"] == 7
    assert named["pmid"] == "1"


def test_mesh_is_indexed_but_never_merged_into_the_abstract(tmp_path):
    """MeSH descriptors are NLM's expert subject labels and so the strongest
    retrieval signal available -- which is exactly why they must stay in their
    own column. `mechanism_recall.py` exists because folding MeSH into the
    matched text makes a keyword measurement circular: the keyword gets scored
    against a label that already encodes the answer."""
    p = tmp_path / "s.jsonl.gz"
    _shard(p, [{"pmid": "1", "title": "A paper", "abstract": "body text",
                "mesh": ["Ferroptosis", "Neoplasms"]}])
    c = _db(tmp_path)
    ri.index_shard(c, p, set(), verbose=False)
    row = c.execute("SELECT title, abstract, mesh FROM doc").fetchone()
    assert "Ferroptosis" not in row[0] and "Ferroptosis" not in row[1], (
        "a MeSH descriptor leaked into the title or abstract, which makes any "
        "keyword measurement over those fields circular")
    assert "Ferroptosis" in row[2]
    # And it must be reachable when asked for deliberately.
    assert c.execute(
        "SELECT COUNT(*) FROM doc WHERE doc MATCH 'mesh:Ferroptosis'").fetchone()[0] == 1
    # ...but NOT when the query is restricted to the text fields.
    assert c.execute(
        "SELECT COUNT(*) FROM doc WHERE doc MATCH '{title abstract}:Ferroptosis'"
    ).fetchone()[0] == 0
    c.close()


def test_the_keyword_report_order_is_established_by_the_renderer():
    """Compensating guard for `mechanism_keyword_precision`'s exemption, for
    the same reason as the evaluation's: it reads the gitignored index, so the
    freshness gate cannot regenerate it, and exemption removes it from that
    file's LIVE list."""
    import mechanism_keyword_precision as mk

    def spec(order):
        return {
            "docs_indexed": 100, "seed": 1, "depth": 10, "aboutness_floor": 0.3,
            "candidates": [
                {"phrase": p, "control": False, "pool": 5, "inspected": 5,
                 "microbe_in_title": 0.8, "both_in_title": 0.6,
                 "verdict": "admit", "example_titles": []}
                for p in order
            ],
            "sample": [], "admitted": list(order), "refused": [],
            "control_separates": True, "worst_real": 0.6, "best_control": 0.1,
        }

    first = mk.render(spec(["alpha", "beta"]))
    second = mk.render(spec(["beta", "alpha"]))
    assert first != second, (
        "render() ignored the candidate order it was given, so this guard is "
        "vacuous")
    assert first.index("alpha") < first.index("beta")
    assert second.index("beta") < second.index("alpha")


def test_the_keyword_report_refuses_to_endorse_a_measure_that_does_not_separate():
    """The controls exist to calibrate the measure. If a generic word scores
    like a specific phrase the column means nothing, and the page must say so
    rather than printing a verdict table under a confident heading. This repo
    shipped a static endorsement once (`manuscript_vs_census.py`) and forcing
    the comparison to fail produced a document contradicting its own opening
    line."""
    import mechanism_keyword_precision as mk

    def spec(separates):
        return {
            "docs_indexed": 100, "seed": 1, "depth": 10, "aboutness_floor": 0.3,
            "candidates": [], "sample": [], "admitted": ["x"], "refused": [],
            "control_separates": separates,
            "worst_real": 0.1 if not separates else 0.6,
            "best_control": 0.5 if not separates else 0.1,
        }

    good = mk.render(spec(True))
    bad = mk.render(spec(False))
    assert "should be acted on" not in good
    assert "should be acted on" in bad, (
        "the page endorses its own table even when the controls show the "
        "measure does not discriminate")


def test_a_query_cannot_reach_mesh_by_default(tmp_path):
    """Zeroing the MeSH bm25 weight is NOT sufficient, and assuming it was is
    the defect this catches.

    FTS5 ranking and FTS5 matching are separate. A weight of 0 removes the
    column from the SCORE, but an unrestricted query still MATCHES on it, so
    the POOL grows and a record whose only occurrence of a term is in its
    expert label counts as a hit. Measured on two synthetic rows: the bare
    query matched the MeSH-only row, the restricted one did not.

    This matters because `mechanism_keyword_precision` scores a candidate
    keyword against the corpus, and MeSH is NLM's answer to the question that
    measurement is asking.
    """
    p = tmp_path / "s.jsonl.gz"
    _shard(p, [
        {"pmid": "1", "title": "Engineered bacteria in tumours",
         "abstract": "body", "mesh": ["Ferroptosis", "Neoplasms"]},
        # The trap: `ferroptosis` appears ONLY in the expert label.
        {"pmid": "2", "title": "An unrelated paper", "abstract": "nothing",
         "mesh": ["Ferroptosis"]},
    ])
    c = _db(tmp_path)
    ri.index_shard(c, p, set(), verbose=False)

    unrestricted = c.execute(
        "SELECT COUNT(*) FROM doc WHERE doc MATCH ?", ('"ferroptosis"',)).fetchone()[0]
    assert unrestricted == 2, (
        "an unrestricted query no longer reaches MeSH, so this guard is "
        "vacuous and the restriction below proves nothing")

    # The module's own query builders must not be reachable that way.
    assert rs.count(c, rs.phrase_query(["ferroptosis"])) == 0, (
        "phrase_query matched a record whose only occurrence is its MeSH "
        "label; the query is not column-restricted")
    assert rs.count(c, rs.bag_of_words("ferroptosis")) == 0, (
        "bag_of_words is not column-restricted")
    # And the restriction must not break ordinary matching.
    assert rs.count(c, rs.phrase_query(["engineered bacteria"])) == 1
    # MeSH stays reachable when a caller asks for it ON PURPOSE.
    assert rs.count(c, rs.phrase_query(["ferroptosis"], fields=("mesh",))) == 2
    c.close()


def test_the_mesh_weight_is_zero_in_every_measurement_arm():
    """A positional weight tuple shorter than the indexed columns leaves the
    rest at FTS5's default of 1.0, which is how MeSH silently entered both
    arms' ranking when the column was added. Both arms must name it."""
    indexed = [c for c in ri.COLUMNS if c in ("title", "abstract", "mesh")]
    for name, w in (("W_UNIFORM", rs.W_UNIFORM), ("W_TITLE_HEAVY", rs.W_TITLE_HEAVY)):
        assert len(w) == len(indexed), (
            f"{name} has {len(w)} weights against {len(indexed)} indexed "
            f"columns {indexed}; the unnamed ones default to 1.0")
        assert w[indexed.index("mesh")] == 0.0, f"{name} gives MeSH weight {w}"


def test_the_verdict_cannot_claim_a_win_while_coverage_falls():
    """`check the column beside it`, as a guard.

    The first version of this report compared MRR alone and printed "arm B
    reaches the known paper where arm A does not" while the column immediately
    beside it recorded arm B reaching FOUR gold papers against arm A's five.
    A ranking win and a coverage loss are different facts.
    """
    import corpus_retrieval_eval as ev

    def spec(a_mrr, b_mrr, a_found, b_found):
        return {
            "docs_indexed": 100, "queries": 1, "gold_total": 5, "depth": 50,
            "at_k": 10,
            "aggregate": {
                "a": {"mrr": a_mrr, "recall_at_k": 0.5, "gold_found": a_found,
                      "median_pool": 1000},
                "b": {"mrr": b_mrr, "recall_at_k": 0.7, "gold_found": b_found,
                      "median_pool": 30},
            },
            "per_query": [{
                "id": "q", "question": "x", "concepts": ["c"], "note": "",
                "gold": ["1"],
                "a": {"pool": 1000, "found": a_found, "gold_total": 5,
                      "ranks": {}, "top3": []},
                "b": {"pool": 30, "found": b_found, "gold_total": 5,
                      "ranks": {}, "top3": []},
            }],
            "provenance": "p", "known_bias": "b",
        }

    # Ranking up, coverage DOWN: must be reported as a trade, not a win.
    trade = ev.render(spec(0.357, 0.533, 5, 4))
    assert "A trade, not a win" in trade, (
        "a ranking gain with a coverage loss was reported as a win")
    assert "FEWER gold papers" in trade

    # Ranking up, coverage up: a clean win is allowed to say so.
    win = ev.render(spec(0.357, 0.533, 4, 5))
    assert "A trade, not a win" not in win
    assert "wins on both axes" in win

    # Ranking down: the negative result must survive.
    loss = ev.render(spec(0.533, 0.357, 5, 4))
    assert "does not beat" in loss


def test_every_candidate_is_matchable_by_the_concept_pattern():
    """A hand-written detector that does not know a word reports 0% for that
    word, and 0% reads in the table as "the control worked".

    That happened: `helicobacter pylori` scored 0% microbe-in-title over
    13,739 records because `CONCEPT_RE` omitted `helicobacter`. The error ran
    in the direction that flattered the measure, which is the direction an
    instrument must never be wrong in.

    Every candidate's own organism word must be reachable by the pattern, or
    its score is a statement about the pattern rather than the literature.
    """
    import mechanism_keyword_precision as mk

    unreachable = []
    for cand in mk.CANDIDATES:
        # The phrase itself must be matchable, since a title containing the
        # phrase is the easiest possible case.
        if not mk.CONCEPT_RE.search(cand["phrase"]):
            unreachable.append(cand["phrase"])
    assert not unreachable, (
        f"CONCEPT_RE cannot match these candidates' own words: {unreachable}. "
        "Each would score 0% for the detector's reasons, not the corpus's.")


def test_the_concept_pattern_still_rejects_something():
    """Widening a pattern until nothing fails makes it useless. If it matches
    arbitrary non-microbial text, the `microbe in title` column is just a
    count of titles."""
    import mechanism_keyword_precision as mk
    for text in ("Radiotherapy dose escalation in glioma",
                 "A randomised trial of pembrolizumab",
                 "Ferroptosis and lipid peroxidation"):
        assert not mk.CONCEPT_RE.search(text), (
            f"CONCEPT_RE matched {text!r}, which names no microbe")
