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

    That generator reads a 20 GB gitignored index built from shards on external
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
            "gold_in_index": len(qs), "gold_missing_from_index": [],
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
            "gold_in_index": 1, "gold_missing_from_index": [],
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
            "sampling": "uniform random draw from the whole pool",
            "candidates": [
                {"phrase": p, "control": False, "pool": 5, "inspected": 5,
                 "microbe_in_title": 0.8, "both_in_title": 0.6,
                 "verdict": "admit", "example_titles": []}
                for p in order
            ],
            "sample": [], "admitted": list(order), "refused": [],
            "control_separates": True, "worst_real": 0.6, "best_control": 0.1,
            "separation_margin": 0.5, "separation_se": 0.03,
            "separation_exceeds_noise": True,
            "separation_worst_real": "alpha", "separation_worst_real_n": 200,
            "separation_best_control": "ctrl", "separation_best_control_n": 200,
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
            "sampling": "uniform random draw from the whole pool",
            "candidates": [], "sample": [], "admitted": ["x"], "refused": [],
            "control_separates": separates,
            "worst_real": 0.1 if not separates else 0.6,
            "best_control": 0.5 if not separates else 0.1,
            "separation_margin": -0.4 if not separates else 0.5,
            "separation_se": 0.03,
            "separation_exceeds_noise": bool(separates),
            "separation_worst_real": "r", "separation_worst_real_n": 200,
            "separation_best_control": "c", "separation_best_control_n": 200,
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
            "at_k": 10, "gold_in_index": 5, "gold_missing_from_index": [],
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


def test_the_keyword_measure_does_not_inspect_a_title_ranked_head(tmp_path):
    """The circularity that inflated every row of the first keyword table.

    `aboutness` asks what share of matches name the concept in the TITLE. The
    first version handed it `search(c, m, limit=DEPTH)`, whose default weights
    are W_TITLE_HEAVY, so bm25 had already selected the documents matching in
    the title -- and a candidate phrase contains its own concept word, so the
    ranking was partly producing the number. Measured on the real index the
    inflation reached 43% against 24% for `helicobacter pylori`, and it scaled
    with pool size, so it hit the CONTROLS hardest and pushed the verdict
    toward refusing.

    The fixture separates the two populations: twenty records carry the phrase
    in the title, twenty carry it only in the abstract. A title-weighted head
    of ten is therefore all title-matches and scores 100%; a uniform draw must
    reach the abstract-only half and score lower.
    """
    import random

    import mechanism_keyword_precision as mk

    c = _db(tmp_path)
    rows = []
    for i in range(20):
        rows.append(("Engineered bacteria attack the tumor %d" % i,
                     "An unrelated abstract body.", "", f"T{i:03}"))
    for i in range(20):
        rows.append(("A study of delivery vectors %d" % i,
                     "Here we use engineered bacteria against a tumor.", "",
                     f"A{i:03}"))
    for title, abstract, mesh, pmid in rows:
        c.execute(
            "INSERT INTO doc (title, abstract, mesh, pmid, pmcid, doi, journal,"
            " year, has_text, cited_by) VALUES (?,?,?,?,'','','','',0,-1)",
            (title, abstract, mesh, pmid))
    c.commit()

    m = rs.phrase_query(["engineered bacteria"])
    assert rs.count(c, m) == 40, "fixture must put every record in the pool"

    ranked = rs.search(c, m, rs.W_TITLE_HEAVY, limit=10)
    assert all(r["pmid"].startswith("T") for r in ranked), (
        "fixture is not exercising the defect: the title-weighted head should "
        "be entirely title-matches")
    assert mk.aboutness(ranked)["both_in_title"] == 1.0

    drawn = mk.inspect_pool(c, m, 10, random.Random("fixed"))
    assert len(drawn) == 10
    assert any(r["pmid"].startswith("A") for r in drawn), (
        "a uniform draw over a 20/20 pool reached no abstract-only record, so "
        "it is still sampling the ranked head")
    assert mk.aboutness(drawn)["both_in_title"] < 1.0, (
        "the uniform draw scores the same as the title-ranked head, which is "
        "the circularity this test exists to prevent")
    c.close()


def test_the_keyword_sample_is_drawn_from_the_same_uniform_set(tmp_path):
    """The judging sample must not be the most title-like matches either.

    A human reading a sample to decide precision is deciding it for the
    KEYWORD, over its pool. Handing them the bm25 head would hand them the
    easiest records and call the result the keyword's precision -- the same
    defect as the measure, one step further downstream where it is harder to
    see.
    """
    import random

    import mechanism_keyword_precision as mk

    c = _db(tmp_path)
    for i in range(30):
        where = ("Engineered bacteria in tumor %d" % i, "Body text.") if i < 15 \
            else ("Vector study %d" % i, "We use engineered bacteria in tumor work.")
        c.execute(
            "INSERT INTO doc (title, abstract, mesh, pmid, pmcid, doi, journal,"
            " year, has_text, cited_by) VALUES (?,?,'',?,'','','','',0,-1)",
            (where[0], where[1], f"P{i:03}"))
    c.commit()
    m = rs.phrase_query(["engineered bacteria"])
    drawn = mk.inspect_pool(c, m, 12, random.Random("fixed"))
    titled = sum(1 for r in drawn if r["title"].lower().startswith("engineered"))
    assert 0 < titled < 12, (
        f"a uniform draw of 12 from a 15/15 pool gave {titled} title-matches; "
        "an all-or-nothing split means the draw is following the ranking")
    c.close()


def test_the_evaluation_refuses_to_call_a_corpus_gap_a_ranking_result():
    """0.000 on both arms is what the WRONG CORPUS looks like.

    This already happened: indexing the expansion crawl alone put none of the
    gold in the index, both arms scored zero, and the page rendered a verdict
    about RANKING. The crawl skips anything `is_held()` knows, so it is the
    census's complement and cannot hold a paper the census had. The report has
    to distinguish "ranked badly" from "never indexed", and it can only do that
    if membership is measured rather than taken from the gold file's prose.
    """
    import corpus_retrieval_eval as ev

    def spec(missing):
        return {
            "docs_indexed": 7_727_217, "queries": 1, "gold_total": 2,
            "depth": 50, "at_k": 10,
            "gold_in_index": 2 - len(missing),
            "gold_missing_from_index": missing,
            "aggregate": {
                "a": {"mrr": 0.0, "recall_at_k": 0.0, "gold_found": 0, "median_pool": 9},
                "b": {"mrr": 0.0, "recall_at_k": 0.0, "gold_found": 0, "median_pool": 9}},
            "per_query": [{
                "id": "q", "question": "q?", "concepts": ["x"], "note": "",
                "gold": ["111", "222"], "gold_absent_from_index": missing,
                "a": {"pool": 9, "found": 0, "gold_total": 2, "ranks": {}, "top3": []},
                "b": {"pool": 9, "found": 0, "gold_total": 2, "ranks": {}, "top3": []}}],
            "provenance": "prov.", "known_bias": "bias."}

    red = ev.render(spec(["111", "222"]))
    assert "not a retrieval result" in red.lower(), (
        "an index holding none of the gold rendered as an ordinary result")
    assert red.index("not a retrieval result") < red.index("| MRR |"), (
        "the warning must come BEFORE the score table, or it is a footnote on "
        "a number the reader has already believed")
    assert "complement" in red, (
        "the page should name the known cause, which is the crawl being the "
        "census's complement by construction")

    ok = ev.render(spec([]))
    assert "not a retrieval result" not in ok.lower(), (
        "a complete index must not carry the corpus-gap warning")
    assert "present in the index" in ok


def test_the_keyword_report_will_not_call_a_gap_inside_the_noise_a_separation():
    """`control_separates` is a bare `>` and would call a one-record gap a win.

    This matters because it actually happened in the other direction: removing
    the title-ranking circularity moved the verdict from "does not
    discriminate" to a `>` that holds by seven points between a 27% estimate
    over 168 records and a 20% estimate over 200. A page that announces a
    calibrated measure on that is how an instrument gets trusted for the wrong
    reason, so the headline has to read the margin against its own standard
    error and not just the inequality.
    """
    import mechanism_keyword_precision as mk

    def spec(worst, best, n, exceeds):
        return {
            "docs_indexed": 100, "seed": 1, "depth": 200, "aboutness_floor": 0.3,
            "sampling": "uniform random draw from the whole pool",
            "candidates": [], "sample": [], "admitted": [], "refused": [],
            "control_separates": True,
            "worst_real": worst, "best_control": best,
            "separation_margin": worst - best,
            "separation_se": 0.04,
            "separation_exceeds_noise": exceeds,
            "separation_worst_real": "real", "separation_worst_real_n": n,
            "separation_best_control": "ctrl", "separation_best_control_n": n,
        }

    thin = mk.render(spec(0.27, 0.20, 200, False))
    assert "Not decidably" in thin, (
        "a gap inside its own standard error was announced as discrimination")
    assert "**Yes.**" not in thin

    wide = mk.render(spec(0.90, 0.10, 200, True))
    assert "**Yes.**" in wide, (
        "a gap of 80 points over n=200 must still be reportable as a "
        "separation, or this guard just suppresses every positive result")
    assert "Not decidably" not in wide


def test_the_separation_standard_error_falls_as_the_sample_grows():
    """A noise bound that ignores n is decoration.

    The whole point of computing an SE is that the same observed gap means
    different things over 20 records and over 2,000. If `_separation` returned
    a constant, the guard above would pass while measuring nothing.
    """
    import mechanism_keyword_precision as mk

    def rows(n):
        real = [{"phrase": "r", "both_in_title": 0.30, "inspected": n}]
        ctrl = [{"phrase": "c", "both_in_title": 0.20, "inspected": n}]
        return mk._separation(real, ctrl)

    small, large = rows(20), rows(2000)
    assert small["separation_margin"] == pytest.approx(large["separation_margin"]), \
        "the fixture must hold the observed gap fixed and vary only n"
    assert small["separation_se"] > large["separation_se"], (
        "the standard error did not fall as the sample grew, so it is not a "
        "function of n")
    assert not small["separation_exceeds_noise"], (
        "a 10-point gap over 20 records is well inside the noise")
    assert large["separation_exceeds_noise"], (
        "the same 10-point gap over 2,000 records is real and must be "
        "reportable, or the bound rejects everything")


def test_a_textless_stub_does_not_burn_the_key_for_the_real_paper(tmp_path):
    """A record that is DROPPED must not be recorded as indexed.

    `seen` means "this paper is in the index". The first version added the key
    the moment it read it -- one line above the check that drops a record with
    neither title nor abstract -- so an ID-only stub claimed the key, was
    discarded, and the same paper arriving later WITH text was skipped as a
    duplicate. It could not heal either: the stub's key never reached the
    `seen` table, so a resumed run would accept the record, but the shard is
    marked scanned and never read again.

    Measured on the real corpus: 21 of 17,584 sampled census records carry an
    identity key and no indexable text (0.12% of that slice), so the stubs
    that claim a key exist. This fixture makes the loss they enable explicit.
    """
    stub = {"pmid": "40000001"}                      # keyed, no title, no abstract
    real = {"pmid": "40000001",
            "title": "Engineered bacteria for pancreatic cancer",
            "abstract": "The paper that answers the question."}
    _shard(tmp_path / "a-stub.jsonl.gz", [stub])
    _shard(tmp_path / "b-real.jsonl.gz", [real])

    c = _db(tmp_path)
    seen: set[str] = set()
    # Scan order matters: the stub is read first, exactly as sorted shard
    # order would deliver it.
    ri.index_shard(c, tmp_path / "a-stub.jsonl.gz", seen, verbose=False)
    ri.index_shard(c, tmp_path / "b-real.jsonl.gz", seen, verbose=False)

    titles = [r[0] for r in c.execute("SELECT title FROM doc")]
    assert titles == ["Engineered bacteria for pancreatic cancer"], (
        f"the real paper was lost behind a textless stub sharing its PMID: "
        f"indexed titles were {titles}")
    assert c.execute("SELECT COUNT(*) FROM seen").fetchone()[0] == 1, (
        "the key should be claimed exactly once, by the record that was "
        "actually indexed")
    c.close()


def test_the_stub_itself_is_still_not_indexed(tmp_path):
    """The complement, so the fix above cannot be 'index everything'.

    Relaxing the dedup claim must not relax the text requirement: a record
    with no title and no abstract has nothing for FTS5 to match and indexing
    it could only add noise to every pool.
    """
    _shard(tmp_path / "only-stubs.jsonl.gz",
           [{"pmid": "40000002"}, {"pmid": "40000003", "title": "", "abstract": ""}])
    c = _db(tmp_path)
    added, bad = ri.index_shard(c, tmp_path / "only-stubs.jsonl.gz", set(), verbose=False)
    assert added == 0, "a textless record was indexed"
    assert bad == 0, "a textless record is not a PARSE failure and must not be counted as one"
    assert c.execute("SELECT COUNT(*) FROM seen").fetchone()[0] == 0, (
        "a record that was never indexed must not claim a dedup key")
    c.close()
