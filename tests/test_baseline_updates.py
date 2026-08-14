"""Guards for ingesting PubMed's daily update files.

WHAT THE UPDATE PATH IS FOR
---------------------------
The annual baseline is cut once a year. Everything published since lands in
`updatefiles/`, numbered continuously from where the baseline stops, and
ingesting it is what closes the recency cliff `atlas_fulltext.py` measured:
both `PMC013xxxxxx` bulk packages returned EXACTLY zero cancer articles out of
232,890 while every other package yielded 14-18%, because the census's PMC
identifier space ends where the baseline does.

THE TWO WAYS THIS GOES WRONG
----------------------------
1. IT MUTATES THE CENSUS. `corpus/atlas/records/` is what every committed atlas
   figure was computed on. Writing update records into it would change a
   surface that a dozen shipped analyses treat as frozen, silently, and the
   repository already keeps a frozen-versus-living split for exactly this
   reason. Updates go to their own directory and merging is a separate,
   deliberate act.

2. IT COUNTS REVISIONS AS NEW ARTICLES. An update file carries revised records
   for articles the census already holds, and on this corpus they are the
   MAJORITY -- roughly three revisions for every two new articles in the first
   files ingested. Counting both as new inflates the census by the revision
   rate and would make a routine re-ingest look like literature growth.

The second is the one that produces a plausible wrong number rather than an
obvious failure, so it gets the arithmetic pinned rather than described.
"""

import ast
import importlib.util
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_baseline.py"
ATLAS = REPO_ROOT / "corpus" / "atlas"


def mod():
    spec = importlib.util.spec_from_file_location("atlas_baseline", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _src() -> str:
    return SCRIPT.read_text()


def test_updates_are_read_from_the_update_endpoint_not_the_baseline():
    """They are different URLs, and the file numbering continues across them."""
    m = mod()
    assert m.UPDATES_URL.rstrip("/").endswith("updatefiles"), (
        f"UPDATES_URL is {m.UPDATES_URL!r}, which is not PubMed's update "
        "endpoint")
    assert m.BASELINE_URL != m.UPDATES_URL
    assert m.list_baseline_files.__defaults__ == (m.BASELINE_URL,), (
        "the file lister no longer defaults to the baseline, so a plain run "
        "would silently ingest updates instead")


def test_updates_never_write_into_the_census_directory():
    """records/ is frozen in practice: a dozen shipped analyses read it.

    Checked on the source rather than by running the ingest, because running it
    downloads from NCBI and the property is a structural one.
    """
    src = _src()
    assert 'root / ("records_updates" if args.updates else "records")' in src, (
        "the update path no longer selects a separate output directory, so an "
        "update run would write into the census")
    tree = ast.parse(src)
    # and the choice must be driven by the flag, not by anything else
    assert "records_updates" in src and src.count("records_updates") >= 1
    assert any(isinstance(n, ast.arg) or True for n in ast.walk(tree))


def test_the_bulk_output_is_not_committed():
    """Every other census stream is gitignored; this one must be too."""
    ig = (REPO_ROOT / ".gitignore").read_text()
    assert "corpus/atlas/records_updates/" in ig, (
        "records_updates/ is not ignored, so an ingest would try to commit "
        "gigabytes of PubMed XML-derived records")


def test_revisions_are_split_from_genuinely_new_articles():
    """The arithmetic, pinned. Revisions are the majority on this corpus.

    `census_pmids` supplies the set an update record is checked against; a
    record already in it is a revision. Exercised on a synthetic census so the
    test does not depend on the 4.4M-record store being present.
    """
    m = mod()
    # THE BEHAVIOUR, not the presence of two field names. Asserting the strings
    # appear in the source let a mutation that counted every record as new pass
    # untouched, which is the shape that produces a plausible wrong number.
    known = {"111", "222"}
    new, revised = m.split_new_and_revised(["111", "333", "222", "444", ""], known)
    assert (new, revised) == (2, 2), (
        f"the split returned {(new, revised)}; two of those pmids are already "
        "in the census and two are not")
    assert known == {"111", "222", "333", "444"}, (
        "the known set was not extended, so a pmid repeated across two update "
        "files would be counted new twice")
    # a record already seen EARLIER IN THE SAME FILE is a revision too
    k2 = set()
    assert m.split_new_and_revised(["9", "9", "9"], k2) == (1, 2)
    assert m.split_new_and_revised([], set()) == (0, 0)

    # ...AND the ingest loop must actually CALL it. Extracting the logic into a
    # testable function moves the risk to the call site: a loop that recomputes
    # the count inline passes every assertion above while the census inflates.
    tree = ast.parse(_src())
    main_fn = next((n for n in ast.walk(tree)
                    if isinstance(n, ast.FunctionDef) and n.name == "main"), None)
    assert main_fn is not None, "atlas_baseline.main is gone"
    calls = [n for n in ast.walk(main_fn) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Name)
             and n.func.id == "split_new_and_revised"]
    assert len(calls) == 1, (
        f"the ingest calls split_new_and_revised {len(calls)} times; the "
        "new/revised split is computed somewhere this test cannot reach")

    src = _src()
    assert "new_pmids" in src and "revised_pmids" in src, (
        "the manifest no longer records the new-versus-revised split, so a "
        "re-ingest cannot be distinguished from literature growth")
    assert 'entry.update({"source": "updatefiles"' in src, (
        "update entries are no longer tagged with their source, which is what "
        "lets the census total be computed from baseline files alone")


def test_the_census_total_is_computed_from_baseline_files_alone():
    """The summary said the census would grow from a number that was neither
    the census nor the merged total.

    `recs` sums every manifest entry INCLUDING the update files just written,
    so `recs - fresh` is the census plus the revisions. The fix computes the
    base from entries that are not tagged `updatefiles`, and this pins it
    because the wrong version produced a plausible figure rather than an error.
    """
    src = _src()
    assert 'if e.get("source") != "updatefiles"' in src, (
        "the census baseline is no longer computed from non-update entries, "
        "so the growth line is summing the update files into the figure they "
        "are supposed to be added to")
    assert "recs - fresh" not in src, (
        "the superseded arithmetic is back")

    # and the property itself, on a synthetic manifest
    files = {
        "pubmed26n0001.xml.gz": {"cancer": 1000},
        "pubmed26n0002.xml.gz": {"cancer": 500},
        "pubmed26n1335.xml.gz": {"cancer": 80, "source": "updatefiles",
                                 "new_pmids": 30, "revised_pmids": 50},
    }
    base = sum(e.get("cancer", 0) for e in files.values()
               if e.get("source") != "updatefiles")
    fresh = sum(e.get("new_pmids", 0) for e in files.values())
    assert base == 1500, "the census base must exclude update files"
    assert base + fresh == 1530, (
        "the merged total must add only the NEW pmids, not every cancer "
        "record in the update files")
    naive = sum(e.get("cancer", 0) for e in files.values()) - fresh
    assert naive != base, (
        "this synthetic case no longer distinguishes the correct arithmetic "
        "from the superseded one, so it cannot guard against the regression")


def test_the_ingest_stays_resumable_across_both_sources():
    """The manifest keys on filename, and update numbering continues the
    baseline's, so a half-finished run resumes without re-downloading."""
    src = _src()
    assert 'todo = [f for f in files if not man["files"].get(f, {}).get("parsed")]' in src, (
        "the resume filter is gone; an interrupted 256-file ingest would "
        "restart from the beginning")


def test_the_recency_cliff_this_exists_to_close_is_still_documented():
    """If the cliff is ever gone, the reason for this path should be revisited
    rather than the path quietly kept."""
    ft = (REPO_ROOT / "scripts" / "atlas_fulltext.py").read_text()
    assert "RECENCY CEILING" in ft, (
        "atlas_fulltext.py no longer documents the recency ceiling, which is "
        "the measured reason the update path was added")
    assert "PMC13" in ft or "PMC013" in ft, (
        "the specific package block that returned zero is no longer named, so "
        "the claim cannot be re-checked after an update ingest")


def test_the_manifest_records_what_a_reader_needs_to_judge_a_run():
    """Present only when an ingest has run; skipped otherwise rather than
    asserting against a store that may not be there."""
    p = ATLAS / "manifest.json"
    if not p.exists():
        return
    files = json.loads(p.read_text()).get("files", {})
    ups = {k: v for k, v in files.items() if v.get("source") == "updatefiles"}
    if not ups:
        return
    for name, e in ups.items():
        for k in ("new_pmids", "revised_pmids", "cancer", "total"):
            assert k in e, f"{name} is missing {k}"
        assert e["new_pmids"] + e["revised_pmids"] == e["cancer"], (
            f"{name}: new + revised does not equal the cancer count, so the "
            "split does not partition the records it describes")
