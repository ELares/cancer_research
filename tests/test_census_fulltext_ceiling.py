"""Guards for the full-text ceiling measurement.

This artifact exists to STOP work, so its guards protect against the failure
mode of an optimistic bound: a ceiling that reads higher than it is, or two
skews that quietly go missing, would license exactly the project the
measurement was run to bound.

The bound is deliberately the most OPTIMISTIC version of the answer -- a PMC
identifier means a record is in PMC, not that its text is in the open-access
subset held here, nor that it parsed. If the real ceiling can only be lower,
then a decision not to build survives the uncertainty and a decision to build
would not.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-fulltext-ceiling.json"
MD = REPO / "analysis/census-fulltext-ceiling.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_ceiling_recomputes_from_its_own_counts(d):
    u = next(r for r in d["rows"] if r["class"] == "undetermined")
    assert d["ceiling_records"] == u["reachable"]
    assert d["unreachable_records"] == u["total"] - u["reachable"]
    assert d["ceiling_share_of_undetermined"] == pytest.approx(
        100 * u["reachable"] / u["total"], abs=0.05)
    assert d["ceiling_share_of_census"] == pytest.approx(
        100 * u["reachable"] / d["census"], abs=0.05)
    # The unreachable majority is the finding; if it ever becomes a minority
    # the conclusion changes and the prose must be rewritten, not kept.
    assert d["unreachable_records"] > d["ceiling_records"], (
        "most undetermined records are now reachable, which inverts this "
        "analysis's conclusion -- rewrite it rather than leaving the prose")


def test_both_skews_are_derived_and_present(d):
    ds, es = d["design_skew"], d["era_skew"]
    rated = [r for r in d["rows"] if r["class"] != "undetermined" and r["rate"]]
    assert ds["highest"]["rate"] == max(r["rate"] for r in rated)
    assert ds["lowest"]["rate"] == min(r["rate"] for r in rated)
    assert ds["fold"] == pytest.approx(ds["highest"]["rate"] / ds["lowest"]["rate"], abs=0.05)
    assert es["fold"] == pytest.approx(es["since"]["rate"] / es["before"]["rate"], abs=0.05)
    md = MD.read_text()
    assert f"a factor of {ds['fold']}" in md
    assert f"a factor of {es['fold']}" in md


def test_the_era_shift_is_stated_as_a_year_difference(d):
    """A fold-change in a rate is abstract; "the recovered set is 11 years
    younger than the pile" is the same fact in a form a reader can act on."""
    assert d["median_year_reachable"] > d["median_year_pile"]
    md = MD.read_text()
    shift = d["median_year_reachable"] - d["median_year_pile"]
    assert f"{shift} years later" in md


def test_it_refuses_to_license_a_correction(d):
    """The whole point. A distribution recovered from the reachable fifth
    describes the readable literature, and publishing it as a correction would
    swap one population for another without saying so."""
    md = MD.read_text()
    assert "never merged into an NLM-labelled one" in md
    assert "most optimistic version of the answer" in md
    for overclaim in ("would fix", "closes the gap", "resolves the undetermined"):
        assert overclaim not in md.lower()


def test_the_scan_spreads_across_eras_not_a_prefix(scanner, tmp_path, monkeypatch):
    """Shards are chronological. A prefix would sample one era and destroy the
    era measurement this script exists for -- a mistake this repo has already
    made once on a different scan."""
    import gzip

    monkeypatch.setattr(scanner, "RECORDS", tmp_path)
    for name, years in (("c", [2020, 2021]), ("a", [1980, 1981]), ("b", [2000])):
        with gzip.open(tmp_path / f"{name}.jsonl.gz", "wt") as fh:
            for year in years:
                fh.write(json.dumps({"year": year}) + "\n")
    result = scanner.scan(stride=2)
    assert result["census"] == 4
    assert set(result["undetermined_by_year"]) == {"1980", "1981", "2020", "2021"}


def test_the_manuscript_states_the_ceiling_where_it_states_the_gap(d):
    """The 44.5% undetermined share and its bound belong together. Reporting
    the gap without the bound invites a reader to assume it is closable."""
    txt = " ".join(MANUSCRIPT.read_text().split())
    if "44.5%" not in txt:
        pytest.skip("the manuscript no longer reports the undetermined share")
    assert f"{d['ceiling_share_of_undetermined']}%" in txt, (
        "the manuscript reports the design-label gap without the measured "
        "ceiling on closing it")
    assert f"{d['unreachable_records']:,}" in txt or "four-fifths" in txt


@pytest.fixture
def scanner(monkeypatch):
    import importlib.util

    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        "fulltext_ceiling_sparse_test", REPO / "scripts/census_fulltext_ceiling.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("totals,reachable,years", [
    ({"undetermined": 1}, {}, {}),
    ({"trial": 1}, {}, {}),
    ({"trial": 1, "animal": 1}, {}, {}),
    ({"trial": 1, "animal": 1}, {"trial": 1}, {}),
    ({"undetermined": 1}, {"undetermined": 1}, {"1990": [1, 1]}),
    ({"undetermined": 1}, {}, {"2020": [1, 0]}),
    ({"undetermined": 2}, {"undetermined": 1}, {"1990": [1, 0], "2020": [1, 1]}),
    ({"undetermined": 2}, {"undetermined": 1}, {"1990": [1, 1], "2020": [1, 0]}),
], ids=["undated-unclassified", "no-undetermined", "zero-reachability",
        "one-zero-design-rate", "only-earlier-era", "only-later-era",
        "earlier-rate-zero", "later-rate-zero"])
def test_sparse_populations_have_observed_counts_and_unavailable_comparisons(
        scanner, totals, reachable, years):
    result = scanner.assemble({
        "census": sum(totals.values()), "by_class_total": totals,
        "by_class_reachable": reachable, "undetermined_by_year": years,
    })
    markdown = scanner.render(result)

    assert result["ceiling_records"] == reachable.get("undetermined", 0)
    assert "None" not in markdown
    assert "a fifth" not in markdown and "four-fifths" not in markdown
    assert "this one is larger" not in markdown
    if not totals.get("undetermined"):
        assert result["ceiling_share_of_undetermined"] is None
        assert result["ceiling_share_of_census"] == 0.0
        assert "No undetermined records" in markdown
    if len([c for c in totals if c != "undetermined"]) < 2:
        assert result["design_skew"]["highest"] is None
        assert "fewer than two" in markdown
    else:
        assert result["design_skew"]["lowest"]["rate"] == 0.0
        assert result["design_skew"]["fold"] is None
        assert "lowest rate is zero" in markdown
    if not years:
        assert result["era_skew"]["before_share_of_undetermined"] is None
        assert "median-year shift is unavailable" in markdown
    if years == {"1990": [1, 1], "2020": [1, 0]}:
        assert result["era_skew"]["fold"] == 0.0
        assert "15 years earlier" in markdown


def test_reassembly_preserves_published_fulltext_counts(scanner, d):
    assert scanner.assemble(d) == d
