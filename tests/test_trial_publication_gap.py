"""Guards for the trial-publication measurement.

This is the one analysis in the project that claims something about studies
that do NOT appear in the literature, so the ways it can be wrong are
asymmetric: over-counting publication understates the gap (safe), while
under-counting it manufactures a finding (not). Every guard here is about
keeping the error on the safe side and keeping the caveats attached to the
number rather than in someone's memory.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import trial_publication_gap as g  # noqa: E402


def _t(nct, *, status="COMPLETED", stype="INTERVENTIONAL", year=2015, pmids=()):
    return {"nct_id": nct, "status": status, "study_type": stype,
            "completion": f"{year}-06-01", "pmids": list(pmids), "phases": []}


def test_only_completed_interventional_trials_are_counted():
    """A terminated or withdrawn trial may have nothing to report, and counting
    it as unpublished charges the literature for a study that never produced a
    result. An observational study has different reporting norms, so mixing
    them compares two populations under one label."""
    assert g.eligible(_t("A"))
    assert not g.eligible(_t("B", status="TERMINATED"))
    assert not g.eligible(_t("C", status="WITHDRAWN"))
    assert not g.eligible(_t("D", status="RECRUITING"))
    assert not g.eligible(_t("E", stype="OBSERVATIONAL"))


def test_a_trial_with_no_completion_date_is_excluded():
    """Time is load-bearing: without a completion year a trial cannot be placed
    in an age band, and pooling it reports lag as gap."""
    t = _t("F")
    t["completion"] = None
    assert not g.eligible(t)
    t["completion"] = "not-a-date"
    assert not g.eligible(t)


def test_the_two_measures_are_reported_separately():
    """THE ERROR THIS PAGE EXISTS TO AVOID. Reading the registry's own
    references field as a publication rate measures REGISTRY UPKEEP: a sponsor
    who never came back leaves a published trial looking unpublished. The
    independent literature check must never be collapsed into it."""
    # The two totals must DIFFER, or swapping one measure for the other
    # changes nothing and the test passes on a collapsed implementation --
    # which is exactly what happened to the first version of this fixture.
    rows = [{"nct": "A", "year": 2015, "phases": [], "enrollment": 1,
             "registry_pmids": 0, "epmc_mentions": 3},
            {"nct": "B", "year": 2015, "phases": [], "enrollment": 1,
             "registry_pmids": 0, "epmc_mentions": 5},
            {"nct": "C", "year": 2015, "phases": [], "enrollment": 1,
             "registry_pmids": 2, "epmc_mentions": 0},
            {"nct": "D", "year": 2015, "phases": [], "enrollment": 1,
             "registry_pmids": 0, "epmc_mentions": 0}]
    s = g.summarise({"rows": rows})
    assert s["registry_says_published"] == 1
    assert s["literature_mentions_it"] == 2, (
        "the literature check is reading the registry field; the two measures "
        "have been collapsed into one")
    assert s["literature_only"] == 2, "trials found only in the literature"
    assert s["registry_only"] == 1, "a trial found only in the registry"
    assert s["neither"] == 1


def test_the_page_states_that_a_mention_over_counts_publication():
    """A paper citing a trial number is not reporting that trial, so the
    measure over-counts publication and any gap is a LOWER bound. If that
    sentence goes, the number reads as stronger than it is."""
    # Enough 2015-2019 rows for the reliable era to exist: the LOWER-bound
    # paragraph is attached to that headline, so a fixture without it renders
    # a page missing the caveat and the test would pass for the wrong reason.
    rows = [{"nct": f"A{i}", "year": 2017, "phases": [], "enrollment": 1,
             "registry_pmids": 0, "epmc_mentions": 0} for i in range(12)]
    md = g.render(g.assemble({"eligible": 10, "sampled": len(rows), "seed": 1,
                                    "rows": rows}))
    assert "LOWER bound" in md
    assert "not necessarily REPORTING" in md
    assert "Registration is not running" in md
    assert "A gap is not misconduct" in md


def test_the_headline_comes_from_the_era_where_the_measure_works():
    """THE CORRECTION THIS PAGE IS BUILT AROUND.

    Citing a registration number in a paper is a modern convention -- it became
    an ICMJE condition of publication in 2005 and spread gradually after. The
    measured mention rate climbs 32% (2000-04) to 83% (2015-19), and an earlier
    draft read that as trials becoming more likely to be published. It is the
    opposite kind of fact: old trials are invisible to this measure whether or
    not they were published, so their low rate measures the CONVENTION.

    The headline must therefore come from the era where the convention holds,
    not from the whole sample.
    """
    rows = ([{"nct": f"OLD{i}", "year": 2002, "phases": [], "enrollment": 1,
              "registry_pmids": 0, "epmc_mentions": 0} for i in range(40)]
            + [{"nct": f"MID{i}", "year": 2017, "phases": [], "enrollment": 1,
                "registry_pmids": 0, "epmc_mentions": 1} for i in range(40)])
    s = g.summarise({"rows": rows})
    assert s["reliable_era"] == "2015-2019"
    assert s["reliable_era_stats"]["pct"] == 100.0
    assert s["reliable_era_gap_pct"] == 0.0, (
        "the headline is being computed over the whole sample, so the old "
        "trials' citation-convention artifact is contaminating it")
    # the pooled figure is far worse, which is exactly the trap
    assert s["literature_pct"] == 50.0


def test_the_page_names_the_artifact_rather_than_reporting_it_as_a_finding():
    """A reader who sees the era table without the explanation concludes that
    old trials went unpublished. The explanation is the finding."""
    rows = ([{"nct": f"O{i}", "year": 2002, "phases": [], "enrollment": 1,
              "registry_pmids": 0, "epmc_mentions": 0} for i in range(12)]
            + [{"nct": f"M{i}", "year": 2017, "phases": [], "enrollment": 1,
                "registry_pmids": 1, "epmc_mentions": 1} for i in range(12)])
    md = g.render(g.assemble({"eligible": 100, "sampled": 24, "seed": 1,
                                    "rows": rows}))
    assert "measurement artifact" in md
    assert "PRINT THE REGISTRATION NUMBER" in md
    assert "withdrawn" in md, "the earlier reading is not retracted on the page"


def test_the_sample_is_seeded_so_the_figure_reproduces():
    """An unseeded sample gives a different number every run, and a claim about
    missing evidence that moves when you re-run it is not a measurement."""
    import inspect
    src = inspect.getsource(g.measure)
    assert "random.Random(SEED)" in src
    assert isinstance(g.SEED, int)


def test_every_number_on_the_page_is_derived():
    """A hand-written figure beside a derived one is this repository's most
    persistent defect."""
    rows = [{"nct": "A", "year": 2015, "phases": [], "enrollment": 1,
             "registry_pmids": 1, "epmc_mentions": 1}] * 7
    md = g.render(g.assemble({"eligible": 4321, "sampled": 7, "seed": 99,
                                    "rows": rows}))
    assert "4,321" in md and "seed 99" in md
    assert "600" not in md, "a default sample size is baked into the prose"


# --- Order-independence, reinstated offline ---------------------------------
#
# `trial_publication_gap` is EXEMPT in tests/test_artifact_freshness.py: its
# input is a live six-minute API sample over shards CI does not have.
# Exemption removes it from LIVE, which is what that file's order gate
# iterates -- so the gate does NOT run on it, and an exemption reason claiming
# otherwise would be false. I have now written that false claim three times
# (corpus_expand_report, fulltext_ceiling, and this), which is why the guard is
# reinstated here rather than asserted there.

def _shuffled(obj, rng):
    if isinstance(obj, dict):
        items = list(obj.items())
        rng.shuffle(items)
        return {k: _shuffled(v, rng) for k, v in items}
    if isinstance(obj, list):
        return [_shuffled(v, rng) for v in obj]
    return obj


def _payload():
    rows = []
    for i in range(30):
        rows.append({"nct": f"NCT{i:08d}", "year": 2000 + (i % 25),
                     "phases": ["PHASE2"], "enrollment": 10 + i,
                     "registry_pmids": i % 3, "epmc_mentions": (i * 7) % 5})
    return g.assemble({"eligible": 5000, "sampled": len(rows), "seed": 1,
                       "rows": rows})


def test_the_page_order_is_established_by_the_renderer():
    """Shuffling every dict must not change a character of the output."""
    import random
    base = g.render(_payload())
    rng = random.Random(4242)
    for _ in range(20):
        assert g.render(_shuffled(_payload(), rng)) == base


def test_the_era_table_order_survives_a_reordered_list():
    """`_shuffled` preserves list order, so reversing `rows` is the one
    transformation that exposes an era table inherited from row order."""
    p = _payload()
    base = g.render(p)
    rev = _payload()
    rev["rows"] = list(reversed(rev["rows"]))
    rev["summary"] = g.summarise(rev)
    assert g.render(rev) == base
