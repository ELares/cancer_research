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


# --- The headline is a proportion over one stratum, not over the sample ---

def _rows(n_mentioned, n_total, year=2017):
    """A reliable-era stratum with a known mention rate."""
    out = []
    for i in range(n_total):
        out.append({"nct": f"NCT{i:08}", "year": year, "phases": [],
                    "enrollment": 10, "registry_pmids": 0,
                    "epmc_mentions": 1 if i < n_mentioned else 0})
    return {"eligible": n_total, "sampled": n_total, "seed": 1, "rows": out}


def test_the_headline_gap_carries_its_interval():
    """16.8% to a tenth of a percent, from 149 trials, is precision the
    stratum cannot support.

    The gap is measured over ONE completion era -- 149 of the 600 sampled --
    so the sample size that matters is the stratum's. `config.py` already
    states a precision with its n and its CI; a claim about MISSING MEDICAL
    EVIDENCE should not be held to a looser standard than a keyword audit.
    """
    d = g.assemble(_rows(124, 149))
    s = d["summary"]
    lo, hi = s["reliable_era_gap_ci"]
    assert lo < s["reliable_era_gap_pct"] < hi, (
        f"the interval {lo}-{hi} does not bracket the point estimate "
        f"{s['reliable_era_gap_pct']}")
    md = g.render(d)
    assert "95% CI" in md, "the headline states a gap with no interval"
    assert f"n={s['reliable_era_stats']['trials']}" in md, (
        "the interval is printed without the n it came from, which is the "
        "half that makes it readable")


def test_a_smaller_stratum_widens_the_interval():
    """An interval that ignores n is decoration.

    The same observed rate over 28 trials and over 1,400 means very different
    things, and the early eras in this table really are that small.
    """
    narrow = g.assemble(_rows(1400, 1680))["summary"]["reliable_era_gap_ci"]
    wide = g.assemble(_rows(23, 28))["summary"]["reliable_era_gap_ci"]
    assert (wide[1] - wide[0]) > (narrow[1] - narrow[0]) * 3, (
        f"interval over 28 trials is {wide}, over 1,680 is {narrow}; it is "
        "not a function of the sample size")


def test_an_absent_trial_share_is_refused_not_rendered():
    """The recurring mechanical failure, and it used to surface as arithmetic.

    The trial shards live on an external share. Unmounted, `load_trials()`
    returns [] and every share is None, and the first version then raised
    `TypeError: unsupported operand type(s) for -: 'int' and 'NoneType'` from
    inside `render` -- which says nothing about a mount. The freshness
    exemption even described this case as rendering "an empty page", a
    behaviour nobody had checked.
    """
    import pytest

    with pytest.raises(g.NoTrialsError) as e:
        g.assemble({"eligible": 0, "sampled": 0, "seed": 1, "rows": []})
    msg = str(e.value)
    assert "mounted" in msg or "mount" in msg, (
        f"the refusal does not name the likely cause: {msg}")
    assert "FERRO_TRIALS_OUT" in msg, (
        "the refusal does not name the override that fixes it")


def test_the_registry_shortfall_share_is_derived():
    """It read "a third of the sample" in hand-written prose, beside a count.

    34.3% happens to be about a third, which is exactly why this kind of
    sentence survives: it is right until the sample changes and then it is
    quietly wrong.
    """
    d = g.assemble(_rows(60, 100))
    d["rows"][0]["registry_pmids"] = 0
    d["summary"] = g.summarise(d)
    md = g.render(d)
    share = round(100 * d["summary"]["literature_only"] / d["summary"]["sampled"], 1)
    assert f"{share}% of the sample" in md, (
        f"the page does not print the derived share {share}%")
    assert "a third of the sample" not in md, (
        "the hand-written fraction is back")


def test_render_only_rebuilds_the_page_from_stored_rows(tmp_path, monkeypatch):
    """Fixing a sentence must not cost 600 live requests against a public API.

    That cost is a standing reason not to fix sentences. `summarise` and
    `render` are pure functions of the stored rows, so the page rebuilds
    offline -- and the figures cannot drift while doing it, which a re-run
    against a live index could not promise.
    """
    import json
    import subprocess

    d = _rows(124, 149)
    out_json = tmp_path / "gap.json"
    out_md = tmp_path / "gap.md"
    out_json.write_text(json.dumps(g.assemble(d), indent=2, sort_keys=True) + "\n")

    monkeypatch.setattr(g, "OUT_JSON", out_json)
    monkeypatch.setattr(g, "OUT_MD", out_md)

    def no_network(*a, **k):
        raise AssertionError("--render-only reached the network")

    monkeypatch.setattr(g, "measure", no_network)
    monkeypatch.setattr(sys, "argv", ["trial_publication_gap.py", "--render-only"])
    assert g.main() == 0
    assert "95% CI" in out_md.read_text()
