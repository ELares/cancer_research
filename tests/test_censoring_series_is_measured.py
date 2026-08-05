"""The censoring artifact must be MEASURED, not remembered.

`atlas-replication.md`'s methodology section rests on a comparison: scoring
cohorts on whether they were EVER replicated shows a steep decline, and that
decline is the observation window shrinking rather than science changing. It is
one of the census's headline findings, and `census_findings.py` re-states it.

For most of its life that comparison was two numbers ("about 60% in 1950 to
17.5% in 2020") written into prose from a development run whose output was never
stored. Neither document could check them, and neither would have noticed them
going stale -- while the repo's own findings page lists "the replication
collapse was my own censoring artifact" as an established result.

These guards pin the property, not the numbers: the series exists in the
artifact, the two documents agree with it, and it still shows what it claims to.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "analysis" / "atlas-replication.json"
DOC = REPO_ROOT / "analysis" / "atlas-replication.md"
FINDINGS = REPO_ROOT / "analysis" / "census-findings.md"


def _raw():
    return json.loads(RAW.read_text())


def _complete(d):
    """Cohorts old enough to have completed the follow-up window."""
    return [r for r in d["cohorts_censored_ever"]
            if r["year"] + d["quiet_years"] <= d["latest_year"]]


def test_the_censored_series_is_stored_not_remembered():
    d = _raw()
    assert d.get("cohorts_censored_ever"), (
        "the ever-replicated series the methodology section argues from is not in "
        "the artifact, so nothing can check it")
    assert len(_complete(d)) >= 10, "too few complete cohorts to support the claim"


def _prose():
    """The censoring PARAGRAPH, not the whole document.

    Checking substring-anywhere let the prose claim 99.9% while the test stayed
    green, because the generator also prints the same cohort into the table
    below it. The endpoints have to be checked where they are asserted.
    """
    txt = DOC.read_text()
    start = txt.index("Scoring each cohort on whether it was *ever* replicated")
    end = txt.index("| first asserted |", start)
    return txt[start:end]


def test_the_document_quotes_the_series_it_stores():
    """Both endpoints in the prose must come from the artifact."""
    d, prose = _raw(), _prose()
    rows = _complete(d)
    assert f"{100*rows[0]['rate']:.1f}% for pairs first asserted in {rows[0]['year']}" \
        in prose, f"the {rows[0]['year']} endpoint in the prose is not the artifact's"
    assert f"{100*rows[-1]['rate']:.1f}% for {rows[-1]['year']}" in prose, (
        f"the {rows[-1]['year']} endpoint in the prose is not the artifact's")


def test_the_derived_arithmetic_is_right():
    """The two quantities that carry the argument were the two unguarded ones.

    "a fall of N points across M years" and "has had N years" are what make the
    censoring case; both survived every earlier assertion here, so a fall of
    99.9 points across 3 years would have shipped.
    """
    d, prose = _raw(), _prose()
    rows = _complete(d)
    lo, hi = rows[0], rows[-1]
    assert f"a fall of {100*(lo['rate'] - hi['rate']):.1f} points" in prose
    assert f"across {hi['year'] - lo['year']} years" in prose
    assert f"has had {d['latest_year'] - lo['year']} years" in prose
    assert f"has had {d['latest_year'] - hi['year']}," in prose


def test_the_artifact_is_internally_consistent():
    """rate must be replicated/pairs, and the >=200 filter must hold."""
    d = _raw()
    for key in ("cohorts", "cohorts_censored_ever"):
        for r in d[key]:
            assert r["pairs"] >= 200, f"{key} {r['year']} is below the reporting floor"
            assert abs(r["rate"] - r["replicated"] / r["pairs"]) < 1e-9, (
                f"{key} {r['year']}: rate does not equal replicated/pairs")
            assert r["replicated"] <= r["pairs"]


def test_the_findings_page_quotes_the_artifact_not_the_other_document():
    """census-findings.md restated this from atlas-replication.md's PROSE."""
    d = _raw()
    rows = _complete(d)
    txt = FINDINGS.read_text()
    assert f"{100*rows[0]['rate']:.1f}% for {rows[0]['year']}" in txt, (
        "the findings page quotes a censored rate the artifact does not hold")
    assert f"{100*rows[-1]['rate']:.1f}% for {rows[-1]['year']}" in txt


def test_the_comparison_still_shows_what_it_claims():
    """If the censored series stopped declining, the argument would be gone.

    This is the assertion that makes the others worth having: the section
    claims a steep apparent decline that the equal-window series does not
    reproduce. Both halves are checked, so a change that flattened the censored
    series -- or steepened the windowed one -- would surface here rather than
    leaving the prose asserting a contrast that no longer exists.
    """
    d = _raw()
    cen = _complete(d)
    assert cen[0]["rate"] - cen[-1]["rate"] > 0.20, (
        "the censored series no longer shows the steep decline the section is "
        f"built on ({100*cen[0]['rate']:.1f}% -> {100*cen[-1]['rate']:.1f}%)")
    # NOT "the windowed series falls less steeply" -- a review proved that
    # cannot fail. Within-W implies ever, so the censored rate dominates
    # everywhere, and at the newest complete cohort the two are IDENTICAL by
    # construction, so the difference in falls reduces to the difference at the
    # old end and is non-negative for free.
    #
    # What is worth pinning is the thing that identity implies and the prose now
    # states: the two series MEET at the recent end. If they ever stopped
    # meeting, the window arithmetic would be wrong.
    win = {r["year"]: r for r in d["cohorts"]}
    newest = cen[-1]
    assert newest["year"] in win, "the newest complete cohort is missing from the window series"
    assert abs(win[newest["year"]]["rate"] - newest["rate"]) < 1e-9, (
        f"the two series disagree at {newest['year']}, where a pair can only be "
        "replicated inside its own window -- the window arithmetic is wrong")
    # And the separation that carries the argument is at the OLD end.
    assert cen[0]["rate"] - win[cen[0]["year"]]["rate"] > 0.20, (
        "the censored and windowed measures no longer separate at the old end, "
        "which is where the entire contrast lives")


def test_the_degenerate_tail_is_labelled_when_it_is_shown():
    """The final cohort has had ~no time; showing its rate needs the caveat."""
    d, txt = _raw(), DOC.read_text()
    tail = d["cohorts_censored_ever"][-1]
    if tail["year"] + d["quiet_years"] <= d["latest_year"]:
        return                                  # not degenerate, no caveat needed
    # Keyed on the SENTENCE, not on the rate string. Gating on "is this rate
    # present" made the test return early in exactly the case it exists for --
    # a stale document whose quoted rate no longer matches the artifact.
    m = re.search(r"reaches ([\d.]+)% at (\d+), which is not evidence of anything: "
                  r"that cohort has had (\d+) years", txt)
    assert m, ("the final cohort's rate is shown without the caveat that it has had "
               "no time to be replicated")
    assert abs(float(m.group(1)) - 100 * tail["rate"]) < 0.06, (
        f"the document says {m.group(1)}%, the artifact says {100*tail['rate']:.1f}%")
    assert int(m.group(2)) == tail["year"]
    assert int(m.group(3)) == d["latest_year"] - tail["year"]
