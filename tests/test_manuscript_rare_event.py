"""Section 5.3's rare-event paragraph must match the sweep it describes.

WHY THIS FILE EXISTS
--------------------
The paragraph quotes about a dozen numbers -- two bounds, a rate, two interval
factors, a death count, a percentage error -- and every one of them is a
hand-written copy of something `analysis/rare-event-findings.json` computes. That
is the shape this repo has been burned by repeatedly: the analysis survives
probing and the sentence drawn from it does not, because the sentence has no
numbers attached to falsify it once the artifact moves.

So every figure is RECOMPUTED here from the committed artifact and asserted to
appear in the paragraph. Nothing is restated.

WHAT IT DOES WHEN THE SWEEP GROWS
---------------------------------
It fails. That is deliberate and is the main point. A 1e11 tier is running; when
those rows land, the largest-n row changes, the bounds tighten by another decade
and this test goes red until the prose is updated to match. A guard that silently
tolerated the new data would leave the manuscript quoting a ten-billion-cell
bound while the repository held a hundred-billion-cell one -- which is precisely
the failure it is here to prevent.

THE ANCHOR MUST FAIL LOUDLY
---------------------------
If the paragraph is reworded past recognition the locator raises rather than
returning an empty string, because a guard that no-ops on a missing anchor
protects nothing and looks identical to a passing one.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MS = REPO_ROOT / "article" / "drafts" / "v1.md"
FINDINGS = REPO_ROOT / "analysis" / "rare-event-findings.json"
SWEEP = REPO_ROOT / "analysis" / "rare-event-sweep.jsonl"

ANCHOR = "**Resolving the zeros.**"


def paragraph() -> str:
    txt = MS.read_text()
    if ANCHOR not in txt:
        raise AssertionError(
            f"{ANCHOR!r} is not in v1.md. The rare-event paragraph was renamed or "
            "removed; this guard cannot check prose it cannot find, and silently "
            "passing would be worse than failing.")
    i = txt.index(ANCHOR)
    return txt[i:txt.index("\n\n", i)]


def findings() -> dict:
    return json.loads(FINDINGS.read_text())


def rows() -> list:
    return [json.loads(l) for l in SWEEP.read_text().splitlines() if l.strip()]


def _sci(x: float) -> str:
    """Format as the paragraph does, e.g. 3.69e-10 / 7.134e-07."""
    return f"{x:.2e}", f"{x:.3e}"


def test_the_zero_event_bounds_match_the_artifact():
    """Both Glycolytic bounds, recomputed, must appear in the prose."""
    para, f = paragraph(), findings()
    for key in ("Glycolytic/Control", "Glycolytic/RSL3"):
        v = f[key]
        assert v["n_dead_at_largest_n"] == 0, (
            f"{key} now has events; it is no longer a zero-event bound and the "
            "paragraph's framing of it as resolution-limited is wrong")
        two, three = _sci(v["bound_or_rate"])
        assert two in para or three in para, (
            f"{key}'s bound {v['bound_or_rate']:.3e} is not quoted in the "
            f"paragraph (looked for {two!r} and {three!r})")


def test_the_largest_n_in_the_prose_matches_the_sweep():
    """The prose says 'ten billion'; the sweep must actually reach 1e10.

    When the 1e11 tier lands this fails, which is the intended behaviour.
    """
    f = findings()
    largest = max(v["largest_n"] for v in f.values())
    words = {int(1e6): "one million", int(1e7): "ten million",
             int(1e8): "a hundred million", int(1e9): "one billion",
             int(1e10): "ten billion", int(1e11): "a hundred billion"}
    assert largest in words, f"unexpected largest n {largest:.0e}"
    para = paragraph()
    assert words[largest] in para, (
        f"the sweep now reaches {largest:.0e} ({words[largest]}) but the paragraph "
        f"does not say so. Update Section 5.3 -- the manuscript is quoting a "
        f"smaller run than the repository holds.")
    assert f"{largest:,}" in para, (
        f"the paragraph should state the sweep's top sample size as {largest:,}")


def test_the_resolved_rate_and_its_interval_match():
    para, v = paragraph(), findings()["PersisterNrf2/Control"]
    assert v["classification"] == "resolved"
    assert f"{v['n_dead_at_largest_n']:,}" in para, "death count not quoted"
    two, three = _sci(v["bound_or_rate"])
    assert three in para or two in para, "resolved rate not quoted"
    lo, hi = v["poisson_ci"]
    assert f"{lo:.2e}" in para and f"{hi:.2e}" in para, "Poisson interval not quoted"
    assert f"{hi/lo:.2f}" in para, (
        f"the interval-width factor {hi/lo:.2f} is not quoted; it is the number "
        "that says what the extra cells actually bought")


def test_the_small_sample_error_is_stated_and_correct():
    """'the single-event estimate was 40% high' must be arithmetic, not a guess."""
    para = paragraph()
    p = sorted((r for r in rows() if r["phenotype"] == "PersisterNrf2"),
               key=lambda r: r["n_cells"])
    first, last = p[0], p[-1]
    drift = abs(first["death_rate"] - last["death_rate"]) / last["death_rate"]
    stated = re.search(r"was (\d+)% high", para)
    assert stated, "the paragraph no longer states the small-sample error"
    assert abs(int(stated.group(1)) - round(drift * 100)) <= 1, (
        f"paragraph says {stated.group(1)}% but the artifact gives "
        f"{drift*100:.0f}% ({first['death_rate']:.3e} -> {last['death_rate']:.3e})")
    # and the claim that the small interval still covered the truth
    assert first["poisson_ci_low"] <= last["death_rate"] <= first["poisson_ci_high"], (
        "the paragraph claims the single-event INTERVAL contained the resolved "
        "value; on this data it does not, so that sentence is now false")


def test_the_selectivity_claim_uses_the_real_reciprocals():
    """'one part in 2.7 billion rather than one part in 271,000'."""
    para = paragraph()
    f = findings()["Glycolytic/RSL3"]
    small = next(r for r in rows()
                 if r["phenotype"] == "Glycolytic" and r["treatment"] == "RSL3"
                 and r["n_cells"] == 10**6)
    now = 1.0 / f["bound_or_rate"]
    before = 1.0 / small["poisson_ci_high"]
    assert f"{now/1e9:.1f} billion" in para, (
        f"expected 'one part in {now/1e9:.1f} billion' from the artifact")
    assert f"{round(before, -3):,.0f}".replace(",", ",") in para or \
           f"{int(round(before/1000)*1000):,}" in para, (
        f"expected the prior bound's reciprocal ~{before:,.0f}")


def test_the_three_stated_limits_are_all_present():
    """The caveats are the load-bearing part; none may be dropped in an edit."""
    para = paragraph().lower()
    for phrase, why in (
        ("zero calibration targets", "the parameters are uncalibrated"),
        ("never interact", "independent cells are not a tumor"),
        ("nested", "the samples are not independent across n"),
    ):
        assert phrase in para, f"the paragraph dropped the caveat that {why}"


def test_the_bound_sentence_quotes_the_LARGEST_n_bounds_positionally():
    """Parse the actual sentence, because "appears somewhere" is not enough.

    TWO WEAKER VERSIONS OF THIS FAILED A MUTATION, both proven rather than
    suspected. Corrupting the first bound from 3.69e-10 to 3.69e-09 survived:

      1. `assert bound in para` passed, because the paragraph quotes the same
         value TWICE (the two Glycolytic conditions share a bound) and the
         untouched second copy satisfied the assertion;
      2. an exhaustive "every figure must appear somewhere in the artifact"
         check ALSO passed, because 3.69e-09 is a perfectly legitimate value --
         it is the bound at n=1e9. The artifact contains a bound at every scale,
         so membership proves nothing about whether the right one was quoted.

    The question is positional: does the sentence about the largest run quote
    the largest run's bounds? So the sentence is parsed and its captures are
    compared to the artifact, in order.
    """
    para, f = paragraph(), findings()
    m = re.search(
        r"upper bound falls to ([\d.]+e-\d+) and ([\d.]+e-\d+) respectively, "
        r"from ([\d.]+e-\d+) at\s+the million-cell point", para)
    assert m, ("could not parse the bound sentence in Section 5.3; it was "
               "reworded and this guard can no longer check it")
    got_ctrl, got_rsl3, got_small = m.groups()

    exp_ctrl = f"{f['Glycolytic/Control']['bound_or_rate']:.2e}"
    exp_rsl3 = f"{f['Glycolytic/RSL3']['bound_or_rate']:.2e}"
    small = [r for r in rows() if r["phenotype"] == "Glycolytic"
             and r["n_cells"] == 10**6]
    exp_small = f"{small[0]['poisson_ci_high']:.2e}"

    assert got_ctrl == exp_ctrl, (
        f"Glycolytic+Control bound: prose says {got_ctrl}, artifact says {exp_ctrl}")
    assert got_rsl3 == exp_rsl3, (
        f"Glycolytic+RSL3 bound: prose says {got_rsl3}, artifact says {exp_rsl3}")
    assert got_small == exp_small, (
        f"million-cell bound: prose says {got_small}, artifact says {exp_small}")


def test_the_resolved_sentence_quotes_its_own_numbers_positionally():
    """Same treatment for the resolved condition's rate and interval."""
    para, v = paragraph(), findings()["PersisterNrf2/Control"]
    m = re.search(r"a rate of ([\d.]+e-\d+) with an interval\s+spanning a factor "
                  r"of\s+([\d.]+) \(\[([\d.]+e-\d+), ([\d.]+e-\d+)\]\)", para)
    assert m, "could not parse the resolved-rate sentence in Section 5.3"
    rate, factor, lo, hi = m.groups()
    lo_e, hi_e = v["poisson_ci"]
    assert rate == f"{v['bound_or_rate']:.3e}", f"rate {rate} vs {v['bound_or_rate']:.3e}"
    assert lo == f"{lo_e:.2e}" and hi == f"{hi_e:.2e}", "interval endpoints differ"
    assert factor == f"{hi_e/lo_e:.2f}", f"width factor {factor} vs {hi_e/lo_e:.2f}"


def test_the_paragraph_does_not_conflate_the_two_zero_event_bounds():
    """It quoted 3.689/n and called it the rule of three, which is 3/n.

    They differ by 23%. Both are proportional to 1/n, so the shape argument was
    unaffected and the error survived three readings of the figure -- including
    two where the markers were visibly sitting above their own reference line.
    """
    import math
    para = paragraph()
    quoted = re.findall(r"3\.69e-\d+", para)
    assert quoted, "the paragraph no longer quotes the zero-event bound"
    # If it quotes 3.689/n it must not also claim to be tracking 3/n.
    bad = re.search(r"tracks 3/n|is the rule of three|rule of three and\s+therefore", para)
    assert not bad, (
        f"the paragraph quotes {quoted[0]} (which is -ln(0.025)/n = 3.689/n) while "
        f"calling it the rule of three (3/n). They differ by 23%.")
    # and the distinction has to be stated, not merely avoided
    assert "3.689/n" in para, "the paragraph should name the constant it quotes"
    assert "rule of three" in para, (
        "the paragraph should say explicitly that this is NOT the rule of three, "
        "since that is the number a reader will assume")
    n = 10**11
    assert abs((-math.log(0.025) / n) / (3.0 / n) - 1.2296) < 1e-3
