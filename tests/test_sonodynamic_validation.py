"""Guards for the sonodynamic frequency-optimum validation.

The load-bearing risk here is not that the model is wrong -- one of its three
claims IS wrong and the page says so. It is that the refutation quietly turns
into a confirmation on some later regeneration, because a refuted claim is the
uncomfortable one and the document has to keep carrying it.

So these guards pin the REFUTATION as hard as the confirmations, and pin the
things that would make either hollow: that the page is regenerated rather than
committed stale, that the closed form is checked against a scan rather than
against itself, and that the two mismatches the page refuses to reason across
are still stated.
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_sonodynamic.py"
MD = REPO / "analysis" / "calibration" / "sonodynamic-validation.md"
JSON = REPO / "analysis" / "calibration" / "sonodynamic-validation.json"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "sonodynamic.rs"
PARAMS = REPO / "simulations" / "ferroptosis-core" / "src" / "params.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_page_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, (
            "sonodynamic-validation.md is stale; run scripts/validate_sonodynamic.py")
        assert JSON.read_text() == before_json, "the JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_python_closed_form_is_the_rust_one():
    """A drift guard, parsed from source rather than restated.

    Two implementations of one formula is the arrangement this repository
    uses to make a validator independent; it only works while they agree.
    """
    src = RUST.read_text()
    assert "10.0 / denom" in src, "the Rust closed form has changed shape"
    assert "alpha_db_cm_mhz * std::f64::consts::LN_10 * depth_cm" in src, (
        "the Rust denominator has changed; the Python re-implementation in "
        "validate_sonodynamic.py must move with it")
    # And the attenuation constant the page quotes comes from params.rs, not
    # from a literal in the validator.
    alpha = float(re.search(r"sdt_alpha:\s*([0-9.]+)", PARAMS.read_text()).group(1))
    assert _d()["alpha_db_cm_mhz_from_params_rs"] == alpha


def test_the_closed_form_is_checked_against_a_scan_and_not_against_itself():
    """The check that makes claim 1 a derivation.

    A guard whose expected value comes from the thing it is guarding proves
    nothing; here the expectation is an independent brute scan.
    """
    for row in _d()["interior"]:
        assert row["agrees"], f"closed form and scan disagree at {row['depth_cm']} cm"
        assert row["beats_low_end"] and row["beats_high_end"], (
            f"the optimum at {row['depth_cm']} cm is not interior -- a "
            "monotonic model would satisfy the closed form and still be wrong")
        assert row["scanned_mhz"] > 0


def test_the_depth_scaling_is_still_reported_as_refuted():
    """The uncomfortable claim, pinned.

    If the comparator's own optimum ever stops being flat, or the model's
    stops sliding, this must be re-derived rather than left saying REFUTED --
    and equally, a regeneration must not quietly upgrade it.
    """
    d = _d()
    assert d["claim_depth_scaling"] == "REFUTED"
    assert d["model_depth_span_ratio"] > 2.0, (
        "the model no longer predicts a large depth swing, so the refutation "
        "no longer describes what the model does")
    assert all(v < 1.6 for v in d["reported_depth_span_ratios"].values()), (
        "the comparator's optimum is no longer flat with depth")
    md = MD.read_text()
    assert "no value of" in md and "flat" in md, (
        "the page no longer states that this is a disagreement rather than a "
        "tuning gap")
    assert "near-field" in md.lower(), (
        "the page no longer names the missing term, which is the only useful "
        "thing a refutation leaves behind")


def test_the_two_confirmed_claims_are_confirmed_for_a_reason():
    d = _d()
    assert d["claim_interior_optimum"] == "CONFIRMED"
    assert d["claim_falls_with_attenuation"] == "CONFIRMED"
    # Not just the model's own direction: the comparator must independently
    # show it too, or claim 2 is the model agreeing with itself.
    assert d["reported_optimum_falls_with_attenuation"], (
        "the comparator no longer shows the attenuation direction, so claim 2 "
        "rests on the model alone")
    for row in d["attenuation"]:
        assert row["falls"]
        assert abs(row["ratio"] - 2.0) < 1e-6, (
            "doubling alpha no longer halves the optimum exactly; that exact "
            "factor is the sharper half of the claim")


def test_the_page_refuses_the_numerical_comparison_and_says_why():
    md = MD.read_text()
    for phrase in ("No numerical agreement", "opposite sign", "Np/m/MHz",
                   "factor of two"):
        assert phrase in md, (
            f"the page no longer states the mismatch {phrase!r}; a numeric "
            "match across either of them would be uninterpretable")
    assert d_verdict() == "PARTIAL"


def d_verdict():
    return _d()["verdict"]


def test_the_comparator_rows_are_the_published_ones():
    """The anchor cannot drift into being convenient.

    These six values are what PMID 26233216 reports; if any of them is edited
    the claim changes and the edit must be deliberate.
    """
    rows = {(e["alpha_np_m_mhz"], e["depth_mm"]): e["reported_khz"]
            for e in _d()["ellens"]}
    assert rows == {
        (5.0, 50): 750, (5.0, 100): 750, (5.0, 150): 750,
        (10.0, 50): 750, (10.0, 100): 500, (10.0, 150): 500,
    }


def test_the_model_answers_inside_the_band_the_comparator_scanned():
    """Weaker than a fit and stated as such: a model whose optimum fell
    outside 250-1500 kHz everywhere would not be answering the same question,
    and 'the directions agree' would mean nothing."""
    band = _d()["band"]
    assert band and all(r["in_scanned_band"] for r in band)
    # Straddling: the model must be above the comparator at the shallow end
    # and below it at the deep end, which is what a sliding optimum against a
    # flat one MEANS. Both on one side would be an offset, not a scaling gap.
    khz = {r["depth_mm"]: r["model_khz"] for r in band}
    assert khz[50] > 750 > khz[150], (
        f"the model no longer straddles the comparator's 750 kHz: {khz}")


def test_the_engine_module_carries_the_threshold_distinction():
    """SDT's limit is a threshold, which is what makes it structurally unlike
    the dose-response arms. If `cavitates` ever becomes a gradient the whole
    section's argument changes."""
    src = RUST.read_text()
    assert "pub fn cavitates(" in src
    assert re.search(r"index >= threshold", src), (
        "cavitation is no longer a threshold comparison")
    assert "MI_DIAGNOSTIC_CAP" in src
    assert "REGULATORY LIMIT, NOT A MEASUREMENT" in src, (
        "the diagnostic cap no longer says it is a regulatory ceiling rather "
        "than the pressure at which tissue cavitates")


def test_the_closed_form_reproduces_by_hand():
    """One arithmetic check the guards do not take from the artifact."""
    # f* = 10 / (alpha * ln10 * z); at alpha=1, z=1 that is 10/ln(10).
    from importlib import util
    spec = util.spec_from_file_location("vs", SCRIPT)
    m = util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert abs(m.optimal_frequency_mhz(1.0, 1.0) - 10.0 / math.log(10.0)) < 1e-12
    assert math.isinf(m.optimal_frequency_mhz(0.0, 1.0))


def test_the_depth_limit_table_marks_both_degenerate_rows():
    """Two rows in this table are not measurements, and both were shipped as
    measurements once before being caught.

    A row whose applicator cannot cavitate at the SURFACE fails everywhere,
    and rendering its zero as a shallow limit is the defect the oncolytic and
    ablation sections each fixed. A row whose limit equals the SCAN BOUND
    never found a limit at all, and 30.00 cm read as an extraordinarily deep
    applicator when it measures how far the scan went. Each marker must have a
    live example, because a rule nothing triggers is a rule nobody can check.
    """
    d = _d()
    rows = d["depth_limits"]
    assert rows, "the depth-limit table is empty"
    scan = d["depth_scan_limit_cm"]
    for r in rows:
        assert r["depth_limit_cm"] >= 0.0
        assert r["unbounded_in_scan"] == (r["depth_limit_cm"] >= scan - 1e-9)
        if r["total_failure"]:
            assert r["depth_limit_cm"] == 0.0, (
                "an applicator that fails at the surface reports a non-zero "
                "depth, which is the exact defect the marker exists for")
    assert any(r["total_failure"] for r in rows), (
        "no row exercises the TOTAL FAILURE marker")
    assert any(r["unbounded_in_scan"] for r in rows), (
        "no row exercises the beyond-the-scan marker")
    assert any(not r["total_failure"] and not r["unbounded_in_scan"]
               for r in rows), (
        "every row is degenerate, so the table measures nothing")
    md = MD.read_text()
    assert "*fails everywhere*" in md and "*beyond the scan*" in md, (
        "the page no longer distinguishes the two degenerate rows in its own "
        "rendering, whatever the JSON says")


def test_the_depth_limit_is_monotone_in_the_things_it_should_be():
    """A stronger applicator cannot reach less deep, and a higher device cap
    cannot hurt -- the cap binds only at shallow depth where the unconstrained
    optimum is high, so raising it can add reach or do nothing."""
    rows = {(r["index_at_reference"], r["max_frequency_mhz"]): r
            for r in _d()["depth_limits"]}
    strengths = sorted({k[0] for k in rows})
    caps = sorted({k[1] for k in rows})
    for cap in caps:
        for a, b in zip(strengths, strengths[1:]):
            assert rows[(b, cap)]["depth_limit_cm"] >= rows[(a, cap)]["depth_limit_cm"], (
                f"a stronger applicator reaches less deep at cap {cap}")
    for st in strengths:
        for a, b in zip(caps, caps[1:]):
            assert rows[(st, b)]["depth_limit_cm"] >= rows[(st, a)]["depth_limit_cm"], (
                f"a higher device cap reduced reach at strength {st}")


def test_the_device_cap_exists_because_the_model_diverges_without_it():
    """The cap is a physical property and not a fudge, and the page has to
    keep saying which. A reader who takes it for a tuning knob would read the
    depth limits as tunable, when what it removes is an idealisation."""
    src = RUST.read_text()
    assert "pub fn usable_frequency_mhz(" in src
    assert "THE CAP IS NOT A FUDGE AND THE MODEL IS WRONG WITHOUT IT" in src
    md = MD.read_text()
    assert "unbounded at shallow depth" in md
    assert "REQUIRED rather than" in md, (
        "the page no longer says the device cap is required rather than "
        "assumed")


def test_the_figure_reads_the_curve_rather_than_recomputing_it():
    """Two implementations of one formula in one repository is a drift this
    project has already had to fix, and a figure disagreeing with the page
    beside it is the worst place for it."""
    d = _d()
    curves = d["curves"]
    assert curves, "the artifact carries no curve for the figure to read"
    gen = (REPO / "scripts" / "generate_conceptual_diagrams.py").read_text()
    start = gen.index("def fig43_sonodynamic_frequency")
    body = gen[start:gen.index("\ndef ", start + 10)]
    assert 'd["curves"]' in body, (
        "fig43 no longer reads the committed curve")
    assert 'd["cavitation_threshold"]' in body, (
        "fig43 hardcodes a threshold instead of reading the artifact's")
    # Every curve must start at a FINITE surface value; that is what the cap
    # buys, and an infinite one would mean the divergence is back.
    for c in curves:
        assert c["points"][0][0] == 0.0
        assert 0.0 < c["points"][0][1] < 1e3, (
            f"curve at strength {c['index_at_reference']} is not finite at "
            "the surface; the device cap is not being applied")
        ys = [pt[1] for pt in c["points"]]
        assert all(b <= a + 1e-12 for a, b in zip(ys, ys[1:])), (
            "the delivered index does not fall monotonically with depth")


def test_the_chapter_table_is_the_artifact_table():
    """Section 6.13 reprints the depth-limit grid by hand.

    A claim in one artifact about the contents of another is this
    repository's recurring defect class, and the same PR had a per-arm line
    count go stale within one commit. Every cell is checked, including the two
    that are MARKERS rather than numbers -- those are the ones a hand-edit
    would most plausibly "tidy" into a value, which is precisely the reading
    the markers exist to prevent.
    """
    md = (REPO / "article" / "drafts" / "v1.md").read_text()
    i = md.index("### 6.13 Photodynamic")
    section = md[i:md.index("### 6.14")]
    assert "Deepest focus still clearing the cavitation threshold" in section, (
        "the chapter no longer reprints the depth-limit table")
    for r in _d()["depth_limits"]:
        if r["total_failure"]:
            expected = "fails everywhere"
        elif r["unbounded_in_scan"]:
            expected = "beyond the 30 cm scan"
        else:
            expected = f"{r['depth_limit_cm']:.2f} cm"
        assert expected in section, (
            f"the chapter's table is missing {expected!r} for applicator "
            f"{r['index_at_reference']} at a {r['max_frequency_mhz']} MHz cap")
    # Every row of the artifact must appear, so a row cannot be dropped
    # silently -- one was, in the first draft of this table.
    body = section[section.index("| applicator |"):]
    body = body[:body.index("\n\n")]
    assert len(body.strip().splitlines()) == 2 + len(
        {r["index_at_reference"] for r in _d()["depth_limits"]}), (
        "the chapter's table has a different number of rows from the "
        "artifact's applicator strengths")
