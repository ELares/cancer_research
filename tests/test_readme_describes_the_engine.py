"""The front door's description of the simulation engine, bound to the engine.

README.md said the worked implementations were "ferroptosis/RSL3 biochemistry
and PDT/SDT depth physics -- there are no others" long after eight more arms
had landed, and named radiotherapy as its worked example of a modality with
"no query, no mechanism tag, no engine term" after `Treatment::Radiation`
shipped. Both were true when written. Neither was guarded, so neither moved.

That is this repository's most-repeated defect -- a claim in one artifact about
another -- and the front door is the worst place for it, because it is what a
reader meets first and the last thing anyone re-reads. The numbers the README
now uses to describe the engine are derived here from the artifacts that
produce them.

NOT a check that the prose is good. It cannot tell whether the description is
fair, only whether its figures are the live ones.
"""
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
README = REPO / "README.md"
ANALYSIS = REPO / "analysis"


WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
         7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven",
         12: "twelve", 13: "thirteen"}


def _j(name):
    return json.loads((ANALYSIS / name).read_text())


def _says(txt: str, n: int, template: str) -> bool:
    """Accept a count written as a digit OR spelled out.

    Prose reads better with words and artifacts carry digits; forcing one on
    the other makes the guard an editor rather than a check.
    """
    return any(template.format(n=form) in txt
               for form in {str(n), WORDS.get(n, str(n))})


def test_the_readme_states_the_live_binary_and_arm_counts():
    txt = README.read_text()
    binaries = len([p for p in (REPO / "simulations").glob("sim-*") if p.is_dir()])
    assert f"**{binaries} Rust simulation binaries**" in txt, (
        f"the README does not say {binaries} binaries; the tree holds that many")
    panel = _j("modality-panel.json")
    arms = [a["arm"] for a in panel["arms"]]
    treatments = [a for a in arms if a != "Control"]
    assert _says(txt, len(treatments),
                 "**{n} treatment arms plus an untreated control**"), (
        f"the README does not state the live arm count ({len(treatments)} plus "
        "a control); `sim-modality-panel` runs "
        f"{sorted(arms)}")


def test_the_readme_quotes_the_measured_depth_gap():
    """The honest half. If the arms ever catch up, this sentence must be
    rewritten rather than left standing -- and if they shrink, likewise."""
    txt = README.read_text()
    d = _j("modality-module-depth.json")
    assert (f"{d['dedicated_modules']} modules and {d['dedicated_code_lines']} lines"
            in txt), "the README does not quote the live dedicated-module size"
    assert (f"{d['ferroptosis_engine_modules']} modules and\n  "
            f"{d['engine_code_lines']:,}" in txt
            or f"{d['ferroptosis_engine_modules']} modules and "
               f"{d['engine_code_lines']:,}" in txt), (
        "the README does not quote the live ferroptosis-engine size")
    lo, hi = round(d["line_ratio_narrow"]), round(d["line_ratio_wide"])
    assert f"roughly {lo} to {hi} times smaller" in txt, (
        f"the README says something other than the live {lo}-to-{hi}x gap "
        f"(exact ratios {d['line_ratio_narrow']} and {d['line_ratio_wide']})")
    assert d["engine_code_lines"] > d["dedicated_code_lines"], (
        "the modality arms now exceed the ferroptosis engine; the README's "
        "'smaller' framing has to be re-derived, not adjusted")


def test_the_readme_findings_quote_the_panel_and_the_barriers():
    txt = README.read_text()
    panel = _j("modality-panel.json")
    arms = {a["arm"]: a["kill_fraction"] for a in panel["arms"]}
    adc, sdt = arms["AntibodyDrugConjugate"], arms["SDT"]
    assert f"kills **{adc * 100:.1f}%** where sonodynamic therapy kills **{sdt * 100:.1f}%**" in txt, (
        f"the README's delivery finding does not quote the live pair "
        f"({adc * 100:.1f}% vs {sdt * 100:.1f}%)")
    assert adc < sdt, "the delivery argument depends on the ADC being smaller"
    ab = panel["adoptive_barriers"]
    collapse = ab["leukaemia_kill_fraction"] / ab["solid_tumour_kill_fraction"]
    assert f"collapses {collapse:,.0f}-fold" in txt, (
        f"the README does not quote the live {collapse:,.0f}-fold CAR-T collapse")
    preds = _j("modality-predictions.json")["P10"]["reach_of_negative_pool"]
    hi_k, lo_k = max(preds, key=float), min(preds, key=float)
    assert (f"falls from {float(preds[hi_k]):.0%} to {float(preds[lo_k]):.1%}"
            in txt), "the README does not quote P10's live reach decline"


def test_the_readme_scope_claims_match_the_audit():
    """The two sentences that were FALSE before this update, pinned so they
    cannot go false again in the other direction."""
    txt = README.read_text()
    sa = _j("scope-audit.json")
    cov = _j("modality-coverage.json")
    silent = sorted(s.replace(".rs", "") for s in sa["engine_modules"]["silent"])
    assert _says(txt, len(silent), "{n} modules ("), (
        f"the README does not state the live silent-module count ({len(silent)})")
    for m in silent:
        assert f"`{m}`" in txt, f"the README does not name the silent module {m}"
    absent = [m for m in cov["mechanisms"] if not m.get("code_modules")]
    assert f"**{len(absent)} now have no engine representation at all**" in txt, (
        f"the README does not state the live absent count ({len(absent)} of "
        f"{len(cov['mechanisms'])})")
    assert "Radiotherapy was this section's worked example" in txt, (
        "the README no longer records that its radiotherapy example was "
        "retracted; a reader cannot tell the claim changed")
    cal = _j("modality-calibration.json")
    vc = cal["verdict_counts"]
    assert (f"{vc['ADMISSIBLE']} of the {len(cal['arms'])} arms" in txt), (
        "the README does not quote the live calibration verdict split")


def test_the_readme_does_not_re_assert_the_retracted_claims():
    """Both sentences this update replaced, refused by name."""
    txt = README.read_text()
    for gone in ("there are no others",
                 "no query, no mechanism tag, no engine term",
                 "**12 Rust simulation binaries**",
                 "for **ferroptosis and physical-ROS therapies specifically**"):
        assert gone not in txt, (
            f"the README asserts the retracted claim {gone!r} again")
