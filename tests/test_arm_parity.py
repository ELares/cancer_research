"""Guards for the per-arm parity table.

The table exists to say how far each treatment arm is from the ferroptosis arm.
It will be read while the gap is being closed, which makes it exactly the kind
of page that flatters its author if nothing checks it: every axis is a number
this project controls, and every one of them can be moved by editing prose
instead of code.

So these guards check the three things that would make it useless:

  * that it is FRESH -- regenerated from the crate and the artifacts, not a
    committed snapshot that agreed with the code six weeks ago (this repository
    has shipped that defect six times, see tests/test_generator_guards_pin_source.py)
  * that the ATTRIBUTION RULE does what the page says it does, which is where
    the first version was wrong in the direction that flattered the arms
  * that an arm the engine does NOT have cannot quietly disappear from the page
    by being built, or by being deleted from a list
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "arm_parity.py"
MD = REPO / "analysis" / "arm-parity.md"
JSON = REPO / "analysis" / "arm-parity.json"
CELL_RS = REPO / "simulations" / "ferroptosis-core" / "src" / "cell.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_page_is_what_the_script_produces_now():
    """A regenerate-and-diff gate, not an artifact-reads-artifact one.

    Every other check here reads the committed JSON. If the generator changed
    and nobody re-ran it, they would all pass against a stale file -- the
    artifact and its guards go stale together, which is the failure mode this
    repository keeps rediscovering.
    """
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, (
            "analysis/arm-parity.md is stale: re-running the generator changes "
            "it. Run scripts/arm_parity.py and commit the result.")
        assert JSON.read_text() == before_json, (
            "analysis/arm-parity.json is stale; run scripts/arm_parity.py")
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_every_arm_in_the_table_is_a_treatment_the_engine_has():
    """Derived from the enum, so an arm cannot be listed into existence."""
    body = re.search(r"pub enum Treatment\s*\{(.*?)\n\}", CELL_RS.read_text(),
                     re.S).group(1)
    variants = {v.strip().rstrip(",") for v in body.splitlines()
                if v.strip() and not v.strip().startswith(("//", "#["))}
    listed = {a["arm"] for a in _d()["arms"]}
    assert listed <= variants, (
        f"the parity table lists arms the engine does not have: "
        f"{sorted(listed - variants)}")
    missing = variants - listed - {"Control"}
    assert not missing, (
        f"the engine has arms the parity table does not measure: "
        f"{sorted(missing)}. An arm that is not in the table is an arm whose "
        "distance from parity nobody is tracking.")


def test_the_comparator_is_still_the_largest_arm():
    """The success condition, stated as a guard.

    If an arm ever exceeds the ferroptosis engine, this page's whole framing --
    a distance measured TOWARDS a comparator -- stops being the right shape and
    has to be re-derived rather than adjusted.
    """
    d = _d()
    base = next(r for r in d["arms"] if r["arm"] == d["comparator"])
    bigger = [r["arm"] for r in d["arms"]
              if r["arm"] != d["comparator"] and r["code_lines"] > base["code_lines"]]
    assert not bigger, (
        f"{bigger} now exceed the comparator. That is the campaign succeeding, "
        "and it means the page must be rewritten around a different frame -- "
        "not that this guard should be relaxed.")


def test_the_book_column_credits_whole_chapters_to_the_arm_they_are_about():
    """The attribution rule, checked at the case that caught it out.

    Reading only ### sections credited the COMPARATOR with zero book words,
    because Chapter 5 is titled "The Ferroptosis Engine" and none of its
    section headings repeat the word. Zero for ferroptosis and two thousand
    for checkpoint blockade is not a small error: it is the table's own
    headline claim, inverted.
    """
    d = _d()
    base = next(r for r in d["arms"] if r["arm"] == d["comparator"])
    assert base["book_words"] > 1000, (
        "the comparator is credited with almost no manuscript words, which "
        "means the attribution rule is reading sections and not chapters")
    assert any("Chapter 5" in s for s in base["book_sections"]), (
        "Chapter 5 (The Ferroptosis Engine) is not attributed to the "
        "ferroptosis arm; the chapter-level rule is not firing")
    assert d["manuscript_words_attributed"] < d["manuscript_words_in_sections"], (
        "every word in the manuscript is attributed to an arm, which cannot be "
        "true -- the unattributed remainder is what makes the column honest")


def test_an_arm_without_a_taxonomy_row_is_not_reported_as_zero():
    """`no row` and `0` are different claims, and this repository has already
    published the second when it meant the first."""
    d = _d()
    for r in d["arms"]:
        if r["mechanism"] is None:
            assert r["census"] is None, (
                f"{r['arm']} has no taxonomy row but carries a census count")
    md = MD.read_text()
    assert "no row" in md, (
        "the literature table no longer distinguishes an unmeasurable arm from "
        "one with no articles")


def test_the_arms_the_engine_lacks_are_named_and_still_lacking():
    """A parity page listing only what exists measures progress and hides
    scope. Each absent arm is checked to be genuinely absent, so the list
    cannot outlive the gap it describes."""
    d = _d()
    names = {a["arm"] for a in d["absent_arms"]}
    assert "Cytotoxic chemotherapy" in names, (
        "chemotherapy has dropped off the absent list; it is the modality most "
        "patients receive and the sharpest thing this engine cannot express")
    body = re.search(r"pub enum Treatment\s*\{(.*?)\n\}", CELL_RS.read_text(),
                     re.S).group(1).lower()
    for token, arm in (("chemo", "Cytotoxic chemotherapy"),
                       ("hormone", "Hormone therapy"),
                       ("radioligand", "Radioligand therapy")):
        if token in body:
            assert arm not in names, (
                f"the engine now has a {token} variant but the page still "
                f"lists {arm} as absent -- move it into ARMS and measure it")


def test_the_page_states_what_it_cannot_measure():
    md = MD.read_text()
    for phrase in ("does not measure", "lower bound", "attributed by HEADING"):
        assert phrase.lower() in md.lower(), (
            f"the page no longer states its own limit: {phrase!r}")
