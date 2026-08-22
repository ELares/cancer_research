"""Guards for the precedent measurement on this project's keystone experiment.

The P1 protocol exists to be handed to a collaborator. This measures how much
comparative literature each of its reagents has, and finds the two arms 17.7x
apart -- abundant precedent for choosing an RSL3 dose, almost none for choosing
an iFSP1 one.

TWO HAZARDS, in opposite directions.

A count of mentions is not evidence about quality, and a thin literature is not
a reason to avoid a reagent: a compound is rare when it is NEW for the same
reason it is rare when it is poor. If this analysis ever reads as an argument
against the experiment, it has been misused.

And the finding must reach the protocol itself, not sit in an analysis
directory. A collaborator reads the protocol; a caveat they never see is a
caveat that does not exist.

OFFLINE: reads only committed artifacts.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-protocol-precedent.json"
MD = REPO / "analysis/census-protocol-precedent.md"
PROTOCOL = REPO / "analysis/p1-wetlab-protocol.md"
THIN = 50


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_reagents_measured_are_the_ones_the_protocol_names(d):
    """Read from the protocol, not typed into the analysis.

    If the protocol switches reagents, the measurement must follow it rather
    than keep reporting on compounds nobody will use.
    """
    txt = PROTOCOL.read_text().lower()
    for name in d["protocol_named"]:
        assert re.search(rf"\b{re.escape(name.lower())}\b", txt), (
            f"{name} is marked as named in the protocol and is not in it")
    # Both arms must have at least one named reagent, or the asymmetry
    # comparison is between an arm and nothing.
    for g in d["groups"]:
        if g["role"].endswith(("(arm 1)", "(arm 2)")):
            assert any(r["in_protocol"] for r in g["rows"]), (
                f"{g['role']} has no protocol-named reagent")


def test_the_asymmetry_recomputes_from_the_arm_totals(d):
    assert d["arm1_precedent"] == next(
        g["protocol_total"] for g in d["groups"] if g["role"].endswith("(arm 1)"))
    assert d["arm2_precedent"] == next(
        g["protocol_total"] for g in d["groups"] if g["role"].endswith("(arm 2)"))
    assert d["arm_asymmetry"] == pytest.approx(
        d["arm1_precedent"] / d["arm2_precedent"], abs=0.05)
    # The thin side must be derived, not asserted.
    thin_is_arm2 = d["arm_asymmetry"] > 1
    assert d["thin_arm"].endswith("(arm 2)") == thin_is_arm2


def test_the_finding_reached_the_protocol_itself(d):
    """A collaborator reads the protocol. A caveat they never see does not exist."""
    txt = " ".join(PROTOCOL.read_text().split())
    assert str(d["arm_asymmetry"]) in txt, (
        "the protocol does not state the asymmetry between its two arms")
    assert f"{d['arm1_precedent']:,}" in txt and f"{d['arm2_precedent']:,}" in txt
    # And it must say what to DO about it, not merely note it.
    assert "single-agent pre-run matters far more" in txt or \
           "Budget for it" in txt, (
        "the protocol states the asymmetry without a consequence, which leaves "
        "a collaborator no action to take")
    assert "harder to interpret" in txt, (
        "the protocol does not say that a negative result on the thin arm is "
        "harder to interpret, which is the consequence that matters most")


def test_it_does_not_argue_against_the_experiment(d):
    """A compound is rare when it is NEW for the same reason it is rare when it
    is poor, and a count cannot tell those apart."""
    # WHITESPACE-NORMALISED. The phrase wraps across a line break in the
    # protocol, and a raw substring search misses it -- the wrapped-prose blind
    # spot this project has already recorded, hitting a guard written to
    # enforce a caveat.
    for text in (MD.read_text(), PROTOCOL.read_text()):
        flat = " ".join(text.split())
        assert "rare in the literature when it is new for the same reason" in flat
    md = " ".join(MD.read_text().split())
    assert "does not mean the experiment is wrong" in md
    for overclaim in ("should not be used", "abandon", "unsuitable reagent",
                      "the experiment is not viable"):
        assert overclaim not in md.lower()


def test_thin_reagents_are_flagged_only_where_the_protocol_names_them(d):
    """Flagging every rare compound would fill the table with noise: most
    compounds in any field are rare, and only the ones this protocol asks a
    collaborator to buy matter here."""
    md = " ".join(MD.read_text().split())
    for g in d["groups"]:
        for r in g["rows"]:
            assert r["thin"] == (r["articles"] < THIN)
            if r["thin"] and r["in_protocol"]:
                assert f"| {r['reagent']} *" in md, (
                    f"{r['reagent']} is thin and protocol-named but unmarked")


def test_the_undercount_direction_is_stated(d):
    """Title-and-abstract counting misses Methods-only mentions. That biases
    every row the same way, which is why the between-row comparison survives --
    and saying so is what stops a reader treating a count as a census of use."""
    md = " ".join(MD.read_text().split())
    assert "never named in the abstract is undercounted" in md
    assert "biases every row the same way" in md
