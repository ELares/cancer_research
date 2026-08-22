"""Guards for the sizing of the literatures MeSH cannot name.

These nine rows are the weakest measurements in the project and they look like
the strongest: precise counts in a table, with no visible mark of the fact that
every other modality here can be cross-checked against a second instrument and
these cannot. The cross-check is unavailable BY CONSTRUCTION -- the absence of a
descriptor is what puts a modality on this list -- so the guards exist to keep
that asymmetry stated rather than implied.

The second hazard is the opposite of the first: having produced numbers, it
becomes tempting to drop the phrase they qualify. "Not measurable" is still
correct and still needs saying; what changes is that it means "not measurable
BY DESCRIPTOR", which is a claim about MeSH rather than about the field.

OFFLINE: reads only committed artifacts.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-unnamed-modalities.json"
MD = REPO / "analysis/census-unnamed-modalities.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
# Pinned literal, not read from the artifact.
VINTAGE_YEAR = 1990


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_every_row_genuinely_lacks_a_descriptor(d):
    """The list must be derived from the descriptor map, not curated.

    A modality that acquires a descriptor should leave this table and become
    cross-checkable; one hand-maintained list would keep it here reporting a
    single-arm number when two are available.
    """
    import yaml

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    has_descriptor = {k.lower() for k, v in mp.items() if v["descriptors"]}
    for r in d["rows"]:
        assert r["modality"] not in has_descriptor, (
            f"{r['modality']} now has a MeSH descriptor, so it can be "
            "cross-checked and does not belong in a single-arm table")


def test_the_missing_cross_check_is_stated_as_structural(d):
    """The asymmetry that makes these rows weaker than the ones beside them."""
    md = MD.read_text()
    assert "BY CONSTRUCTION" in md
    assert "no second arm" in md
    assert "LESS certain" in md, (
        "the report does not say these figures are less certain than the "
        "descriptor-backed ones, which is how a precise-looking count reads "
        "unless told otherwise")
    assert "82.5%" in md, (
        "the single instrument's measured precision is not stated")


def test_the_vintage_flag_is_derived_and_its_mechanism_named(d):
    """A modality whose literature predates its own existence is matching an
    older use of its words."""
    for r in d["rows"]:
        assert r["vintage_suspect"] == bool(
            r["pre_vintage_share"] and r["pre_vintage_share"] >= 10)
        if r["first_year"]:
            assert r["first_year"] <= (r["median_year"] or r["first_year"])
    assert sorted(d["suspect"]) == sorted(
        r["modality"] for r in d["rows"] if r["vintage_suspect"])
    md = MD.read_text()
    if d["suspect"]:
        assert "liechtenstein" in md.lower(), (
            "the substring collision is not shown, so the flag reads as a rule "
            "rather than an observed failure")
        assert "copper ionophore" in md
        assert "upper bounds" in md


def test_flagged_rows_are_published_not_filtered(d):
    """A date cut-off would be a second judgement layered on the first."""
    for m in d["suspect"]:
        row = next(r for r in d["rows"] if r["modality"] == m)
        assert row["articles"] > 0, "a flagged row was filtered out"
        assert f"| {m} *" in MD.read_text(), (
            f"{m} is flagged in the data but not marked in the table")
    md = MD.read_text()
    assert "rather than filtered" in md


def test_totals_recompute(d):
    assert d["total_articles"] == sum(r["articles"] for r in d["rows"])
    assert d["total_trials"] == sum(r["trials"] for r in d["rows"])
    for r in d["rows"]:
        if r["articles"]:
            assert r["trial_share"] == pytest.approx(
                100 * r["trials"] / r["articles"], abs=0.02)
        assert r["late_phase"] <= r["trials"], (
            f"{r['modality']} has more late-phase records than trials")


def test_the_manuscript_states_what_not_measurable_costs():
    """The phrase appears repeatedly; at least once it must carry a size.

    Used bare, it invites the reading it was meant to prevent -- that there is
    nothing there. Radioligand therapy alone has hundreds of indexed trials and
    no descriptor.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    d = json.loads(JSON.read_text())
    rl = next(r for r in d["rows"] if r["modality"] == "radioligand-therapy")
    assert f"{d['total_articles']:,} articles" in txt, (
        "the manuscript uses 'not measurable' without ever saying how much "
        "literature it refers to")
    assert f"{rl['trials']:,}" in txt or f"{rl['articles']:,}" in txt
    assert "not measurable by descriptor" in txt.lower(), (
        "the manuscript does not qualify the phrase, so it reads as a claim "
        "about the field rather than about MeSH")
