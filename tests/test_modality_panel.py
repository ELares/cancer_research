"""Guards for `analysis/modality-panel.{md,json}` and Figure 31.

The panel runs every treatment arm against the same tumour, which makes it the
first artifact in this repository where a non-ferroptosis arm produces a number
a reader can compare -- and the first that invites being read as a ranking of
therapies, which it is not.

So the guards are about three things: that the arms really are comparable
(same tumour, same seed), that the refusal is stated and cannot be quietly
dropped, and that the structural claim -- most arms no longer route through
the ferroptosis engine -- is measured rather than asserted.
"""
import ast
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/modality-panel.md"
JSON_ = REPO / "analysis/modality-panel.json"
BINARY = REPO / "simulations/sim-modality-panel/src/main.rs"
CELL_RS = REPO / "simulations/ferroptosis-core/src/cell.rs"
DIAGRAMS = REPO / "scripts/generate_conceptual_diagrams.py"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


@pytest.fixture(scope="module")
def md():
    return MD.read_text()


def test_every_treatment_variant_appears_in_the_panel(d):
    """The panel is only 'every applicable arm' if it is.

    Derived from the enum rather than from a list here, so a new variant that
    is not run fails instead of being quietly omitted -- which is exactly how
    an 'applicable' count drifts away from what a run can select.
    """
    src = CELL_RS.read_text()
    body = re.search(r"pub enum Treatment\s*\{(.*?)\n\}", src, re.S).group(1)
    variants = [v.strip().rstrip(",") for v in body.splitlines()
                if v.strip() and not v.strip().startswith(("//", "#["))]
    ran = {a["arm"] for a in d["arms"]}
    missing = set(variants) - ran
    assert not missing, (
        f"{sorted(missing)} are `Treatment` variants the panel does not run, "
        "so the report's 'every applicable arm' is false")
    assert ran <= set(variants), f"the panel runs arms that are not variants: "
    assert len(ran) == len(variants) == d["n_arms"]


def test_the_arms_are_actually_comparable(d):
    """One tumour, one seed. If the arms saw different tumours the table would
    be comparing noise."""
    assert d["seed"] is not None and d["n_cells"] > 1000
    src = BINARY.read_text()
    assert "a.seed.wrapping_add" in src, (
        "the arms no longer derive from a shared seed, so a difference "
        "between two rows is not necessarily the arm")
    # And the binary must SAY so, since the comparability is the premise.
    assert "SAME tumour" in src


def test_the_refusal_is_stated_in_every_place_it_can_be_read(d, md):
    """A table of kill fractions is read as a ranking unless it stops you.

    Three surfaces carry the refusal -- the JSON a machine reads, the report a
    person reads, and the figure a person skims -- because dropping it from
    any one of them leaves a route to the wrong reading.
    """
    assert "not_a_ranking" in d, "the JSON key that carries the refusal is gone"
    assert "not claims about clinical efficacy" in d["not_a_ranking"], (
        "the JSON's refusal no longer says what it refuses")
    assert "CALIBRATION_STATUS" in d["not_a_ranking"]
    assert "It is not a ranking of therapies" in md
    fig = DIAGRAMS.read_text()
    assert "NOT a ranking" in fig
    # Every row must carry its own calibration tier, so the caveat survives a
    # reader who only looks at the table.
    for a in d["arms"]:
        assert a["calibration"].strip(), a["arm"]
        assert a["route"].strip() and a["limited_by"].strip(), a["arm"]


def test_the_structural_claim_is_counted_not_asserted(d, md):
    """'Most arms no longer go through the ferroptosis engine' is the answer
    to the criticism this campaign began from, so it is a count."""
    ferro = [a for a in d["arms"] if "ferroptosis engine" in a["route"]]
    assert d["n_ferroptosis_routed"] == len(ferro)
    assert d["n_other_routes"] == d["n_arms"] - len(ferro)
    assert d["n_other_routes"] > d["n_ferroptosis_routed"], (
        f"only {d['n_other_routes']} of {d['n_arms']} arms route outside the "
        "ferroptosis engine; the report's structural claim is no longer true")
    assert len(d["distinct_routes"]) >= 5, d["distinct_routes"]
    assert f"{d['n_arms']} arms, {len(d['distinct_routes'])} distinct routes" in md


def test_the_delivery_finding_is_real_and_the_payloads_match(d, md):
    """The ADC row is the sharpest thing in the panel and rests on the two
    arms sharing a payload. If they ever stop sharing it, the comparison is
    between two different drugs and says nothing about delivery."""
    by = {a["arm"]: a for a in d["arms"]}
    adc, sdt = by["AntibodyDrugConjugate"], by["SDT"]
    assert adc["kill_fraction"] < sdt["kill_fraction"], (
        "the ADC no longer underperforms SDT; the delivery paragraph is stale")
    assert sdt["kill_fraction"] / max(adc["kill_fraction"], 1e-12) > 5.0, (
        "the delivery gap has shrunk below a factor of five and the report "
        "calls it the sharpest row")
    src = BINARY.read_text()
    assert "params.sdt_ros * avail" in src, (
        "the ADC arm no longer uses SDT's exogenous-ROS parameter, so the two "
        "rows do not share a payload and the comparison is between drugs")
    assert "binding-site barrier" in md


def test_the_sdt_pdt_coincidence_is_explained_rather_than_hidden(d, md):
    """Two identical rows in a comparison table look like a bug. They are a
    property of a depth-free panel, and the report has to say so."""
    by = {a["arm"]: a for a in d["arms"]}
    if abs(by["SDT"]["kill_fraction"] - by["PDT"]["kill_fraction"]) < 1e-9:
        assert "depth-free" in md.lower() or "DEPTH-FREE" in md
        assert "depth-reach-comparison.md" in md, (
            "the report notes the coincidence without pointing at the artifact "
            "where the difference actually lives")
    else:
        assert True  # they diverged; nothing to explain


def test_the_figure_types_no_number_and_groups_by_route():
    """Same rule as Figure 30: a diagram summarising a measurement is where a
    stale number hides."""
    tree = ast.parse(DIAGRAMS.read_text())
    fn = next((n for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef) and n.name == "fig31_modality_panel"), None)
    assert fn, "fig31_modality_panel is gone"
    src = ast.get_source_segment(DIAGRAMS.read_text(), fn)
    assert "modality-panel.json" in src
    big = [n.value for n in ast.walk(fn)
           if isinstance(n, ast.Constant) and isinstance(n.value, int)
           and not isinstance(n.value, bool) and n.value > 1_000]
    assert not big, f"typed integers on the figure: {big}"
    # The grouping IS the figure -- a plain bar chart would invite the ranking
    # reading the report refuses.
    assert "route_key" in src and "what kills" in src
