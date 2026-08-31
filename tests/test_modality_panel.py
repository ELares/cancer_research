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


def test_the_adoptive_barrier_section_exists_and_its_factors_multiply_out():
    """The PR's headline deliverable had ZERO coverage.

    No test referenced `adoptive_barriers` at all, and the renderer guards it
    with `if ab:` -- so deleting the whole block from the JSON made the "One
    construct, two diseases" section vanish silently with the full suite green,
    while every other section of the same report was guarded.

    The identity is the load-bearing part. An earlier version computed the
    antigen factor as `collapse / (delivery * persistence)`, a RESIDUAL that
    absorbed whatever the other two did not explain -- so "the factors multiply
    to the collapse" was true by construction and could not fail. All three are
    read from independently measured quantities now, which means the product
    CAN disagree with the collapse, and this is where that would surface.
    """
    d = json.loads(JSON_.read_text())
    ab = d.get("adoptive_barriers")
    assert ab, "the adoptive-barrier section is gone from the panel output"
    for k in ("leukaemia_kill_fraction", "solid_tumour_kill_fraction",
              "solid_tumour_kill_fraction_before_ceiling",
              "delivery_efficiency_solid", "persistence_at_run_end_solid",
              "antigen_ceiling_solid", "antigen_ceiling_binds"):
        assert k in ab, f"the panel no longer reports {k}"

    collapse = ab["leukaemia_kill_fraction"] / ab["solid_tumour_kill_fraction"]
    product = ((1.0 / ab["delivery_efficiency_solid"])
               * (1.0 / ab["persistence_at_run_end_solid"])
               * (ab["solid_tumour_kill_fraction_before_ceiling"]
                  / ab["solid_tumour_kill_fraction"]))
    assert abs(product / collapse - 1.0) < 0.01, (
        f"the three factors multiply to {product:,.0f}x against a {collapse:,.0f}x "
        "collapse, so the page's decomposition is incomplete")

    # The ceiling is a CAP, so its factor is 1 exactly when it does not fire.
    # Checking the flag against the arithmetic keeps the binary and the report
    # from disagreeing about which happened.
    fired = ab["antigen_ceiling_binds"]
    ceiling_x = (ab["solid_tumour_kill_fraction_before_ceiling"]
                 / ab["solid_tumour_kill_fraction"])
    assert fired == (ceiling_x > 1.000_001), (
        f"the binary says the ceiling binds={fired} while its own numbers give "
        f"a factor of {ceiling_x}")
    assert ab["leukaemia_kill_fraction"] > ab["solid_tumour_kill_fraction"], (
        "the solid tumour is no longer the harder case, which would invert the "
        "whole section")


def test_the_rendered_page_states_whichever_ceiling_verdict_is_true():
    """The verdict is prose about a computed quantity, so it must be computed.

    A reviewer set the antigen fraction to 1e-9, which made the ceiling the ONLY
    thing determining the kill, and the page went on saying it "does not bind"
    while the headline moved by three orders of magnitude.
    """
    d = json.loads(JSON_.read_text())
    ab = d["adoptive_barriers"]
    md = MD.read_text()
    assert "## One construct, two diseases" in md
    if ab["antigen_ceiling_binds"]:
        assert "**The antigen ceiling binds here**" in md
        assert "does not bind here" not in md
    else:
        assert "**The antigen ceiling does not bind here**" in md
    # and the collapse the prose quotes is the one the JSON supports
    collapse = ab["leukaemia_kill_fraction"] / ab["solid_tumour_kill_fraction"]
    assert f"{collapse:,.0f}-fold collapse" in md
