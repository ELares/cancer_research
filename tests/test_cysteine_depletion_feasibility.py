"""Guards for the bacterial cysteine-depletion feasibility scan.

PMID 39696900 delivers cyst(e)inase from engineered bacteria. The engine already
carries `erastin_xc_inhib`/`transsulfuration_floor`, so the only genuinely new
piece is the map from a residual CONCENTRATION to a supply FRACTION, anchored on
a measured transporter Km.

Three things here would be easy to break silently and each is tied to a
recomputed quantity rather than to prose:

1. The map must not degenerate to the identity. "Depletion equals inhibition" is
   the intuition the whole page refutes, and it is exactly what a well-meaning
   simplification of `relative_supply` would produce.
2. The marginal return CROSSES one. An earlier draft of the generator claimed it
   "only reaches 1 asymptotically" while the table it rendered on the same page
   showed 1.30, so the numeric property is guarded, not the sentence.
3. The engine-structure claims are re-read from the Rust source here, not from
   the generator's JSON. A guard that reads only a committed artifact cannot fail
   when the engine moves underneath it -- the defect
   `tests/test_generator_guards_pin_source.py` exists to block.
"""

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from cysteine_depletion_feasibility import (  # noqa: E402
    MissingInputError,
    depletion_for_inhibition,
    inhibition_from_depletion,
    marginal_crossover,
    marginal_return,
    read_structural_facts,
    relative_supply,
    saturation,
)

CAL = REPO_ROOT / "analysis" / "calibration"
RAW = CAL / "cysteine-depletion-feasibility.json"
DOC = CAL / "cysteine-depletion-feasibility.md"
CONSTS = CAL / "cysteine_depletion_constants.csv"
POOLS = CAL / "cysteine_depletion_pools.csv"
CORE = REPO_ROOT / "simulations" / "ferroptosis-core" / "src"


def _data():
    return json.loads(RAW.read_text())


def test_the_reference_concentration_is_derived_from_the_formulation_not_typed():
    """65.0 mg/L of a 313.0 g/mol salt is 207.7 uM. Recomputed from the CSV."""
    rows = {r["constant"]: r for r in csv.DictReader(CONSTS.open())}
    mass = float(rows["rpmi1640_cystine_2hcl_mass"]["value"])
    mw = float(rows["cystine_2hcl_molar_mass"]["value"])
    expected = mass / mw * 1000.0
    assert abs(_data()["reference_uM"] - expected) < 1e-9
    # Sanity on the magnitude, so a unit slip cannot pass the equality above.
    assert 150.0 < expected < 300.0


def test_the_map_is_not_the_identity_which_is_the_whole_finding():
    """Depletion must buy LESS than its share wherever substrate remains.

    If `relative_supply` were ever simplified to a linear ratio, inhibition would
    equal depletion and this page would say nothing. Checked across the range,
    not at one point.
    """
    km, s0 = 179.4, 207.7
    for d in [0.1, 0.25, 0.5, 0.75, 0.9]:
        inhib = inhibition_from_depletion(d, s0, km)
        assert inhib < d, f"depletion {d} must over-promise, got inhibition {inhib}"
    # And the gap is substantial, not a rounding artifact: half the substrate
    # removed is worth well under half the supply.
    assert inhibition_from_depletion(0.5, s0, km) < 0.4


def test_the_inverse_round_trips_through_the_closed_form():
    km, s0 = 179.4, 207.7
    for target in [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]:
        d = depletion_for_inhibition(target, s0, km)
        assert abs(inhibition_from_depletion(d, s0, km) - target) < 1e-9


def test_a_deep_supply_cut_needs_near_complete_depletion():
    """The headline: 90% supply cut costs ~95% depletion, and >90% at any Km here."""
    s0 = 207.7
    for km in [89.7, 152.7, 179.4, 206.1, 358.8]:
        d = depletion_for_inhibition(0.9, s0, km)
        assert d > 0.9, f"Km {km} gave depletion {d}"
    assert abs(depletion_for_inhibition(0.9, s0, 179.4) - 0.951) < 0.002


def test_the_marginal_return_crosses_one_rather_than_approaching_it():
    """The retracted claim was that it tops out at 1. It ends at (Km+S0)/Km."""
    km, s0 = 179.4, 207.7
    start = marginal_return(0.0, s0, km)
    end = marginal_return(1.0, s0, km)
    assert start < 1.0, "the first substrate removed must under-deliver"
    assert end > 1.0, "the last substrate removed must over-deliver"
    # Independently: start is Km/(Km+S0) and end is its reciprocal.
    assert abs(start - km / (km + s0)) < 1e-12
    assert abs(end - (km + s0) / km) < 1e-12
    cross = marginal_crossover(s0, km)
    assert 0.0 < cross < 1.0
    assert marginal_return(cross - 0.05, s0, km) < 1.0 < marginal_return(cross + 0.05, s0, km)


def test_saturation_is_zero_at_zero_and_monotone():
    km = 179.4
    assert saturation(0.0, km) == 0.0
    assert saturation(-5.0, km) == 0.0
    vals = [saturation(s, km) for s in [1, 10, 100, 1000, 10000]]
    assert vals == sorted(vals)
    assert vals[-1] < 1.0


def test_relative_supply_is_one_at_the_reference():
    assert abs(relative_supply(207.7, 207.7, 179.4) - 1.0) < 1e-12


def test_the_structural_claims_are_read_from_rust_not_from_the_json():
    """Re-parse the engine. This is the guard that can fail when the engine moves."""
    facts = read_structural_facts()
    # The claim the page leans on: the escape is global, so a hypoxia-graded
    # escape is not expressible. If someone adds a per-cell floor, this fails and
    # the page must be rewritten -- which is the point.
    assert facts["floor_is_params_field"] is True
    assert facts["floor_is_cell_field"] is False
    assert facts["escape_can_vary_per_cell"] is False
    assert facts["nrf2_is_cell_field"] is True
    # Every glutathione resynthesis site must be gated by the cystine factor, or
    # NRF2 could rescue a complete block and the topology argument collapses.
    assert facts["gsh_resynthesis_sites"] > 0
    assert facts["resynthesis_sites_gated_by_cystine"] == facts["gsh_resynthesis_sites"]
    assert facts["nrf2_cannot_rescue_complete_block"] is True
    # And the committed page must agree with the live source.
    assert _data()["structure"] == facts


def test_the_cystine_factor_formula_still_has_the_shape_the_map_assumes():
    """The map treats `xc_cystine_factor` as a supply multiplier in [0,1].

    Parsed from source: a change to the formula invalidates the whole page, so it
    must not be possible to change it and leave this green.
    """
    src = (CORE / "params.rs").read_text()
    m = re.search(r"pub fn xc_cystine_factor\(&self\)[^{]*\{(.*?)\n    \}", src, re.S)
    assert m, "xc_cystine_factor not found in params.rs"
    body = m.group(1)
    # PER FIELD, not a substring search over the whole body. A single
    # `"clamp(0.0, 1.0)" in body` check passed a mutation that unclamped one of
    # the two operands, because the other still carried a clamp -- the
    # fix-one-site-claim-all-sites shape.
    for field in ("erastin_xc_inhib", "transsulfuration_floor"):
        assert re.search(rf"self\.{field}\s*\.clamp\(0\.0,\s*1\.0\)", body), (
            f"{field} must be clamped into [0,1] for the concentration map to hold")


def test_the_page_reports_both_pool_numbers_so_neither_stands_alone():
    """The correction is that a pool of 1 measured a phrase pair."""
    d = _data()
    pool = {r["key"]: r["pool"] for r in d["pool_rows"]}
    assert pool["phrase_pair_ferro_engbact"] == 1
    assert pool["cystine_import_ferro"] > 1000
    doc = DOC.read_text()
    assert f"{pool['cystine_import_ferro']:,}" in doc
    assert "vocabulary" in doc.lower()
    # And it must NOT publish a multiplier over a pool of one.
    assert "3404x" not in doc and "3,404x" not in doc


def test_the_pool_table_keeps_the_csv_order_not_an_alphabetical_one():
    """`sort_keys=True` reordering a stored dict is on record in this repo.

    It made `--render-only` produce a different document from a fresh run. The
    artifact stores ordered pairs, so the two agree; this pins the ordering
    itself rather than only the agreement, because agreement would also hold if
    BOTH orders became alphabetical.
    """
    csv_order = [r["key"] for r in csv.DictReader(POOLS.open())]
    stored_order = [r["key"] for r in _data()["pool_rows"]]
    assert stored_order == csv_order
    assert stored_order != sorted(stored_order), (
        "the CSV order happens to be alphabetical, so this guard cannot "
        "distinguish it from a sort; reorder the CSV to keep the check live")


def test_no_layer_is_claimed_to_have_been_built():
    doc = DOC.read_text().lower()
    assert "no layer" in doc
    # The page must not promote the dose magnitude to anchored.
    assert "unanchored" in doc


def test_it_refuses_rather_than_rendering_from_absent_inputs(tmp_path):
    import cysteine_depletion_feasibility as m
    with pytest.raises(MissingInputError):
        m.load_constants(tmp_path / "nope.csv")
    with pytest.raises(MissingInputError):
        m.load_pools(tmp_path / "nope.csv")
    # A partial file is refused too, not silently half-rendered.
    partial = tmp_path / "partial.csv"
    partial.write_text("constant,value,unit,sd,source_pmid,source_doi,source_note\n"
                       "xct_cystine_km,179.4,uM,26.7,34287957,x,y\n")
    with pytest.raises(MissingInputError):
        m.load_constants(partial)


def test_every_committed_constant_carries_a_source_or_a_stated_reason():
    for r in csv.DictReader(CONSTS.open()):
        assert r["source_note"].strip(), f"{r['constant']} has no provenance note"
        assert r["unit"].strip(), f"{r['constant']} has no unit"


def test_render_only_reproduces_the_page_without_recomputing():
    before = DOC.read_text()
    r = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "cysteine_depletion_feasibility.py"),
         "--render-only"],
        capture_output=True, text=True, cwd=REPO_ROOT)
    assert r.returncode == 0, r.stderr
    assert DOC.read_text() == before


def test_every_landmark_paper_names_the_pool_it_was_actually_found_in():
    """An earlier draft asserted PMID 41704298 was among the `cyst(e)inase` records.

    It is not: it sits in the cysteine-depletion, bacteria and transsulfuration
    pools. Membership is committed data now, so the claim is checkable offline,
    and this pins the specific error so it cannot come back as prose.
    """
    landmarks = {lm["pmid"]: lm["pools"] for lm in _data()["landmarks"]}
    assert "cysteinase_any" not in landmarks["41704298"], (
        "41704298 is not in the cyst(e)inase pool; that was the retracted claim")
    # The papers the page does place there must really be there.
    assert "cysteinase_any" in landmarks["27869804"]
    assert "cysteinase_any" in landmarks["32241947"]
    # The target paper and the Km paper must each be reachable at all.
    assert landmarks["39696900"], "the target paper matched no pool"
    assert "xct_kinetics" in landmarks["34287957"]
    # Every pool named by a landmark must be a pool the page actually reports.
    known = {r["key"] for r in _data()["pool_rows"]}
    for pmid, pools in landmarks.items():
        assert pools, f"{pmid} matched no pool"
        unknown = set(pools) - known
        assert not unknown, f"{pmid} cites pools absent from the table: {unknown}"


def test_the_two_bacterial_papers_claim_is_derived_from_membership():
    """The page says the bacteria pool holds two papers. Recount it."""
    landmarks = _data()["landmarks"]
    bact = [lm["pmid"] for lm in landmarks if "bacteria_cysteine_ferro" in lm["pools"]]
    assert len(bact) == 2, f"expected two bacterial landmarks, got {bact}"
    doc = DOC.read_text()
    assert f"holds {len(bact)} bacterial papers" in doc
    for pmid in bact:
        assert pmid in doc


def test_a_landmark_with_no_pool_is_refused_rather_than_cited(tmp_path):
    import cysteine_depletion_feasibility as m
    bad = tmp_path / "landmarks.csv"
    bad.write_text("pmid,found_in_pools,note\n39696900,,the target paper\n")
    with pytest.raises(MissingInputError):
        m.load_landmarks(bad)
    empty = tmp_path / "empty.csv"
    empty.write_text("pmid,found_in_pools,note\n")
    with pytest.raises(MissingInputError):
        m.load_landmarks(empty)


def test_the_prose_quotes_the_committed_constants_rather_than_literals():
    """"Derived number, handwritten prose" is a recorded defect class here.

    The renderer once typed `65.0 mg/L` and `313.0 g/mol` beside a derived
    concentration, so editing the CSV would have moved the figure and left the
    sentence behind. Read from the CSV, so a change that is not regenerated fails.
    """
    rows = {r["constant"]: r for r in csv.DictReader(CONSTS.open())}
    doc = DOC.read_text()
    for key in ("rpmi1640_cystine_2hcl_mass", "cystine_2hcl_molar_mass"):
        assert rows[key]["value"] in doc, (
            f"{key} = {rows[key]['value']} is not quoted in the page; the prose "
            "has drifted from the committed constant")
    pmid = rows["xct_cystine_km"]["source_pmid"]
    assert f"PMID {pmid}" in doc, "the Km's source PMID is not quoted in the page"
    # And the Km itself, which the whole map turns on.
    assert rows["xct_cystine_km"]["value"] in doc
    assert rows["xct_cystine_km"]["sd"] in doc


def test_the_map_ships_its_bounds_not_just_its_arithmetic():
    """The three assumptions are where the map breaks before the numbers do.

    Guarded by the SUBSTANCE rather than the heading: the parallel-acquisition
    route is the bound that would most change a reader's use of the map, and a
    heading-only check would survive its deletion.
    """
    doc = DOC.read_text()
    # The parallel route an enzyme in the medium cannot reach.
    assert "35287706" in doc, "the macropinocytosis bound is not cited"
    assert "macropinocytosis" in doc.lower()
    # The proportionality step that lets a transporter rate stand in for supply.
    assert "proportional to cystine uptake" in doc
    # The initial-rate-vs-steady-state limit.
    assert "initial-rate" in doc.lower()


def test_the_headline_is_marked_in_vitro_only_because_the_gap_is_regime_specific():
    """The depletion-vs-supply gap exists only ABOVE the Km.

    Below it the saturation curve is nearly linear and depletion maps almost
    one-to-one, so presenting the culture figure as a general property would
    overclaim. Recomputed here rather than read, and the low reference must be
    flagged as not-free-cystine so it is never promoted to an in-vivo anchor.
    """
    d = _data()
    g = d["regime"]
    km = d["km_uM"]
    assert g["low_reference_uM"] < km, "the low reference must sit below the Km"
    assert d["reference_uM"] > km, "the medium reference must sit above the Km"
    # Independently: half-depletion buys more at the low reference than at medium.
    low = inhibition_from_depletion(0.5, g["low_reference_uM"], km)
    high = inhibition_from_depletion(0.5, d["reference_uM"], km)
    assert low > high, "the gap must shrink below the Km, or the qualification is wrong"
    assert abs(g["inhibition_at_half_depletion_low"] - low) < 1e-12
    assert abs(g["inhibition_at_half_depletion_reference"] - high) < 1e-12
    # The low reference must never be presented as free cystine.
    assert g["low_reference_is_free_cystine"] is False
    doc = DOC.read_text()
    assert "in-vitro only" in doc
    assert "not free cystine" in doc.lower() or "NOT free cystine" in doc
