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
    # This asserted `floor_is_cell_field is False` and `escape_can_vary_per_cell
    # is False`, and #916 added the per-cell floor, so it fired -- which is what
    # it was for: the page's "not expressible" section was obsolete and had to be
    # rewritten. It now pins the AGREEMENT between the page and the source plus
    # the invariants that must hold either way, rather than one side of a fact
    # that was always going to move. The page itself branches on
    # `escape_can_vary_per_cell`, so it reads correctly in both states.
    assert facts["floor_is_params_field"] is True
    assert facts["escape_can_vary_per_cell"] is facts["floor_is_cell_field"], (
        "the two must be the same question; a per-cell floor IS the per-cell escape")
    assert facts["nrf2_is_cell_field"] is True
    # Every glutathione resynthesis site must be gated by the cystine factor, or
    # NRF2 could rescue a complete block and the topology argument collapses.
    assert facts["gsh_resynthesis_sites"] > 0
    assert facts["resynthesis_sites_gated_by_cystine"] == facts["gsh_resynthesis_sites"]
    assert facts["nrf2_cannot_rescue_complete_block"] is True
    # And the committed page must agree with the live source.
    assert _data()["structure"] == facts


def test_the_spatial_consumer_claim_is_measured_not_asserted(tmp_path, monkeypatch):
    """The sentence this page used to publish was true when written and false the
    moment #918 landed, and nothing failed, because it was prose: "`erastin_xc_inhib`
    is reachable in no simulation binary at all, and sim-tme-3d applies no parameter
    overrides". Both halves are detectors now. A value-level assertion cannot tell a
    READ from a literal that happens to be right today, so this points them at
    doctored sources and requires the answers to MOVE.
    """
    import cysteine_depletion_feasibility as m

    live = m.read_structural_facts()
    assert live["arm_is_spatially_reachable"] is True
    assert live["floor_has_spatial_consumer"] is True

    # A binary doing neither must answer False for both, or the detectors are
    # being satisfied by something other than what they name.
    blind = tmp_path / "blind.rs"
    blind.write_text("fn main() {\n    let mut params = Params::default();\n}\n")
    monkeypatch.setattr(m, "TME3D_RS", blind)
    out = m.read_structural_facts()
    assert out["arm_is_spatially_reachable"] is False
    assert out["floor_has_spatial_consumer"] is False

    # A write that exists only inside the crate's test module is NOT a consumer.
    # This repo has already shipped a detector whose single match was an
    # `assert!` message in a `#[cfg(test)]` block.
    test_only = tmp_path / "test_only.rs"
    test_only.write_text(
        "fn main() {}\n\n#[cfg(test)]\nmod tests {\n"
        "    fn ab() {\n"
        "        apply_param_overrides(&mut params, ov);\n"
        "        gc.cell.transsulfuration_floor = Some(0.5);\n"
        "    }\n}\n"
    )
    monkeypatch.setattr(m, "TME3D_RS", test_only)
    out = m.read_structural_facts()
    assert out["floor_has_spatial_consumer"] is False, (
        "an A/B test setting a parameter is not a spatial consumer")
    assert out["arm_is_spatially_reachable"] is False


def test_the_page_states_the_reachability_it_measured():
    """Page and source cannot drift: the two fact lines are rendered from the
    detectors, and while the arm IS reachable the retracted sentence must be gone.
    """
    doc = DOC.read_text()
    facts = read_structural_facts()
    assert (f"the #502 arm is driveable in sim-tme-3d: "
            f"{facts['arm_is_spatially_reachable']}") in doc
    assert (f"something writes the per-cell floor in sim-tme-3d: "
            f"{facts['floor_has_spatial_consumer']}") in doc
    if facts["arm_is_spatially_reachable"]:
        assert "reachable in no simulation binary" not in doc, (
            "the arm is reachable; the page must not still say it is not")
        assert "applies no parameter overrides" not in doc


def test_the_cystine_factor_formula_still_has_the_shape_the_map_assumes():
    """The map treats `xc_cystine_factor` as a supply multiplier in [0,1].

    Parsed from source: a change to the formula invalidates the whole page, so it
    must not be possible to change it and leave this green.
    """
    src = (CORE / "params.rs").read_text()
    # `xc_cystine_factor_for`, because #916 made the no-argument form delegate to
    # it: anchoring on the old name found a one-line body and reported both
    # clamps missing. The arithmetic lives wherever the clamps are, so the guard
    # follows the arithmetic rather than the name it used to have.
    m = re.search(r"pub fn xc_cystine_factor_for\(&self[^{]*\{(.*?)\n    \}", src, re.S)
    assert m, "xc_cystine_factor_for not found in params.rs"
    body = m.group(1)
    # PER FIELD, not a substring search over the whole body. A single
    # `"clamp(0.0, 1.0)" in body` check passed a mutation that unclamped one of
    # the two operands, because the other still carried a clamp -- the
    # fix-one-site-claim-all-sites shape.
    assert re.search(r"self\.erastin_xc_inhib\s*\.clamp\(0\.0,\s*1\.0\)", body), (
        "erastin_xc_inhib must be clamped into [0,1] for the concentration map to hold")
    # The floor reaches the clamp through the per-cell fallback now, so pin the
    # fallback AND the clamp rather than one literal expression.
    assert "unwrap_or(self.transsulfuration_floor)" in body, (
        "the per-cell floor no longer falls back to the global field")
    assert re.search(r"\.clamp\(0\.0,\s*1\.0\)", body.split("unwrap_or")[1]), (
        "the resolved floor is not clamped into [0,1]")
    # And the no-argument form must still exist and delegate, or every existing
    # caller silently changed meaning.
    assert re.search(r"pub fn xc_cystine_factor\(&self\)[^{]*\{\s*"
                     r"self\.xc_cystine_factor_for\(None\)", src), (
        "xc_cystine_factor no longer delegates to the per-cell form")


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


def test_this_generator_hand_writes_no_quantity_into_its_prose():
    """Run the repo's hand-written-figure detector against this generator.

    `tests/test_no_handwritten_figures_in_generators.py` is ENFORCED on five
    files by its own design and says so, so it passing says nothing about this
    one. Enrolling this file there would mean loosening its `IDENT_TOKEN` to
    swallow bare eight-digit numbers, which would blind it for the five files it
    does protect -- tightening or widening a scan breaks the other side. So the
    detector is run here instead, and every surviving hit must be an IDENTIFIER
    this page genuinely cites rather than a quantity.

    EXACT, not a floor: an unrecognised hit fails, so a newly hand-written
    figure cannot hide behind the known false positives.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "_hwf", REPO_ROOT / "tests" / "test_no_handwritten_figures_in_generators.py")
    hwf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hwf)
    gen = REPO_ROOT / "scripts" / "cysteine_depletion_feasibility.py"

    # Identifiers this page cites: the PMIDs it names and the issues it builds on.
    # Not a style allowlist -- each is a reference that cannot "go stale" the way
    # a count can, because it denotes one fixed document. The issue numbers earn
    # the same treatment as the PMIDs for the same reason, and they need naming
    # here because the detector reads `#916 library` as a count while `(#916)`
    # reads as a reference -- a difference in the following word, not in kind.
    identifiers = {"39696900", "41704298", "27869804", "32241947",
                   "34287957", "41422256", "31819185", "35287706", "5054467",
                   "502", "916", "918", "920"}
    # Two non-numeric shapes the detector reads as counts and which are not:
    # "a tension" (the article, not a quantity) and the salt name "2HCl".
    non_quantities = {"tension", "2HCl"}

    unexplained = []
    for _lineno, _why, frag, _txt in hwf.offending(gen):
        token = frag.split()[0].rstrip(",.;:")
        if token in identifiers or token in non_quantities:
            continue
        if any(w in non_quantities for w in frag.split()):
            continue
        unexplained.append(frag)
    assert not unexplained, (
        "these look like hand-written quantities in the generator's prose, so "
        f"the sentence can outlive the measurement: {unexplained}")


def test_the_detector_import_actually_found_something_to_triage():
    """A triage over an empty list would pass for the wrong reason.

    If the detector stopped returning hits -- a broken import, a changed API --
    the test above would go green while checking nothing. This pins that it is
    still reading this file and still finding the identifier mentions.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "_hwf2", REPO_ROOT / "tests" / "test_no_handwritten_figures_in_generators.py")
    hwf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hwf)
    hits = hwf.offending(REPO_ROOT / "scripts" / "cysteine_depletion_feasibility.py")
    assert len(hits) >= 5, (
        f"the detector returned only {len(hits)} hits on a page citing eight "
        "PMIDs; it has probably stopped scanning this file")


def test_the_single_cell_line_anchor_is_shown_to_be_the_whole_literature():
    """"One paper" is a measurement here, not a confession of a shallow search.

    Of the records naming cystine uptake only two report kinetic constants, and
    one of those is a 1972 non-tumour study -- so in cancer cell lines there is
    exactly one. Recomputed from the committed pools, because this is what turns
    "measure it in a PDAC line" from tidying-up into a novel experiment.
    """
    pool = {r["key"]: r["pool"] for r in _data()["pool_rows"]}
    assert pool["cystine_uptake_kinetic"] == 2, (
        "the scarcity claim rests on this count; if it moved, the prose must too")
    # The claim is only interesting if the denominator is large.
    assert pool["cystine_uptake_any"] > 100
    assert pool["cystine_uptake_kinetic"] < pool["cystine_uptake_any"] / 50
    # And that no pancreatic kinetic measurement exists is the load-bearing half.
    assert pool["cystine_uptake_pancreatic"] > 0, (
        "a zero here would mean the query is broken, not that none exist")
    doc = DOC.read_text()
    assert str(pool["cystine_uptake_any"]) in doc
    assert "5054467" in doc, "the non-tumour kinetic paper must be named, not elided"


def test_the_map_is_stated_to_set_only_the_inhibition_not_the_floor():
    """The map equals the engine's supply factor exactly only at a zero floor.

    `xc_cystine_factor() = 1 - erastin_xc_inhib * (1 - transsulfuration_floor)`,
    so a reader who sets both from this page's figures would double-count the
    escape. Recomputed here so the arithmetic behind the warning is checked, not
    just its presence.
    """
    def factor(inhib, floor):
        return 1.0 - inhib * (1.0 - floor)

    inhib = inhibition_from_depletion(0.5, 207.7, 179.4)
    # At a zero floor the map's remaining supply IS the engine's factor.
    assert abs(factor(inhib, 0.0) - (1.0 - inhib)) < 1e-12
    # With a floor the engine retains MORE supply, so the page must warn.
    assert factor(inhib, 0.3) > 1.0 - inhib
    doc = DOC.read_text()
    assert "erastin_xc_inhib` alone" in doc
    assert "count the floor twice" in doc
