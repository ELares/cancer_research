"""Guards for the penetration access-vs-depth measurement (#926).

Every guard here drives `assemble()` on a DOCTORED copy of the committed raw
measurements rather than reading the stored derived fields. A guard that only
reads the artifact cannot fail when the generator changes, because the artifact
is committed and does not move until someone re-runs it -- the shape
`tests/test_generator_guards_pin_source.py` exists to catch.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import penetration_access_vs_depth as mod  # noqa: E402

ARTIFACT = REPO / "analysis" / "calibration" / "penetration-access-vs-depth.json"
pytestmark = pytest.mark.skipif(
    not ARTIFACT.exists(),
    reason="penetration-access-vs-depth.json not committed")


def _raw() -> dict:
    """The stored measurements with the derived block removed, so `assemble()`
    must rebuild every published field rather than echoing one."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    d.pop("derived", None)
    return d


def _derived(raw: dict) -> dict:
    return mod.assemble(raw)["derived"]


def test_the_wall_observable_does_not_respond_to_the_penetration_length():
    """The headline sits at r = 0, where exp(-r/lambda) is 1 for every lambda.

    Non-vacuity matters more than the equality here: two drugs with the SAME
    lambda would also agree, and that would say nothing. So the lambdas are
    asserted to differ first.
    """
    raw = _raw()
    lam = raw["penetration_length_um"]
    assert len(lam) >= 2, "need two drugs for this to be a comparison at all"
    assert len(set(lam.values())) == len(lam), (
        f"the drugs must differ in penetration length or the agreement below "
        f"is trivial: {lam}")
    v = _derived(raw)
    assert v["lambda_ratio"] > 1.5, (
        f"lambda ratio {v['lambda_ratio']} is too small for this to be a "
        "meaningful span")
    assert v["lambda_cancels_at_the_wall"], v["wall_identical_across_drugs"]


def test_a_wall_kill_that_moved_with_lambda_would_be_reported_as_such():
    """The cancellation flag is WIRED, not a constant: perturb one drug's wall
    kill and the verdict must stop claiming lambda cancels."""
    raw = _raw()
    drug = sorted(raw["wall_death_rate"])[0]
    tissue = sorted(raw["wall_death_rate"][drug])[0]
    raw["wall_death_rate"][drug][tissue] += 0.01
    v = _derived(raw)
    assert not v["lambda_cancels_at_the_wall"]
    assert v["verdict"] == "observable_depends_on_penetration_length"


def test_the_three_headline_numbers_are_one_dose_response():
    """All three tissues share the RNG stream at bin 0, so the headline triple is
    one function at three concentrations -- verified by reproducing each number
    exactly from the single-cell binding."""
    raw = _raw()
    c0 = raw["vessel_wall_concentration"]
    assert len(set(c0.values())) == len(c0), (
        f"the tissues must differ in vessel-wall concentration: {c0}")
    v = _derived(raw)
    assert v["headline_is_one_dose_response"], v["wall_reproduced_by_binding"]
    for tissue, r in v["wall_reproduced_by_binding"].items():
        assert r["binding"] == r["binary"], (
            f"{tissue}: binding {r['binding']} != binary {r['binary']}")


def test_the_ordering_claim_is_entailed_by_monotonicity():
    """The ordering holds in every draw BECAUSE the dose-response is monotone,
    which is a property of the kill switch and not of the tissue comparison."""
    raw = _raw()
    prior = raw["prior_ordering"]
    assert prior["ordering_preserved"] == prior["n_draws"], prior
    v = _derived(raw)
    assert v["dose_response_monotone"]
    assert v["ordering_entailed_not_measured"]
    assert v["verdict"] == "direction_robust_by_construction_magnitude_not"


def test_a_non_monotone_dose_response_would_stop_the_entailment_claim():
    """`ordering_entailed_not_measured` must track the measured monotonicity.

    Without this the flag could be hardcoded True and every sentence resting on
    it would be unfalsifiable."""
    raw = _raw()
    grid = raw["dose_response"]["grid"]
    hi = max(grid, key=lambda k: float(k))
    grid[hi] = 0.0  # a dip at the top end
    v = _derived(raw)
    assert not v["dose_response_monotone"]
    assert not v["ordering_entailed_not_measured"]
    assert v["verdict"] == "direction_not_entailed"


def test_the_cns_figure_is_mostly_its_own_no_drug_control():
    """Recomputed from the raw numbers, not read back: the drug-attributable
    part of the CNS headline must be a small fraction of the headline itself."""
    raw = _raw()
    v = _derived(raw)
    cns = v["access_order"][-1]
    a = v["drug_attributable"][cns]
    ctrl = raw["dose_response"]["control"]
    drug_only = raw["wall_death_rate"][sorted(raw["penetration_length_um"],
                                              key=lambda d: -raw["penetration_length_um"][d])[0]][cns] - ctrl
    assert a["control"] == ctrl
    assert a["drug_attributable"] == pytest.approx(drug_only, abs=1e-12)
    assert ctrl > 0, "a zero control would make this section vacuous"
    assert a["drug_attributable"] < a["control"], (
        "the claim is that most of the CNS number is the control; if the drug "
        f"contributes more ({a['drug_attributable']}) than the control "
        f"({a['control']}) the prose must be rewritten")
    # and the well-vascularized end must NOT have that property, or the
    # statement is about the model's baseline rather than about the gradient.
    well = v["drug_attributable"][v["access_order"][0]]
    assert well["drug_attributable"] > well["control"]


def test_access_and_distance_rank_the_tissues_differently():
    """The published ordering is the ACCESS one. If the two axes agreed, nothing
    would turn on which the headline measures."""
    raw = _raw()
    v = _derived(raw)
    assert not v["access_and_distance_agree"]
    poorly = [t for t in v["access_order"] if "poorly" in t.lower()][0]
    cns = [t for t in v["access_order"] if "cns" in t.lower()][0]
    half = raw["half_distance_um"]
    c0 = raw["vessel_wall_concentration"]
    assert c0[poorly] > c0[cns], "access ranks poorly above CNS"
    assert half[poorly] > half[cns], (
        "distance ranks CNS above poorly -- that disagreement is the finding")


def test_the_deep_crossover_is_checked_against_the_engine_not_asserted():
    """The crossover is derived from the stored presets, so on its own it would
    be a guard computing its own expectation. It is validated by requiring the
    two shipped lambdas to straddle it AND the measured deepest-cell
    concentrations to invert on exactly the predicted side."""
    raw = _raw()
    v = _derived(raw)
    assert v["crossover_bracketed_by_shipped_lambdas"], (
        f"crossover {v['deep_crossover_lambda_um']} not confirmed by the "
        f"engine: {v['deep_ordering_inverts']} at {raw['penetration_length_um']}")
    assert len(set(v["deep_ordering_inverts"].values())) == 2, (
        "both arms landing on the same side makes the bracket vacuous")
    lo, hi = mod.TANNOCK_BAND_UM
    assert v["crossover_inside_tannock_band"] == (
        lo <= v["deep_crossover_lambda_um"] <= hi)


def test_a_crossover_that_did_not_predict_the_inversion_is_refused():
    """Flip one measured deepest-cell concentration so the engine disagrees with
    the derived crossover; the confirmation flag must drop."""
    raw = _raw()
    v0 = _derived(copy.deepcopy(raw))
    assert v0["crossover_bracketed_by_shipped_lambdas"]
    poorly = [t for t in v0["access_order"] if "poorly" in t.lower()][0]
    cns = [t for t in v0["access_order"] if "cns" in t.lower()][0]
    for drug in raw["deepest_bin"]:
        raw["deepest_bin"][drug][cns]["concentration"] = (
            raw["deepest_bin"][drug][poorly]["concentration"] + 1.0)
    v = _derived(raw)
    assert not v["crossover_bracketed_by_shipped_lambdas"]


def test_render_only_reassembles_rather_than_echoing():
    """A derived field that is echoed on re-render makes every guard on it
    inert, which is how a sibling generator shipped a stale table."""
    stored = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert mod.assemble(stored) == stored, (
        "assemble() is not idempotent on its own output, so --render-only "
        "cannot reproduce a fresh run")
    assert mod.assemble(_raw())["derived"] == stored["derived"]


def test_the_report_states_the_claim_it_does_not_make():
    """The ordering claim is robust-by-construction, NOT false. A page that
    dropped that distinction would be overstating the finding."""
    raw = _raw()
    md = mod.render(mod.assemble(raw))
    assert "by construction" in md
    assert "not false" in md.lower()
    # and it must not dress the concentration inversion as a kill inversion
    assert "inversion is in the concentration, not the kill" in md


SCRIPT = REPO / "scripts" / "penetration_access_vs_depth.py"
TRANSPORT_RS = REPO / "simulations" / "ferroptosis-core" / "src" / "drug_transport.rs"
TISSUE_PK_RS = REPO / "simulations" / "sim-tissue-pk" / "src" / "main.rs"


def _src() -> str:
    """The generator's own text. Guards that read only the committed artifact
    cannot fail when the generator changes, because the artifact does not move
    until someone re-runs it."""
    return SCRIPT.read_text(encoding="utf-8")


def test_the_generator_derives_its_figures_rather_than_stating_them():
    """A committed artifact does not move when its generator is edited, so the
    guards above could all pass against a renderer that had been replaced by
    literals. Pin the source: none of the published figures may appear in it."""
    src = _src()
    body = src.split('"""', 2)[-1]  # the module docstring may quote figures
    for literal in ("0.121", "0.026", "0.018", "50.98", "0.0658", "14.7"):
        assert literal not in body, (
            f"{literal!r} is hardcoded in the generator body; every published "
            "figure must come from the measurement")
    # and the observable must stay the one the published interval used
    assert "extract_tissue_pk_observables" in src, (
        "the headline observable must be read through the repo's own extractor "
        "or it can drift away from the interval it qualifies")
    # the render-only path must RE-ASSEMBLE, not echo
    assert "assemble(json.loads(OUT_JSON.read_text()))" in src


def test_the_engine_still_has_the_structure_this_page_describes():
    """Both central claims are structural facts about Rust code, so they are
    re-parsed here rather than asserted in prose. If the engine changes, this
    page must be rewritten instead of going quietly stale."""
    transport = TRANSPORT_RS.read_text(encoding="utf-8")
    # The wall value is lambda-free: c0 is built from the permeability alone and
    # lambda enters only through exp(-r/lambda), which is 1 at r = 0.
    assert "let c0 = drug.vessel_wall_conc * tissue.vascular_permeability;" in transport, (
        "`concentration_at_distance` no longer builds the vessel-wall value "
        "from permeability alone, so section 1 of the page is out of date")
    assert "(c0 * (-r_um / lambda).exp()).min(1.0)" in transport, (
        "the exponential no longer carries r_um/lambda, so the claim that "
        "lambda cancels at r = 0 must be re-derived")

    pk = TISSUE_PK_RS.read_text(encoding="utf-8")
    seed_lines = [ln for ln in pk.splitlines() if "let cell_seed" in ln]
    assert len(seed_lines) == 1, seed_lines
    assert "tissue" not in seed_lines[0], (
        "`cell_seed` now depends on the tissue, so the three tissues no longer "
        "share an RNG stream and the entailment argument does not hold")
    for token in ("bin_idx", "seed"):
        assert token in seed_lines[0], seed_lines[0]
