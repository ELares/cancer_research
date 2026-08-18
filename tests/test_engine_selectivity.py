"""Guards for the modality selectivity contrast (#728).

WHAT THIS REPORTS
-----------------
Kill rate for every phenotype under every treatment the engine models, and the
tumour-versus-CAF contrast the engine has always been able to produce and has
never published. Two properties fall out of it that nothing else in the repo
records:

  SDT AND PDT ARE BIT-IDENTICAL in the single-cell path, because `sdt_ros` and
  `pdt_ros` share a default and nothing else differs. They are distinguished
  only by depth physics, in another module.

  RSL3 KILLS EXACTLY ZERO NON-TUMOUR CELLS. The project's load-bearing
  selectivity assumption is encoded in the `Stromal` parameters, so a ratio here
  restates the assumption instead of testing it.

THE ERROR THIS MUST NOT MAKE
-----------------------------
Calling any of it a therapeutic index. `Stromal` models cancer-associated
fibroblasts: tumour-resident, recruited by the tumour, parameterised to model
shielding. Not healthy tissue. An earlier version of #728 proposed converting it
and that is a category error -- so the guard checks the report keeps saying so,
in both the artifact and the renderer.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "engine_selectivity.py"
MD = REPO_ROOT / "analysis" / "engine-selectivity.md"
JSON_OUT = REPO_ROOT / "analysis" / "engine-selectivity.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def test_it_is_not_called_a_therapeutic_index():
    """Stromal is CAFs; dividing by it is not selectivity against normal tissue."""
    md, src = MD.read_text(), SCRIPT.read_text()
    assert "not a therapeutic index" in md.lower(), (
        "the report no longer disclaims being a therapeutic index")
    assert "cancer-associated fibroblasts" in md and \
        "cancer-associated fibroblasts" in src.lower().replace("_", " ") or \
        "CANCER-ASSOCIATED FIBROBLASTS" in src, (
        "the report no longer says what Stromal actually models")
    # In the RENDERED section, not merely somewhere in the file: the phrase
    # also appears in the module docstring, so a file-wide check passes even
    # when the sentence is deleted from the report the reader sees.
    assert "category error" in md, (
        "the REPORT no longer records that converting CAFs into a "
        "normal-tissue proxy was the withdrawn proposal")
    renderer = src[src.index("def render("):]
    assert "category error" in renderer, (
        "the renderer no longer emits the category-error note; the docstring "
        "mentioning it is not what a reader of the report sees")


def test_a_control_baseline_is_present_and_low():
    """A contrast against a population that dies on its own is not about treatment."""
    d = _doc()
    assert "Control" in d["treatments"], "no untreated baseline was run"
    c = d["treatments"]["Control"]["_contrast"]
    assert c["tumour_max"] < 0.5, (
        f"untreated tumour death is {100*c['tumour_max']:.1f}%; every contrast "
        "would be measuring spontaneous death")
    src = SCRIPT.read_text()
    assert 'if ctrl["tumour_max"] > 0.5:' in src, (
        "the baseline sanity check is gone from the generator")


def test_the_contrast_is_reported_both_ways():
    """One flattering phenotype should not carry a selectivity claim."""
    d, md = _doc(), MD.read_text()
    for t, row in d["treatments"].items():
        c = row["_contrast"]
        assert c["tumour_min"] <= c["tumour_max"]
        if c["ratio_best_case"] and c["ratio_worst_case"]:
            assert c["ratio_worst_case"] <= c["ratio_best_case"], (
                f"{t}: worst case exceeds best case")
    assert "best-case tumour:CAF" in md and "worst-case tumour:CAF" in md, (
        "only one side of the contrast is reported, so a modality can be "
        "flattered by its most-killed phenotype")


def test_the_sdt_pdt_identity_is_detected_not_assumed():
    """It is derived from the numbers, and must be reported while it holds."""
    d, md = _doc(), MD.read_text()
    sdt, pdt = d["treatments"]["SDT"], d["treatments"]["PDT"]
    same = all(abs(sdt[p]["death_rate"] - pdt[p]["death_rate"]) < 1e-12
               for p in sdt if not p.startswith("_"))
    if same:
        assert "same modality here" in md, (
            "SDT and PDT are bit-identical and the report does not say so, so a "
            "reader would take a single-cell contrast between them as evidence "
            "about two different therapies")
    else:
        assert "same modality here" not in md, (
            "the report claims SDT and PDT are identical when they no longer "
            "are; the defaults must have diverged")
    src = SCRIPT.read_text()
    assert "identical = all(" in src, (
        "the identity is no longer computed from the death rates; asserting it "
        "in prose would make it a claim rather than an observation")


def test_the_zero_denominator_is_flagged_as_a_tell():
    """An exactly-zero non-tumour kill is what an assumed answer looks like."""
    d, md = _doc(), MD.read_text()
    rsl3 = d["treatments"]["RSL3"]["Stromal"]["death_rate"]
    if rsl3 == 0.0:
        # THIS GUARD PINNED THE RETRACTED READING. It required the phrase
        # "true by construction", i.e. the fingerprint inference a control in
        # the same table refutes. What must be flagged is the UNDEFINED ratio
        # and the assumption-restating problem, neither of which needs it.
        assert "zero denominator" in md, (
            "RSL3 kills exactly zero non-tumour cells and the report does not "
            "flag the undefined ratio")
        assert "restate the assumption rather than test it" in md, (
            "the report no longer says a ratio against these parameters "
            "restates the selectivity assumption instead of testing it")
        assert d["treatments"]["RSL3"]["_contrast"]["ratio_best_case"] is None, (
            "a ratio was computed against a zero denominator")
    # and the CODE must still refuse it. The artifact is static, so a generator
    # that started dividing by zero would leave the assertion above passing.
    # BOTH ratios, counted. A substring check passes while one of the two is
    # changed to divide by zero, because the other still carries the phrase.
    n_guarded = SCRIPT.read_text().count("if caf > 0 else None")
    assert n_guarded == 2, (
        f"{n_guarded} of the 2 ratios guard against a zero denominator; a "
        "ratio against zero non-tumour deaths would be published as a number")
    src = SCRIPT.read_text()
    assert 'if rsl3_caf == 0.0:' in src, (
        "the zero-denominator branch is gone from the renderer")


def test_the_missing_normal_tissue_phenotype_stays_named():
    """The gap is the point; without it this is just a table."""
    md = MD.read_text()
    assert "ACSL4-low" in md, (
        "the report no longer says what a real normal-tissue parameter set "
        "would need, so the absence stops being actionable")
    assert "cannot disagree with the project about selectivity" in md or \
        "cannot currently check" in md, (
        "the report no longer states that the assumption is uncheckable inside "
        "the model that assumes it")


def test_the_uncalibrated_status_is_carried():
    md = MD.read_text()
    assert "Uncalibrated" in md and "CALIBRATION_STATUS" in md, (
        "these are default parameters and the report must say so, or the "
        "numbers read as measurements of biology")


def test_the_zero_is_not_read_as_a_fingerprint():
    """A control in the same row refutes the causal reading.

    The page said an exactly-zero non-tumour kill "is what a parameter set
    chosen to produce it looks like". RSL3 also kills exactly zero of a TUMOUR
    phenotype, printed in the same row of the same table, so the exactness is
    a property of the RSL3 path against a resistant parameterisation rather
    than evidence about how `Stromal` was chosen.
    """
    d, md = _doc(), MD.read_text()
    rsl3 = d["treatments"]["RSL3"]
    tumour_zeros = sorted(ph for ph, row in rsl3.items()
                          if not ph.startswith("_") and ph != "Stromal"
                          and row["death_rate"] == 0.0)
    if rsl3["Stromal"]["death_rate"] != 0.0:
        return
    assert tumour_zeros, (
        "no tumour phenotype has an exact zero any more, so the control that "
        "refutes the fingerprint reading is gone; re-check the argument")
    for ph in tumour_zeros:
        assert f"`{ph}`" in md, (
            f"`{ph}` is a tumour phenotype RSL3 also kills exactly zero of, "
            "and the report does not name it as the control")
    assert "That inference is withdrawn" in md, (
        "the fingerprint inference is stated without its withdrawal")
    # the CLAIMING form must not return
    for m in re.finditer(r"chosen to produce it looks like", md):
        w = md[max(0, m.start() - 300):m.end() + 300]
        assert re.search(r"earlier version|withdrawn|refutes", w, re.I), (
            "the report reads the exact zero as a fingerprint again, with a "
            "tumour phenotype at exactly zero in the same table")
    # and what survives must still be stated
    assert "restate the assumption rather than test it" in md


def test_the_seed_spread_uses_genuinely_disjoint_seeds():
    """`sim_batch` draws cell i from seed+2i, so runs must be spaced by 2n.

    Nothing said so, and the published ratios carried no interval -- anyone
    re-running at a neighbouring seed would have confirmed them by
    construction.
    """
    d, md = _doc(), MD.read_text()
    n = d["n_per_cell_type"]
    assert d.get("seed_stride") == 2 * n, (
        f"the seed stride is {d.get('seed_stride')} for n={n:,}; runs spaced "
        f"by less than {2*n:,} share cells and are not independent")
    sp = d.get("seed_spread") or {}
    assert sp, "the multi-seed spread is gone; the ratios carry no interval"
    for t, s in sp.items():
        assert s["n_seeds"] >= 4, f"{t}: only {s['n_seeds']} seeds"
        assert s["best_lo"] <= s["best_hi"] and s["worst_lo"] <= s["worst_hi"]
        assert s["worst_lo"] <= s["best_hi"], f"{t}: worst exceeds best"
        assert f"{s['best_lo']:.2f}x - {s['best_hi']:.2f}x" in md, (
            f"{t}'s measured best-case range is not rendered")
    assert f"spaced by {2*n:,}" in md
    # the SOURCE property, verified against the Rust rather than asserted
    lib = (REPO_ROOT / "simulations" / "ferroptosis-python" / "src" / "lib.rs")
    src = lib.read_text()
    assert "seed.wrapping_add((i as u64) * 2)" in src, (
        "sim_batch no longer derives its per-cell seed as seed + 2i, so the "
        "stride this report states is no longer the right one")


def test_the_point_estimates_are_not_presented_alone():
    """A single-seed figure quoted to a decimal, with a measured range beside it."""
    d, md = _doc(), MD.read_text()
    sp = d.get("seed_spread") or {}
    for t, s in sp.items():
        pt = d["treatments"][t]["_contrast"].get("ratio_best_case")
        if pt is None:
            continue
        assert "single-seed point estimates" in md, (
            "the point ratios are shown without saying they are single-seed")
        # and the range must actually contain, or bracket, the point
        assert s["best_lo"] - 1e-6 <= pt <= s["best_hi"] + 1e-6, (
            f"{t}: the published point ratio {pt:.2f} lies outside the "
            f"measured range [{s['best_lo']:.2f}, {s['best_hi']:.2f}], so one "
            "of the two was not computed from the same model")
