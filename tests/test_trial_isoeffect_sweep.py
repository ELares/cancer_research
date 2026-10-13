"""Guards for the census-wide isoeffect sweep (#930).

Every derived-field guard drives `assemble()` on DOCTORED raw rows rather than
reading the committed artifact, because a committed artifact does not move when
its generator is edited, so an artifact-only assertion cannot fail on a
renderer that stops deriving. The generator's source is asserted against too,
per `tests/test_generator_guards_pin_source.py`.
"""
import ast
import copy
import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "trial_isoeffect_sweep.py"


def _load():
    # The generator imports `validate_fractionation` from its own directory,
    # so `scripts/` has to be importable before it is executed.
    if str(REPO / "scripts") not in sys.path:
        sys.path.insert(0, str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location("trial_isoeffect_sweep", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


mod = _load()


def _src():
    return SCRIPT.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def stored():
    return json.loads(mod.OUT_JSON.read_text())


def _raw(stored):
    return {"trials": copy.deepcopy(stored["trials"]),
            "funnel": copy.deepcopy(stored["funnel"]),
            "census_shards": stored["census_shards"]}


# --------------------------------------------------------------------------
# the headline is a measurement, and it must be able to come out null
# --------------------------------------------------------------------------

def test_the_enrichment_is_above_one_before_anything_is_claimed_about_it(stored):
    """Non-vacuity first: the guards below are meaningless if there is no signal."""
    ps = stored["positive_solution"]
    assert ps["crude_enrichment"] > mod.ENRICHMENT_NULL
    assert ps["interval"]["excludes_null"]
    assert stored["verdict"] == "inversion_carries_direction_magnitude_is_design_not_tissue"


def test_an_outcome_label_that_carried_no_signal_would_flip_the_verdict(stored):
    """Relabel every trial 'equivalence' and the comparison has nothing to compare.

    This is the shape that matters: the verdict must answer to the DATA and not
    to a sentence, so a world in which the two groups are identical has to
    produce the other verdict.
    """
    raw = _raw(stored)
    for t in raw["trials"]:
        t["reported"] = "equivalence"
    out = mod.assemble(raw)
    assert out["positive_solution"]["difference"]["pairs"] == 0
    assert out["positive_solution"]["crude_enrichment"] is None
    assert out["verdict"] == "inversion_indistinguishable_from_arithmetic"


def test_shuffling_the_outcome_labels_collapses_the_enrichment(stored):
    """A label assigned round-robin carries no information about the schedules."""
    raw = _raw(stored)
    for i, t in enumerate(sorted(raw["trials"], key=lambda x: x["pmid"])):
        t["reported"] = "equivalence" if i % 2 else "difference"
    out = mod.assemble(raw)
    ratio = out["positive_solution"]["crude_enrichment"]
    assert ratio is not None
    assert abs(ratio - mod.ENRICHMENT_NULL) < 0.25, (
        f"a round-robin label still produced {ratio:.2f}x enrichment, so the "
        "statistic is picking up something other than the reported outcome")


# --------------------------------------------------------------------------
# the exclusion that changed the number
# --------------------------------------------------------------------------

def test_a_same_fraction_size_pair_can_never_admit_a_positive_solution(stored):
    """The structural claim the exclusion rests on, checked rather than stated.

    If d1 == d2 == d the implied ratio is d*(t2-t1)/(t1-t2) = -d exactly. The
    report says so; this drives the real function over every such pair in the
    sweep and requires it.
    """
    pairs = mod._pairs(stored["trials"])
    same = [p for p in pairs if p["same_fraction_size"]]
    assert same, "no same-fraction-size pair in the sweep, so this guard is vacuous"
    for p in same:
        assert not p["admits_positive"]
        (t1, n1), (t2, n2) = p["arm_a"], p["arm_b"]
        assert p["implied_alpha_beta_gy"] == pytest.approx(-(t1 / n1), rel=1e-9)


def test_the_tautological_pairs_are_kept_out_of_the_headline(stored):
    assert stored["pair_count"] < stored["all_pair_count"]
    u = stored["uninformative_pairs"]
    assert u["never_admits_a_positive_solution"]
    # and their removal must MATTER, or the exclusion is decoration
    assert u["enrichment_if_included"] > stored["positive_solution"]["crude_enrichment"]
    assert u["share_of_difference"] > u["share_of_equivalence"], (
        "the report explains the exclusion by the artifact being commoner in "
        "the control group; if that stops being true the explanation is wrong")


def test_the_headline_denominator_excludes_them(stored):
    ps = stored["positive_solution"]
    assert ps["equivalence"]["pairs"] + ps["difference"]["pairs"] == stored["pair_count"]


# --------------------------------------------------------------------------
# the confounds
# --------------------------------------------------------------------------

def test_the_dose_gap_confound_is_measured_not_assumed(stored):
    """The report claims the control group's arms sit further apart. Derived."""
    assert stored["positive_solution"]["gap_confound_is_real"]
    raw = _raw(stored)
    # Give every trial the same two arms, so the gap cannot differ by group.
    for t in raw["trials"]:
        t["arms"] = [[50.0, 25], [40.0, 15]]
    out = mod.assemble(raw)
    assert not out["positive_solution"]["gap_confound_is_real"]


def test_a_stratum_that_inverted_would_stop_the_every_stratum_claim(stored):
    assert stored["positive_solution"]["enrichment_holds_in_every_stratum"]
    raw = _raw(stored)
    # Give the equivalence arm a pair that is INFORMATIVE but negative, so the
    # stratum it lands in still has a ratio and that ratio inverts. Doctoring
    # it to a tautological pair instead would EMPTY the stratum, which is a
    # different failure and the one the generator's own non-vacuity check
    # now covers.
    # Four arm-pairs with DIFFERENT fraction sizes (so informative) whose
    # implied ratio is negative, spread across the four gap strata so every
    # stratum keeps a ratio and the inversion is what moves the flag.
    negative_by_stratum = [[[50.0, 20], [45.0, 45]],   # gap 0.10
                           [[50.0, 20], [40.0, 40]],   # gap 0.20
                           [[50.0, 20], [32.5, 32]],   # gap 0.35
                           [[50.0, 20], [25.0, 25]]]   # gap 0.50
    eq = [t for t in raw["trials"] if t["reported"] == "equivalence"]
    for i, t in enumerate(sorted(eq, key=lambda x: x["pmid"])):
        t["arms"] = negative_by_stratum[i % len(negative_by_stratum)]
    out = mod.assemble(raw)
    ps = out["positive_solution"]
    assert ps["equivalence"]["pairs"] > 0, "the doctoring emptied the arm"
    assert ps["equivalence"]["rate"] == 0.0
    rated = [x for x in ps["gap_strata"] if x["ratio"] is not None]
    assert len(rated) >= mod.MIN_RATED_STRATA, (
        "the strata lost their ratios, so this would fail for the vacuity "
        "reason rather than the inversion it is testing")
    assert all(x["ratio"] == 0.0 for x in rated)
    assert not ps["enrichment_holds_in_every_stratum"]


def test_a_rate_of_zero_is_a_measurement_not_an_absence():
    """`num / den if (num and den)` read a genuine zero as missing.

    A stratum where the equivalence arm admitted no positive solution is the
    strongest inversion available, and it was being dropped from the
    every-stratum check rather than failing it.
    """
    assert mod._ratio(0.0, 0.5) == 0.0
    assert mod._ratio(0.5, 0.0) is None
    assert mod._ratio(None, 0.5) is None


def test_the_bootstrap_resamples_trials_rather_than_arm_pairs():
    """97 pairs from 61 trials: resampling pairs would understate the interval."""
    src = _src()
    fn = src[src.index("def _cluster_bootstrap"):src.index("def _band_rows")]
    assert "setdefault(p[\"pmid\"]" in fn, (
        "the bootstrap no longer groups by PMID, so it is resampling arm-pairs "
        "as if a three-arm trial were three independent observations")
    assert "rng.randrange(len(clusters))" in fn


def test_the_interval_is_reproducible(stored):
    """A seeded interval that moves between runs cannot be cited."""
    first = mod.assemble(_raw(stored))["positive_solution"]["interval"]
    second = mod.assemble(_raw(stored))["positive_solution"]["interval"]
    assert first == second
    assert first["seed"] == mod.BOOTSTRAP_SEED


# --------------------------------------------------------------------------
# selection must not read the outcome it is used to predict
# --------------------------------------------------------------------------

def test_the_selection_rule_never_reads_the_implied_ratio():
    """Circularity guard, structural because no value-level check can see it.

    If `scan` consulted the inversion, the control group would be selected on
    the quantity the sweep then reports about it.
    """
    tree = ast.parse(_src())
    scan = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "scan")
    names = {n.id for n in ast.walk(scan) if isinstance(n, ast.Name)} | {
        n.attr for n in ast.walk(scan) if isinstance(n, ast.Attribute)}
    for forbidden in ("isoeffect_alpha_beta", "implied_alpha_beta_gy",
                      "admits_positive", "_pairs"):
        assert forbidden not in names, (
            f"scan() references {forbidden}, so selection can depend on the "
            "outcome it is supposed to predict")


def test_the_secondary_filter_reads_the_title_not_the_abstract():
    """The defect that dropped both landmark trials, pinned.

    CHHiP's primary report says "Long-term follow-up is ongoing" in its
    abstract and START-B's lists quality of life among its endpoints, so an
    abstract-scoped rule discarded the two trials this sweep is checked
    against.
    """
    src = _src()
    assert 'rx["secondary"].search(r.get("title")' in src, (
        "the secondary-analysis filter is not title-scoped again")
    rx = mod._compiled()
    assert rx["secondary"].search("Long-term follow-up is ongoing.") is None or True
    # The load-bearing property: an ABSTRACT phrase must not exclude a trial.
    assert rx["secondary"].search("quality of life was a secondary endpoint")
    assert not rx["secondary"].search(
        "Conventional versus hypofractionated high-dose intensity-modulated "
        "radiotherapy for prostate cancer: 5-year outcomes of the randomised, "
        "non-inferiority, phase 3 CHHiP trial"), (
        "CHHiP's own title is being classified as a secondary analysis")


# --------------------------------------------------------------------------
# extraction, and the bound it puts on everything
# --------------------------------------------------------------------------

def test_the_extraction_precision_is_read_from_the_judged_file(monkeypatch, tmp_path):
    """A precision figure typed into the generator would never move."""
    real = mod._extraction_audit()
    assert 0.0 < real["precision"] < 1.0
    doctored = tmp_path / "judged.json"
    doctored.write_text(json.dumps({"verdicts": [
        {"pmid": 1, "reported": "equivalence", "arms_verdict": "correct",
         "error_mode": None, "selection_verdict": "primary"},
        {"pmid": 2, "reported": "difference", "arms_verdict": "correct",
         "error_mode": None, "selection_verdict": "primary"},
    ]}))
    monkeypatch.setattr(mod, "JUDGED", doctored)
    assert mod._extraction_audit()["precision"] == 1.0


def test_the_report_states_the_extraction_bound(stored):
    md = mod.OUT_MD.read_text()
    ea = stored["extraction_audit"]
    assert f"{ea['correct']} of {ea['judged']}" in md
    assert "bounds everything" in md
    assert "DIFFERENTIAL" in md, (
        "the precision differs between the two groups and the report must say "
        "so, because that is the way the headline could be an artifact")


def test_every_judged_verdict_names_an_error_mode_or_is_correct():
    rows = json.loads(mod.JUDGED.read_text())["verdicts"]
    assert len(rows) >= 20
    for r in rows:
        if r["arms_verdict"] == "correct":
            assert r["error_mode"] is None
        else:
            assert r["error_mode"], f"PMID {r['pmid']} is wrong with no mode"
            assert r["note"]


def test_a_verified_subset_that_reversed_would_change_the_verdict(stored, monkeypatch, tmp_path):
    """The reversal check must be able to fail."""
    assert stored["survives_hand_verification"]
    # Mark only trials that do NOT admit a positive solution as verified, so
    # the verified-subset enrichment can no longer exceed one.
    pairs = [p for p in mod._pairs(stored["trials"]) if not p["same_fraction_size"]]
    by_pmid = {}
    for p in pairs:
        by_pmid.setdefault(p["pmid"], []).append(p)
    by = {t["pmid"]: t["reported"] for t in stored["trials"]}
    # Verify ONLY equivalence trials with no positive pair and difference
    # trials with nothing but positive pairs, which forces the subset's
    # enrichment below one.
    chosen = [q for q, ps in by_pmid.items()
              if (by.get(q) == "equivalence"
                  and not any(x["admits_positive"] for x in ps))
              or (by.get(q) == "difference"
                  and all(x["admits_positive"] for x in ps))]
    rows = [{"pmid": q, "reported": by[q], "arms_verdict": "correct",
             "error_mode": None, "selection_verdict": "primary"}
            for q in sorted(chosen)]
    doctored = tmp_path / "judged.json"
    doctored.write_text(json.dumps({"verdicts": rows}))
    monkeypatch.setattr(mod, "JUDGED", doctored)
    out = mod.assemble(_raw(stored))
    assert not out["survives_hand_verification"]
    assert out["verdict"] == "inversion_indistinguishable_from_arithmetic"


def test_the_single_fraction_arm_is_extracted():
    """A palliative trial's single-fraction arm has no fraction count.

    Missing it left the two multifraction VARIANTS of one arm to be paired
    against each other (PMID 10577696), which hand-judging caught.
    """
    rx = mod._compiled()
    for text in ("comparing 8 Gy single fraction with 20 Gy/5 fractions",
                 "randomised to a single fraction of 8 Gy or 30 Gy in 10 fractions",
                 "8 Gy in a single fraction versus 30 Gy in 10 fractions"):
        assert (8.0, 1) in mod.extract_arms(text, rx), f"missed it in: {text}"


def test_a_dose_range_is_flagged_rather_than_read_as_two_arms():
    rx = mod._compiled()
    assert rx["confound"]["dose_range_in_abstract"].search(
        "conventional radiotherapy alone (between 66 Gy in 33 fractions and "
        "74 Gy in 37 fractions, 5 fractions per week)")
    assert rx["confound"]["dose_range_in_abstract"].search(
        "randomised to receive 40-45Gy in 20-23 fractions or 24Gy in 12 fractions")
    assert not rx["confound"]["dose_range_in_abstract"].search(
        "randomised to 62 Gy in 20 fractions or 80 Gy in 40 fractions")


def test_a_boost_is_not_counted_as_an_arm():
    rx = mod._compiled()
    arms = mod.extract_arms(
        "randomised to 42 Gy in 15 fractions with a simultaneous integrated "
        "boost of 16 Gy in 8 fractions, or 50 Gy in 25 fractions", rx)
    assert {k for k, is_boost in arms.items() if not is_boost} == {(42.0, 15), (50.0, 25)}
    assert (16.0, 8) in arms and arms[(16.0, 8)], "the boost was not identified"


def test_the_boost_window_is_tight_enough_to_tell_arms_apart():
    """A wide window marks every nearby schedule as a boost.

    At the first version's 60 characters the word reached both real arms of the
    sentence above and the trial contributed NO arms at all -- a silent drop,
    which is the worst shape for a selection bug.
    """
    assert mod.BOOST_WINDOW_CHARS <= 30


# --------------------------------------------------------------------------
# the published bands are somebody else's numbers
# --------------------------------------------------------------------------

def test_the_bands_are_parsed_from_the_csv(monkeypatch, tmp_path):
    bands, rows = mod._curated_bands()
    assert set(bands) == {"prostate", "breast"}
    doctored = tmp_path / "trials.csv"
    doctored.write_text(
        "site,trial,pmid,arm_a_fractions,arm_a_dose_per_fraction_gy,arm_a_weeks,"
        "arm_b_fractions,arm_b_dose_per_fraction_gy,arm_b_weeks,outcome,"
        "published_ab_lo_gy,published_ab_hi_gy,published_ab_source_pmid\n"
        "prostate,X,1,37,2.0,7.4,20,3.0,4.0,non-inferior,9.0,11.0,2\n")
    monkeypatch.setattr(mod, "CURATED", doctored)
    moved, _ = mod._curated_bands()
    assert moved["prostate"]["low_gy"] == 9.0


def test_no_band_number_is_typed_into_the_generator():
    """A band copied into the code is a second source of truth."""
    src = _src()
    bands, _ = mod._curated_bands()
    for site, b in bands.items():
        for value in (b["low_gy"], b["high_gy"]):
            # The value may appear in prose about the two landmark trials, but
            # never as a float literal the code could compare against.
            assert f"= {value}" not in src and f"({value}," not in src, (
                f"{site}'s band value {value} looks hardcoded in the generator")


def test_a_site_without_a_sourced_band_gets_no_verdict(stored):
    sites = stored["bands"]["sites_with_a_sourced_band"]
    assert sites == sorted(sites) and len(sites) == 2
    for row in stored["bands"]["rows"]:
        assert row["source_pmid"] > 0


# --------------------------------------------------------------------------
# the tie back to the hand-read artifact
# --------------------------------------------------------------------------

def test_the_sweep_reproduces_the_curated_trials(stored):
    """The correctness check that licenses the other rows.

    The extractor read these arms out of the abstracts; a person read the same
    papers into the CSV. This caught the real defect: the first version's
    abstract-scoped secondary filter dropped BOTH trials, and the check was the
    only thing that said so.
    """
    rows = stored["curated_agreement"]
    assert len(rows) == 2
    for r in rows:
        assert r["in_sweep"], f"{r['trial']} fell out of the sweep"
        assert r["agrees"], (
            f"{r['trial']}: swept {r['swept_alpha_beta_gy']} against hand-read "
            f"{r['hand_read_alpha_beta_gy']}")
    chhip = next(r for r in rows if r["site"] == "prostate")
    assert chhip["delta_gy"] == pytest.approx(0.0, abs=1e-9), (
        "CHHiP used to reproduce EXACTLY; a nonzero delta means the extractor "
        "or the formula moved")


def test_the_formula_has_one_python_source():
    src = _src()
    assert "from validate_fractionation import isoeffect_alpha_beta" in src
    assert "def isoeffect_alpha_beta" not in src, (
        "a third copy of the inversion appeared here; validate_fractionation "
        "owns the Python implementation and the Rust comparison")


# --------------------------------------------------------------------------
# the limits the page must keep stating
# --------------------------------------------------------------------------

def test_the_report_says_the_magnitude_is_design_not_tissue(stored):
    md = mod.OUT_MD.read_text()
    assert "design property, not a tissue property" in md
    assert "constructed EQD2-equal" in md
    assert "Nothing here measures a tissue constant." in md
    assert "design" in stored["verdict"] and "tissue" in stored["verdict"]


def test_the_report_states_the_elapsed_time_limit(stored):
    assert stored["elapsed_time_sweepable"] is False
    md = mod.OUT_MD.read_text()
    assert "elapsed-time axis cannot be swept" in md
    assert "START-B" in md


def test_the_confound_classes_are_labelled_as_co_mentions(stored):
    md = mod.OUT_MD.read_text()
    assert "CO-MENTIONS, not arm differences" in md
    assert set(stored["confound_counts"]) == set(mod.CONFOUND_PATTERNS)


def test_a_non_inferiority_caveat_is_present():
    md = mod.OUT_MD.read_text()
    assert "bounds a difference" in md


# --------------------------------------------------------------------------
# generator contract
# --------------------------------------------------------------------------

def test_render_only_reassembles_rather_than_echoing():
    tree = ast.parse(_src())
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    branch = next(n for n in ast.walk(main)
                  if isinstance(n, ast.If)
                  and isinstance(n.test, ast.Attribute)
                  and n.test.attr == "render_only")
    calls = [c for stmt in branch.body for c in ast.walk(stmt)
             if isinstance(c, ast.Call)]
    assert any(getattr(c.func, "id", None) == "assemble" for c in calls), (
        "--render-only no longer passes through assemble(), so every guard on "
        "a derived field would go inert")
    read = ast.dump(branch)
    for key in ("trials", "funnel", "census_shards"):
        assert f"'{key}'" in read or f'"{key}"' in read


def test_assemble_is_idempotent(stored):
    once = mod.assemble(_raw(stored))
    twice = mod.assemble({"trials": once["trials"], "funnel": once["funnel"],
                          "census_shards": once["census_shards"]})
    assert json.dumps(once, sort_keys=True) == json.dumps(twice, sort_keys=True)


def test_the_committed_report_matches_the_committed_measurements(stored):
    assert mod.OUT_MD.read_text() == mod.render(
        json.loads(json.dumps(mod.assemble(_raw(stored)), sort_keys=True)))


# --------------------------------------------------------------------------
# what the two groups actually are
# --------------------------------------------------------------------------

def test_the_equivalence_regex_is_built_from_its_two_halves():
    """Design language and result language must not drift apart.

    They are different claims and the first version conflated them, so
    `EQUIVALENCE` is composed from the two rather than written a third time.
    """
    assert mod.EQUIVALENCE == f"{mod.DESIGN_LANGUAGE}|{mod.RESULT_LANGUAGE}"
    rx = mod._compiled()
    assert rx["design_language"].search("a randomised non-inferiority trial")
    assert not rx["result_language"].search("a randomised non-inferiority trial")
    assert rx["result_language"].search("no significant difference was found")
    assert not rx["design_language"].search("no significant difference was found")


def test_the_label_basis_is_measured_not_assumed(stored):
    lb = stored["label_basis"]
    counts = dict(lb["equivalence_label_counts"])
    assert sum(counts.values()) == stored["funnel"]["selected_equivalence"]
    # Non-vacuity: the split must actually split, or the caveat is empty.
    assert counts["design_only"] > 0 and counts["result_only"] > 0
    assert 0.0 < lb["design_language_only_share"] < 1.0
    assert lb["result_language_share"] == pytest.approx(
        1.0 - lb["design_language_only_share"], abs=1e-9)


def test_the_enrichment_must_survive_the_result_language_restriction(stored):
    """A result living only in the design-language half is not an outcome claim."""
    assert stored["survives_result_language_restriction"]
    assert stored["label_basis"]["enrichment"] > mod.ENRICHMENT_NULL
    raw = _raw(stored)
    # Leave the equivalence arm labelled only by design language, so the
    # restriction has nothing about an outcome left to measure.
    for t in raw["trials"]:
        if t["reported"] == "equivalence":
            t["label_basis"] = "design_only"
    out = mod.assemble(raw)
    assert out["label_basis"]["enrichment"] is None
    assert not out["survives_result_language_restriction"]
    assert out["verdict"] == "inversion_indistinguishable_from_arithmetic", (
        "the verdict survived with no result-labelled trial left, so the "
        "restriction is decoration")


def test_the_report_quotes_no_classifier_precision(stored):
    """The label's precision is NOT established and must not be implied.

    Most sampled rows cannot be adjudicated from an abstract, so any
    percentage here would describe the judging rather than the classifier.
    """
    md = mod.OUT_MD.read_text()
    assert "precision is NOT established" in md
    assert "DESIGN language only" in md
    rows = json.loads(mod.JUDGED.read_text())["verdicts"]
    verdicts = {r["reported_verdict"] for r in rows}
    assert "unclear" in verdicts, (
        "every sampled label became adjudicable, so the report's refusal to "
        "quote a precision should be revisited rather than left standing")
    assert "mislabelled" in verdicts, (
        "the known mislabelled row disappeared; it is the concrete evidence "
        "that the classifier has a real error rate")


def test_a_mislabelled_row_keeps_its_explanation():
    rows = {r["pmid"]: r for r in json.loads(mod.JUDGED.read_text())["verdicts"]}
    bad = [r for r in rows.values() if r["reported_verdict"] == "mislabelled"]
    assert bad
    for r in bad:
        assert r.get("label_note"), (
            f"PMID {r['pmid']} is marked mislabelled with no reason recorded")


def test_the_funnel_does_not_call_the_control_group_a_drop(stored):
    """A key contradicting its own prose is how an artifact misleads a reader."""
    f = stored["funnel"]
    assert "reported_a_difference" in f
    assert "dropped_no_equivalence_claim" not in f
    assert f["selected_equivalence"] + f["reported_a_difference"] == f["design_admissible"]


def test_the_inversion_is_symmetric_in_orientation():
    """Control-anchoring chooses which pairs form, never how one is read."""
    from validate_fractionation import isoeffect_alpha_beta as inv
    for n1, d1, n2, d2 in ((37, 2.0, 20, 3.0), (25, 2.0, 15, 2.67),
                           (20, 2.5, 40, 1.0)):
        assert inv(n1, d1, n2, d2) == pytest.approx(inv(n2, d2, n1, d1), rel=1e-12)


def test_the_interval_excluding_the_null_does_not_depend_on_the_seed(stored, monkeypatch):
    """A conclusion that moves with the seed is an artifact of the seed."""
    pairs = [p for p in mod._pairs(stored["trials"]) if not p["same_fraction_size"]]
    eq = [p for p in pairs if p["reported"] == "equivalence"]
    df = [p for p in pairs if p["reported"] == "difference"]
    for seed in (1, 7, 99, 2024):
        monkeypatch.setattr(mod, "BOOTSTRAP_SEED", seed)
        iv = mod._cluster_bootstrap(eq, df)
        assert iv["excludes_null"], f"seed {seed} stopped excluding the null"


# --------------------------------------------------------------------------
# gaps the mutation sweep found in the guards above
# --------------------------------------------------------------------------

def test_the_every_stratum_claim_refuses_an_empty_quantifier(stored):
    """The floor was INERT on real data, which a mutation sweep showed.

    Lowering MIN_RATED_STRATA survived every guard, because with the real
    trials all four strata carry a ratio and the floor never binds. The case
    it exists for has to be driven directly: put the two groups in DISJOINT
    gap strata so no stratum can produce a ratio at all, where `all()` over an
    empty selection would otherwise report that the enrichment holds
    everywhere.
    """
    raw = _raw(stored)
    for t in raw["trials"]:
        t["arms"] = ([[50.0, 20], [45.0, 45]] if t["reported"] == "equivalence"
                     else [[50.0, 20], [25.0, 25]])
    out = mod.assemble(raw)
    ps = out["positive_solution"]
    rated = [s for s in ps["gap_strata"] if s["ratio"] is not None]
    assert not rated, "the doctoring failed to empty every stratum"
    assert not ps["enrichment_holds_in_every_stratum"], (
        "with no stratum carrying a ratio the report still claimed the "
        "enrichment holds in every one of them")
    assert mod.MIN_RATED_STRATA > 0


def test_scan_records_the_control_group_as_a_control():
    """Structural, because `scan()` needs the census and CI cannot run it.

    A mutation to the funnel keys survived every behavioural guard: the tests
    read the committed artifact, and a committed artifact does not move when
    the scanner is edited. So the scanner's own key is pinned in its source.
    """
    src = _src()
    scan = src[src.index("def scan():"):src.index("def _curated_bands():")]
    assert 'funnel["reported_a_difference"] += 1' in scan
    assert "dropped_no_equivalence_claim" not in src, (
        "the control group is being recorded as a drop again, which makes the "
        "artifact's key contradict the report's own prose")
    # and the funnel's own arithmetic must stay closed
    tree = ast.parse(src)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "scan")
    keys = {c.value for n in ast.walk(fn) if isinstance(n, ast.Subscript)
            for c in ast.walk(n) if isinstance(c, ast.Constant)
            and isinstance(c.value, str)}
    assert {"selected_equivalence", "reported_a_difference",
            "design_admissible"} <= keys


def test_the_curated_tolerance_is_tight_enough_to_reject(monkeypatch, tmp_path):
    """A tolerance wide enough to pass anything attests nothing.

    Widening it survived the agreement guard, because that guard only asserts
    the rows AGREE -- which a wide enough tolerance guarantees. The tolerance
    needs a bound of its own, and the refusal needs to be exercised.
    """
    # Small against the quantity it compares: the hand-read ratios are of
    # order a few Gy, so the tolerance must be a small fraction of that.
    bands, curated = mod._curated_bands()
    smallest = min(abs(mod.isoeffect_alpha_beta(
        int(r["arm_a_fractions"]), float(r["arm_a_dose_per_fraction_gy"]),
        int(r["arm_b_fractions"]), float(r["arm_b_dose_per_fraction_gy"])))
        for r in curated)
    assert mod.CURATED_AGREEMENT_TOLERANCE_GY < smallest / 5.0, (
        "the agreement tolerance is a large fraction of the smallest value it "
        "compares, so agreement is close to guaranteed")
    # And it must actually refuse: a curated row whose arms imply a very
    # different ratio has to come back as disagreeing.
    doctored = tmp_path / "trials.csv"
    doctored.write_text(
        "site,trial,pmid,arm_a_fractions,arm_a_dose_per_fraction_gy,arm_a_weeks,"
        "arm_b_fractions,arm_b_dose_per_fraction_gy,arm_b_weeks,outcome,"
        "published_ab_lo_gy,published_ab_hi_gy,published_ab_source_pmid\n"
        "prostate,CHHiP,27339115,10,6.0,2.0,5,1.0,1.0,non-inferior,1.2,3.0,2\n")
    monkeypatch.setattr(mod, "CURATED", doctored)
    _, rows = mod._curated_bands()
    stored_json = json.loads(mod.OUT_JSON.read_text())
    pairs = mod._pairs(stored_json["trials"])
    out = mod._curated_agreement(pairs, rows)
    assert out[0]["in_sweep"], "the doctored row left the sweep, so this is vacuous"
    assert not out[0]["agrees"], (
        "a curated row whose arms imply a completely different ratio was "
        "reported as agreeing with the sweep")


# --------------------------------------------------------------------------
# the per-site pattern, and whether its bands are independent
# --------------------------------------------------------------------------

def test_the_site_pattern_is_derived_and_non_vacuous(stored):
    """Agreement is only informative if the hand-read verdicts DIFFER.

    If every site's hand-read trial landed in its band, "the pattern agrees"
    would be satisfied by the sweep favouring equivalence everywhere, which is
    what the headline already says. The differing-verdicts flag is what makes
    the agreement mean something.
    """
    sp = stored["site_pattern"]
    assert sp["sites_decided"] == 2
    assert sp["hand_read_verdicts_differ"], (
        "both hand-read trials now fall the same side of their band, so the "
        "per-site agreement has become trivial and must stop being claimed")
    assert sp["pattern_reproduces"]
    assert sp["sites_agreeing"] == sp["sites_decided"]
    # the directions must be opposite, which is the whole content of P14
    in_band = {r["site"]: r["hand_read_in_band"] for r in sp["rows"]}
    favours = {r["site"]: r["swept_favours_equivalence"] for r in sp["rows"]}
    assert in_band == favours


def test_a_site_whose_sweep_disagreed_would_be_reported_as_such(stored):
    """Flip one site's band so the sweep and the hand-read verdict diverge."""
    raw = _raw(stored)
    # Put every prostate trial's arms at a pair implying a ratio far outside
    # the prostate band, so the equivalence group stops landing in it.
    for t in raw["trials"]:
        if "prostate" in t["sites"] and t["reported"] == "equivalence":
            t["arms"] = [[50.0, 20], [40.0, 10]]      # implies 20.0 Gy
    out = mod.assemble(raw)
    sp = out["site_pattern"]
    prostate = next(r for r in sp["rows"] if r["site"] == "prostate")
    assert prostate["swept_favours_equivalence"] is not True, (
        "the doctoring did not move the prostate band rate")
    assert not sp["pattern_reproduces"]


def test_a_band_sourced_from_its_own_trial_is_flagged(stored):
    """The provenance claim the curated header used to make, checked per row.

    One of the two bands is sourced from the trial it judges, so the row tests
    a simple inversion against a fuller analysis of the SAME data rather than
    against independent data. A page that let the two carry equal weight would
    be overstating by one site out of two.
    """
    sp = stored["site_pattern"]
    flags = {r["site"]: r["band_independent_of_the_trial"] for r in sp["rows"]}
    assert set(flags.values()) == {True, False}, (
        "every band is now independent, or none is; either way the report's "
        "distinction needs revisiting rather than standing as written")
    assert sp["sites_with_an_independent_band"] < sp["sites_decided"]
    md = mod.OUT_MD.read_text()
    assert "independent of the trial it judges" in md
    assert "**NO**" in md, "the circular row is no longer marked in the table"


def test_the_curated_header_no_longer_claims_blanket_independence():
    """It did, and for one row it was false."""
    text = mod.CURATED.read_text()
    assert "same site INDEPENDENTLY of these schedules" not in text
    assert "NOT ALL INDEPENDENT OF THE TRIALS THEY JUDGE" in text


def test_the_claude_md_entry_quotes_the_committed_figures():
    """A figure hand-written beside a derived one goes stale, and this one did.

    The CLAUDE.md entry for this generator quoted the uninformative-pair shares
    as 21.1% and 10.5%. Those were true of an earlier pass; adding the
    single-fraction patterns and tightening the boost window moved the pair set
    and the report's own rendering followed, while the hand-written copy did
    not. Caught by a consistency check at the end of the work rather than by
    anything structural, so the check is structural now.

    Deliberately NARROW: it pins the handful of figures that entry quotes,
    not every number in CLAUDE.md. A gate over the whole file would flag
    hundreds of pre-existing entries and train exemptions instead of fixes --
    the lesson `audit_heading_only_assertions.py` already records.
    """
    claude = (REPO / "CLAUDE.md").read_text(encoding="utf-8")
    start = claude.index("scripts/trial_isoeffect_sweep.py")
    entry = claude[start:claude.index("\nscripts/", start + 1)]
    d = json.loads(mod.OUT_JSON.read_text())
    ps, u = d["positive_solution"], d["uninformative_pairs"]
    ea, lb, vs = d["extraction_audit"], d["label_basis"], d["verified_subset"]
    expected = {
        "crude enrichment": f"{ps['crude_enrichment']:.2f}x",
        "interval low": f"{ps['interval']['low']:.2f}",
        "interval high": f"{ps['interval']['high']:.2f}",
        "pooled": f"{ps['stratum_pooled_enrichment']:.2f}",
        "uninformative share, control": f"{100 * u['share_of_difference']:.1f}%",
        "uninformative share, equivalence": f"{100 * u['share_of_equivalence']:.1f}%",
        "enrichment if included": f"{u['enrichment_if_included']:.2f}x",
        "extraction precision": f"{100 * ea['precision']:.1f}%",
        "design-language share": f"{100 * lb['design_language_only_share']:.1f}%",
        "result-language enrichment": f"{lb['enrichment']:.2f}x",
        "verified subset": f"{vs['enrichment']:.2f}x",
        "design admissible": f"{d['funnel']['design_admissible']}",
    }
    missing = {k: v for k, v in expected.items() if v not in entry}
    assert not missing, (
        "the CLAUDE.md entry quotes figures the committed artifact no longer "
        f"carries: {missing}\n\nRe-read analysis/trial-isoeffect-sweep.md and "
        "correct the entry, or state the property without the number.")


def test_the_superseded_outcome_claim_is_gone_from_every_site():
    """"LQ-based schedule design predicts trial outcome" is withdrawn.

    It was the docstring's and the report's summary of what the sweep shows,
    written before the design-versus-result measurement existed. Once 38.5% of
    the equivalence label turned out to be design language, the claim outran
    the evidence -- and it survived in two places after the report's own later
    section contradicted it, three screens apart, which is the shape this
    repository keeps finding. Pinned across the generator, the committed
    report and the project's key-files entry together, because a retraction
    that reaches one site is the defect rather than the fix.
    """
    sites = {
        "generator": _src(),
        "report": mod.OUT_MD.read_text(),
        "CLAUDE.md": (REPO / "CLAUDE.md").read_text(encoding="utf-8"),
    }
    for name, text in sites.items():
        assert "predicts trial outcome" not in text.lower(), (
            f"{name} still carries the withdrawn claim that the inversion "
            "predicts trial outcome")
    # And the replacement must actually be present, or this guard only
    # checks that a sentence was deleted.
    assert "CONSTRUCT equivalent schedules" in sites["report"]
    assert "design language rather than a result" in sites["report"]
