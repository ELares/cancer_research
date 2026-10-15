"""Guards for `scripts/prcc_under_admissibility.py` (#944).

Every guard here DRIVES `assemble()` on a doctored input. Reading the committed
artifact instead is the defect this repository has recorded at least seven times
across three pull requests: a derived field's guard that reads the generated
file cannot fail when the generator changes, because the file does not move
until someone re-runs it.

The fixture is a small analytic world rather than a random one, so what each
guard expects is a consequence of the fixture's construction rather than of a
number someone read off a run. Its `published` half is COMPUTED from its own
rows, which makes the reproduction check self-consistent by construction -- so
the guards that matter doctor one side of it and require the comparison to be
withheld.
"""

import ast
import copy
import re
import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import prcc_under_admissibility as P  # noqa: E402
import intervals_under_admissibility as SIB  # noqa: E402
import run_prcc as RP  # noqa: E402
import sobol_sensitivity as SOBOL  # noqa: E402

SCRIPT = REPO / "scripts" / "prcc_under_admissibility.py"

NAMES = list(RP.PARAM_RANGES)
N_FIXTURE = 600  # comfortably above the conditioned floor, so that floor's own guard drives it

# `"nan" in md` matches the word "dominant". Word boundaries, in one place, so
# the two call sites cannot drift apart.
_NOT_A_NUMBER = re.compile(r"\b(?:nan|inf|-inf)\b", re.I)


def _no_unrenderable_number(md):
    hit = _NOT_A_NUMBER.search(md)
    assert hit is None, f"rendered page carries {hit.group(0)!r} at {hit.start()}"


def _frac(X, name):
    lo, hi = RP.PARAM_RANGES[name]
    return (X[:, NAMES.index(name)] - lo) / (hi - lo)


def _fixture_rows():
    """An analytic world with the three shapes the page has to handle.

    * the untreated arms rise with `gsh_scav_efficiency`, so admissibility
      selects on that margin and on nothing else -- which is also what the real
      data does, and is what makes the circularity check inert on it;
    * `Persister/SDT` is SATURATED on the full design and de-saturates under
      conditioning, the shape the manuscript's bound is quoted at;
    * `Glycolytic/RSL3` is constant on the admissible subset, so it has no
      conditioned ranking at all.
    """
    X = RP.latin_hypercube_sample(N_FIXTURE, list(RP.PARAM_RANGES.values()), seed=7)
    gsh = _frac(X, "gsh_scav_efficiency")
    sdt = _frac(X, "sdt_ros")
    lpp = _frac(X, "lp_propagation")
    lpr = _frac(X, "lp_rate")
    gpx = _frac(X, "gpx4_rate")
    conditions = [P._label(c) for c in P._all_conditions()]
    rows = []
    for i in range(N_FIXTURE):
        # EVERY condition must depend on at least two parameters. A condition
        # that is an exact monotone function of ONE makes its ranked column
        # perfectly rank-correlated with that input, and `compute_prcc` inverts
        # the rank-correlation matrix -- so the fixture's first version raised
        # `LinAlgError: Singular matrix` rather than testing anything.
        worst = 0.05 * gsh[i] + 0.012 * lpp[i]
        vals = {
            "Persister/SDT": 1.0 - 0.001 * (1 - sdt[i]) - 0.02 * (1 - gsh[i]),
            "Persister/RSL3": 0.2 + 0.5 * lpp[i] - 0.2 * gpx[i],
            # Strongly convex, so the full design carries a
            # Pearson-against-Spearman divergence and the escalation flag set
            # is NON-EMPTY -- without one, a guard asserting that the right
            # number of flags on the wrong parameters is not a reproduction
            # compares two empty lists and passes for that reason. Convex in
            # `lp_rate` rather than in a DRUG parameter: the circularity
            # fixture makes admissibility track `sdt_ros`, and a ninth power of
            # a parameter restricted to its low end collapses to a constant, so
            # that column became a function of ONE parameter on the admissible
            # subset and the estimator went singular. Fourth time for that
            # trap, in a new disguise.
            "Glycolytic/SDT": 0.5 + 0.3 * lpr[i] ** 9 - 0.1 * sdt[i],
            "Glycolytic/RSL3": (0.0 if worst < 0.02
                               else 0.001 * lpp[i] + 0.0005 * gpx[i]),
            "Persister/Control": worst,
            "PersisterNrf2/SDT": 0.9 - 0.05 * gpx[i] + 0.02 * sdt[i],
            "Glycolytic/Control": 0.5 * worst,
            "PersisterNrf2/Control": 0.9 * worst,
        }
        row = {"draw": str(i)}
        for j, nm in enumerate(NAMES):
            row[nm] = f"{X[i, j]:.12g}"
        for label in conditions:
            row[label] = f"{vals[label]:.12g}"
        for cells in P.RESOLUTION_CELLS:
            for seed in P.RESOLUTION_SEEDS:
                # Coarser cell counts admit MORE draws here, the direction the
                # real measurement found; the seed moves the boundary slightly.
                nudge = (0.004 if cells == min(P.RESOLUTION_CELLS) else 0.0)
                nudge += 0.0005 * (seed % 3)
                row[f"worst_untreated_n{cells}_seed{seed}"] = f"{max(worst - nudge, 0.0):.12g}"
        rows.append(row)
    return rows, X


def _published_from(rows, X):
    """The published artifact, computed from the fixture's own rows."""
    conditions = []
    for c in RP.CONDITIONS:
        label = P._label(c)
        y = np.array([float(r[label]) for r in rows], float)
        coeff = {nm: v for nm, (v, _p) in zip(NAMES, P.compute_prcc(X, y))}
        conditions.append({
            "phenotype": c[0], "treatment": c[1], "output": "death_rate",
            "mean_death_rate": float(y.mean()), "std_death_rate": float(y.std()),
            "prcc": coeff, "p_values": {nm: 0.0 for nm in NAMES},
            "ranked_by_abs_prcc": sorted(NAMES, key=lambda nm: -abs(coeff[nm])),
        })
    # The escalation flags, computed the same way the generator does, so the
    # fixture's published half is SELF-CONSISTENT. Leaving them empty made the
    # baseline fail to reproduce the moment a fixture condition started
    # flagging, which is the fixture lying rather than the code being wrong.
    from scipy.stats import pearsonr, spearmanr

    flags = []
    for c in RP.CONDITIONS:
        y = np.array([float(r[P._label(c)]) for r in rows], float)
        if y.std() < P.DEGENERATE_STD:
            continue
        for j, nm in enumerate(NAMES):
            div = abs(float(pearsonr(X[:, j], y)[0])
                      - float(spearmanr(X[:, j], y)[0]))
            if div > 0.2:
                flags.append({"condition": f"{c[0]} \u00d7 {c[1]}",
                              "parameter": nm, "divergence": round(div, 3)})
    return {
        "metadata": {"n_samples": len(rows), "n_cells_per_condition": 1000,
                     "seed": 42, "n_parameters": len(NAMES), "parameters": NAMES,
                     "parameter_ranges": {nm: list(RP.PARAM_RANGES[nm])
                                          for nm in NAMES}},
        "conditions": conditions,
        "sobol_decision": {"warranted": len(flags) > 2,
                           "divergence_threshold": 0.2,
                           "non_monotonic_flags": flags, "n_flags": len(flags),
                           "criterion": "x", "reason": "y"},
    }


def _sobol_rows():
    """Flags whose straddling counts are a property of the fixture: the first
    parameter straddles on every row, the LAST TWO on none.

    Two rather than one, matching the real measurement -- with a single
    zero-straddling parameter there is no order to publish, and a mutation
    reverting that list to stored dict order survived for exactly that reason.
    A fixture that cannot express the defect is not a test of it.
    """
    sob = [p[0] for p in SOBOL.PARAMS]
    rows = []
    for r in range(64):
        row = {"row": str(r), "A": "1", "B": "0" if r % 2 else "1"}
        for i, nm in enumerate(sob):
            if i == 0:
                row[f"AB_{nm}"] = "0"          # straddles every row
            elif i >= len(sob) - 2:
                row[f"AB_{nm}"] = "1"          # straddles none
            else:
                row[f"AB_{nm}"] = "0" if r % 4 == 0 else "1"
        rows.append(row)
    return rows


def _manuscript():
    return {
        "sdt": {"condition": "Persister/SDT", "dominant": "sdt_ros",
                "dominant_quoted": 0.63, "others_below": 0.16,
                "sentence": "For Persister x SDT, a single parameter dominates"},
        "rsl3": {"condition": "Persister/RSL3",
                 "named": [("lp_propagation", 0.92), ("gpx4_rate", -0.87)],
                 "count_above": 2, "above": 0.4,
                 "sentence": "For Persister x RSL3, the LP cascade parameters"},
    }


@pytest.fixture
def raw():
    rows, X = _fixture_rows()
    return {"rows": rows, "sobol_rows": _sobol_rows(),
            "sobol_control": [{"condition": P._label(c), "low": 0.0, "high": 0.0,
                               "identical": True} for c in SIB.CONTROL_CONDITIONS],
            "published": _published_from(rows, X), "manuscript": _manuscript()}


def _design_patch(monkeypatch, raw):
    """`assemble` reads the design from the committed artifact; point it at the
    fixture's instead, so a guard never silently mixes the two."""
    monkeypatch.setattr(P, "_design", lambda: {
        "n_samples": raw["published"]["metadata"]["n_samples"],
        "n_cells": raw["published"]["metadata"]["n_cells_per_condition"],
        "seed": raw["published"]["metadata"]["seed"],
        "parameters": list(raw["published"]["metadata"]["parameters"])})
    monkeypatch.setattr(P, "_published", lambda: raw["published"])


# --------------------------------------------------------------------------
# Reproduction
# --------------------------------------------------------------------------
def test_the_fixture_reproduces_so_the_doctored_cases_below_mean_something(
        monkeypatch, raw):
    _design_patch(monkeypatch, raw)
    d = P.assemble(raw)
    assert d["reproduction"]["reproduces"], d["reproduction"]
    assert d["reproduction"]["n_coefficients"] == len(RP.CONDITIONS) * len(NAMES)


def test_a_design_that_does_not_reproduce_withholds_the_comparison(
        monkeypatch, raw):
    """A page comparing against a re-run rather than against what is published
    is not correcting anything, so the state has to be visible in the prose."""
    _design_patch(monkeypatch, raw)
    bad = copy.deepcopy(raw)
    row = bad["published"]["conditions"][0]
    row["prcc"]["lp_rate"] = row["prcc"]["lp_rate"] + 0.3
    d = P.assemble(bad)
    assert not d["reproduction"]["reproduces"]
    assert d["reproduction"]["max_coefficient_deviation"] > P.REPRODUCTION_TOLERANCE
    assert "withheld" in P.render(d)


def test_the_reproduction_bound_sits_above_the_float_noise_and_below_a_readable_difference():
    """Two-directional, because a tolerance loose enough to accept a real
    discrepancy and one tight enough to reject a reproduction are the same
    defect mirrored. Both sides are read from the committed run rather than
    compared against a number typed here."""
    measured = json.loads((REPO / "analysis" / "prcc-under-admissibility.json")
                          .read_text())["reproduction"]["max_coefficient_deviation"]
    assert measured <= P.REPRODUCTION_TOLERANCE, (
        "the committed run does not reproduce under its own bound")
    smallest_published = min(
        abs(v) for c in json.loads((REPO / "analysis" / "prcc-results.json").read_text())
        ["conditions"] for v in c["prcc"].values() if v != 0.0)
    assert P.REPRODUCTION_TOLERANCE < smallest_published / 1000, (
        "the bound is not far below the smallest published coefficient, so it "
        "could accept a difference a reader would act on")


# --------------------------------------------------------------------------
# The design constants that cannot be imported
# --------------------------------------------------------------------------
def test_the_draw_seed_offset_is_read_from_the_published_generator():
    assert P._pinned_seed_offset() == P.DRAW_SEED_OFFSET, (
        "run_prcc's per-draw seed moved; every figure on the page moves with it")


def test_a_draw_seed_this_parser_cannot_read_raises_rather_than_falling_back(
        monkeypatch, tmp_path):
    """Falling back to the pinned value is how a design change stops being
    noticed, so the shapes are driven one at a time."""
    (tmp_path / "scripts").mkdir()
    src = (REPO / "scripts" / "run_prcc.py").read_text()
    for broken in ("seed=i + offset_from_somewhere", "seed=i"):
        doctored = src.replace(f"seed=i + {P.DRAW_SEED_OFFSET}", broken, 1)
        assert doctored != src
        (tmp_path / "scripts" / "run_prcc.py").write_text(doctored)
        monkeypatch.setattr(P, "REPO", tmp_path)
        with pytest.raises(ValueError):
            P._pinned_seed_offset()
    (tmp_path / "scripts" / "run_prcc.py").write_text(
        src.replace("result = fc.sim_batch(", "result = other(", 1))
    monkeypatch.setattr(P, "REPO", tmp_path)
    with pytest.raises(ValueError):
        P._pinned_seed_offset()


def test_the_constant_output_rule_is_the_published_generator_s_own():
    assert P._pinned_degenerate_std() == P.DEGENERATE_STD


def test_this_page_judges_admissibility_exactly_as_its_siblings_do():
    """Two definitions of admissible that can drift apart is the whole reason
    the constraint is imported rather than restated."""
    assert P.CONTROL_CONDITIONS is SIB.CONTROL_CONDITIONS
    assert P.DRUG_PARAMETERS is SIB.DRUG_PARAMETERS
    assert P.STRONG_SELECTION == SIB.STRONG_SELECTION
    src = SCRIPT.read_text()
    tree = ast.parse(src)
    defined = {t.id for n in tree.body if isinstance(n, ast.Assign)
               for t in n.targets if isinstance(t, ast.Name)}
    for forbidden in ("CONTROL_CONDITIONS", "BASELINE_MAX", "DRUG_PARAMETERS",
                      "PARAM_RANGES", "CONDITIONS", "STRONG_SELECTION"):
        assert forbidden not in defined, (
            f"{forbidden} is defined here as well as imported, so the two can "
            f"disagree")


# --------------------------------------------------------------------------
# The manuscript's figures
# --------------------------------------------------------------------------
def test_the_manuscript_figures_are_parsed_from_the_manuscript():
    claims = P.manuscript_claims()
    body = (REPO / "article" / "drafts" / "v1.md").read_text()
    for half in ("sdt", "rsl3"):
        assert claims[half]["sentence"] in body
    assert claims["sdt"]["dominant"] in NAMES
    assert [nm for nm, _v in claims["rsl3"]["named"]] == [
        nm for nm, _v in claims["rsl3"]["named"] if nm in NAMES]


def test_moving_the_manuscript_bound_moves_the_verdict(monkeypatch, raw):
    """The decisive guard against a bound hardcoded here: the fixture's SDT
    others sit at known magnitudes, so a bound below them must fail the claim
    and a bound above them must pass it, with NOTHING else changed."""
    _design_patch(monkeypatch, raw)
    base = P.assemble(raw)
    others = [abs(v) for nm, v in
              base["per_condition"]["Persister/SDT"]["admissible"]["coefficients"].items()
              if nm != "sdt_ros"]
    high, low = max(others), min(others)
    loose = copy.deepcopy(raw)
    loose["manuscript"]["sdt"]["others_below"] = high * 2
    tight = copy.deepcopy(raw)
    tight["manuscript"]["sdt"]["others_below"] = low / 2
    assert P.assemble(loose)["scored"]["sdt_bound"]["holds_admissible"] is True
    assert P.assemble(tight)["scored"]["sdt_bound"]["holds_admissible"] is False


def test_moving_the_manuscript_count_moves_the_verdict(monkeypatch, raw):
    _design_patch(monkeypatch, raw)
    base = P.assemble(raw)
    measured = base["scored"]["rsl3_count"]["admissible"]
    agree = copy.deepcopy(raw)
    agree["manuscript"]["rsl3"]["count_above"] = measured
    disagree = copy.deepcopy(raw)
    disagree["manuscript"]["rsl3"]["count_above"] = measured + 1
    assert P.assemble(agree)["scored"]["rsl3_count"]["holds_admissible"] is True
    assert P.assemble(disagree)["scored"]["rsl3_count"]["holds_admissible"] is False
    # And the BOUND the count is taken over, which the version above left
    # unguarded: a mutation hardcoding it SURVIVED, because moving only the
    # claimed count never consults it. Every parameter must count when the bar
    # is on the floor and none when it is above every coefficient.
    low, high = copy.deepcopy(raw), copy.deepcopy(raw)
    low["manuscript"]["rsl3"]["above"] = 0.0
    high["manuscript"]["rsl3"]["above"] = 1.0
    assert P.assemble(low)["scored"]["rsl3_count"]["admissible"] == len(NAMES)
    assert P.assemble(high)["scored"]["rsl3_count"]["admissible"] == 0


def test_a_manuscript_sentence_that_no_longer_parses_raises(monkeypatch, tmp_path):
    """The page scores those two figures, so it must break rather than carry a
    stale quotation of a sentence someone edited."""
    body = (REPO / "article" / "drafts" / "v1.md").read_text()
    claims = P.manuscript_claims()
    for half in ("sdt", "rsl3"):
        doctored = tmp_path / f"v1-{half}.md"
        doctored.write_text(body.replace(claims[half]["sentence"], "edited away", 1))
        monkeypatch.setattr(P, "MANUSCRIPT", doctored)
        with pytest.raises(ValueError):
            P.manuscript_claims()


# --------------------------------------------------------------------------
# Conditioning
# --------------------------------------------------------------------------
def test_a_condition_constant_under_conditioning_publishes_no_ranking(
        monkeypatch, raw):
    """A rank partial correlation against a constant is undefined, and the one
    thing that must never happen is a ranking sorted on NaN."""
    _design_patch(monkeypatch, raw)
    d = P.assemble(raw)
    c = d["per_condition"]["Glycolytic/RSL3"]
    assert c["degenerate_under_conditioning"] is True
    assert c["withheld_for_too_few_draws"] is False, (
        "withheld for the wrong reason: this condition is constant, not thin")
    assert c["admissible"]["ranking"] is None
    assert c["admissible"]["coefficients"] is None
    md = P.render(d)
    assert "undefined" in md
    _no_unrenderable_number(md)


def test_no_rendered_figure_is_nan_when_a_condition_degenerates(monkeypatch, raw):
    _design_patch(monkeypatch, raw)
    flat = copy.deepcopy(raw)
    # Constant on the ADMISSIBLE subset only. Flattening the whole column makes
    # the column constant on the full design too, and `compute_prcc` inverts a
    # rank-correlation matrix -- so the published half could not have been
    # produced either, and the fixture raises instead of testing the page.
    for r in flat["rows"]:
        if float(r["Persister/Control"]) < P.BASELINE_MAX:
            r["Persister/RSL3"] = "0.4"
    flat["published"] = _published_from(
        flat["rows"], np.array([[float(r[nm]) for nm in NAMES] for r in flat["rows"]]))
    d = P.assemble(flat)
    assert d["per_condition"]["Persister/RSL3"]["degenerate_under_conditioning"]
    _no_unrenderable_number(P.render(d))


def test_the_conditioned_coefficients_are_recomputed_not_echoed(monkeypatch, raw):
    """Moving one admissible draw's outcome must move the conditioned
    coefficient and leave the full-design one alone only if the draw is
    inadmissible -- so the two halves cannot be the same read."""
    _design_patch(monkeypatch, raw)
    base = P.assemble(raw)
    moved = copy.deepcopy(raw)
    before = base["per_condition"]["Persister/RSL3"]["admissible"]["coefficients"]
    for r in moved["rows"]:
        if float(r["Persister/Control"]) < P.BASELINE_MAX:
            r["Persister/RSL3"] = str(1.0 - float(r["Persister/RSL3"]))
    after = P.assemble(moved)["per_condition"]["Persister/RSL3"]["admissible"]["coefficients"]
    assert any(abs(after[nm] - before[nm]) > 0.1 for nm in NAMES)


def test_the_arm_conditioned_on_is_marked_as_a_diagnostic(monkeypatch, raw):
    """Selecting draws on an outcome and then reporting that outcome's
    sensitivity is selection on the dependent variable."""
    _design_patch(monkeypatch, raw)
    d = P.assemble(raw)
    marked = [k for k, v in d["per_condition"].items() if v["conditioned_on"]]
    assert marked == ["Persister/Control"], marked
    assert "the arm conditioned on" in P.render(d)
    assert "dependent variable" in P.render(d)


def test_one_dead_untreated_arm_is_enough_to_make_a_draw_inadmissible(
        monkeypatch, raw):
    """The constraint is on the WORST untreated arm, not on an average or a
    best case. A draw whose Persister control is dead while the other two
    phenotypes survive must be inadmissible -- the rule the manuscript states
    and the one the sibling pages apply."""
    _design_patch(monkeypatch, raw)
    doctored = copy.deepcopy(raw)
    arms = [P._label(c) for c in SIB.CONTROL_CONDITIONS]
    dead = "Persister/Control"
    assert dead in arms
    for r in doctored["rows"]:
        lpp = float(r["lp_propagation"])
        gpx = float(r["gpx4_rate"])
        # TWO parameters, and varying. Varying because a constant column is one
        # the published artifact could not hold and `assemble` now says so; two
        # because a column that is an exact monotone function of ONE makes the
        # rank-correlation matrix the estimator inverts singular. That trap has
        # now caught three fixtures in this file, which is why it is spelled
        # out here rather than left to the next one.
        r[dead] = f"{0.5 + 0.001 * lpp + 0.0005 * gpx:.12g}"
        for other in arms:
            if other != dead:
                r[other] = f"{0.0001 * lpp + 0.00005 * gpx:.12g}"
    d = P.assemble(doctored)
    assert d["admissibility"]["n_admissible"] == 0, (
        "a draw with a dead untreated phenotype is being admitted, so "
        "admissibility is not judged on the worst arm")
    assert d["conditioned_floor"]["bound"] is True
    md = P.render(d)
    assert "No conditioned ranking is published at all" in md
    _no_unrenderable_number(md)


def test_too_few_admissible_draws_withholds_every_conditioned_ranking(
        monkeypatch, raw):
    """INERT on the real data -- 841 draws survive against a floor well below
    that, asserted below -- so it is driven by a fixture, an inert floor
    reachable only that way being a defect class this repository has on
    record."""
    _design_patch(monkeypatch, raw)
    live = json.loads((REPO / "analysis" / "prcc-under-admissibility.json").read_text())
    assert live["conditioned_floor"]["bound"] is False, (
        "the floor binds on the real data; this guard's premise is gone")
    assert live["conditioned_floor"]["n_admissible"] > live["conditioned_floor"]["floor"]
    d0 = P.assemble(raw)
    assert d0["conditioned_floor"]["bound"] is False
    assert all(c["admissible"] is not None or c["degenerate_under_conditioning"]
               for c in d0["per_condition"].values())
    # Just under the floor: every conditioned ranking is withheld, and the
    # reason is distinguishable from a constant output.
    doctored = copy.deepcopy(raw)
    floor = d0["conditioned_floor"]["floor"]
    kept = 0
    for r in doctored["rows"]:
        if float(r["Persister/Control"]) < P.BASELINE_MAX and kept < floor - 1:
            kept += 1
            continue
        lpp = float(r["lp_propagation"])
        for c in SIB.CONTROL_CONDITIONS:
            r[P._label(c)] = f"{0.5 + 0.001 * lpp:.12g}"
    d = P.assemble(doctored)
    assert d["admissibility"]["n_admissible"] == floor - 1
    assert d["conditioned_floor"]["bound"] is True
    for label, c in d["per_condition"].items():
        # The MOMENTS survive -- a subset of any size has a mean -- and only the
        # RANKING is withheld. Asserting the whole half is absent would be
        # asserting the absence-as-zero defect this split exists to remove.
        assert c["withheld_for_too_few_draws"] is True, label
        assert c["admissible"]["coefficients"] is None, label
        assert c["admissible"]["ranking"] is None, label
        assert c["admissible"]["mean"] is not None, label


def test_the_floor_is_the_larger_of_the_imported_and_the_arithmetic_one():
    """A floor chosen here would be a second opinion about a question the
    sibling pages already answer, and an arithmetic floor alone would publish a
    ranking with a single degree of freedom."""
    assert P._conditioned_floor(11) == max(SIB.MIN_ADMISSIBLE_DRAWS, 13)
    assert P._conditioned_floor(10_000) == 10_002, (
        "the arithmetic floor stops binding above the imported one")


# --------------------------------------------------------------------------
# The circularity falsifier -- inert on the real data, so driven here
# --------------------------------------------------------------------------
def test_a_selected_drug_parameter_breaks_the_circularity_check(monkeypatch, raw):
    """Measured, no drug parameter is selected on. A check nothing in the data
    can fail has to be driven by a fixture or it is decoration."""
    _design_patch(monkeypatch, raw)
    assert P.assemble(raw)["circularity_clear"] is True
    drug = P.DRUG_PARAMETERS[0]
    doctored = copy.deepcopy(raw)
    lo, hi = RP.PARAM_RANGES[drug]
    lo2, hi2 = RP.PARAM_RANGES["fsp1_rate"]
    for r in doctored["rows"]:
        # Admissibility now tracks the DRUG parameter instead of a defence.
        # The second term is not decoration: a control arm that is an exact
        # monotone function of ONE parameter makes the published half's
        # rank-correlation matrix singular, so the fixture would raise.
        u = (float(r[drug]) - lo) / (hi - lo)
        v = (float(r["fsp1_rate"]) - lo2) / (hi2 - lo2)
        for c in SIB.CONTROL_CONDITIONS:
            r[P._label(c)] = f"{0.05 * u + 0.004 * v:.12g}"
    d = P.assemble(doctored)
    assert d["circularity_clear"] is False
    assert d["selection"][drug]["strongly_selected"] is True
    assert "NOT clear" in P.render(d)


def test_the_artifact_never_carries_a_nan_or_an_infinity(monkeypatch, raw):
    """`json.dumps` writes NaN as the bare token `NaN`, which is not JSON.

    So a NaN here does not merely read badly on the page -- it makes the
    committed artifact unparseable by a strict reader, and every figure derived
    from an empty slice produces one. Driven at the state that generates them.
    """
    _design_patch(monkeypatch, raw)
    empty = copy.deepcopy(raw)
    for r in empty["rows"]:
        lpp, gpx = float(r["lp_propagation"]), float(r["gpx4_rate"])
        for c in SIB.CONTROL_CONDITIONS:
            r[P._label(c)] = f"{0.5 + 0.001 * lpp + 0.0005 * gpx:.12g}"
    for case in (raw, empty):
        d = P.assemble(case)
        text = json.dumps(d, indent=2, sort_keys=True)
        for token in ("NaN", "Infinity", "-Infinity"):
            assert token not in text, f"{token} in the artifact"
        json.loads(text, parse_constant=lambda c: (_ for _ in ()).throw(
            AssertionError(f"non-JSON constant {c!r} in the artifact")))
        _no_unrenderable_number(P.render(d))


def test_a_withheld_falsifier_is_not_reported_as_passing(monkeypatch, raw):
    """A selection that kept nothing has not selected on a drug parameter, and
    reading that as the check passing would be reading an absence as a
    result -- the defect class this repository already has on record."""
    _design_patch(monkeypatch, raw)
    assert P.assemble(raw)["circularity_clear"] is True
    empty = copy.deepcopy(raw)
    for r in empty["rows"]:
        lpp, gpx = float(r["lp_propagation"]), float(r["gpx4_rate"])
        for c in SIB.CONTROL_CONDITIONS:
            r[P._label(c)] = f"{0.5 + 0.001 * lpp + 0.0005 * gpx:.12g}"
    d = P.assemble(empty)
    assert d["circularity_clear"] is None, (
        "a falsifier with nothing to test is reported as a verdict")
    md = P.render(d)
    assert "WITHHELD" in md
    assert "**Clear.**" not in md


# --------------------------------------------------------------------------
# Stratification
# --------------------------------------------------------------------------
def test_the_stratification_loss_is_measured_in_both_directions(monkeypatch, raw):
    """A filtered sample that happens to stay uniform must NOT be reported as
    having lost stratification, or the measurement is just a restatement of the
    admissible count."""
    _design_patch(monkeypatch, raw)
    d = P.assemble(raw)
    n = d["admissibility"]["n_draws"]
    for nm, v in d["stratification"].items():
        assert v["full"] < 5.0 / n, (
            f"{nm}: the full design is not uniform to a few draw-widths, so "
            f"the comparison against the admissible subset means nothing")
    # The bar is DERIVED from the admissible count, so halving the surviving
    # draws must move it rather than leave a chosen constant in place.
    null = d["stratification_null"]["admissible"]
    expected = P.KS_COEFFICIENT / math.sqrt(d["admissibility"]["n_admissible"])
    assert abs(null - expected) < 1e-12
    for nm, v in d["stratification"].items():
        assert abs(v["times_null"] - v["admissible"] / null) < 1e-9, nm
    assert d["stratification"]["gsh_scav_efficiency"]["admissible"] > \
        d["stratification"][P.DRUG_PARAMETERS[0]]["admissible"], (
        "the margin admissibility selects on is not the most deformed one")
    # The negative case: conditioning that removes NOTHING must not be
    # reported as having deformed anything. Selecting every other row instead
    # would compare a bar of a few hundredths against the sampling noise of a
    # random half, which is the same order -- so the test would be measuring
    # its own fixture's luck.
    spread = copy.deepcopy(raw)
    for r in spread["rows"]:
        for c in SIB.CONTROL_CONDITIONS:
            # Scaled under the constraint rather than flattened to zero: a
            # constant column is one the published artifact could not hold, and
            # `assemble` now says so instead of admitting it.
            r[P._label(c)] = f"{float(r[P._label(c)]) * 0.1:.12g}"
    d2 = P.assemble(spread)
    assert d2["admissibility"]["n_admissible"] == d2["admissibility"]["n_draws"]
    assert d2["stratification_lost"] == [], d2["stratification_lost"]
    # And the null itself must track the subset size, not sit at a constant.
    assert d2["stratification_null"]["admissible"] < null, (
        "the same-size null did not tighten when more draws survived, so it is "
        "not derived from the admissible count")


# --------------------------------------------------------------------------
# The escalation decision published in the same artifact
# --------------------------------------------------------------------------
def test_the_escalation_recomputation_ignores_the_untreated_arms(monkeypatch, raw):
    """A first pass recomputed this over the untreated arms as well and got one
    flag more than the artifact publishes, which would have been a reproduction
    failure reported as a conditioning effect."""
    _design_patch(monkeypatch, raw)
    base = P.assemble(raw)["escalation"]["full"]["n_flags"]
    doctored = copy.deepcopy(raw)
    untreated = [P._label(c) for c in SIB.CONTROL_CONDITIONS
                 if P._label(c) not in [P._label(x) for x in RP.CONDITIONS]]
    assert untreated, "the fixture has no untreated arm outside the design"
    X = np.array([[float(r[nm]) for nm in NAMES] for r in raw["rows"]], float)
    gsh = _frac(X, "gsh_scav_efficiency")
    # ONE dimension: a strongly convex (so Pearson and Spearman diverge) but
    # still sub-`Persister/Control` arm, which leaves the max -- and therefore
    # every admissibility verdict -- bit-identical. Doctoring it upward would
    # change which draws are admissible and the guard would be testing two
    # things at once.
    for i, r in enumerate(doctored["rows"]):
        val = float(r["Persister/Control"]) * 0.9 * gsh[i] ** 11
        for label in untreated:
            r[label] = f"{val:.12g}"
    d = P.assemble(doctored)
    assert (d["admissibility"]["n_admissible"]
            == P.assemble(raw)["admissibility"]["n_admissible"])
    # Non-vacuity: the doctored arm really would flag if it were consulted.
    from scipy.stats import pearsonr, spearmanr
    y = np.array([float(r[untreated[0]]) for r in doctored["rows"]], float)
    j = NAMES.index("gsh_scav_efficiency")
    divergence = abs(float(pearsonr(X[:, j], y)[0])
                     - float(spearmanr(X[:, j], y)[0]))
    threshold = raw["published"]["sobol_decision"]["divergence_threshold"]
    assert divergence > threshold, (
        f"the doctored untreated arm does not diverge ({divergence:.3f}), so "
        f"this guard would pass whatever the recomputation consulted")
    assert d["escalation"]["full"]["n_flags"] == base


def test_the_escalation_decision_must_reproduce_before_it_is_read(monkeypatch, raw):
    _design_patch(monkeypatch, raw)
    assert P.assemble(raw)["escalation"]["reproduces"] is True
    bad = copy.deepcopy(raw)
    bad["published"]["sobol_decision"]["n_flags"] = 99
    assert P.assemble(bad)["escalation"]["reproduces"] is False
    # And the harder case, which a count comparison cannot see: the RIGHT
    # NUMBER of flags on the WRONG parameters is not a reproduction.
    live = P.assemble(raw)["escalation"]["full"]["flags"]
    swapped = copy.deepcopy(raw)
    swapped["published"]["sobol_decision"]["n_flags"] = len(live)
    swapped["published"]["sobol_decision"]["non_monotonic_flags"] = [
        {"condition": f["condition"].replace("/", " \u00d7 "),
         "parameter": "fenton_rate" if f["parameter"] != "fenton_rate"
                      else "lp_rate",
         "divergence": f["divergence"]} for f in live]
    d = P.assemble(swapped)
    assert d["escalation"]["full"]["n_flags"] == d["escalation"]["published"]["n_flags"]
    assert d["escalation"]["reproduces"] is False, (
        "the right number of flags on the wrong parameters counts as a "
        "reproduction")
    # The label spelling must NOT be what breaks it: the two artifacts write a
    # condition differently and the comparison normalises for that.
    same = copy.deepcopy(raw)
    same["published"]["sobol_decision"]["n_flags"] = len(live)
    same["published"]["sobol_decision"]["non_monotonic_flags"] = [
        {"condition": f["condition"].replace("/", " \u00d7 "),
         "parameter": f["parameter"], "divergence": f["divergence"]}
        for f in live]
    assert P.assemble(same)["escalation"]["reproduces"] is True


def test_a_conditioned_escalation_flip_is_reported_as_one(monkeypatch, raw):
    """The published decision is one flag from flipping, so the page must be
    able to SAY it flipped. As measured it does not, which is why this is
    driven rather than observed."""
    _design_patch(monkeypatch, raw)
    assert P.assemble(raw)["escalation"]["survives_conditioning"] is True
    doctored = copy.deepcopy(raw)
    rng = np.random.default_rng(11)
    for r in doctored["rows"]:
        if float(r["Persister/Control"]) < P.BASELINE_MAX:
            u = float(_frac(np.array([[float(r[nm]) for nm in NAMES]]),
                            "lp_propagation")[0])
            r["Persister/RSL3"] = f"{u ** 11 * rng.uniform(0.9, 1.1):.12g}"
    d = P.assemble(doctored)
    assert d["escalation"]["admissible"]["n_flags"] > \
        d["escalation"]["full"]["n_flags"]


# --------------------------------------------------------------------------
# The Sobol design
# --------------------------------------------------------------------------
def test_the_straddling_fractions_are_recomputed_from_the_committed_flags(
        monkeypatch, raw):
    _design_patch(monkeypatch, raw)
    sob = [p[0] for p in SOBOL.PARAMS]
    d = P.assemble(raw)
    assert d["sobol_design"]["straddle"][sob[0]]["straddling_fraction"] == 1.0
    assert d["sobol_design"]["straddle"][sob[-1]]["straddling"] == 0
    moved = copy.deepcopy(raw)
    for r in moved["sobol_rows"]:
        r[f"AB_{sob[0]}"] = r["A"]
    d2 = P.assemble(moved)
    assert d2["sobol_design"]["straddle"][sob[0]]["straddling"] == 0, (
        "the straddling count is carried rather than recomputed from the flags")


def test_a_parameter_whose_swaps_never_straddle_is_reported_as_zero(
        monkeypatch, raw):
    """The zero is the load-bearing end of the range: it is what makes
    filtering a comparison between a stripped parameter and an intact one."""
    _design_patch(monkeypatch, raw)
    d = P.assemble(raw)
    assert d["sobol_design"]["straddle_min"] == 0.0
    assert d["sobol_design"]["straddle_max"] == 1.0
    md = P.render(d)
    assert "of none" in md


def test_the_excluded_parameter_s_inertness_can_fail_and_then_withholds(
        monkeypatch, raw):
    """The Sobol design drops one parameter, so the untreated arms are judged
    at its default. That is only sound if it is inert untreated."""
    _design_patch(monkeypatch, raw)
    assert P.assemble(raw)["sobol_design"]["sdt_ros_inert"] is True
    bad = copy.deepcopy(raw)
    bad["sobol_control"][0] = {"condition": "x", "low": 0.0, "high": 0.3,
                               "identical": False}
    d = P.assemble(bad)
    assert d["sobol_design"]["sdt_ros_inert"] is False
    assert "NOT inert" in P.render(d)


# --------------------------------------------------------------------------
# The resolution band
# --------------------------------------------------------------------------
def test_the_resolution_band_is_derived_from_the_committed_columns(
        monkeypatch, raw):
    _design_patch(monkeypatch, raw)
    d = P.assemble(raw)
    assert set(d["resolution"]["band"]) == set(P.RESOLUTION_CELLS)
    assert d["resolution"]["falls_with_cell_count"] is True
    # Drop one convention: the band must shrink rather than report a count it
    # no longer has.
    fewer = copy.deepcopy(raw)
    dropped = min(P.RESOLUTION_CELLS)
    for r in fewer["rows"]:
        for seed in P.RESOLUTION_SEEDS:
            r.pop(f"worst_untreated_n{dropped}_seed{seed}")
    d2 = P.assemble(fewer)
    assert dropped not in d2["resolution"]["band"]
    # Invert the direction: the page must stop claiming it falls.
    inverted = copy.deepcopy(raw)
    for r in inverted["rows"]:
        for seed in P.RESOLUTION_SEEDS:
            r[f"worst_untreated_n{dropped}_seed{seed}"] = "0.9"
    assert P.assemble(inverted)["resolution"]["falls_with_cell_count"] is False


def _reversed_keys(o):
    if isinstance(o, dict):
        return {k: _reversed_keys(v) for k, v in reversed(list(o.items()))}
    if isinstance(o, list):
        return [_reversed_keys(v) for v in o]
    return o


def test_every_published_order_is_established_by_the_renderer(monkeypatch, raw):
    """Reverse the key order of every dict `assemble` produces; the page must
    be byte-identical.

    THE GENERIC GATE CANNOT SEE THIS GENERATOR, which is why the property needs
    its own guard. `tests/test_artifact_freshness.py` shuffles the stored
    artifact and re-assembles -- but `_split_stored` re-reads the committed
    CSVs (deliberately, so no derived field is ever echoed), so the shuffle
    never reaches `assemble` and the renderer is handed the same key order
    every time. Its sibling `headline_intervals_under_admissibility` shares
    that property.

    The committed-markdown check does exercise insertion order against
    `sort_keys=True`, and a mutation reverting one list to stored dict order
    SURVIVED both: for the pair in question the two orders coincide. A
    coincidence is not a property.
    """
    _design_patch(monkeypatch, raw)
    d = P.assemble(raw)
    assert P.render(_reversed_keys(d)) == P.render(d), (
        "the page changes when the derived dicts' key order does, so at least "
        "one published order is inherited from iteration rather than "
        "established by an explicit sort in render()")


# --------------------------------------------------------------------------
# Restraint
# --------------------------------------------------------------------------
def test_the_data_conditioned_headline_count_stays_zero():
    """Conditioning removes draws and adds no data. The refusal breaks if the
    identifiability report ever says otherwise while this page stands."""
    d = json.loads((REPO / "analysis" / "identifiability-report.json").read_text())
    assert d["data_constrained_in_production"] == 0


def test_the_page_does_not_call_the_admissible_subset_a_design_or_a_posterior(
        monkeypatch, raw):
    _design_patch(monkeypatch, raw)
    md = P.render(P.assemble(raw))
    assert "not a Latin hypercube" in md
    assert "not inference" in md
    low = md.lower()
    assert "posterior" not in low or "not a posterior" in low


def test_the_page_quotes_the_manuscript_sentences_it_scores(monkeypatch, raw):
    """Quoting them is what lets a reader check the scoring; verifying they are
    still in the manuscript is what stops the quotation outliving the text."""
    body = (REPO / "article" / "drafts" / "v1.md").read_text()
    claims = P.manuscript_claims()
    live = P.assemble(P.load())
    md = P.render(live)
    for half in ("sdt", "rsl3"):
        assert claims[half]["sentence"] in body
        assert claims[half]["sentence"] in md


def test_the_committed_page_reports_the_finding_it_exists_for():
    """Read from the committed artifact deliberately: these are properties of
    the REAL measurement, not of a fixture, and they are the reason the page
    exists. Everything mechanical about the derivation is driven above."""
    d = json.loads((REPO / "analysis" / "prcc-under-admissibility.json").read_text())
    assert d["reproduction"]["reproduces"] is True
    adm = d["admissibility"]
    assert adm["published_control_mean"] > d["design"]["baseline_max"], (
        "the published Control row no longer exceeds the constraint, which is "
        "this page's premise")
    assert 0.0 < adm["fraction"] < 1.0
    assert d["circularity_clear"] is True, (
        "a drug parameter is now selected on, so the comparison this page "
        "makes has become circular and the prose must change")
    assert d["sobol_design"]["sdt_ros_inert"] is True
