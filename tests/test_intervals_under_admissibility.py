"""Guards for the admissibility-conditioned intervals (#936).

The load-bearing property is not the conditioned numbers -- it is that the
FULL-PRIOR half reproduces the committed report. Without that, the page compares
a correction against a re-run rather than against what is published, and every
claim it scores is scored against the wrong baseline. Everything else protects
the derivation: one sample partitioned rather than two sampled, each published
claim scored on the property it actually rests on, and a threshold parsed from
the sentence that states it.
"""

from __future__ import annotations

import ast
import copy
import csv
import gzip
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load("intervals_under_admissibility")
SCRIPT = SCRIPTS / "intervals_under_admissibility.py"


@pytest.fixture(scope="module")
def raw():
    return M.load()


@pytest.fixture(scope="module")
def committed():
    return json.loads(M.OUT_JSON.read_text())


def _doctor(raw):
    return {"rows": [dict(r) for r in raw["rows"]], "published": raw["published"]}


# ---------------------------------------------------------------- freshness

def test_the_committed_json_is_what_the_generator_produces(raw, committed):
    assert M.assemble(raw) == committed


def test_the_committed_markdown_is_what_the_renderer_produces(committed):
    assert M.OUT_MD.read_text() == M.render(committed)


def test_the_render_only_path_reassembles():
    tree = ast.parse(SCRIPT.read_text())
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    branch = next(n for n in ast.walk(main)
                  if isinstance(n, ast.If) and isinstance(n.test, ast.Attribute)
                  and n.test.attr == "render_only")
    assert any(isinstance(c, ast.Call) and getattr(c.func, "id", None) == "assemble"
               for stmt in branch.body for c in ast.walk(stmt))


# ------------------------------------------- the published design is one source

def test_the_design_constants_are_imported_not_restated():
    """Restating any of these is how a reproduction quietly stops reproducing."""
    src = SCRIPT.read_text()
    assert "from uncertainty_intervals import" in src
    for name in ("N_SAMPLES", "RNG_SEED", "SIM_SEED", "N_CELLS", "QUANTILES",
                 "CONDITIONS", "PARAMS"):
        assert not re.search(rf"^{name}\s*=", src, re.M), (
            f"{name} is redefined here instead of imported")
    import uncertainty_intervals as U
    for name in ("N_SAMPLES", "RNG_SEED", "SIM_SEED", "N_CELLS"):
        assert getattr(M, name) == getattr(U, name), name
    assert list(M.CONDITIONS) == list(U.CONDITIONS)
    assert list(M.QUANTILES) == list(U.QUANTILES)


def test_the_constraint_has_one_python_source():
    src = SCRIPT.read_text()
    assert "from headline_at_fitted import BASELINE_MAX" in src
    assert not re.search(r"^BASELINE_MAX\s*=", src, re.M)
    from headline_at_fitted import BASELINE_MAX
    assert M.BASELINE_MAX == BASELINE_MAX


def test_the_committed_draws_are_the_published_designs_own_draws(raw):
    """The strongest check that this is a reproduction: the parameter columns
    must equal the matrix the published generator's seed produces."""
    expected = M.draw_matrix()
    names = [p[0] for p in M.PARAMS]
    assert len(raw["rows"]) == expected.shape[0] == M.N_SAMPLES
    for i in (0, 1, len(raw["rows"]) // 2, len(raw["rows"]) - 1):
        for j, name in enumerate(names):
            assert float(raw["rows"][i][name]) == pytest.approx(
                expected[i][j], rel=1e-9), (i, name)


def test_the_draws_file_carries_every_condition(raw):
    header = set(raw["rows"][0])
    for condition in list(M.CONDITIONS) + list(M.CONTROL_CONDITIONS):
        assert M._label(condition) in header, condition
    for name, _, _ in M.PARAMS:
        assert name in header, name


# -------------------------------------------------- reproduction is load-bearing

def test_the_full_prior_half_reproduces_the_published_table(committed):
    assert committed["reproduces_published"] is True
    treated = committed["design"]["treated_conditions"]
    assert len(committed["conditions_reproduced"]) == len(treated)
    for label in treated:
        c = committed["conditions"][label]
        assert c["published"] is not None, label
        assert c["reproduces_published"] is True, label
        for key in ("median", "low", "high"):
            assert abs(c["full_prior"][key] - c["published"][key]) \
                <= M.REPRODUCTION_TOLERANCE, (label, key)


def test_a_doctored_draws_file_fails_reproduction(raw):
    """Non-vacuity: the check must be able to notice that the numbers differ."""
    d = _doctor(raw)
    label = M._label(M.CONDITIONS[0])
    for row in d["rows"]:
        row[label] = "0.000000"
    after = M.assemble(d)
    assert after["conditions"][label]["reproduces_published"] is False
    assert after["reproduces_published"] is False
    assert after["verdict"] == \
        "cannot_compare_the_published_numbers_were_not_reproduced"


def test_the_reproduction_tolerance_is_the_published_tables_own_rounding(committed):
    """DERIVED from the table's printing precision, not picked.

    A hand-picked 0.011 was 22x looser than the largest real deviation, which
    would have accepted a genuine one-point discrepancy while reporting a
    reproduction. It must also be LARGE enough that the real deviations pass, or
    the page would refuse to compare against a table it does reproduce.
    """
    expected = M.ROUNDING_SLACK * 0.5 * 10 ** (-M.PUBLISHED_DECIMALS) / 100.0
    assert M.REPRODUCTION_TOLERANCE == pytest.approx(expected)
    assert 0.0 < M.REPRODUCTION_TOLERANCE <= 0.002
    worst = 0.0
    for label in committed["design"]["treated_conditions"]:
        c = committed["conditions"][label]
        for key in ("median", "low", "high"):
            worst = max(worst, abs(c["full_prior"][key] - c["published"][key]))
    assert worst <= M.REPRODUCTION_TOLERANCE, (
        f"largest real deviation {worst} exceeds the derived tolerance")
    # And the tolerance must not be wildly larger than what it admits, which is
    # the defect it replaced.
    assert M.REPRODUCTION_TOLERANCE <= 4.0 * max(worst, 1e-9), (
        "the tolerance is far looser than the agreement it certifies")


def test_the_verdict_is_withheld_when_reproduction_fails(raw):
    d = _doctor(raw)
    for row in d["rows"]:
        for condition in M.CONDITIONS:
            row[M._label(condition)] = "0.500000"
    after = M.assemble(d)
    assert after["reproduces_published"] is False
    assert "cannot_compare" in after["verdict"]


# ------------------------------------------ one sample partitioned, not two drawn

def test_the_partition_uses_the_same_draws_for_both_halves(raw, committed):
    """Drawing a second sample would confound conditioning with sampling noise.
    Structural: both halves must be selected from one row list."""
    src = SCRIPT.read_text()
    tree = ast.parse(src)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "assemble")
    body = ast.get_source_segment(src, fn)
    assert "zip(rows, admissible)" in body, (
        "the admissible half is not selected from the same rows")
    assert "default_rng" not in body, "assemble draws its own sample"
    # And the counts must be consistent with a partition of one sample.
    adm = committed["admissibility"]
    for label in committed["design"]["treated_conditions"]:
        c = committed["conditions"][label]
        assert c["full_prior"]["n"] == adm["n_draws"], label
        assert c["admissible_only"]["n"] == adm["n_admissible"], label


def test_the_admissible_fraction_is_derived(raw, committed):
    assert committed["admissibility"]["fraction"] == pytest.approx(
        committed["admissibility"]["n_admissible"] / committed["admissibility"]["n_draws"])
    d = _doctor(raw)
    for row in d["rows"]:
        for condition in M.CONTROL_CONDITIONS:
            row[M._label(condition)] = "0.900000"
    after = M.assemble(d)
    assert after["admissibility"]["n_admissible"] == 0
    assert after["admissibility"]["reportable"] is False
    # The point of doctoring: a hard-coded fraction that happens to equal the
    # real one passes every check above. Here it cannot.
    assert after["admissibility"]["fraction"] == pytest.approx(0.0), (
        "the admissible fraction did not follow the draws, so it is not derived")


def test_no_interval_is_quoted_on_too_few_admissible_draws(raw):
    original = M.MIN_ADMISSIBLE_DRAWS
    try:
        M.MIN_ADMISSIBLE_DRAWS = 10 ** 9
        after = M.assemble(raw)
    finally:
        M.MIN_ADMISSIBLE_DRAWS = original
    for label in after["design"]["treated_conditions"]:
        assert after["conditions"][label]["admissible_only"] is None, label


def test_the_admissible_floor_is_below_the_real_count(committed):
    assert 0 < M.MIN_ADMISSIBLE_DRAWS <= committed["admissibility"]["n_admissible"]


def test_quantiles_use_the_published_method(raw):
    """A hand-rolled percentile index does not reproduce `np.quantile`'s
    interpolation, which is what produced the published numbers."""
    assert "np.quantile" in SCRIPT.read_text()
    values = [0.0, 0.25, 0.5, 0.75, 1.0]
    got = M._quantiles(values)
    low, median, high = np.quantile(np.asarray(values, float), M.QUANTILES)
    assert (got["low"], got["median"], got["high"]) == pytest.approx(
        (float(low), float(median), float(high)))


# ---------------------------------------------------- the published claims

def test_every_quoted_claim_is_still_in_the_committed_report(raw):
    """Scoring a sentence nobody publishes any more is worse than not scoring
    it: the page would read as a correction to a live claim."""
    published = raw["published"]
    for fragment, label, prop in M.published_claims():
        assert fragment in published, (
            f"this fragment is no longer in uncertainty-intervals-report.md, so "
            f"the claim being scored is not the one published: {fragment!r}")
        assert label in [M._label(c) for c in M.CONDITIONS], label
        assert prop in ("still_the_widest", "still_wide", "still_a_floor",
                        "still_above_its_quoted_floor"), prop


def test_the_claims_live_in_committed_data_not_in_the_generator():
    """A quoted claim legitimately contains its author's figures; a generator
    may not hand-write one. Keeping the quotations in data is what makes the two
    rules compatible -- and it lets a reviewer audit the list without reading
    Python."""
    assert M.CLAIMS_CSV.exists()
    src = SCRIPT.read_text()
    for fragment, _, _ in M.published_claims():
        assert fragment not in src, (
            f"this claim is quoted in the generator source, where its figures "
            f"are indistinguishable from the generator asserting them: {fragment!r}")
    assert M.CLAIMS_CSV.stat().st_size < 8192


def test_the_self_quoted_threshold_is_parsed_from_the_claim():
    """Typing the number beside a sentence that states it is two copies, and the
    sentence is the one that gets edited."""
    assert M.quoted_floor("sits >98% across the entire prior") == pytest.approx(0.98)
    assert M.quoted_floor("above >42.5% of them") == pytest.approx(0.425)
    with pytest.raises(ValueError):
        M.quoted_floor("a claim quoting no percentage of itself")
    # And it is CALLED where the claim is scored. A substring check matched the
    # `def quoted_floor(fragment)` line, so it passed while the call site had
    # been replaced by a typed constant -- a guard satisfied by the thing it is
    # meant to verify merely existing.
    src = SCRIPT.read_text()
    tree = ast.parse(src)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "assemble")
    called = {getattr(c.func, "id", None) for c in ast.walk(fn)
              if isinstance(c, ast.Call)}
    assert "quoted_floor" in called, (
        "assemble does not call quoted_floor, so the threshold is typed rather "
        "than parsed from the claim it scores")


def test_the_claim_verdicts_are_derived(raw, committed):
    """As measured, three are refuted. The guard drives the other direction:
    make every condition narrow and the width-based claims must flip."""
    assert len(committed["claims_refuted"]) > 0
    assert len(committed["claims_survived"]) > 0
    d = _doctor(raw)
    for row in d["rows"]:
        for condition in M.CONDITIONS:
            row[M._label(condition)] = "0.500000"
        for condition in M.CONTROL_CONDITIONS:
            row[M._label(condition)] = "0.000000"
    after = M.assemble(d)
    for claim in after["published_claims"]:
        if claim["property"] == "still_wide":
            assert claim["survives"] is False, claim["claim"]


def test_the_widest_interval_claim_reports_its_full_prior_margin(raw, committed):
    """A claim true by a hair and then reversed is a different story from one
    that was never close, and quoting only the winner hides which.

    Driven as well as read: this guard first read only the committed artifact,
    so deleting the margin from the basis string tripped nothing but freshness.
    That is the third time in three PRs I have written a derived-field guard
    against the file instead of against the derivation.
    """
    for source in (committed, M.assemble(raw)):
        claim = next(c for c in source["published_claims"]
                     if c["property"] == "still_the_widest")
        assert "full prior:" in claim["basis"], claim["basis"]
        assert "margin" in claim["basis"], claim["basis"]
        assert "admissible:" in claim["basis"], claim["basis"]


def test_the_wide_bar_is_a_declared_threshold(committed):
    assert 0.0 < M.WIDE_INTERVAL < 1.0
    # It must actually discriminate: some condition on each side of it.
    widths = [committed["conditions"][c]["admissible_only"]["width"]
              for c in committed["design"]["treated_conditions"]
              if committed["conditions"][c]["admissible_only"]]
    assert any(w >= M.WIDE_INTERVAL for w in widths)
    assert any(w < M.WIDE_INTERVAL for w in widths), (
        "every interval clears the wide bar, so it separates nothing")


# ------------------------------------------------- the per-condition directions

def test_the_two_directions_are_both_reported(committed):
    """One mechanism, two opposite-looking corrections. Reporting either alone
    would flatter whichever was looked at first."""
    assert committed["upper_tails_collapsed"], "no collapsed upper tail"
    assert committed["lower_tails_extended"], "no extended lower tail"
    assert not (set(committed["upper_tails_collapsed"])
                & set(committed["lower_tails_extended"])), \
        "a condition cannot be in both"
    md = M.render(committed)
    assert "upper tails collapse" in md
    assert "lower tails extend" in md


def test_a_condition_at_a_floor_or_ceiling_is_reported_unmoved(committed):
    assert committed["conditions_unmoved"], (
        "no condition is unmoved, so the per-condition framing is untested")
    for label in committed["conditions_unmoved"]:
        c = committed["conditions"][label]
        assert c["median_moved_materially"] is False, label
        assert not c["upper_tail_collapsed"], label
        assert not c["lower_tail_extended"], label


def test_the_direction_flags_are_derived(raw):
    d = _doctor(raw)
    label = M._label(M.CONDITIONS[0])
    # Identical in both halves -> nothing can have moved.
    for row in d["rows"]:
        row[label] = "0.300000"
    after = M.assemble(d)["conditions"][label]
    assert after["median_shift"] == pytest.approx(0.0)
    assert after["median_moved_materially"] is False
    assert after["upper_tail_collapsed"] is False
    assert after["lower_tail_extended"] is False


def test_the_material_shift_bar_is_declared_and_binds(committed):
    assert 0.0 < M.MATERIAL_SHIFT < 0.5
    shifts = [abs(committed["conditions"][c]["median_shift"])
              for c in committed["design"]["treated_conditions"]
              if committed["conditions"][c].get("median_shift") is not None]
    assert any(s >= M.MATERIAL_SHIFT for s in shifts)
    assert any(s < M.MATERIAL_SHIFT for s in shifts), (
        "every median moved, so the bar separates nothing")


# ------------------------------------------------------------- the refusals

def test_the_page_refuses_to_call_conditioning_inference(committed):
    md = M.render(committed)
    assert "Conditioning is not inference" in md
    assert "not a posterior" in md
    assert "data-conditioned-headline count is" in md
    report = json.loads((REPO / "analysis" / "identifiability-report.json").read_text())
    assert report.get("data_constrained_in_production") in (0, None), (
        "the identifiability report now claims a data-conditioned headline, so "
        "this page's refusal has to be revisited rather than left standing")


def test_the_page_names_the_spatial_intervals_as_outstanding(committed):
    """Leaving them silently uncorrected beside a corrected page would imply
    they had been checked."""
    md = M.render(committed)
    assert "are NOT" in md and "corrected here" in md
    assert "outstanding rather than applied" in md
    for sibling in ("headline-uncertainty-report.md",):
        assert (REPO / "analysis" / sibling).exists(), sibling


def test_the_page_states_the_phenotype_subset(committed):
    md = M.render(committed)
    assert "upper bound on admissibility" in md
    assert len(committed["design"]["untreated_conditions"]) == \
        len(M.CONTROL_CONDITIONS)


# ------------------------------------------------------------- housekeeping

def test_a_blank_field_is_not_read_as_zero(raw):
    """Every condition column must be present and numeric for every draw; a
    blank silently read as zero would make a draw look admissible."""
    for row in raw["rows"]:
        for condition in list(M.CONDITIONS) + list(M.CONTROL_CONDITIONS):
            value = row[M._label(condition)]
            assert value not in (None, ""), condition
            assert 0.0 <= float(value) <= 1.0


def test_ordered_tables_survive_the_sorted_dump(committed):
    assert isinstance(committed["published_claims"], list)
    assert isinstance(committed["design"]["treated_conditions"], list)
    assert [c["claim"] for c in committed["published_claims"]] == \
        [f for f, _, _ in M.published_claims()]


def test_the_json_dump_options_are_stable():
    assert "json.dumps(data, indent=2, sort_keys=True)" in SCRIPT.read_text()


def test_the_draws_file_is_small_enough_to_commit():
    assert M.DRAWS_CSV_GZ.exists()
    assert M.DRAWS_CSV_GZ.stat().st_size < 2 * 1024 * 1024


def test_the_module_docstring_quotes_no_measurement():
    doc = re.sub(r"#\d+", "", M.__doc__)
    standalone = re.findall(r"(?<![A-Za-z0-9_])\d[\d.,%]*", doc)
    assert not standalone, f"the docstring carries numbers: {standalone}"


def test_the_generator_reads_no_network(committed):
    src = SCRIPT.read_text()
    for forbidden in ("urllib", "requests", "corpus/atlas", "subprocess"):
        assert forbidden not in src, forbidden
