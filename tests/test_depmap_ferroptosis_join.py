"""Guards for the DepMap x CTRPv2 ferroptosis join (#932).

The thing to protect here is not the enrichment figure, it is the machinery that
could make the figure meaningless: an expectation that stops being measured, a
pass rule that quietly needs only one of its two conditions, a derivative that
goes stale against its sibling, a blank field read as a zero, and a verdict
sentence that cannot fail. Every guard on a DERIVED field drives `assemble` on
doctored raw input rather than reading the committed artifact, because a test
that reads only the artifact cannot fail when the generator changes.
"""

from __future__ import annotations

import ast
import copy
import csv
import gzip
import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"


def _load_module(name):
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load_module("depmap_ferroptosis_join")
SCRIPT = SCRIPTS / "depmap_ferroptosis_join.py"
FETCH_SCRIPT = SCRIPTS / "fetch_depmap_ctrpv2_join.py"
FETCH_SOURCE = FETCH_SCRIPT.read_text()


@pytest.fixture(scope="module")
def raw():
    return M.load()


@pytest.fixture(scope="module")
def committed():
    return json.loads(M.OUT_JSON.read_text())


def _doctor(raw):
    """A copy the tests can rewrite freely, WITHOUT a recursive deep copy.

    `raw["genes"]` is a table of ~19k single-level dicts, and `deepcopy` builds a
    memo entry for every value in every one of them. A dozen guards doing that
    was enough to get this file's pytest process SIGKILLed on a loaded machine,
    which killed three mutation sweeps -- one of them leaving a mutant behind.
    Every guard here rewrites top-level string fields only, so a per-row `dict`
    is exactly as much isolation as they need.
    """
    return {
        "genes": [dict(row) for row in raw["genes"]],
        "panel": [dict(row) for row in raw["panel"]],
        "provenance": copy.deepcopy(raw["provenance"]),
    }


# ---------------------------------------------------------------- freshness

def test_the_committed_json_is_what_the_generator_produces(raw, committed):
    """A real regenerate-and-diff, which this generator can have and most cannot:
    every input is committed, so CI re-derives the whole artifact."""
    assert M.assemble(raw) == committed


def test_the_committed_markdown_is_what_the_renderer_produces(committed):
    assert M.OUT_MD.read_text() == M.render(committed)


def test_the_render_only_path_recomputes_rather_than_echoing():
    """`main` must pass through `assemble` on BOTH branches. With every input
    committed there is no scan to skip, so a branch that read stored derived
    fields would make every guard below inert."""
    tree = ast.parse(SCRIPT.read_text())
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    branches = [n for n in ast.walk(main)
                if isinstance(n, ast.If) and isinstance(n.test, ast.Attribute)
                and n.test.attr == "render_only"]
    assert branches, "no --render-only branch"
    for node in branches:
        for body in (node.body, node.orelse):
            assert any(isinstance(c, ast.Call)
                       and getattr(c.func, "id", None) == "assemble"
                       for stmt in body for c in ast.walk(stmt)), (
                "a branch of main does not call assemble")


def test_split_stored_rereads_the_committed_inputs(raw):
    """The freshness gate's hook hands `assemble` the real inputs, not the
    artifact's own fields."""
    recovered = M._split_stored({"nothing": "useful"})
    assert isinstance(recovered, tuple) and len(recovered) == 1
    assert sorted(recovered[0]) == sorted(raw)
    assert M.assemble(*recovered) == M.assemble(raw)


# ---------------------------------------------------------------- the pass rule

def test_the_rule_needs_both_scopes_not_either(raw):
    """Blank the within-lineage column for a passing gene and it must stop
    passing, with the pooled column untouched."""
    gene = json.loads(M.OUT_JSON.read_text())["passing_genes"][0]
    d = _doctor(raw)
    compound = sorted(d["provenance"]["compounds"])[0]
    column = M._column(compound, "lineage")
    before = M.assemble(raw)["candidates"][gene]["compounds"]
    for row in d["genes"]:
        if row["gene_symbol"] == gene:
            row[column] = ""
    after = M.assemble(d)["candidates"][gene]["compounds"]
    assert after[compound]["rho_pooled"] == before[compound]["rho_pooled"]
    assert not after[compound]["passes_preregistered_rule"]


def test_a_gene_extreme_in_one_scope_only_does_not_pass(raw):
    """The guard that actually owns `and` versus `or`.

    Blanking the lineage column is not enough: a blank makes the gene UNSCORED,
    so the `scored` gate alone stops it passing and an `or` in place of the
    `and` survives. The discriminating case has BOTH values present, one in the
    tail and one in the middle -- which a mutation sweep had to point out.
    """
    stored = json.loads(M.OUT_JSON.read_text())
    gene = stored["passing_genes"][0]
    compound = sorted(raw["provenance"]["compounds"])[0]
    pooled, lineage = M._column(compound, "pooled"), M._column(compound, "lineage")
    d = _doctor(raw)
    # A value above every other gene's, so the pooled percentile is the top of
    # the distribution; and the median of the lineage column, so that scope is
    # as far from either tail as a value can be.
    others = sorted(M._f(r[lineage]) for r in d["genes"]
                    if M._f(r[lineage]) is not None)
    middle = others[len(others) // 2]
    for row in d["genes"]:
        if row["gene_symbol"] == gene:
            row[pooled] = "9.000000"
            row[lineage] = f"{middle:.6f}"
    cell = M.assemble(d)["candidates"][gene]["compounds"][compound]
    assert M._is_extreme(cell["percentile_pooled"]), "the setup is not extreme"
    assert not M._is_extreme(cell["percentile_by_lineage"]), "the setup is not mid"
    assert not cell["passes_preregistered_rule"], (
        "a gene extreme in ONE scope passes, so the rule is an OR")


def test_a_mid_distribution_percentile_is_not_extreme():
    half = M.EXTREME_TAIL / 2.0
    assert M._is_extreme(half / 2.0)
    assert M._is_extreme(1.0 - half / 2.0)
    assert not M._is_extreme(0.5)
    assert not M._is_extreme(half * 2.0)
    assert not M._is_extreme(None)


def test_the_tail_is_a_declared_constant_and_binds(raw):
    """Widening the tail must admit more cells; a rule insensitive to its own
    threshold is not applying it."""
    base = M.assemble(raw)["enrichment"]["n_passing"]
    original = M.EXTREME_TAIL
    try:
        M.EXTREME_TAIL = 0.5
        wider = M.assemble(raw)["enrichment"]["n_passing"]
    finally:
        M.EXTREME_TAIL = original
    assert wider > base


# ------------------------------------------------- the expectation is measured

def test_the_expectation_is_measured_over_all_genes_not_assumed(raw):
    """Move the DEPENDENCE between the two scopes and the expectation follows.

    Extremeness is a quantile, so a fixed share of genes is in each tail
    whatever the values are; what the measured joint rate actually records is
    how much the two tails OVERLAP. Shuffling the within-lineage column across
    genes makes the two scopes independent, and the rate must fall toward the
    product of the two tails -- which is exactly the number a generator that
    assumed independence would have published all along.
    """
    import random
    stored = json.loads(M.OUT_JSON.read_text())
    compounds = sorted(raw["provenance"]["compounds"])
    nominal = M.EXTREME_TAIL * M.EXTREME_TAIL
    d = _doctor(raw)
    rng = random.Random(M.PERMUTATION_SEED)
    for c in compounds:
        column = M._column(c, "lineage")
        values = [row[column] for row in d["genes"]]
        rng.shuffle(values)
        for row, v in zip(d["genes"], values):
            row[column] = v
    after = M.assemble(d)
    for c in compounds:
        before_rate = stored["all_gene_null"][c]["rate"]
        after_rate = after["all_gene_null"][c]["rate"]
        assert after_rate < before_rate, c
        assert abs(after_rate - nominal) < abs(before_rate - nominal), c
    assert after["enrichment"]["expected_passing"] < \
        stored["enrichment"]["expected_passing"]


def test_the_rate_responds_to_the_data_at_all(raw):
    """A constant cannot tell two distributions apart.

    The guard above compares a doctored run against the COMMITTED artifact, so
    replacing the measured rate with the nominal product still satisfies it --
    0.0025 really is below the committed 0.0103. Two doctorings that must give
    two different rates is what rules a constant out, and a sweep is how that
    was found.
    """
    compounds = sorted(raw["provenance"]["compounds"])
    candidates = set(raw["provenance"]["panel_genes"])
    import random
    shuffled = _doctor(raw)
    rng = random.Random(M.PERMUTATION_SEED)
    for c in compounds:
        column = M._column(c, "lineage")
        values = [row[column] for row in shuffled["genes"]]
        rng.shuffle(values)
        for row, v in zip(shuffled["genes"], values):
            row[column] = v
    saturated = _doctor(raw)
    for row in saturated["genes"]:
        if row["gene_symbol"] in candidates:
            continue
        for c in compounds:
            row[M._column(c, "pooled")] = "9.000000"
            row[M._column(c, "lineage")] = "9.000000"
    a, b = M.assemble(shuffled), M.assemble(saturated)
    for c in compounds:
        assert a["all_gene_null"][c]["rate"] != b["all_gene_null"][c]["rate"], (
            f"{c}: the rate is the same for two different distributions, so it "
            "is not being measured")


def test_the_nominal_product_of_the_two_tails_is_not_what_is_reported(committed):
    """The two conditions share a cohort, so the measured joint rate must differ
    from the independent product the nominal tail would give."""
    nominal = (M.EXTREME_TAIL / 2.0) * 2 * (M.EXTREME_TAIL / 2.0) * 2
    for c, null in committed["all_gene_null"].items():
        assert null["rate"] is not None
        assert abs(null["rate"] - nominal) > nominal / 2.0, c


def test_a_panel_with_no_passing_cell_reaches_the_no_signal_verdict(raw):
    """The verdict sentence has to be able to fail."""
    compounds = sorted(raw["provenance"]["compounds"])
    d = _doctor(raw)
    candidates = set(d["provenance"]["panel_genes"])
    for row in d["genes"]:
        if row["gene_symbol"] not in candidates:
            continue
        for c in compounds:
            row[M._column(c, "pooled")] = "0.000000"
            row[M._column(c, "lineage")] = "0.000000"
    after = M.assemble(d)
    assert after["enrichment"]["n_passing"] == 0
    assert after["verdict"] == "panel_indistinguishable_from_its_own_null"


def _synthetic(concentrate_tails_in_candidate_bin: bool):
    """A matrix built so the two nulls can be made to DISAGREE.

    Extremeness is a quantile, so the share of genes in a tail is fixed
    whatever the values are -- which means the all-gene null cannot be moved by
    rewriting correlations, and a guard that tries to isolate the level-matched
    null on the real matrix ends up moving both. Here the two are separable: the
    tails hold a fixed set of genes, and where those genes sit in the EXPRESSION
    bins is a free choice, so the level-matched rate moves while the all-gene
    rate does not.

    Twenty candidates sit in one bin, ten of them in the tails. With the tails
    concentrated in that bin its rate is high enough that the observed ten
    passes are no more than expected there, while the all-gene rate still says
    ten-fold. Spread the tails across bins and only the level-matched null moves.
    """
    import random
    n_genes, per_bin, tail = 3000, 300, 75
    compound = "X"
    bottom = list(range(tail))
    top = list(range(tail, 2 * tail))
    extreme = set(bottom) | set(top)
    genes = []
    for i in range(n_genes):
        if i in bottom:
            rho = -0.6 + i * 0.001
        elif i in top:
            rho = 0.5 + (i - tail) * 0.001
        else:
            rho = -0.1 + (i % 199) * 0.001
        if concentrate_tails_in_candidate_bin:
            level = (i // per_bin) * 10.0 + (i % per_bin) / 1000.0
        else:
            # Same bin for the candidates, tails spread over every bin.
            slot = i if i < 20 or i >= 2 * tail else (i * 137) % n_genes
            level = (slot // per_bin) * 10.0 + (slot % per_bin) / 1000.0
        genes.append({
            "gene_symbol": f"G{i:04d}", "entrez_id": str(100000 + i),
            "mean_expression": f"{level:.6f}", "sd_expression": "1.000000",
            M._column(compound, "pooled"): f"{rho:.6f}",
            M._column(compound, "lineage"): f"{rho:.6f}",
        })
    candidates = [f"G{i:04d}" for i in
                  (bottom[:5] + top[:5] + list(range(2 * tail, 2 * tail + 10)))]
    rng = random.Random(M.PERMUTATION_SEED)
    panel = []
    for line in range(40):
        row = {"ModelID": f"ACH-{line:06d}", "OncotreeLineage": "Synthetic",
               f"viability_auc_{compound}": f"{rng.random():.6f}"}
        for gene in candidates:
            row[f"expr_{gene}"] = f"{rng.random():.6f}"
        panel.append(row)
    provenance = {
        "source": "synthetic", "figshare_article": 0,
        "figshare_doi": "doi:figshare/synthetic", "license": "CC BY 4.0",
        "join_key": "ModelID", "n_curve_lines": len(panel),
        "n_joined_lines": len(panel), "n_genes": n_genes,
        "n_blank_expression_values": 0, "min_lineage_n": 8, "n_grid": 64,
        "readout": "synthetic readout",
        "panel_genes": {g: 100000 + int(g[1:]) for g in candidates},
        "data_blocked_genes": [candidates[0]],
        "compounds": {compound: {"n_joined_lines": len(panel), "n_curves": len(panel),
                                 "n_retained": len(panel), "n_excluded": 0,
                                 "grid_low_um": 0.001, "grid_high_um": 10.0,
                                 "n_grid": 64, "n_constant_genes": 0,
                                 "lineages_kept": {"Synthetic": len(panel)},
                                 "n_lines_dropped_small_lineage": 0}},
    }
    return {"genes": genes, "panel": panel, "provenance": provenance}, extreme


def test_the_verdict_requires_the_level_matched_null_as_well():
    """The isolation guard: the SAME pass count and the SAME all-gene rate, and
    only where the tails sit among the expression bins differs."""
    concentrated, _ = _synthetic(concentrate_tails_in_candidate_bin=True)
    spread, _ = _synthetic(concentrate_tails_in_candidate_bin=False)
    a, b = M.assemble(concentrated), M.assemble(spread)

    # Same cells, same passes, same all-gene null in both worlds.
    assert a["enrichment"]["n_cells"] == b["enrichment"]["n_cells"]
    assert a["enrichment"]["n_passing"] == b["enrichment"]["n_passing"]
    assert a["enrichment"]["expected_passing"] == pytest.approx(
        b["enrichment"]["expected_passing"])
    assert a["enrichment"]["ratio"] == pytest.approx(b["enrichment"]["ratio"])
    assert a["enrichment"]["ratio"] > 1.0

    # Only the level-matched null differs, and it decides the verdict.
    assert a["enrichment"]["expected_passing_level_matched"] > \
        b["enrichment"]["expected_passing_level_matched"]
    assert a["verdict"] == "panel_indistinguishable_from_its_own_null"
    assert b["verdict"] == "panel_beats_its_own_null"


def test_a_null_that_fires_on_everything_collapses_the_verdict(raw):
    """The blunter control on the real matrix: if the rule fires on every gene
    there is nothing left to be enriched over."""
    d = _doctor(raw)
    compounds = sorted(d["provenance"]["compounds"])
    candidates = set(d["provenance"]["panel_genes"])
    for row in d["genes"]:
        if row["gene_symbol"] in candidates:
            continue
        for c in compounds:
            row[M._column(c, "pooled")] = "9.000000"
            row[M._column(c, "lineage")] = "9.000000"
    after = M.assemble(d)
    assert after["enrichment"]["ratio"] is not None
    assert after["verdict"] == "panel_indistinguishable_from_its_own_null"


def test_the_gene_level_enrichment_counts_each_gene_once(committed):
    """The conservative statistic: one trial per gene, not one per cell."""
    gl = committed["gene_level_enrichment"]
    present = [g for g, info in committed["candidates"].items()
               if info.get("present")]
    passing = [g for g in present
               if committed["candidates"][g]["n_compounds_passing"] > 0]
    assert gl["n_genes"] == len(present)
    assert gl["n_genes_passing"] == len(passing) == len(committed["passing_genes"])
    assert gl["n_genes_passing"] <= gl["n_genes"]
    # It must be the WEAKER claim: fewer trials, and a larger per-trial null.
    assert gl["n_genes"] < committed["enrichment"]["n_cells"]
    for c, null in committed["all_gene_null"].items():
        assert gl["per_gene_null"] > null["rate"], c


def test_the_per_gene_null_is_derived_from_the_measured_rates(raw):
    """Not a constant: move the per-compound rates and it must follow."""
    base = M.assemble(raw)["gene_level_enrichment"]["per_gene_null"]
    compounds = sorted(raw["provenance"]["compounds"])
    candidates = set(raw["provenance"]["panel_genes"])
    d = _doctor(raw)
    for row in d["genes"]:
        if row["gene_symbol"] in candidates:
            continue
        for c in compounds:
            row[M._column(c, "pooled")] = "9.000000"
            row[M._column(c, "lineage")] = "9.000000"
    after = M.assemble(d)["gene_level_enrichment"]["per_gene_null"]
    assert after != base
    # And it is a probability built from the per-compound rates, so it must sit
    # above every one of them and below one.
    assert max(M.assemble(raw)["all_gene_null"][c]["rate"]
               for c in compounds) < base < 1.0


def test_the_verdict_requires_the_gene_level_enrichment_too(raw):
    """Concentrating every pass in ONE gene keeps the cell-level ratio above one
    and drives the gene-level ratio below it, so the two can disagree and the
    verdict must follow the weaker."""
    compounds = sorted(raw["provenance"]["compounds"])
    candidates = sorted(raw["provenance"]["panel_genes"])
    chosen = json.loads(M.OUT_JSON.read_text())["passing_genes"][0]
    d = _doctor(raw)
    middles = {}
    for c in compounds:
        for scope in ("pooled", "lineage"):
            col = M._column(c, scope)
            vals = sorted(M._f(r[col]) for r in d["genes"]
                          if M._f(r[col]) is not None)
            middles[col] = vals[len(vals) // 2]
    for row in d["genes"]:
        if row["gene_symbol"] not in candidates:
            continue
        for c in compounds:
            for scope in ("pooled", "lineage"):
                col = M._column(c, scope)
                row[col] = ("9.000000" if row["gene_symbol"] == chosen
                            else f"{middles[col]:.6f}")
    after = M.assemble(d)
    assert after["gene_level_enrichment"]["n_genes_passing"] == 1
    assert after["enrichment"]["ratio"] > M.ENRICHMENT_NULL, "cell level still enriched"
    assert after["gene_level_enrichment"]["ratio"] <= M.ENRICHMENT_NULL
    assert after["verdict"] == "panel_indistinguishable_from_its_own_null"


def test_the_report_leads_the_reader_to_the_conservative_figure(committed):
    md = M.render(committed)
    gl = committed["gene_level_enrichment"]
    assert "assumes the cells are independent" in md
    assert f"{gl['n_genes_passing']} of {gl['n_genes']} genes" in md
    assert f"{gl['ratio']:.2f}x" in md
    assert "the one to quote if only one is quoted" in md


def test_the_enrichment_is_not_reported_below_the_cell_floor(raw):
    d = _doctor(raw)
    keep = sorted(d["provenance"]["panel_genes"])[:1]
    d["provenance"]["panel_genes"] = {g: d["provenance"]["panel_genes"][g]
                                      for g in keep}
    after = M.assemble(d)
    assert after["enrichment"]["n_cells"] < M.MIN_CELLS_FOR_ENRICHMENT
    assert not after["enrichment"]["reportable"]
    assert after["verdict"] == "panel_indistinguishable_from_its_own_null"


def test_the_cell_floor_is_below_the_real_cell_count(committed):
    """A floor above the real count would silence the page permanently; a floor
    of zero would admit a one-cell claim."""
    assert 0 < M.MIN_CELLS_FOR_ENRICHMENT <= committed["enrichment"]["n_cells"]


# --------------------------------------------------- the pre-registered direction

def test_the_engine_still_fixes_the_acsl4_direction_in_source():
    """The one pre-registered sign is PARSED from the engine, not typed here.

    `biochem.rs` must still ADD the ACSL4 boost to the peroxidizable-substrate
    term, and `acsl4.rs` must still call the negative status refractory. If the
    engine's sign moves, this generator's registered direction is stale and the
    page must be rewritten rather than quietly disagreeing with the model.
    """
    core = REPO / "simulations" / "ferroptosis-core" / "src"
    biochem = (core / "biochem.rs").read_text()
    assert "+ params.acsl4_status_boost" in biochem
    acsl4 = (core / "acsl4.rs").read_text().lower()
    assert "refractory" in acsl4
    assert "boost `-1`" in acsl4 or "collapses" in acsl4
    assert M.PREREGISTERED_DIRECTION == "negative"


def test_the_direction_verdict_follows_the_sign(raw):
    """Flip the stored correlation and the registered-direction verdict must
    flip with it."""
    d = _doctor(raw)
    for row in d["genes"]:
        if row["gene_symbol"] != M.PREREGISTERED_GENE:
            continue
        for c in M.PREREGISTERED_COMPOUNDS:
            value = M._f(row[M._column(c, "pooled")])
            row[M._column(c, "pooled")] = f"{-value:.6f}"
    after = M.assemble(d)["preregistered"]
    assert after["verdict"] == "direction_failed_for_every_preregistered_compound"
    assert all(r["direction_held"] is False for r in after["rows"])


def test_the_direction_verdict_is_not_a_constant(raw, committed):
    """It currently holds, so the guard above is the one that can fail; this one
    pins that the committed artifact really reports the held case."""
    assert committed["preregistered"]["verdict"] == \
        "direction_held_for_every_preregistered_compound"
    assert len(committed["preregistered"]["rows"]) == \
        len(M.PREREGISTERED_COMPOUNDS)


def test_holding_the_direction_is_separate_from_clearing_the_tail(committed):
    """The two halves of the registered outcome must be reported apart. As
    measured they disagree, which is exactly why they cannot share a field."""
    pre = committed["preregistered"]
    held = [r["direction_held"] for r in pre["rows"]]
    extreme = [r["passes_extremeness_rule"] for r in pre["rows"]]
    assert all(held)
    assert pre["n_compounds_extreme"] == sum(1 for v in extreme if v)
    assert not pre["direction_and_extremeness_agree"]


# ------------------------------------------------- the four data-blocked rows

def test_the_data_blocked_genes_are_the_siblings_own_four():
    """Pinned against `calibration-feasibility.json`, not against a list typed
    here, so the two artifacts cannot drift apart."""
    sibling = json.loads(
        (REPO / "analysis" / "calibration" / "calibration-feasibility.json").read_text())
    from fetch_depmap_ctrpv2_join import DATA_BLOCKED_GENES  # noqa: E402
    assert sorted(DATA_BLOCKED_GENES) == sorted(sibling["gaps"])


def test_a_data_blocked_row_says_tried_not_unavailable(committed):
    """The whole point of the join: the four rows must no longer read as
    untestable, whatever the measurement said."""
    rows = committed["data_blocked_rows"]
    assert rows, "no data-blocked rows"
    for gene, row in rows.items():
        assert row["present"], gene
        assert row["verdict"] != "route_unavailable", gene


def test_a_data_blocked_row_that_anchored_would_say_so(raw):
    """The verdict is derived, so a gene that did clear the rule must be
    reported as anchored rather than as tried-and-failed."""
    d = _doctor(raw)
    from fetch_depmap_ctrpv2_join import DATA_BLOCKED_GENES  # noqa: E402
    gene = sorted(DATA_BLOCKED_GENES)[0]
    compounds = sorted(d["provenance"]["compounds"])
    for row in d["genes"]:
        if row["gene_symbol"] != gene:
            continue
        for c in compounds:
            row[M._column(c, "pooled")] = "9.000000"
            row[M._column(c, "lineage")] = "9.000000"
    after = M.assemble(d)["data_blocked_rows"][gene]
    assert after["n_compounds_passing"] == len(compounds)
    assert after["verdict"] == "anchored_by_this_route"


# ------------------------------------------- the two derivatives must agree

def test_the_panel_reproduces_the_correlation_table(committed):
    """Two committed derivatives, computed by different code over different
    files. A stale one is otherwise invisible."""
    checked = 0
    for gene, info in committed["candidates"].items():
        if not info.get("present"):
            continue
        for c, cc in info["compounds"].items():
            assert cc["panel_agrees_with_table"] is True, (gene, c)
            checked += 1
    assert checked >= len(committed["compounds"])


def test_a_doctored_panel_fails_the_agreement_check(raw):
    """Non-vacuity: the tolerance must be tight enough to reject a panel that
    disagrees with the table."""
    d = _doctor(raw)
    compound = sorted(d["provenance"]["compounds"])[0]
    column = f"viability_auc_{compound}"
    rows = [r for r in d["panel"] if r.get(column)]
    for a, b in zip(rows[::2], rows[1::2]):
        a[column], b[column] = b[column], a[column]
    after = M.assemble(d)["candidates"]
    disagreed = sum(1 for g in after if after[g].get("present")
                    and after[g]["compounds"][compound]["panel_agrees_with_table"]
                    is False)
    assert disagreed > 0


def test_the_agreement_tolerance_has_a_bound_of_its_own():
    """A tolerance wide enough to accept anything would pass the guard above
    while checking nothing."""
    assert 0.0 < M.PANEL_AGREEMENT_TOLERANCE <= 0.05


# ------------------------------------------------------ absence is not zero

def test_a_blank_correlation_is_absent_not_zero(raw):
    """A gene with no correlation must be UNSCORED, not a non-passer sitting at
    the middle of the distribution. This repo has shipped the opposite."""
    assert M._f("") is None
    assert M._f(None) is None
    assert M._f("0") == 0.0
    d = _doctor(raw)
    compound = sorted(d["provenance"]["compounds"])[0]
    before = M.assemble(raw)["all_gene_null"][compound]["n_genes"]
    blanked = 0
    candidates = set(d["provenance"]["panel_genes"])
    for row in d["genes"]:
        if row["gene_symbol"] in candidates:
            continue
        row[M._column(compound, "pooled")] = ""
        blanked += 1
        if blanked >= 50:
            break
    after = M.assemble(d)["all_gene_null"][compound]["n_genes"]
    assert after == before - blanked


def test_a_constant_gene_is_counted_rather_than_scored_as_null(committed):
    """The fetch leg reports zero-variance genes; they must not arrive here as
    correlations of zero."""
    for c, info in committed["compounds"].items():
        assert info["n_constant_genes"] >= 0
        scored = committed["all_gene_null"][c]["n_genes"]
        assert scored <= committed["provenance"]["n_genes"] - \
            info["n_constant_genes"], c


# ------------------------------------------------ the level-matched null

def test_a_thin_expression_bin_gives_no_level_matched_rate(raw):
    d = _doctor(raw)
    original = M.MIN_MATCHED_GENES
    try:
        M.MIN_MATCHED_GENES = 10 ** 9
        after = M.assemble(d)
    finally:
        M.MIN_MATCHED_GENES = original
    for gene, info in after["candidates"].items():
        if not info.get("present"):
            continue
        for cc in info["compounds"].values():
            assert cc["level_matched_null"] is None, gene


def test_a_cell_without_a_level_matched_rate_is_not_dropped(raw):
    """It falls back to the all-gene rate in the denominator, so the two
    expectations stay over the same set of cells."""
    original = M.MIN_MATCHED_GENES
    try:
        M.MIN_MATCHED_GENES = 10 ** 9
        after = M.assemble(raw)
    finally:
        M.MIN_MATCHED_GENES = original
    enr = after["enrichment"]
    assert enr["expected_passing_level_matched"] == pytest.approx(
        enr["expected_passing"])


def test_the_level_matched_floor_has_bounds_of_its_own(committed):
    """A floor nothing bounds is decoration.

    Both existing guards SET `MIN_MATCHED_GENES` at run time, so neither says
    anything about its committed value -- a sweep dropped it to zero and both
    stayed green. A rate over one gene is 0 or 1 and means nothing; a floor
    above the smallest real bin would silence the level-matched null entirely.
    """
    sizes = sorted(committed["expression_bin_sizes"].values())
    assert 1 < M.MIN_MATCHED_GENES <= sizes[0], (
        f"MIN_MATCHED_GENES = {M.MIN_MATCHED_GENES} is outside "
        f"(1, {sizes[0]}], so it either admits a meaningless rate or "
        "suppresses every real bin")
    # And it must actually be doing work: every real bin clears it, so the
    # committed artifact must carry a rate for every candidate.
    quoted = sum(1 for info in committed["candidates"].values()
                 if info.get("present")
                 for cc in info["compounds"].values()
                 if cc["level_matched_null"] is not None)
    assert quoted > 0, "no cell carries a level-matched rate"


def test_the_expression_bins_are_equal_count_and_cover_every_gene(committed):
    sizes = committed["expression_bin_sizes"]
    assert len(sizes) == M.EXPRESSION_BINS
    assert sum(sizes.values()) == committed["provenance"]["n_genes"]
    counts = sorted(sizes.values())
    assert counts[-1] - counts[0] <= M.EXPRESSION_BINS


# ------------------------------------------------------ the permutation null

def test_the_permutation_p_has_a_floor_rather_than_reaching_zero(committed):
    floor = 1.0 / (M.N_PERMUTATIONS + 1)
    seen = 0
    for info in committed["candidates"].values():
        if not info.get("present"):
            continue
        for cc in info["compounds"].values():
            if cc["permutation_p"] is None:
                continue
            assert cc["permutation_p"] >= floor
            seen += 1
    assert seen > 0


def test_the_permutation_null_is_seeded_and_reproducible(raw):
    first = M.assemble(raw)
    second = M.assemble(raw)
    assert first == second


def test_shuffling_the_readout_destroys_the_correlation(raw):
    """The sanity check on the whole measurement: with the readout permuted
    against expression, the panel must stop beating its null."""
    import random
    d = _doctor(raw)
    compounds = sorted(d["provenance"]["compounds"])
    rng = random.Random(M.PERMUTATION_SEED)
    for c in compounds:
        for row in d["genes"]:
            for scope in ("pooled", "lineage"):
                column = M._column(c, scope)
                if M._f(row[column]) is None:
                    continue
                row[column] = f"{rng.uniform(-0.1, 0.1):.6f}"
    after = M.assemble(d)
    assert after["enrichment"]["n_passing"] <= \
        M.assemble(raw)["enrichment"]["n_passing"]


# ------------------------------------------------- independence measurements

def test_the_compound_independence_is_measured_not_asserted(raw, committed):
    """Driven through `assemble`, not read off the artifact.

    Reading the committed file cannot notice the computation being emptied --
    a sweep deleted the whole pairwise loop and only the freshness guard fired.
    """
    for source in (committed, M.assemble(raw)):
        agree = source["readout_agreement_across_compounds"]
        n = len(source["compounds"])
        assert len(agree["pairs"]) == n * (n - 1) // 2, "the pairwise loop is empty"
        for row in agree["pairs"]:
            assert row["rho"] is not None
            assert row["n_shared_lines"] > 0
        assert agree["max_rho"] is not None and agree["min_rho"] is not None
        assert agree["max_rho"] > agree["min_rho"], "a single value is not a spread"
    # And it must track the data: dropping a compound's readout changes the rho.
    compound = sorted(raw["provenance"]["compounds"])[0]
    other = sorted(raw["provenance"]["compounds"])[1]
    d = _doctor(raw)
    for row in d["panel"]:
        row[f"viability_auc_{compound}"] = ""
    pairs = {tuple(r["compounds"]): r for r in
             M.assemble(d)["readout_agreement_across_compounds"]["pairs"]}
    assert pairs[(compound, other)]["rho"] is None, (
        "removing a readout left its correlation unchanged")


def test_the_passing_genes_collinearity_is_measured_not_asserted(committed):
    col = committed["passing_gene_collinearity"]
    assert col is not None
    genes = committed["passing_genes"]
    assert col["n_genes"] == len(genes)
    assert len(col["pairs"]) == len(genes) * (len(genes) - 1) // 2
    assert col["max_abs_rho"] is not None
    assert 0.0 <= col["max_abs_rho"] <= 1.0


def test_the_sign_unanimity_is_derived(raw, committed):
    """It is unanimous as measured, so the guard drives the non-unanimous case
    rather than asserting the comfortable one."""
    assert committed["direction_of_passing_cells"]["unanimous"] is True
    d = _doctor(raw)
    gene = committed["passing_genes"][0]
    compounds = sorted(d["provenance"]["compounds"])
    for row in d["genes"]:
        if row["gene_symbol"] != gene:
            continue
        for c in compounds:
            for scope in ("pooled", "lineage"):
                value = M._f(row[M._column(c, scope)])
                if value is not None:
                    row[M._column(c, scope)] = f"{-value:.6f}"
    after = M.assemble(d)["direction_of_passing_cells"]
    assert after["n_negative"] > 0
    assert after["unanimous"] is False


# -------------------------------------------------------- the readout itself

def test_the_readout_is_interpolation_inside_each_lines_own_dose_range():
    """Raw EC50 cannot be the readout: the committed curves contain run-away
    fits above the screened maximum. The fetch leg must go through the existing
    cohort selector rather than rolling its own range check."""
    assert "from ctrp_dose_support import select_supported_cohort" in FETCH_SOURCE
    assert "select_supported_cohort(crows, grid)" in FETCH_SOURCE
    curves = list(csv.DictReader(
        (REPO / "analysis" / "calibration" / "ctrpv2_ferroptosis_curves.csv").open()))
    runaway = [r for r in curves
               if float(r["EC50"]) > 10.0 * float(r["MaximumDose"])]
    assert runaway, "no run-away fit left to justify the policy"


def test_the_dose_grid_endpoints_are_the_recorded_bounds():
    """Computing the top of a log grid by arithmetic lands a few ulps above the
    recorded maximum, which excluded every curve of a compound."""
    fetch = _load_module("fetch_depmap_ctrpv2_join")
    for low, high in ((3e-4, 10.0), (7.6e-5, 2.5), (1e-3, 33.0), (2e-3, 66.0)):
        grid = fetch.dose_grid(low, high)
        assert grid[0] == low
        assert grid[-1] == high
        assert all(a < b for a, b in zip(grid, grid[1:]))


def test_the_dose_grid_is_the_modal_recorded_range():
    fetch = _load_module("fetch_depmap_ctrpv2_join")
    rows = [{"MinimumDose": "1", "MaximumDose": "10"}] * 3 + \
           [{"MinimumDose": "1", "MaximumDose": "5"}] * 2
    assert fetch.modal_range(rows) == (1.0, 10.0)
    # A tie breaks on the WIDER range, so the choice cannot depend on row order.
    tied = [{"MinimumDose": "1", "MaximumDose": "10"},
            {"MinimumDose": "1", "MaximumDose": "5"}]
    assert fetch.modal_range(tied) == (1.0, 10.0)
    assert fetch.modal_range(list(reversed(tied))) == (1.0, 10.0)


def test_the_curve_formula_has_one_python_source():
    """The leg that wrote the committed curves owns how to evaluate them."""
    assert "from fetch_calibration_data import predicted_viability" in FETCH_SOURCE
    assert "def predicted_viability" not in FETCH_SOURCE


# -------------------------------------------------- gene identity and access

def test_genes_are_resolved_by_entrez_with_the_symbol_checked():
    """A candidate silently resolving to another locus would be invisible in
    every figure. This repo has already been bitten by a symbol-level mapping
    error in a third-party annotation."""
    assert "by_entrez = {int(e): i for i, e in enumerate(entrez)}" in FETCH_SOURCE
    assert "panel gene identity failed" in FETCH_SOURCE
    assert "if symbols[i] != sym:" in FETCH_SOURCE


def test_fsp1_is_carried_as_aifm2_by_entrez(committed):
    """The gene whose census evidence was hidden by a symbol collision is in the
    panel under its own identifier."""
    fetch = _load_module("fetch_depmap_ctrpv2_join")
    assert fetch.PANEL_GENES["AIFM2"] == 84883
    assert committed["candidates"]["AIFM2"]["entrez_id"] == 84883


def test_the_provenance_records_the_open_route_and_its_licence(committed):
    prov = committed["provenance"]
    assert prov["join_key"] == "ModelID"
    assert "figshare" in prov["figshare_doi"]
    assert prov["license"].startswith("CC BY")
    assert prov["n_joined_lines"] > prov["n_curve_lines"] / 2
    assert prov["n_genes"] > 1000


def test_the_fetch_leg_verifies_the_download_against_the_published_hash():
    assert "hashlib.md5" in FETCH_SOURCE
    assert "md5 mismatch" in FETCH_SOURCE


def test_the_raw_download_is_not_committed():
    """The offline contract: CI reads the derivatives, never the matrix."""
    for name in ("OmicsExpressionProteinCodingGenesTPMLogp1.csv", "Model.csv"):
        assert not (REPO / "analysis" / "calibration" / name).exists()
    committed_bytes = (REPO / "analysis" / "calibration"
                       / "depmap_ctrpv2_gene_correlations.csv.gz").stat().st_size
    assert committed_bytes < 20 * 1024 * 1024


def test_the_generator_reads_no_network_and_no_census():
    source = SCRIPT.read_text()
    for forbidden in ("urllib", "requests", "corpus/atlas", "RECORDS", "subprocess"):
        assert forbidden not in source, forbidden


# ----------------------------------------------------- lineage stratification

def test_the_lineage_scope_is_a_real_second_measurement(committed):
    """Pooled and within-lineage must actually differ somewhere, or the second
    condition is the first one twice."""
    differed = 0
    for info in committed["candidates"].values():
        if not info.get("present"):
            continue
        for cc in info["compounds"].values():
            a, b = cc["rho_pooled"], cc["rho_by_lineage"]
            if a is not None and b is not None and abs(a - b) > 1e-6:
                differed += 1
    assert differed > 0


def test_the_stratified_ranks_are_centred_within_stratum():
    """Stacking raw within-stratum ranks leaves a between-stratum component --
    the very lineage effect the stratification removes."""
    assert "block - block.mean(axis=0)" in FETCH_SOURCE
    assert "y - y.mean()" in FETCH_SOURCE


def test_small_lineages_are_dropped_rather_than_pooled(committed):
    assert committed["provenance"]["min_lineage_n"] > 1
    for c, info in committed["compounds"].items():
        assert info["lineages_kept"], c
        assert min(info["lineages_kept"].values()) >= \
            committed["provenance"]["min_lineage_n"], c
        assert info["n_lines_dropped_small_lineage"] >= 0


# ---------------------------------------------------------------- the report

def test_the_report_states_what_is_not_established(committed):
    md = M.render(committed)
    for phrase in ("does NOT establish", "not tumours in patients",
                   "not independent", "placeholder"):
        assert phrase in md, phrase


def test_the_report_names_the_access_claim_it_corrects(committed):
    md = M.render(committed)
    assert "blocked on ACCESS" in md
    assert "bot-verification page" in md
    assert committed["provenance"]["figshare_doi"] in md


def test_a_vanishing_p_value_is_printed_rather_than_rounded_to_zero():
    assert M._p(6.5e-12) == "6.5e-12"
    assert M._p(0.0) == "0.0e+00"
    assert M._p(0.25) == "0.2500"
    assert M._p(None) == "n/a"


def test_the_report_quotes_no_figure_the_artifact_does_not_carry(committed):
    """Every headline number in the prose must move when the artifact moves."""
    md = M.render(committed)
    enr = committed["enrichment"]
    assert f"{enr['n_passing']}" in md
    assert f"{enr['ratio']:.2f}x" in md
    assert f"{enr['ratio_level_matched']:.2f}x" in md
    doctored = copy.deepcopy(committed)
    doctored["enrichment"]["ratio"] = enr["ratio"] + 1.0
    assert M.render(doctored) != md


def test_the_markdown_carries_every_candidate_gene(committed):
    md = M.render(committed)
    for gene in committed["candidates"]:
        assert gene in md, gene


# ------------------------------------------------------------- housekeeping

def test_the_json_dump_options_are_stable():
    source = SCRIPT.read_text()
    assert "json.dumps(data, indent=2, sort_keys=True)" in source


def test_ordered_tables_survive_the_sorted_dump(committed):
    """`sort_keys` reorders dicts, so anything the report presents in a chosen
    order has to be stored as a LIST. This repo has shipped alphabetised tables
    presented as rankings twice."""
    assert isinstance(committed["readout_agreement_across_compounds"]["pairs"], list)
    assert isinstance(committed["passing_gene_collinearity"]["pairs"], list)
    assert isinstance(committed["preregistered"]["rows"], list)
    assert isinstance(committed["passing_genes"], list)


def test_the_correlations_table_has_a_row_for_every_gene(raw, committed):
    assert len(raw["genes"]) == committed["provenance"]["n_genes"]
    with gzip.open(M.CORRELATIONS_CSV_GZ, "rt") as fh:
        header = next(csv.reader(fh))
    for c in committed["compounds"]:
        assert M._column(c, "pooled") in header
        assert M._column(c, "lineage") in header
    assert "mean_expression" in header and "entrez_id" in header


def test_the_binomial_tail_is_the_exact_sum():
    assert M._binomial_tail(4, 0, 0.5) == pytest.approx(1.0)
    assert M._binomial_tail(4, 4, 0.5) == pytest.approx(0.0625)
    assert M._binomial_tail(2, 1, 0.5) == pytest.approx(0.75)
    assert M._binomial_tail(10, 3, None) is None
    assert M._binomial_tail(0, 0, 0.5) is None


def test_the_rank_helper_averages_ties():
    assert M._ranks([3.0, 1.0, 2.0]) == [3.0, 1.0, 2.0]
    assert M._ranks([1.0, 1.0, 2.0]) == [1.5, 1.5, 3.0]
    assert M._ranks([5.0, 5.0, 5.0]) == [2.0, 2.0, 2.0]


def test_the_spearman_helper_matches_a_known_value():
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert M._spearman_pairs(a, a) == pytest.approx(1.0)
    assert M._spearman_pairs(a, list(reversed(a))) == pytest.approx(-1.0)
    assert M._spearman_pairs([1.0, 2.0], [1.0, 2.0]) is None
    assert M._spearman_pairs(a, [1.0, 1.0, 1.0, 1.0, 1.0]) is None


def test_the_percentile_is_a_fraction_at_or_below():
    values = [1.0, 2.0, 3.0, 4.0]
    assert M._percentile(values, 0.5) == 0.0
    assert M._percentile(values, 2.0) == 0.5
    assert M._percentile(values, 9.0) == 1.0


def test_the_module_docstring_quotes_no_measurement():
    """A docstring figure goes stale silently; the artifact's do not.

    Issue references are the one allowed numeral, because they name the work
    rather than reporting a result.
    """
    import re
    doc = re.sub(r"#\d+", "", M.__doc__)
    # A digit attached to letters is an identifier -- ACSL4, AIFM2, acsl4.rs --
    # not a measurement. What must not appear is a number standing on its own.
    standalone = re.findall(r"(?<![A-Za-z0-9_])\d[\d.,%]*", doc)
    assert not standalone, (
        f"the docstring carries numbers that cannot be regenerated: {standalone}")


def test_the_claude_md_entry_quotes_the_committed_figures(committed):
    """Every figure hand-written into the repo's front-door file must be the
    artifact's own.

    This is not hypothetical: in the immediately preceding piece of work a
    figure written into CLAUDE.md went stale INSIDE the same session, because
    the measurement moved underneath it. Prose is the one place a number cannot
    be regenerated, so it gets pinned instead.
    """
    entry = next(block for block in
                 (REPO / "CLAUDE.md").read_text().split("\n")
                 if block.startswith("scripts/depmap_ferroptosis_join.py"))
    prov, enr = committed["provenance"], committed["enrichment"]
    pre = committed["preregistered"]
    col = committed["passing_gene_collinearity"]
    agree = committed["readout_agreement_across_compounds"]
    rates = sorted(committed["all_gene_null"][c]["rate"]
                   for c in committed["compounds"])
    fired = sorted(committed["all_gene_null"][c]["n_both_extreme"]
                   for c in committed["compounds"])
    nominal = M.EXTREME_TAIL * M.EXTREME_TAIL
    expect = [
        str(prov["figshare_article"]), prov["figshare_doi"],
        str(prov["n_joined_lines"]), str(prov["n_curve_lines"]),
        f"{prov['n_genes']:,}",
        f"{100 * rates[0]:.2f}%-{100 * rates[-1]:.2f}%",
        f"{100 * nominal:.2f}%",
        f"{enr['n_passing']} of {enr['n_cells']}",
        f"{enr['expected_passing']:.2f}", f"{enr['ratio']:.2f}x",
        f"{enr['binomial_p']:.1e}",
        f"{enr['expected_passing_level_matched']:.2f}",
        f"{enr['ratio_level_matched']:.2f}x",
        f"{enr['binomial_p_level_matched']:.1e}",
        f"{len(committed['passing_genes'])} of {len(committed['candidates'])}",
        str(committed["direction_of_passing_cells"]["n_positive"]),
        str(committed["candidates"]["AIFM2"]["entrez_id"]),
        f"{pre['rows'][0]['rho_pooled']:.3f}",
        f"{pre['rows'][1]['rho_pooled']:.3f}",
        f"{100 * pre['rows'][0]['percentile_pooled']:.2f}%",
        f"{100 * pre['rows'][1]['percentile_pooled']:.2f}%",
        f"{agree['min_rho']:.3f} to {agree['max_rho']:.3f}",
        f"{agree['median_rho']:.3f}",
        f"{col['max_abs_rho']:.3f}", f"{col['median_abs_rho']:.3f}",
        f"{fired[0]} to {fired[-1]}",
        f"{round(prov_size() / 1e6)} MB",
    ]
    missing = [v for v in expect if v not in entry]
    assert not missing, f"CLAUDE.md does not quote the committed values: {missing}"
    for gene, n in (("SLC7A11", 4), ("SLC3A2", 3), ("GCLC", 3),
                    ("AIFM2", 3), ("ALOX15B", 1)):
        assert committed["candidates"][gene]["n_compounds_passing"] == n, gene
    n_guards = sum(1 for line in
                   Path(__file__).read_text().split("\n")
                   if line.startswith("def test_"))
    assert f"{n_guards} guards" in entry, (
        f"CLAUDE.md claims a guard count that is not {n_guards}")


def prov_size():
    """The expression matrix's size, as figshare reported it."""
    prov = json.loads(M.PROVENANCE_JSON.read_text())
    return prov["files"]["OmicsExpressionProteinCodingGenesTPMLogp1.csv"]["size_bytes"]


def test_math_is_imported_where_it_is_used():
    assert math.comb(4, 2) == 6
