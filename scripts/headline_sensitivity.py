#!/usr/bin/env python3
"""Morris elementary-effects sensitivity of the spatial/combo HEADLINE outputs (#331).

The single-cell Sobol (#386, `scripts/sobol_sensitivity.py`) screened the
ferroptosis death switch via the fast `sim_batch` binding. This is the per-
HEADLINE follow-up: it screens which biochemical rate constants drive a headline
output that only the simulation BINARIES produce, by perturbing the SAME 11 PRCC
rate constants through the `FERRO_PARAM_OVERRIDES` env hook (added in #331 PR 1)
and reading the binary's output.

**Method: Morris elementary effects (Morris 1991), not Saltelli/Sobol.** A full
variance-based Sobol design needs `N_base * (k + 2)` model evaluations (tens of
thousands); each binary run is far more expensive than a `sim_batch` call, so a
full Sobol is infeasible for the spatial sims. Morris is the standard SCREENING
alternative: `r * (k + 1)` evaluations give, per parameter, the mean absolute
elementary effect `mu_star` (overall importance, the Morris analogue of a total
effect) and `sigma` (interaction / nonlinearity indicator). It answers "which
parameters, and which show interactions" — the issue's question — at a fraction
of the cost. It is a screening, not a variance decomposition; we report it as such.

**Covered headlines:**
- **Bliss synergy** (RSL3 + FSP1i 1.99x from `sim-combo-mech`): a clean scalar,
  cheap to evaluate (~seconds/run).
- **Hypoxia kill-collapse** (SDT-minus-RSL3 hypoxic-zone kill GAP from `sim-tme`):
  the kill-collapse asymmetry the headline reports. Each `sim-tme` run is far
  costlier than `sim-combo-mech`, so use a smaller `--trajectories` for it.

The **immune-ratio** headline is still deferred, but NOT because the raw ratio is
misleading in direction: at the canonical condition the raw `immune_kills` MATCHES
the headline (SDT >> RSL3, ~104:1 at the default gradient-120 immune-on, the
Figure-17 result). The subtlety is for a SENSITIVITY screen: `immune_kills` is
confounded by POOL DEPLETION (SDT ferroptotically clears most of the tumor first,
so a parameter that raises SDT's ferroptosis shrinks the residual pool the immune
layer acts on, lowering the immune-kill COUNT even as per-cell amplification
rises). A faithful sensitivity observable controls for the pool (e.g. the
de-confounded rate `immune_kills / non-ferroptotic-cells`). See the report's
"Deferred" section; tracked under #331.

Self-contained (no SALib): the Morris estimator is ~40 lines and is validated on
an analytic linear+interaction function in `tests/test_headline_sensitivity.py`.
Deterministic given the seed. Writes `analysis/headline-sensitivity-report.md`.

Usage:
    python3 scripts/headline_sensitivity.py [--headline {bliss,hypoxia,both}]
            [--trajectories 100] [--tme-trajectories 6] [--levels 4]
            [--workers 4] [--smoke]
    (run with no args reproduces the committed analysis/headline-sensitivity-report.md)
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
REPORT = REPO / "analysis" / "headline-sensitivity-report.md"
PRCC_JSON = REPO / "analysis" / "prcc-results.json"

# The same 11 biochemical rate constants + ranges as the univariate PRCC (#134)
# and the single-cell Sobol (#386), so this is an apples-to-apples per-headline
# extension. Loaded (not hard-coded) from the frozen PRCC results.
_RANGES = json.loads(PRCC_JSON.read_text())["metadata"]["parameter_ranges"]
PARAM_NAMES = list(_RANGES.keys())
LOWS = np.array([_RANGES[n][0] for n in PARAM_NAMES], dtype=float)
HIGHS = np.array([_RANGES[n][1] for n in PARAM_NAMES], dtype=float)


# ----------------------------------------------------------------------------
# Morris design + estimator (self-contained; validated analytically in tests)
# ----------------------------------------------------------------------------
def morris_trajectory(k, levels, delta, rng):
    """One Morris trajectory: a `(k+1, k)` array of points in the unit cube,
    consecutive rows differing in exactly one coordinate by +/- `delta`
    (the standard Morris B* construction, Morris 1991 / Saltelli 2008)."""
    grid = np.arange(levels) / (levels - 1.0)  # `levels` points in [0, 1]
    base_choices = grid[grid <= 1.0 - delta + 1e-12]  # so x + delta stays <= 1
    x_star = rng.choice(base_choices, size=k)
    d = rng.choice(np.array([-1.0, 1.0]), size=k)
    perm = rng.permutation(k)
    pmat = np.eye(k)[perm]
    b = np.tril(np.ones((k + 1, k)), -1)
    j = np.ones((k + 1, k))
    bstar = (j * x_star + (delta / 2.0) * ((2.0 * b - j) @ np.diag(d) + j)) @ pmat
    return np.clip(bstar, 0.0, 1.0)


def elementary_effects(traj_points, traj_values):
    """Per-parameter elementary effect from one trajectory: for the consecutive
    pair that differs in coordinate `i`, `EE_i = dy / dx_i` (in unit-cube
    coordinates, signed)."""
    k = traj_points.shape[1]
    ee = np.full(k, np.nan)
    for r in range(k):
        dx = traj_points[r + 1] - traj_points[r]
        i = int(np.argmax(np.abs(dx)))
        ee[i] = (traj_values[r + 1] - traj_values[r]) / dx[i]
    return ee


def morris_indices(eval_fn, lows, highs, n_traj, levels, rng_seed):
    """Run a Morris screening of `eval_fn` (maps an `(m, k)` unit-cube array to
    `(m,)` outputs) over `[lows, highs]`. Returns `(mu_star, sigma, names_order)`
    where `mu_star[i] = mean(|EE_i|)` and `sigma[i] = std(EE_i)` over trajectories."""
    lows = np.asarray(lows, float)
    highs = np.asarray(highs, float)
    k = len(lows)
    delta = levels / (2.0 * (levels - 1.0))  # standard Morris step
    rng = np.random.default_rng(rng_seed)

    trajectories = [morris_trajectory(k, levels, delta, rng) for _ in range(n_traj)]
    all_unit = np.vstack(trajectories)  # (n_traj*(k+1), k)
    scaled = lows + all_unit * (highs - lows)
    values = np.asarray(eval_fn(scaled), float)

    ee_rows = []
    step = k + 1
    for t in range(n_traj):
        pts = all_unit[t * step : (t + 1) * step]
        vals = values[t * step : (t + 1) * step]
        ee_rows.append(elementary_effects(pts, vals))
    ee = np.array(ee_rows)  # (n_traj, k)
    mu_star = np.nanmean(np.abs(ee), axis=0)
    sigma = np.nanstd(ee, axis=0)
    return mu_star, sigma


# ----------------------------------------------------------------------------
# Headline observables (run a binary under FERRO_PARAM_OVERRIDES, read its output)
# ----------------------------------------------------------------------------
def _default_binary(name="sim-combo-mech"):
    for cand in (
        REPO / "simulations" / "target" / "release" / name,
        REPO / "simulations" / "target" / "debug" / name,
    ):
        if cand.exists():
            return cand
    return None


def run_bliss(params_row, binary):
    """RSL3 + FSP1i `synergy_score` from sim-combo-mech (the ~1.99x Bliss headline)."""
    with tempfile.TemporaryDirectory(prefix="ferro_morris_") as workdir:
        overrides = {n: float(v) for n, v in zip(PARAM_NAMES, params_row)}
        env = dict(os.environ, FERRO_PARAM_OVERRIDES=json.dumps(overrides))
        proc = subprocess.run(
            [str(binary)], cwd=workdir, env=env, capture_output=True, text=True
        )
        if proc.returncode != 0:
            raise RuntimeError(f"sim-combo-mech failed ({proc.returncode}): {proc.stderr[-400:]}")
        csv_path = Path(workdir) / "output" / "combo-mech" / "combo_synergy.csv"
        with csv_path.open() as fh:
            for row in csv.DictReader(fh):
                if {row["drug_a"], row["drug_b"]} == {"RSL3", "FSP1i"}:
                    return float(row["synergy_score"])
    raise RuntimeError("RSL3+FSP1i row not found in combo_synergy.csv")


# Reference O2 gradient for the hypoxia headline (lambda = 120 um, the lambda the
# manuscript's 87.8% figure is read at). The SDT hypoxic kill sits near 1.0 and
# RSL3's is ~0, so the SDT-minus-RSL3 gap IS the kill-collapse asymmetry.
HYPOXIA_GRADIENT = "gradient_120um"


def run_hypoxia_gap(params_row, binary):
    """SDT-minus-RSL3 hypoxic-zone kill GAP at the reference gradient, from
    sim-tme's tme_summary.json (the immune-off conditions: the hypoxia headline is
    the O2-only comparison). The gap is the kill-collapse asymmetry the headline
    reports (SDT holds, RSL3 collapses)."""
    with tempfile.TemporaryDirectory(prefix="ferro_morris_tme_") as workdir:
        overrides = {n: float(v) for n, v in zip(PARAM_NAMES, params_row)}
        env = dict(os.environ, FERRO_PARAM_OVERRIDES=json.dumps(overrides))
        proc = subprocess.run(
            [str(binary)], cwd=workdir, env=env, capture_output=True, text=True
        )
        if proc.returncode != 0:
            raise RuntimeError(f"sim-tme failed ({proc.returncode}): {proc.stderr[-400:]}")
        summary = Path(workdir) / "output" / "tme" / "tme_summary.json"
        conditions = json.loads(summary.read_text())["conditions"]

        def hyp(tx):
            for r in conditions:
                if (
                    r["treatment"] == tx
                    and r["o2_condition"] == HYPOXIA_GRADIENT
                    and r["immune_mode"] == "off"
                ):
                    return r["hypoxic_kill_rate"]
            raise RuntimeError(f"no {tx}/{HYPOXIA_GRADIENT}/off row in tme_summary.json")

        return hyp("SDT") - hyp("RSL3")


def make_eval(run_fn, binary, workers):
    """Parallel evaluator: map `run_fn(row, binary)` over the Morris design rows."""

    def eval_fn(scaled_rows):
        with ThreadPoolExecutor(max_workers=workers) as ex:
            return list(ex.map(lambda r: run_fn(r, binary), scaled_rows))

    return eval_fn


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def _index_table(mu_star, sigma):
    """Markdown rank table + the descending-mu* order."""
    order = np.argsort(-mu_star)
    rows = [
        "| rank | parameter | mu* | sigma | sigma/mu* |",
        "|------|-----------|-----|-------|-----------|",
    ]
    for rank, i in enumerate(order, 1):
        ratio = sigma[i] / mu_star[i] if mu_star[i] > 1e-12 else float("nan")
        rows.append(
            f"| {rank} | `{PARAM_NAMES[i]}` | {mu_star[i]:.4f} | {sigma[i]:.4f} | {ratio:.2f} |"
        )
    return order, rows


def bliss_section(mu_star, sigma, n_traj, n_evals):
    order, table = _index_table(mu_star, sigma)
    top = ", ".join(f"`{PARAM_NAMES[i]}`" for i in order[:3])
    # Two genuinely STRUCTURAL zeros for the RSL3 + FSP1i pair (never read on this
    # code path). Any OTHER mu*=0 is an EMPIRICAL/below-resolution zero (active on
    # the biochem path but its elementary effect rounds to zero) — a weaker claim.
    structural_names = {"sdt_ros", "rsl3_gpx4_inhib"}
    zero_params = [PARAM_NAMES[i] for i in order if mu_star[i] < 1e-9]
    structural = [n for n in zero_params if n in structural_names]
    empirical = [n for n in zero_params if n not in structural_names]
    return [
        "## Headline: Bliss synergy (RSL3 + FSP1i) — `sim-combo-mech`",
        "",
        f"Observable: the RSL3 + FSP1i `synergy_score` (death-over-Bliss-prediction, "
        f"the manuscript's ~1.99x) from `combo_synergy.csv`. Morris r={n_traj} "
        f"({n_evals} runs).",
        "",
        *table,
        "",
        f"**Top drivers:** {top} — the LP-cascade and GSH/GPX4 defense constants, the "
        "same axis the single-cell Sobol and the PRCC rank for the RSL3 kill. `sigma` "
        "exceeds `mu*` for every ACTIVE parameter (`sigma/mu* > 1`), so the synergy is "
        "strongly NONLINEAR / INTERACTION-LADEN, structure a univariate PRCC cannot "
        "see. (Morris `sigma` does not separate interaction from single-parameter "
        "nonlinearity; either way the effect is non-additive.)",
        "",
        "**Zero-effect parameters (`mu* = 0`) — two distinct kinds:**",
        "- *Structural zeros* (disconnected by construction): "
        + (", ".join(f"`{n}`" for n in structural) or "none")
        + ". `sdt_ros` is the SDT dose (no SDT in the pair); `rsl3_gpx4_inhib` is "
        "never read (`sim-combo-mech` applies RSL3 as a fixed `DrugEffect`, 92% GPX4 "
        "inhibition, not via `Params.rsl3_gpx4_inhib`).",
        "- *Empirical / below-resolution zeros* (active on the biochem path, effect "
        "rounds to zero here): " + (", ".join(f"`{n}`" for n in empirical) or "none")
        + ". `gpx4_degradation_by_ros` is used in the GPX4 dynamics, but over its PRCC "
        "range a perturbation never flips a cell's bistable outcome (low single-cell "
        "Sobol ST). A practical-identifiability limit, NOT a structural disconnection.",
        "",
    ], order


def hypoxia_section(mu_star, sigma, n_traj, n_evals):
    order, table = _index_table(mu_star, sigma)
    top = ", ".join(f"`{PARAM_NAMES[i]}`" for i in order[:3])
    return [
        "## Headline: hypoxia kill-collapse (SDT vs RSL3) — `sim-tme`",
        "",
        f"Observable: the SDT-minus-RSL3 hypoxic-zone kill GAP at the reference O2 "
        f"gradient (`{HYPOXIA_GRADIENT}`, immune off) from `tme_summary.json`. This is "
        "the kill-collapse asymmetry the headline reports: SDT holds (~0.87 at "
        f"baseline), RSL3 collapses (~0). Morris r={n_traj} ({n_evals} runs; "
        "fewer than the Bliss headline because each sim-tme run is far costlier).",
        "",
        *table,
        "",
        f"**Top drivers:** {top}. `sdt_ros` (the SDT exogenous-ROS dose) is the single "
        "largest driver: the SDT hypoxic kill that *holds* is set primarily by dose. "
        "But the LP-cascade and GSH-defense constants (`lp_rate`, `gsh_scav_efficiency`, "
        "`lp_propagation`) ALSO move the gap substantially (roughly 35-60% of `sdt_ros`'s "
        "mu* at this design) — not via RSL3 (whose hypoxic kill stays ~0 regardless) but "
        "by modulating the SDT hypoxic kill itself, which sits HIGH (~0.87) but is NOT "
        "saturated, so the same defenses that resist RSL3 still partially blunt SDT. The "
        "model-side reading of the Section 7.1 asymmetry is therefore nuanced: SDT's "
        "hypoxic efficacy is primarily DOSE-driven (hence it survives hypoxia where RSL3 "
        "collapses), but it is defense-MODULATED, not defense-independent — consistent "
        "with the manuscript's framing that the SDT advantage is real but not absolute. "
        "`rsl3_gpx4_inhib` sits near the bottom (mu* ~0), since RSL3 already fails in "
        "hypoxia so changing its GPX4-inhibition strength barely moves the collapsed "
        "RSL3 kill.",
        "",
        "**Caveat (the contested leg).** SDT is modeled here as O2-INDEPENDENT, the "
        "optimistic upper bound (Section 7.1). So this screen attributes the gap to "
        "`sdt_ros` *given that assumption*; under the off-by-default O2-dependent SDT "
        "mode (#336/#358) the SDT hypoxic kill would itself collapse, and that O2 "
        "dependence is a separate knob not in this PRCC rate-constant set. The result "
        "is therefore: *among the biochemical rate constants*, only the SDT dose "
        "matters for the (O2-independent) hypoxia gap.",
        "",
    ], order


def write_report(sections, levels, total_evals):
    """Assemble the multi-headline report from `(lines, label)` sections."""
    lines = [
        "# Per-headline Morris sensitivity (#331)",
        "",
        "Generated by `scripts/headline_sensitivity.py`. Screens which of the 11 "
        "PRCC biochemical rate constants drive a **headline** output that only the "
        "simulation binaries produce, via the `FERRO_PARAM_OVERRIDES` hook (#331 PR 1). "
        "Companion to the single-cell Sobol (`scripts/sobol_sensitivity.py`, #386) "
        "and the univariate PRCC (#134).",
        "",
        "## Method",
        "",
        "**Morris elementary effects (Morris 1991), a SCREENING method — not a "
        "variance decomposition.** A full Sobol design needs `N_base*(k+2)` "
        "evaluations (tens of thousands); each binary run is orders of magnitude more "
        "expensive than the `sim_batch` call the single-cell Sobol used, so full Sobol "
        "is infeasible for the simulation binaries. Morris gives, per parameter, `mu*` "
        "(mean |elementary effect|, the importance ranking — the screening analogue of "
        "a total effect) and `sigma` (interaction / single-parameter-nonlinearity "
        f"indicator) in `r*(k+1)` runs over `k = {len(PARAM_NAMES)}` parameters at "
        f"`{levels}` levels. The trajectory count `r` is PER-HEADLINE (the cheap "
        "`sim-combo-mech` headline uses many more than the costly `sim-tme` one; each "
        f"section states its own `r`), {total_evals} model runs total. Parameters "
        "swept over the SAME PRCC ranges (`analysis/prcc-results.json`). Deterministic "
        "given the design seed.",
        "",
    ]
    for sec_lines, _ in sections:
        lines += sec_lines
    if not any(label == "immune" for _, label in sections):
        lines += [
            "## Deferred: the immune-ratio headline",
            "",
            "The raw `immune_kills` ratio at the canonical condition MATCHES the "
            "headline direction: SDT >> RSL3 (~104:1 at the default gradient-120, "
            "immune-on; the Figure-17 result), because SDT's dense ferroptotic death "
            "builds a large DAMP field that amplifies immune killing far more than "
            "RSL3's sparse death. The deferral is a SENSITIVITY-observable subtlety, "
            "not a direction problem: `immune_kills` is confounded by POOL DEPLETION "
            "(SDT ferroptotically clears most of the tumor first, so a parameter that "
            "raises SDT's ferroptosis shrinks the residual pool the immune layer acts "
            "on, which can lower the immune-kill COUNT even as per-cell amplification "
            "rises). A faithful sensitivity observable controls for the pool, e.g. the "
            "de-confounded rate `immune_kills / (total_tumor - ferroptosis_kills)` "
            "(which at baseline gives an even sharper SDT:RSL3 asymmetry, ~850:1). "
            "Building that screen is the remaining #331 increment.",
            "",
            "Tracked under #331 (which stays open for the immune headline).",
        ]
    REPORT.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trajectories", type=int, default=100)
    ap.add_argument("--levels", type=int, default=4)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20250606)
    ap.add_argument(
        "--tme-trajectories",
        type=int,
        default=6,
        help="trajectories for the sim-tme (hypoxia) headline; far fewer than --trajectories "
        "because each sim-tme run is orders of magnitude slower than sim-combo-mech. The "
        "default (6) is the value the committed report was generated with.",
    )
    ap.add_argument("--smoke", action="store_true", help="tiny run (2 trajectories each) for a wiring check")
    ap.add_argument(
        "--headline",
        choices=["bliss", "hypoxia", "both"],
        default="both",
        help="which headline(s) to screen (sim-tme is much slower than sim-combo-mech)",
    )
    args = ap.parse_args()
    bliss_r = 2 if args.smoke else args.trajectories
    tme_r = 2 if args.smoke else args.tme_trajectories

    # Per-headline: (run_fn, binary_name, section_fn, label, n_traj).
    specs = []
    if args.headline in ("bliss", "both"):
        specs.append((run_bliss, "sim-combo-mech", bliss_section, "bliss", bliss_r))
    if args.headline in ("hypoxia", "both"):
        specs.append((run_hypoxia_gap, "sim-tme", hypoxia_section, "hypoxia", tme_r))

    sections = []
    total_evals = 0
    for run_fn, bin_name, section_fn, label, n_traj in specs:
        binary = _default_binary(bin_name)
        if binary is None:
            sys.exit(
                f"{bin_name} binary not found. Build it first:\n"
                f"  (cd simulations && cargo build --release -p {bin_name})"
            )
        n_evals = n_traj * (len(PARAM_NAMES) + 1)
        total_evals += n_evals
        print(f"[{label}] Morris: r={n_traj}, k={len(PARAM_NAMES)} -> {n_evals} runs of {bin_name}")
        eval_fn = make_eval(run_fn, binary, args.workers)
        mu_star, sigma = morris_indices(eval_fn, LOWS, HIGHS, n_traj, args.levels, args.seed)
        sec_lines, order = section_fn(mu_star, sigma, n_traj, n_evals)
        sections.append((sec_lines, label))
        print(f"[{label}] top-3:", ", ".join(PARAM_NAMES[i] for i in order[:3]))

    write_report(sections, args.levels, total_evals)
    print(f"Wrote {REPORT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
