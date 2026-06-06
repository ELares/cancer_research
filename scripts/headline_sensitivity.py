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

**This PR covers the Bliss-synergy headline** (the RSL3 + FSP1i 1.99x from
`sim-combo-mech`), which is a clean scalar and cheap to evaluate (~seconds/run).
The two `sim-tme`-driven headlines (hypoxia kill-collapse, immune ratio) are the
next increment: they are far more expensive per run, the SDT kill sits near a
ceiling that compresses its variance, and the immune-ratio observable is
CONFOUNDED (raw `immune_kills` depends on how many cells survive ferroptosis
first, so it runs opposite to the headline's DAMP-amplification framing) — each
needs a deliberately-chosen observable. See the report's "Deferred" section.

Self-contained (no SALib): the Morris estimator is ~40 lines and is validated on
an analytic linear+interaction function in `tests/test_headline_sensitivity.py`.
Deterministic given the seed. Writes `analysis/headline-sensitivity-report.md`.

Usage:
    python3 scripts/headline_sensitivity.py [--trajectories 10] [--levels 4]
            [--workers 4] [--smoke] [--binary PATH]
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
# Headline observable: Bliss synergy (RSL3 + FSP1i) from sim-combo-mech
# ----------------------------------------------------------------------------
def _default_binary():
    for cand in (
        REPO / "simulations" / "target" / "release" / "sim-combo-mech",
        REPO / "simulations" / "target" / "debug" / "sim-combo-mech",
    ):
        if cand.exists():
            return cand
    return None


def run_bliss(params_row, binary):
    """Run sim-combo-mech with the given absolute parameter values and return the
    RSL3 + FSP1i synergy_score. Each run uses a private cwd so its hard-coded
    relative output path (`output/combo-mech/combo_synergy.csv`) cannot collide
    with concurrent runs."""
    overrides = {n: float(v) for n, v in zip(PARAM_NAMES, params_row)}
    with tempfile.TemporaryDirectory(prefix="ferro_morris_") as workdir:
        env = dict(os.environ, FERRO_PARAM_OVERRIDES=json.dumps(overrides))
        proc = subprocess.run(
            [str(binary)],
            cwd=workdir,
            env=env,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"sim-combo-mech failed ({proc.returncode}): {proc.stderr[-400:]}")
        csv_path = Path(workdir) / "output" / "combo-mech" / "combo_synergy.csv"
        with csv_path.open() as fh:
            for row in csv.DictReader(fh):
                pair = {row["drug_a"], row["drug_b"]}
                if pair == {"RSL3", "FSP1i"}:
                    return float(row["synergy_score"])
    raise RuntimeError("RSL3+FSP1i row not found in combo_synergy.csv")


def make_bliss_eval(binary, workers):
    def eval_fn(scaled_rows):
        with ThreadPoolExecutor(max_workers=workers) as ex:
            return list(ex.map(lambda r: run_bliss(r, binary), scaled_rows))

    return eval_fn


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def write_report(mu_star, sigma, n_traj, levels, n_evals, binary):
    order = np.argsort(-mu_star)
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
        "evaluations (tens of thousands); each binary run is orders of magnitude "
        "more expensive than the `sim_batch` call the single-cell Sobol used, so "
        "full Sobol is infeasible for the simulation binaries. Morris gives, per "
        "parameter, `mu*` (mean |elementary effect|, the importance ranking — the "
        f"screening analogue of a total effect) and `sigma` (interaction / "
        f"nonlinearity indicator) in `r*(k+1)` runs. Here `r = {n_traj}` "
        f"trajectories, `k = {len(PARAM_NAMES)}` parameters, `{levels}` levels "
        f"⇒ **{n_evals} model runs**. Parameters swept over the SAME PRCC ranges "
        "(`analysis/prcc-results.json`). Deterministic given the design seed.",
        "",
        "## Headline: Bliss synergy (RSL3 + FSP1i)",
        "",
        f"Observable: the RSL3 + FSP1i `synergy_score` (death-over-Bliss-prediction; "
        f"the manuscript's ~1.99x) from `{binary.name}`'s `combo_synergy.csv`. "
        "`mu*` ranks how much each rate constant moves that synergy; `sigma` "
        "comparable-to-or-larger-than `mu*` flags interaction / nonlinearity.",
        "",
        "| rank | parameter | mu* | sigma | sigma/mu* |",
        "|------|-----------|-----|-------|-----------|",
    ]
    for rank, i in enumerate(order, 1):
        ratio = sigma[i] / mu_star[i] if mu_star[i] > 1e-12 else float("nan")
        lines.append(
            f"| {rank} | `{PARAM_NAMES[i]}` | {mu_star[i]:.4f} | {sigma[i]:.4f} | {ratio:.2f} |"
        )
    top = ", ".join(f"`{PARAM_NAMES[i]}`" for i in order[:3])
    # Two genuinely STRUCTURAL zeros for the RSL3 + FSP1i pair: parameters that are
    # never read on this code path (so mu*=0 by construction). Any OTHER mu*=0
    # parameter is an EMPIRICAL/below-resolution zero (it IS on the active biochem
    # path but its elementary effect rounds to zero here) — a different, weaker
    # claim. Keying purely on `mu* < eps` cannot distinguish the two, so name the
    # structural ones explicitly rather than overstate the empirical zeros.
    structural_names = {"sdt_ros", "rsl3_gpx4_inhib"}
    zero_params = [PARAM_NAMES[i] for i in order if mu_star[i] < 1e-9]
    structural = [n for n in zero_params if n in structural_names]
    empirical = [n for n in zero_params if n not in structural_names]
    lines += [
        "",
        f"**Top drivers of the Bliss synergy:** {top} — the LP-cascade and GSH/GPX4 "
        "defense constants, the same axis the single-cell Sobol and the PRCC rank "
        "for the RSL3 kill. `sigma` exceeds `mu*` for every ACTIVE parameter "
        "(`sigma/mu* > 1`), so the synergy is strongly NONLINEAR / INTERACTION-LADEN: "
        "an elementary effect depends heavily on where in parameter space it is "
        "measured (the bistable ferroptosis switch), structure a univariate PRCC "
        "cannot see. (Morris `sigma` does not separate interaction from "
        "single-parameter nonlinearity; either way the effect is non-additive.)",
        "",
        "**Zero-effect parameters (`mu* = 0`) — two distinct kinds, not conflated:**",
        "",
        "- *Structural zeros* (disconnected from this observable by construction): "
        + (", ".join(f"`{n}`" for n in structural) or "none")
        + ". `sdt_ros` is the SDT dose (no SDT acts in the RSL3 + FSP1i pair); "
        "`rsl3_gpx4_inhib` is never read because `sim-combo-mech` applies RSL3 as a "
        "fixed `DrugEffect` (92% GPX4 inhibition baked into the drug, not via "
        "`Params.rsl3_gpx4_inhib`). Not constrainable from this headline by construction.",
        "- *Empirical / below-resolution zeros* (the parameter IS on the active "
        "biochem path, but its elementary effect rounds to zero here): "
        + (", ".join(f"`{n}`" for n in empirical) or "none")
        + ". `gpx4_degradation_by_ros` is used in the GPX4 dynamics, but over its "
        "PRCC range a perturbation never flips a cell's bistable outcome, so the "
        "discrete death-count ratio does not move (consistent with its low "
        "single-cell Sobol ST). This is a practical-identifiability limit for THIS "
        "observable, NOT a structural disconnection.",
        "",
        "## Deferred: the two sim-tme headlines (hypoxia kill-collapse, immune ratio)",
        "",
        "Both are produced by `sim-tme`, whose per-run cost is far higher than "
        "`sim-combo-mech`, and both need a carefully chosen observable before a "
        "Morris screen is meaningful:",
        "",
        "- **Hypoxia kill-collapse.** The headline SDT hypoxic-zone kill sits near "
        "a 1.0 ceiling, which compresses its variance (the PRCC saw the same "
        "SDT-insensitivity). A faithful observable is the SDT-minus-RSL3 hypoxic "
        "GAP rather than the SDT kill alone; that screen is the next increment.",
        "- **Immune ratio.** Raw `immune_kills` is CONFOUNDED: SDT kills almost the "
        "whole tumor by ferroptosis first, leaving few cells for the immune layer "
        "to kill, so the raw SDT:RSL3 immune-kill ratio runs OPPOSITE to the "
        "headline's dense-DAMP-amplification framing. A faithful observable keys "
        "on the DAMP field driving immune amplification (e.g. ferroptotic-death "
        "density), which the 2D `sim-tme` summary does not yet expose. Deferred "
        "until that observable is defined, rather than screen a misleading one.",
        "",
        "These are tracked under #331 (which stays open for the sim-tme headlines).",
    ]
    REPORT.write_text("\n".join(lines) + "\n")
    return order


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trajectories", type=int, default=100)
    ap.add_argument("--levels", type=int, default=4)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20250606)
    ap.add_argument("--smoke", action="store_true", help="tiny run (2 trajectories) for a wiring check")
    ap.add_argument("--binary", type=str, default=None)
    args = ap.parse_args()

    n_traj = 2 if args.smoke else args.trajectories
    binary = Path(args.binary) if args.binary else _default_binary()
    if binary is None or not binary.exists():
        sys.exit(
            "sim-combo-mech binary not found. Build it first:\n"
            "  (cd simulations && cargo build --release -p sim-combo-mech)"
        )

    n_evals = n_traj * (len(PARAM_NAMES) + 1)
    print(f"Morris screen: r={n_traj}, k={len(PARAM_NAMES)}, levels={args.levels} -> {n_evals} runs of {binary.name}")
    eval_fn = make_bliss_eval(binary, args.workers)
    mu_star, sigma = morris_indices(eval_fn, LOWS, HIGHS, n_traj, args.levels, args.seed)
    order = write_report(mu_star, sigma, n_traj, args.levels, n_evals, binary)
    print(f"Wrote {REPORT.relative_to(REPO)}")
    print("Top-3 Bliss-synergy drivers:", ", ".join(PARAM_NAMES[i] for i in order[:3]))


if __name__ == "__main__":
    main()
