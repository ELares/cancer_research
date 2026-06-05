#!/usr/bin/env python3
"""Prior-predictive uncertainty intervals on the single-cell kill switch (#332).

The manuscript's headline simulation numbers are reported as POINT estimates from
a single parameter set, when nearly every biochemical rate constant is an
uncertain placeholder (see `simulations/calibration/CALIBRATION_STATUS.md`,
`parameter_provenance.md`). A contribution-grade model should report INTERVALS.

This is the foundational, in-repo-tractable piece of that: a prior-predictive
uncertainty propagation over the single-cell ferroptosis death rate (the
`ferroptosis_core` `sim_batch` switch the spatial headline claims sit on top of).
We place a UNIFORM PRIOR over each biochemical rate constant across its PRCC range
(`analysis/prcc-results.json`, the same plausible ranges the PRCC #134 / Sobol
#331 used), Monte-Carlo sample the joint prior, run the forward model for each
phenotype x treatment condition, and report the resulting credible interval
(2.5 / 50 / 97.5 percentiles) on the death rate alongside the default-parameter
point estimate.

HONEST SCOPE — this is a PRIOR-PREDICTIVE distribution, NOT a data-conditioned
Bayesian POSTERIOR. We have no in-repo kill-rate calibration dataset to condition
on, so this reports how much the death rate MOVES under the documented parameter
uncertainty (forward uncertainty), not a calibrated posterior. Turning the
uniform priors into a posterior needs the external GDSC/DepMap/PRISM drug-response
data (issue #330) fed through an ABC/likelihood step; that is the remaining #332
work. Likewise this covers the SINGLE-CELL outputs; propagating the same priors
through the SPATIAL headline outputs (Bliss synergy, hypoxia collapse, immune
ratio, penetration) needs a spatial-sim ensemble harness and is deferred (the
same single-cell-first scoping as #331).

Reproducible: fixed numpy + simulation seeds; reads only the committed binding and
`prcc-results.json`. Writes `analysis/uncertainty-intervals-report.md`.
"""

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "simulations" / "ferroptosis-python"))

# The compiled `ferroptosis_core` extension is imported LAZILY so the pure-numpy
# estimator (`prior_predictive_intervals`) and its analytic validation can run
# without it (the Python CI does not build the extension). Only the model-specific
# functions need it; they call `_fc()`.
_FC = None


def _fc():
    global _FC
    if _FC is None:
        import ferroptosis_core

        _FC = ferroptosis_core
    return _FC


REPORT = REPO / "analysis" / "uncertainty-intervals-report.md"

N_CELLS = 2000  # cells per sim_batch (Monte Carlo) evaluation
SIM_SEED = 42  # fixed across samples so the death rate is a deterministic function
# of the parameters (the per-cell RNG is reproducible), isolating the PARAMETER
# uncertainty from sim_batch Monte Carlo noise.
N_SAMPLES = 2000  # prior draws (joint uniform over the PRCC ranges)
RNG_SEED = 20250605

# The phenotype x treatment cells whose death rates underlie the manuscript's
# single-cell headline comparisons (Figure 7 / §5: Persister vs Glycolytic, RSL3
# vs SDT). OXPHOS is the RSL3-resistant counterpoint.
CONDITIONS = [
    ("Persister", "RSL3"),
    ("Persister", "SDT"),
    ("Glycolytic", "RSL3"),
    ("Glycolytic", "SDT"),
    ("OXPHOS", "RSL3"),
    ("OXPHOS", "SDT"),
]

# Uniform priors over the SAME biochemical rate constants + ranges as the PRCC
# (#134) and the Sobol screen (#331), loaded from the committed PRCC results so
# the three analyses share one parameter-range source of truth. Unlike the Sobol
# (which excluded the SDT-only `sdt_ros` for its RSL3 observable), this propagates
# ALL of them, since the condition set includes SDT, where `sdt_ros` is the
# dominant driver.
_PRCC_RANGES = json.loads((REPO / "analysis" / "prcc-results.json").read_text())[
    "metadata"
]["parameter_ranges"]
PARAMS = [(name, lo, hi) for name, (lo, hi) in _PRCC_RANGES.items()]


def prior_predictive_intervals(eval_fn, lows, highs, n_samples, rng_seed, quantiles):
    """Generic prior-predictive propagation: draw `n_samples` joint-uniform
    parameter vectors over the box `[lows, highs]`, push them through the
    vector-valued model `eval_fn` (maps an `(n_samples, k)` array to an
    `(n_samples, n_outputs)` array), and return, per output, the requested
    `quantiles` plus the mean and std. Estimator-only, so it can be validated on
    an analytic distribution. Returns `(stats, samples_out)` where `stats` is an
    `(n_outputs, len(quantiles) + 2)` array (quantiles..., mean, std)."""
    lows = np.asarray(lows, float)
    highs = np.asarray(highs, float)
    k = lows.shape[0]
    rng = np.random.default_rng(rng_seed)
    draws = lows + rng.random((n_samples, k)) * (highs - lows)
    out = np.atleast_2d(eval_fn(draws))
    if out.shape[0] == n_samples:  # eval returned (n_samples, n_outputs)
        out = out.T  # -> (n_outputs, n_samples)
    qs = np.quantile(out, quantiles, axis=1).T  # (n_outputs, len(quantiles))
    mean = out.mean(axis=1, keepdims=True)
    std = out.std(axis=1, keepdims=True)
    stats = np.hstack([qs, mean, std])
    return stats, out


def evaluate(sample_rows):
    """Death rate for every (sample, condition): returns an
    `(n_samples, n_conditions)` array. Each row is one joint-uniform parameter
    draw (ABSOLUTE values over the PRCC ranges), applied as `sim_batch` overrides;
    each column is one phenotype x treatment condition."""
    names = [p[0] for p in PARAMS]
    out = np.empty((len(sample_rows), len(CONDITIONS)))
    for i, row in enumerate(sample_rows):
        overrides = {n: float(row[j]) for j, n in enumerate(names)}
        for c, (pheno, treat) in enumerate(CONDITIONS):
            res = _fc().sim_batch(pheno, treat, n=N_CELLS, seed=SIM_SEED, **overrides)
            out[i, c] = res["death_rate"]
    return out


def point_estimates():
    """The default-parameter death rate for each condition (the manuscript's
    point-estimate operating point), for side-by-side comparison with the
    interval."""
    return [
        _fc().sim_batch(pheno, treat, n=N_CELLS, seed=SIM_SEED)["death_rate"]
        for (pheno, treat) in CONDITIONS
    ]


QUANTILES = [0.025, 0.5, 0.975]


def run():
    lows = np.array([p[1] for p in PARAMS], float)
    highs = np.array([p[2] for p in PARAMS], float)
    stats, _ = prior_predictive_intervals(
        evaluate, lows, highs, N_SAMPLES, RNG_SEED, QUANTILES
    )
    points = point_estimates()
    return stats, points


def fmt_pct(x):
    return f"{100 * x:.1f}%"


def write_report(stats, points):
    lines = []
    lines.append("# Prior-predictive uncertainty intervals on the kill switch (#332)\n")
    lines.append(
        "Generated by `scripts/uncertainty_intervals.py`. Companion to the univariate "
        "PRCC (#134) and the variance-based Sobol screen (#331): where those rank "
        "WHICH parameters matter, this reports HOW MUCH the single-cell death rate "
        "moves under the documented parameter uncertainty.\n"
    )
    lines.append(
        f"**Method.** Uniform priors over the {len(PARAMS)} biochemical rate constants, "
        f"each across its PRCC range (`analysis/prcc-results.json`); "
        f"{N_SAMPLES} joint draws propagated through `sim_batch` "
        f"(n={N_CELLS}, seed={SIM_SEED}) for each phenotype x treatment condition; "
        f"per-condition 2.5 / 50 / 97.5 percentiles of the resulting death-rate "
        f"distribution. Fixed seeds (numpy {RNG_SEED}, sim {SIM_SEED}).\n"
    )
    lines.append(
        "**This is a PRIOR-PREDICTIVE distribution, not a data-conditioned posterior.** "
        "With no in-repo kill-rate calibration dataset, this is the FORWARD uncertainty "
        "(how much the death rate moves given the parameter uncertainty), not a "
        "calibrated Bayesian posterior. A true posterior needs the external "
        "GDSC/DepMap/PRISM drug-response data (#330) through an ABC/likelihood step. "
        "It also covers the SINGLE-CELL switch only; the SPATIAL headline outputs "
        "(Bliss synergy, hypoxia collapse, immune ratio, penetration) need a "
        "spatial-sim ensemble and are deferred (same single-cell-first scoping as "
        "#331).\n"
    )
    lines.append("## Death-rate intervals\n")
    lines.append(
        "`point` = default-parameter death rate (the manuscript's operating point); "
        "`median`, `2.5%`, `97.5%` = percentiles of the prior-predictive "
        "distribution; `width` = 97.5% - 2.5% (the size of the credible interval).\n"
    )
    lines.append("| Phenotype | Treatment | point | median | 2.5% | 97.5% | width |")
    lines.append("|---|---|--:|--:|--:|--:|--:|")
    for (pheno, treat), st, pt in zip(CONDITIONS, stats, points):
        lo, med, hi = st[0], st[1], st[2]
        lines.append(
            f"| {pheno} | {treat} | {fmt_pct(pt)} | {fmt_pct(med)} | "
            f"{fmt_pct(lo)} | {fmt_pct(hi)} | {fmt_pct(hi - lo)} |"
        )
    lines.append("")
    lines.append("## Reading the intervals\n")
    lines.append(
        "- The point estimates are NOT the centre of the prior-predictive "
        "distribution: the manuscript runs the default parameter set, whereas these "
        "intervals span the full plausible range, so a wide interval means the "
        "headline number is sensitive to the (uncalibrated) parameter choice and "
        "should be read as order-of-magnitude.\n"
        "- **The RSL3 death rates carry the widest intervals.** The bistable "
        "switch's tipping point sits inside the swept range, so the Persister x RSL3 "
        "death rate spans nearly the full [0, 1] under the parameter uncertainty: "
        "the manuscript's RSL3 point estimates are the LEAST constrained and most "
        "order-of-magnitude. The Glycolytic x RSL3 floor (~0% across the entire "
        "prior) is the robust exception, and OXPHOS x RSL3 stays low-but-wide.\n"
        "- **SDT death is near-ceiling only for the high-basal-ROS phenotypes.** "
        "Persister and OXPHOS x SDT sit near 1 across the prior (a genuine ceiling), "
        "but the Glycolytic x SDT interval is WIDE because the lower end of the "
        "`sdt_ros` prior does not overwhelm glycolytic defenses, so even an SDT "
        "advantage is not uniformly parameter-robust.\n"
        "- Because several intervals are wide and OVERLAP, the qualitative ordering "
        "the manuscript leans on (Persister vs Glycolytic, SDT vs RSL3) should be "
        "read against interval OVERLAP, not point estimates alone.\n"
    )
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT}")


def main():
    stats, points = run()
    write_report(stats, points)
    for (pheno, treat), st, pt in zip(CONDITIONS, stats, points):
        print(
            f"  {pheno:10s} {treat:5s} point={fmt_pct(pt)} "
            f"median={fmt_pct(st[1])} [{fmt_pct(st[0])}, {fmt_pct(st[2])}]"
        )


if __name__ == "__main__":
    main()
