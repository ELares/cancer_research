#!/usr/bin/env python3
"""Joint multi-inducer ABC posterior for the single-cell switch (#500).

#330 anchored the switch to a SINGLE inducer (RSL3/GPX4i). #332 produced an ABC
posterior but on that one inducer, and showed the in-vivo PRCC priors and the
in-vitro data are disjoint. #502 added System Xc-/erastin to the core. This script
closes #500: it conditions the SHARED single-cell switch jointly on a multi-inducer
panel (a GPX4 inhibitor AND a System Xc- inhibitor at once) and reports a POSTERIOR
with credible intervals, not another tuned point.

WHAT IS JOINTLY FIT
-------------------
A single shared cascade drives BOTH inducers; only the per-mechanism potency
differs:
  * shared:   lp_propagation, lp_rate, gpx4_rate, gsh_scav_efficiency
  * RSL3:     k_um           (GPX4i dose -> rsl3_gpx4_inhib = dose/(dose+k_um))
  * erastin:  k_erastin, hill (System Xc- dose -> erastin_xc_inhib =
                               dose^h/(dose^h+k_erastin^h), run under Control)

ABC: draw from in-vitro-spanning uniform priors, simulate BOTH dose-response
curves per draw, score by the JOINT distance (RSL3 RMSE vs ML162 + erastin RMSE
vs CTRPv2 erastin), accept the closest fraction -> the joint posterior. Marginal
2.5/50/97.5 credible intervals per parameter.

POSTERIOR-PREDICTIVE CHECK
--------------------------
The accepted draws predict a HELD-OUT GPX4 inhibitor (ML210, never used in the
distance), reported as a posterior-predictive RMSE distribution + coverage of the
held-out points inside the 95% predictive band. This is generalization, not
training fit.

HONESTY / SCOPE (the load-bearing caveat)
-----------------------------------------
This is an IN-VITRO joint posterior. It puts credible intervals on the IN-VITRO
single-cell switch. It does NOT condition the in-vivo / spatial manuscript
headlines: the in-vivo PRCC priors and the in-vitro data are DISJOINT (#332,
re-confirmed here for the joint fit), so the in-vivo/spatial headline numbers stay
PRIOR-predictive (the existing headline_uncertainty intervals). The manuscript
switch magnitudes therefore carry the in-vitro joint posterior intervals where the
claim is in-vitro, and stay prior-predictive where the claim is in-vivo/spatial,
and this script states which is which. Conditioning the in-vivo headlines needs an
in-vivo ferroptosis dataset that does not exist publicly.

Run (needs the compiled `ferroptosis_core` extension; not run in CI):
  python3 scripts/abc_joint_posterior.py
Writes analysis/calibration/joint-posterior.{md,json} + joint-posterior-predictive.png.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import calibrate_erastin as ce  # noqa: E402
import calibrate_kill_switch as ck  # noqa: E402

OUT_MD = REPO_ROOT / "analysis" / "calibration" / "joint-posterior.md"
OUT_JSON = REPO_ROOT / "analysis" / "calibration" / "joint-posterior.json"
OUT_PNG = REPO_ROOT / "analysis" / "calibration" / "joint-posterior-predictive.png"

# In-vitro-spanning uniform priors (name, lo, hi). The cascade ranges match the
# #332 in-vitro-spanning priors; the erastin potency ranges match the #502 grid.
PRIORS = (
    ("lp_propagation", 0.05, 1.0),
    ("lp_rate", 0.03, 1.0),
    ("gpx4_rate", 0.15, 0.6),
    ("gsh_scav_efficiency", 0.2, 0.9),
    ("k_um", 0.1, 4.0),       # RSL3 (GPX4i) potency
    ("k_erastin", 3.0, 50.0),  # erastin (System Xc-) potency
    ("hill", 1.0, 6.0),        # erastin Hill slope
)
SHARED = ("lp_propagation", "lp_rate", "gpx4_rate", "gsh_scav_efficiency")

# The in-vivo PRCC ranges the in-vitro posterior supersedes (the disjunction).
INVIVO_PRCC = {"lp_propagation": (0.05, 0.2), "lp_rate": (0.03, 0.12)}

N_DRAWS = 1500
ACCEPT_FRAC = 0.02
RNG_SEED = 12345
SIM_N = 2000
SIM_SEED = 42
QUANTILES = (2.5, 50.0, 97.5)


def _shared(params):
    return {k: params[k] for k in SHARED}


def model_rsl3(doses, params, n=SIM_N, seed=SIM_SEED):
    fc = ck._fc()
    out = []
    for d in doses:
        inhib = ck.dose_to_inhib(d, params["k_um"])
        death = fc.sim_batch(ck.PHENOTYPE, "RSL3", n=n, seed=seed,
                             rsl3_gpx4_inhib=inhib, **_shared(params))["death_rate"]
        out.append(1.0 - death)
    return out


def model_erastin(doses, params, n=SIM_N, seed=SIM_SEED):
    fc = ck._fc()
    out = []
    for d in doses:
        inhib = ce.dose_to_inhib(d, params["k_erastin"], params["hill"])
        death = fc.sim_batch(ce.PHENOTYPE, "Control", n=n, seed=seed,
                             erastin_xc_inhib=inhib, **_shared(params))["death_rate"]
        out.append(1.0 - death)
    return out


def run(args):
    curves = ck.load_curves()
    rsl3_doses = list(ck.DOSE_GRID_UM)
    erastin_doses = list(ce.DOSE_GRID_UM)

    emp_rsl3 = ck.empirical_median_viability(curves[ck.FIT_COMPOUND], rsl3_doses)        # ML162
    emp_erastin = ck.empirical_median_viability(curves[ce.COMPOUND], erastin_doses)      # erastin
    emp_heldout = ck.empirical_median_viability(curves[ck.HELDOUT_GPX4I], rsl3_doses)    # ML210

    rng = np.random.default_rng(RNG_SEED)
    names = [p[0] for p in PRIORS]
    lows = np.array([p[1] for p in PRIORS])
    highs = np.array([p[2] for p in PRIORS])
    draws = rng.uniform(lows, highs, size=(args.n_draws, len(PRIORS)))

    distances = np.empty(args.n_draws)
    rsl3_models = []
    erastin_models = []
    for i in range(args.n_draws):
        params = dict(zip(names, draws[i]))
        m_rsl3 = model_rsl3(rsl3_doses, params)
        m_erastin = model_erastin(erastin_doses, params)
        rsl3_models.append(m_rsl3)
        erastin_models.append(m_erastin)
        # Joint distance: equally-weighted sum of the two per-curve RMSEs.
        distances[i] = ck.rmse(m_rsl3, emp_rsl3) + ck.rmse(m_erastin, emp_erastin)

    n_accept = max(10, int(args.n_draws * ACCEPT_FRAC))
    accept_idx = np.argsort(distances)[:n_accept]
    posterior = draws[accept_idx]
    eps = float(distances[accept_idx].max())

    post = {}
    for j, name in enumerate(names):
        q = np.percentile(posterior[:, j], QUANTILES)
        width = float(highs[j] - lows[j])
        post_width = float(q[2] - q[0])
        post[name] = {
            "q2_5": round(float(q[0]), 4),
            "median": round(float(q[1]), 4),
            "q97_5": round(float(q[2]), 4),
            # Fraction of the prior width the posterior still occupies: ~1 ==
            # unconstrained by the data, small == well constrained.
            "posterior_width_frac_of_prior": round(post_width / width, 3),
        }

    # Disjunction with the in-vivo PRCC ranges (the #332 finding, re-checked jointly).
    disjoint = {}
    for name, (lo, hi) in INVIVO_PRCC.items():
        disjoint[name] = {
            "invivo_prcc_range": [lo, hi],
            "posterior_q2_5": post[name]["q2_5"],
            "posterior_median": post[name]["median"],
            "entire_95pct_posterior_above_invivo_max": post[name]["q2_5"] > hi,
        }

    # Posterior-predictive on held-out ML210: each accepted draw predicts ML210.
    acc_rsl3 = np.array([rsl3_models[i] for i in accept_idx])          # training curve
    heldout_pred = np.array([model_rsl3(rsl3_doses, dict(zip(names, posterior[k])))
                             for k in range(len(posterior))])
    pp_band = np.percentile(heldout_pred, [2.5, 50.0, 97.5], axis=0)
    emp_h = np.array(emp_heldout)
    inside = (emp_h >= pp_band[0]) & (emp_h <= pp_band[2])
    pp_rmses = [ck.rmse(heldout_pred[k], emp_heldout) for k in range(len(posterior))]

    # Which shared params are unconstrained (posterior ~ prior).
    unconstrained = [n for n in names if post[n]["posterior_width_frac_of_prior"] >= 0.6]

    result = {
        "n_draws": args.n_draws,
        "n_accepted": int(n_accept),
        "accept_frac": ACCEPT_FRAC,
        "epsilon_joint_distance": round(eps, 4),
        "inducer_panel": {
            "fit": [ck.FIT_COMPOUND, ce.COMPOUND],
            "heldout": ck.HELDOUT_GPX4I,
        },
        "posterior": post,
        "unconstrained_params": unconstrained,
        "disjunction_with_invivo_priors": disjoint,
        "heldout_posterior_predictive": {
            "compound": ck.HELDOUT_GPX4I,
            "coverage_inside_95pct_band": f"{int(inside.sum())}/{len(emp_heldout)}",
            "median_pp_rmse": round(float(np.median(pp_rmses)), 4),
            "band_q2_5": [round(float(v), 4) for v in pp_band[0]],
            "band_median": [round(float(v), 4) for v in pp_band[1]],
            "band_q97_5": [round(float(v), 4) for v in pp_band[2]],
            "empirical": [round(v, 4) for v in emp_heldout],
        },
        "curves": {
            "rsl3_doses_um": rsl3_doses,
            "erastin_doses_um": erastin_doses,
            "empirical_rsl3_ml162": [round(v, 4) for v in emp_rsl3],
            "empirical_erastin": [round(v, 4) for v in emp_erastin],
            "posterior_median_rsl3": [round(float(v), 4) for v in np.percentile(acc_rsl3, 50, axis=0)],
            "posterior_median_erastin": [
                round(float(v), 4)
                for v in np.percentile(np.array([erastin_models[i] for i in accept_idx]), 50, axis=0)
            ],
        },
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    _plot(result, accept_idx, rsl3_models, erastin_models, pp_band)
    write_report(result)
    print(f"accepted {n_accept}/{args.n_draws}; eps={eps:.4f}")
    print(f"held-out ML210 coverage {result['heldout_posterior_predictive']['coverage_inside_95pct_band']}, "
          f"median PP RMSE {result['heldout_posterior_predictive']['median_pp_rmse']}")
    print(f"unconstrained: {unconstrained}")
    print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)} + {OUT_MD.relative_to(REPO_ROOT)} + {OUT_PNG.name}")
    return result


def _plot(r, accept_idx, rsl3_models, erastin_models, pp_band):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = r["curves"]
    rsl3_doses = c["rsl3_doses_um"]
    erastin_doses = c["erastin_doses_um"]
    acc_rsl3 = np.array([rsl3_models[i] for i in accept_idx])
    acc_er = np.array([erastin_models[i] for i in accept_idx])
    r_band = np.percentile(acc_rsl3, [2.5, 50, 97.5], axis=0)
    e_band = np.percentile(acc_er, [2.5, 50, 97.5], axis=0)
    hp = r["heldout_posterior_predictive"]

    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    for a in ax:
        a.set_xscale("log")
        a.set_ylim(0, 1.1)
        a.set_xlabel("dose (µM)")
        a.set_ylabel("viability")

    ax[0].fill_between(rsl3_doses, r_band[0], r_band[2], alpha=0.25, color="C0", label="95% posterior")
    ax[0].plot(rsl3_doses, r_band[1], "-", color="C0", label="posterior median")
    ax[0].plot(rsl3_doses, c["empirical_rsl3_ml162"], "ko", label="ML162 (fit)")
    ax[0].set_title("RSL3 / GPX4i (ML162): fit")
    ax[0].legend(fontsize=7)

    ax[1].fill_between(erastin_doses, e_band[0], e_band[2], alpha=0.25, color="C1", label="95% posterior")
    ax[1].plot(erastin_doses, e_band[1], "-", color="C1", label="posterior median")
    ax[1].plot(erastin_doses, c["empirical_erastin"], "ks", label="erastin (fit)")
    ax[1].set_title("erastin / System Xc-: fit")
    ax[1].legend(fontsize=7)

    ax[2].fill_between(rsl3_doses, pp_band[0], pp_band[2], alpha=0.25, color="C2", label="95% predictive")
    ax[2].plot(rsl3_doses, pp_band[1], "-", color="C2", label="predictive median")
    ax[2].plot(rsl3_doses, hp["empirical"], "k^", label=f"{hp['compound']} (held-out)")
    ax[2].set_title(f"{hp['compound']}: held-out posterior-predictive")
    ax[2].legend(fontsize=7)

    fig.suptitle("Joint multi-inducer in-vitro posterior (#500): RSL3 + erastin", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=120)
    plt.close(fig)


def write_report(r):
    p = r["posterior"]
    d = r["disjunction_with_invivo_priors"]
    hp = r["heldout_posterior_predictive"]

    def ptab():
        lines = ["| parameter | 2.5% | median | 97.5% | width (frac of prior) |", "|---|---|---|---|---|"]
        for name, v in p.items():
            lines.append(f"| `{name}` | {v['q2_5']} | {v['median']} | {v['q97_5']} | {v['posterior_width_frac_of_prior']} |")
        return "\n".join(lines)

    disj_lines = "\n".join(
        f"- `{n}`: in-vivo PRCC range {v['invivo_prcc_range']}, in-vitro posterior median {v['posterior_median']} "
        f"(2.5% = {v['posterior_q2_5']}); entire 95% posterior above the in-vivo max: "
        f"**{v['entire_95pct_posterior_above_invivo_max']}**."
        for n, v in d.items()
    )

    md = f"""# Joint multi-inducer in-vitro posterior (#500)

Generated by `scripts/abc_joint_posterior.py` (needs the compiled `ferroptosis_core`
extension; not run in CI). Builds on the #330 GPX4i fit, the #502 System Xc-/erastin
mechanism, and the #332 single-inducer ABC.

## What this is

The shared single-cell switch (`lp_propagation`, `lp_rate`, `gpx4_rate`,
`gsh_scav_efficiency`) is conditioned JOINTLY on **two inducer mechanisms at once**
(a GPX4 inhibitor, ML162, and a System Xc- inhibitor, erastin), with per-mechanism
potencies (`k_um` for RSL3; `k_erastin`, `hill` for erastin). ABC over {r['n_draws']}
in-vitro-spanning prior draws, accepting the closest {r['accept_frac']:.0%}
({r['n_accepted']} draws) by the joint distance (ML162 RMSE + erastin RMSE),
epsilon = {r['epsilon_joint_distance']}.

## Joint posterior (credible intervals, not a point)

{ptab()}

The `width (frac of prior)` column is how much of the prior range the posterior
still spans: ~1.0 means the data barely constrains it, small means well
constrained. **Unconstrained (>= 0.6 of prior width):** {", ".join(f"`{x}`" for x in r['unconstrained_params']) or "none"}.
These are the parameters the in-vitro dose-response panel does not identify (e.g.
the GSH/GPX4 axis is partly degenerate with the LP cascade, consistent with the
PRCC/Sobol identifiability findings); they are reported as intervals, not points.

## Held-out posterior-predictive ({hp['compound']}, never used in the fit)

The accepted draws predict a held-out GPX4 inhibitor ({hp['compound']}):
**coverage {hp['coverage_inside_95pct_band']} of held-out points inside the 95%
predictive band**, median posterior-predictive RMSE **{hp['median_pp_rmse']}**. This
is generalization to an unseen inducer of the same class, not training fit. See
`joint-posterior-predictive.png` (right panel).

## The load-bearing caveat: in-vitro only, disjoint from the in-vivo priors

The joint posterior is **in-vitro**. Re-checking the #332 disjunction for the joint
fit:

{disj_lines}

So the in-vitro joint posterior lies ABOVE the in-vivo PRCC priors used for the
spatial/headline prior-predictive intervals. The consequence is unchanged from
#332 and is the honest scope of #500:

- **In-vitro switch claims** carry these joint-posterior credible intervals.
- **In-vivo / spatial headline numbers** (hypoxia asymmetry, Bliss synergy,
  penetration gap, immune ratio) **cannot** be conditioned on in-vitro data and
  stay **prior-predictive** (the existing `headline_uncertainty.py` /
  `uncertainty_intervals.py` intervals, `analysis/identifiability-report.md`).
  Conditioning them needs an in-vivo ferroptosis dataset that does not exist
  publicly.

This is the posterior #500 asked for (a real multi-inducer posterior with
held-out generalization), reported with exactly the scope the data supports.
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-draws", type=int, default=N_DRAWS)
    args = ap.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
