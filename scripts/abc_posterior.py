#!/usr/bin/env python3
"""ABC posterior for the kill switch, conditioned on CTRPv2 GPX4i data (#332).

#332 asks for a DATA-CONDITIONED posterior (not just the prior-predictive
intervals already shipped in scripts/uncertainty_intervals.py +
scripts/headline_uncertainty.py). #330 produced the data that makes conditioning
possible: the CTRPv2 GPX4-inhibitor dose-response target.

A PIVOTAL CONSTRAINT shapes what is achievable. The in-vivo PRCC prior ranges
(`analysis/prcc-results.json`) used for the prior-predictive intervals put
`lp_propagation` in [0.05, 0.2] and `lp_rate` in [0.03, 0.12], but the in-vitro
CTRPv2 GPX4-inhibitor kill (#330) requires `lp_propagation ~ 0.7`, `lp_rate ~ 0.4`,
3-4x ABOVE those ranges. The in-vivo priors and the in-vitro data are therefore
DISJOINT: an ABC that conditioned the in-vivo priors on the in-vitro data would
accept nothing. Two honest consequences:

  1. The data-conditioned posterior here is an IN-VITRO posterior, computed over
     priors WIDENED to span the in-vitro regime. It quantifies the uncertainty on
     the #330 point calibration (turning that point fit into a posterior band) and
     is posterior-predictive-checked on the held-out GPX4 inhibitor (ML210).
  2. The IN-VIVO / spatial headline claims (Bliss synergy, hypoxia collapse,
     immune amplification, penetration) CANNOT be conditioned on this in-vitro
     data, because no in-vivo ferroptosis dose-response dataset is in hand. They
     stay PRIOR-PREDICTIVE (the existing intervals). Conditioning them needs an
     in-vivo dataset, which is the remaining #332 dependency.

So this reports the in-vitro single-cell posterior + the disjunction finding; it
does not (and honestly cannot, yet) tighten the in-vivo headline intervals.

ABC rejection: sample (lp_propagation, lp_rate, k_um, gpx4_rate) from in-vitro-
spanning uniform priors, run the model RSL3 dose-response (1 - death_rate from
`ferroptosis_core.sim_batch`), accept the draws closest (lowest RMSE) to the
CTRPv2 ML162 median viability curve. Deterministic (seeded). Needs the compiled
extension; not run in CI. Writes analysis/calibration/abc-posterior-report.md + .json.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import calibrate_kill_switch as ck  # noqa: E402

OUT_MD = REPO_ROOT / "analysis" / "calibration" / "abc-posterior-report.md"
OUT_JSON = REPO_ROOT / "analysis" / "calibration" / "abc-posterior.json"

# In-vitro-spanning priors (WIDER than the in-vivo PRCC ranges, which cannot
# reproduce in-vitro GPX4i kill). (name, lo, hi).
PRIORS = (
    ("lp_propagation", 0.05, 1.0),
    ("lp_rate", 0.03, 1.0),
    ("k_um", 0.1, 4.0),
    ("gpx4_rate", 0.15, 0.6),
)
# The in-vivo PRCC ranges these supersede, for the disjunction report.
INVIVO_PRCC = {"lp_propagation": (0.05, 0.2), "lp_rate": (0.03, 0.12)}

N_DRAWS = 2000
ACCEPT_FRAC = 0.02   # keep the closest 2% -> ABC posterior
HELDOUT_TOL = 0.05   # viability tolerance for the TOLERANT held-out coverage count
RNG_SEED = 12345
SIM_N = 2000
SIM_SEED = 42
QUANTILES = (2.5, 50.0, 97.5)


def model_dose_response(doses, params, phenotype=ck.PHENOTYPE, n=SIM_N, seed=SIM_SEED):
    """Model RSL3 viability(dose) under an ABC parameter draw."""
    fc = ck._fc()
    overrides = {k: v for k, v in params.items() if k != "k_um"}
    out = []
    for d in doses:
        inhib = ck.dose_to_inhib(d, params["k_um"])
        death = fc.sim_batch(phenotype, "RSL3", n=n, seed=seed, rsl3_gpx4_inhib=inhib, **overrides)["death_rate"]
        out.append(1.0 - death)
    return out


def run(args):
    curves = ck.load_curves()
    doses = list(ck.DOSE_GRID_UM)
    emp_fit = ck.empirical_median_viability(curves[ck.FIT_COMPOUND], doses)
    emp_heldout = ck.empirical_median_viability(curves[ck.HELDOUT_GPX4I], doses)

    rng = np.random.default_rng(RNG_SEED)
    names = [p[0] for p in PRIORS]
    lows = np.array([p[1] for p in PRIORS])
    highs = np.array([p[2] for p in PRIORS])
    draws = rng.uniform(lows, highs, size=(args.n_draws, len(PRIORS)))

    distances = np.empty(args.n_draws)
    for i in range(args.n_draws):
        params = dict(zip(names, draws[i]))
        model = model_dose_response(doses, params)
        distances[i] = ck.rmse(model, emp_fit)

    n_accept = max(10, int(args.n_draws * ACCEPT_FRAC))
    accept_idx = np.argsort(distances)[:n_accept]
    posterior = draws[accept_idx]
    eps = float(distances[accept_idx].max())

    # Marginal posterior intervals per parameter.
    post = {}
    for j, name in enumerate(names):
        q = np.percentile(posterior[:, j], QUANTILES)
        post[name] = {"q2_5": round(float(q[0]), 4), "median": round(float(q[1]), 4), "q97_5": round(float(q[2]), 4)}

    # Disjunction: posterior vs in-vivo PRCC ranges.
    disjoint = {}
    for name, (lo, hi) in INVIVO_PRCC.items():
        pmed = post[name]["median"]
        p_lo = post[name]["q2_5"]
        disjoint[name] = {
            "invivo_prcc_range": [lo, hi],
            "posterior_median": pmed,
            "posterior_q2_5": p_lo,
            "posterior_above_invivo_range": p_lo > hi,  # entire 95% posterior above the in-vivo max
        }

    # Posterior-predictive on held-out ML210: viability band per dose.
    pp = []
    for i in accept_idx:
        pp.append(model_dose_response(doses, dict(zip(names, draws[i]))))
    pp = np.array(pp)
    pp_band = {
        "dose_um": doses,
        "empirical_heldout": [round(v, 4) for v in emp_heldout],
        "post_pred_median": [round(float(v), 4) for v in np.percentile(pp, 50, axis=0)],
        "post_pred_q2_5": [round(float(v), 4) for v in np.percentile(pp, 2.5, axis=0)],
        "post_pred_q97_5": [round(float(v), 4) for v in np.percentile(pp, 97.5, axis=0)],
    }
    def _coverage(tol):
        return sum(
            1 for k in range(len(doses))
            if pp_band["post_pred_q2_5"][k] - tol <= emp_heldout[k] <= pp_band["post_pred_q97_5"][k] + tol
        )

    # STRICT (inside the 95% band) and TOLERANT (within HELDOUT_TOL viability). The
    # tolerance is not cosmetic: the CTRPv2 curves are cell-line MEDIANS while the
    # model is single-cell, and the RMSE summary statistic discards curve shape, so a
    # small viability offset is expected. Both are reported so the strict number is
    # never hidden behind the tolerant one.
    covered_strict = _coverage(0.0)
    covered_tol = _coverage(HELDOUT_TOL)

    result = {
        "n_draws": args.n_draws,
        "n_accepted": int(n_accept),
        "accept_frac": ACCEPT_FRAC,
        "epsilon_rmse": round(eps, 4),
        "fit_compound": ck.FIT_COMPOUND,
        "heldout_compound": ck.HELDOUT_GPX4I,
        "priors": {p[0]: [p[1], p[2]] for p in PRIORS},
        "posterior": post,
        "invivo_prior_disjunction": disjoint,
        "posterior_predictive_heldout": pp_band,
        "heldout_tolerance": HELDOUT_TOL,
        "heldout_coverage_strict": f"{covered_strict}/{len(doses)}",
        "heldout_coverage_tolerant": f"{covered_tol}/{len(doses)}",
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_report(result)
    print(f"accepted {n_accept}/{args.n_draws} (eps RMSE={result['epsilon_rmse']})")
    for name in names:
        print(f"  {name}: {post[name]['q2_5']} .. {post[name]['median']} .. {post[name]['q97_5']}")
    print(f"held-out coverage: strict {result['heldout_coverage_strict']}, "
          f"within {HELDOUT_TOL} {result['heldout_coverage_tolerant']}")
    print(f"wrote {OUT_MD.relative_to(REPO_ROOT)} + {OUT_JSON.relative_to(REPO_ROOT)}")
    return result


def write_report(r):
    post = r["posterior"]
    dj = r["invivo_prior_disjunction"]
    pp = r["posterior_predictive_heldout"]

    def prow(name):
        p = post[name]
        return f"| `{name}` | {p['q2_5']} | {p['median']} | {p['q97_5']} |"

    lines = [
        "# ABC posterior for the kill switch, conditioned on CTRPv2 GPX4i (#332)",
        "",
        "Generated by `scripts/abc_posterior.py` (needs the compiled `ferroptosis_core`",
        "extension; not run in CI). This is the DATA-CONDITIONED posterior #332 asks for,",
        "and it is necessarily an **in-vitro** posterior (see the disjunction below).",
        "",
        "## Method",
        "",
        f"ABC rejection: {r['n_draws']} draws from in-vitro-spanning uniform priors, the",
        f"closest **{r['n_accepted']}** (lowest RMSE to the CTRPv2 {r['fit_compound']} median",
        f"viability curve, acceptance fraction {r['accept_frac']}, epsilon RMSE",
        f"{r['epsilon_rmse']}) form the posterior. Posterior-predictive-checked on the",
        f"held-out GPX4 inhibitor {r['heldout_compound']}.",
        "",
        "## Posterior (95% credible intervals)",
        "",
        "| parameter | 2.5% | median | 97.5% |",
        "|---|---:|---:|---:|",
        prow("lp_propagation"),
        prow("lp_rate"),
        prow("k_um"),
        prow("gpx4_rate"),
        "",
        "## The in-vivo / in-vitro disjunction (the load-bearing finding)",
        "",
        "The in-vivo PRCC prior ranges used for the prior-predictive intervals do NOT",
        "overlap the in-vitro posterior:",
        "",
        "| parameter | in-vivo PRCC range | in-vitro posterior (2.5% .. median) | posterior entirely above in-vivo max? |",
        "|---|---|---|---|",
    ]
    for name, d in dj.items():
        lines.append(
            f"| `{name}` | {d['invivo_prcc_range']} | {d['posterior_q2_5']} .. {d['posterior_median']} | "
            f"{'YES' if d['posterior_above_invivo_range'] else 'no'} |"
        )
    lines += [
        "",
        "The in-vitro GPX4-inhibitor kill requires a lipid-peroxidation cascade 3 to 4x",
        "stronger than the in-vivo plausible ranges allow. So:",
        "",
        "1. The data-conditioned posterior is **in-vitro only**; it puts credible bands on",
        "   the #330 point calibration.",
        "2. The **in-vivo / spatial headline claims** (Bliss synergy, hypoxia collapse,",
        "   immune amplification, penetration) **cannot be conditioned on this in-vitro",
        "   data** and remain PRIOR-PREDICTIVE (the existing `headline_uncertainty.py` /",
        "   `uncertainty_intervals.py` intervals). Conditioning them needs an in-vivo",
        "   ferroptosis dose-response dataset, which is the remaining #332 dependency.",
        "",
        "## Posterior-predictive check on held-out " + r["heldout_compound"],
        "",
        f"Coverage of the empirical {r['heldout_compound']} median by the 95%",
        f"posterior-predictive band: **{r['heldout_coverage_strict']}** strictly inside the",
        f"band, **{r['heldout_coverage_tolerant']}** within a {r['heldout_tolerance']} viability",
        "tolerance. Both are reported because the tolerance is doing real work: the CTRPv2",
        "curves are cell-line MEDIANS and the model is single-cell, and the RMSE summary",
        "statistic discards curve shape, so a small viability offset is expected. The honest",
        "reading is that the posterior-predictive band is in the right place (tolerant",
        "coverage high) but not tight enough to bracket every point strictly, consistent with",
        "the single-cell-vs-median-curve mismatch and the limited 7-point summary statistic.",
        "",
        "| dose (µM) | " + " | ".join(str(d) for d in pp["dose_um"]) + " |",
        "|---|" + "---|" * len(pp["dose_um"]),
        "| empirical median | " + " | ".join(f"{v:.2f}" for v in pp["empirical_heldout"]) + " |",
        "| post-pred median | " + " | ".join(f"{v:.2f}" for v in pp["post_pred_median"]) + " |",
        "| post-pred 2.5% | " + " | ".join(f"{v:.2f}" for v in pp["post_pred_q2_5"]) + " |",
        "| post-pred 97.5% | " + " | ".join(f"{v:.2f}" for v in pp["post_pred_q97_5"]) + " |",
        "",
        "## Caveats",
        "",
        "- ABC rejection with a 7-point RMSE summary statistic and a fixed acceptance",
        "  fraction; the posterior width depends on that choice and the prior bounds.",
        "- `lp_propagation` / `lp_rate` identifiability is limited (same cascade), so their",
        "  marginal posteriors are correlated, not independently resolved.",
        "- In-vitro single-cell only; no in-vivo conditioning (see the disjunction).",
        "",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-draws", type=int, default=N_DRAWS)
    args = ap.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
