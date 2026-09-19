#!/usr/bin/env python3
"""Joint multi-inducer ABC calibration of the single-cell switch.

Fit shared cascade parameters to supported CTRPv2 ML162 and erastin fitted-curve
medians. Accept prior draws only within 1.10 times the reference vector's joint
RMSE. ML210 is held out by compound in the same screen, with overlapping cell
lines; it is not independent assay validation. Bands describe parameter draws,
not experimental replicate uncertainty.

Fewer than MIN_POSTERIOR accepted draws produce an underpowered diagnostic:
quantiles in JSON are retained for audit, but reports and plots suppress posterior
inference. The corrected 40,000-draw run accepted four draws. It cannot support
joint credible intervals or parameter-information claims. The spatial headlines
remain prior-predictive, with no validated transfer from this in-vitro target.

Run with the compiled ferroptosis_core extension:
  python scripts/abc_joint_posterior.py --n-draws 40000
Rebuild prose (and the underpowered status plot) without sampling:
  python scripts/abc_joint_posterior.py --render-only
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
ACCEPT_FRAC = 0.02  # retained only for the artifact's historical record

# Acceptance tolerance, as a multiple of a distance known to be reachable.
TOLERANCE_FACTOR = 1.10
# Below this the posterior quantiles are not worth reporting as a posterior.
MIN_POSTERIOR = 20
# The vector this repository already holds: the #330 CTRPv2 cascade fit plus
# #502's shared-switch erastin parameters, with the other two at Params::default.
# It is the BAR, never a sample -- it is not added to the accepted set.
REFERENCE_VECTOR = {"lp_propagation": 0.7, "lp_rate": 0.4, "gpx4_rate": 0.30,
                    "gsh_scav_efficiency": 0.5, "k_um": 0.25,
                    "k_erastin": 3.0, "hill": 6.0}
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
    curves, source = ck.load_target_data()
    rsl3_doses = list(ck.DOSE_GRID_UM)
    erastin_doses = list(ce.DOSE_GRID_UM)

    emp_rsl3, rsl3_support = ck.empirical_target(curves[ck.FIT_COMPOUND], rsl3_doses)
    emp_erastin, erastin_support = ck.empirical_target(curves[ce.COMPOUND], erastin_doses)
    emp_heldout, heldout_support = ck.empirical_target(curves[ck.HELDOUT_GPX4I], rsl3_doses)

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

    # ACCEPTANCE IS A TOLERANCE, NOT A QUANTILE.
    #
    # This used to be `n_accept = n_draws * 0.02` -- keep the best 2% whatever
    # they score -- which gave the run no floor. Epsilon was an OUTPUT, whatever
    # the last accepted draw happened to hit, and nothing had to meet it. The
    # consequence was measured: the returned posterior's median fitted WORSE than
    # a vector already committed in this repository (0.2413 against 0.2202) while
    # reporting epsilon 0.35, and 0 of 300 fresh uniform draws beat that vector.
    # The accepted set was a shell between the achievable and the arbitrary.
    # Diagnosis: `analysis/calibration/abc-acceptance-diagnostic.md`.
    #
    # So acceptance is now anchored to a distance that is demonstrably reachable:
    # the reference vector this repository already holds, scored by the same
    # distance in the same run. A draw is accepted only if it is within
    # TOLERANCE_FACTOR of that, which makes epsilon a criterion.
    #
    # Using a hand-tuned reference to gate an inference is worth being explicit
    # about. It does not feed the posterior -- the reference is never added to
    # the accepted set -- it only sets the bar. The alternative is a bar set by
    # whatever the sampler happened to draw, which is what produced the defect.
    reference = (ck.rmse(model_rsl3(rsl3_doses, REFERENCE_VECTOR), emp_rsl3)
                 + ck.rmse(model_erastin(erastin_doses, REFERENCE_VECTOR), emp_erastin))
    eps = float(reference * TOLERANCE_FACTOR)
    accept_idx = np.flatnonzero(distances <= eps)
    accept_idx = accept_idx[np.argsort(distances[accept_idx])]

    # An UNDER-POWERED run says so instead of padding itself back to a quota.
    # Silently topping up to 2% is exactly how a posterior of inferior vectors
    # came to be reported as a posterior.
    underpowered = len(accept_idx) < MIN_POSTERIOR
    if underpowered:
        print(f"UNDER-POWERED: {len(accept_idx)} of {args.n_draws} draws met the "
              f"tolerance {eps:.4f} (reference {reference:.4f} x {TOLERANCE_FACTOR}); "
              f"{MIN_POSTERIOR} are needed for a usable posterior. Reporting the "
              f"shortfall rather than widening the tolerance to fill the quota.",
              file=sys.stderr)
    if len(accept_idx) < 2:
        raise SystemExit(
            f"only {len(accept_idx)} draw(s) met the tolerance; nothing to "
            f"summarise. Increase --n-draws (the qualifying rate here is about "
            f"{len(accept_idx) / max(args.n_draws, 1):.2e} per draw).")
    posterior = draws[accept_idx]
    n_accept = len(accept_idx)

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

    # Which params are unconstrained -- judged against a NULL, not a constant.
    #
    # This was `>= 0.6`, a bare threshold that ignored how many draws were
    # accepted. That is the thing which decides what an uninformative posterior
    # looks like: with 30 accepted draws the 2.5-97.5 span of samples drawn from
    # the prior and nothing else still covers ~0.90 of the prior width, because
    # 30 points rarely reach the corners. So 0.6 sat far below anything noise
    # produces and the flag fired on WELL-constrained parameters -- it labelled
    # lp_propagation, lp_rate and gpx4_rate unconstrained while they sat at the
    # 0th percentile of that null. Full treatment in
    # `scripts/abc_posterior_information.py`.
    null = np.percentile(
        np.random.default_rng(0).uniform(0, 1, size=(20000, int(n_accept))),
        [2.5, 97.5], axis=1)
    null_widths = null[1] - null[0]
    null_p5 = float(np.percentile(null_widths, 5))
    unconstrained = [n for n in names
                     if post[n]["posterior_width_frac_of_prior"] > null_p5]

    result = {
        "target_source": source,
        "target_support": {ck.FIT_COMPOUND: rsl3_support, ce.COMPOUND: erastin_support,
                           ck.HELDOUT_GPX4I: heldout_support},
        "n_draws": args.n_draws,
        "n_accepted": int(n_accept),
        "accept_frac": round(n_accept / args.n_draws, 5),
        "acceptance_rule": "tolerance",
        "reference_distance": round(float(reference), 4),
        "tolerance_factor": TOLERANCE_FACTOR,
        "underpowered": bool(underpowered),
        "min_posterior": MIN_POSTERIOR,
        "uninformative_null_p5_width": round(null_p5, 3),
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
    if posterior_is_usable(result):
        print(f"held-out ML210 coverage {result['heldout_posterior_predictive']['coverage_inside_95pct_band']}, "
              f"median PP RMSE {result['heldout_posterior_predictive']['median_pp_rmse']}")
        print(f"unconstrained: {unconstrained}")
    else:
        print("UNDERPOWERED: joint intervals, information classifications, and "
              "generalization claims withheld; JSON summaries are audit-only diagnostics.")
    print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)} + {OUT_MD.relative_to(REPO_ROOT)} + {OUT_PNG.name}")
    return result


def posterior_is_usable(r):
    """Require the sample-size guard even if an input flag is stale."""
    return not r["underpowered"] and r["n_accepted"] >= r["min_posterior"]


def _plot_underpowered(r):
    """Render targets and sampling status without implying credible bands."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = r["curves"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    panels = (
        ("ML162: fit target", c["rsl3_doses_um"], c["empirical_rsl3_ml162"]),
        ("Erastin: fit target", c["erastin_doses_um"], c["empirical_erastin"]),
        ("ML210: held-out compound", c["rsl3_doses_um"],
         r["heldout_posterior_predictive"]["empirical"]),
    )
    for ax, (title, doses, target) in zip(axes, panels):
        ax.semilogx(doses, target, "ko-", label="Supported fitted-curve median")
        ax.set(xlabel="Dose (µM)", ylabel="Viability", ylim=(0, 1.1), title=title)
        ax.legend(fontsize=7)
    fig.suptitle(
        f"UNDERPOWERED: {r['n_accepted']} / {r['n_draws']:,} draws accepted; "
        f"minimum {r['min_posterior']} required\n"
        "Joint credible intervals and predictive bands withheld", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=120)
    plt.close(fig)


def _plot(r, accept_idx, rsl3_models, erastin_models, pp_band):
    if not posterior_is_usable(r):
        _plot_underpowered(r)
        return
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

    ax[0].fill_between(rsl3_doses, r_band[0], r_band[2], alpha=0.25, color="C0", label="95% parameter-draw band")
    ax[0].plot(rsl3_doses, r_band[1], "-", color="C0", label="posterior median")
    ax[0].plot(rsl3_doses, c["empirical_rsl3_ml162"], "ko", label="ML162 (fit)")
    ax[0].set_title("RSL3 / GPX4i (ML162): fit")
    ax[0].legend(fontsize=7)

    ax[1].fill_between(erastin_doses, e_band[0], e_band[2], alpha=0.25, color="C1", label="95% parameter-draw band")
    ax[1].plot(erastin_doses, e_band[1], "-", color="C1", label="posterior median")
    ax[1].plot(erastin_doses, c["empirical_erastin"], "ks", label="erastin (fit)")
    ax[1].set_title("erastin / System Xc-: fit")
    ax[1].legend(fontsize=7)

    ax[2].fill_between(rsl3_doses, pp_band[0], pp_band[2], alpha=0.25, color="C2", label="95% parameter-draw band")
    ax[2].plot(rsl3_doses, pp_band[1], "-", color="C2", label="predictive median")
    ax[2].plot(rsl3_doses, hp["empirical"], "k^", label=f"{hp['compound']} (held-out)")
    ax[2].set_title(f"{hp['compound']}: held-out posterior-predictive")
    ax[2].legend(fontsize=7)

    fig.suptitle("Joint multi-inducer in-vitro posterior (#500): RSL3 + erastin", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=120)
    plt.close(fig)


def _write_underpowered_report(r):
    md = f"""# Joint multi-inducer calibration: underpowered

Generated by `scripts/abc_joint_posterior.py` from `joint-posterior.json`.

**Underpowered:** {r['n_accepted']} of {r['n_draws']:,} prior draws met the
acceptance criterion; the minimum is {r['min_posterior']}. Joint credible
intervals, parameter-information classifications, and held-out generalization
claims are withheld. The JSON quantiles and coverage are diagnostic summaries
of this small accepted set, not a usable posterior.

## Dose support

{ck.support_markdown(r['target_support'])}

The unsupported 100 µM erastin target has been removed. The corrected erastin
grid ends at 30 µM. Priors, reference vector, tolerance factor, and minimum
accepted count are unchanged.

## Acceptance and history

The joint distance is ML162 RMSE plus erastin RMSE. The reference distance
**{r['reference_distance']}** times **{r['tolerance_factor']}** sets epsilon to
**{r['epsilon_joint_distance']}**. The accepted fraction is {r['accept_frac']:.2%};
that fraction is an outcome, not a quota.

An earlier version accepted a fixed 2% of draws. Its reference distance 0.2202
and median-vector distance 0.2413 used the original targets, including the
unsupported 100 µM value; they are not directly comparable to this run.
Acceptance is now a TOLERANCE. The former fixed 0.6 threshold is no longer used
for the information-width check. For this small set, its sample-size null P5
width is **{r['uninformative_null_p5_width']}**, retained only as a diagnostic;
no parameter is classified as informed or unconstrained here. See
`abc-acceptance-diagnostic.md` and `abc-information-content.md`.

## What can be concluded

This run did not establish a usable joint posterior under the documented
sampling budget and acceptance rule. It does not distinguish insufficient
sampling from model or target mismatch. More efficient sampling should be
assessed against the same fixed target and criterion before interpreting
intervals. A larger draw count alone does not establish biological validity.

The plot shows only supported fitted-curve targets and this sampling status.
A parameter-draw band omits uncertainty in experimental outcomes and fitted
curves. Dose-wise targets are not independent observations; ML210 shares cell
lines with the training compound within the same screen. Independent assay
validation remains pending.

The separately refitted GPX4 point fit and single-inducer ABC remain available
in `kill-switch-calibration.md` and `abc-posterior-report.md`. They do not supply
a joint posterior. In-vivo/spatial headline intervals remain prior-predictive;
this run establishes no transfer of in-vitro parameters to those outputs.
"""
    OUT_MD.write_text(md, encoding="utf-8")


def write_report(r):
    if not posterior_is_usable(r):
        _write_underpowered_report(r)
        return
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
    separation = (
        "Both reported 95% intervals lie above the corresponding PRCC ranges."
        if all(v["entire_95pct_posterior_above_invivo_max"] for v in d.values())
        else "At least one reported 95% interval does not lie entirely above its PRCC range."
    )
    sampling_note = (
        f"**Underpowered:** only {r['n_accepted']} draws met the criterion, below "
        f"the {r['min_posterior']}-draw reporting threshold. These quantiles are unstable."
        if r["underpowered"] else
        f"The accepted set meets the {r['min_posterior']}-draw minimum; this is a "
        "sampling check, not evidence of biological validity."
    )

    md = f"""# Joint multi-inducer in-vitro posterior (#500)

Generated by `scripts/abc_joint_posterior.py` (needs the compiled `ferroptosis_core`
extension; not run in CI). Builds on the #330 GPX4i fit, the #502 System Xc-/erastin
mechanism, and the #332 single-inducer ABC.

## Dose support

{ck.support_markdown(r['target_support'])}

Erastin uses the original grid through 30 µM. The unsupported 100 µM point
has been removed. The reference distance and acceptance threshold are recomputed
on these corrected targets, with the same priors and tolerance factor.

## History of this run's acceptance rule

An earlier version accepted a fixed 2% of draws, so epsilon was an OUTPUT --
whatever the last accepted draw scored -- and nothing had to meet it. That run
(1,500 draws, 30 accepted, epsilon 0.35) returned a posterior whose median fitted
WORSE than a vector already committed in this repository: 0.2413 against 0.2202
on its original target set, which included the unsupported 100 µM point.
Those historical distances are not directly comparable with distances on the
corrected targets above. Acceptance is now a TOLERANCE anchored to the
reachable reference, so epsilon is a criterion; the run reports a shortfall
rather than padding itself back to a quota. History and current-target diagnostics:
`analysis/calibration/abc-acceptance-diagnostic.md`.

This note lives in the generator rather than in the document, because the
corrections written directly into the generated markdown were erased the next
time it was generated -- which is the obvious failure in hindsight and was not
obvious at the time.

## What this is

The shared single-cell switch (`lp_propagation`, `lp_rate`, `gpx4_rate`,
`gsh_scav_efficiency`) is conditioned JOINTLY on **two inducer mechanisms at once**
(a GPX4 inhibitor, ML162, and a System Xc- inhibitor, erastin), with per-mechanism
potencies (`k_um` for RSL3; `k_erastin`, `hill` for erastin). ABC over {r['n_draws']}
in-vitro-spanning prior draws. A draw is accepted when its joint distance
(ML162 RMSE + erastin RMSE) is at most **{r['epsilon_joint_distance']}**:
the reference distance **{r['reference_distance']}** multiplied by
**{r['tolerance_factor']}**. This accepted **{r['n_accepted']} draws**
({r['accept_frac']:.2%}); that fraction is an outcome, not a quota.

{sampling_note}

## Joint posterior (credible intervals, not a point)

{ptab()}

The `width (frac of prior)` column measures marginal interval width relative
to the prior range. Its reference is the width obtained by sampling the prior
alone with the same number of accepted draws. The 5th percentile of that null
is **{r['uninformative_null_p5_width']}**. Parameters whose width exceeds this
threshold show no detected contraction by this criterion:
**{", ".join(f"`{x}`" for x in r['unconstrained_params']) or "none"}**.
The former fixed 0.6 threshold is no longer used. Marginal contraction measures
information under this ABC design; it does not establish unique mechanistic
identification or a good fit.

## Held-out posterior-predictive ({hp['compound']}, never used in the fit)

The accepted parameter draws produce model viability curves compared with the
held-out compound {hp['compound']}: **{hp['coverage_inside_95pct_band']} fitted-curve
target values** fall inside the central 95% parameter-draw band, with median
held-out RMSE **{hp['median_pp_rmse']}**. Simulation seed and population size are
fixed. The band omits uncertainty in the fitted CTRPv2 curves, variation across
cell lines, and a measurement-error model; it is not a calibrated 95% interval
for experimental outcomes. The dose-wise checks share fitted curves and are
not independent validation observations. ML210 was held out by compound within
the same screen, with overlapping cell lines. See
`joint-posterior-predictive.png` (right panel).

## The load-bearing caveat: in-vitro only, disjoint from the in-vivo priors

The joint posterior is **in-vitro**. Re-checking the #332 disjunction for the joint
fit:

{disj_lines}

{separation} These PRCC ranges are sensitivity bands around chosen defaults;
separation from them is not independent validation of an in-vivo regime.
The calibration scope remains:

- **In-vitro switch claims** carry these joint-posterior credible intervals.
- **In-vivo / spatial headline numbers** (hypoxia asymmetry, Bliss synergy,
  penetration gap, immune ratio) **cannot** be conditioned on in-vitro data and
  stay **prior-predictive** (the existing `headline_uncertainty.py` /
  `uncertainty_intervals.py` intervals, `analysis/identifiability-report.md`).
  The current workflow has no validated transfer from this fitted in-vitro
  observable to those spatial headline parameters.
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-draws", type=int, default=N_DRAWS)
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the markdown from the committed JSON without "
                         "re-running the sampling (a 40k-draw run is ~14 minutes)")
    args = ap.parse_args()
    if args.render_only:
        # The prose lives in this file, so a wording change should not require
        # re-running the inference. Precedent: scripts/atlas_ambiguity.py.
        result = json.loads(OUT_JSON.read_text())
        write_report(result)
        if not posterior_is_usable(result):
            _plot_underpowered(result)
        print(f"re-rendered {OUT_MD.relative_to(REPO_ROOT)} from "
              f"{OUT_JSON.relative_to(REPO_ROOT)} (no sampling)")
        return 0
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
