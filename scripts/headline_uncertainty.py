#!/usr/bin/env python3
"""Prior-predictive uncertainty intervals for the Bliss and sim-tme headlines (#332).

`scripts/uncertainty_intervals.py` reports prior-predictive credible intervals for
the SINGLE-CELL kill switch (via `sim_batch`). Its docstring names the tractable
follow-up: propagate the SAME uniform priors through the SPATIAL/COMBO headline
outputs that the simulation BINARIES produce. This script does that, sampling the
joint uniform prior over the same 11 PRCC rate constants and running a binary
under the `FERRO_PARAM_OVERRIDES` hook (#331) for each draw, for two headline
families:

- `--headline bliss` (default): the RSL3 + FSP1i Bliss synergy from
  `sim-combo-mech` (the manuscript's ~1.99x). Cheap (~seconds/run), so a full
  Monte-Carlo is feasible (default 300 draws).
- `--headline sim-tme`: the two spatial `sim-tme` headlines extracted from one
  run per draw — the hypoxia kill-collapse GAP (SDT minus RSL3 hypoxic-zone kill)
  and SDT's pool-de-confounded immune kill rate. `sim-tme` costs ~4 min/run, so
  this uses a smaller default ensemble (100 draws) and the 2.5/97.5 TAILS are
  read cautiously (median + spread, not the exact bounds).

It reuses the binary-invocation + parameter-range machinery from
`headline_sensitivity.py` (`run_bliss`, `run_sim_tme_observables`, `PARAM_NAMES`,
`LOWS`, `HIGHS`, `_default_binary`) so the priors, override mapping, and
observables are identical to the Morris sensitivity screen.

This is PRIOR-PREDICTIVE (forward) uncertainty over UNIFORM priors, NOT a
data-conditioned Bayesian/ABC posterior — turning the priors into a posterior
needs the external GDSC/DepMap drug-response data (#330). It captures PARAMETER
uncertainty only, not STRUCTURAL uncertainty. The penetration headline
(`drug_transport` Krogh) is the remaining deferred extension.

Usage:
    python scripts/headline_uncertainty.py [--headline bliss] [--samples 300] [--workers 8]
    python scripts/headline_uncertainty.py --headline sim-tme [--tme-samples 100] [--workers 8]

Deterministic given the sample count + the fixed seed. Writes
analysis/headline-uncertainty-report.md (bliss) or
analysis/headline-uncertainty-tme-report.md (sim-tme).
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from headline_sensitivity import (  # noqa: E402  (path insert above)
    HIGHS,
    LOWS,
    PARAM_NAMES,
    PENETRATION_TISSUES,
    _default_binary,
    extract_tissue_pk_observables,
    extract_tme_observables,
    read_bliss_synergy,
    run_bliss,
    run_sim_tme_observables,
    run_tissue_pk_observables,
)

REPORT = REPO / "analysis" / "headline-uncertainty-report.md"
TME_REPORT = REPO / "analysis" / "headline-uncertainty-tme-report.md"
PEN_REPORT = REPO / "analysis" / "headline-uncertainty-penetration-report.md"
SEED = 332  # deterministic prior draws (issue number; fixed for reproducibility)
DEFAULT_SAMPLES = 300
# sim-tme is far costlier than sim-combo-mech (~4 min/run single-threaded), so its
# prior-predictive sample is much smaller; the wider-tail caveat is documented.
DEFAULT_TME_SAMPLES = 100
# sim-tissue-pk is cheap (~4 s/run), so the penetration headline gets a full sample.
DEFAULT_PEN_SAMPLES = 300


def sample_prior(n_samples, seed=SEED):
    """n_samples joint-uniform draws over the 11 PRCC parameter ranges (one row
    per draw, columns in PARAM_NAMES order). Deterministic given (n_samples, seed)."""
    rng = np.random.default_rng(seed)
    return rng.uniform(LOWS, HIGHS, size=(n_samples, len(PARAM_NAMES)))


def prior_predictive_bliss(n_samples, workers, binary):
    """Run sim-combo-mech under each prior draw and partition the outcomes into
    three categories, returning (finite_synergy_array, n_failed, n_undefined):

    - finite synergy: a defined `synergy_score`, summarized into the interval;
    - n_undefined: the binary ran (exit 0) but emitted a NON-FINITE synergy.
      sim-combo-mech writes NaN when the Bliss baseline `<= 0.001`, i.e. both
      single agents kill ~0% (a ferroptosis-resistant corner of the prior), so
      the synergy RATIO is mathematically undefined. These must be dropped from
      the percentiles (a single NaN poisons `np.percentile`) and counted
      separately rather than silently passed through as "successes";
    - n_failed: the subprocess itself errored.
    """
    draws = sample_prior(n_samples)

    def _one(row):
        try:
            return run_bliss(row, binary)  # may be a finite float or NaN
        except Exception:
            return None  # subprocess / parse failure

    with ThreadPoolExecutor(max_workers=workers) as ex:  # subprocess releases the GIL
        results = list(ex.map(_one, draws))
    return _partition(results)


def _partition(results):
    """Split run outcomes into (finite_synergy_array, n_failed, n_undefined):
    `None` = subprocess/parse failure; a non-finite float (NaN/inf) = undefined
    Bliss (binary ran but the ratio is undefined); a finite float = a valid
    synergy. Pure, so the NaN-dropping contract is unit-tested without the
    binary — the bug this guards is a single NaN poisoning every percentile."""
    finite = [v for v in results if v is not None and np.isfinite(v)]
    n_failed = sum(1 for v in results if v is None)
    n_undefined = sum(1 for v in results if v is not None and not np.isfinite(v))
    return np.array(finite, dtype=float), n_failed, n_undefined


def _pctiles(values):
    return {
        "n": int(values.size),
        "p2_5": float(np.percentile(values, 2.5)),
        "median": float(np.percentile(values, 50)),
        "p97_5": float(np.percentile(values, 97.5)),
        "min": float(values.min()),
        "max": float(values.max()),
        "mean": float(values.mean()),
    }


# ----------------------------------------------------------------------------
# sim-tme headlines (hypoxia kill-collapse gap + de-confounded immune rate)
# ----------------------------------------------------------------------------
def prior_predictive_tme(n_samples, workers, binary):
    """Run sim-tme under each prior draw and extract BOTH headline observables
    from the one (costly) run, via `headline_sensitivity.run_sim_tme_observables`:
    `hypoxia` (the SDT-minus-RSL3 hypoxic-zone kill GAP) and `immune` (SDT's
    pool-de-confounded immune kill rate, bounded [0,1]). Returns
    `(hyp_array, imm_array, n_failed)`. A draw whose sim-tme run errors is dropped
    (counted in n_failed); a non-finite observable is also dropped defensively
    (the Bliss-review lesson — one NaN poisons np.percentile), though the
    de-confounded immune rate is floored to avoid div-by-zero at the source."""
    draws = sample_prior(n_samples)

    def _one(row):
        try:
            return run_sim_tme_observables(row, binary)  # {"hypoxia":.., "immune":..}
        except Exception:
            return None

    with ThreadPoolExecutor(max_workers=workers) as ex:  # subprocess releases the GIL
        results = list(ex.map(_one, draws))
    return _partition_tme(results)


def _partition_tme(results):
    """Split sim-tme outcomes into (hyp_finite_array, imm_finite_array, n_failed).
    Each element of `results` is a `{"hypoxia":float, "immune":float}` dict or
    `None` (run failure). Per-observable non-finite values are dropped defensively
    (the Bliss-review lesson — one NaN poisons np.percentile), independently for
    each observable. Pure, so the NaN/failure handling is unit-tested without the
    costly binary."""
    n_failed = sum(1 for r in results if r is None)
    hyp = np.array(
        [r["hypoxia"] for r in results if r is not None and np.isfinite(r["hypoxia"])],
        dtype=float,
    )
    imm = np.array(
        [r["immune"] for r in results if r is not None and np.isfinite(r["immune"])],
        dtype=float,
    )
    return hyp, imm, n_failed


def _default_tme(binary):
    """The unperturbed sim-tme observables (no FERRO_PARAM_OVERRIDES) — the
    manuscript point estimates the intervals bracket."""
    with tempfile.TemporaryDirectory(prefix="ferro_pp_tme_default_") as workdir:
        env = dict(os.environ)
        env.pop("FERRO_PARAM_OVERRIDES", None)
        proc = subprocess.run([str(binary)], cwd=workdir, env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"default sim-tme failed: {proc.stderr[-300:]}")
        summary = Path(workdir) / "output" / "tme" / "tme_summary.json"
        conditions = json.loads(summary.read_text())["conditions"]
        return extract_tme_observables(conditions)


def write_tme_report(hyp_stats, imm_stats, default_obs, n_failed, n_samples):
    lines = []
    lines.append("# Prior-predictive uncertainty: sim-tme headlines (#332)\n")
    lines.append(
        "Generated by `scripts/headline_uncertainty.py --headline sim-tme`. Propagates "
        "the same joint UNIFORM prior over the 11 PRCC rate constants (the ranges the "
        "PRCC #134 / Sobol #331 / Morris screen use) through the two spatial `sim-tme` "
        "headlines, extracted from one (costly, ~4 min) run per draw via "
        "`headline_sensitivity.run_sim_tme_observables`:\n"
        "- **hypoxia**: the SDT-minus-RSL3 hypoxic-zone kill GAP (the kill-collapse "
        "asymmetry — SDT holds, RSL3 collapses);\n"
        "- **immune**: SDT's pool-de-confounded immune kill rate "
        "`immune_kills / (total_tumor - ferroptosis_kills)` (bounded [0,1]).\n"
    )
    lines.append("## Method\n")
    lines.append(
        f"- {n_samples} Monte-Carlo draws from the joint uniform prior (seed {SEED}, "
        f"deterministic); {hyp_stats['n']} runs informed the hypoxia interval, "
        f"{imm_stats['n']} the immune interval"
        + (f"; {n_failed} draws dropped (sim-tme run failure)" if n_failed else "")
        + ".\n"
        "- **Small-sample caveat:** sim-tme costs ~4 min/run, so this uses far fewer "
        f"draws ({n_samples}) than the single-cell ensemble (2000); the 2.5/97.5 TAILS are "
        "correspondingly less stable than the single-cell or Bliss intervals — read the "
        "median + the broad spread, not the exact bounds.\n"
        "- **Prior-predictive (forward) uncertainty over UNIFORM priors — NOT a "
        "data-conditioned posterior** (blocked on the #330 GDSC data); PARAMETER, not "
        "STRUCTURAL, uncertainty.\n"
    )
    lines.append("## Hypoxia kill-collapse gap (SDT − RSL3 hypoxic-zone kill)\n")
    lines.append("| quantity | value |")
    lines.append("|---|---|")
    lines.append(f"| default point estimate | {default_obs['hypoxia']:.3f} |")
    lines.append(f"| prior-predictive median | {hyp_stats['median']:.3f} |")
    lines.append(f"| 95% prior-predictive interval | [{hyp_stats['p2_5']:.3f}, {hyp_stats['p97_5']:.3f}] |")
    lines.append(f"| full range (min, max) | [{hyp_stats['min']:.3f}, {hyp_stats['max']:.3f}] |")
    lines.append("")
    lines.append("## Immune de-confounded kill rate (SDT, immune on)\n")
    lines.append("| quantity | value |")
    lines.append("|---|---|")
    lines.append(f"| default point estimate | {default_obs['immune']:.3f} |")
    lines.append(f"| prior-predictive median | {imm_stats['median']:.3f} |")
    lines.append(f"| 95% prior-predictive interval | [{imm_stats['p2_5']:.3f}, {imm_stats['p97_5']:.3f}] |")
    lines.append(f"| full range (min, max) | [{imm_stats['min']:.3f}, {imm_stats['max']:.3f}] |")
    lines.append("")
    lines.append("## Scope\n")
    lines.append(
        "Completes the spatial-headline prior-predictive started with the Bliss synergy "
        "(`analysis/headline-uncertainty-report.md`); the single-cell kill-switch intervals "
        "are in `analysis/uncertainty-intervals-report.md`. The data-conditioned posterior "
        "(ABC) for all of these is blocked on #330. The penetration headline (`drug_transport` "
        "Krogh) is a pure-Python ensemble, a remaining smaller extension.\n"
    )
    TME_REPORT.write_text("\n".join(lines), encoding="utf-8")


# ----------------------------------------------------------------------------
# Penetration headline (per-tissue vessel-wall RSL3 kill, sim-tissue-pk)
# ----------------------------------------------------------------------------
def prior_predictive_penetration(n_samples, workers, binary):
    """Run sim-tissue-pk under each prior draw and collect the RSL3-like
    vessel-wall death rate for each of the three tissue scenarios via
    `headline_sensitivity.run_tissue_pk_observables`. Returns
    `({tissue_key: finite_array}, n_failed)`. The transport params are held at
    their fixed (uncalibrated) per-tissue presets; the perturbed 11 PRCC biochem
    rates drive the kill conversion, so this is the SAME prior as the other
    headlines applied at each tissue's fixed transport scenario."""
    draws = sample_prior(n_samples)

    def _one(row):
        try:
            return run_tissue_pk_observables(row, binary)  # {tissue_key: dr}
        except Exception:
            return None

    with ThreadPoolExecutor(max_workers=workers) as ex:  # subprocess releases the GIL
        results = list(ex.map(_one, draws))
    return _partition_penetration(results)


def _partition_penetration(results):
    """Split per-tissue penetration outcomes into ({tissue_key: finite_array},
    n_failed). Each element of `results` is a `{tissue_key: death_rate}` dict or
    `None` (run failure). Per-tissue non-finite values are dropped defensively
    (the Bliss-review NaN lesson). Pure, so unit-tested without the binary."""
    keys = [k for k, _ in PENETRATION_TISSUES]
    n_failed = sum(1 for r in results if r is None)
    per_tissue = {
        k: np.array(
            [r[k] for r in results if r is not None and np.isfinite(r.get(k, float("nan")))],
            dtype=float,
        )
        for k in keys
    }
    return per_tissue, n_failed


def _default_penetration(binary):
    """The unperturbed per-tissue vessel-wall death rates (no FERRO_PARAM_OVERRIDES)
    — the manuscript point estimates (12.1% / 2.6% / 1.8%) the intervals bracket."""
    with tempfile.TemporaryDirectory(prefix="ferro_pp_pen_default_") as workdir:
        env = dict(os.environ)
        env.pop("FERRO_PARAM_OVERRIDES", None)
        proc = subprocess.run([str(binary)], cwd=workdir, env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"default sim-tissue-pk failed: {proc.stderr[-300:]}")
        summary = Path(workdir) / "output" / "tissue-pk" / "tissue_pk_summary.json"
        return extract_tissue_pk_observables(json.loads(summary.read_text()))


_TISSUE_LABELS = {
    "well_vascularized": "Well-vascularized epithelial",
    "poorly_vascularized": "Poorly-vascularized (pancreatic)",
    "cns_bbb": "Neuroectodermal (CNS / BBB)",
}


def write_penetration_report(per_tissue_stats, default_obs, n_failed, n_samples):
    lines = []
    lines.append("# Prior-predictive uncertainty: drug-penetration headline (#332)\n")
    lines.append(
        "Generated by `scripts/headline_uncertainty.py --headline penetration`. Propagates "
        "the same joint UNIFORM prior over the 11 PRCC rate constants through the "
        "tissue-specific vessel-wall RSL3 kill (`sim-tissue-pk`, the manuscript's "
        "40% in 2D culture dropping to 12.1% / 2.6% / 1.8% at the vessel wall across three "
        "tissue scenarios). Each radial-bin kill is `sim_cell` under RSL3 with GPX4 "
        "inhibition scaled by the local Krogh concentration, so the 11 biochem rates drive "
        "the kill while the transport params (penetration length, vascular permeability, "
        "ECM tortuosity) are held at their FIXED, uncalibrated per-tissue presets.\n"
    )
    lines.append("## Method\n")
    lines.append(
        f"- {n_samples} Monte-Carlo draws from the joint uniform prior (seed {SEED}, "
        f"deterministic)"
        + (f"; {n_failed} draws dropped (sim-tissue-pk run failure)" if n_failed else "")
        + ".\n"
        "- The interval here is the BIOCHEM-rate uncertainty at each tissue's fixed "
        "transport scenario; the transport parameters themselves are uncalibrated "
        "placeholders (their own prior is NOT swept — there is no PRCC-style range table "
        "for them — so the across-tissue ORDERING, not the absolute per-tissue magnitude, "
        "is the robust reading).\n"
        "- **Prior-predictive (forward) uncertainty over UNIFORM priors — NOT a "
        "data-conditioned posterior** (blocked on the #330 GDSC data); PARAMETER, not "
        "STRUCTURAL, uncertainty.\n"
    )
    lines.append("## Vessel-wall RSL3 death rate, by tissue\n")
    lines.append("| tissue | default point | median | 95% prior-predictive interval | full range |")
    lines.append("|---|---|---|---|---|")
    for key, _ in PENETRATION_TISSUES:
        s = per_tissue_stats[key]
        lines.append(
            f"| {_TISSUE_LABELS[key]} | {default_obs[key]:.3f} | {s['median']:.3f} | "
            f"[{s['p2_5']:.3f}, {s['p97_5']:.3f}] | [{s['min']:.3f}, {s['max']:.3f}] |"
        )
    lines.append("")
    lines.append("## Scope\n")
    lines.append(
        "Completes the prior-predictive propagation of #332's headline claims over the 11 "
        "PRCC biochem rate constants (single-cell kill switch, Bliss synergy, the two "
        "sim-tme headlines, and now penetration). The remaining piece is the "
        "data-conditioned ABC posterior, blocked on the external GDSC/DepMap/PRISM "
        "drug-response data (#330); turning these uniform priors into a posterior, and "
        "sweeping the transport parameters' own (currently undocumented) ranges, are the "
        "open extensions.\n"
    )
    PEN_REPORT.write_text("\n".join(lines), encoding="utf-8")


def write_report(stats, default_synergy, n_failed, n_undefined, n_samples):
    lines = []
    lines.append("# Prior-predictive uncertainty: Bliss synergy headline (#332)\n")
    lines.append(
        "Generated by `scripts/headline_uncertainty.py`. Propagates a joint UNIFORM "
        "prior over the 11 PRCC rate constants (`analysis/prcc-results.json`, the same "
        "ranges the PRCC #134 / Sobol #331 / Morris screen use) through the RSL3 + FSP1i "
        "**Bliss synergy** headline (`synergy_score` from `sim-combo-mech`, the "
        "manuscript's ~1.99x), via the `FERRO_PARAM_OVERRIDES` hook. Reuses "
        "`headline_sensitivity.run_bliss` so the prior + observable match the sensitivity "
        "screen exactly.\n"
    )
    lines.append("## Method\n")
    lines.append(
        f"- {n_samples} Monte-Carlo draws from the joint uniform prior (seed {SEED}, "
        f"deterministic); {stats['n']} yielded a defined synergy and form the interval, "
        f"{n_undefined} undefined-Bliss (dropped), {n_failed} binary failures (dropped).\n"
        "- **Undefined-Bliss draws** are a ferroptosis-resistant corner of the prior where "
        "BOTH single agents kill ~0%, so the Bliss baseline is ~0 and the synergy ratio "
        "(death / Bliss) is mathematically undefined; sim-combo-mech emits NaN there "
        "(`bliss <= 0.001`). These are dropped from the percentiles (one NaN would poison "
        "`np.percentile`) and counted separately rather than passed through as successes — "
        "their existence is itself an honest note on the headline's robustness.\n"
        "- **Prior-predictive (forward) uncertainty over UNIFORM priors — NOT a "
        "data-conditioned Bayesian/ABC posterior.** Turning the priors into a posterior "
        "needs the external GDSC/DepMap/PRISM drug-response data (#330). Captures "
        "PARAMETER uncertainty only, not STRUCTURAL uncertainty (model form, the fixed "
        "RSL3 `DrugEffect`, the Bliss independence assumption).\n"
    )
    lines.append("## Result\n")
    lines.append("| quantity | RSL3 + FSP1i synergy_score |")
    lines.append("|---|---|")
    lines.append(f"| default point estimate (no perturbation) | {default_synergy:.3f} |")
    lines.append(f"| prior-predictive median | {stats['median']:.3f} |")
    lines.append(f"| 95% prior-predictive interval | [{stats['p2_5']:.3f}, {stats['p97_5']:.3f}] |")
    lines.append(f"| full range (min, max) | [{stats['min']:.3f}, {stats['max']:.3f}] |")
    lines.append("")
    lines.append(
        f"The headline Bliss synergy is reported in the manuscript as a point estimate "
        f"(~1.99x; the default run here gives {default_synergy:.2f}x). Under prior "
        f"uncertainty over the 11 rate constants it spans a 95% interval of "
        f"**[{stats['p2_5']:.2f}x, {stats['p97_5']:.2f}x]** (median {stats['median']:.2f}x). "
        "A synergy > 1 means the RSL3 + FSP1i combination kills more than the Bliss "
        "independence prediction; the interval shows whether that synergy conclusion is "
        "robust to the uncertain rate constants or driven by the specific placeholder "
        "values.\n"
    )
    lines.append("## Scope\n")
    lines.append(
        "Bliss is the first headline propagated because `sim-combo-mech` is the cheapest "
        "binary. The costlier `sim-tme` headlines (hypoxia kill-collapse, the immune "
        "ratio) are the next extension of this prior-predictive harness; the single-cell "
        "kill-switch intervals are in `analysis/uncertainty-intervals-report.md`. The "
        "data-conditioned posterior (ABC) for all of these is blocked on #330.\n"
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def _run_bliss(args):
    binary = _default_binary("sim-combo-mech")
    if binary is None:
        raise SystemExit(
            "sim-combo-mech binary not found; build it: "
            "cd simulations && cargo build --release -p sim-combo-mech"
        )
    # Default point estimate: the no-override (FERRO_PARAM_OVERRIDES-unset) run,
    # i.e. the manuscript ~1.99x the interval brackets.
    default_synergy = _default_bliss(binary)
    values, n_failed, n_undefined = prior_predictive_bliss(args.samples, args.workers, binary)
    if values.size == 0:
        raise SystemExit("no draw produced a defined synergy; nothing to report")
    stats = _pctiles(values)
    write_report(stats, default_synergy, n_failed, n_undefined, args.samples)
    print(
        f"Bliss prior-predictive: default={default_synergy:.3f}  "
        f"median={stats['median']:.3f}  95% CI=[{stats['p2_5']:.3f}, {stats['p97_5']:.3f}]  "
        f"(n={stats['n']}, {n_undefined} undefined-Bliss, {n_failed} failed)"
    )
    print(f"wrote {REPORT.relative_to(REPO)}")


def _run_tme(args):
    binary = _default_binary("sim-tme")
    if binary is None:
        raise SystemExit(
            "sim-tme binary not found; build it: "
            "cd simulations && cargo build --release -p sim-tme"
        )
    n = args.tme_samples
    default_obs = _default_tme(binary)
    hyp, imm, n_failed = prior_predictive_tme(n, args.workers, binary)
    if hyp.size == 0 or imm.size == 0:
        raise SystemExit("no sim-tme draw produced a usable observable; nothing to report")
    hyp_stats, imm_stats = _pctiles(hyp), _pctiles(imm)
    write_tme_report(hyp_stats, imm_stats, default_obs, n_failed, n)
    print(
        f"sim-tme prior-predictive (n={n}, {n_failed} failed):\n"
        f"  hypoxia gap : default={default_obs['hypoxia']:.3f}  median={hyp_stats['median']:.3f}  "
        f"95% CI=[{hyp_stats['p2_5']:.3f}, {hyp_stats['p97_5']:.3f}]\n"
        f"  immune rate : default={default_obs['immune']:.3f}  median={imm_stats['median']:.3f}  "
        f"95% CI=[{imm_stats['p2_5']:.3f}, {imm_stats['p97_5']:.3f}]"
    )
    print(f"wrote {TME_REPORT.relative_to(REPO)}")


def _run_penetration(args):
    binary = _default_binary("sim-tissue-pk")
    if binary is None:
        raise SystemExit(
            "sim-tissue-pk binary not found; build it: "
            "cd simulations && cargo build --release -p sim-tissue-pk"
        )
    n = args.pen_samples
    default_obs = _default_penetration(binary)
    per_tissue, n_failed = prior_predictive_penetration(n, args.workers, binary)
    if any(per_tissue[k].size == 0 for k, _ in PENETRATION_TISSUES):
        raise SystemExit("a tissue scenario produced no usable draw; nothing to report")
    per_tissue_stats = {k: _pctiles(per_tissue[k]) for k, _ in PENETRATION_TISSUES}
    write_penetration_report(per_tissue_stats, default_obs, n_failed, n)
    print(f"penetration prior-predictive (n={n}, {n_failed} failed):")
    for key, _ in PENETRATION_TISSUES:
        s = per_tissue_stats[key]
        print(
            f"  {key:20s}: default={default_obs[key]:.3f}  median={s['median']:.3f}  "
            f"95% CI=[{s['p2_5']:.3f}, {s['p97_5']:.3f}]"
        )
    print(f"wrote {PEN_REPORT.relative_to(REPO)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--headline",
        choices=["bliss", "sim-tme", "penetration"],
        default="bliss",
        help="bliss = RSL3+FSP1i synergy (fast); sim-tme = hypoxia gap + immune rate "
        "(~4 min/run); penetration = per-tissue vessel-wall RSL3 kill (sim-tissue-pk, ~4 s/run)",
    )
    ap.add_argument("--samples", type=int, default=DEFAULT_SAMPLES, help="Bliss prior draws")
    ap.add_argument(
        "--tme-samples", type=int, default=DEFAULT_TME_SAMPLES, help="sim-tme prior draws (costly)"
    )
    ap.add_argument(
        "--pen-samples", type=int, default=DEFAULT_PEN_SAMPLES, help="penetration prior draws"
    )
    ap.add_argument("--workers", type=int, default=min(8, (os.cpu_count() or 4)))
    args = ap.parse_args()

    if args.headline == "bliss":
        _run_bliss(args)
    elif args.headline == "sim-tme":
        _run_tme(args)
    else:
        _run_penetration(args)


def _default_bliss(binary):
    """The unperturbed Bliss synergy (no FERRO_PARAM_OVERRIDES) — the manuscript
    point estimate the interval brackets. Reuses the shared CSV reader so the
    no-override path and run_bliss's override path never drift."""
    with tempfile.TemporaryDirectory(prefix="ferro_pp_default_") as workdir:
        env = dict(os.environ)
        env.pop("FERRO_PARAM_OVERRIDES", None)  # ensure a clean default run
        proc = subprocess.run([str(binary)], cwd=workdir, env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"default sim-combo-mech failed: {proc.stderr[-300:]}")
        return read_bliss_synergy(workdir)


if __name__ == "__main__":
    main()
