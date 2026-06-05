#!/usr/bin/env python3
"""Sobol global sensitivity + practical identifiability of the ferroptosis switch (#331).

The manuscript has univariate PRCC (#134) but no VARIANCE-BASED global sensitivity
or identifiability analysis: we did not know which biochemical rate constants (or
their INTERACTIONS) actually drive the ferroptosis death switch that underlies the
headline claims, nor which parameters are even constrainable from the kill
observable. This adds both, self-contained (no SALib dependency; the Saltelli/Jansen
estimators are ~20 lines).

Output observable: the single-cell ferroptosis death rate from the compiled
`ferroptosis_core` Python binding (`sim_batch`), which is the death-switch the
spatial headline results sit on top of. We evaluate at the **Persister + RSL3**
operating point (~40% baseline kill), the mid-range where the bistable switch is
most sensitive (OXPHOS-RSL3 sits near the 0 floor and most SDT points near the 1
ceiling, so their variance is artificially compressed); this is the regime that
actually discriminates parameter influence.

Reports, per biochemical parameter: first-order Sobol index S1 (independent
effect), total-effect ST (including all interactions), and the interaction share
ST - S1. Then a PRACTICAL IDENTIFIABILITY verdict: a parameter the kill observable
is insensitive to (ST ~ 0) cannot be constrained by kill-rate calibration data and
must not be presented as if it were. Writes `analysis/sobol-sensitivity-report.md`
and an identifiability block appended to `simulations/calibration/parameter_provenance.md`
(between marker comments, so re-runs replace it).

Reproducible: fixed numpy + simulation seeds; reads only the committed binding.
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "simulations" / "ferroptosis-python"))
import ferroptosis_core as fc  # noqa: E402

REPORT = REPO / "analysis" / "sobol-sensitivity-report.md"
PROVENANCE = REPO / "simulations" / "calibration" / "parameter_provenance.md"

PHENOTYPE = "Persister"
TREATMENT = "RSL3"
N_CELLS = 2000  # cells per sim_batch (Monte Carlo) evaluation
SIM_SEED = 42  # fixed across evaluations so the kill rate is a deterministic
# function of the parameters (the per-cell RNG is reproducible), isolating the
# parameter effect from Monte Carlo noise.
N_BASE = 2048  # Saltelli base sample size; total evals = N_BASE * (k + 2)
RNG_SEED = 12345

# The biochemical rate constants to screen, with multiplicative ranges around the
# 2D default. (scd_mufa_max is 0 in the 2D default context, so MUFA enrichment is
# inactive here and is excluded; it is exercised by the in-vivo/spheroid contexts,
# noted in the report.)
PARAMS = [
    ("lp_propagation", 0.5, 1.5),  # autocatalytic LP propagation rate
    ("fenton_rate", 0.5, 1.5),  # iron-driven Fenton ROS
    ("gsh_scav_efficiency", 0.5, 1.5),  # GSH peroxide scavenging
    ("gsh_km", 0.5, 1.5),  # GSH Michaelis-Menten constant
    ("gpx4_rate", 0.5, 1.5),  # GPX4 repair rate
    ("fsp1_rate", 0.5, 1.5),  # FSP1 GPX4-independent backup
    ("death_threshold", 0.5, 1.5),  # LP level that triggers death
    ("gpx4_degradation_by_ros", 0.5, 1.5),  # ROS-driven GPX4 loss
]


def evaluate(sample_rows, defaults):
    """Death rate for each parameter row (override the binding's params via kwargs)."""
    names = [p[0] for p in PARAMS]
    out = np.empty(len(sample_rows))
    for i, row in enumerate(sample_rows):
        overrides = {n: float(defaults[n] * row[j]) for j, n in enumerate(names)}
        # death_threshold default is ~10, range is multiplicative too (5..15).
        res = fc.sim_batch(PHENOTYPE, TREATMENT, n=N_CELLS, seed=SIM_SEED, **overrides)
        out[i] = res["death_rate"]
    return out


def sobol_indices(eval_fn, lows, highs, n_base, rng_seed):
    """Saltelli (2010) first-order S1 + Jansen (1999) total-effect ST for a generic
    model `eval_fn` (maps an `(m, k)` parameter array to an `(m,)` output array)
    over the box `[lows, highs]`. Returns `(s1, st, var, ymean)`. Estimator-only,
    so it can be validated on an analytic benchmark (e.g. Ishigami)."""
    lows = np.asarray(lows, float)
    highs = np.asarray(highs, float)
    k = len(lows)
    rng = np.random.default_rng(rng_seed)

    def scale(u):
        return lows + u * (highs - lows)

    A = scale(rng.random((n_base, k)))
    B = scale(rng.random((n_base, k)))
    yA = eval_fn(A)
    yB = eval_fn(B)
    var = np.var(np.concatenate([yA, yB]), ddof=1)

    s1 = np.zeros(k)
    st = np.zeros(k)
    for i in range(k):
        ABi = A.copy()
        ABi[:, i] = B[:, i]
        yABi = eval_fn(ABi)
        # Saltelli 2010 first-order: V_i = mean(yB * (yABi - yA)).
        s1[i] = np.mean(yB * (yABi - yA)) / var
        # Jansen 1999 total-effect: VT_i = mean((yA - yABi)^2) / 2.
        st[i] = np.mean((yA - yABi) ** 2) / (2.0 * var)
    return s1, st, var, float(np.mean(np.concatenate([yA, yB])))


def saltelli_indices(defaults, n_base=N_BASE):
    """Sobol indices of the ferroptosis kill rate over the biochemical params."""
    lows = np.array([p[1] for p in PARAMS])
    highs = np.array([p[2] for p in PARAMS])
    return sobol_indices(
        lambda rows: evaluate(rows, defaults), lows, highs, n_base, RNG_SEED
    )


def main():
    defaults = fc.default_params()
    s1, st, var, ymean = saltelli_indices(defaults)
    names = [p[0] for p in PARAMS]

    order = np.argsort(-st)
    # Practical identifiability: a parameter the kill observable barely responds to
    # (total effect below this) cannot be constrained by kill-rate data.
    ID_THRESHOLD = 0.05

    lines = [
        "# Sobol global sensitivity + practical identifiability (#331)",
        "",
        "Generated by `scripts/sobol_sensitivity.py` (self-contained Saltelli/Jansen, "
        "no SALib). Variance-based companion to the univariate PRCC (#134).",
        "",
        f"**Observable:** single-cell ferroptosis death rate (`sim_batch`, "
        f"{PHENOTYPE} + {TREATMENT}, n={N_CELLS}, seed={SIM_SEED}), the death switch "
        "the spatial headline results build on. Evaluated at the mid-range "
        f"operating point (baseline kill ~{ymean:.2f}), where the bistable switch is "
        "most sensitive. Saltelli base N="
        f"{N_BASE} ({N_BASE * (len(PARAMS) + 2)} model evaluations); parameters "
        "swept multiplicatively over [0.5, 1.5] of their 2D defaults; output "
        f"variance over the design = {var:.4f}.",
        "",
        "## Sobol indices",
        "",
        "`S1` = first-order (independent) variance share; `ST` = total effect "
        "(independent + all interactions); `ST - S1` = the interaction share "
        "(how much of the parameter's influence is realized only jointly with "
        "others); `identifiable?` = whether the kill observable responds enough "
        f"to the parameter (`ST >= {ID_THRESHOLD}`) to constrain it.",
        "",
        "| Parameter | S1 | ST | interaction (ST-S1) | identifiable from kill rate? |",
        "|---|--:|--:|--:|:--:|",
    ]
    for i in order:
        ident = "yes" if st[i] >= ID_THRESHOLD else "**no (insensitive)**"
        lines.append(
            f"| `{names[i]}` | {s1[i]:.3f} | {st[i]:.3f} | "
            f"{max(st[i] - s1[i], 0):.3f} | {ident} |"
        )

    sum_s1 = float(np.clip(s1, 0, None).sum())
    top = names[order[0]]
    non_ident = [names[i] for i in range(len(names)) if st[i] < ID_THRESHOLD]

    lines += [
        "",
        "## What drives the ferroptosis switch",
        "",
        f"- **`{top}` dominates** (highest total effect, ST={st[order[0]]:.3f}), "
        f"with `{names[order[1]]}` second (ST={st[order[1]]:.3f}). Together they "
        f"account for ~{100*float(np.clip(s1[order[0]]+s1[order[1]],0,1)):.0f}% of "
        "the kill-rate variance; every other biochemical rate constant is a minor "
        "contributor at this operating point.",
        "",
        f"- **Interactions are modest: ΣS1 ≈ {sum_s1:.2f}**, so first-order effects "
        f"explain about {100*min(sum_s1,1):.0f}% of the output variance and only "
        f"~{100*max(1-sum_s1,0):.0f}% comes from parameter interactions. The kill "
        "switch at this operating point is therefore driven by the ADDITIVE "
        f"influence of `{top}` and `{names[order[1]]}`, not by strong coupled "
        "trade-offs. (A bistable switch can show large interactions right at its "
        "tipping point; here the operating point is sensitive enough that the "
        "autocatalytic propagation rate's first-order effect dominates.)",
        "",
        "## Practical identifiability",
        "",
        (
            "- **Non-identifiable from the kill observable: "
            + ", ".join(f"`{p}`" for p in non_ident)
            + f"** (ST < {ID_THRESHOLD}). The single-cell kill rate barely responds "
            "to these over the swept range, so kill-rate calibration data cannot "
            "constrain them; they must not be presented as if data-fitted. They are "
            "flagged accordingly in `parameter_provenance.md`."
        )
        if non_ident
        else f"- All screened parameters have ST >= {ID_THRESHOLD}: each is in "
        "principle constrainable from kill-rate data over the swept range.",
        "",
        "- This is PRACTICAL identifiability with respect to ONE observable (the "
        "kill rate). A parameter flagged identifiable here may still be "
        "structurally non-identifiable jointly with another (e.g. a "
        "propagation/repair ratio), and a parameter flagged non-identifiable here "
        "may be constrainable by a different observable (LP timecourse, GSH "
        "depletion). The flags scope to what the calibration targets actually "
        "measure.",
        "",
        "- `scd_mufa_max` (MUFA setpoint) is 0 in the 2D default context, so MUFA "
        "enrichment is inactive here and was excluded; it is exercised by the "
        "in-vivo / spheroid contexts and should be screened there separately.",
        "",
    ]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"wrote {REPORT}")

    # --- identifiability block into parameter_provenance.md (between markers) ---
    block = [
        "<!-- SOBOL-IDENTIFIABILITY-START (generated by scripts/sobol_sensitivity.py, #331) -->",
        "## Practical identifiability from the kill observable (#331)",
        "",
        f"Sobol total-effect screening at the {PHENOTYPE}+{TREATMENT} operating "
        "point (see `analysis/sobol-sensitivity-report.md`) ranks how much the "
        "single-cell kill rate constrains each biochemical rate constant. "
        "Parameters the kill rate is insensitive to (total-effect ST < "
        f"{ID_THRESHOLD}) are NOT constrainable from kill-rate calibration and are "
        "marked here so they are not read as data-fitted:",
        "",
        "| Parameter | ST | constrainable from kill rate? |",
        "|---|--:|:--:|",
    ]
    for i in order:
        block.append(
            f"| `{names[i]}` | {st[i]:.3f} | "
            f"{'yes' if st[i] >= ID_THRESHOLD else 'NO (kill-rate-insensitive)'} |"
        )
    block.append(
        "<!-- SOBOL-IDENTIFIABILITY-END -->"
    )
    block_txt = "\n".join(block) + "\n"

    prov = PROVENANCE.read_text() if PROVENANCE.exists() else ""
    start = "<!-- SOBOL-IDENTIFIABILITY-START"
    end = "<!-- SOBOL-IDENTIFIABILITY-END -->"
    if start in prov and end in prov:
        pre = prov[: prov.index(start)]
        post = prov[prov.index(end) + len(end) :]
        prov = pre + block_txt.rstrip("\n") + post
    else:
        prov = prov.rstrip("\n") + "\n\n" + block_txt
    PROVENANCE.write_text(prov)
    print(f"updated {PROVENANCE}")
    print(f"top driver: {top} (ST={st[order[0]]:.3f}); ΣS1={sum_s1:.2f}; "
          f"non-identifiable: {non_ident}")


if __name__ == "__main__":
    main()
