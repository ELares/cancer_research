#!/usr/bin/env python3
"""PRCC global sensitivity analysis for ferroptosis simulation parameters.

Uses Latin Hypercube Sampling to explore the 11-dimensional parameter space
and computes Partial Rank Correlation Coefficients (PRCC) against death_rate
for key phenotype × treatment conditions.

Requires: ferroptosis_core Python bindings
  cd simulations && maturin develop -m ferroptosis-python/Cargo.toml --release

Usage:
  python scripts/run_prcc.py                  # 2000 samples (default)
  python scripts/run_prcc.py --samples 500    # quick test run

Reference: Marino et al., J Theor Biol 254(3):178-196, 2008.
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

try:
    import ferroptosis_core as fc
except ImportError:
    print("ERROR: ferroptosis_core Python bindings not installed.")
    print("Run: cd simulations && maturin develop -m ferroptosis-python/Cargo.toml --release")
    raise SystemExit(1)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_FILE = PROJECT_ROOT / "analysis" / "prcc-results.json"

# 11 parameters matching the existing local sensitivity analysis (sim-original).
# Ranges: ±50% of default, capped at biologically meaningful bounds.
PARAM_RANGES = {
    "fenton_rate":              (0.01, 0.04),
    "gsh_scav_efficiency":      (0.25, 1.0),
    "lp_rate":                  (0.03, 0.12),
    "lp_propagation":           (0.05, 0.20),
    "gpx4_rate":                (0.15, 0.60),
    "fsp1_rate":                (0.04, 0.16),
    "nrf2_gsh_rate":            (0.0125, 0.05),
    "gpx4_degradation_by_ros":  (0.001, 0.004),
    "death_threshold":          (5.0, 20.0),
    "sdt_ros":                  (2.5, 10.0),
    "rsl3_gpx4_inhib":          (0.80, 0.99),
}

CONDITIONS = [
    ("Persister", "SDT"),
    ("Persister", "RSL3"),
    ("Glycolytic", "SDT"),
    ("Glycolytic", "RSL3"),
    ("Persister", "Control"),
    ("PersisterNrf2", "SDT"),
]


def latin_hypercube_sample(n, ranges, seed=42):
    """Generate n Latin Hypercube Samples across len(ranges) parameters."""
    rng = np.random.RandomState(seed)
    p = len(ranges)
    samples = np.zeros((n, p))
    for j, (lo, hi) in enumerate(ranges):
        # Stratified: one point per stratum, randomly placed within stratum
        strata = (np.arange(n) + rng.uniform(size=n)) / n
        rng.shuffle(strata)
        samples[:, j] = lo + strata * (hi - lo)
    return samples


def compute_prcc(inputs, output):
    """Compute PRCC for each input column against the output vector.

    Uses the rank correlation matrix inversion method
    (Blower & Dowlatabadi, 1994; Marino et al., 2008).

    Returns list of (prcc, p_value) tuples, one per input column.
    """
    n, p = inputs.shape
    data = np.column_stack([inputs, output])
    ranked = np.column_stack([rankdata(data[:, j]) for j in range(p + 1)])
    R = np.corrcoef(ranked.T)
    R_inv = np.linalg.inv(R)

    results = []
    for i in range(p):
        prcc = -R_inv[i, p] / math.sqrt(R_inv[i, i] * R_inv[p, p])
        # t-test for significance (df = n - 2 - p)
        df = n - 2 - p
        if abs(prcc) < 1.0 and df > 0:
            t_stat = abs(prcc) * math.sqrt(df / (1.0 - prcc * prcc))
            # Normal approximation for df > 1000
            p_value = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(t_stat / math.sqrt(2.0))))
        else:
            p_value = 0.0
        results.append((float(prcc), float(p_value)))
    return results


def main():
    parser = argparse.ArgumentParser(description="PRCC global sensitivity analysis.")
    parser.add_argument("--samples", type=int, default=2000, help="Number of LHS samples (default: 2000)")
    parser.add_argument("--cells", type=int, default=1000, help="Cells per condition per sample (default: 1000)")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed (default: 42)")
    args = parser.parse_args()

    param_names = list(PARAM_RANGES.keys())
    ranges = list(PARAM_RANGES.values())

    print(f"PRCC Analysis: {args.samples} LHS samples × {len(CONDITIONS)} conditions × {args.cells} cells/condition")
    print(f"Parameters: {len(param_names)}, Seed: {args.seed}\n")

    X = latin_hypercube_sample(args.samples, ranges, seed=args.seed)

    all_results = {
        "metadata": {
            "n_samples": args.samples,
            "n_cells_per_condition": args.cells,
            "seed": args.seed,
            "n_parameters": len(param_names),
            "parameters": param_names,
            "parameter_ranges": {name: list(r) for name, r in PARAM_RANGES.items()},
        },
        "conditions": [],
    }

    from scipy.stats import pearsonr, spearmanr

    # Store raw death_rates per condition for Sobol decision check after PRCC
    raw_death_rates = {}

    for phenotype, treatment in CONDITIONS:
        print(f"  {phenotype} × {treatment}...", end="", flush=True)
        death_rates = np.zeros(args.samples)
        for i in range(args.samples):
            kwargs = {param_names[j]: float(X[i, j]) for j in range(len(param_names))}
            result = fc.sim_batch(phenotype, treatment, n=args.cells, seed=i + 100, **kwargs)
            death_rates[i] = result["death_rate"]

        raw_death_rates[(phenotype, treatment)] = death_rates
        prcc_results = compute_prcc(X, death_rates)

        condition_data = {
            "phenotype": phenotype,
            "treatment": treatment,
            "output": "death_rate",
            "mean_death_rate": float(np.mean(death_rates)),
            "std_death_rate": float(np.std(death_rates)),
            "prcc": {param_names[j]: prcc_results[j][0] for j in range(len(param_names))},
            "p_values": {param_names[j]: prcc_results[j][1] for j in range(len(param_names))},
            "ranked_by_abs_prcc": sorted(
                param_names, key=lambda p: abs(prcc_results[param_names.index(p)][0]), reverse=True
            ),
        }
        all_results["conditions"].append(condition_data)

        top3 = condition_data["ranked_by_abs_prcc"][:3]
        top3_str = ", ".join(f"{p}={condition_data['prcc'][p]:+.3f}" for p in top3)
        print(f" done. Top 3: {top3_str}")

    # Sobol decision criterion: compare Pearson (linear) vs Spearman (monotonic)
    # correlations for each parameter × condition. Large divergence indicates
    # non-monotonic relationships where PRCC (which assumes monotonicity) may
    # be unreliable, warranting Sobol analysis.
    print("\n  Checking Sobol escalation criterion (Pearson vs Spearman divergence)...")
    DIVERGENCE_THRESHOLD = 0.2
    non_monotonic_flags = []

    for (phenotype, treatment), death_rates in raw_death_rates.items():
        if np.std(death_rates) < 1e-10:
            continue  # constant output — skip
        for j, pname in enumerate(param_names):
            r_pearson, _ = pearsonr(X[:, j], death_rates)
            r_spearman, _ = spearmanr(X[:, j], death_rates)
            divergence = abs(r_pearson - r_spearman)
            if divergence > DIVERGENCE_THRESHOLD:
                non_monotonic_flags.append({
                    "condition": f"{phenotype} × {treatment}",
                    "parameter": pname,
                    "pearson": round(r_pearson, 3),
                    "spearman": round(r_spearman, 3),
                    "divergence": round(divergence, 3),
                })

    sobol_warranted = len(non_monotonic_flags) > 2
    all_results["sobol_decision"] = {
        "warranted": sobol_warranted,
        "criterion": f"|Pearson − Spearman| > {DIVERGENCE_THRESHOLD} on >2 parameter×condition pairs",
        "divergence_threshold": DIVERGENCE_THRESHOLD,
        "non_monotonic_flags": non_monotonic_flags,
        "n_flags": len(non_monotonic_flags),
        "reason": (
            f"Found {len(non_monotonic_flags)} parameter×condition pairs with "
            f"|Pearson − Spearman| > {DIVERGENCE_THRESHOLD}. "
            + ("Stage 2 (Sobol) warranted." if sobol_warranted
               else "Below the >2 threshold for Sobol escalation.")
        ),
    }
    print(f"  Flags: {len(non_monotonic_flags)}, Sobol warranted: {sobol_warranted}")

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(json.dumps(all_results, indent=2) + "\n")
    print(f"\nResults written to {OUTPUT_FILE}")
    print(f"Sobol warranted: {sobol_warranted}")


if __name__ == "__main__":
    main()
