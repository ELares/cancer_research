#!/usr/bin/env python3
"""
Ferroptosis simulation parameter sweep and visualization.

Demonstrates the ferroptosis-core Python API with:
1. Basic simulation — single cell and batch
2. RSL3 dose-response curve — GPX4 inhibition sweep
3. Phenotype comparison — heatmap of death rates
4. 2D vs in-vivo context — MUFA protection effect
5. Combination exploration — dual-pathway depletion

Requirements: ferroptosis_core (build with maturin), matplotlib, numpy

Run:
    python parameter_sweep.py              # generates all plots
    python parameter_sweep.py --no-plots   # print tables only

Build the extension first:
    cd simulations
    python -m venv .venv && source .venv/bin/activate   # virtualenv required
    pip install maturin matplotlib numpy
    maturin develop -m ferroptosis-python/Cargo.toml --release
"""

import argparse
import sys
import time

try:
    import ferroptosis_core as fc
except ImportError:
    print("ERROR: ferroptosis_core not installed.")
    print("Build it first: maturin develop -m ferroptosis-python/Cargo.toml --release")
    sys.exit(1)

# ============================================================
# 1. Basics
# ============================================================

def demo_basics():
    print("=" * 60)
    print("1. BASICS")
    print("=" * 60)

    # Default parameters
    params = fc.default_params()
    print(f"\nDefault params ({len(params)} parameters):")
    for k in sorted(params)[:5]:
        print(f"  {k} = {params[k]}")
    print(f"  ... and {len(params) - 5} more\n")

    # Single cell
    result = fc.sim_cell("Persister", "RSL3", seed=42)
    print(f"Single cell (Persister + RSL3):")
    print(f"  Dead: {result['dead']}, LP: {result['lp']:.2f}, GSH: {result['gsh']:.2f}, GPX4: {result['gpx4']:.2f}\n")

    # Batch simulation with timing
    t0 = time.time()
    stats = fc.sim_batch("Persister", "RSL3", n=10000, seed=42)
    elapsed = time.time() - t0
    print(f"Batch (10,000 Persister + RSL3, {elapsed:.3f}s):")
    print(f"  Death rate: {stats['death_rate']:.1%} [{stats['ci_low']:.1%}, {stats['ci_high']:.1%}]")
    print(f"  Mean LP: {stats['mean_lp']:.2f}, GSH: {stats['mean_gsh']:.2f}, GPX4: {stats['mean_gpx4']:.2f}\n")


# ============================================================
# 2. RSL3 Dose-Response
# ============================================================

def demo_dose_response(plot=True):
    print("=" * 60)
    print("2. RSL3 DOSE-RESPONSE (GPX4 inhibition sweep)")
    print("=" * 60)

    inhibitions = [i / 20 for i in range(21)]  # 0.0 to 1.0 in 0.05 steps
    death_rates = []
    ci_lows = []
    ci_highs = []

    for inhib in inhibitions:
        stats = fc.sim_batch("Persister", "RSL3", n=1000, seed=42,
                             rsl3_gpx4_inhib=inhib)
        death_rates.append(stats["death_rate"])
        ci_lows.append(stats["ci_low"])
        ci_highs.append(stats["ci_high"])

    print(f"\n{'Inhibition':>12} {'Death Rate':>12} {'95% CI':>20}")
    print("-" * 46)
    for i in range(0, len(inhibitions), 4):
        inhib = inhibitions[i]
        dr = death_rates[i]
        lo, hi = ci_lows[i], ci_highs[i]
        print(f"  {inhib:>10.2f} {dr:>11.1%} [{lo:.1%}, {hi:.1%}]")

    if plot:
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(inhibitions, [d * 100 for d in death_rates], 'b-o', markersize=4)
            ax.fill_between(inhibitions,
                            [lo * 100 for lo in ci_lows],
                            [hi * 100 for hi in ci_highs],
                            alpha=0.2, color='blue')
            ax.set_xlabel("RSL3 GPX4 Inhibition Strength")
            ax.set_ylabel("Persister Death Rate (%)")
            ax.set_title("RSL3 Dose-Response on Persister Cells\n(n=1000 per point)")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 100)
            ax.grid(True, alpha=0.3)
            fig.savefig("dose_response.png", dpi=150, bbox_inches="tight")
            print(f"\n  Plot saved: dose_response.png")
            plt.close()
        except ImportError:
            print("\n  (matplotlib not available — skipping plot)")
    print()


# ============================================================
# 3. Phenotype Comparison
# ============================================================

def demo_phenotype_comparison(plot=True):
    print("=" * 60)
    print("3. PHENOTYPE × TREATMENT COMPARISON")
    print("=" * 60)

    phenotypes = ["Glycolytic", "OXPHOS", "Persister", "PersisterNrf2", "Stromal"]
    treatments = ["Control", "RSL3", "SDT", "PDT"]

    matrix = []
    print(f"\n{'':>15}", end="")
    for tx in treatments:
        print(f"{tx:>10}", end="")
    print()
    print("-" * 55)

    for pheno in phenotypes:
        row = []
        print(f"{pheno:>15}", end="")
        for tx in treatments:
            stats = fc.sim_batch(pheno, tx, n=1000, seed=42)
            row.append(stats["death_rate"])
            print(f"{stats['death_rate']:>9.1%}", end="")
        matrix.append(row)
        print()

    if plot:
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            fig, ax = plt.subplots(figsize=(8, 5))
            data = np.array(matrix) * 100
            im = ax.imshow(data, cmap="YlOrRd", aspect="auto", vmin=0, vmax=100)
            ax.set_xticks(range(len(treatments)))
            ax.set_xticklabels(treatments)
            ax.set_yticks(range(len(phenotypes)))
            ax.set_yticklabels(phenotypes)
            for i in range(len(phenotypes)):
                for j in range(len(treatments)):
                    color = "white" if data[i, j] > 50 else "black"
                    ax.text(j, i, f"{data[i, j]:.0f}%", ha="center", va="center",
                            color=color, fontsize=10)
            ax.set_title("Ferroptosis Death Rate by Phenotype × Treatment\n(n=1000 per condition)")
            plt.colorbar(im, ax=ax, label="Death Rate (%)")
            fig.savefig("phenotype_heatmap.png", dpi=150, bbox_inches="tight")
            print(f"\n  Plot saved: phenotype_heatmap.png")
            plt.close()
        except ImportError:
            print("\n  (matplotlib not available — skipping plot)")
    print()


# ============================================================
# 4. 2D vs In-Vivo (MUFA Protection)
# ============================================================

def demo_invivo(plot=True):
    print("=" * 60)
    print("4. 2D vs IN-VIVO CONTEXT (MUFA protection)")
    print("=" * 60)

    mufa_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    rsl3_rates = []
    sdt_rates = []

    for mufa in mufa_values:
        rsl3 = fc.sim_batch("Persister", "RSL3", n=1000, seed=42,
                            scd_mufa_rate=0.01, scd_mufa_max=0.5,
                            initial_mufa_protection=mufa, scd_mufa_decay=0.005)
        sdt = fc.sim_batch("Persister", "SDT", n=1000, seed=42,
                           scd_mufa_rate=0.01, scd_mufa_max=0.5,
                           initial_mufa_protection=mufa, scd_mufa_decay=0.005)
        rsl3_rates.append(rsl3["death_rate"])
        sdt_rates.append(sdt["death_rate"])

    print(f"\n{'MUFA':>8} {'RSL3':>10} {'SDT':>10} {'Protection':>12}")
    print("-" * 42)
    for i, mufa in enumerate(mufa_values):
        prot = rsl3_rates[0] / rsl3_rates[i] if rsl3_rates[i] > 0 else float('inf')
        print(f"  {mufa:>6.1f} {rsl3_rates[i]:>9.1%} {sdt_rates[i]:>9.1%} {prot:>11.1f}×")

    if plot:
        try:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(mufa_values, [r * 100 for r in rsl3_rates], 'r-o', label="RSL3", markersize=5)
            ax.plot(mufa_values, [r * 100 for r in sdt_rates], 'b-s', label="SDT", markersize=5)
            ax.set_xlabel("Initial MUFA Protection Level")
            ax.set_ylabel("Persister Death Rate (%)")
            ax.set_title("MUFA Protection: RSL3 Collapses, SDT Resists\n(n=1000 per point)")
            ax.legend()
            ax.set_xlim(0, 0.6)
            ax.set_ylim(0, 105)
            ax.grid(True, alpha=0.3)
            fig.savefig("mufa_protection.png", dpi=150, bbox_inches="tight")
            print(f"\n  Plot saved: mufa_protection.png")
            plt.close()
        except ImportError:
            print("\n  (matplotlib not available — skipping plot)")
    print()


# ============================================================
# 5. Combination Exploration
# ============================================================

def demo_combinations():
    print("=" * 60)
    print("5. COMBINATION EXPLORATION (dual-pathway depletion)")
    print("=" * 60)

    # Single agents
    ctrl = fc.sim_batch("Persister", "Control", n=1000, seed=42)
    rsl3 = fc.sim_batch("Persister", "RSL3", n=1000, seed=42)

    # FSP1 inhibition via reduced fsp1_rate (simulates FSP1i effect)
    fsp1i = fc.sim_batch("Persister", "RSL3", n=1000, seed=42,
                         rsl3_gpx4_inhib=0.0, fsp1_rate=0.01)  # FSP1i only

    # Combination: RSL3 + reduced FSP1
    combo = fc.sim_batch("Persister", "RSL3", n=1000, seed=42,
                         fsp1_rate=0.01)  # RSL3 + FSP1i

    # Bliss prediction
    bliss = rsl3["death_rate"] + fsp1i["death_rate"] - rsl3["death_rate"] * fsp1i["death_rate"]
    synergy = combo["death_rate"] / bliss if bliss > 0.001 else float('nan')

    print(f"\n  Control:      {ctrl['death_rate']:.1%}")
    print(f"  RSL3 alone:   {rsl3['death_rate']:.1%}")
    print(f"  FSP1i alone:  {fsp1i['death_rate']:.1%}")
    print(f"  RSL3 + FSP1i: {combo['death_rate']:.1%}")
    print(f"  Bliss prediction: {bliss:.1%}")
    print(f"  Synergy score: {synergy:.2f}× {'(SYNERGISTIC)' if synergy > 1.1 else '(~additive)'}")

    # Pathway traces
    print(f"\n  Pathway traces (mean final values):")
    print(f"  {'Condition':>15} {'GPX4':>8} {'GSH':>8} {'LP':>8}")
    print(f"  {'-'*41}")
    for name, s in [("Control", ctrl), ("RSL3", rsl3), ("FSP1i", fsp1i), ("RSL3+FSP1i", combo)]:
        print(f"  {name:>15} {s['mean_gpx4']:>8.3f} {s['mean_gsh']:>8.3f} {s['mean_lp']:>8.2f}")
    print()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Ferroptosis simulation examples")
    parser.add_argument("--no-plots", action="store_true", help="Skip plot generation")
    args = parser.parse_args()

    plot = not args.no_plots

    print(f"ferroptosis-core Python bindings")
    print(f"Parameters: {len(fc.default_params())} biochemistry rate constants")
    print(f"Phenotypes: Glycolytic, OXPHOS, Persister, PersisterNrf2, Stromal")
    print(f"Treatments: Control, RSL3, SDT, PDT")
    print(f"Contexts: 2d (default), invivo (SCD1/MUFA protection)\n")

    demo_basics()
    demo_dose_response(plot=plot)
    demo_phenotype_comparison(plot=plot)
    demo_invivo(plot=plot)
    demo_combinations()

    print("=" * 60)
    print("All examples complete.")
    if plot:
        print("Plots saved: dose_response.png, phenotype_heatmap.png, mufa_protection.png")
    print("=" * 60)


if __name__ == "__main__":
    main()
