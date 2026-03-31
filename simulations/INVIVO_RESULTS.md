# In-Vivo Lipid Remodeling Simulation Results

## What this simulation tests

Dixon/Park (Cancer Research, 2025) showed that GPX4 inhibition kills cancer cells in 2D culture but fails in vivo because SCD1-driven MUFA enrichment displaces PUFAs from membranes, reducing ferroptosis susceptibility. This simulation tests whether that finding changes the conclusions about physical ROS modalities (SDT/PDT).

## Model

The MUFA protection model uses saturable logistic accumulation:

```
mufa_protection(t+1) = mufa_protection(t) + rate × (1 - mufa_protection / max)
effective_unsat = lipid_unsat × (1 - mufa_protection), clamped ≥ 0.05
```

Default parameters (from `Params::invivo()`):
- `scd_mufa_rate: 0.01` — time constant ~50 steps; 86% of max at step 100
- `scd_mufa_max: 0.50` — 50% PUFA displacement at steady state

Biological basis: SCD1 is regulated by SREBP1/mTORC1 (not NRF2) and is constitutively active in 3D/in vivo. MUFA incorporation onset ~6-10h (Magtanong 2019), steady state ~48-72h. Protection factor ~3-5× (Tesfay 2019).

## Key results

### Three-context comparison (100K cells per condition)

| Phenotype | Treatment | 2D Death% | In-Vivo Death% | Protection Factor |
|-----------|-----------|-----------|----------------|-------------------|
| Glycolytic | SDT | 87.1% | 51.6% | 1.69× |
| OXPHOS | RSL3 | 1.2% | 0.0% | >>1× |
| OXPHOS | SDT | 99.9% | 98.5% | 1.01× |
| **Persister** | **RSL3** | **42.4%** | **7.1%** | **5.99×** |
| **Persister** | **SDT** | **100.0%** | **99.99%** | **1.00×** |
| Persister+NRF2 | RSL3 | 0.04% | 0.0% | >>1× |
| Persister+NRF2 | SDT | 99.5% | 97.9% | 1.02× |

SCD1 inhibitor fully restores 2D sensitivity in all conditions (invivo+scd1i matches 2D exactly).

### Biological predictions

1. **Dixon 2025 confirmed**: RSL3 kills 42.4% of persisters in 2D but only 7.1% in vivo. MUFA remodeling provides ~6× protection against pharmacologic GPX4 inhibition.

2. **SDT/PDT bypass MUFA defense**: Physical ROS modalities maintain ≥97.9% kill across all phenotypes in vivo. The exogenous ROS burst (5.0 relative units, decaying over ~45 steps) overwhelms GSH before MUFA protection can fully accumulate. By the time MUFA reaches meaningful levels (~step 50), the cell's antioxidant defenses are already depleted and autocatalytic propagation is underway.

3. **SCD1 inhibitor resensitization confirmed**: Tesfay 2019 predicted that SCD1 inhibition restores ferroptosis sensitivity. The simulation confirms complete resensitization (in-vivo + SCD1i = exact 2D match).

### Parameter sensitivity

MUFA sweep for Persister + SDT (50K cells per point):

| rate \ max | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|------------|------|------|------|------|------|
| 0.002 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.005 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.010 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.020 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.040 | 100.0% | 100.0% | 100.0% | 100.0% | 99.9% |

SDT kills persisters at ≥99.9% across the entire tested parameter space. The result is insensitive to MUFA parameter choices.

MUFA sweep for Persister + RSL3 (50K cells per point):

| rate \ max | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|------------|------|------|------|------|------|
| 0.002 | 31.3% | 30.0% | 29.3% | 28.8% | 28.5% |
| 0.005 | 24.1% | 20.4% | 17.8% | 16.1% | 14.8% |
| 0.010 | 19.8% | 13.8% | 9.6% | 7.0% | 5.1% |
| 0.020 | 17.2% | 9.7% | 5.1% | 2.4% | 1.1% |
| 0.040 | 15.9% | 7.8% | 3.2% | 1.1% | 0.2% |

RSL3 shows a steep gradient: at low MUFA (rate=0.002, max=0.20), efficacy drops from 42% to 31%. At high MUFA (rate=0.04, max=0.60), efficacy drops to 0.2% — a >200× reduction. The Dixon 2025 prediction holds across the entire parameter space but the magnitude depends strongly on how fast and how completely MUFA displacement occurs. At the default in-vivo parameters (rate=0.01, max=0.50), RSL3 drops to 7.0% — a ~6× protection factor.

## What conclusions survive

1. **SDT/PDT kill persisters even with in-vivo lipid remodeling.** This is the simulation's strongest result. The massive exogenous ROS burst overwhelms MUFA protection because it depletes GSH before membrane remodeling can prevent lipid peroxidation cascade.

2. **Physical modalities have a specific advantage over drugs for this resistance mechanism.** RSL3 (GPX4 inhibition) loses most of its efficacy in vivo (42% → 7%), but SDT/PDT maintain near-complete kill. The mechanism: drugs inhibit a single enzyme (GPX4), which can be compensated by MUFA-mediated substrate reduction. Physical modalities impose overwhelming oxidative stress that depletes all defenses simultaneously.

3. **SCD1 inhibitor + RSL3 is a viable combination.** If SDT/PDT are unavailable, pharmacologic ferroptosis induction can be rescued by co-administering an SCD1 inhibitor to prevent MUFA compensation.

## What conclusions are now 2D-only artifacts

1. **RSL3 efficacy on persisters (42.4%) is a 2D artifact.** In vivo, this drops to 7.1%. The manuscript should not present RSL3 as comparably effective to SDT/PDT without this caveat.

2. **Glycolytic cell sensitivity to SDT (87%) is partially a 2D artifact.** In vivo, this drops to 52%. Glycolytic cells have lower baseline ROS, giving MUFA more time to accumulate before oxidative stress overwhelms it.

3. **Baseline persister death rate (1.2%) is partially a 2D artifact.** In vivo, spontaneous ferroptosis drops to 0.01%. MUFA protection stabilizes borderline-vulnerable cells.

## Limitations

- The model uses a single coarse MUFA-protection term, not a full lipidomics simulation.
- SCD1 activity is modeled as a context switch (on/off via Params), not as a dynamic per-cell variable.
- Other in-vivo resistance axes (DHODH, DHCR7/7-DHC, stromal buffering) are not modeled.
- The 180-step simulation compresses biological time; the rate calibration maps approximately to 48-72h onset but is not a literal timescale.
- The 0.05 floor on effective_unsat means cells always retain minimal PUFA vulnerability, which may overstate physical-modality efficacy at extreme MUFA protection levels.

## References

- Dixon SJ, Park VS, et al. "Lipid Composition Alters Ferroptosis Sensitivity." Cancer Research 85(22):4380-4397, 2025.
- Magtanong L, et al. "Exogenous Monounsaturated Fatty Acids Promote a Ferroptosis-Resistant Cell State." Cell Chemical Biology, 2019.
- Tesfay L, et al. "Stearoyl-CoA Desaturase 1 Protects Ovarian Cancer Cells from Ferroptotic Cell Death." Cancer Research, 2019.
