# In-Vivo Lipid Remodeling Simulation Results

## What this simulation tests

Dixon/Park (Cancer Research, 2025) showed that GPX4 inhibition kills cancer cells in 2D culture but fails in vivo because SCD1-driven MUFA enrichment displaces PUFAs from membranes, reducing ferroptosis susceptibility. This simulation tests whether that finding changes the conclusions about physical ROS modalities (modeled here as a shared exogenous ROS burst — SDT and PDT use identical parameters in this binary).

## Model

The MUFA protection model uses saturable logistic accumulation with natural lipid-turnover decay:

```
growth = rate × (1 - M / max)
decay  = decay_rate × M
M(t+1) = M(t) + growth - decay
effective_unsat = lipid_unsat × (1 - M), clamped ≥ 0.05
```

Steady state (with SCD1 active): `M_ss = rate × max / (rate + decay × max)`

Default parameters (from `Params::invivo()`):
- `scd_mufa_rate: 0.01` — accumulation rate
- `scd_mufa_max: 0.50` — maximum PUFA displacement
- `scd_mufa_decay: 0.005` — natural lipid turnover
- `initial_mufa_protection: 0.40` — steady state (established in-vivo remodeling)

Cells start with pre-accumulated protection representing established tumors, not freshly seeded 2D culture.

Biological basis: SCD1 is regulated by SREBP1/mTORC1 (not NRF2) and is constitutively active in 3D/in vivo. MUFA incorporation onset ~6-10h (Magtanong 2019), steady state ~48-72h. Protection factor ~3-5× (Tesfay 2019). Membrane phospholipid half-life ~24-48h drives the decay term.

## Three contexts

1. **2D**: `initial_mufa_protection=0, rate=0, decay=0` — standard in-vitro conditions, no MUFA remodeling
2. **In-vivo**: `initial_mufa_protection=0.40, rate=0.01, decay=0.005` — established tumor with active SCD1 maintaining MUFA at steady state
3. **In-vivo + SCD1 inhibitor**: `initial_mufa_protection=0.40, rate=0, decay=0.005` — SCD1 blocked, existing MUFA decays via natural lipid turnover

The SCD1i context is NOT identical to 2D. Cells start with pre-existing MUFA that gradually depletes, producing intermediate results.

## Key results

### Three-context comparison (100K cells per condition)

| Phenotype | Treatment | 2D Death% | In-Vivo Death% | SCD1i Death% | Protection Factor |
|-----------|-----------|-----------|----------------|-------------|-------------------|
| Glycolytic | Exo. ROS | 87.1% | 12.9% | 26.2% | 6.76× |
| OXPHOS | RSL3 | 1.2% | 0.0% | 0.01% | >>1× |
| OXPHOS | Exo. ROS | 99.9% | 90.5% | 96.0% | 1.10× |
| **Persister** | **RSL3** | **42.4%** | **2.3%** | **7.2%** | **18.6×** |
| **Persister** | **Exo. ROS** | **100.0%** | **99.98%** | **100.0%** | **1.00×** |
| Persister+NRF2 | RSL3 | 0.04% | 0.0% | 0.0% | >>1× |
| Persister+NRF2 | Exo. ROS | 99.5% | 90.5% | 94.4% | 1.10× |

Note: "Exo. ROS" = exogenous ROS modality. SDT and PDT are modeled identically in this binary (shared `sdt_ros`/`pdt_ros` = 5.0). Independent conclusions about SDT vs PDT cannot be drawn from this simulation.

### Biological predictions

1. **Dixon 2025 confirmed**: RSL3 kills 42.4% of persisters in 2D but only 2.3% in vivo. Pre-accumulated MUFA provides ~19× protection against pharmacologic GPX4 inhibition.

2. **Exogenous ROS still effective on persisters**: Physical ROS modalities maintain 99.98% kill on persisters even with established MUFA defense. The exogenous ROS burst overwhelms pre-accumulated MUFA protection. However, for other phenotypes the effect is substantial — glycolytic cells drop from 87% to 13%, and NRF2-compensated persisters drop from 99.5% to 90.5%.

3. **SCD1 inhibitor partially resensitizes**: SCD1i produces intermediate results (RSL3 on persisters: 2.3% → 7.2%), not full restoration to 2D levels (42.4%). This is because existing MUFA decays gradually — it is not instantly depleted. Full resensitization would require waiting for complete lipid turnover.

### Parameter sensitivity

MUFA sweep for Persister + SDT (50K cells per point):

| rate \ max | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|------------|------|------|------|------|------|
| 0.002 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.005 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.010 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.020 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 0.040 | 100.0% | 100.0% | 100.0% | 100.0% | 99.9% |

Exogenous ROS kills persisters at ≥99.9% across the entire tested parameter space. The result is insensitive to MUFA parameter choices.

MUFA sweep for Persister + RSL3 (50K cells per point):

| rate \ max | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|------------|------|------|------|------|------|
| 0.002 | 31.3% | 30.0% | 29.3% | 28.8% | 28.5% |
| 0.005 | 24.1% | 20.4% | 17.8% | 16.1% | 14.8% |
| 0.010 | 19.8% | 13.8% | 9.6% | 7.0% | 5.1% |
| 0.020 | 17.2% | 9.7% | 5.1% | 2.4% | 1.1% |
| 0.040 | 15.9% | 7.8% | 3.2% | 1.1% | 0.2% |

RSL3 shows a steep gradient: at low MUFA (rate=0.002, max=0.20), efficacy drops from 42% to 31%. At high MUFA (rate=0.04, max=0.60), efficacy drops to 0.2% — a >200× reduction. Note: the sweep starts cells at mufa_protection=0 (onset scenario), not at steady state. The main comparison uses steady-state initial conditions.

## What conclusions survive

1. **Exogenous ROS kills persisters even with established in-vivo lipid remodeling.** This is the simulation's strongest result. For persisters specifically, pre-accumulated MUFA barely dents exogenous ROS efficacy (99.98% vs 100%).

2. **Physical modalities have a specific advantage over drugs for this resistance mechanism.** RSL3 loses most of its efficacy (42% → 2.3%), but exogenous ROS maintains near-complete kill on persisters. The mechanism: drugs inhibit a single enzyme (GPX4), which MUFA-mediated substrate reduction can compensate. Exogenous ROS imposes overwhelming oxidative stress that depletes all defenses simultaneously.

3. **The advantage is strongest for persisters, weaker for other phenotypes.** Glycolytic cells (87% → 13%) and NRF2-compensated persisters (99.5% → 90.5%) are substantially affected by established MUFA. Only the FSP1-downregulated persister phenotype retains near-complete vulnerability.

## What conclusions are now 2D-only artifacts

1. **RSL3 efficacy on persisters (42.4%) is a 2D artifact.** In vivo, this drops to 2.3%.

2. **Exogenous ROS efficacy on glycolytic cells (87%) is largely a 2D artifact.** In vivo, this drops to 13%.

3. **The near-100% kill across all phenotypes under SDT/PDT is partially a 2D artifact.** In vivo, only the persister phenotype retains near-complete kill. Other phenotypes show 10-30% survival.

## Limitations

- SDT and PDT are modeled identically (shared exogenous ROS parameter). Independent claims about SDT vs PDT are not supported by this simulation.
- The decay model is first-order (constant fractional turnover). Real lipid dynamics are more complex.
- The MUFA parameter sweep starts cells at mufa_protection=0 rather than steady state, so it represents onset dynamics rather than established tumors.
- Other in-vivo resistance axes (DHODH, DHCR7/7-DHC, stromal buffering) are not modeled.
- The 0.05 floor on effective_unsat means cells always retain minimal PUFA vulnerability.

## References

- Dixon SJ, Park VS, et al. "Lipid Composition Alters Ferroptosis Sensitivity." Cancer Research 85(22):4380-4397, 2025.
- Magtanong L, et al. "Exogenous Monounsaturated Fatty Acids Promote a Ferroptosis-Resistant Cell State." Cell Chemical Biology, 2019.
- Tesfay L, et al. "Stearoyl-CoA Desaturase 1 Protects Ovarian Cancer Cells from Ferroptotic Cell Death." Cancer Research, 2019.
