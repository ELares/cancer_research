# Resistant-State Map

First-pass scaffold for analyzing the corpus by resistant state rather than by modality alone.

These state assignments are keyword-derived heuristics. They are intended to support prioritization and literature review, not to assert that a paper experimentally validated a state transition.

Current resistant-state coverage in the index: 10/4830 records (0.2%).


## Resistant States Tracked

- **epigenetically-plastic**
- **nrf2-compensated-ferroptosis-resistant**
- **oxphos-dependent-persister**
- **slc7a11-high-disulfidptosis-prone**
- **stromal-sheltered-immune-excluded**
- **therapy-induced-senescence**

## State × Mechanism Counts

| Resistant State | Top linked mechanisms | Tagged articles |
|---|---|---|
| **epigenetically-plastic** | none | 0 |
| **nrf2-compensated-ferroptosis-resistant** | none | 0 |
| **oxphos-dependent-persister** | metabolic-targeting (2), epigenetic (1), bioelectric (1) | 2 |
| **slc7a11-high-disulfidptosis-prone** | nanoparticle (1), sonodynamic (1) | 1 |
| **stromal-sheltered-immune-excluded** | immunotherapy (7), bispecific-antibody (1), oncolytic-virus (1), metabolic-targeting (1), car-t (1) | 7 |
| **therapy-induced-senescence** | none | 0 |

## Interpretation

- The repo should use these states as the primary decision layer when comparing interventions.
- Physical ROS modalities should be framed as best-matched to OXPHOS-dependent, ferroptosis-prone persisters rather than as a universal answer to resistance.
- Senescence, stromal sheltering, and NRF2/SLC7A11 compensation should be treated as parallel escape states, not edge cases.