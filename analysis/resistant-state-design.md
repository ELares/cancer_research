# Resistant-State Design

First-pass scaffold for analyzing the repo by resistant state rather than by modality alone.

## Why this layer exists

The current repo is strongest when it explains why a residual state becomes vulnerable, not when it argues that a single modality is universally important. The discussion around SDT/PDT/CAP is therefore best organized as:

- resistant state first
- matched intervention second

## Initial resistant states to track

### OXPHOS-dependent persister

- Hallmarks: oxidative phosphorylation, mitochondrial ROS, lipid-peroxidation-prone membranes, persister survival after frontline therapy
- Candidate matched interventions: SDT, PDT, CAP, GPX4 inhibition, FSP1 inhibition, HDAC pre-sensitization

### NRF2-compensated ferroptosis-resistant state

- Hallmarks: NRF2 activation, high GSH, GPX4/FSP1 buffering, antioxidant compensation
- Candidate matched interventions: epigenetic reprogramming, backup-defense suppression, combination ferroptosis targeting

### SLC7A11-high / disulfidptosis-relevant state

- Hallmarks: xCT/SLC7A11 dependence, cystine import, ferroptosis escape with a possible disulfidptosis liability
- Candidate matched interventions: glucose restriction logic, metabolic stress, disulfidptosis-oriented follow-up

### Therapy-induced senescence

- Hallmarks: growth arrest, SASP, BCL-2/BCL-XL dependence, immune-modulatory secretome
- Candidate matched interventions: senolytics, SASP suppression, immune remodeling

### Stromal-sheltered / immune-excluded state

- Hallmarks: CAFs, extracellular matrix, immune exclusion, physical shielding, delivery failure
- Candidate matched interventions: stromal remodeling, ECM-targeting, locoregional therapy combinations

### Epigenetically plastic state

- Hallmarks: chromatin-state switching, KDM5/EZH2/HDAC involvement, reversible drug tolerance
- Candidate matched interventions: HDAC/KDM5/EZH2 targeting, state-collapse combinations

## Implications for this repo

- PDT/SDT should remain in scope, but as interventions matched to OXPHOS-skewed, ferroptosis-prone residual disease.
- CAP belongs in the same modality-class discussion and should no longer be absent from the framing.
- Radioligands, cell therapies, vaccines, and stromal biology are not alternative curiosities; they are comparison lanes against the same residual-disease problem.

## Implementation note

The tagging pipeline now has a `resistant_states` axis in code, implemented as composite rules rather than single-keyword OR matches, but the full corpus has not yet been re-tagged in this PR. Regenerating that layer should happen in a follow-up once the taxonomy changes are reviewed.
