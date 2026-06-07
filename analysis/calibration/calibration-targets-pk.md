# Ferroptosis-inducer PK: the #334 calibration target

This is the measured-data target for the tumor-PK layer, which until now has been
honest-but-uncalibrated scaffolding (`simulations/calibration/CALIBRATION_STATUS.md`:
"RSL3 pharmacokinetics, Uncalibrated, order-of-magnitude estimates"). It is the
**data file** for #334; the model anchoring is `scripts/calibrate_pk.py` (which
writes `pk-calibration.md` + `.json`).

Committed artifact:
- `analysis/calibration/pk_measured_data.csv`: published summary PK (NCA metrics and
  popPK parameters) for two ferroptosis-relevant drugs, in long format
  (`drug, species, population, route, dose, compartment, parameter, value, unit,
  source_pmid, note`). These are transcribed published table values with PMIDs, not
  digitized curves, so the file is the reproducible target (no raw download).

## Why these two drugs

The canonical tool compounds the biochem layer uses (RSL3, ML162, ML210) and the
prototype system-xc- inducer erastin have **no usable in-vivo PK**: erastin's poor
metabolic stability and solubility are exactly why a more stable analog was
engineered. So the model's RSL3 concentrations are in-vitro / abstract, not
pharmacokinetically grounded, and that gap is real, not an oversight.

The two drugs with extractable, ferroptosis-relevant in-vivo PK are:

| drug | role | why it anchors the model |
|---|---|---|
| **IKE** (imidazole ketone erastin) | system-xc- inhibitor, engineered for in-vivo stability | the ONLY public ferroptosis-specific dataset with a PAIRED plasma + tumor concentration-time course; the mechanism (system-xc- to ferroptosis) is the model's biology |
| **sorafenib** | clinical kinase inhibitor that also induces ferroptosis via system-xc- | a published human population-PK model gives a complete, human-scale parameter set for a forward sanity check |

- **IKE**: Zhang et al. 2019, Cell Chem Biol, PMID 30799221 (mouse; IP/IV/PO route
  comparison plus a paired plasma+tumor distribution study in NCG SUDHL6 xenografts).
- **sorafenib**: Jain et al. 2011, Br J Clin Pharmacol, PMID 21392074 (human
  population PK: one-compartment + transit absorption + enterohepatic recycling).

## The numbers (IKE, the primary anchor)

The distribution study (NCG SUDHL6 xenograft, IP 50 mg/kg) reports paired plasma and
tumor NCA, the rows the fit uses:

| compartment | Tmax (h) | Cmax (ng/mL) | terminal t1/2 (h) | AUC (ng*h/mL) |
|---|---:|---:|---:|---:|
| plasma | 1.35 | 5185 | 1.83 | 10926 |
| tumor  | 3.30 | 2516 | 3.50 | 9857 |

What this tells the model:
1. **Tissue:plasma partition** Kp = AUC_tumor / AUC_plasma = **0.90**. The presets'
   `partition_coeff` is currently 0.15 to 0.5 (estimated); the measured value is
   nearly 1, i.e. tumor exposure of this system-xc- inducer almost matches plasma.
2. **Plasma to tumor delay** of ~2 h (Tmax 1.35 to 3.30 h) and **slower tumor
   clearance** (terminal t1/2 3.50 vs 1.83 h): the signature of a real distributional
   tissue compartment.
3. **A structural constraint**: the plasma AUC/Cmax ratio (2.11 h) is below the hard
   1-compartment floor `e * Tmax` (3.67 h), proving the disposition is
   multi-compartment (a fast distribution phase). See the fit report.

The IKE route-comparison rows (NOD/SCID cohort: IP Cmax 19515, IV 11384, PO 5203
ng/mL) are committed for context (they show IP > IV > PO exposure and an absolute
bioavailability reference) but are a DIFFERENT cohort, so the paired distribution
study is the one fit.

## sorafenib popPK (the human-scale anchor)

| parameter | value | meaning |
|---|---:|---|
| CL/F | 8.13 L/h | apparent clearance |
| V/F | 213 L | apparent volume (80 kg reference) |
| MTT | 1.98 h | mean absorption transit time |
| k_EHC | 0.857 1/h | enterohepatic recycling rate |
| F_ent | 0.498 | fraction entering enterohepatic recycling |

These give a derived terminal half-life of ~18 h and a predicted 400 mg twice-daily
steady-state average concentration `Cavg,ss = Dose/(CL/F * tau)` of ~4.1 mg/L,
within the clinically reported ~3 to 5 mg/L range, a forward consistency check on the
model's PK scale at human dose.

## Caveats the fit leg must handle (do NOT skip these)

This file is the target, not a calibration. Turning it into an anchored PK model
requires honesty about three limits, all carried into the fit report:

1. **Identifiability**: four summary NCA numbers per compartment do not uniquely
   determine a multi-compartment plasma model. The fit anchors the robust,
   identifiable quantities (partition, half-lives, peak level/timing, delay
   direction) and reports the rest (AUC, tumor Tmax magnitude) as flagged residuals,
   rather than over-fitting an unidentifiable model.
2. **Scope**: this anchors the PK structure + partition + delay for two drugs. It
   does NOT recalibrate the per-tumor-type presets (breast/pancreatic/GBM/melanoma/
   sarcoma): no public per-tumor-type measured PK exists for a ferroptosis inducer,
   so those stay documented estimates.
3. **Species / form**: IKE is a mouse anchor with summary NCA targets (the
   per-timepoint curve is in the cited supplement); sorafenib is human but a forward
   check, not a per-timepoint fit. The tool-compound (RSL3/ML) PK gap is flagged, not
   fabricated.

Until the fit lands, the tumor-PK layer stays labelled uncalibrated in
`CALIBRATION_STATUS.md`.
