# Diagnostic-to-Therapy Matching Audit

First-pass extraction of diagnostic → targetable feature → intervention chains. Matching requires the intervention link plus at least one of (diagnostic, feature).

- Articles with at least one diagnostic-therapy link: **129** / 4830
- Chains evaluated: 6

## Chain Counts

| Chain | Articles | Top cancer types | Top evidence levels |
|---|---|---|---|
| **psma-imaging-to-radioligand** | 4 | prostate (3), kidney (1), breast (1) | preclinical-invivo (1) |
| **sstr-imaging-to-prrt** | 0 | . | . |
| **pdl1-ihc-to-checkpoint** | 13 | lung (6), bladder (2), melanoma (2) | phase2-clinical (4), phase1-clinical (2), clinical-other (1) |
| **tmb-msi-to-immunotherapy** | 33 | colorectal (15), lung (4), melanoma (3) | preclinical-invivo (5), clinical-other (4), phase2-clinical (2) |
| **neoantigen-profiling-to-mrna-vaccine** | 79 | melanoma (18), pancreatic (16), lung (8) | preclinical-invivo (12), phase1-clinical (6), preclinical-invitro (5) |
| **oncolytic-susceptibility-to-virotherapy** | 1 | melanoma (1) | preclinical-invivo (1) |

## Example Papers


### psma-imaging-to-radioligand

*Diagnostic:* psma pet, psma imaging, psma scan | *Feature:* psma expression, psma-positive, psma positive | *Intervention:* 177lu-psma, lu-psma, psma radioligand

- **PMID 35750683** (2022, 242 cites) — unclassified — immunotherapy, synthetic-lethality — *Targeting signaling pathways in prostate cancer: mechanisms and clinical trials.*
- **PMID 31640747** (2019, 44 cites) — preclinical-invivo — untagged — *Targeting of prostate-specific membrane antigen for radio-ligand therapy of triple-negative breast cancer.*
- **PMID 38302933** (2024, 16 cites) — unclassified — immunotherapy, radioligand-therapy — *A multicentric, single arm, open-label, phase I/II study evaluating PSMA targeted radionuclide therapy in adult patients*

### pdl1-ihc-to-checkpoint

*Diagnostic:* pd-l1 immunohistochemistry, pd-l1 ihc, pd-l1 staining | *Feature:* pd-l1 positive, pd-l1 high, pd-l1 expression | *Intervention:* pembrolizumab, nivolumab, atezolizumab

- **PMID 29347993** (2018, 328 cites) — phase2-clinical — untagged — *Updated efficacy of avelumab in patients with previously treated metastatic Merkel cell carcinoma after ≥1 year of follo*
- **PMID 33203644** (2021, 118 cites) — phase2-clinical — epigenetic, immunotherapy — *Entinostat plus Pembrolizumab in Patients with Metastatic NSCLC Previously Treated with Anti-PD-(L)1 Therapy.*
- **PMID 39722028** (2024, 49 cites) — unclassified — antibody-drug-conjugate, car-t, immunotherapy — *Current and future immunotherapy for breast cancer.*

### tmb-msi-to-immunotherapy

*Diagnostic:* tumor mutational burden, tmb-high, tmb-h | *Feature:* tmb-high, tmb-h, msi-h | *Intervention:* pembrolizumab, nivolumab, checkpoint inhibitor

- **PMID 37492581** (2023, 80 cites) — unclassified — immunotherapy — *Immune escape and resistance to immunotherapy in mismatch repair deficient tumors.*
- **PMID 39123202** (2024, 72 cites) — unclassified — car-t, immunotherapy — *Recent developments in immunotherapy for gastrointestinal tract cancers.*
- **PMID 37655661** (2023, 57 cites) — unclassified — immunotherapy — *Linked CD4+/CD8+ T cell neoantigen vaccination overcomes immune checkpoint blockade resistance and enables tumor regress*

### neoantigen-profiling-to-mrna-vaccine

*Diagnostic:* neoantigen prediction, neoantigen profiling, neoantigen identification | *Feature:* neoantigen, neo-antigen, tumor-specific antigen | *Intervention:* mrna vaccine, mrna cancer vaccine, personalized vaccine

- **PMID 37165196** (2023, 1153 cites) — phase1-clinical — immunotherapy, mRNA-vaccine, nanoparticle — *Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer.*
- **PMID 33473101** (2021, 257 cites) — clinical-other — immunotherapy, mRNA-vaccine — *Personalized neoantigen pulsed dendritic cell vaccine for advanced lung cancer.*
- **PMID 38584166** (2024, 189 cites) — phase2-clinical — immunotherapy, mRNA-vaccine — *Personalized neoantigen vaccine and pembrolizumab in advanced hepatocellular carcinoma: a phase 1/2 trial.*

### oncolytic-susceptibility-to-virotherapy

*Diagnostic:* viral receptor expression, nectin-1 expression, cd46 expression | *Feature:* viral entry receptor, nectin-1, cd46 | *Intervention:* t-vec, talimogene, oncolytic herpes

- **PMID 34205379** (2021, 13 cites) — preclinical-invivo — crispr, oncolytic-virus — *Nectin-1 Expression Correlates with the Susceptibility of Malignant Melanoma to Oncolytic Herpes Simplex Virus In Vitro *

## Interpretation

- This is a first-pass pilot covering 6 diagnostic-therapy chains across 4 modalities (radioligands, checkpoint selection, mRNA vaccines, oncolytic viruses).
- The matching rule (intervention required + at least one other link) is conservative; papers that discuss only a diagnostic or only a therapy without the chain are excluded.
- Chain counts depend on keyword coverage and should not be read as exhaustive. Papers using non-standard terminology for diagnostics or interventions may be missed.