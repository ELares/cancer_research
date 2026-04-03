# Evidence Gold-Set Evaluation

Manual labels are present for **100** sampled records from `analysis/evidence-gold-set-v1.csv`, with label assignments stored in `analysis/evidence-gold-labels-v1.csv`.

- Sampling design: 10 `predicted-tagged` and 10 `predicted-untagged` rows for each of `immunotherapy`, `mRNA-vaccine`, `electrochemical-therapy`, `ttfields`, and `synthetic-lethality`.
- All sampled rows currently have manual labels.

## Overall Metrics

- Exact-label accuracy: **46/100 (46.0%)**
- Binary evidence-detection precision: **96.0%** (48/50)
- Binary evidence-detection recall: **55.2%** (48/87)
- Binary evidence-detection F1: **0.701**
- Gold positive rows: **87**; predicted positive rows: **50**

- Exact-label scoring treats an empty heuristic prediction as equivalent to `none-applicable`, because both represent an intentional no-evidence assignment.
## Per-Label Metrics

| Label | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| **phase3-clinical** | 2 | 2 | 0 | 50.0% | 100.0% | 0.667 |
| **phase2-clinical** | 0 | 1 | 0 | 0.0% | 0.0% | 0.000 |
| **phase1-clinical** | 3 | 2 | 0 | 60.0% | 100.0% | 0.750 |
| **clinical-other** | 3 | 1 | 10 | 75.0% | 23.1% | 0.353 |
| **preclinical-invivo** | 19 | 4 | 10 | 82.6% | 65.5% | 0.731 |
| **preclinical-invitro** | 7 | 5 | 13 | 58.3% | 35.0% | 0.438 |
| **theoretical** | 1 | 0 | 19 | 100.0% | 5.0% | 0.095 |
| **none-applicable** | 11 | 39 | 2 | 22.0% | 84.6% | 0.349 |

## Per-Mechanism Exact Accuracy

| Mechanism | Labeled rows | Exact accuracy | Predicted positive | Gold positive |
|---|---|---|---|---|
| **electrochemical-therapy** | 20 | 55.0% | 10 | 17 |
| **immunotherapy** | 20 | 55.0% | 10 | 16 |
| **mRNA-vaccine** | 20 | 35.0% | 10 | 18 |
| **synthetic-lethality** | 20 | 35.0% | 10 | 20 |
| **ttfields** | 20 | 50.0% | 10 | 16 |

## Most Common Confusions

- **theoretical -> none-applicable**: 15  
  Example PMIDs: 41660849 (immunotherapy), 35265614 (mRNA-vaccine), 35818395 (mRNA-vaccine)
- **clinical-other -> none-applicable**: 10  
  Example PMIDs: 39830952 (immunotherapy), 40753395 (immunotherapy), 41313664 (immunotherapy)
- **preclinical-invitro -> none-applicable**: 9  
  Example PMIDs: 21878233 (electrochemical-therapy), 36671876 (electrochemical-therapy), 30101194 (synthetic-lethality)
- **preclinical-invivo -> none-applicable**: 5  
  Example PMIDs: 37655661 (immunotherapy), 38858600 (immunotherapy), 37428918 (mRNA-vaccine)
- **preclinical-invivo -> preclinical-invitro**: 3  
  Example PMIDs: 39679828 (immunotherapy), 40069686 (immunotherapy), 39380383 (mRNA-vaccine)
- **theoretical -> preclinical-invivo**: 2  
  Example PMIDs: 40406146 (immunotherapy), 38043609 (synthetic-lethality)
- **theoretical -> preclinical-invitro**: 2  
  Example PMIDs: 39403328 (mRNA-vaccine), 39305483 (synthetic-lethality)
- **preclinical-invitro -> preclinical-invivo**: 2  
  Example PMIDs: 27693939 (electrochemical-therapy), 36358594 (ttfields)
- **none-applicable -> phase2-clinical**: 1  
  Example PMIDs: 35890409 (mRNA-vaccine)
- **preclinical-invivo -> phase1-clinical**: 1  
  Example PMIDs: 40519325 (mRNA-vaccine)
- **preclinical-invivo -> clinical-other**: 1  
  Example PMIDs: 30071778 (electrochemical-therapy)
- **preclinical-invitro -> phase3-clinical**: 1  
  Example PMIDs: 29284495 (ttfields)

## Interpretation

- The current evidence tagger behaves like a conservative detector: it rarely assigns evidence to rows manually labeled `none-applicable`, but it misses a large share of valid evidence-bearing records.
- The largest blind spots in this sample are unlabeled `theoretical`, `clinical-other`, and `preclinical-invitro` studies. That lines up with earlier qualitative audit notes.
- The gold set supports using coverage-aware manuscript language. The current heuristic is much more reliable for `if tagged, usually real` than for `if untagged, probably absent`.
