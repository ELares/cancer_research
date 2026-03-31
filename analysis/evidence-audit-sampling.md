# Evidence Audit Sampling Notes

This note documents a manual spot-check of uncategorized primary-study-like records after the taxonomy rerun and evidence-coverage audit refresh in PR #17.

## Sample reviewed

| PMID | Mechanism(s) | Current status | Manual readout |
|---|---|---|---|
| 27298410 | immunotherapy, oncolytic-virus | uncategorized | Explicit phase Ib/IV clinical trial wording. This is a classifier blind spot and should be tagged as clinical evidence. |
| 33473101 | immunotherapy, mRNA-vaccine | uncategorized | Investigator-initiated single-arm pilot clinical study. Real clinical evidence, but not phase-labeled. This is a schema blind spot, not just a missing keyword. |
| 34861036 | mRNA-vaccine | uncategorized | Generic clinical trial / vaccine-response study in patients, without phase labeling. Schema blind spot. |
| 35661819 | immunotherapy, mRNA-vaccine | uncategorized | Single-patient clinical report / case-style translational report. Schema blind spot. |
| 25442132 | electrochemical-therapy | uncategorized | Consensus/standards paper. Should remain outside the evidence tiers. |
| 33080774 | ttfields | uncategorized | Mechanistic lab study with clear preclinical orientation, but not announced with simple `in vitro` keywords. Likely classifier miss. |
| 23095807 | ttfields | uncategorized | Case-report / small clinical series. Schema blind spot under the current six-tier model. |
| 35444283 | crispr, synthetic-lethality | uncategorized | Preclinical screen/model paper with tumor-regression language. Likely classifier miss. |

## What this sample suggests

- Some residual misses are true classifier misses and can be fixed with better metadata/text handling.
  Examples: explicit phase-designation in `pub_types`; preclinical model papers that do not literally say `in vivo` or `in vitro`.
- A larger share of the remaining uncategorized pool appears to be outside the current six-tier schema.
  Examples: `Clinical Trial` without phase labeling, feasibility studies, pilot studies, retrospective outcome studies, and case reports.
- Consensus/standards documents also remain intentionally outside the evidence tiers and should not be treated as tagging failures.

## Practical implication

The next meaningful decision is not just "add more keywords." It is whether the repo wants to keep the current six-tier evidence model, with coverage-aware caveats, or introduce a broader `clinical-other` bucket for non-phase patient studies and case-series evidence.
