# Evidence Labeling Adjudication Protocol

Decision rules for resolving disagreements between primary and second raters. Commit this document BEFORE any second-rater labeling begins.

## Scope

Applies to the 60-article second-rater overlap subset (24% of the 250-article v2 gold set, 6 per mechanism). Both raters must have read and agreed to the labeling guidelines (evidence-labeling-guidelines-v2.md) before labeling.

## Decision Rules

### 1. Agreement
Both raters assign the same tier → **label is that tier**. No further action needed.

### 2. Adjacent Disagreement
Raters assign tiers one level apart in the ordinal hierarchy (e.g., preclinical-invivo vs preclinical-invitro, or phase1-clinical vs clinical-other).

→ **Primary rater's label stands** unless the second rater provides a specific textual justification citing a passage in the article. If a cited passage supports the second rater's tier, relabel accordingly and note in `gold_notes`.

### 3. Non-Adjacent Disagreement
Raters assign tiers 2+ levels apart (e.g., phase2-clinical vs preclinical-invitro, or theoretical vs preclinical-invivo).

→ **Flag for review.** Both raters re-read the abstract and methods section independently, then discuss. If agreement is reached, use that tier. If still no agreement, assign the **more conservative (lower) tier** and document the disagreement in `gold_notes`.

### 4. Evidence vs None-Applicable
One rater assigns an evidence tier, the other assigns none-applicable.

→ Re-read the abstract and methods. Apply this test: **does the article report ANY primary research data** (even a single experiment, case report, or computational result)? If yes → assign the appropriate evidence tier. If the article only discusses, reviews, or synthesizes other work → assign none-applicable.

### 5. Tie-Breaking (Last Resort)
If discussion fails to resolve after Steps 2-4, the **more conservative (lower) tier** wins. Document the disagreement in `gold_notes` with both raters' tiers and reasoning.

## Ordinal Hierarchy (for adjacency determination)

From highest to lowest:
1. phase3-clinical
2. phase2-clinical
3. phase1-clinical
4. clinical-other
5. preclinical-invivo
6. preclinical-invitro
7. theoretical
8. none-applicable

Adjacent = 1 step apart. Non-adjacent = 2+ steps apart.

## Test-Retest Fallback

If a second rater is unavailable, the primary rater relabels the 60-article overlap subset after a **minimum 7-day delay**. This measures test-retest reliability rather than inter-rater agreement. Report as "test-retest kappa" rather than "inter-rater kappa" and note the single-rater limitation.

## Metrics to Compute

After both raters complete the overlap subset:
- Cohen's kappa (unweighted) — strict categorical agreement
- Cohen's kappa (linear weighted) — ordinal-aware, penalizes distant disagreements less
- Per-tier agreement rate (only for tiers with N ≥ 5 in the overlap)
- Confusion matrix (primary × second rater)

Run: `python scripts/compute_inter_rater_agreement.py`
