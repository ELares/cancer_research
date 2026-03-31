# Taxonomy Rerun Notes

This note summarizes what changed after re-running the tag, index, and analysis pipeline under the tightened taxonomy introduced after PR #16.

## What materially changed

- The new intervention taxonomy is now live in the corpus and index.
- The newly surfaced intervention families are present but still relatively small:
  - `radioligand-therapy`: 52 tagged full-text articles
  - `phagocytosis-checkpoint`: 28
  - `targeted-protein-degradation`: 19
  - `cold-atmospheric-plasma`: 3
- A separate `biology_processes` layer now captures broad biology that would have been too noisy as top-level mechanisms:
  - `tme-stroma`: 1050
  - `autophagy`: 85
  - `senescence-sasp`: 35
  - `cuproptosis`: 7
  - `disulfidptosis`: 2

## What did not materially change

- The top of the mechanism-cancer matrix is still dominated by the same large axes:
  - immunotherapy-lung
  - immunotherapy-melanoma
  - TTFields-glioblastoma
  - immunotherapy-breast
  - immunotherapy-pancreatic
- The highest-priority zero-count gaps remain concentrated in older mechanism/cancer pairs rather than being overturned by the newly added mechanisms.

## Main interpretation

The tightened taxonomy corrected a structural problem in the repo:

- intervention tags are no longer inflated by broad biology/process terms
- resistant-state tags are conservative instead of noisy

The cost of that correction is sparsity. The resistant-state layer currently matches only 10/4830 full-text articles, which is safer than false-positive over-tagging but too thin to support strong corpus-level resistant-state conclusions yet.

## Known taxonomy and query artifacts

- Broad biology terms such as `tme-stroma`, `autophagy`, and `senescence-sasp` were previously liable to act like top-level intervention families. Separating them into `biology_processes` reduced the risk of mistaking background cancer biology coverage for therapeutic modality depth.
- The repo's historical gap lists should not be read as if they were stable under ontology changes. Splitting intervention tags from process tags preserves the high-level dominance of immunotherapy, TTFields, and established device/drug classes, but it weakens confidence in any claim that depended on heterogeneous category counts alone.
- Gap-analysis totals now use unique article counts for mechanisms and cancer types rather than co-tag pair counts. This is the more defensible denominator, but it does change which rows cross the reporting thresholds compared with earlier gap-analysis outputs.
- The resistant-state layer was previously vulnerable to false-positive keyword matches. Tightening it into composite rules revealed that much of the earlier apparent state coverage was likely taxonomy inflation rather than validated resistant-state evidence.
- Some apparent zero-count gaps are already known search artifacts rather than real absences. The clearest confirmed example remains `synthetic-lethality × myeloma`, where independent PubMed verification shows relevant PARP-related literature despite the corpus-level non-detection.
- `radioligand-therapy` remains underrepresented in the fetched full-text corpus relative to the real clinical field. Landmark clinical evidence such as the VISION trial (`PMID 34161051`) is not currently present in the local archive, so `0 phase-labeled trial evidence detected` for radioligands should be treated as a corpus-coverage artifact rather than a biologically meaningful absence.
- See `analysis/landmark-corpus-gaps.md` for a short, manually curated shortlist of missing full-text papers that are important enough to change how absence claims should be phrased.
- The newly added intervention families do matter, but the rerun suggests they are additive rather than revolutionary at the current corpus scale: `radioligand-therapy` (52), `phagocytosis-checkpoint` (28), `targeted-protein-degradation` (19), and `cold-atmospheric-plasma` (3) do not overturn the top-level structure of the mechanism-cancer matrix.

## Implication for issue #10

This rerun completes the first half of the issue:

- expanded taxonomy/query scaffolding is now applied to the corpus
- regenerated analysis outputs now reflect the revised taxonomy

What still remains is the interpretation pass:

- determine which previously reported gaps were true search artifacts
- decide whether the resistant-state rules should be broadened carefully, based on manual validation rather than keyword loosening
