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

## Implication for issue #10

This rerun completes the first half of the issue:

- expanded taxonomy/query scaffolding is now applied to the corpus
- regenerated analysis outputs now reflect the revised taxonomy

What still remains is the interpretation pass:

- determine which previously reported gaps were true search artifacts
- decide whether the resistant-state rules should be broadened carefully, based on manual validation rather than keyword loosening
