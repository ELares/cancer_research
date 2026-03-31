# CLAUDE.md

## Author

Ezequiel Lares

## Repo Goal

Maintain this repository as an open cancer-research synthesis workspace that compares therapeutic mechanisms, evidence depth, resistant-state biology, pathway targets, and simulation-based hypotheses without assuming a single answer in advance.

The repo began around a PDT/SDT-persister-ferroptosis thesis. That thesis is still worth evaluating, but it is no longer the only organizing frame. The repo should stay open to alternative modalities, targets, pathways, and interpretations when the corpus or external literature supports them.

## Working Principle

Do not let the current manuscript title or earlier repo framing hard-code the conclusion.

When reviewing or extending the work:

- treat PDT/SDT as one candidate lane, not the default winner
- prefer resistant-state and evidence-quality questions over modality loyalty
- distinguish intervention classes from broad biology/process terms
- treat corpus-level non-detection as provisional unless externally verified
- surface known artifacts, missing landmark papers, and tagging limitations early

## Current Workstreams

- manuscript drafting and revision
- corpus fetching, enrichment, tagging, and indexing
- evidence-tier audits and coverage caveats
- taxonomy and search refinement
- pathway-target and resistant-state analysis
- simulation work around ferroptosis and escape mechanisms
- broader strategy review of alternative therapies and biological bottlenecks

## Current Repo State

- local full-text corpus: 4,830 records
- abstract-only archive: 5,584 records
- mechanism taxonomy, evidence tiers, pathway-targets, biology-process tags, and resistant-state scaffolding are all active
- evidence tagging is improved but still incomplete
- some landmark papers are known to be missing from the local full-text archive

## What To Optimize For

- claims that are traceable and caveated
- taxonomy choices that do not inflate conclusions
- docs and manuscript language that reflect uncertainty honestly
- outputs that help compare alternatives fairly
- maintainable scripts and reproducible reruns

## What To Avoid

- assuming the repo exists only to defend PDT/SDT
- writing stronger absence claims than the evidence model supports
- confusing patient-study signal with phase-labeled trial maturity
- treating broad process coverage as proof of therapeutic depth
- letting historical framing survive after the underlying analysis changed

## Key Files

```text
README.md                                 repo-level framing and entry points
article/drafts/v1.{md,tex,pdf}            manuscript drafts
analysis/evidence-coverage-audit.md       evidence-tier coverage and guardrails
analysis/taxonomy-rerun-notes.md          taxonomy/query caveats after reruns
analysis/pathway-target-audit.md          pathway-target coverage
analysis/landmark-corpus-gaps.md          known missing papers that distort claims
corpus/INDEX.jsonl                        master index
scripts/                                  Python pipeline
simulations/                              Rust simulation work
tags/                                     precomputed tag indexes
```

## Search Conventions

Prefer fast repo-native inspection first:

```bash
rg "term" scripts analysis article
rg --files corpus/by-pmid | head
sed -n '1,120p' analysis/evidence-coverage-audit.md
```

## Writing Conventions

- every strong claim should be traceable to the corpus, analysis outputs, or external verification
- separate generated outputs from handwritten interpretation notes
- use coverage-aware language such as `not detected in the local keyword-derived analysis` where appropriate
- if a known taxonomy artifact or corpus gap applies, mention it directly rather than burying it
- keep the repo open to thesis revision rather than optimizing for rhetorical neatness
