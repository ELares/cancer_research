# CLAUDE.md

## Author

Ezequiel Lares

## Project Goal

Produce a publishable perspective article proposing that the drug-tolerant persister ferroptosis field should evaluate physical ROS-generating modalities (PDT, SDT) as spatially selective alternatives to pharmacologic ferroptosis inducers. Includes Monte Carlo biochemical simulation as computational validation.

## Central Hypothesis

The persister-ferroptosis field searches only among drugs. Physical ROS modalities (PDT 355 articles, SDT 121) trigger ferroptosis + ICD but haven't been proposed as a persister-targeting class. SDT extends to deep tumors where PDT can't reach.

**What we claim**: The modality-class framing is absent from the literature. SDT as depth-extended PDT for persister ferroptosis is a new connection.
**What we don't claim**: Any individual component is new, or that physical modalities will outperform drugs.

## Corpus

10,413 articles, 1,668 journals, 2015-2026. PubMed (8,220) + Semantic Scholar (2,193). 19 mechanisms, 22 cancer types.

## Simulation

Rust-based Monte Carlo: 16M cells, 4 phenotypes × 4 treatments. Features autocatalytic LP propagation (GSH-gated), dynamic GPX4, FSP1 pathway. All baselines <2%, sensitivity 22/22 holds. Located in `simulations/`.

## Article

`article/drafts/v1.pdf` — 34 pages, 7 figures, 2 tables, 114 references. All citations verified.

## Key Files

```
article/drafts/v1.{md,tex,pdf}    # The manuscript
article/figures/fig*.{pdf,png}      # 7 figures
article/references/bibliography.bib # BibTeX
corpus/by-pmid/{PMID}.md           # 10,413 articles
corpus/INDEX.jsonl                   # Master index
tags/by-mechanism/*.txt              # PMID lists
simulations/src/main.rs              # Rust simulation
simulations/simulation_results.json  # Results
scripts/                             # Python pipeline
books/                               # Reference textbooks (LFS)
```

## Search

```bash
cat tags/by-mechanism/sonodynamic.txt
grep "ferroptosis" corpus/by-pmid/*.md
grep "GPX4" corpus/by-pmid/*.md
```

## Conventions

- Every claim traceable to PMID
- Known findings acknowledged as known
- Nanosonosensitizer confound explicitly flagged
- Failed trials cited (CheckMate-498, BIND-014, Pexa-Vec)
- 114/114 references verified
- Simulation validated: all baselines <2%, sensitivity 22/22
