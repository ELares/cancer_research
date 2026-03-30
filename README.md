# Beyond Pharmacologic Ferroptosis Inducers: Physical ROS Modalities for Drug-Tolerant Persister Cells

A cross-literature analysis of 4,830 full-text cancer articles, plus a separate archive of 5,584 abstract-only records, with Monte Carlo biochemical simulation proposing that the persister-ferroptosis field should evaluate physical modalities (PDT, SDT) as spatially selective alternatives to pharmacologic ferroptosis inducers.

## The Idea

Drug-tolerant persister cells are ferroptosis-sensitive (established, PMID:41481741). The field is searching for clinical ferroptosis inducers using *only pharmacologic agents* (RSL3, erastin). Meanwhile, PDT (355 ferroptosis articles) and SDT (121) trigger ferroptosis through ROS + produce immunogenic cell death — advantages systemic drugs lack.

**The proposal**: Physical ROS modalities should be evaluated as persister-targeting tools offering (1) spatial selectivity and (2) ICD for immunotherapy synergy. SDT extends this to deep tumors where PDT can't reach.

**What's novel**: The modality-class question — should the persister field look beyond drugs? — has not been systematically framed despite 355 PDT and 121 SDT ferroptosis papers existing independently.

**Key caveat**: PDT has 40 years of development without demonstrating robust ICD-immune synergy in randomized trials.

## Computational Simulation

A Rust-based Monte Carlo simulation of the ferroptosis cascade (16M cells across 16 conditions) validates the biochemical plausibility:

| Phenotype | Control | RSL3 | SDT/PDT |
|-----------|---------|------|---------|
| Glycolytic | 0.00% | 0.00% | 87.2% |
| OXPHOS | 0.04% | 1.1% | 99.9% |
| **Persister (FSP1↓)** | 1.2% | **42.5%** | **100.0%** |
| Persister+NRF2 | 0.00% | 0.05% | 99.5% |

Key features: autocatalytic LP propagation gated by GSH/GPX4, dynamic GPX4 regulation, FSP1 as GPX4-independent repair pathway. All baselines <2%. Sensitivity: 22/22 holds (100%). RSL3 kills persisters (42.5%) — matching published biology (PMID:41481741). NRF2 rescues from RSL3 but NOT from SDT.

## Article

**Author:** Ezequiel Lares

- `article/drafts/v1.md` — Markdown draft (~11,800 words, 114 references)
- `article/drafts/v1.tex` — LaTeX with proper tables and embedded figures
- `article/drafts/v1.pdf` — Compiled PDF manuscript
- `article/references/bibliography.bib` — BibTeX bibliography (114 entries)

10+ review rounds including adversarial peer review, falsification analysis, novelty assessment, and simulation validation.

## Figures

| Figure | Content |
|--------|---------|
| Fig 1 | Publication trends 2015-2025 |
| Fig 2 | 19×22 mechanism-cancer heatmap |
| Fig 3 | Ferroptosis engagement comparison (χ²=97.3, p<10⁻²³) |
| Fig 4 | Molecular pathway overlap across modalities |
| Fig 5 | Literature disconnect between communities |
| Fig 6 | SDT ferroptosis-ICD chain evidence depth |
| Fig 7 | Monte Carlo simulation results (1M cells/condition) |

## Repository Structure

```
article/drafts/          # Manuscript (md, tex, pdf)
article/figures/         # 7 publication-quality figures (pdf + png)
article/references/      # BibTeX bibliography
corpus/by-pmid/          # 4,830 full-text articles with YAML frontmatter
corpus/abstracts/by-pmid/ # 5,584 abstract-only records
tags/                    # Pre-computed indexes
scripts/                 # Python pipeline (fetch, enrich, tag, analyze, figures)
simulations/             # Rust Monte Carlo ferroptosis simulation
analysis/                # Hypothesis documents and data analysis
books/                   # Reference textbooks (LFS)
plans/                   # Research plans
```

## Reproduction

```bash
# Corpus
pip install -r requirements.txt && cp .env.example .env
cd scripts/
python fetch_articles.py --query-file queries.txt --max 500
python fetch_semantic_scholar.py --mode search --max 200
python enrich_metadata.py && python tag_articles.py
python build_index.py && python analyze_corpus.py
python generate_figures.py

# Simulation
cd ../simulations && cargo build --release
./target/release/ferroptosis-sim

# PDF
cd ../article/drafts && pdflatex v1.tex && bibtex v1 && pdflatex v1.tex && pdflatex v1.tex
```

## License

MIT License. See [LICENSE](LICENSE) for details.
