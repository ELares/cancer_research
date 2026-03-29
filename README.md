# Physical ROS-Generating Modalities as Spatially Selective Ferroptosis Inducers for Drug-Tolerant Persister Cells

A cross-literature analysis of 10,413 cancer articles proposing that the persister-ferroptosis field should evaluate physical modalities (PDT, SDT) as spatially selective alternatives to pharmacologic ferroptosis inducers.

## The Idea

Drug-tolerant persister cells are ferroptosis-sensitive (established, PMID:41481741). The field is searching for clinical ferroptosis inducers using *only pharmacologic agents* (RSL3, erastin). Meanwhile, PDT (355 ferroptosis articles) and SDT (121) trigger ferroptosis through ROS + produce immunogenic cell death — advantages systemic drugs lack.

**The proposal**: Physical ROS modalities should be evaluated as persister-targeting tools offering (1) spatial selectivity and (2) ICD for immunotherapy synergy. SDT extends this to deep tumors where PDT can't reach.

**What's novel**: The modality-class question — should the persister field look beyond drugs? — has not been systematically framed despite 355 PDT and 121 SDT ferroptosis papers existing independently.

**What's known**: Persister ferroptosis (44 papers), PDT-ferroptosis-ICD (67 papers), SDT-ferroptosis (121 papers). The individual components are published; the cross-community connection is absent.

**Key caveat**: PDT has 40 years of development without demonstrating robust ICD-immune synergy in randomized trials. The pharmacologic-to-physical translation may not improve outcomes.

## Decisive Experiment

Compare SDT vs RSL3 vs erastin for killing persister cells (using existing models from PMID:41481741). Measure both cell death AND ICD markers (calreticulin, HMGB1, STING). If SDT produces equivalent killing with superior ICD, the hypothesis is supported.

## Article

**Author:** Ezequiel Lares

- `article/drafts/v1.md` — Markdown draft (~11,200 words, 114 verified references)
- `article/drafts/v1.tex` — LaTeX version with `\cite{}` references
- `article/references/bibliography.bib` — BibTeX bibliography (114 entries)

10+ review rounds including adversarial peer review, falsification analysis, and novelty assessment with honest prior-art acknowledgment.

## Corpus

10,413 articles from 1,668 journals (2015-2026). Sources: PubMed (8,220) + Semantic Scholar (2,193). Enriched with OpenAlex, PubTator3, and iCite.

## Repository Structure

```
article/drafts/v1.md           # The manuscript
analysis/                      # Hypothesis documents + data analysis
corpus/by-pmid/                # 10,413 articles with YAML frontmatter
tags/                          # Pre-computed indexes
scripts/                       # Reproducible Python pipeline
plans/                         # Research plans
```

## Reproduction

```bash
pip install -r requirements.txt && cp .env.example .env
cd scripts/
python fetch_articles.py --query-file queries.txt --max 500
python fetch_semantic_scholar.py --mode search --max 200
python enrich_metadata.py && python tag_articles.py
python build_index.py && python analyze_corpus.py
python verify_references.py
```

## Review History

9 rounds: scientific integrity → citation verification (112/112 clean) → adversarial peer review → falsification → hypothesis distillation → novelty assessment with honest acknowledgment of prior art.

## License

Research use. Pre-publication — do not redistribute without permission.
