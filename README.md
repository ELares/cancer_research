# Sonodynamic Therapy as a Physically-Triggered Ferroptosis Inducer for Post-Resistance Cancer Therapy

A cross-literature analysis of 10,413 cancer research articles identifying SDT as uniquely positioned among physical modalities to exploit the ferroptosis vulnerability of therapy-resistant persister cells.

## The Hypothesis

Recent work (Higuchi et al., *Science Advances* 2026) demonstrates that therapy-resistant cancer persister cells are sensitized to ferroptosis. The question our analysis addresses: **which therapeutic modality best exploits this vulnerability?**

Pharmacologic ferroptosis inducers (RSL3, erastin) lack spatial selectivity. Our corpus-wide analysis identifies **sonodynamic therapy (SDT)** as a candidate alternative — the only physical modality that engages ferroptosis at scale (39 articles, vs 0-1 for TTFields/HIFU/IRE), with the added advantage of local ROS delivery via focused ultrasound and generation of immunogenic cell death that could synergize with checkpoint immunotherapy.

**What's genuinely novel** (verified):
1. The cross-modality comparison on ferroptosis engagement has not been published
2. The connection of persister-ferroptosis biology to SDT specifically has not been made

**What's already known** (and we say so):
- Persister cells are ferroptosis-sensitive (PMID:41481741, Science Advances 2026)
- SDT triggers ferroptosis via ROS/GSH depletion (dozens of nanoparticle papers)
- Resistance creates vulnerabilities (Gatenby adaptive therapy, collateral sensitivity)

**Key caveat**: Most SDT-ferroptosis data involves engineered nanosonosensitizers designed to deplete GSH, not inherent SDT properties. This confound must be resolved experimentally.

## Decisive Experiment

Compare SDT vs RSL3 vs erastin for killing persister cells (using existing models from PMID:41481741). Measure both cell death AND ICD markers (calreticulin, HMGB1, STING). If SDT produces equivalent killing with superior ICD, the hypothesis is supported.

## Article

`article/drafts/v1.md` — ~11,000 words, 112 verified references, underwent 9 review rounds including adversarial peer review and falsification analysis.

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
