# Depth-Extended Ferroptosis-ICD Therapy: Can SDT Reach Tumors That PDT Cannot?

A cross-literature analysis of 10,413 cancer research articles examining whether sonodynamic therapy can extend photodynamic therapy's ferroptosis-ICD mechanism to deep-seated tumors.

## The Hypothesis

PDT already exploits persister-cell ferroptosis vulnerability through ROS-mediated ferroptosis and immunogenic cell death (355 PubMed ferroptosis articles, 67 ferroptosis+ICD). But PDT is limited to superficial tumors (light penetrates millimeters).

SDT uses ultrasound instead of light to trigger the same mechanism (121 PubMed ferroptosis articles). Ultrasound penetrates centimeters. For deep tumors (pancreatic, hepatic, pelvic), SDT may deliver ferroptosis-ICD where PDT physically cannot.

**This is an incremental advance, not a paradigmatic one.** The mechanism is shared; the advantage is depth.

**What's genuinely novel**: The connection of the emerging persister-ferroptosis biology to SDT as a depth-extended alternative to PDT for deep tumors.

**What's already known**: Persister ferroptosis sensitivity (PMID:41481741), PDT-ferroptosis-ICD (355 papers), SDT-ferroptosis (121 papers), resistance tradeoffs (Gatenby, collateral sensitivity).

**Key caveat**: PDT has been in development 40 years without demonstrating robust ICD-immune synergy clinically. SDT must outperform this baseline.

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
