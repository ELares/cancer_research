# CLAUDE.md

## Project Goal

This repository produces a publishable hypothesis-driven research article from systematic analysis of the cancer therapy literature. The goal is not to summarize oncology — it is to identify non-obvious, evidence-grounded therapeutic hypotheses that could change what scientists test.

## Central Hypothesis

**Therapy-resistant tumors that switch to oxidative phosphorylation (OXPHOS) become selectively vulnerable to sonodynamic therapy (SDT)-triggered ferroptosis, creating a rational post-resistance therapeutic sequence.**

This rests on two novel findings from corpus-wide analysis:
1. SDT is the only physical modality engaging ferroptosis at scale (39 articles vs 0-1 for TTFields/HIFU/IRE) — a comparison nobody has published
2. OXPHOS-resistance (61 articles) and SDT-ferroptosis (39 articles) are bridged by only 4 papers in 10,413 — a major blind spot

## Corpus

- **10,413 articles** from 1,668 journals (2015-2026)
- Sources: PubMed (8,220 via E-utilities) + Semantic Scholar (2,193)
- Enriched with: OpenAlex (OA/citations), PubTator3 (gene/drug NER), iCite (impact metrics)
- 19 therapeutic mechanisms, 22 cancer types, 6 evidence tiers

## Repository Structure

```
corpus/by-pmid/{PMID}.md    # Articles with YAML frontmatter
tags/by-mechanism/*.txt      # Pre-computed PMID lists per mechanism
tags/by-cancer-type/*.txt    # Pre-computed PMID lists per cancer type
corpus/INDEX.jsonl           # Master index (one JSON line per article)
article/drafts/v1.md         # The manuscript (~14,700 words, 98 refs)
analysis/                    # Hypothesis docs and data analysis
scripts/                     # Reproducible pipeline (Python 3.11+)
plans/                       # Research plans and status docs
```

## How to Search the Corpus

```bash
# By mechanism
cat tags/by-mechanism/sonodynamic.txt

# By gene/drug
grep "GPX4" corpus/by-pmid/*.md

# Full text search
grep -l "ferroptosis" corpus/by-pmid/*.md

# By cancer type
cat tags/by-cancer-type/glioblastoma.txt

# High-impact articles (iCite RCR > 10)
grep "icite_rcr: [1-9][0-9]" corpus/by-pmid/*.md

# Specific article
cat corpus/by-pmid/35199647.md

# Master index
cat corpus/INDEX.jsonl | head -5
```

## Scripts

All scripts run from `scripts/` directory:

| Script | Purpose |
|--------|---------|
| `fetch_articles.py` | PubMed search + OpenAlex + PMC full text |
| `fetch_semantic_scholar.py` | S2 search, citation discovery, TLDR enrichment |
| `enrich_metadata.py` | PubTator3 gene/drug NER + iCite citation metrics |
| `tag_articles.py` | Auto-tag by mechanism, cancer type, evidence level |
| `build_index.py` | Rebuild INDEX.jsonl master index |
| `analyze_corpus.py` | Generate all analysis files |
| `verify_references.py` | Check article references against corpus data |

## Article File Format

Each article in `corpus/by-pmid/` has YAML frontmatter:

```yaml
pmid: "35199647"
doi: 10.1172/JCI149258
title: "Tumor Treating Fields dually activate STING..."
journal: "The Journal of clinical investigation"
year: 2022
mechanisms: [immunotherapy, ttfields]
cancer_types: [glioblastoma, melanoma]
evidence_level: preclinical-invivo
genes: [AIM2, STING, cGAS]
drugs: []
icite_rcr: 8.23
```

## Key Design Decisions

- **Automated keyword tagging** (not manual curation) — disclosed honestly in article methods
- **Word-boundary matching** for short keywords (<=4 chars) to prevent false positives
- **`resilient_get()`** with 2 retries and exponential backoff for all HTTP calls
- **Rate limiters per API** in config.py
- **PMID-keyed flat files** for deterministic paths and fast Grep/Glob access
- **Tag index files** for O(1) mechanism/cancer lookups
- **JSONL index** for fast programmatic corpus queries

## What the Article Claims (and Doesn't)

**Claims**:
- SDT is quantitatively unique among physical modalities on ferroptosis (comparison not published elsewhere)
- OXPHOS-resistance and SDT-ferroptosis are connected by only 4 papers (blind spot identified)
- A specific therapeutic sequence (OXPHOS detection → sub-ablative SDT → checkpoint immunotherapy) is rationally grounded but untested

**Does not claim**:
- SDT cures cancer
- The resistance tradeoff is universal or proven
- Preclinical data predicts clinical success (Pexa-Vec precedent cited)
- SDT's ferroptosis engagement is inherent vs nanosonosensitizer-dependent (explicitly flagged)

## Conventions

- Every claim must be traceable to a cited PMID
- Distinguish established evidence from speculation
- All references verified against corpus (98/98 clean)
- Failed trials cited with PMIDs as counterweight
- The article underwent 8 review rounds (see README)
