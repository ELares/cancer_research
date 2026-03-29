# CLAUDE.md

## Project Goal

Produce a publishable hypothesis article proposing sonodynamic therapy as a physically-triggered ferroptosis inducer for post-resistance cancer therapy.

## Central Hypothesis

Therapy-resistant cancer persister cells are ferroptosis-sensitive (established by others, PMID:41481741). SDT is the only physical modality engaging ferroptosis at scale (39 articles vs 0-1 for TTFields/HIFU/IRE). SDT may be superior to pharmacologic ferroptosis inducers (RSL3, erastin) because it delivers ROS locally and triggers immunogenic cell death for immunotherapy synergy.

**What we claim**: The cross-modality comparison and the SDT-persister connection are new.
**What we don't claim**: We did not discover that resistance creates ferroptosis vulnerability.

## Corpus

10,413 articles, 1,668 journals, 2015-2026. PubMed (8,220) + Semantic Scholar (2,193). 19 mechanisms, 22 cancer types.

## Article

`article/drafts/v1.md` — ~11,000 words, 112 references (all verified), 9 review rounds.

## Key Files

```
corpus/by-pmid/{PMID}.md       # Articles with YAML frontmatter
corpus/INDEX.jsonl              # Master index
tags/by-mechanism/*.txt         # PMID lists per mechanism
article/drafts/v1.md            # The manuscript
analysis/                       # Hypothesis docs and analysis
scripts/                        # Python pipeline
```

## Search

```bash
cat tags/by-mechanism/sonodynamic.txt          # SDT articles
grep "ferroptosis" corpus/by-pmid/*.md         # Ferroptosis mentions
grep "GPX4" corpus/by-pmid/*.md                # GPX4 gene
cat corpus/INDEX.jsonl | head -5               # Index preview
```

## Scripts

| Script | Purpose |
|--------|---------|
| `fetch_articles.py` | PubMed search + OpenAlex + PMC |
| `fetch_semantic_scholar.py` | S2 search + TLDR |
| `enrich_metadata.py` | PubTator3 + iCite |
| `tag_articles.py` | Auto-tag mechanism/cancer/evidence |
| `build_index.py` | Rebuild INDEX.jsonl |
| `analyze_corpus.py` | Generate analysis files |
| `verify_references.py` | Check refs against corpus |

## Conventions

- Every claim traceable to PMID
- Known findings acknowledged as known (ferroptosis vulnerability, adaptive therapy)
- Engineered nanosonosensitizer confound explicitly flagged
- Failed trials cited (CheckMate-498, BIND-014, Pexa-Vec)
- 112/112 references verified against corpus or PubMed
