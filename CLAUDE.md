# CLAUDE.md

## Project Goal

Produce a publishable hypothesis article proposing sonodynamic therapy as a physically-triggered ferroptosis inducer for post-resistance cancer therapy.

## Central Hypothesis

PDT exploits persister-cell ferroptosis via ROS-ICD (355 PubMed articles) but is depth-limited. SDT uses the same mechanism via ultrasound (121 articles) and can reach deep tumors PDT cannot. The question: can SDT extend PDT's ferroptosis-ICD approach to deep-seated tumors?

**What we claim**: SDT as a depth-extended PDT alternative for persister ferroptosis in deep tumors is a connection not made in the literature.
**What we don't claim**: SDT is unique or superior to PDT. PDT dominates on every metric. SDT's only advantage is tissue penetration depth.

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
