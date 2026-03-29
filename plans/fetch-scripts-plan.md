# Fetch Scripts Implementation Plan

## Architecture

Four scripts, each standalone but composable:

```
fetch_articles.py  →  enrich_metadata.py  →  tag_articles.py  →  build_index.py
   (search +            (PubTator +           (auto-tag by          (rebuild
    download)            iCite +               mechanism +           INDEX.jsonl +
                         OpenAlex OA)          cancer type)          tag files)
```

Plus a shared `config.py` for API keys, rate limiting, and paths.

## Script 1: fetch_articles.py

**Purpose**: Search PubMed, download abstracts + full text, save as markdown

**Flow**:
1. Accept search query (MeSH term or keyword) + max results
2. PubMed esearch → list of PMIDs
3. PubMed efetch → metadata + abstracts (batch of 200)
4. OpenAlex → OA status + PDF URLs (batch of 50 DOIs)
5. PMC BioC → full text for OA articles
6. Save each article as `corpus/by-pmid/{PMID}.md`
7. Append to `corpus/by-doi/DOI_LOOKUP.jsonl`

**Rate limits to respect**:
- NCBI: 10 req/s with key, 3 req/s without
- OpenAlex: ~10 req/s with polite email
- PMC BioC: ~3 req/s (be conservative)

## Script 2: enrich_metadata.py

**Purpose**: Add annotations to existing articles

**Flow**:
1. Scan `corpus/by-pmid/` for articles missing enrichment
2. PubTator3 → genes, diseases, chemicals, mutations (batch of 100 PMIDs)
3. NIH iCite → citation metrics (batch of 200 PMIDs)
4. Update YAML frontmatter in each article file

## Script 3: tag_articles.py

**Purpose**: Auto-tag articles and build tag index files

**Flow**:
1. Read all articles in corpus
2. Match mechanisms and cancer types from title/abstract/MeSH
3. Write tag files in `tags/by-mechanism/`, `tags/by-cancer-type/`, etc.

## Script 4: build_index.py

**Purpose**: Rebuild master INDEX.jsonl from all corpus files

**Flow**:
1. Glob all `corpus/by-pmid/*.md`
2. Parse YAML frontmatter
3. Write one JSON line per article to `corpus/INDEX.jsonl`

## Dependencies

- `requests` — HTTP calls
- `python-dotenv` — Load .env
- `pyyaml` — Parse/write YAML frontmatter
- `lxml` — Parse PubMed XML responses
- `tqdm` — Progress bars
