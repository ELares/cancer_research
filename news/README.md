# News Sources

This directory stores vetted non-peer-reviewed sources (news articles, institutional announcements, expert commentary) that complement the PubMed-indexed corpus.

## Framework

All news sources are evaluated against the criteria in [`analysis/news-source-criteria.md`](../analysis/news-source-criteria.md):

- **Tier 1** (institutional): NIH/NCI, WHO, FDA, journal news sections — cite as evidence
- **Tier 2** (journalism): STAT News, Reuters Health, Endpoints — cite as context with verification
- **Tier 3** (commentary): Expert blogs, advocacy organizations — cite as expert opinion only

## Structure

```
news/
├── README.md           ← this file
├── by-source/          ← articles organized by source domain (populated by scripts/fetch_news.py)
│   ├── cancer-gov/
│   ├── statnews-com/
│   └── ...
└── examples/           ← manually-processed example articles from framework validation
    ├── tier1/
    ├── tier2/
    └── tier3/
```

## Pipeline (Issue #99)

Automated article fetching, claim extraction, verification, and scoring is tracked in issue #99. The criteria document (this issue, #98) defines the rules; the pipeline automates them.
