# Convergent Therapeutic Mechanisms in Cancer

A systematic analysis of 8,220 cancer research articles across 19 therapeutic mechanisms and 22 cancer types, producing a publishable review article on multi-mechanism convergence in oncology.

## What This Is

This repository contains everything needed to reproduce a large-scale literature analysis of the cancer therapy landscape and a peer-review-ready manuscript arguing that the next frontier lies at convergence zones where mechanistically orthogonal therapies intersect.

**Article**: `article/drafts/v1.md` (~12,700 words, 95 verified references, 7 figure placeholders)

**Target journals**: *Trends in Cancer*, *Nature Reviews Cancer*

## Key Findings

- **8,220 articles** analyzed from 1,507 journals (2015-2026)
- **37% reference 2+ mechanisms**, but subsample analysis estimates only **~9% are experimental combination studies**
- **14x publication growth** from 134 articles (2015) to 1,911 (2025)
- **47 unexplored mechanism-combination pairs** identified (independently verified)
- **Translational bottlenecks**: sonodynamic therapy (503 articles, limited clinical evidence), novel nanoparticle platforms (1,000+ articles, few beyond Phase 1 despite approved nanomedicines like Abraxane/Doxil)
- **Physical-immune interface** identified as the least-explored convergence territory, though clinical synergy evidence remains absent
- **6 failed trials** (CheckMate-498, KEYNOTE-361, BIND-014, Pexa-Vec, CAR-T solid tumor failures, TTFields compliance) discussed as counterweight

## Repository Structure

```
cancer_cure/
├── README.md                  # This file
├── CLAUDE.md                  # AI assistant project guide
├── requirements.txt           # Python dependencies
├── .env.example               # API key template
│
├── article/                   # The manuscript
│   ├── drafts/v1.md           # Current draft (12,700 words, 95 refs)
│   └── references/            # Verified reference list
│
├── corpus/                    # 8,220 research articles
│   ├── INDEX.jsonl            # Master index (one JSON line per article)
│   ├── by-pmid/               # 8,220 markdown files with YAML frontmatter
│   └── by-doi/                # DOI → PMID lookup
│
├── tags/                      # Pre-computed tag indexes
│   ├── by-mechanism/          # 19 mechanism files (PMID lists)
│   ├── by-cancer-type/        # 22 cancer type files
│   ├── by-evidence-level/     # 6 evidence tier files
│   └── by-journal/            # 1,507 journal files
│
├── analysis/                  # Data-driven analysis outputs
│   ├── mechanism-matrix.md    # 19x22 mechanism-cancer cross-tabulation
│   ├── convergence-map.md     # Multi-mechanism co-occurrence patterns
│   ├── gap-analysis.md        # Zero-publication gaps + independent verification
│   ├── evidence-tiers.md      # Clinical evidence level per mechanism
│   ├── key-findings.md        # Top 100 articles by iCite impact
│   └── timeline.md            # 2015-2026 breakthrough timeline
│
├── plans/                     # Research plans and methodology docs
│   ├── project-status-and-next-steps.md
│   ├── research-sources-and-directory-plan.md
│   └── fetch-scripts-plan.md
│
└── scripts/                   # Reproducible pipeline
    ├── fetch_articles.py      # PubMed search + OpenAlex + PMC full text
    ├── enrich_metadata.py     # PubTator3 annotations + iCite metrics
    ├── tag_articles.py        # Auto-tag mechanism/cancer/evidence
    ├── build_index.py         # Rebuild INDEX.jsonl
    ├── analyze_corpus.py      # Generate all analysis files
    ├── verify_references.py   # Check article refs against corpus
    ├── config.py              # API keys, rate limiters, keyword dicts
    ├── article_io.py          # Shared file I/O utilities
    └── queries.txt            # 22 PubMed search queries
```

## How to Reproduce

### 1. Setup

```bash
pip install -r requirements.txt
cp .env.example .env
# Add your API keys to .env (OpenAlex, optionally NCBI, Semantic Scholar)
```

### 2. Build the Corpus

```bash
cd scripts/

# Fetch articles (22 queries × 500 max each, ~35 min)
python fetch_articles.py --query-file queries.txt --max 500

# Enrich with gene/drug annotations + citation metrics (~3 min)
python enrich_metadata.py

# Auto-tag by mechanism, cancer type, evidence level (~30 sec)
python tag_articles.py

# Build master index (~15 sec)
python build_index.py
```

### 3. Generate Analysis

```bash
python analyze_corpus.py
```

### 4. Verify References

```bash
python verify_references.py
```

## Article Review History

The manuscript underwent 6 rounds of increasingly adversarial review:

| Round | Focus | Corrections |
|-------|-------|-------------|
| 1. First-pass review | Scientific integrity, methodology | 16 fixes: false method claims, mechanism errors, journal mismatches |
| 2. Citation spot-check | 20 PMIDs verified against corpus | 7 fixes: wrong authors, wrong journals, claim mischaracterization |
| 3. Full reference audit | All 91 references verified | 25 fixes: 4 wrong authors, 15 wrong journals, 1 wrong year |
| 4. Adversarial peer review | Break the paper mentally | 8 flaws addressed: inflated convergence rate, unverified gaps, missing approved therapies, missing failed trials |
| 5. Citation & evidence audit | Every claim vs its citation | Failed trials cited, citation mismatches fixed, evaluative language softened |
| 6. Falsification review | Try to disprove the hypothesis | 6 assumptions strengthened: ICD translation gap, rational omissions, multiplicative compliance, taxonomy sensitivity |

**Final state**: 95 references (all verified), 12,700 words, honest methodology, balanced thesis with counterarguments.

## Corpus Details

Each article in `corpus/by-pmid/` is a markdown file with structured YAML frontmatter:

```yaml
pmid: "35199647"
doi: 10.1172/JCI149258
title: "Tumor Treating Fields dually activate STING and AIM2..."
authors: [Chen Dongjiang, ...]
journal: "The Journal of clinical investigation"
year: 2022
mechanisms: [immunotherapy, ttfields]
cancer_types: [glioblastoma, melanoma]
evidence_level: preclinical-invivo
genes: [AIM2, STING, cGAS]
icite_rcr: 8.23
```

**Search the corpus**:
```
# By mechanism
cat tags/by-mechanism/ttfields.txt

# By gene
grep "BRAF" corpus/by-pmid/*.md

# Full text search
grep -l "Warburg effect" corpus/by-pmid/*.md

# High-impact articles
grep "icite_rcr: [1-9][0-9]" corpus/by-pmid/*.md
```

## APIs Used

| API | Purpose | Key Required |
|-----|---------|-------------|
| PubMed E-utilities | Article search + metadata | Free (optional key for 10 req/s) |
| OpenAlex | OA status, citations, topics | Free (email for polite pool) |
| PMC BioC | Full-text download | Free |
| PubTator3 | Gene/disease/drug annotations | Free |
| NIH iCite | Citation impact metrics (RCR) | Free |

## License

Research use. Article draft is pre-publication — do not redistribute without permission.
