# Project Status & Next Steps — Cancer Cure Research

**Last updated**: 2026-03-28
**Phase**: Corpus built. Analysis & article writing next.

---

## What Has Been Done

### 1. Repository Setup

- `CLAUDE.md` — Project guide with research scope, questions, conventions
- `.env` / `.env.example` — API keys (OpenAlex, pending Semantic Scholar, NCBI, CORE)
- `.gitignore` — Protects .env, __pycache__, venv
- `requirements.txt` — Python deps: requests, python-dotenv, pyyaml, lxml, tqdm

### 2. Journal & API Verification (64 journals, 14 APIs)

Full results in `plans/research-sources-and-directory-plan.md`.

**Journals verified**: 64 across 7 categories (oncology, general medicine, bioelectricity, nanotechnology, immunotherapy, gene therapy, metabolism). 4 URL corrections applied (Lancet Oncology, 3x Taylor & Francis).

**APIs tested live**:
- **Working (no key)**: PubMed E-utilities, PMC FTP, PMC BioC, PubTator3, ClinicalTrials.gov, OpenAlex, CrossRef, Europe PMC, bioRxiv, NIH iCite
- **Working (key needed)**: Semantic Scholar (free key, pending), CORE (free key, not registered)
- **Working (caveats)**: Unpaywall (rejects @example.com emails; OpenAlex subsumes it)

### 3. Fetch Pipeline — 4 Scripts Built & Tested

All scripts are in `scripts/`. Run from that directory.

| Script | Purpose | Runtime (8K articles) |
|--------|---------|----------------------|
| `fetch_articles.py` | PubMed search + metadata + OpenAlex OA + PMC full text | ~35 min |
| `enrich_metadata.py` | PubTator3 gene/drug annotations + NIH iCite metrics | ~3 min |
| `tag_articles.py` | Auto-tag mechanism, cancer type, evidence level | ~30 sec |
| `build_index.py` | Rebuild INDEX.jsonl master index | ~15 sec |

**Supporting modules**:
- `config.py` — API keys, rate limiters, base URLs, keyword dictionaries, `resilient_get()` with retry
- `article_io.py` — Shared `load_article()`, `save_article()`, `load_frontmatter()` functions

**Usage**:
```bash
cd scripts/

# Full pipeline
python fetch_articles.py --query-file queries.txt --max 500
python enrich_metadata.py
python tag_articles.py
python build_index.py

# Single query
python fetch_articles.py "cancer immunotherapy" --max 100

# Re-enrich specific articles
python enrich_metadata.py --pmids 12345 67890 --force

# Dry-run tagging
python tag_articles.py --dry-run
```

**Key design decisions**:
- Articles stored as PMID-keyed markdown files with YAML frontmatter (`corpus/by-pmid/{PMID}.md`)
- Tag index files in `tags/by-mechanism/`, `tags/by-cancer-type/`, etc. — one PMID per line
- JSONL master index at `corpus/INDEX.jsonl` — one JSON object per line, sorted by year desc
- Word-boundary matching for short keywords (<=4 chars) to prevent false positives
- `resilient_get()` with 2 retries and exponential backoff for all HTTP calls
- Rate limiters per API (NCBI 9/s, OpenAlex 9/s, PMC BioC 2/s, PubTator 2/s, iCite 5/s)

### 4. Corpus — 8,220 Articles

**Search queries used** (in `scripts/queries.txt`):
1. Tumor Treating Fields + cancer
2. Alternating electric fields + neoplasms
3. Sonodynamic therapy + cancer
4. High intensity focused ultrasound + cancer
5. Pulsed electromagnetic field + neoplasms
6. Electrochemical therapy + cancer
7. Electrolytic ablation + tumor
8. Bioelectric + cancer
9. Membrane potential + cancer + therapy
10. Irreversible electroporation + cancer
11. CRISPR + cancer + therapy (2020-2026)
12. CAR-T + cancer (2020-2026)
13. Antibody-drug conjugate + cancer (2023-2026)
14. Bispecific antibody + cancer (2023-2026)
15. mRNA vaccine + cancer (2022-2026)
16. Synthetic lethality + cancer (2020-2026)
17. Oncolytic virus + cancer (2020-2026)
18. Cancer metabolism + warburg (2020-2026)
19. Epigenetic therapy + cancer (2022-2026)
20. Microbiome + cancer + immunotherapy (2022-2026)
21. Combination immunotherapy + cancer (2023-2026)
22. Nanoparticle + drug delivery + cancer (2023-2026)

**Corpus stats**:
| Metric | Count |
|--------|-------|
| Total articles | 8,220 |
| Open access | 4,818 (59%) |
| With PMC full text | ~3,500+ |
| With mechanism tags | 7,912 (96%) |
| With cancer type tags | 5,165 (63%) |
| With evidence level | 3,344 (41%) |
| With iCite RCR | 5,981 (73%) |
| Journals represented | 1,507 |
| Year range | 1950–2026 |

**Mechanism distribution** (articles):
- immunotherapy: 2,982
- nanoparticle: 1,231
- bioelectric: 821
- CAR-T: 776
- antibody-drug-conjugate: 549
- CRISPR: 548
- oncolytic virus: 548
- synthetic lethality: 515
- sonodynamic: 503
- electrochemical therapy: 503
- TTFields: 502
- mRNA vaccine: 486
- HIFU: 440
- bispecific antibody: 434
- epigenetic: 263
- metabolic targeting: 250
- microbiome: 141
- frequency therapy: 110
- electrolysis: 61

**Cancer type distribution** (articles):
- breast: 875, lung: 685, glioblastoma: 661, prostate: 559, pancreatic: 555
- colorectal: 442, liver: 366, melanoma: 320, ovarian: 295, lymphoma: 231
- leukemia: 199, gastric: 162, myeloma: 137, head-and-neck: 135, bladder: 128
- sarcoma: 117, kidney: 116, cervical: 114, mesothelioma: 73, neuroblastoma: 33, esophageal: 27, thyroid: 19

**Evidence level distribution**:
- preclinical-invivo: 1,689
- preclinical-invitro: 1,060
- phase3-clinical: 204
- theoretical: 136
- phase1-clinical: 130
- phase2-clinical: 125

### 5. Article File Format

Each file in `corpus/by-pmid/` has YAML frontmatter + markdown body:

```yaml
---
pmid: "35199647"
doi: 10.1172/JCI149258
pmcid: PMC9012294
title: "Tumor Treating Fields dually activate STING and AIM2..."
authors: [Chen Dongjiang, Le Son B, ...]
journal: "The Journal of clinical investigation"
year: 2022
is_oa: true
oa_status: gold
cited_by_count: 120
mesh_terms: [Glioblastoma, Inflammasomes, ...]
mechanisms: [immunotherapy, ttfields]
cancer_types: [glioblastoma, melanoma, mesothelioma]
evidence_level: preclinical-invivo
genes: [AIM2, STING, cGAS, caspase 1]
drugs: []
icite_rcr: 8.23
icite_percentile: 96.9
---

# Title

## Abstract
...

## Full Text
...
```

### 6. How to Search the Corpus

```
# By mechanism
Read tags/by-mechanism/ttfields.txt
Grep 'mechanisms:.*ttfields' corpus/by-pmid/

# By cancer type
Read tags/by-cancer-type/glioblastoma.txt

# By gene/drug
Grep "BRAF" corpus/by-pmid/

# Full text search
Grep "Warburg effect" corpus/by-pmid/

# By evidence level
Read tags/by-evidence-level/phase3-clinical.txt

# High-impact articles (iCite RCR > 10)
Grep "icite_rcr: [1-9][0-9]" corpus/by-pmid/

# Specific article
Read corpus/by-pmid/35199647.md

# Master index
Read corpus/INDEX.jsonl
```

---

## What Needs to Be Done Next

### Phase 1: Analysis (build in `analysis/`)

**Step 1: Mechanism-Cancer Matrix**
- Create `analysis/mechanism-matrix.md`
- For each mechanism × cancer type pair, summarize: number of articles, highest evidence level, key findings, top-cited papers
- This becomes a core figure in the article

**Step 2: Convergence Map**
- Create `analysis/convergence-map.md`
- Identify articles that span multiple mechanisms (e.g., TTFields + immunotherapy, nanoparticle + CRISPR)
- Map which combinations have evidence and which are unexplored
- Look for synergistic effects reported in multi-mechanism studies

**Step 3: Gap Analysis**
- Create `analysis/gap-analysis.md`
- Identify mechanism × cancer type pairs with zero or very few articles
- Identify promising mechanisms stuck in preclinical that could advance
- Find patterns: which mechanisms work across many cancer types vs. niche ones

**Step 4: Evidence Tier Summary**
- Create `analysis/evidence-tiers.md`
- For each mechanism, list the highest-evidence findings (Phase 3 > Phase 2 > Phase 1 > preclinical)
- Highlight which mechanisms have clinical validation vs. still theoretical

**Step 5: Key Findings Extraction**
- For the top 50-100 highest-impact articles (by iCite RCR), extract:
  - Primary finding
  - Mechanism of action
  - Efficacy data (response rates, survival, etc.)
  - Limitations
- Store in `analysis/key-findings.md`

**Step 6: Timeline of Breakthroughs**
- Create `analysis/timeline.md`
- Map major milestones chronologically: FDA approvals, landmark trials, first-in-class discoveries
- Focus on 2020-2026 for recency

### Phase 2: Article Writing (build in `article/`)

**Target journal**: Perspective/Review article for *Trends in Cancer*, *Nature Reviews Cancer*, or *Cancer Research* (review track).

**Proposed article structure**:

```
article/drafts/v1.md

1. Abstract (250 words)
   - Context: fragmented cancer research landscape
   - Thesis: convergent mechanisms offer curative potential
   - Key findings from our analysis
   - Implications

2. Introduction
   - Current state of cancer therapy
   - The case for looking beyond single-mechanism approaches
   - Scope of this review

3. Methods
   - Systematic search strategy (PubMed, PMC, OpenAlex)
   - Inclusion/exclusion criteria
   - Categorization framework (mechanisms, cancer types, evidence levels)
   - Analysis approach

4. Results / Analysis
   4.1 Physical destruction mechanisms (TTFields, HIFU, sonodynamic, frequency, electrolysis)
   4.2 Bioelectric modulation (membrane potential, electroporation)
   4.3 Immune-based approaches (immunotherapy, CAR-T, bispecific, mRNA vaccines, oncolytic)
   4.4 Precision molecular approaches (CRISPR, synthetic lethality, ADCs, epigenetic)
   4.5 Metabolic & microenvironment approaches (metabolic targeting, microbiome, nanoparticle delivery)
   4.6 Cross-mechanism convergence patterns
   4.7 The mechanism-cancer matrix: what works where

5. Discussion
   5.1 Underexplored intersections with high potential
   5.2 Barriers to clinical translation
   5.3 The case for combinatorial/multi-modal strategies
   5.4 Frequency-based and bioelectric approaches: an overlooked frontier

6. Conclusion
   - Summary of convergent opportunities
   - Call to action for multi-disciplinary research

7. Figures
   - Fig 1: Mechanism-cancer matrix heatmap
   - Fig 2: Evidence level pyramid per mechanism
   - Fig 3: Convergence network diagram
   - Fig 4: Timeline of breakthroughs 2020-2026

8. References (BibTeX → journal format)
```

**Article conventions**:
- Every claim cited with PMID
- Distinguish established vs. emerging vs. theoretical
- Use precise scientific terminology
- Maintain objectivity
- Target ~8,000-12,000 words for a comprehensive review

### Phase 3: Polish & Submit

- Generate BibTeX bibliography from cited PMIDs
- Create publication-quality figures
- Format per target journal guidelines
- Write cover letter
- Submit

---

## Pending Items / Blockers

| Item | Status | Action Needed |
|------|--------|---------------|
| Semantic Scholar API key | Pending approval | Add to .env when received |
| NCBI API key | Not registered | Register at ncbi.nlm.nih.gov/account for 10 req/s |
| CORE API key | Not registered | Register at core.ac.uk/services/api (optional) |
| Corpus expansion | Optional | Could add more queries for specific gaps found in analysis |
| Preprint coverage | Not done | Could add Europe PMC / bioRxiv queries for cutting-edge work |

---

## Repository Structure (Final)

```
cancer_cure/
├── CLAUDE.md                          # Project guide
├── .env                               # API keys (gitignored)
├── .env.example                       # Template
├── .gitignore
├── requirements.txt                   # Python dependencies
│
├── plans/
│   ├── research-sources-and-directory-plan.md   # Journal list, APIs, directory design
│   ├── fetch-scripts-plan.md                    # Script architecture
│   └── project-status-and-next-steps.md         # THIS FILE
│
├── corpus/
│   ├── INDEX.jsonl                    # Master index (8,220 entries)
│   ├── by-pmid/                       # 8,220 article files
│   │   ├── 26670971.md
│   │   ├── 27668386.md
│   │   └── ... (8,220 files)
│   └── by-doi/
│       └── DOI_LOOKUP.jsonl           # DOI → PMID mapping
│
├── tags/
│   ├── by-mechanism/                  # 19 mechanism tag files
│   │   ├── immunotherapy.txt (2,982 PMIDs)
│   │   ├── nanoparticle.txt (1,231 PMIDs)
│   │   └── ...
│   ├── by-cancer-type/                # 22 cancer type tag files
│   │   ├── breast.txt (875 PMIDs)
│   │   ├── lung.txt (685 PMIDs)
│   │   └── ...
│   ├── by-evidence-level/             # 6 evidence level tag files
│   └── by-journal/                    # 1,507 journal tag files
│
├── analysis/                          # TO BE BUILT — Phase 1
│   ├── mechanism-matrix.md
│   ├── convergence-map.md
│   ├── gap-analysis.md
│   ├── evidence-tiers.md
│   ├── key-findings.md
│   ├── timeline.md
│   └── notes/
│
├── article/                           # TO BE BUILT — Phase 2
│   ├── drafts/
│   │   └── v1.md
│   ├── figures/
│   ├── references/
│   │   └── bibliography.bib
│   └── supplementary/
│
└── scripts/
    ├── config.py                      # Shared config, rate limiters, keywords
    ├── article_io.py                  # Shared load/save functions
    ├── fetch_articles.py              # PubMed + OpenAlex + PMC fetcher
    ├── enrich_metadata.py             # PubTator + iCite enrichment
    ├── tag_articles.py                # Auto-tagger + index builder
    ├── build_index.py                 # INDEX.jsonl rebuilder
    └── queries.txt                    # 22 PubMed search queries
```

---

## Key Technical Notes for Future Context

1. **Python 3.11+** required (uses `dict | None` type hints)
2. All scripts run from `scripts/` directory (`cd scripts/ && python fetch_articles.py ...`)
3. Rate limiters are global singletons in `config.py` — safe for single-process use only
4. The `resilient_get()` function retries on 5xx and connection errors (2 retries, exponential backoff)
5. Keyword matching uses `\b` word boundaries for keywords <=4 chars to prevent substring false positives (e.g., "all" no longer matches "overall")
6. Dangerous short keywords were removed: "ire", "all", "adc", "bite", "rct", "brca"
7. OpenAlex is our primary OA discovery tool (replaces Unpaywall)
8. PMC BioC API returns full text as structured JSON — we extract text passages excluding REF/SUPPL/TABLE/FIG sections
9. PubTator3 provides gene/disease/chemical/mutation NER annotations per PMID
10. iCite provides Relative Citation Ratio (field-normalized impact) and NIH percentile
11. The `oa_url` field is normalized to empty string (not null) for consistent Grep behavior
12. Enrichment script respects `--skip-pubtator` / `--skip-icite` flags in its filter logic
