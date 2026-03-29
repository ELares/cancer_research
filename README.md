# The Obligate Tradeoff of Resistance

A systematic analysis of 10,413 cancer research articles identifying a non-obvious therapeutic connection: therapy-resistant tumors that switch to oxidative phosphorylation become selectively vulnerable to sonodynamic therapy-triggered ferroptosis.

## The Core Hypothesis

**Resistance escape routes impose obligate biophysical costs.** When tumors survive therapy by switching to OXPHOS (oxidative phosphorylation), they acquire high iron demand, elevated mitochondrial ROS, and lipid-rich membranes — precisely the preconditions for ferroptosis. Sonodynamic therapy (SDT) is the only physical modality that triggers ferroptosis at scale (39 articles) through a traceable ROS → GSH depletion → GPX4 inactivation chain. Yet only 4 papers in 10,413 bridge these two literatures.

**What's genuinely novel** (not what's already known):
1. **A comparison that hasn't been published**: No paper compares physical modalities on ferroptosis engagement. SDT has 39 ferroptosis articles; TTFields, HIFU, and IRE have 0-1 each.
2. **A therapeutic connection that hasn't been made**: OXPHOS-resistance (61 articles) and SDT-ferroptosis (39 articles) overlap in only 4 papers — a 0.04% bridge between two large, separate literatures.

**What's already known** (and we say so):
- Ferroptosis as a cancer vulnerability (2,045-citation review exists)
- SDT triggers ferroptosis via ROS and GSH depletion (dozens of papers)
- Immunogenic cell death following ferroptosis (established)

**Article**: `article/drafts/v1.md` (~14,700 words, 98 verified references, 8 figure placeholders)

**Target journals**: *Trends in Cancer*, *Nature Reviews Cancer*, *Cancer Discovery*

## Corpus

10,413 articles from 1,668 journals (2015-2026), sourced from PubMed (8,220) and Semantic Scholar (2,193).

| Metric | Count |
|--------|-------|
| Total articles | 10,413 |
| Open access | 6,395 (61%) |
| With mechanism tags | 9,684 (93%) |
| With iCite citation metrics | 7,905 (76%) |
| Mechanisms tracked | 19 |
| Cancer types tracked | 22 |
| Phase 3-tier articles | 273 |

## Repository Structure

```
cancer_cure/
├── article/drafts/v1.md           # The manuscript (14,700 words, 98 refs)
├── article/references/            # Verified reference list
│
├── analysis/                      # Data-driven findings
│   ├── hypothesis-sdt-ferroptosis-icd.md    # SDT ferroptosis-ICD hypothesis
│   ├── principle-resistance-tradeoff.md     # The resistance tradeoff principle
│   ├── distilled-hypotheses-final.md        # What survived ruthless scrutiny
│   ├── deep-pattern-analysis.md             # 4 candidate breakthroughs ranked
│   ├── mechanism-matrix.md                  # 19×22 cross-tabulation
│   ├── convergence-map.md                   # Multi-mechanism co-occurrence
│   ├── gap-analysis.md                      # Zero-publication gaps
│   ├── evidence-tiers.md                    # Clinical evidence by mechanism
│   ├── key-findings.md                      # Top 100 articles by impact
│   └── timeline.md                          # 2015-2026 breakthroughs
│
├── corpus/                        # 10,413 articles
│   ├── INDEX.jsonl                # Master index
│   ├── by-pmid/                   # Markdown files with YAML frontmatter
│   └── by-doi/                    # DOI → PMID lookup
│
├── tags/                          # Pre-computed indexes
│   ├── by-mechanism/              # 19 mechanism files
│   ├── by-cancer-type/            # 22 cancer type files
│   ├── by-evidence-level/         # 6 evidence tiers
│   └── by-journal/                # 1,668 journal files
│
├── scripts/                       # Reproducible pipeline
│   ├── fetch_articles.py          # PubMed search + OpenAlex + PMC full text
│   ├── fetch_semantic_scholar.py  # S2 search + citation discovery + TLDR
│   ├── enrich_metadata.py         # PubTator3 + iCite metrics
│   ├── tag_articles.py            # Auto-tag mechanism/cancer/evidence
│   ├── build_index.py             # Rebuild INDEX.jsonl
│   ├── analyze_corpus.py          # Generate analysis files
│   ├── verify_references.py       # Check article refs against corpus
│   ├── config.py                  # API keys, rate limiters, keywords
│   ├── article_io.py              # Shared I/O utilities
│   └── queries.txt                # PubMed search queries
│
├── plans/                         # Research plans and status
├── CLAUDE.md                      # AI assistant project guide
├── requirements.txt               # Python dependencies
└── .env.example                   # API key template
```

## How to Reproduce

```bash
pip install -r requirements.txt
cp .env.example .env  # Add API keys

cd scripts/
python fetch_articles.py --query-file queries.txt --max 500   # PubMed corpus
python fetch_semantic_scholar.py --mode search --max 200       # S2 expansion
python enrich_metadata.py                                      # Annotations + citations
python tag_articles.py                                         # Auto-tag
python build_index.py                                          # Master index
python analyze_corpus.py                                       # Analysis files
python verify_references.py                                    # Reference check
```

## Review History

The manuscript underwent 8 review rounds:

| Round | Focus | Key Changes |
|-------|-------|-------------|
| 1 | Scientific integrity | 16 fixes: false method claims, mechanism errors |
| 2 | Citation spot-check (20 PMIDs) | 7 fixes: wrong authors, journals |
| 3 | Full reference audit (91 refs) | 25 fixes: hallucinated authors/journals |
| 4 | Adversarial peer review | 8 flaws: inflated convergence, missing trials |
| 5 | Citation & evidence audit | Failed trials cited, speculative claims hedged |
| 6 | Falsification review | 6 assumptions strengthened, ICD gap acknowledged |
| 7 | Hypothesis distillation | Eliminated known ideas, sharpened novel claims |
| 8 | Principle extraction | Resistance tradeoff framework, OXPHOS→ferroptosis |

## APIs Used

| API | Purpose | Key Required |
|-----|---------|-------------|
| PubMed E-utilities | Article search + metadata | Free (optional for 10 req/s) |
| Semantic Scholar | Citation graph + TLDR + broader search | Free key |
| OpenAlex | OA status, citations, topics | Free (email for polite pool) |
| PMC BioC | Full-text download | Free |
| PubTator3 | Gene/disease/drug annotations | Free |
| NIH iCite | Citation impact metrics | Free |

## License

Research use. Article draft is pre-publication — do not redistribute without permission.
