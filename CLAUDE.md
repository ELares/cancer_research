# CLAUDE.md

## Project Goal

This repository is a research project aimed at producing a publishable scientific journal article. The core mission is to systematically analyze the global landscape of cancer research — across all cancer types — and identify promising, underexplored, or convergent therapeutic mechanisms that could lead to novel curative approaches.

## Research Scope

### Primary Research Questions

1. **What are the most significant recent findings across global cancer research?** — Systematic review of high-impact journals, preprints, and clinical trial data across all cancer types.
2. **Can natural frequencies (resonance-based therapies) selectively destroy cancer cells?** — Investigate tumor-treating fields (TTFields), piezoelectric nanoparticles, ultrasound-mediated therapies, and bioelectromagnetic approaches.
3. **Can electrolysis or bioelectric modulation be used as a therapeutic mechanism?** — Explore electrolytic tumor ablation, bioelectric signaling reprogramming, and electrochemical therapy (EChT).
4. **What novel targeted mechanisms are emerging that challenge conventional paradigms?** — Identify approaches beyond traditional chemo/immunotherapy: synthetic lethality, metabolic reprogramming, epigenetic editing, oncolytic viruses, microbiome-mediated therapies, and others.
5. **Where do these mechanisms converge?** — Look for combinatorial or synergistic strategies that could yield breakthrough results.

### Cancer Types

All cancer types are in scope. Cross-cancer pattern analysis is a priority — mechanisms that work across multiple cancer types are of highest interest.

### Key Therapeutic Domains to Investigate

- **Frequency-based therapies**: TTFields (Optune), resonant frequency destruction, pulsed electromagnetic fields (PEMF), focused ultrasound (HIFU/FUS), sonodynamic therapy
- **Electrolysis & bioelectric approaches**: Electrochemical therapy (EChT), bioelectric membrane potential manipulation, iontophoresis for drug delivery, irreversible electroporation (IRE)
- **Targeted & emerging mechanisms**: CRISPR-based gene editing, CAR-T and next-gen cell therapies, antibody-drug conjugates (ADCs), bispecific antibodies, mRNA cancer vaccines, synthetic lethality (PARP inhibitors, etc.), metabolic targeting (Warburg effect exploitation), epigenetic reprogramming, oncolytic virotherapy, nanoparticle-mediated delivery, microbiome modulation
- **Convergent / combinatorial strategies**: Multi-modal approaches that combine the above

## Output: Publishable Research Article

The end product is a scientific review/perspective article suitable for submission to a peer-reviewed journal. The article should:

- Follow standard scientific article structure (Abstract, Introduction, Methods, Results/Analysis, Discussion, Conclusion, References)
- Be rigorously sourced with proper citations (APA or journal-specific format)
- Present original synthesis — not just a literature dump, but novel analysis of patterns, gaps, and opportunities
- Include clear figures/diagrams where appropriate
- Target journals such as: *Nature Reviews Cancer*, *The Lancet Oncology*, *Cancer Research*, *Trends in Cancer*, or similar

## Repository Structure

```
cancer_cure/
├── CLAUDE.md                          # This file (project guide)
├── plans/                             # Research plans and methodology
│
├── corpus/                            # ALL downloaded articles live here
│   ├── INDEX.jsonl                    # Master index (one JSON line per article)
│   ├── by-pmid/                       # Article content as .md with YAML frontmatter
│   │   └── {PMID}.md                  # e.g., 38000001.md
│   └── by-doi/
│       └── DOI_LOOKUP.jsonl           # DOI → PMID mapping
│
├── tags/                              # Pre-computed tag indexes (PMID lists)
│   ├── by-mechanism/                  # e.g., ttfields.txt, car-t.txt, crispr.txt
│   ├── by-cancer-type/                # e.g., glioblastoma.txt, breast.txt
│   ├── by-evidence-level/             # e.g., phase3-clinical.txt, preclinical-invivo.txt
│   └── by-journal/                    # e.g., nature-cancer.txt
│
├── analysis/                          # Cross-cutting synthesis and notes
│   └── notes/
│
├── article/                           # The manuscript
│   ├── drafts/
│   ├── figures/
│   ├── references/                    # bibliography.bib
│   └── supplementary/
│
└── scripts/                           # Fetch, enrich, index, tag automation
```

### Article File Format

Each article in `corpus/by-pmid/` is a markdown file with YAML frontmatter for structured search:

```yaml
---
pmid: "38000001"
doi: "10.1038/..."
title: "..."
authors: ["Smith J", "Lee K"]
journal: "Nature Cancer"
year: 2024
mechanisms: ["ttfields", "immunotherapy"]
cancer_types: ["glioblastoma"]
evidence_level: "phase2-clinical"
mesh_terms: ["..."]
genes: ["PD-L1"]
drugs: ["pembrolizumab"]
date_added: "2026-03-28"
---
```

### How to Search the Corpus

- **By mechanism**: `Grep 'mechanisms:.*ttfields' corpus/by-pmid/` or `Read tags/by-mechanism/ttfields.txt`
- **By cancer type**: `Read tags/by-cancer-type/glioblastoma.txt`
- **By gene/drug**: `Grep "BRAF" corpus/by-pmid/`
- **Full text search**: `Grep "Warburg effect" corpus/by-pmid/`
- **Specific article**: `Read corpus/by-pmid/{PMID}.md`

## Methodology

1. **Literature Collection** — Systematic search of PubMed, Google Scholar, bioRxiv, clinical trial registries (ClinicalTrials.gov), and major oncology journals
2. **Categorization & Tagging** — Organize findings by cancer type, mechanism, stage of research (preclinical/clinical/approved), and efficacy data
3. **Gap Analysis** — Identify underexplored intersections (e.g., frequency therapy + immunotherapy combinations)
4. **Synthesis** — Develop original thesis on convergent curative strategies
5. **Article Writing** — Draft, review, revise to journal submission standards

## Conventions

- All claims must be traceable to a cited source
- Distinguish clearly between: established evidence, emerging evidence, and speculative/theoretical
- Use precise scientific terminology
- Maintain objectivity — present evidence for and against each approach
- Date all research notes so temporal context is preserved

## Important Notes

- This is a serious scientific endeavor aimed at real publication — rigor is paramount
- No unsupported health claims — everything must be evidence-based
- Prioritize recent research (2020-present) but include foundational/landmark studies where relevant
- Consider both efficacy AND safety/feasibility of proposed mechanisms
