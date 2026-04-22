# Book Outline

**Working title:** Mapping Cancer Therapy Mechanisms, Evidence Gaps, and Resistant-State Hypotheses

**Target:** ~50,000 words (current: ~16,000). 5 Parts, 11 chapters, 3 appendices.

---

## Part I: Why This Exists

*NEW content. Issue #101.*

### Chapter 1: The Problem (~3,500 words)

**Scope:** The human reality of cancer: survival statistics (sourced from WHO, NCI, IARC), why treatments fail for some patients, the residual disease problem, and what the current therapeutic landscape looks like. Written for non-specialists — no jargon without inline definition.

**Non-scope:** This chapter does NOT make therapeutic recommendations, does not introduce ferroptosis or simulation details, and does not advocate for any specific modality. It motivates the problem; it does not propose solutions.

**Source material:** README.md narrative voice, WHO cancer fact sheet (news-source-criteria.md example 1), NCI statistics, landmark papers on drug resistance and persister cells.

### Chapter 2: What We Set Out To Do (~2,500 words)

**Scope:** The open-science mission, why computational approaches complement experiments, what Monte Carlo simulation means in plain language, what corpus analysis is, and an honest framing of what we can and cannot contribute. Reader guide for different audiences (researchers, students, clinicians, engineers).

**Non-scope:** Does not present results or findings. Does not describe the simulation architecture in detail (that's Chapter 5).

**Source material:** README.md, CLAUDE.md guiding principles.

---

## Part II: What We Found

*MIGRATED from current v1.md sections 1-3 + 4.1-4.2. Expansion owned by issues #102, #104.*

### Chapter 3: The Literature Landscape (~5,000 words current, expandable to ~7,000)

**Scope:** Introduction to the cross-literature analysis, search strategy, corpus construction, categorization framework, analytical approach, limitations, corpus overview, all 19 mechanism categories (physical, immune, precision molecular, delivery), cross-mechanism convergence, the mechanism-cancer matrix, and interpretive framing of what the corpus supports.

**Non-scope:** Does not present simulation results (Part III). Does not present evidence-tier scoring (Chapter 4). Does not make clinical recommendations.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 1. Introduction | 3.1 |
| 2.1 Search Strategy | 3.2 |
| 2.2 Inclusion/Exclusion | 3.3 |
| 2.3 Categorization | 3.4 |
| 2.4 Analytical Approach | 3.5 |
| 2.5 Limitations | 3.6 |
| 3.1 Corpus Overview | 3.7 |
| 3.2 Physical Mechanisms | 3.8 |
| 3.3 Bioelectric | 3.9 |
| 3.4 Immune-Based | 3.10 |
| 3.5 Precision Molecular | 3.11 |
| 3.6 Delivery/Microenvironment | 3.12 |
| 3.7 Cross-Mechanism Convergence | 3.13 |
| 3.8 Mechanism-Cancer Matrix | 3.14 |
| 4.1 Main Contribution | 3.15 |
| 4.2 Corpus Support for PDT/SDT/CAP | 3.16 |

**Source material:** corpus/INDEX.jsonl, analysis/convergence-map.md, analysis/mechanism-matrix.md, analysis/taxonomy-rerun-notes.md.

### Chapter 4: Evidence Quality and Gaps (~3,000 words current, expandable to ~5,000)

**Scope:** Evidence maturity scoring, weighted evidence methodology, gold-set evaluation (46% exact, 96% precision, 55% recall), tissue-of-origin layer, diagnostic-therapy matching, evidence coverage caveats, and known corpus gaps.

**Non-scope:** Does not discuss simulation findings. Does not make claims about which therapies are "best" — only which have deeper or shallower evidence bases.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 3.9 Evidence Maturity and Weighted Scoring | 4.1 |

**Source material:** analysis/evidence-coverage-audit.md, analysis/evidence-gold-eval.md, analysis/tissue-evidence-summary.md, analysis/diagnostic-therapy-audit.md, analysis/landmark-corpus-gaps.md, analysis/weighted-evidence-summary.md.

---

## Part III: What The Simulations Show

*MIGRATED from current v1.md sections 4.3-4.6. Expansion owned by issue #102.*

### Chapter 5: The Ferroptosis Engine (~3,000 words)

**Scope:** How the ferroptosis simulation works, explained for non-programmers. What each parameter represents biologically. The cell model (GSH, GPX4, FSP1, lipid peroxidation, death threshold). Why Monte Carlo. Key calibration decisions and their sources. What the model CAN and CANNOT represent.

**Non-scope:** Does not present TME results (Chapter 6) or drug combination results (Chapter 7). Does not claim the model is a faithful replica of biology — it is a hypothesis-exploration tool.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 4.3 Computational Simulation (overview, lines 309-341) | 5.1 |

**Source material:** simulations/ferroptosis-core/src/biochem.rs, simulations/ferroptosis-core/src/params.rs, analysis/hypothesis-sdt-ferroptosis-icd.md.

### Chapter 6: Drug Combinations and Penetration (~4,000 words)

**Scope:** Spatial tumor simulation (depth-dependent kill), vulnerability window timing, drug combination synergy (Bliss independence), and the in-vitro-to-in-vivo gap from tissue-specific pharmacokinetics.

**Non-scope:** Does not claim SDT is clinically superior. Presents the simulation findings and their implications, but frames the modality question as open.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 4.3.1 Spatial Tumor Simulation | 6.1 |
| 4.3.2 Vulnerability Window | 6.2 |
| 4.3.3 Combination Modeling | 6.3 |

**Source material:** simulations/sim-spatial, sim-vuln-window, sim-combo-mech, sim-tumor-pk output.

### Chapter 7: Resistance Mechanisms (~6,000 words)

**Scope:** Four TME resistance/amplification mechanisms, each with: plain-language opening, simulation finding, "what this means," "how to test this," and honest limitations. The key insight: three resistance mechanisms (hypoxia, stromal shielding, acidic pH) selectively penalize pharmacologic ferroptosis inducers while leaving physical modalities unaffected. One amplification mechanism (immune coupling) favors dense spatial kill patterns.

**Non-scope:** Does not overstate simulation confidence. Each finding is framed as "the model predicts X, subject to [listed assumptions]." Does not propose clinical protocols.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 4.3.4 Oxygen Gradients and Hypoxia | 7.1 |
| 4.3.5 Immune Coupling | 7.2 |
| 4.3.6 Stromal Shielding | 7.3 |
| 4.3.7 Acidic pH | 7.4 |

**Cross-reference update:** Internal references "Section 4.3.4/4.3.5/4.3.6" become "Section 7.1/7.2/7.3."

**Source material:** simulations/sim-tme output, analysis/principle-resistance-tradeoff.md.

### Chapter 8: Counterarguments and Failure Modes (~3,500 words)

**Scope:** Honest assessment of what could go wrong. Other directions the repo was underweighting (radioligand therapy, immune strategies, pathway-target layers, stromal biology). Historical precedents (PDT's 40-year translational plateau). The nanosonosensitizer confound. Evolutionary escape. Where the simulation architecture is structurally too simple.

**Non-scope:** This chapter must not soften counterarguments with "but actually our approach handles this." If the counterargument is valid, say so. The chapter's job is intellectual honesty, not defense.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 4.4 Other Directions | 8.1 |
| 4.5 The Modality Question | 8.2 |
| 4.6 Counterarguments, Precedents | 8.3 |

**Source material:** analysis/gap-analysis.md, analysis/pathway-target-audit.md.

---

## Part IV: What Should Happen Next

*PARTIALLY NEW. Issues #103 (Ch 9-10), #105 (news integration sidebars).*

### Chapter 9: Research Roadmap (~6,000 words)

**Scope:** Current Section 4.7 expanded into experiment-level detail. Each proposed experiment includes: hypothesis, model system, expected outcome if right, expected outcome if wrong, confidence tier (high/medium/low based on how well-calibrated the underlying simulation parameters are). Ordered by confidence and feasibility.

**Non-scope:** Does not claim to be a grant application or a clinical protocol. Experimental suggestions are hypothesis-generating guidance, not validated assay designs. Timelines and costs are rough illustrations, not commitments.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 4.7 Research Priorities | 9.1 |

**Source material:** simulation outputs, analysis/distilled-hypotheses-final.md.

### Chapter 10: The Broader Landscape (~3,500 words)

**Scope:** Current Section 5 (Conclusion) expanded with field context. Where the current therapeutic landscape is heading (immunotherapy plateau, CAR-T solid tumor challenge, radioligand rise, AI in drug discovery). How this work fits. What's missing from everyone's approach. Competing hypotheses worth taking seriously.

**Non-scope:** Does not make funding recommendations or predict winners. Presents the landscape honestly.

**Migrated content:**
| Current section | New location |
|----------------|-------------|
| 5. Conclusion | 10.1 |

**Source material:** analysis/deep-pattern-analysis.md, news sources (after #99 pipeline).

### Chapter 11: How To Contribute (~2,000 words)

**Scope:** Practical guide for researchers (use ferroptosis-core, challenge parameters), engineers (extend simulations, Python bindings, PhysiCell), students (thesis project ideas), and clinicians (what to watch for in upcoming trials). Links to code, data, and the MIT license.

**Non-scope:** Does not repeat technical content from earlier chapters. Points to them.

**Source material:** README.md, ferroptosis-core documentation, ferroptosis-python, ferroptosis-ffi.

---

## Part V: References and Tools

*NEW content. Issue #104.*

### Appendix A: Parameter Documentation (~3,000 words)

Every simulation parameter: name, default value, biological meaning, literature source or "estimated," confidence level, sensitivity at 2x/0.5x.

### Appendix B: Reproduction Guide (~1,500 words)

Step-by-step: clone, build, run analysis, run simulations, compile PDF.

### Appendix C: Glossary (~3,500 words)

80-120 terms, plain-language definitions, cross-referenced to chapters where used.

---

## Word Budget Summary

| Part | Chapters | Current words | Target words |
|------|----------|--------------|-------------|
| I: Why This Exists | 1-2 | 0 | 6,000 |
| II: What We Found | 3-4 | ~8,000 | 12,000 |
| III: Simulations | 5-8 | ~7,000 | 16,500 |
| IV: What's Next | 9-11 | ~1,000 | 11,500 |
| V: References/Tools | Apps A-C | 0 | 8,000 |
| **Total** | **11 + 3 apps** | **~16,000** | **~54,000** |

## Chapter → Issue Mapping

| Chapter | Owner Issue | Depends On |
|---------|-----------|------------|
| 1-2 | #101 | #100 (this outline) |
| 3-4 | #102, #104 | #100 |
| 5-8 | #102 | #100 |
| 9-10 | #103 | #102 (findings must exist to reference) |
| 11 | #103 | #100 |
| Apps A-C | #104 | #102 (content must exist to audit) |
| News sidebars | #105 | #99, #101, #102 |
