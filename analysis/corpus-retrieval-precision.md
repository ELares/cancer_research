# Retrieval, not corpus size, was the bottleneck

Measured over **7,727,217 indexed titles and abstracts**, 5 mechanism questions, 10 known-relevant papers.

## The result

| | arm A: bag of words | arm B: phrase + title weight |
|---|--:|--:|
| MRR | 0.351 | 0.483 |
| recall@10 | 50.0% | 70.0% |
| gold papers reached (of 10, depth 50) | 5 | 4 |
| median pool searched | 6,638,629 | 71 |

**A trade, not a win, and the direction matters more than the averages.** Arm B ranks better (MRR 0.483 against 0.351, recall@10 70% against 50%) while reaching FEWER gold papers overall (4 against 5).

The per-question table below shows why, and it is not a tuning problem. Where the phrase matches the literature's own vocabulary, arm B finds the paper in a pool of tens of records against arm A's millions. Where it does not, arm B returns **nothing at all** -- its pool simply does not contain the answer. The failure mode is total rather than graceful, which is the opposite shape to arm A's (always finds something, usually buried).

The clearest case is `bacteria-immune-tagging`: arm B's pool of 289 records holds 0 of 5 gold papers, because the phrases ['engineered bacteria', 'tumor'] are not the words those papers use. So the concrete next step is vocabulary, not weighting.

## Per question

### bacteria-immune-tagging

> engineered bacteria that tag tumors for immune recognition

The mechanism the owner described: bacteria mark a tumour so white cells attack it.

| arm | pool | gold reached | best rank |
|---|--:|--:|--:|
| A | 6,196,053 | 1/5 | 23 |
| B | 289 | 0/5 | — |

Arm A's top 3:

- `40157514` Metabolic reprogramming of tumor microenviroment by engineered bacteri
- `42439475` Artificial Host-Guest Recognition Directs Glycometabolically Engineere
- `33425777` Human Short Peptidoglycan Recognition Protein PGLYRP1/Tag-7/PGRP-S Inh

Arm B's top 3:

- `38623680` Engineered Bacteria Detect Tumor DNA: Evaluation and Prospect.
- `37429881` Author Correction: Modular-designed engineered bacteria for precision 
- `38688906` Author Correction: Spatiotemporal control of engineered bacteria to ex

### ultrasound-gated-bacteria

> ultrasound controlled bacteria for cancer immunotherapy

Intersects this repo's own sonodynamic thesis.

| arm | pool | gold reached | best rank |
|---|--:|--:|--:|
| A | 5,877,824 | 1/1 | 6 |
| B | 29 | 1/1 | 6 |

Arm A's top 3:

- `35918309` Spatiotemporal control of engineered bacteria to express interferon-γ 
- `40688674` Sonogenetics regulation of intelligent engineered bacteria to reverse 
- `38688906` Author Correction: Spatiotemporal control of engineered bacteria to ex

Arm B's top 3:

- `38688906` Author Correction: Spatiotemporal control of engineered bacteria to ex
- `40104641` Ultrasound-responsive engineered bacteria mediated specific controlled
- `38640931` Ultrasound-visible engineered bacteria for tumor chemo-immunotherapy.

### bacteria-fenton-iron

> bacteria generating iron catalyzed fenton reaction in tumors

Fenton chemistry is iron-catalysed ROS, the ferroptosis mechanism this repo models.

| arm | pool | gold reached | best rank |
|---|--:|--:|--:|
| A | 6,771,944 | 1/1 | 22 |
| B | 72 | 1/1 | 4 |

Arm A's top 3:

- `35027953` Iron-doxorubicin prodrug loaded liposome nanogenerator programs multim
- `34547257` Photothermal nanozyme-ignited Fenton reaction-independent ferroptosis 
- `27855492` Sodium ascorbate kills Candida albicans in vitro via iron-catalyzed Fe

Arm B's top 3:

- `38098651` Long-term continuous degradation of carbon nanotubes by a bacteria-dri
- `34147772` Identifying driving factors of humic acid formation during rice straw 
- `38138725` A Synergistic Antibacterial Study of Copper-Doped Polydopamine on Ti3C

### microbiome-gpx4

> gut microbiome metabolite regulating GPX4 and ferroptosis

GPX4 is the gene the ferroptosis engine models; a microbial metabolite altering it is a causal bridge.

| arm | pool | gold reached | best rank |
|---|--:|--:|--:|
| A | 6,638,629 | 1/2 | 2 |
| B | 71 | 1/2 | 1 |

Arm A's top 3:

- `37752615` The gut metabolite 3-hydroxyphenylacetic acid rejuvenates spermatogeni
- `33779497` The gut microbiota metabolite capsiate promotes Gpx4 expression by act
- `42180552` Interactions between the gut microbiome and ferroptosis in degenerativ

Arm B's top 3:

- `33779497` The gut microbiota metabolite capsiate promotes Gpx4 expression by act
- `40876612` Ginsenoside Rb3 modulates gut microbiota to alleviate cerebral inflamm
- `40856206` Dietary hispidulin ameliorated chemically-induced colitis by inhibitin

### ferroptosis-macrophage-bacteria

> ferroptosis in macrophages defending against intracellular bacteria

Ferroptosis as an immune effector mechanism rather than a tumour-cell death mode.

| arm | pool | gold reached | best rank |
|---|--:|--:|--:|
| A | 6,654,762 | 1/1 | 1 |
| B | 16 | 1/1 | 1 |

Arm A's top 3:

- `35265210` Ferroptotic stress promotes macrophages against intracellular bacteria
- `25289573` Are inflammatory phagocytes responsible for resistance to facultative 
- `36211962` Heme oxygenase-1 modulates ferroptosis by fine-tuning levels of intrac

Arm B's top 3:

- `35265210` Ferroptotic stress promotes macrophages against intracellular bacteria
- `38059412` Endogenous chemicals guard health through inhibiting ferroptotic cell 
- `42039264` Targeted intestinal barrier repair via probiotic-derived engineered ou

## What this does and does not establish

Gold PMIDs were identified by field-restricted phrase search against the LIVE Europe PMC API (TITLE_ABS phrase queries, CITED desc), and each title was read to confirm the paper is genuinely about the mechanism named. They were NOT selected by running the local index, so neither arm of the evaluation chose its own answers. Membership in the corpus was verified independently against corpus/atlas/identity.sqlite: 10 of the 11 papers first identified are held; PMID 37358255 is not, and was dropped rather than counted as a retrieval failure.

**Known bias.** Because the gold was found by phrase search, every gold paper is guaranteed to contain its mechanism phrase. That favours arm B, which searches for phrases. The evaluation therefore measures whether phrasing plus field weighting can recover a paper the corpus demonstrably holds -- known-item retrieval -- and does NOT claim to measure precision over unjudged results.

Generated by `scripts/corpus_retrieval_eval.py`.
