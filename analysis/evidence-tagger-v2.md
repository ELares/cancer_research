# Evidence tagger v2: how the accuracy gain was obtained (#TAGGER-V2)

Off-by-default layer (`FERRO_EVIDENCE_V2=1`) that roughly halves the evidence-tier
error rate. The frozen corpus, `corpus/INDEX.jsonl`, and every manuscript number
are unchanged with the flag off, guarded by
`tests/test_evidence_v2.py::test_flag_off_is_byte_identical_to_frozen_corpus`.

Regenerate every number below with:

```bash
python scripts/evaluate_evidence_v2.py     # writes analysis/evidence-v2-eval.md
```

## Headline

| evaluation set | n | v1 exact | v2 exact | error reduction |
|---|---|---|---|---|
| HELDOUT (never inspected during development) | 170 | 51.2% | **80.6%** | **2.52x** |
| CONSENSUS (both annotators agree) | 77 | 51.9% | **75.3%** | **1.95x** |
| DEV (v1 human labels, tuned on) | 100 | 46.0% | 63.0% | 1.46x |

Binary evidence-detection F1 rises from 0.705 to 0.966 on the held-out set,
with precision holding at 96.0% (v1: 100%) and recall rising from 54.4% to 97.3%.

The conservative claim is the CONSENSUS row: **error halved, 1.95x**. The
HELDOUT row is the generalization test but shares a methodological bias with its
labels (see caveats).

## What was wrong

The production tagger reads title + MeSH + PubTator entity strings + abstract.
`include_full_text` defaults to `False`, so 216 MB of stored full text was never
read, and binary recall sat at 55%.

Naively enabling full text lifts recall to 96.6% but drops precision from 96.0%
to 89.4%, and the damage is one-directional: preclinical articles get promoted
into clinical tiers. Measured on the v1 gold set, enabling whole full text fixed
20 predictions and broke 14. The cause is that Introduction and Discussion
sections describe **other people's** studies. PMID 40700574 is a mouse study
whose Discussion contains "The most recent mRNA vaccine to undergo a phase 3
clinical trial ..." — that one sentence promoted it to `phase3-clinical`.

The project's own labeling guideline already states the correct rule: take the
tier from the article's primary research, "not in cited references or background
discussion". v2 implements that rule mechanically.

## The five changes, each traced to a measured error cluster

1. **Section-scoped full text** (`scripts/evidence_sections.py`). The
   `## Full Text` blob is split into sections and each is classified SELF
   (Methods, Results, Trial design, ...), CITED (Introduction, Discussion,
   Conclusion, ...) or DROP (References, Funding, Ethics, ...). Only SELF is
   read. Headings survive as bare blank-line-surrounded lines; this shape finds
   headings in 98.7% of a 600-article random sample, and unrecognised headings
   inherit the enclosing class so Methods subsections ("Calcein assay") stay SELF
   without enumeration. Median SELF fraction of full text: 46%.

2. **Expanded `theoretical` vocabulary.** `theoretical` scored 5% recall (1 of
   20) with **zero** keyword hits on 12 of 12 inspected misses. The five shipped
   terms describe mechanistic modelling; the dominant genre in this corpus is
   public-cohort bioinformatics — "Identification of Tumor Antigens and Immune
   Subtypes in Lung Adenocarcinoma", "Multi-seed searching algorithm for codon
   optimization". None say "simulation" or "in silico". Added TCGA/GEO mining,
   signature construction, docking, ML, and device-dosimetry vocabulary.
   Recall 5% -> 35% (dev), 0% -> 54.5% (held-out).

3. **In-vivo reagent guard** (`INVIVO_REAGENT_ANTIPATTERNS`). Bare `murine` and
   `animal model` match reagent and cell-line prose. Measured: PMID 39557959
   (in-vitro) was promoted by "murine lewis lung carcinoma **cell line** ll/2".
   An in-vivo keyword hit is discounted when every occurrence sits inside a
   reagent context window.

4. **Opinion pub-type veto.** Editorials, comments and letters were the largest
   single `none-applicable` leak once full text is read — an editorial
   discussing a phase III trial inherited that trial's tier. #346 already
   established this veto for the MeSH branch ("a commentary is never primary
   evidence"); v2 applies it to the whole decision. This single change moved the
   held-out error reduction from 2.24x to 2.52x.

5. **Theoretical-dominance rule.** v1 evaluates tiers in a fixed order with
   `theoretical` last and first-match-wins, so a bioinformatics paper that merely
   names a cell line is tagged `preclinical-invitro`. v2 lets `theoretical` win
   when the computational signal is strong (>= 3 distinct markers) and the
   wet-lab signal is incidental (no in-vivo hit, <= 1 in-vitro hit). This does
   not override the guideline rule that real wet-lab work outranks in-silico
   work.

### Changes that were tried and rejected

Six `clinical-other` phrases were added and then removed after measurement,
because they fired more on preclinical than clinical records:

| phrase | false | true |
|---|---|---|
| `we treated` | 7 | 1 |
| `institutional review board` | 7 | 5 |
| `informed consent was obtained` | 5 | 1 |
| `were recruited` | 2 | 1 |

`we treated` matches "we treated cells with ..."; animal protocols and
tissue-donor consent statements carry IRB and consent language. Keeping them
cost 13 records to `preclinical-invivo -> clinical-other` confusion.

A "no Methods section implies review" rule was also tried and rejected: measured
on the v1 gold set it is only 26% precise as a lone `none-applicable` predictor
(54% of none-applicable records lack a Methods section, but so do 25-35% of
genuine primary-research records). `evidence_sections.has_methods_section()` is
retained as a diagnostic with that caveat documented.

## Caveats that bound these numbers

**Annotator agreement is 77%.** The v1 human labels and the independent
full-text relabel agree on 77 of 100 records. Disagreement concentrates on
`theoretical` and the in-vitro/in-vivo boundary — the same classes the tagger
finds hardest. Neither column is ground truth, and ~77% is the practical
ceiling for this task. A meaningful fraction of what looks like tagger error is
annotator disagreement: PMID 41554743 is labeled `preclinical-invitro` by v1 but
its Methods describe orthotopic tumor-bearing mice.

**The held-out labels are not human labels.** They come from an LLM relabeling
pass that read Methods and Results and applied the committed guidelines, blocked
from reading any gold file or the tagger's own `evidence_level` field. They are
independent of both the human labels and the tagger, but they are machine
labels and `analysis/evidence-gold-set-v3-fulltext.csv` says so in its header.

**Shared methodological bias.** The held-out labels and the v2 tagger both read
full text under the same "ignore Introduction/Discussion" rule, so the 2.52x
held-out figure is an upper estimate. The 1.95x consensus figure — where a human
annotator who never saw that rule independently agreed — is the defensible one.

**Small strata.** `phase2-clinical` and `phase3-clinical` have 0-3 records per
evaluation set, so their per-tier recall is not measurable here. The
Wilson intervals in `analysis/evidence-v2-eval.md` are wide for every tier.

## What is still wrong

Largest remaining held-out error clusters under v2:

- 10 x `preclinical-invitro` -> `preclinical-invivo`. Partly genuine annotator
  disagreement (papers that do both), partly residual reagent leakage.
- 5 x `none-applicable` -> `preclinical-invivo`. Narrative reviews with no
  review pub-type and no review-shaped title.
- `phase2` / `phase3` remain driven almost entirely by NLM publication types;
  the text path for them is unmeasured at this sample size.

The honest next step is a larger, dual-annotated gold set with a reported
agreement statistic, not further keyword tuning against 100 records.

## Promotion checklist

This layer is deliberately NOT applied to the production corpus. Promoting it
would re-tag `corpus/INDEX.jsonl` and move published numbers, so it requires:

1. a decision to re-freeze the corpus and re-run every dependent analysis;
2. updating the manuscript's 46% / 96% / 55% figures and their caveats;
3. re-running `scripts/oa_bias_analysis.py` and the evidence-tier figures;
4. recording the re-tag in `analysis/provenance.jsonl`.
