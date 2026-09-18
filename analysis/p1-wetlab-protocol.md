# P1 wet-lab protocol: GPX4i + FSP1i synergy in an FSP1-low persister line (#496)

This is the detailed, runnable protocol for **the keystone** falsifiable
prediction in the program, P1 (preregistered in `PREREGISTRATION.md`, experiment
brief E4). It is written so a wet-lab collaborator can execute it in isolation,
score it against a pre-stated numeric threshold, and report a clean confirm or
refute.

> **The keystone designation is `PREREGISTRATION.md`'s, not this document's.**
> This protocol used to call P1 "the highest-leverage falsifiable prediction in
> the program" on its own authority, while the preregistration named a different
> experiment. That disagreement is resolved (#619): the preregistration is the
> registering document and makes the call, on the basis that E4 is the most
> DECISIVE experiment rather than the most novel. E1, testing the hypoxia leg,
> remains the most novel test and is sequenced second. If that designation
> changes, it changes there and this file follows.

> **Status: this is a protocol, not a result.** No wet-lab work has been run by
> the project (single non-domain author plus AI). Executing it requires a
> collaborator with a tissue-culture lab. The point of writing it in full is to
> make the test cheap to adopt and impossible to re-interpret after the fact.

## Why P1 first

P1 depends on the **fewest contested assumptions** of the eight predictions. It
does not rely on the contested SDT oxygen-dependence (P4), the poorly characterized
RSL3 pKa (P7), spatial penetration physics (P2, P8), or the immune-coupling layer
(P5). The existence of **GPX4 and FSP1 as parallel ferroptosis-defense pathways** is
grounded in published biology (Doll 2019, PMID 31634900; Bersuker 2019, PMID
31634899; Hangauer 2017 persister GPX4-dependence, PMID 29088702). That foundation
motivates the experiment but does not guarantee synergy at every dose or in
every cell state. One dose matrix tests the interaction for the chosen exposures
and persister system; it does not establish or refute pathway independence on
its own.

## Model prediction being tested

- **Directional claim:** GPX4 inhibition + FSP1 inhibition is **synergistic** in
  FSP1-low persister-enriched cells.
- **Quantitative model output:** observed-to-Bliss **ratio about 1.99x**, not a
  Bliss excess. The committed 300-draw prior-predictive report gives a 95%
  interval of **[1.000, 5.242]** and a full sampled range of **[0.953, 7.800]**
  (`analysis/headline-uncertainty-report.md`). The rounded lower bound includes
  the additive null and some draws are sub-additive, so supra-additivity is not
  guaranteed. Improving on each single agent and exceeding Bliss independence
  are separate claims. These priors are not conditioned on experimental data.
- **Pre-stated falsification threshold:** Chou-Talalay combination index (CI)
  **greater than 0.8** at the matched-effect dose ratio,
  **or** measured combined kill at or below the Bliss-independence prediction within
  assay error.

These are the existing P1 decision thresholds; this clarification does not amend
`PREREGISTRATION.md`. Its original wording calls the 1.99x ratio an "excess" and
interprets the rounded lower interval bound as robust supra-additivity; the
artifact above does not support those interpretations. CI > 0.8 is the project's
chosen failure threshold, not the general definition of additivity: values
between 0.8 and 1 may still indicate weaker synergy under Chou-Talalay
([Chou 2010](https://pubmed.ncbi.nlm.nih.gov/20068163/)).

## Materials

- **Cell line.** A line with **low baseline FSP1 (AIFM2) expression** in which a
  drug-tolerant persister state can be enriched. Candidates: a persister-derived
  subline (e.g. EGFR-driven NSCLC persisters in the Hangauer 2017 system) or any
  line where FSP1-low status is confirmed by Western blot / qPCR before use.
  Confirm FSP1 status; do not assume it.
- **Persister enrichment.** Standard drug-tolerant-persister induction for the
  chosen line (e.g. sustained high-dose targeted therapy for ~9 to 14 days until a
  slow-cycling drug-tolerant population remains), per the line's established
  protocol.
- **Drugs.** GPX4 inhibitor: **RSL3** (or ML162, with its own single-agent response
  established). FSP1 inhibitor: **iFSP1**. Brequinar may be examined in a separate
  mechanistic comparison, not as an interchangeable selective FSP1 reagent (see
  target-attribution note below). Ferroptosis-pathway confirmation: **C11-BODIPY 581/591**
  (lipid peroxidation). Rescue control: **ferrostatin-1** (Fer-1) and
  **liproxstatin-1**.

**Target attribution.** Brequinar targets DHODH, the mitochondrial defense studied
by [Mao et al. (2021)](https://doi.org/10.1038/s41586-021-03539-7).
[Mishima et al. (2023)](https://doi.org/10.1038/s41586-023-06269-0) showed that
ferroptosis sensitization at high inhibitor concentrations can instead reflect
FSP1 inhibition; the [authors' reply](https://doi.org/10.1038/s41586-023-06270-7)
discusses the dependence on context. Interpret a brequinar comparison with
concentration-specific target engagement and DHODH/FSP1 genetic controls. A drug
response alone cannot distinguish those mechanisms, and a brequinar result
cannot replace the GPX4i-by-iFSP1 matrix specified for P1.

### How much precedent each arm has

The two arms are named symmetrically above and are not symmetric in the
literature, which is worth knowing before the bench rather than after.
Measured over the 13,346 census articles carrying the
`Ferroptosis` descriptor (`analysis/census-protocol-precedent.md`):

| reagent | articles naming it | role here |
|---|--:|---|
| RSL3 | 609 | arm 1, primary |
| ML162 | 12 | arm 1, offered surrogate |
| iFSP1 | 23 | arm 2, primary |
| brequinar | 12 | separate DHODH / off-target comparison |

The historical grouping pools 621 reagent mentions for arm 1 and
35 for its backup-defense comparison list — a factor of 17.7. The latter
combines iFSP1 and brequinar; it is not a count of selective FSP1-inhibitor
evidence or a claim that the two reagents are interchangeable. Three consequences
for whoever runs this:

1. **Dose-finding is asymmetric.** There is abundant published guidance for an
   RSL3 EC50 in a given line and almost none for iFSP1, so the single-agent
   pre-run matters far more on arm 2 than on arm 1. Budget for it.
2. **The offered surrogates are thinner than the primaries.** ML162 appears in
   roughly a fiftieth as many articles as RSL3, so a result obtained with it
   has correspondingly less to be compared against. Prefer RSL3 unless
   the chosen model system justifies the substitution, and say which was used.
3. **A negative result on arm 2 is harder to interpret.** With this little
   comparative work, a null could mean the hypothesis is wrong or that the
   dose or the FSP1-low status was wrong, and the protocol's pre-stated
   falsification threshold cannot distinguish them on its own.

None of this argues against the experiment. FSP1 inhibitors are recent — iFSP1
comes from the 2019 papers this protocol builds on — and a compound is rare in
the literature when it is new for the same reason it is rare when it is poor. A
count cannot tell those apart. What it can do is say where the thin ice is.

## Design

- **Dose matrix:** a full RSL3 (or ML162) by iFSP1 checkerboard, at least 6 by 6,
  spanning roughly 0.25x to 4x each single-agent EC50 (anchor the EC50s in a
  single-agent pre-run on the same persister-enriched cells).
- **Arms:** single-agent RSL3 series, single-agent iFSP1 series, the full
  combination matrix, vehicle, and a Fer-1 co-treatment of the most synergistic
  combination well (must rescue if the death is ferroptotic).
- **Replicates:** at least 3 biological replicates (independent persister
  inductions / passages), each with at least 3 technical replicates. Power the
  design to resolve a CI of 0.7 versus 0.9 (straddling P1's decision threshold)
  at the matrix center.

## Readouts

1. **Viability / death:** a live/dead or ATP-viability readout at a fixed endpoint
   (e.g. 24 to 48 h), used to compute the combination index, Bliss ratio and
   Bliss excess.
2. **Pathway confirmation:** C11-BODIPY lipid-peroxidation signal (flow or imaging)
   at the synergistic well, and **Fer-1 rescue** of that well (ferroptosis-specific
   death, not generic cytotoxicity).
3. **FSP1 status:** Western blot / qPCR confirming the line is FSP1-low at the time
   of the assay (the prediction is conditioned on it).

## Analysis

- Compute the **Chou-Talalay combination index** (CompuSyn or equivalent) at the
  matched-effect dose ratio. Separately report Bliss expected kill,
  `E = A + B - A*B` for fractional single-agent kills A and B, the **Bliss ratio**
  `observed / E`, and the **Bliss excess** `observed - E`. The ratio is undefined
  when E is zero; do not report an infinite synergy score.
- For the manuscript's illustrative rates (A = 0.40, B = 0.037, observed = 0.841),
  E = 0.4222, the ratio is approximately 1.99, and excess is approximately 0.419
  (41.9 percentage points). Chou-Talalay CI cannot be derived from this ratio;
  it needs the dose-response data.
- Pre-registered decision: **confirm** if CI is at or below 0.8 (synergy) with the
  combined effect above the Bliss-independence prediction; **refute** if CI is
  greater than 0.8 or the combined effect is at or below independence within assay
  error.
- Report the raw matrix, the CI surface, and the Fer-1 rescue, regardless of
  outcome.

## Expected result and what a refutation means

- **If confirmed:** the tested combination meets P1's synergy criterion in this
  cell state and dose region, supporting the combination case study. This does
  not independently establish biochemical pathway independence or clinical
  efficacy.
- **If refuted:** P1's specified synergy prediction fails for the tested system.
  The result does not by itself show that the two defenses are not independent.
  Check target engagement, response saturation, exposure timing and alternative
  defense pathways before assigning a mechanism. Report the outcome as
  prominently as a confirmation (the preregistration honesty clause).

## Cost and timeline

Low. A single dose-matrix viability experiment with a C11-BODIPY confirmation and a
Fer-1 rescue is standard 2D tissue culture, executable in a few weeks by one
person once the persister line is in hand. This is deliberately the cheapest
high-leverage test in the menu (`PREREGISTRATION.md` Part 2, E4).

## What this would change

This is the first test of a model prediction against biology the model was **not
built on**. A single confirmed-or-refuted result converts the headline numbers from
the model predicting its own behavior into a claim that has survived (or failed)
contact with data. That is the single biggest credibility step available to the
project (issue #496).
