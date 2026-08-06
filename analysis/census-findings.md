# What the census established (#CENSUS-FINDINGS)

Assembled by `scripts/census_findings.py` from the committed JSON of each
analysis, so it cannot drift from the measurements it summarises.

The question this answers: going from a 4,830-article keyword corpus to a
4,403,994-article MeSH census, what changed about what this project
believes?

---

## 1. The corpus was a 213-fold non-uniform sample

Measured with the same MeSH labels on both corpora, so descriptor scope
cancels. The frozen corpus captures
**41.9%** of the census's `mrna-vaccine` articles and
**0.20%** of its `epigenetic` -- a **213-fold** spread.

Any relative prevalence computed on the frozen corpus inherits that. It
is not a small uniform sample of the literature; it is a wildly uneven
one.

*What changed:* claims of the form "research is concentrated here" are
now measurable rather than assumed.
*Source:* `atlas-landscape.md`

## 2. The manuscript understated its own headline, and over-broadened another

**Volume.** Pharmacological to physical runs **9.1 : 1** by the
manuscript's keyword method and **17.6 : 1** on the census, because the
corpus over-samples physical modalities 3.3x -- exactly what a corpus
built from queries about them would do. The claim survives and was
understated by about half.

**Maturity.** "Physical modalities remain comparatively preclinical" does
not hold as a class. HIFU is **7.10%** clinical against CAR-T's
**6.64%**, both on precise descriptors, and HIFU and sonodynamic differ
by 1.6x. The defensible claim is narrower: sonodynamic therapy
specifically is early -- which is the mechanism this work rests on.

*What changed:* one claim strengthened, one narrowed. Both are in the
manuscript.
*Source:* `atlas-landscape.md`

## 3. The thesis sits on roughly thirty papers

The corpus contains no ferroptosis query and no PDT query, so it could
not measure this. The census can:

* ferroptosis-indexed cancer articles: **13,346**, growing 400 (2020) to 4,496 (2025)
* x drug resistance: **479**
* x photodynamic therapy: **177**
* x sonodynamic therapy: **32**

The resistance leg is supported by a literature. The sonodynamic leg,
the thesis's central mechanism, is supported by roughly thirty papers,
and the cautionary precedent (PDT) is better established than the thing
it is a precedent for.

*What changed:* the simulation work is carrying more of the argument
than the citation count suggested, and now says so.
*Source:* `atlas-thesis-position.md`

## 4. The two planning documents disagree about the keystone

P1 (persister/resistance) sits on **479**
ferroptosis articles; P4 (hypoxia), which `PREREGISTRATION.md` calls
the keystone, sits on **64**. The P1 protocol
calls P1 the highest-leverage prediction. Neither designation cited
evidence.

Not a quality ranking -- a sparse leg is where novelty lives. The
asymmetry that matters is that a negative P4 is ambiguous (mechanism
wrong, or experiment not yet worked out) while a negative P1 is simply
a negative result.

*What changed:* the choice can now be made deliberately.
*Source:* `atlas-prediction-position.md`

## 5. What the literature says to model next

| gene | articles | engine handle |
|---|---|---|
| HO-1 | 1,583 | none |
| P53 | 1,497 | none |
| TRANSFERRIN RECEPTOR | 1,482 | none |
| TNF-ALPHA | 1,013 | none |

None of them became a layer. When four of these were checked for
a calibration target (#616), the route that partially anchored ACSL4 --
cBioPortal within-cohort z-scores -- turned out to recover the normal
distribution at the z<-1 cut for every gene tested, so that cut carries
no gene-specific signal. At the deeper z<-2 cut TP53 does separate
(above the normal expectation in 31 of 32 cancer types), but it bounds
a prevalence rather than a dose-response, so it is recorded as a weak
anchor and no layer was written. Read the table as where the
literature's attention and the available data fail to overlap, not as
a backlog.

*Source:* `atlas-model-gaps.md`, `calibration-feasibility.md`

---

## What the census did NOT support

Reported here because a findings page that only lists wins is marketing.

**A repair this project made to its own co-mention layer made it worse.** #617 measured the layer at 55.5% precision (recomputed at 51% once every stratum was hand-judged rather than assumed), traced it to a filter that exempted multi-word forms, and replaced the single-token test with two measured filters. Re-measured on a fresh sample after the rebuild, precision FELL to 42%. The false positives had been multi-word because the multi-word channel was the unfiltered one; closing it moved the pressure to the channel just opened, and the top offenders are now bare English words. The replacement filters were then measured and do not separate true matches from false ones at all.
**It was then repaired for real, and the fix is measured at 88%** (#628). What separates true matches from false ones is not how much support a form has but whether the form is a NAME of the entity it resolves to, checked against NLM and NCBI rather than against the corpus. `treatment` is not a name of any descriptor; `xCT` is a name of SLC7A11. A blind panel of three judges who never saw the first verdicts put it at 88%, and the hostile bound, resolving every borderline case against the layer, is 80%.
So the standing finding is not that the layer is broken. It is that this project shipped a filter justified by an error distribution the filter itself changed, did not notice for two issues, and needed a fresh sample drawn after the rebuild to see it. The cost of the real fix is 28% of true matches, paid entirely on MeSH terms (35%) and not at all on genes (0%).
*Source:* `comention-regression.md`, `comention-authority-result.md`

**Literature-based discovery does not work as built.** The shipped ABC
ranking scores 10.1% precision@20 against
15.5% for ranking the same candidates by
popularity, and no standard link predictor beats that baseline either.
The candidate SET is genuinely informative, so that half stands. The obvious second half -- a bad RANKER -- does not follow, and two later measurements say why. Ordering the methods by how hub-selecting each one is reproduces the precision leaderboard exactly (rank correlation 1.00), so the metric rewards NOT correcting for candidate degree; and blending any seed-specific signal into a degree-only prior adds nothing measurable at any weight. On this metric a degree-correcting ranker and a bad one cannot be told apart.
*Source:* `atlas-discovery-eval.md`, `atlas-discovery-degree-bias.md`, `atlas-discovery-headroom.md`

**Replication looked like it was collapsing, and was not.** Scoring cohorts on whether they were EVER replicated gives 60.1% for 1950 to 13.4% for 2021; that is the observation window shrinking, not science changing, since the older cohort has had decades to acquire a second paper and the newer one had 5 years. On an equal 5-year window from each pair's own first assertion the decline is modest, and the recent end is an upper bound because of MeSH indexing lag.
*Source:* `atlas-replication.md`

**Most of what the census cannot corroborate is about what has been STUDIED, not about what is true.** Of the 20 simulation-module claims, 11 are corroborated by at least one asserting article and 9 by none. Read flat, that looks like 9 unsupported claims. It is not: a pair can only be asserted if BOTH its entities are written about, and the weaker entity's partner count across these claims runs from 6 to 2,792. Every claim that HAS support has a weaker entity of at least 108 partners, and 7 of the 9 zeros fall below that (Spearman rho = 0.86, and the association survives dropping every GPX4 pair).
**But which zeros are 'genuine' cannot be identified, and the source document deliberately names none.** The line is a sample minimum set by a one-article row; 45% of all asserted pairs in the graph sit below it, so it is not a detectability limit; and running the same procedure on the pair-level co-mention column inverts the correlation and returns a disjoint pair of exceptions. The finding is that a zero is a poor guide to a claim's truth -- not that any particular zero is interesting.
*Source:* `atlas-module-support.md`

**The entity collisions are not as bad as containment suggests.**
50.8% of relation rows touch a contested identifier, but only
1.35% rest on an uncorroborated one.
*Source:* `atlas-ambiguity-impact.md`

---

## Integrity problems found along the way

* **3 module citations pointed at unrelated papers** -- a Nature news item on fetal-tissue policy, a Theriogenology paper on embryo vitrification, and a PMID that does not resolve. Corrected. (`atlas-citation-audit.md`)
* **55.7% of the news pipeline's "verified" links** shared no content word with the claim they verified. Root cause: a claim yielding the single search term `Seven` matched 835,973 records and the five newest were accepted. The linker is fixed and re-run, and the same measurement now reads 1.9% -- but by WITHDRAWAL, not repair: 30 of the 44 verifications were dropped outright, and 19 of the 33 surviving pairs clear the bar on oncology boilerplate alone. (`news-verification-audit.md`)
* **`FSP1` resolved to a spastic-paraplegia gene**, leaving the manuscript's headline GPX4+FSP1 claim with zero typed relations. Blocklist now covers 346 measured sense collisions. (`atlas-ambiguity.md`)
* **The manuscript's mechanism count was 19 where the index carries 23** -- an undocumented 20-article threshold presented as coverage, and the four it hid were mostly physical modalities.

None of these were computational errors. Every one was a true statement describing something narrower than what it was used for.

---

## Every layer now carries a bound

* co-mention precision: **88% measured** on the layer as shipped (blind panel 88%, hostile bound 80%), superseding both the 60.1%-89.4% corroboration bound and the 42% the layer measured before the authority filter was turned on
* contradictions: ambiguity inflates the flag rate 1.45x
* emergence: 99.0% precision, 99.6% recall
* FSP1 disambiguation: 97.4%, with 75% of corrections extrapolated and that extrapolation independently tested

---

*Regenerate with `python scripts/census_findings.py`. Every figure is read
from the JSON of the analysis that produced it.*
