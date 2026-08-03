# Mission

## The old mission

> Analyse thousands of cancer research articles and run biochemical simulations
> to propose new approaches for targeting drug-resistant tumour cells.

That mission has been substantially executed. It produced a 4,830-article
corpus, a 46,700-word manuscript, an 11-binary simulation suite, and a
calibration ledger. It also hit its ceiling, and the ceiling is measurable:

**4,830 articles is 0.11% of the cancer literature.** PubMed indexes 4,276,707
articles under the MeSH tree C04 (Neoplasms). The corpus was assembled from 19
hand-written keyword queries with a 500-record-per-query cap and inconsistent
date windows, and it contains no photodynamic-therapy query and no ferroptosis
query at all — the two topics the project's own simulation half is built on.

So every statement of the form "research is concentrated here" or "this is a
gap" has been, in truth, a statement about a retrieval design. That is the
honest reading of the last two years of work, and it is the reason for a new
mission rather than more of the old one.

## The new mission

**Map the entire recorded history of cancer research into a single normalized,
queryable evidence base — and mine it for what no one has had the vantage point
to see.**

Not a sample. The census. Every article humanity has indexed on cancer, with its
entities resolved to stable identifiers, its claims attributable to a sentence,
and its evidence tier measured rather than assumed. Then: find the patterns,
the contradictions, the replications that never happened, and the connections
that exist across literatures nobody reads together.

The simulation engine stops being the centrepiece and becomes the instrument:
where the literature implies a mechanism, the engine tests whether the mechanism
can carry the weight put on it.

## What that means concretely

| | Now | Target |
|---|---|---|
| Articles | 4,830 | 4,276,707 |
| Selection | 19 keyword queries, 500-record cap | MeSH tree C04, the complete census |
| Full text | 4,830 (98.7% open-access, a biased slice) | ~919,000 open-access, with the bias measured |
| Entities | unnormalized surface strings (3,260 spellings of ~1,200 genes) | NCBI Gene, MeSH, ChEBI, Cellosaurus identifiers |
| Relations | none | typed, scored, sentence-attributed |
| Evidence tiers | 57.8% of records untagged | measured, with a reported annotator ceiling |

## What this cannot do, stated plainly

This section is not modesty. It is the part that makes the rest credible.

**We cannot obtain all full text.** Roughly 21% of cancer articles are in PMC
open access. The remainder is paywalled, and neither scraping it nor
redistributing it is legal or intended. Full-text analyses therefore run on an
open-access subset whose bias must be measured and reported every time, not
waved at. Abstracts and MeSH indexing, by contrast, we can have completely.

**A literature census does not cure anything.** It finds what has been missed,
what has been contradicted, what was never replicated, and what has never been
connected. Those are hypotheses. Every one of them still has to survive a
laboratory, and this project has no laboratory.

**Scale is not rigour.** Four million badly-tagged articles are worse than five
thousand well-tagged ones, because the errors become invisible. Every layer in
the atlas ships with a measured error rate or it does not ship. The evidence
tagger rebuild is the template: it reports 2.59x error reduction *and* the 77%
annotator-agreement ceiling that bounds it.

**"Solving cancer" is not a thing one repository does.** The honest ambition is
to build the map that other people — with labs, with patients, with funding —
can navigate, and to give it away. That is what the mission means. Anything
stronger is marketing, and this project's only real asset is that it has never
been marketing.

## How the work is ordered

1. **Acquire the census.** *Done, and it holds.* The first pass took 4,203,236
   cancer articles from 39,994,988 scanned — 98.3% of what PubMed itself reports
   for `neoplasms[mh]`, so the local MeSH filter reproduces NLM's tree
   explosion. Alongside it: **783,271** articles recovered that MeSH has not
   indexed yet (the un-indexed share climbs to 37.6% in the most recent baseline
   files), **520,143** open-access full texts on external storage, and
   **7,951,325** typed relations over 1,603,105 PMIDs.

   Two corrections the census forced on itself, both instructive:

   * A tree-C04 definition **misses foundational mechanism papers**. Both
     founding FSP1 papers (Doll 2019, Bersuker 2019, *Nature*) are absent from
     it, and 62% of all ferroptosis papers sit outside C04. The definition now
     includes nine adjacent experimental-context descriptors, with the C04 core
     kept separable.
   * Neither a keyword corpus nor a MeSH census is a superset of the other. The
     old 4,830-article corpus contains ~223 cancer papers the census misses, and
     47 that are not cancer at all.
2. **Normalize.** Identifiers, not strings. A gene is `NCBIGene:2879`, not
   twelve spellings of GPX4.
3. **Measure every layer.** Precision and recall with intervals, against labels
   whose own agreement is reported.
4. **Mine.** Contradiction detection, replication tracking, temporal emergence,
   and literature-based discovery across disconnected fields.
5. **Test.** Where the map implies a mechanism, the simulation engine tries to
   break it.
6. **Give it away.** Every layer citable, reproducible, and MIT-licensed.

## What stays the same

The tenets in `CLAUDE.md` do not change, and the first one governs this
document: let the evidence lead. If the census says the project's own thesis was
an artifact of its search design, that is the finding, and it gets published
like any other.
