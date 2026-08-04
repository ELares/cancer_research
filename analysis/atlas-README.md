# The atlas: start here

The atlas is a census of the cancer literature and a set of layers built on top
of it. This page orients someone arriving cold: what exists, what each layer can
and cannot answer, and which caveats are load-bearing rather than decorative.

The short version: **the atlas is a map, and the map has measured holes.** Every
number below ships with the error rate that bounds it, because a census that
hides its error rate is worse than a small corpus that reports one.

## The three corpus surfaces

Do not confuse them. They answer different questions and only one is frozen.

| surface | what it is | mutability |
|---|---|---|
| `corpus/by-pmid/` + `corpus/INDEX.jsonl` | the 4,830-article **frozen** snapshot every manuscript number is computed from | immutable |
| `corpus/living/` | NEW records found after the freeze | append-only, dated |
| `corpus/atlas/` | the **whole cancer literature**, re-derived | rebuilt in place |

Anything that must reproduce a published number reads the frozen index. Anything
that wants coverage reads the atlas. They are not supersets of each other, which
is itself a finding — see *Coverage* below.

## What is in the atlas

| layer | size | built by |
|---|---|---|
| MeSH-indexed cancer articles | **4,403,994** (4,203,236 C04 + 200,758 adjacent) | `atlas_baseline.py` |
| recovered, not yet MeSH-indexed | 783,271 | `atlas_unindexed.py` |
| open-access full texts (on external storage) | 1,100,218 | `atlas_fulltext.py` |
| typed, normalized relations | 7,951,325 over 1,603,105 PMIDs | `atlas_relations.py` |
| queryable relation index | 2,186,309 entity pairs | `atlas_graph.py` |
| full-text sentence co-mentions | 6.2M+ pairs (build in progress) | `atlas_comention.py` |

## The layers, and what each can actually answer

### Census — `atlas_baseline.py`

Ingests the PubMed annual baseline (1,334 gzipped XML files, ~40M records).
E-utilities **cannot** do this: `esearch` over `neoplasms[mh]` reports all
4,276,707 hits but paging caps at `retstart` 9,999 even with WebEnv history.

*Cancer* means MeSH tree **C04** (704 descriptors, from NLM's own SPARQL
endpoint) **union** nine adjacent experimental-context descriptors.

> **Load-bearing caveat.** The adjacent extension exists because a C04-only
> definition misses foundational mechanism papers. Both founding FSP1 papers
> (Doll 2019, Bersuker 2019, *Nature*) are absent from C04 — their only
> tumour-related descriptors are `Cell Line, Tumor` (tree A11) and `Gene
> Expression Regulation, Neoplastic` (G05). **62% of all ferroptosis papers sit
> outside C04.** Records matched only by the extension carry
> `cancer_basis: "adjacent"`, so the C04 core stays separable.

### Recovery — `atlas_unindexed.py`

MeSH indexing lags publication, and the lag is not uniform: the un-indexed share
trends sharply upward with recency -- 0.0% in 1975-era baseline files, 5.2% by
file 409, 18.4% by file 817, 37.6% by file 959 -- though not monotonically, since
individual files fluctuate. A pure MeSH census therefore loses
disproportionately much of the recent literature.

> **Load-bearing caveat.** These are **text-matched**, not MeSH-indexed, and live
> in a separate stream tagged `source: "text-match"`. The matcher's accuracy is
> measured against MeSH truth on indexed articles: **precision 75.7%, recall
> 95.6%**. So roughly one in four recovered records is not really cancer. Quote
> that number wherever a count from this layer is used.

### The full-text recency ceiling

Full text is kept only for articles the census already knows, so the PMC bulk
cannot be matched past whatever the PubMed baseline contains. On the 2026-06-17
bulk both `PMC013xxxxxx` packages returned **exactly zero** cancer articles from
232,890, while every other package yielded 14–18%.

> **A cliff, not a slope.** MeSH indexing lag produces a gradual decline; an
> abrupt 17.7% → 14.7% → 0.0% does not. The census's PMC identifier space stops
> at `PMC128xxxx`, so nothing in the `PMC13` block can match at all — the PMC
> full-text release is simply newer than the PubMed baseline behind the census.
>
> Cost: an estimated **32,000–41,000** cancer full texts (232,890 articles at the
> 13.9–17.7% interquartile yield of the reachable packages). Closing it needs a
> newer PubMed baseline, not a code change. It compounds with the MeSH lag
> already documented under *Recovery*: both bite hardest at the recent end.

### Coverage — `atlas_coverage.py` → `atlas-coverage.md`

Supplies the denominator every manuscript ratio was missing. The frozen corpus
holds **0.086%** of the census. **20.9%** of the census has a PMC identifier, which upper-bounds any
full-text claim (a PMC id is necessary but not sufficient for the text to be in
the redistributable open-access subset -- the actual pull retrieved fewer) — against the frozen corpus's 98.7%
open-access, which is not what the literature looks like.

It also measures the frozen corpus against the census, and the result cuts both
ways: 940 frozen records are simply not yet indexed, 270 are indexed under
non-C04 trees (of which ~223 *are* cancer papers), and 47 are genuinely not
cancer. **45 of those 47 carry no mechanism tag**, so that residue is a
*query-level* defect in `scripts/queries.txt`, distinct from the *tag-level*
precision problem measured separately.

### Relations and the graph — `atlas_relations.py`, `atlas_graph.py`

NCBI's PubTator3 bulk release: entities resolved to NCBI Gene and MeSH
identifiers, and relations over a **fixed** predicate vocabulary (`associate`,
`treat`, `cause`, `inhibit`, `stimulate`, `positive_correlate`,
`negative_correlate`, `cotreat`, `interact`, `compare`, `prevent`,
`drug_interact`).

> **Three load-bearing caveats.**
>
> 1. **There is no biological-process entity type.** Types are Gene, Chemical,
>    Disease, Species, CellLine, Variant. `ferroptosis` is not addressable, so
>    "X induces ferroptosis" **cannot be expressed**. Every mechanistic claim
>    collapses to a gene-gene or gene-chemical pair, which is a real loss.
> 2. **About half of all relations are `associate`**, which is nearer
>    co-mention than knowledge. Always read the per-predicate breakdown.
> 3. **Query by identifier, not symbol.** See the entity audit.

### Co-mention — `atlas_comention.py`

The recall complement to the relation graph, and the reason literature-based
discovery is worth revisiting. PubTator's relations are extracted from
ABSTRACTS, and their edge recall is far too low for anything that reasons from
absence: GPX4 and caspase-3 share 236 PubMed abstracts and **no graph edge**.

This applies PubTator's own normalized vocabulary to the open-access full texts
already on disk, extracting entity co-mention at SENTENCE level. Deterministic,
offline, no LLM. Measured on 6 of 28 shards, so these are lower bounds:

| pair | PubTator edge | full-text co-mention |
|---|---|---|
| GPX4–caspase-3 | 0 | 63 |
| AIFM2–GPX4 | 0 | 69 |
| GPX4–SLC7A11 | 31 | 1,725 |
| ACSL4–GPX4 | 3 | 756 |

> **Load-bearing caveat.** Co-mention carries **no predicate, no direction and no
> polarity**. "X does not inhibit Y" produces the same edge as "X inhibits Y".
> It cannot replace the typed relations; it answers only "does the literature
> discuss this pair at all", which is precisely the question a zero in the
> relation column raises. The alias filter is the other weak point: 87.7% of
> surface forms survive it, and a spurious one poisons every downstream count.
>
> **That filter was about shape, not sense.** Its length rule happens to exclude
> `psa`, `p21`, `p62` and `er`, but 75 measured sense collisions passed it —
> including `cox-2` and `fsp1`, so full-text COX-2 co-mentions were counted
> against mitochondrial cytochrome c oxidase and FSP1 against a
> spastic-paraplegia gene. A sense filter now runs after the shape filter: one
> form is redirected to its measured cancer-domain sense, the other 74 dropped.
> Corrected counts arrive with the next co-mention build.

### Entity audit — `atlas_entity_audit.py` → `atlas-entity-audit.md`

Checks the symbols this project queries against NCBI's own record. **1 mismatch
in 53**, and it is the worst possible one: PubTator3 maps 357 mentions of `FSP1`
and 90 of "ferroptosis suppressor protein 1" to gene **51062 = ATL1 (atlastin
GTPase 1)**, a hereditary spastic paraplegia gene. The real FSP1 is **84883
(AIFM2)**.

That collision sits directly under the manuscript's headline GPX4+FSP1
Bliss-synergy claim. Querying the graph for `FSP1` silently returned ATL1's edges.

> This audit **reports** the collision; it does not fix it, and for a while
> nothing did. The two layers below now do, and they replaced the earlier
> position that reproducing the error was sufficient.

### Entity ambiguity — `atlas_ambiguity.py` → `atlas-ambiguity.md`

Asks the question the audit did not: **how many other FSP1s are there?**
Measured across all three entity types, genes are the outlier by roughly ten
times — 27.7% of gene mentions sit on a contested surface form against 2.0% for
chemicals and 2.7% for diseases, because MeSH is a curated vocabulary with one
preferred term per concept while NCBI lists `FSP1` as an official alias of three
different genes.

That contested share is **not** an error rate. Contested forms split three ways:
species ambiguity (human vs mouse *GAPDH*, benign), hierarchical granularity
(*Glioblastoma* under *Glioma*, lossy but not wrong) and genuine sense collisions
(*EREG* vs *ESR1* for `ER`, damaging). Only the last is an error, and pooling
them would overstate the damage several-fold.

The damaging class reaches the field's most-queried concepts: `ER` resolves to
**epiregulin** rather than the estrogen receptor across 54,293 mentions, `COX-2`
to mitochondrial cytochrome c oxidase rather than PTGS2, `p21` to a
**pseudogene**, and `PC` splits across **prostatic** and **pancreatic**
neoplasms. `atlas_graph.resolve` now returns `None` on these instead of a
plausible-looking wrong entity, so an analysis fails loudly.

### Sense disambiguation — `atlas_disambiguate.py` → `atlas-disambiguation.md`

`FSP1` gets no blanket fix, because there is nothing to fix it *to*: it is an
official alias of AIFM2, S100A4 **and** ATL1, and remapping everything to AIFM2
would corrupt the S100A4 papers, which are the plurality. So the sense is decided
per paper.

Against a gold set of the 242 papers that declare their own expansion — a label
independent of anything the classifier reads — PubTator3 is **36.3%** accurate
and this layer is **96.2%**. The number is only meaningful because every
label-defining phrase is masked out before a feature is read; without that the
measurement would be circular.

**Most corrections are extrapolated, and that is tested.** Accuracy is measured on
papers that declare a sense, but 75% of the 1,191 corrections land on papers that
declare nothing. Publication year checks that population independently: of 175
undeclaring papers corrected to AIFM2, **zero** predate 2019 — against roughly 37
expected if the classifier were assigning it without regard to the biology.

The payoff lands on this project's own headline. GPX4+AIFM2 co-mention papers go
**98 → 257**, and typed GPX4↔AIFM2 relations go **0 → 15** — there were none at
all, because every one had been filed under atlastin. The corrections are
applied when the index is built, so `atlas_module_support` moved from **9 of 20**
corroborated module claims to **10 of 20**.

### Are the curated senses right? — `atlas_domain_sense.py`

`DOMAIN_SENSE` says which sense the literature means for five high-volume
collisions. Those entries were written from domain knowledge, so this checks them
against the corpus by counting papers that declare a sense in their own words.

All five hold: `ER`→ESR1 at 98.7% of 451 declaring papers (epiregulin declared
**zero** times), `COX-2`→PTGS2 541/541, `PSA`→KLK3 528/528, `p62`→SQSTM1 96.2%,
`p21`→CDKN1A 89.6%. **The majority vote is wrong in four of the five** — not
noisy, but wrong in a consistent direction, picking the sense the cancer
literature never means.

One curated claim *was* wrong: the `psa` note asserted the vote returns NPEPPS
when it returns KLK3, and the committed scan already said so. It has been
corrected, and a test now requires every curated sense to be measured.

This also reframes FSP1 as the **unusual** case rather than the typical one. Its
senses are genuinely balanced (110 AIFM2 against 132 S100A4), which is why it
needs a per-paper decision while these five need only a default.

### How much to discount — `atlas_ambiguity_impact.py`

The question a reader actually has. The obvious answer is wrong by ~40×: 50.8% of
relation rows touch a contested identifier, which measures **containment**, not
error — most ESR1 edges come from papers that wrote "estrogen receptor" in full.

The measurement that means something asks whether an ambiguous form was the
*only* route to an assignment. **1.35% of relation rows** rest on one, and that
is an upper bound since the vote is sometimes right.

> **The practical lesson is an asymmetry.** Diffuse damage is ~1%, below the
> extractor's own ~79.6 F1 error, and should not change any aggregate conclusion.
> Damage to a query about one specific ambiguous entity can approach **100%** —
> GPX4–FSP1 had zero correct typed relations. Small in aggregate, total in the
> particular.

### Module support — `atlas_module_support.py` → `atlas-module-support.md`

Each `ferroptosis-core` realism layer was added on the strength of one or two
papers. This asks how many *distinct* cancer articles assert the same entity
relation, and whether the module's own cited PMID is among them. **10 of 20
corroborated** — `SLC7A11–GPX4` 31 articles, `erastin–SLC7A11` 29, `IFNG–SLC7A11`
8 with the cited PMID present, and six others resting on 1–19.

> **The denominator is hand-made.** Those 20 are author-written claims with
> author-chosen proxy entity pairs, covering 19 of roughly 30 library modules.
> "10 of 20" is therefore a statement about that curated list, not a survey of
> the library, and a different choice of proxy pairs would give a different
> fraction.

An absence here is **not** evidence against a mechanism, and the co-mention
layer proved it: **10 of the 11 modules with no asserted relation ARE discussed
in full text.** The clearest case is `fsp1` (AIFM2–GPX4), the parallel pathway
behind the manuscript's headline Bliss-synergy claim — zero relations, 69
full-text co-mentions.

### Contradictions — `atlas_contradictions.py` → `atlas-contradictions.md`

**4,667** pairs asserted in both directions, **6,764** chemicals asserted both to
`treat` and `cause` a disease. Ranked by the *weaker* side, so a 12-vs-9 split
outranks 50-vs-1.

**4 of 20 simulation modules sit on a contested edge** while citing a single
paper that does not mention the disagreement. None of that says a module is
wrong — it says the docs should state which side they took.

> This is a **reading queue, not a verdict**. The predicates are extraction
> outputs (~79.6 F1 on BioRED), the graph does not record which entity is the
> subject so direction of effect is lost, and no context is attached — a
> relation true in one cell line and false in another appears as a contradiction.

### How much of that is real? — `atlas_contradiction_quality.py`

The contradictions layer named extraction error as a caveat and never measured
it. Two failure modes are now bounded, with opposite answers.

**A single paper asserting both directions** would be extraction inconsistency
rather than disagreement between studies. It happens to **1 paper in 115,024**.
That mode is effectively absent, and the conflicts really are between papers.

**Entity ambiguity manufacturing disagreement** is the one that bites. Merging
two entities under one identifier merges two literatures, and two literatures
about different biology will disagree. Pairs involving a measured sense collision
are **1.45×** more likely to be flagged contradictory (95% CI 1.37–1.53).

> That is not a popularity artifact, which it easily could have been: colliding
> identifiers are contested *because* they are heavily mentioned, and more
> assertions means more chance of showing both directions. Stratifying by
> assertion count and pooling (Mantel-Haenszel) leaves it at 1.45× against a
> crude 1.47×, it holds inside every stratum, and it *rises* with assertion count
> — the direction conflation predicts.
>
> It is an association, not an attribution: it does not license subtracting 45%
> of the ambiguous conflicts. Check any conflict involving a blocklisted symbol
> for conflation before reading it as a scientific dispute.

### Emergence — `atlas_emergence.py` → `atlas-emergence.md`

Which relations are *new*. **4,567** pairs have ≥80% of their support since 2021,
and the top of the list is checkable: relatlimab–nivolumab (approved 2022),
belzutifan–VHL (2021), enfortumab vedotin–pembrolizumab (2023), durvalumab–biliary
tract (TOPAZ-1).

> Emergence is not importance. A relation can be new because the entity was only
> recently named, because a technique made it measurable, or because a field is
> faddish. The year range is shown so a reader can tell those apart.
>
> **Measured accuracy** (`atlas_emergence_error.py`). The share comes from the
> 60-PMID sample, so for well-supported pairs it is an estimate that a threshold
> turns into a yes/no. Against the exact share from every dated asserting paper:
> 89.4% of pairs carry no more papers than the sample holds and are exact; on the
> 10.6% genuinely estimated, median share error 0.017 and the decision is 86.4%
> precise at 93.2% recall; across all examined pairs, **99.0% precision, 99.6%
> recall**. Error grows with support (0.013 at 61–120 papers, 0.035 above 2,000),
> so do not read one pair's share as exact when hundreds of papers sit behind it.

### Discovery, and its measured hit rate — `atlas_discovery.py`, `atlas_discovery_eval.py`

Swanson ABC: if A relates to B and B to C, but A and C have never been discussed
together, the A–C link is a hypothesis the literature implies and nobody has
stated. This is the analysis the census was built for, and a 4,830-article corpus
cannot do it — discovery needs both literatures present at once, and the whole
point is that they do not cite each other.

**It has now been measured, and the ranking fails.** A time split rebuilds the
graph as it stood before year Y, ranks the absent candidates, and counts a hit
when the literature first asserts that pair in Y or later. The comparison that
matters is not against random — almost anything beats random on a clustered
co-occurrence graph — but against ranking the *same candidates* by popularity:

| ranking | degree correction | precision@20 |
|---|---|---|
| popularity | none — it *is* degree | **17.2%** |
| raw bridge count | none | 16.5% |
| Adamic-Adar | down-weights hub bridges | 16.4% |
| resource allocation | harder | 14.6% |
| **ABC (shipped)** | divides out candidate degree | **12.3%** |
| Jaccard | normalises by both degrees | 5.0% |
| random | — | 2.6% |

Paired bootstrap for ABC: −0.99 hits of 20, 95% CI [−1.27, −0.71], behind on 97
of 200 seeds, reproduced at split years 2015, 2018 and 2021.

**Nothing beats popularity**, and the rankings order themselves by how hard each
corrects for degree — the harder the correction, the worse it does. So this
cannot be repaired by swapping in a better link predictor; the standard ones were
tried. New edges genuinely do attach preferentially to well-connected entities,
so removing degree removes most of what predicts the next edge.

> **A good candidate generator and a bad ranker.** Both rankings beat random by
> ~6×, so restricting attention to 2-hop bridged entities is genuinely
> informative; the ordering within that set is what fails. `atlas_discovery.py`
> claimed to correct for popularity — measured against what the literature went
> on to say, it does not, and its docstring now says so.
>
> Nor is popularity a *good* ranking. On a graph where well-studied entities keep
> accruing edges, predicting that a famous gene gains another relation is easy and
> not very useful.
>
> **And the evaluation's own target is arguable.** It scores a ranking by whether
> it anticipates the literature, while Swanson-style discovery is *for*
> connections the literature is slow to reach — a genuinely overlooked pair scores
> here as a miss. The narrow claim is the one the module made and failed: it says
> it corrects for popularity, and doing so does not help.

## Two traps that have already bitten

Both produced plausible-looking wrong answers, which is why they are written down.

1. **Baseline files are ordered chronologically.** A partially-rebuilt census
   holds only the *oldest* literature. The coverage report once showed
   ferroptosis and CAR-T at zero, and the emergence detector once found zero
   emerging relations and dated everything to 1965-1997. Both now detect and
   announce the condition rather than reporting through it.
2. **Aggregation keys must be unique.** In `seed_replication.py`, the key
   `(treatment, o2_condition, immune_mode)` collided across three distinct
   matrix blocks, presenting between-condition difference as seed variance.

## Reproducing

```bash
python scripts/atlas_baseline.py              # the census (resumable, ~90 min)
python scripts/atlas_unindexed.py --validate  # measure the text matcher first
python scripts/atlas_unindexed.py             # then recover
python scripts/atlas_fulltext.py              # open-access full text -> NAS
python scripts/atlas_relations.py             # PubTator3 bulk
python scripts/atlas_graph.py --build         # queryable index (~86 s)

python scripts/atlas_coverage.py
python scripts/atlas_entity_audit.py          # needs network
python scripts/atlas_module_support.py
python scripts/atlas_contradictions.py
python scripts/atlas_emergence.py
```

Bulk data is gitignored. `FERRO_ATLAS_ROOT` and `FERRO_ATLAS_FULLTEXT` move it to
external storage. Unit guards that need neither network nor data:
`pytest tests/test_atlas.py`.

## What the atlas is not

It is not a knowledge base and not a source of truth. It is a map of what the
literature *says*, with measured uncertainty about how well it says it. Nothing
in it substitutes for reading the paper, and nothing in it is evidence that a
mechanism is real — only that the field discusses it, how often, in which
direction, and since when.
