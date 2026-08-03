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
| MeSH-indexed cancer articles | 4,203,236 | `atlas_baseline.py` |
| recovered, not yet MeSH-indexed | 783,271 | `atlas_unindexed.py` |
| open-access full texts (on external storage) | 520,143 | `atlas_fulltext.py` |
| typed, normalized relations | 7,951,325 over 1,603,105 PMIDs | `atlas_relations.py` |
| queryable relation index | 2,186,309 entity pairs | `atlas_graph.py` |

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
climbs monotonically with recency, from 0.0% in 1975-era baseline files to
**37.6%** in the most recent. A pure MeSH census therefore loses exactly the
recent literature.

> **Load-bearing caveat.** These are **text-matched**, not MeSH-indexed, and live
> in a separate stream tagged `source: "text-match"`. The matcher's accuracy is
> measured against MeSH truth on indexed articles: **precision 75.7%, recall
> 95.6%**. So roughly one in four recovered records is not really cancer. Quote
> that number wherever a count from this layer is used.

### Coverage — `atlas_coverage.py` → `atlas-coverage.md`

Supplies the denominator every manuscript ratio was missing. The frozen corpus
holds **0.086%** of the census. **20.9%** of the census has a PMC id, which is
the ceiling on any full-text claim — against the frozen corpus's 98.7%
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

### Entity audit — `atlas_entity_audit.py` → `atlas-entity-audit.md`

Checks the symbols this project queries against NCBI's own record. **1 mismatch
in 53**, and it is the worst possible one: PubTator3 maps 357 mentions of `FSP1`
and 90 of "ferroptosis suppressor protein 1" to gene **51062 = ATL1 (atlastin
GTPase 1)**, a hereditary spastic paraplegia gene. The real FSP1 is **84883
(AIFM2)**.

That collision sits directly under the manuscript's headline GPX4+FSP1
Bliss-synergy claim. Querying the graph for `FSP1` silently returns ATL1's edges.

> The resolver deliberately **reproduces** this rather than patching it. Its job
> is to report what the data says; the audit's job is to catch it. Re-run the
> audit before trusting any symbol-based result.

### Module support — `atlas_module_support.py` → `atlas-module-support.md`

Each `ferroptosis-core` realism layer was added on the strength of one or two
papers. This asks how many *distinct* cancer articles assert the same entity
relation, and whether the module's own cited PMID is among them. **9 of 20
corroborated** — `SLC7A11–GPX4` 31 articles, `erastin–SLC7A11` 29, `IFNG–SLC7A11`
8 with the cited PMID present.

An absence here is **not** evidence against a mechanism: it may be inexpressible
(caveat 1 above), or in full text the extractor did not reach.

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

### Emergence — `atlas_emergence.py` → `atlas-emergence.md`

Which relations are *new*. **4,567** pairs have ≥80% of their support since 2021,
and the top of the list is checkable: relatlimab–nivolumab (approved 2022),
belzutifan–VHL (2021), enfortumab vedotin–pembrolizumab (2023), durvalumab–biliary
tract (TOPAZ-1).

> Emergence is not importance. A relation can be new because the entity was only
> recently named, because a technique made it measurable, or because a field is
> faddish. The year range is shown so a reader can tell those apart.

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
