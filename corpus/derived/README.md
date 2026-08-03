# Derived corpus layers

This directory holds **re-derivations of the frozen corpus**: the same PMIDs,
re-tagged by a newer method, kept separate so the frozen snapshot stays
reproducible.

It is the third of three corpus surfaces. The separation is load-bearing.

| | Frozen manuscript corpus | Living review | Derived layers |
|---|---|---|---|
| Location | `corpus/by-pmid/`, `corpus/abstracts/by-pmid/`, `corpus/INDEX.jsonl` | `corpus/living/<date>/` | `corpus/derived/<layer>/` |
| Contains | the snapshot every manuscript number is computed from | NEW records found after the freeze | the SAME records, re-tagged by a newer method |
| Mutability | **immutable** | append-only dated increments | regenerated in place from the frozen corpus |
| Produced by | the one-time corpus build | `scripts/living_review_update.py` | a named build script per layer |

A derived layer never writes to `corpus/INDEX.jsonl`, `corpus/by-pmid/` or
`tags/`. Anything that must reproduce a manuscript number reads the frozen
index; anything that wants the better estimate joins a derived layer on `pmid`.

## Layers

### `evidence-v2/` — evidence tiers from the v2 tagger (#TAGGER-V2)

Built by `scripts/build_evidence_v2_index.py`. One JSON object per line:

```json
{"pmid": "40700574", "evidence_level_v2": "preclinical-invivo",
 "evidence_level_frozen": "", "changed": true}
```

The frozen tagger read only title, MeSH and abstract, so it left 57.8% of
records with no evidence level. The v2 tagger reads the Methods and Results
sections of the stored full text and cuts exact-label error 2.59x on held-out
records (`analysis/evidence-v2-eval.md`). Against the frozen tags it changes
1,058 of 4,830 records (21.9%): 634 gain a tag, 20 lose one, and `theoretical`
rises from 43 to 158.

**These are a better estimate, not a correction.** Annotator agreement on this
task is 77%, and the tiers that move most here — `theoretical`, and the
in-vitro / in-vivo boundary, where 230 records shift — are precisely the ones
two annotators disagree about most. 2,178 records (45%) remain untagged.

Full delta: `analysis/evidence-v2-corpus-delta.md`. Method and caveats:
`analysis/evidence-tagger-v2.md`. Promotion checklist at the end of that file.
