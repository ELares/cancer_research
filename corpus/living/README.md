# Living review (PRISMA-LSR incremental updates, #349)

This directory holds the **living review**: a continuously-updated addendum to the
**frozen** manuscript corpus. It exists so the consolidation does not go stale,
following PRISMA for Living Systematic Reviews (PRISMA-LSR).

## Frozen vs living — the load-bearing separation

| | Frozen manuscript corpus | Living review |
|---|---|---|
| Location | `corpus/by-pmid/`, `corpus/abstracts/by-pmid/`, `corpus/INDEX.jsonl` | `corpus/living/<date>/`, `analysis/living-review/<date>.md` |
| Mutability | **immutable** — the snapshot every manuscript number is computed from | append-only dated increments |
| Produced by | the one-time corpus build (`fetch_articles.py` + `build_index.py` + `tag_articles.py`) | `scripts/living_review_update.py`, run on a schedule |

`scripts/living_review_update.py` re-runs the **same committed mechanism queries**
(`scripts/queries.txt`) against PubMed for a recent publication window, diffs the
results against the frozen corpus PMIDs, tags the **new** records with the same
mechanism keywords the frozen corpus used, and writes:

- `corpus/living/<date>/index.jsonl` — the dated incremental index (metadata +
  mechanism tags for the new records only), and
- `analysis/living-review/<date>.md` — the delta changelog (new records per query,
  new landmark detections).

It **never** writes to the frozen files, so the manuscript's numbers stay
reproducible and the living index is unambiguously an addendum.
(`tests/test_living_review.py` asserts this separation.)

## Cadence

The scheduled GitHub Action `.github/workflows/living-review.yml` runs **monthly**
(and on demand) over the trailing window and uploads the dated index + changelog
as a **workflow artifact** (it does not commit them, so `main` and the frozen
corpus stay clean). The dated outputs are therefore git-ignored; run the script
locally to materialize them:

```bash
python scripts/living_review_update.py --since 2026-01-01            # full: writes the dated index
python scripts/living_review_update.py --since 2026-01-01 --dry-run  # changelog only, no index
```

Promoting a living increment into the frozen corpus (for a future manuscript
revision) is a deliberate, manual step — not something the Action does.
