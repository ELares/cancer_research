# Living review (PRISMA-LSR incremental updates, #349)

This directory holds the **living review**: a continuously-updated addendum to the
**frozen retrieval** used for historical results and method comparisons. It
tracks recent papers matched by that retrieval's committed query set,
following PRISMA for Living Systematic Reviews (PRISMA-LSR).

The current manuscript's literature counts and figures instead use the cancer
census, described in [`analysis/atlas-README.md`](../../analysis/atlas-README.md).
This query-based addendum does not update the census or supply its denominator;
each census analysis declares which census streams it reads.

## Frozen vs living — the load-bearing separation

| | Frozen retrieval | Living review |
|---|---|---|
| Location | `corpus/by-pmid/`, `corpus/abstracts/by-pmid/`, `corpus/INDEX.jsonl` | `corpus/living/<date>/`, `analysis/living-review/<date>.md` |
| Mutability | **immutable** — preserves historical results and method comparisons | dated increments; rerunning on the same date replaces that date's output |
| Produced by | the one-time corpus build (`fetch_articles.py` + `build_index.py` + `tag_articles.py`) | `scripts/living_review_update.py`, run on a schedule |

`scripts/living_review_update.py` re-runs the **same committed mechanism queries**
(`scripts/queries.txt`) against PubMed for a recent publication window, diffs the
results against the frozen corpus PMIDs, tags the **new** records with the
current committed mechanism keywords, and writes:

- `corpus/living/<date>/index.jsonl` — the dated incremental index (metadata +
  mechanism tags for the new records only), and
- `analysis/living-review/<date>.md` — the delta changelog (new records per query,
  new landmark detections).

It **never** writes to the frozen files, so historical results remain
reproducible and the living index is unambiguously an addendum. "New" here means
absent from the frozen retrieval; it does not mean absent from earlier living
increments or from the census.
(`tests/test_living_review.py` asserts this separation.)

`python scripts/corpus_identity_index.py` includes each local dated
`index.jsonl` when rebuilding the download-deduplication index. Its PMID, PMC ID
and DOI keys count as held metadata, without claiming that full text is held.

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

Changing the frozen retrieval or incorporating these records into a census
analysis requires a separate, explicit revision; the Action does neither.
