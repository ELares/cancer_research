# Adjudicating direction candidates from a census

The hypoxia and thesis direction scanners generate keyword candidates. Their
historical classifiers made direction errors, so an unreviewed scan cannot
establish which biological direction leads. New runs do not reuse the committed
title-only adjudication CSVs.

## Preserve and inspect the historical results

```bash
python scripts/census_hypoxia_direction.py --render-only
python scripts/census_thesis_direction.py --render-only
```

These commands render the embedded results in the corresponding
`analysis/census-*-direction.json` without acquiring source data, reading the
adjudication CSV, recalculating its correction, or rewriting the JSON. Historical
CSV labels remain available for inspection. Neither artifact records a complete
candidate inventory with article identifiers, so their historical record linkage
cannot be verified from the committed data. Equal aggregate counts and matching
title prefixes do not supply that missing evidence.

The thesis snapshot's 74.6% correction and 53.8–86.5% range are retained as a
historical conditional calculation. The range propagates the decided exploit
sample's precision interval while treating adjudicated obstacle counts as fixed;
it does not cover all sources of uncertainty or establish that the sample
represents a new population. The hypoxia snapshot's 11 protective and 7 sensitising
labels describe reviewed titles, not a census-wide direction estimate.

## Prepare a verified local review packet

For a complete baseline review, use the packet builder before filling any
decisions. It requires the acquisition `manifest.json` and the indexed
`records/*.jsonl.gz` stream. It reads every shard (`stride=1`), checks the exact
baseline inventory and declared record counts, and hashes the same compressed
bytes it parses. It rejects unreadable, malformed, missing, extra, aliased, or
changing shards. Update, unindexed, and C04-only streams are not merged.

```bash
export FERRO_ATLAS_ROOT=/path/to/atlas
python scripts/census_direction_review.py \
  --atlas-root "$FERRO_ATLAS_ROOT" \
  --bundle-dir /path/outside/repository/direction-review

python scripts/census_direction_review.py \
  --verify-bundle /path/outside/repository/direction-review
```

The destination must not exist, must be outside the repository, and must not
overlap census inputs or public reports. The complete packet is prepared before
publication through a sibling staging directory and rename. It contains:

| Local file | Purpose |
|---|---|
| `intersection-records.jsonl` | Complete union of both descriptor intersections, including both/neither keyword classes. |
| `hypoxia/census-hypoxia-direction.json`, `.md` | Fresh, unreviewed hypoxia report. |
| `hypoxia/candidates.csv` | Every singly classified hypoxia candidate with blank decisions and reasons. |
| `thesis/census-thesis-direction.json`, `.md` | Fresh, unreviewed thesis report. |
| `thesis/candidates.csv` | Every singly classified thesis candidate with blank decisions and reasons. |
| `readiness.json`, `readiness.md` | Snapshot inventory, counts, missing-text coverage, and source/payload fingerprints. |

Article text stays in this local packet. Only the text-free readiness pair is
also written to `analysis/census-direction-review.json` and `.md`. Those two
public writes are not an atomic transaction against disk failure or interruption.
The historical hypoxia and thesis artifacts are not overwritten. Public prose
can be rebuilt offline with `python scripts/census_direction_review.py --render-only`.

The [September 2026 readiness report](../analysis/census-direction-review.md)
reconciles 4,403,994 records across 1,334 shards: 4,203,236 C04 records and
200,758 adjacent records. The scanners add no C04-only restriction. The packet
retains 765 records in the union of the intersections: 291 hypoxia records and
479 thesis records, with five shared records. Their worksheets contain 32 and
266 singly classified candidates respectively. One thesis candidate lacks an
abstract; another missing abstract is in the excluded neither class. All
intersection records have titles. These are keyword workloads, not reviewed
biological findings or evidence of current literature coverage.

Source hashes identify the exact builder, reader and scanner code. For later
verification, use a checkout matching those hashes; the recorded Git commit
provides context, and this published packet was built from a clean checkout.
`--verify-bundle` needs no original census: it checks the fixed file inventory,
source and payload hashes, internal snapshot count consistency, and regenerated
reports, blank worksheets, and cohort fingerprints from the retained records.
It cannot revalidate absent shards, prove the union's completeness, authenticate
the provenance, establish global PMID uniqueness, or assess biological labels.
Repeat the full build with the original census to recheck its inputs. The two
analysis cohorts independently reject duplicate records and repeated nonempty
PMIDs.

Preserve the original packet unchanged. **Copy each candidate worksheet outside
the packet before entering reviewer decisions**; changing its blank CSV makes
packet verification fail. Keep separate initial decisions from independent
readers, reviewer identities, evidence sources and locators, and a documented
disagreement-resolution record. The initial reading can hide the keyword-label
column to reduce anchoring, but the final import must preserve that column and
all identity fields. Record any additional source evidence in the reason or a
sidecar keyed by the record fingerprint; do not replace the exported title or
abstract. The retained both/neither records support a separate missed-candidate
audit, but are not silently added to the current directional denominator.

Independent source review and completed adjudication remain outstanding. The
import workflow below rereads the same census and checks the completed copy
against its current cohort; offline packet verification is not an adjudication
import and does not assign decisions.


## Export a new review worksheet

Set `FERRO_ATLAS_ROOT` to an acquired census containing `records/*.jsonl.gz`.
Use a separate output directory outside the repository so article text stays
local and a new scan does not replace the published historical files:

```bash
python scripts/census_hypoxia_direction.py \
  --output-dir /path/outside/repository/hypoxia-review \
  --export-candidates /path/outside/repository/hypoxia-review/candidates.csv

python scripts/census_thesis_direction.py \
  --output-dir /path/outside/repository/thesis-review \
  --export-candidates /path/outside/repository/thesis-review/candidates.csv
```

The worksheet contains every singly classified candidate's full title and
abstract, PMID when available, original keyword label, record fingerprint, and
cohort fingerprint. Decisions and reasons start blank. The JSON report also
retains the candidate evidence. Records matching both patterns or neither
remain in the intersection counts and cohort fingerprint, but are excluded
from this worksheet's directional denominator.

`--stride N` selects every Nth sorted shard; it does not select every Nth record
or create a random sample. Export and import must use the same selected
intersection. Reordering identical records or moving them to another root does
not invalidate the fingerprints. Changing an article's content, changing the
selected intersection, or changing a keyword classification does. Duplicate
records and repeated nonempty PMIDs within the intersection are rejected.

## Complete and import the review

Fill only `adjudicated` and `reason` for every exported row. Keep all identity,
text, and keyword-label columns unchanged. A reason is required even when a
decision is ambiguous. The allowed decisions are:

| Analysis | Decisions |
|---|---|
| Hypoxia direction | `protects`, `sensitises`, `off-topic`, `ambiguous` |
| Thesis direction | `exploit`, `obstacle`, `ambiguous` |

Document the evidence reviewed, reviewer identities, and any disagreements
alongside the completed worksheet. Independent review and source checks are
needed to assess judgment accuracy; the importer checks identity and coverage,
not the scientific validity of a label.

```bash
python scripts/census_hypoxia_direction.py \
  --output-dir /path/outside/repository/hypoxia-review \
  --adjudication /path/outside/repository/hypoxia-review/completed.csv

python scripts/census_thesis_direction.py \
  --output-dir /path/outside/repository/thesis-review \
  --adjudication /path/outside/repository/thesis-review/completed.csv
```

Every singly classified record must appear exactly once with matching cohort,
record, text, and keyword-label fields. Missing, extra, duplicate, unfinished,
or mismatched rows stop the run before report replacement. Legacy title-only
CSVs are rejected. These imports use direct counts of the complete current
candidate set, without extrapolating a historical sample's precision. Ambiguous
and off-topic decisions are excluded from the directional denominator and
reported explicitly; a cohort with no directional decisions has no directional
share. Complete candidate review still cannot measure recall among records that
the keyword patterns missed.

`--adjudication` and `--export-candidates` are mutually exclusive and cannot be
combined with `--render-only`. A later `--render-only --output-dir ...` renders
the stored report without relabeling it. Missing or empty selected census input
is an error; a readable census with no intersection or no singly classified
candidates is a valid result. An explicit adjudication import with no candidates
is rejected rather than producing a percentage from unrelated rows.

JSON serialization, Markdown rendering, and any worksheet preparation finish
before output replacement. Parsing, coverage, and rendering failures preserve
existing reports. Separate writes are not an atomic transaction against disk
failure or interruption. Fingerprints detect mismatches; they do not authenticate
an adjudicator, validate a biological conclusion, or recover missing historical
record identities.
