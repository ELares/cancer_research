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

## Export a new review worksheet

Set `FERRO_ATLAS_ROOT` to an acquired census containing `records/*.jsonl.gz`.
Use a separate output directory so a new scan does not replace the published
historical files:

```bash
python scripts/census_hypoxia_direction.py \
  --output-dir scratchpad/hypoxia-review \
  --export-candidates scratchpad/hypoxia-review/candidates.csv

python scripts/census_thesis_direction.py \
  --output-dir scratchpad/thesis-review \
  --export-candidates scratchpad/thesis-review/candidates.csv
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
  --output-dir scratchpad/hypoxia-review \
  --adjudication scratchpad/hypoxia-review/completed.csv

python scripts/census_thesis_direction.py \
  --output-dir scratchpad/thesis-review \
  --adjudication scratchpad/thesis-review/completed.csv
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
