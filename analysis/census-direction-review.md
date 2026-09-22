# Direction candidates prepared for independent review

**Unreviewed: no biological decisions have been assigned.** The local packet contains complete blank worksheets for the singly classified hypoxia and thesis candidates. Historical reports and labels remain separate.

Captured 2026-09-22T16:08:52.756432Z with Python 3.14.6. The indexed `records/` snapshot has 4,403,994 parsed records across 1,334 shards; its manifest declares 4,403,994. Every sorted shard was read (`stride=1`). This validates the acquired local snapshot, not its coverage of literature published since acquisition.

The local packet retains 765 complete records in the union of the two descriptor intersections. This public report contains counts and fingerprints without article text.

The indexed stream can include C04 and adjacent records; the direction scanners do not add a C04-only restriction. Snapshot cancer-basis counts: `C04` 4,203,236, `adjacent` 200,758.

## Candidate coverage

| Analysis | Intersection | Singly classified candidates | Both/neither excluded |
|---|---:|---:|---:|
| hypoxia | 291 | 32 | 259 |
| thesis | 479 | 266 | 213 |

All four lexical classes are retained in the intersection. Only singly classified candidates enter the exported worksheets; these counts do not measure biological direction or classifier accuracy. Even complete candidate adjudication cannot establish recall among excluded or unindexed articles.

| Analysis | Lexical class | Records | Missing title | Missing abstract | Missing both |
|---|---|---:|---:|---:|---:|
| hypoxia | protects | 12 | 0 | 0 | 0 |
| hypoxia | sensitises | 20 | 0 | 0 | 0 |
| hypoxia | both | 5 | 0 | 0 | 0 |
| hypoxia | neither | 254 | 0 | 0 | 0 |
| thesis | exploit | 237 | 0 | 1 | 0 |
| thesis | obstacle | 29 | 0 | 0 | 0 |
| thesis | both | 37 | 0 | 0 | 0 |
| thesis | neither | 176 | 0 | 1 | 0 |

Missing text includes empty and whitespace-only fields. Missing both is included in each separate missing-field count. Exported evidence must remain unchanged: record later source checks and any additional evidence in review reasons or a sidecar keyed by the record fingerprint.

## Reproduction and review

Generate a new local packet outside the repository using the same acquired census and code:

```bash
python scripts/census_direction_review.py --atlas-root "$FERRO_ATLAS_ROOT" \
  --bundle-dir /path/outside/repository/direction-review
python scripts/census_direction_review.py \
  --verify-bundle /path/outside/repository/direction-review
python scripts/census_direction_review.py --render-only
```

The packet directory must not already exist. All parsing, hashing, serialization and rendering finish before publication. The packet uses a sibling staging directory and rename; the subsequent public JSON and Markdown writes are not an atomic pair against disk failure or interruption.

Offline verification checks the fixed packet paths and hashes, recorded source-code fingerprints, and all regenerated scanner reports, blank worksheets, cohort fingerprints and readiness counts against the retained intersection records. It does not reread the full census, prove that the retained union is complete, authenticate provenance, or assess scientific labels. Full input verification requires the original shards and a new build.

Independent readers and source checks remain outstanding. Preserve separate reviewer decisions, identities, evidence sources and locators, and documented disagreement resolution. Only completed identity-matched candidate worksheets can be imported using `docs/CENSUS_DIRECTION_ADJUDICATION.md`; historical title-only labels cannot be carried forward.

## Provenance

Source URL recorded by the snapshot: `https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/`.

Source inventory SHA-256: `ebf586c5dd2d276b70fa688c818103454570534b57575cb5c943e917eb7f0224`.

Acquisition manifest SHA-256: `c6d300992fcbb336c6d0b0927442fb7faed88bd940b18bc447d3aff6d3e49484`.

Baseline manifest SHA-256: `6d56da3c5ba6d9b53b457abfe1176effa3ce260b584c3f2fbfe5d8167523814b`.

Git context: `257e9523d5a37b42ddd69b24fac1846c29cd8b92`; working tree dirty: `False`. Source hashes below are authoritative for the code used; the Git head alone may predate uncommitted code.

| Analysis | Cohort SHA-256 |
|---|---|
| hypoxia | `2387ecd78216b8481f953e2f1f841825857f0e6551467a7ebe67c017a1c86cf3` |
| thesis | `aec7e435ddcfa66f7f90edc0726a0b90a7f15dfca2b367df4412297e0bc111df` |

| Source file | SHA-256 |
|---|---|
| `scripts/census_direction_review.py` | `f76dbb656ea3a3cac1014f3097123d4f6570d6919c6ecc28e200a11329d1bcd4` |
| `scripts/census_snapshot.py` | `783aeb63a4f5a9da471fd9be66520cacc62b3b69e41758b1c4bae39b886e9db6` |
| `scripts/census_hypoxia_direction.py` | `8ed7db36ae862167e5a16d64f1f546915b0e3367a13596e86d1286ccfffaae54` |
| `scripts/census_thesis_direction.py` | `02cab7828710b33e132629e49a89283e94b9866f85bab0ca8313af64128c65c1` |
| `scripts/census_adjudication.py` | `e63412490d2342476e2f9cb9a4080999597c98e3c9023f1850e71177b0630579` |
| `scripts/census_input.py` | `517578586090750adb46e1ec6aad2e80ff57d2b793afa5013f6df1dd8d9cd910` |
| `scripts/atlas_baseline.py` | `9ac870db4bb4a89e280ded514bffc4340ae7612a608fbb77ef6261bcd6a542ca` |

| Local packet payload | SHA-256 |
|---|---|
| `intersection-records.jsonl` | `63952a1fc5c59aebc4d261b618b63b7b6855c91299e01d4d16c11aff45a0edc6` |
| `hypoxia/census-hypoxia-direction.json` | `e3e0074dfaa6dfaaed28c1a1a607efd4895ff9ec05b690fab27ba6a724785a1a` |
| `hypoxia/census-hypoxia-direction.md` | `1d238ae9f7b7076b171e126f659fe4177d92881feb44615bfaa33b45a9b566ba` |
| `hypoxia/candidates.csv` | `84040d1a3729ae43e2755f2f688b382c0b95be0a3add31c425d37fe3150aa126` |
| `thesis/census-thesis-direction.json` | `cd7ba12cbf7eb00ee08bf124b85543d3ee05cbb4eefa3fc8c5bde187f5bbf445` |
| `thesis/census-thesis-direction.md` | `7f0849ff1c8c1a26b19d6091f7a49d340919ceb59bcb150e1a7fceda059a6a21` |
| `thesis/candidates.csv` | `48b93863feb07132a3c756eb4ae945e6a245cc685a6510cfffacd261b0a27316` |
