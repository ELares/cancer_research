# Original numerical source bytes for historical joint-sampler archives

The importance and resample/move archives identify their numerical inputs with
SHA-256 hashes. They predate the controlled immune study and must retain those
original identities when current source files change.

The `.source` files here were recovered from tracked Git objects, at the paths
and commits recorded in `manifest.json`. Every recovered byte sequence was
checked against the hash already present in the historical archives. This is
retrospective source preservation, not a claim that this directory existed at
the original experiment freeze. Existing raw archives and generated reports
are unchanged. No simulation was rerun and no new numerical outcomes were
substituted.

The first two came from `cd8374d3a6dafacdc9ed2bb5d612a0efc1943d98`; the immune
source hash also appears in all three importance archives.

`simulations/Cargo.lock` was added when the October dependency batch bumped
rand, serde, clap and pyo3. Seven archives -- three importance and three
resample, plus the offline rebuild check -- identify the lock among their
numerical inputs, so changing it made their declared input unresolvable and
`test_archive_source_hashes_cover_all_frozen_numerical_inputs` failed with
`no registered historical source`. The bytes were recovered from
`4f9f76c48269f134a011190412c4d7d8fce06c19`, the commit that last set them,
rather than from a branch tip that keeps moving. This is the gate working as
designed: a dependency bump is a change to a declared input of a frozen
result, and the result's provenance has to survive it. Nothing about the
archived numbers changed -- what changed is that their inputs are now
preserved where the resolver can still verify them.

`scripts/archived_numerical_sources.py` first checks the current file against
the requested hash. If it differs, the resolver requires an exact registered
pair of original path and hash, reads the corresponding content-addressed
snapshot, and checks its bytes. Missing, unregistered or corrupted historical
sources fail validation. A matching filename or a declared hash alone is not
accepted as evidence.

The original `abc_joint_resample.py` snapshot is the numerical program recorded
by those archives. Today's reader adds historical-source verification while
reconstructing the same reports; it does not relabel its current source as the
original program. The resolver is imported only when offline reconstruction
encounters differing source identities. Matching-source captures still hash
their current numerical inputs before and after sampling, and do not execute
this fallback. The resolver and this registry are dependencies of historical
archive verification, not additional numerical inputs to those old runs.
