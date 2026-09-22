# Original numerical source bytes for historical joint-sampler archives

The importance and resample/move archives identify their numerical inputs with
SHA-256 hashes. They predate the controlled immune study and must retain those
original identities when current source files change.

The two `.source` files here were recovered from the tracked Git object
`cd8374d3a6dafacdc9ed2bb5d612a0efc1943d98`, at the paths in `manifest.json`.
Every recovered byte sequence was checked against the hash already present in
the historical resample archives; the immune source hash also appears in all
three importance archives. This is retrospective source preservation, not a
claim that this directory existed at the original experiment freeze. Existing
raw archives and generated reports are unchanged. No simulation was rerun and
no new numerical outcomes were substituted.

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
