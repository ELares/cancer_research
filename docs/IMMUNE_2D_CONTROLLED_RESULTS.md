# Controlled sources change activation exposure at matched DAMP mass

All 18 conditions and 36 contrasts in the
[frozen protocol](IMMUNE_2D_CONTROLLED_PROTOCOL.md) completed on the first capture
attempt, without retries, dropped conditions or outcome-dependent methods
changes. Every equal-total-mass comparison gave greater accumulated activation
exposure for the declared 256-source pattern than for the 64-source pattern,
while its maximum eligible DAMP was lower. This is a result for these fixed
source and recipient patterns under the existing 2D rules, not a universal
density law or an explanation of the historical SDT:RSL3 immune-kill ratio.

The [generated report](../analysis/immune-2d-controlled.md) publishes every
condition and directed contrast. The [complete derived data](../analysis/immune-2d-controlled.json)
also retain the conditional DAMP means, threshold fractions and alternative
recipient normalization. All quantities below refer to that complete design;
none was selected to replace a prespecified endpoint.

## What the controlled comparisons show

The primary response sums `DAMP/(DAMP+50)` at available recipient sites with
DAMP at least 0.01 during steps 60–179, then divides by the same 3,840 reference
recipients in every condition. It measures accumulated activation steps per
reference recipient. It is neither a probability nor a realized immune-kill
count: sources are imposed pulses and recipients remain fixed, without
biochemistry, random immune draws or endogenous death.

All eight equal-total-mass contrasts are positive. Dividing the same injected
budget among 256 declared source positions increased this response relative to
64 positions, with timing and recipients held fixed. The lower eligible maxima
in those same comparisons show why peak DAMP alone is not the accumulated
exposure endpoint. Source positions and release per source both change in this
comparison. It establishes the effect of these specific budget-partitioning
interventions; it does not isolate a source-count effect while also holding
per-source release fixed.

All eight release-amount contrasts are positive when the budget increases at
fixed source pattern, timing and recipients. The four comparisons with
per-source release fixed at 20 model units are also positive, but they increase
both source count and total injected mass. Those comparisons answer a different
question from the equal-total-mass contrasts. The eligibility floor and transport
cutoff remain active, so these results do not establish general proportionality
or unrestricted linear superposition.

All eight timing contrasts favor release at step 60 over step 30. Step 60 is
the start of the fixed observation window; a step-30 pulse has already undergone
thirty additional transport and clearance updates when that window opens.
The comparison concerns exposure within this imposed window, not a calibrated
biological release schedule or a general advantage of later treatment.

For each source/budget/timing combination, the 960-recipient mask has exactly
one quarter of the full mask's eligible opportunities and, within the frozen
floating-point tolerance, one quarter of its primary exposure. The field is
unchanged, and retained recipients have identical observations in both arms.
Normalizing by each arm's available-recipient count gives the same exposure
within that tolerance. This reflects the particular balanced subset used here;
it is not a rule for arbitrary missing recipients and does not model survival
selection or depletion by killing.

No eligible opportunity reached activation 0.5 or 0.9 in any condition. The
largest eligible DAMP was about **6.26263 model units**, below the threshold 50
for half-maximal activation. Both zero controls have zero field, opportunities
and exposure; their conditional means and maxima remain undefined. These null
threshold results apply to the declared masks, budgets and timing. They do not
prove that other settings cannot reach the thresholds.

## Freeze, compatibility and scope

The complete implementation and protocol were committed before production at
`421603d2a347b8db6087287ed8cba426ba5dfe1c`. The
[archive manifest](../analysis/immune-2d-controlled/manifest.json) records the
source inventory, binary identity, pinned Rust 1.96.0 toolchain, Python version,
platform, capture times and all payload hashes. Its controlled-run log contains
the 18 conditions once each, in their declared order.

After the capture, an isolated rebuild of the archived source reproduced both
the executable hash and the complete controlled observation bytes. The frozen
reader also reconstructed both reports exactly, and its bundled Python tests
passed with the expected archive-absent test skipped. This was a deterministic
reproduction check of already observed conditions, not a retry, new study or
independent replicate.

Before the controlled conditions ran, the same executable reproduced the
historical 33-condition summary with SHA-256
`e04f9699ddca9b84145cdecbcd8e21abb0bc760dc3d909908358d4554a5eee98`
and the canonical 2D observer output byte for byte. The existing golden,
canonical 2D and 3D archives, and twenty-block replication archive are unchanged.
The new archive reconciles injections and clearance, per-step/per-recipient
exposure ledgers, recipient identities, availability invariance and every
declared contrast.

This deterministic design supplies no independent-replicate uncertainty,
confidence interval or biological validation. DAMP remains an uncalibrated
model field, and steps have no new conversion to hours. Its controlled
interventions do not estimate how much each mechanism contributes to the
original treatment contrast. P5 and its experimental falsification threshold
remain unchanged; independent assay mapping and raw biological replicates are
still required.

## Reconstruct the reports without simulations

From the repository root:

```bash
python3 scripts/immune_2d_controlled_report.py
```

This validates the archived hashes and frozen inputs, reconciles the ledgers,
and regenerates both derived reports. It requires no Rust execution. The archive
contains per-step and per-recipient aggregates, not every local field value or
each recipient's eligibility chronology. Reconstruction therefore cannot
independently recover transport or every individual activation evaluation.
FNV fingerprints check consistency between field trajectories; SHA-256 checks
payload integrity. Neither substitutes for a transport replay or biological
validation.

## Replay the frozen source in temporary storage

The following uses the trusted source bundle after archive validation, runs its
synthetic Python tests, and builds its executable with the pinned toolchain.
Git initialization is needed because the source-inventory tests use
`git ls-files`; it does not modify the original repository. The one test that
requires the production archive is skipped in the isolated source snapshot.
Python 3.10 or later, pytest, Git, and Rust 1.96.0 through rustup are required.
Leave every `FERRO_*` environment variable unset.

```bash
(
  set -eu
  repo_root="$(git rev-parse --show-toplevel)"
  replay_dir="$(mktemp -d)"
  source_dir="$replay_dir/source"
  python3 - "$repo_root" <<'PY'
from pathlib import Path
import sys
repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / "scripts"))
import immune_2d_controlled_report as report
report.load_archive(repo / "analysis/immune-2d-controlled")
PY
  mkdir "$source_dir"
  tar -xzf "$repo_root/analysis/immune-2d-controlled/sources.tar.gz" \
    -C "$source_dir"
  git -C "$source_dir" init -q
  git -C "$source_dir" add .
  (
    cd "$source_dir"
    python3 -m pytest -q tests/test_immune_2d_controlled_report.py
  )
  python3 - "$source_dir" "$replay_dir" "$repo_root" <<'PY'
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys

source, replay, repo = map(Path, sys.argv[1:])
assert not any(name.startswith("FERRO_") for name in os.environ)
sys.path.insert(0, str(source / "scripts"))
import immune_measurement_report as common

# The frozen helper selects Cargo's reported executable, including custom targets.
built = common.build_binary(
    ["rustup", "run", "1.96.0", "cargo"], source / "simulations", "sim-tme")
binary = replay / "sim-tme"
shutil.copy2(built, binary)
run = replay / "run"
run.mkdir()
subprocess.run([str(binary), "--immune-controlled"], cwd=run, check=True)
actual = (run / "output/tme/immune_controlled.json").read_bytes()
expected = gzip.decompress(
    (repo / "analysis/immune-2d-controlled/observations.json.gz").read_bytes())
assert actual == expected, "Frozen replay differs from archived observation bytes"
print(f"Frozen observations reproduced byte for byte; replay retained at {replay}")
PY
)
```

This re-executes the already declared controlled conditions; it is not a new
study or an independent replicate. The original archive and reports remain
untouched. The manifest records the capture platform for comparison if a replay
on another toolchain or platform differs; retain and investigate any mismatch
instead of updating the frozen observations to admit it.
