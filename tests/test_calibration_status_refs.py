"""Guard the CALIBRATION_STATUS.md <-> targets.yaml cross-references (#588 audit).

CALIBRATION_STATUS.md calls targets.yaml the "machine-checked authority" and
cites specific target IDs as `targets.yaml: <id>`. Five of those IDs had drifted
to names that do not exist in targets.yaml (gpx4_recovery_rate,
mufa_protection_factor, fsp1_hdac_persister, pdt_depth_attenuation,
sdt_depth_attenuation), so the "machine-checked" claim was hollow. This pins the
invariant: every `targets.yaml: <id>` reference must resolve to a real target id.
"""

import re
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
STATUS = REPO / "simulations" / "calibration" / "CALIBRATION_STATUS.md"
TARGETS = REPO / "simulations" / "calibration" / "targets.yaml"

# `targets.yaml: <id>` (the cross-reference convention used in the Status doc).
# IDs are lower snake_case and may start with a digit (e.g. 3d_immune_sdt_dominates).
_REF = re.compile(r"targets\.yaml:\s*([0-9a-z][0-9a-z_]*)")


def test_every_calibration_status_target_ref_resolves():
    ids = {t["id"] for t in yaml.safe_load(TARGETS.read_text(encoding="utf-8"))["targets"]}
    refs = set(_REF.findall(STATUS.read_text(encoding="utf-8")))
    assert refs, "expected at least one `targets.yaml: <id>` reference in CALIBRATION_STATUS.md"
    dangling = sorted(r for r in refs if r not in ids)
    assert not dangling, (
        f"CALIBRATION_STATUS.md references target IDs absent from targets.yaml: {dangling}. "
        f"Real IDs: {sorted(ids)}"
    )
