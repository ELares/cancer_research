# Revised proposal validation on known synthetic targets

Revision `local-move-v2` uses `resample_move_local` and analytic scores only.

**All fixed synthetic checks passed.**

The [prior study](../../analysis/calibration/proposal-synthetic-validation.json) remains separate; its JSON SHA256 is
`1b045df4b719b1131c9c15256e0a803959c4b57a7df787e17cbf26e311ae3e82`.

The unchanged seven-dimensional fixtures are A = [0, 1/4]^7 and
B = [1/4, 1] × [5/8, 7/8]^6 (mass 1/4096; mode masses 1/4 and 3/4),
plus the boundary slab [0, 1/8] × [0, 1]^6 (mass 1/8). Exact moments and
the same fixed error gates are stored in JSON. Fixture derivations are in the prior report.

Each learned run uses 8,192 independent production attempts.
Pilot particles do not enter the estimates. The fixed pilot plan is:

```json
{
  "islands": 2,
  "kernel_ceiling": 0.5,
  "kernel_floor": 0.005,
  "kernel_neighbors": 32,
  "kernel_scale": 1.5,
  "levels": 14,
  "local_probability": 0.5,
  "local_scale_ceiling": 0.25,
  "local_scale_floor": 0.0025,
  "local_scale_multiplier": 0.5,
  "moves_per_level": 4,
  "neighbor_gap_ratio": 3.0,
  "particles_per_island": 256,
  "retained_fraction": 0.5,
  "singleton_scale": 0.04,
  "uniform_weight": 0.2
}
```

## Learned proposals

| Fixture | Seed | Accepted | ESS | Relative mass error | Largest moment error | Result |
|---|---:|---:|---:|---:|---:|---|
| boundary_slab | 2026092401 | 6094 | 5470.9 | 0.016 | 0.0068 | pass |
| boundary_slab | 2026092402 | 6004 | 5386.5 | 0.005 | 0.0061 | pass |
| boundary_slab | 2026092403 | 6064 | 5287.6 | -0.001 | 0.0072 | pass |
| separated_boxes | 2026092401 | 1333 | 1000.1 | -0.026 | 0.0070 | pass |
| separated_boxes | 2026092402 | 1293 | 831.6 | -0.044 | 0.0223 | pass |
| separated_boxes | 2026092403 | 1275 | 923.3 | -0.036 | 0.0259 | pass |

## Reused negative controls

The three previously inspected proposals centered only on B are retained with their
20% uniform component. These are known failure controls, not fresh validation runs.

| Seed | A hits | B hits | ESS | Usual screens | Known-truth screens |
|---|---:|---:|---:|---|---|
| 2026092001 | 0 | 3730 | 1061.3 | pass | fail |
| 2026092002 | 0 | 3770 | 823.9 | pass | fail |
| 2026092003 | 0 | 3751 | 676.4 | fail | fail |

The reused fixtures informed this revision. Fresh random seeds do not make these
independent validation of general mode discovery. The unchanged oracle precision bound
applies only to its ideal exact-target proposal, not to these learned kernels.
High ESS and agreement cannot exclude an unseen region. These checks establish neither
biological validity nor coverage of unknown biological acceptance regions.
