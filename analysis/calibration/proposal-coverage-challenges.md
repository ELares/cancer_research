# Coverage challenges on new synthetic geometries

**Prespecified checks failed. All planned outcomes are retained.**

The local-move sampler is unchanged. These geometries and all seeds, budgets
and criteria were committed before evaluation. They were not used to develop
the sampler, but were deliberately designed with knowledge of its earlier
failure modes. This is a bounded challenge study, not an unbiased sample of
possible targets or independent biological validation.

See the [frozen protocol](../../docs/COVERAGE_CHALLENGE_PLAN.md) for exact
geometry, analytic masses, moments, partition rules and complete seed lists.
Each learned run uses 29,184 pilot attempts and 8,192 independent production
attempts. No pilot particles contribute to estimates; no model dose calls
are made. All proposal weights use the full normalized mixture.

## Learned proposals

| Geometry | Seed | Accepted | ESS | Relative mass error | Maximum region error | Maximum moment error | Result |
|---|---:|---:|---:|---:|---:|---:|---|
| annular_cylinder | 2026092511 | 685 | 426.8 | -0.0055 | 0.0315 | 0.0328 | pass |
| annular_cylinder | 2026092512 | 727 | 514.8 | +0.0609 | 0.0293 | 0.0290 | pass |
| annular_cylinder | 2026092513 | 654 | 496.3 | +0.0142 | 0.0132 | 0.0131 | pass |
| rotated_box | 2026092501 | 509 | 317.8 | -0.0103 | 0.0466 | 0.0211 | pass |
| rotated_box | 2026092502 | 450 | 221.1 | -0.0418 | 0.0823 | 0.0483 | fail |
| rotated_box | 2026092503 | 399 | 213.3 | -0.0736 | 0.0329 | 0.0390 | fail |
| unequal_balls | 2026092521 | 1526 | 1184.1 | -0.0299 | 0.0087 | 0.0129 | pass |
| unequal_balls | 2026092522 | 1503 | 1162.8 | -0.0049 | 0.0171 | 0.0161 | pass |
| unequal_balls | 2026092523 | 1511 | 1205.3 | -0.0072 | 0.0211 | 0.0165 | pass |

Failed checks:

- rotated_box, seed 2026092502: maximum_weight, region_masses.
- rotated_box, seed 2026092503: maximum_weight.

## Analytic controls

Full-target oracles draw directly from known analytic targets. Support-hole
oracles deliberately omit half of the rotated or annular target, or the
smaller ball. Their densities are correct only on their restricted support.
They are intentionally invalid importance proposals for the full target;
they show that ordinary weight diagnostics cannot establish support.
They are distinct from the earlier defensive-mixture failure controls.
No oracle points or identities enter learned pilots or production.

| Geometry | Seed | Control | ESS | Usual screens | Truth screens |
|---|---:|---|---:|---|---|
| annular_cylinder | 2026092611 | full_target_oracle | 8192.0 | pass | pass |
| annular_cylinder | 2026092612 | full_target_oracle | 8192.0 | pass | pass |
| annular_cylinder | 2026092613 | full_target_oracle | 8192.0 | pass | pass |
| annular_cylinder | 2026092711 | support_hole_oracle | 8192.0 | pass | fail |
| annular_cylinder | 2026092712 | support_hole_oracle | 8192.0 | pass | fail |
| annular_cylinder | 2026092713 | support_hole_oracle | 8192.0 | pass | fail |
| rotated_box | 2026092601 | full_target_oracle | 8192.0 | pass | pass |
| rotated_box | 2026092602 | full_target_oracle | 8192.0 | pass | pass |
| rotated_box | 2026092603 | full_target_oracle | 8192.0 | pass | pass |
| rotated_box | 2026092701 | support_hole_oracle | 8192.0 | pass | fail |
| rotated_box | 2026092702 | support_hole_oracle | 8192.0 | pass | fail |
| rotated_box | 2026092703 | support_hole_oracle | 8192.0 | pass | fail |
| unequal_balls | 2026092621 | full_target_oracle | 8192.0 | pass | pass |
| unequal_balls | 2026092622 | full_target_oracle | 8192.0 | pass | pass |
| unequal_balls | 2026092623 | full_target_oracle | 8192.0 | pass | pass |
| unequal_balls | 2026092721 | support_hole_oracle | 8192.0 | pass | fail |
| unequal_balls | 2026092722 | support_hole_oracle | 8192.0 | pass | fail |
| unequal_balls | 2026092723 | support_hole_oracle | 8192.0 | pass | fail |

## Scope and reproducibility

- all_nine_full_target_oracles_pass: pass.
- all_nine_learned_runs_pass: fail.
- all_nine_support_holes_fail_truth: pass.
- all_nine_support_holes_pass_usual_screens: pass.

Archives retain pilot endpoints and histories, fitted kernels, independent
production points, component identities, densities, scores and decisions.
Reassembly checks source and archive hashes, replays production randomness,
recomputes analytic endpoint scores and reconciles pilot counts. It does not
replay unarchived intermediate pilot trajectories. Numerical assessment
roundoff can vary across platforms; counts and pass/fail decisions remain exact.
Runtime labels describe the archived learned runs; reconstructed analytic
controls execute in the reader's environment rather than the archived runtime.

The controls establish only that the fixed checks can detect specified
support defects and that their accuracy is achievable with analytic oracles.
Passing learned runs cannot exclude a region outside these challenges.
A failed challenge does not alter the historical biological scores or prove
that the biological target shares that geometry. No sampler tuning, pooled
posterior, assay validation, or clinical claim is inferred from this study.
