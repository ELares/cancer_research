# Final-proposal covariance comparison

**Prespecified checks failed. Every planned outcome is retained.**

Only the final independent production proposal changes. Each trio shares an unchanged
29,184-attempt pilot; each arm then makes 8,192 independent production attempts.
Out-of-cube Gaussian draws retain their full mixture densities, have zero target
weight, and remain in the denominator. No boundary clipping or retry is used.

The old three geometries are regression cases already seen during earlier work.
Only shifted_rotated_box is a new prospective geometry, selected with knowledge of
earlier failures. None of these results is independent biological validation.

The unbounded diagonal and correlated arms share centers, multiplicities and marginal
variances; they differ only in off-diagonal covariance. Comparison with the historical
bounded arm also changes boundary handling and kernel marginal scales.

## All proposal outcomes

| Geometry | Seed | Arm | Outside | Accepted | ESS | Max weight | Mass error | Result |
|---|---:|---|---:|---:|---:|---:|---:|---|
| annular_cylinder | 2026092811 | bounded_diagonal | 0 | 729 | 457.2 | 0.00838 | -0.0179 | pass |
| annular_cylinder | 2026092811 | unbounded_diagonal | 879 | 724 | 428.0 | 0.01315 | +0.0784 | pass |
| annular_cylinder | 2026092811 | correlated | 832 | 705 | 387.9 | 0.00945 | -0.0362 | fail |
| annular_cylinder | 2026092812 | bounded_diagonal | 0 | 750 | 539.7 | 0.00579 | +0.0378 | pass |
| annular_cylinder | 2026092812 | unbounded_diagonal | 919 | 651 | 413.5 | 0.01053 | -0.0172 | pass |
| annular_cylinder | 2026092812 | correlated | 966 | 696 | 406.5 | 0.01041 | +0.0453 | pass |
| annular_cylinder | 2026092813 | bounded_diagonal | 0 | 725 | 540.1 | 0.00944 | +0.0402 | pass |
| annular_cylinder | 2026092813 | unbounded_diagonal | 1035 | 590 | 447.2 | 0.00696 | -0.0186 | pass |
| annular_cylinder | 2026092813 | correlated | 1003 | 655 | 418.9 | 0.01286 | +0.0043 | pass |
| rotated_box | 2026092801 | bounded_diagonal | 0 | 365 | 185.8 | 0.03023 | +0.0757 | fail |
| rotated_box | 2026092801 | unbounded_diagonal | 50 | 347 | 200.1 | 0.02206 | -0.0151 | fail |
| rotated_box | 2026092801 | correlated | 36 | 1095 | 502.2 | 0.01405 | +0.0145 | pass |
| rotated_box | 2026092802 | bounded_diagonal | 0 | 573 | 241.8 | 0.02868 | +0.0202 | fail |
| rotated_box | 2026092802 | unbounded_diagonal | 18 | 584 | 235.3 | 0.02604 | +0.1228 | fail |
| rotated_box | 2026092802 | correlated | 9 | 1265 | 456.0 | 0.01728 | -0.0027 | pass |
| rotated_box | 2026092803 | bounded_diagonal | 0 | 580 | 269.6 | 0.01726 | -0.0433 | pass |
| rotated_box | 2026092803 | unbounded_diagonal | 5 | 580 | 261.4 | 0.01875 | -0.1040 | pass |
| rotated_box | 2026092803 | correlated | 5 | 1331 | 546.2 | 0.01773 | -0.0140 | pass |
| shifted_rotated_box | 2026092831 | bounded_diagonal | 0 | 437 | 165.1 | 0.04651 | -0.0049 | fail |
| shifted_rotated_box | 2026092831 | unbounded_diagonal | 456 | 454 | 192.1 | 0.02690 | +0.1557 | fail |
| shifted_rotated_box | 2026092831 | correlated | 361 | 1162 | 532.7 | 0.00933 | -0.0976 | pass |
| shifted_rotated_box | 2026092832 | bounded_diagonal | 0 | 386 | 205.2 | 0.02670 | -0.0805 | fail |
| shifted_rotated_box | 2026092832 | unbounded_diagonal | 307 | 401 | 181.3 | 0.02728 | +0.0625 | fail |
| shifted_rotated_box | 2026092832 | correlated | 244 | 1054 | 581.7 | 0.01189 | -0.0753 | fail |
| shifted_rotated_box | 2026092833 | bounded_diagonal | 0 | 523 | 260.9 | 0.01435 | +0.0008 | fail |
| shifted_rotated_box | 2026092833 | unbounded_diagonal | 419 | 466 | 154.3 | 0.03286 | -0.0122 | fail |
| shifted_rotated_box | 2026092833 | correlated | 328 | 1221 | 535.9 | 0.00973 | -0.0627 | pass |
| unequal_balls | 2026092821 | bounded_diagonal | 0 | 1562 | 974.5 | 0.00700 | -0.0401 | pass |
| unequal_balls | 2026092821 | unbounded_diagonal | 645 | 1465 | 1010.7 | 0.00560 | -0.0130 | pass |
| unequal_balls | 2026092821 | correlated | 633 | 1228 | 778.3 | 0.01187 | -0.0026 | pass |
| unequal_balls | 2026092822 | bounded_diagonal | 0 | 1645 | 1270.2 | 0.00327 | +0.0382 | pass |
| unequal_balls | 2026092822 | unbounded_diagonal | 656 | 1492 | 1177.5 | 0.00320 | +0.0255 | pass |
| unequal_balls | 2026092822 | correlated | 629 | 1209 | 889.8 | 0.00668 | -0.0166 | pass |
| unequal_balls | 2026092823 | bounded_diagonal | 0 | 1388 | 975.6 | 0.00696 | -0.0514 | pass |
| unequal_balls | 2026092823 | unbounded_diagonal | 721 | 1290 | 978.0 | 0.00424 | -0.0320 | pass |
| unequal_balls | 2026092823 | correlated | 692 | 993 | 701.6 | 0.00721 | -0.0326 | pass |

Failed checks:

- annular_cylinder, 2026092811, correlated: moments.
- rotated_box, 2026092801, bounded_diagonal: ess, maximum_weight.
- rotated_box, 2026092801, unbounded_diagonal: maximum_weight.
- rotated_box, 2026092802, bounded_diagonal: maximum_weight.
- rotated_box, 2026092802, unbounded_diagonal: maximum_weight, moments.
- shifted_rotated_box, 2026092831, bounded_diagonal: ess, maximum_weight, moments.
- shifted_rotated_box, 2026092831, unbounded_diagonal: ess, maximum_weight, moments, normalizing_mass.
- shifted_rotated_box, 2026092832, bounded_diagonal: maximum_weight, moments.
- shifted_rotated_box, 2026092832, unbounded_diagonal: ess, maximum_weight, moments.
- shifted_rotated_box, 2026092832, correlated: moments.
- shifted_rotated_box, 2026092833, bounded_diagonal: moments.
- shifted_rotated_box, 2026092833, unbounded_diagonal: ess, maximum_weight, moments.

## Analytic controls

| Geometry | Seed | Kind | Usual screens | Truth screens |
|---|---:|---|---|---|
| annular_cylinder | 2026092911 | full_target_oracle | pass | pass |
| annular_cylinder | 2026092912 | full_target_oracle | pass | pass |
| annular_cylinder | 2026092913 | full_target_oracle | pass | pass |
| annular_cylinder | 2026093011 | support_hole_oracle | pass | fail |
| annular_cylinder | 2026093012 | support_hole_oracle | pass | fail |
| annular_cylinder | 2026093013 | support_hole_oracle | pass | fail |
| rotated_box | 2026092901 | full_target_oracle | pass | pass |
| rotated_box | 2026092902 | full_target_oracle | pass | pass |
| rotated_box | 2026092903 | full_target_oracle | pass | pass |
| rotated_box | 2026093001 | support_hole_oracle | pass | fail |
| rotated_box | 2026093002 | support_hole_oracle | pass | fail |
| rotated_box | 2026093003 | support_hole_oracle | pass | fail |
| shifted_rotated_box | 2026092931 | full_target_oracle | pass | pass |
| shifted_rotated_box | 2026092932 | full_target_oracle | pass | pass |
| shifted_rotated_box | 2026092933 | full_target_oracle | pass | pass |
| shifted_rotated_box | 2026093031 | support_hole_oracle | pass | fail |
| shifted_rotated_box | 2026093032 | support_hole_oracle | pass | fail |
| shifted_rotated_box | 2026093033 | support_hole_oracle | pass | fail |
| unequal_balls | 2026092921 | full_target_oracle | pass | pass |
| unequal_balls | 2026092922 | full_target_oracle | pass | pass |
| unequal_balls | 2026092923 | full_target_oracle | pass | pass |
| unequal_balls | 2026093021 | support_hole_oracle | pass | fail |
| unequal_balls | 2026093022 | support_hole_oracle | pass | fail |
| unequal_balls | 2026093023 | support_hole_oracle | pass | fail |

## Interpretation and replay

- all_correlated_runs_pass: fail.
- all_full_target_oracles_pass: pass.
- all_support_holes_fail_truth: pass.
- all_support_holes_pass_usual: pass.

The all-run rule is a fixed benchmark criterion, not a calibrated confidence statement.
Arm comparisons are descriptive across three pilots per geometry; no superiority test or
pooled posterior is inferred. Normalizer estimates include every attempted draw;
self-normalized region masses and moments retain finite-sample bias.
Finite moment and CDF checks cannot establish a complete distribution or precise tails.

Archives retain shared pilot endpoints/history and every arm’s proposal, points, component
IDs, densities, scores and decisions. Replay verifies source/archive hashes, reconstructs
fits and production randomness, and audits pilot endpoints and counters. Intermediate
pilot paths are not archived or independently replayed. Runtime labels describe archived
learned runs; analytic controls are regenerated in the reader’s environment.
Offline replay verifies source bytes against their recorded hashes, not the Git
history behind the informational commit label. No Git database is needed for replay.

Recorded protocol/source commit: `0b26cbdc9dccd30f7cc1bc7aac25b3b095dde0bc`.
See [the protocol](../../docs/CORRELATED_PROPOSAL_PLAN.md) for all fixed inputs and criteria.

Rebuild without training: `python scripts/proposal_correlated_study.py --render-only`.
