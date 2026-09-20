# Proposal development, 19 September 2026

Eighteen development trials selected `mixed_gap32` for a separate prospective validation. These results do not validate the selected method and do not authorize biological runs. All attempted configurations, including failures, are retained.

The first declaration fixed four configurations before running any of their 12 trials. After all 12 were inspected, a second declaration fixed two follow-up configurations before running their six trials. Both declarations used the same paired development seeds `2026092301`, `2026092302`, and `2026092303`. No reserved validation or biological seeds were used.

Every trial used the separated-box fixture in seven dimensions, the unchanged final distance threshold 1, two islands of 256 particles, 14 levels, four moves per level, and 8,192 separate production attempts. The pilot budget was 29,184 attempted points; independence-MH pre-rejection can avoid some score evaluations. Analytic evaluations used no biological simulator calls. Production weights were the inverse complete frozen-mixture density, with rejected attempts retained in the normalizer denominator.

| Configuration | Seed suffix | Pass | ESS | Estimated A mass | Relative normalizer error | Failed checks |
|---|---:|:---:|---:|---:|---:|---|
| independent_k32 | 2301 | no | 72.39 | 0.1970 | -6.42% | ess, maximum_weight, relative_mcse |
| independent_k32 | 2302 | no | 40.03 | 0.3114 | +9.26% | ess, maximum_weight, relative_mcse |
| independent_k32 | 2303 | no | 93.59 | 0.0989 | -16.19% | ess, maximum_weight, relative_mcse, normalizing_mass, mode_masses, moments |
| independent_k8 | 2301 | no | 270.05 | 0.2284 | +4.42% | maximum_weight |
| independent_k8 | 2302 | no | 1.66 | 0.7754 | +222.43% | ess, maximum_weight, relative_mcse, normalizing_mass, mode_masses, moments |
| independent_k8 | 2303 | no | 1446.99 | 0.0000 | -23.64% | normalizing_mass, every_region_observed, mode_masses, moments |
| mixed_k32 | 2301 | yes | 848.71 | 0.2437 | -4.47% | none |
| mixed_k32 | 2302 | yes | 841.04 | 0.2514 | -5.30% | none |
| mixed_k32 | 2303 | no | 845.12 | 0.0000 | -27.35% | normalizing_mass, every_region_observed, mode_masses, moments |
| mixed_k8 | 2301 | no | 1224.36 | 0.0000 | -24.45% | normalizing_mass, every_region_observed, mode_masses, moments |
| mixed_k8 | 2302 | no | 2.20 | 0.7325 | +167.66% | ess, maximum_weight, relative_mcse, normalizing_mass, mode_masses, moments |
| mixed_k8 | 2303 | no | 412.45 | 0.2418 | -1.26% | maximum_weight |
| mixed_gap32 | 2301 | yes | 972.90 | 0.2519 | +1.40% | none |
| mixed_gap32 | 2302 | yes | 842.57 | 0.2641 | -0.12% | none |
| mixed_gap32 | 2303 | yes | 757.59 | 0.2436 | +4.18% | none |
| staged_gap32 | 2301 | yes | 960.87 | 0.2557 | +4.47% | none |
| staged_gap32 | 2302 | yes | 979.01 | 0.2827 | -3.55% | none |
| staged_gap32 | 2303 | yes | 941.58 | 0.2600 | +6.01% | none |

The known conditional A mass is 0.25 and the total accepted prior mass is 1/4,096. The unchanged checks require at least 20 accepted attempts, ESS at least 200, maximum normalized weight at most 0.02, relative normalizer MCSE at most 0.1, relative normalizer error at most 0.15, mode-mass error at most 0.075, moment error at most 0.05, observations in both regions, and both pilots reaching the final threshold.

`independent_k32` reproduces the original move with 32 unique neighbors; `independent_k8` changes that limit to eight. Both keep the original scale multiplier 1.5 and scale bounds [0.005, 0.5]. Duplicate particles retain their mixture multiplicity but do not count as unique neighbors. All fitted proposals contain 20% uniform mass and equally weighted bounded Gaussian product kernels.

`mixed_k32` and `mixed_k8` choose a local move with probability 0.5 for each particle each sweep. The local Gaussian step is a single shared diagonal vector `clip(0.5 * median(all component scales), 0.0025, 0.25)`, frozen throughout the level. Coordinates are reflected repeatedly using modulo 2. Its symmetric transition needs only the hard-target acceptance test; global moves retain the independence-MH ratio `q(old)/q(new)`. Random-call order is global candidate draws for all particles, the local mask, normal draws for selected local particles, then MH uniforms for all particles. Unused random draws are not extra model evaluations.

`mixed_gap32` additionally sorts up to 32 unique neighbors, including the center, in globally standardized Euclidean distance. It cuts before the first adjacent positive-distance jump strictly greater than threefold, ignoring the zero self-distance. It uses widths 0.04 for a singleton neighborhood; otherwise it applies the unchanged multiplier and bounds. `staged_gap32` uses the same gap rule but performs exactly two local sweeps with the prelocal shared step, refits the global proposal once, and performs two independence sweeps with that proposal frozen.

Both gap configurations passed all three paired development trials. `mixed_gap32` was selected because it adds less scheduling and refitting complexity; three seeds do not establish a performance difference. The failed `mixed_k32` trial retained only 20 distinct A centers, fewer than its 32-neighbor requirement, and its A kernels borrowed distant B centers. The resulting median A coordinate widths reached approximately 0.41–0.46, while passing trials had widths around 0.08–0.09. This documents the geometry problem addressed by the gap rule. It does not prove the new heuristic discovers every region or preserves its mass in other targets.

A fresh-process check against `scripts/resample_move_local.py` reproduced the three selected development runs exactly: final pilot clouds and distances, all original history fields, proposal parameters, production points/component IDs/densities, and assessment dictionaries. `candidate_verification.json` identifies the checked source hash. Four additional scratch checks covered reflected-kernel symmetry, minority-gap widths, singleton fallback, and conditional-uniform moments; their separate seed 2026092399 was not a learned-proposal trial.

A separate final check against numerical-source freeze `91ca37aae10e2f22c19d0c6c4a2923bc4220b61c` again reproduced all three selected development runs exactly. `final_candidate_verification.json` records that commit and the final candidate source hash; `verify_final_candidate.py.txt` preserves a path-normalized snapshot of the verification script. The earlier verification file and its hashes remain unchanged. This confirms implementation agreement with the selected development algorithm, not independent validation.

The two gzip archives preserve the original JSON bytes and all 18 pilot/proposal/result records; gzip time is zero and its filename field is empty. All five published `.txt` script snapshots use path-only normalization: absolute checkout and scratch directory prefixes are replaced by the explicit `CHECKOUT_ROOT` and `DEVELOPMENT_DIR` placeholders. No other script bytes changed. Exact original scripts remain outside the repository. For each script, `metadata.json` retains their `original_sha256` and separately records the published snapshot hash as `sha256`, along with replacement counts. Metadata paths use logical `DEVELOPMENT_DIR` prefixes. Neither gzip archive nor any numerical result was changed by this normalization. The original repository source hashes, both candidate verification records, and scratch mathematical checks remain preserved; original source hashes alone do not describe the in-memory substitutions.

Production can be replayed without executing any historical scratch code. From the repository root, use the bounded proposal stored in each record and the second child of its seed; the first child was reserved for its pilot:

```python
import gzip, json, sys
import numpy as np
sys.path.insert(0, "scripts")
from bounded_proposal import BoundedGaussianMixture
from proposal_synthetic_validation import summarize

folder = "analysis/calibration/proposal-development-20260919"
for name in ("development_results", "followup_results"):
    with gzip.open(f"{folder}/{name}.json.gz", "rt") as handle:
        records = json.load(handle)["results"]
    for record in records:
        proposal = BoundedGaussianMixture.from_dict(record["proposal"])
        stream = np.random.SeedSequence(record["seed"]).spawn(2)[1]
        points, ids = proposal.sample(np.random.default_rng(stream), 8192)
        rebuilt = summarize(points, proposal.log_density(points), record["fixture"])
        rebuilt["pilot_reached_epsilon"] = record["pilot"]["all_islands_reached_epsilon"]
        rebuilt["passed"] &= rebuilt["pilot_reached_epsilon"]
        # Compare decisions exactly and floats with a 1e-12 tolerance across platforms.
```

Production replay verifies weighting and decisions conditional on each archived proposal. It does not independently replay pilot adaptation. Reconstructing a historical pilot also requires its archived scratch algorithm, original source files, and recorded random-call order; bind the explicit path placeholders to local directories before using a normalized snapshot. These `.txt` files are provenance records and are not executed by artifact tests. Paired development results are not independent validation replicates; pilot particles are not independent posterior observations. Fresh, predeclared validation remains required before biological sampling.
