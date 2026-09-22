"""Audit additive RNG seed addresses for the declared 20-block 2D study.

The three treatment arms intentionally retain the engine's existing offsets.
Disjoint block envelopes prevent reuse of a seed argument across these blocks;
they do not prove independence of pseudorandom streams or model outcomes.
This module performs integer arithmetic and deterministic geometry checks only.
It neither samples a generator nor runs the simulator.
"""
import math


BLOCK_COUNT = 20
BASE_SEED = 42
BLOCK_STRIDE = 1 << 32
U64_MAX = (1 << 64) - 1
ARMS = ("Control", "RSL3", "SDT")


def _integer(value, name: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an exact integer >= {minimum}")
    return value


def block_seed(block: int) -> int:
    """Return 42 + block * 2**32 for declared blocks 1 through 20."""
    _integer(block, "block", 1)
    if block > BLOCK_COUNT:
        raise ValueError(f"block must be between 1 and {BLOCK_COUNT}")
    return BASE_SEED + block * BLOCK_STRIDE


def _dimensions(config: dict) -> tuple[int, int, int]:
    rows = _integer(config["grid_rows"], "grid_rows", 1)
    cols = _integer(config["grid_cols"], "grid_cols", 1)
    steps = _integer(config["n_steps"], "n_steps", 1)
    return rows, cols, steps


def rng_namespace(seed: int, config: dict) -> dict:
    """Inclusive envelope for geometry and all potential three-arm seeds.

    Every grid index and step 0..n_steps-1 is included conservatively, even
    when a treatment, dead cell, or immune-eligibility gate skips an RNG call.
    Sparse seed addresses inside the envelope are not asserted to be used.
    """
    _integer(seed, "seed")
    if seed > U64_MAX:
        raise ValueError("seed exceeds u64")
    rows, cols, steps = _dimensions(config)
    rng = config["rng"]
    if rng["scheme"] != "additive_per_cell":
        raise ValueError("unsupported RNG scheme")
    fields = ("treatment_seed_stride", "init_offset", "biochem_offset",
              "biochem_step_stride", "immune_offset", "immune_step_stride")
    offsets = {name: _integer(rng[name], f"rng.{name}") for name in fields}
    arm_offset = (len(ARMS) - 1) * offsets["treatment_seed_stride"]
    max_cell = rows * cols - 1
    maximum = seed + arm_offset + max_cell + max(
        offsets["init_offset"],
        offsets["biochem_offset"] + (steps - 1) * offsets["biochem_step_stride"],
        offsets["immune_offset"] + (steps - 1) * offsets["immune_step_stride"],
    )
    if maximum > U64_MAX:
        raise ValueError("RNG seed namespace would wrap u64")
    return {"minimum": seed, "maximum": maximum}


def _tumor_indices(config: dict) -> set[int]:
    rows, cols, _steps = _dimensions(config)
    size, radius = config["cell_size_um"], config["tumor_radius_um"]
    for name, value in (("cell_size_um", size), ("tumor_radius_um", radius)):
        if (type(value) not in (int, float) or not math.isfinite(value)
                or value <= 0):
            raise ValueError(f"{name} must be finite and positive")
    radius_cells = radius / size
    return {row * cols + col for row in range(rows) for col in range(cols)
            if (row - rows / 2) ** 2 + (col - cols / 2) ** 2 <= radius_cells ** 2}


def audit(config: dict) -> dict:
    """Validate all block envelopes and quantify historical 42/43 address reuse.

    SDT initializes an RNG at every grid index, including stromal cells;
    Control and RSL3 skip that draw. Tumor membership depends on the circle,
    not the geometry RNG. Thus SDT initialization overlap is actual call-site
    overlap, whereas later biochemistry/immune counts are potential overlap
    before the dynamic survival, release and DAMP gates are evaluated.
    """
    roots = [block_seed(block) for block in range(1, BLOCK_COUNT + 1)]
    namespaces = [{"block": block, "seed": seed, **rng_namespace(seed, config)}
                  for block, seed in enumerate(roots, 1)]
    disjoint = all(left["maximum"] < right["minimum"]
                   for left, right in zip(namespaces, namespaces[1:]))
    if not disjoint:
        raise ValueError("declared block RNG namespace envelopes overlap")

    rows, cols, steps = _dimensions(config)
    immune_start = _integer(config["immune_start_step"], "immune_start_step")
    immune_steps = max(steps - immune_start, 0)
    tumor = _tumor_indices(config)
    # 42 + i == 43 + j exactly when i == j + 1. This set calculation
    # includes only indices that are tumor cells in both configurations.
    overlap = sum(index - 1 in tumor for index in tumor)
    biochem_pairs = overlap * steps * len(ARMS)
    immune_pairs = overlap * immune_steps * len(ARMS)
    historical = [rng_namespace(seed, config) for seed in (42, 43)]
    overlap_min = max(item["minimum"] for item in historical)
    overlap_max = min(item["maximum"] for item in historical)
    return {
        "schema_version": 1,
        "block_count": BLOCK_COUNT,
        "base_seed": BASE_SEED,
        "block_stride": BLOCK_STRIDE,
        "roots": roots,
        "namespaces": namespaces,
        "disjoint": disjoint,
        "no_u64_wrap": True,
        "grid_cells": rows * cols,
        "tumor_cells": len(tumor),
        "historical_adjacent": {
            "seeds": [42, 43],
            "envelope_overlap": {
                "minimum": overlap_min, "maximum": overlap_max,
                "integer_addresses": max(overlap_max - overlap_min + 1, 0),
            },
            "sdt_initialization": {
                "all_grid": {"calls_per_run": rows * cols,
                             "shared_seed_addresses": rows * cols - 1},
                "tumor_only": {"calls_per_run": len(tumor),
                               "shared_seed_addresses": overlap},
            },
            "potential_matched_tumor_call_sites": {
                "per_step_per_arm": overlap,
                "biochemistry": {"steps": steps, "arms": len(ARMS),
                                 "matched_call_sites": biochem_pairs},
                "immune": {"steps": immune_steps, "arms": len(ARMS),
                           "matched_call_sites": immune_pairs},
                "total_including_sdt_tumor_initialization": overlap + biochem_pairs + immune_pairs,
            },
        },
        "interpretation": [
            "Namespaces are inclusive envelopes over potential seed arguments, not dense executed streams.",
            "The block envelopes include geometry and all three treatment arms and do not overlap across blocks.",
            "SDT initialization constructs RNGs for every grid cell; Control and RSL3 do not use initialization RNG draws.",
            "Historical 42/43 initialization counts identify identical seed arguments at different cell indices, not identical biological states.",
            "Later call-site counts are potential corresponding-arm/phase/step matches, not observed draws or distinct seed values.",
            "Within-block geometry reuse and treatment/step seed-address relationships are unchanged.",
            "Seed-address reuse alone does not quantify outcome correlation or bias; disjoint envelopes do not prove PRNG independence.",
        ],
    }
