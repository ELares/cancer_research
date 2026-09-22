"""Seed arithmetic checks need no simulator execution or stochastic outcomes."""
import copy
import importlib.util
import json
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "immune_seed_blocks", REPO / "scripts/immune_seed_blocks.py")
seeds = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(seeds)


@pytest.fixture
def config():
    return json.loads((REPO / "scripts/immune_2d_measurement_config_v1.json").read_text())["config"]


@pytest.fixture
def small():
    return {
        "grid_rows": 4, "grid_cols": 4, "cell_size_um": 1.0,
        "tumor_radius_um": 1.8, "n_steps": 3, "immune_start_step": 1,
        "rng": {
            "scheme": "additive_per_cell", "treatment_seed_stride": 300,
            "init_offset": 0, "biochem_offset": 16, "biochem_step_stride": 16,
            "immune_offset": 100, "immune_step_stride": 20,
        },
    }


def test_declared_roots_use_one_based_blocks_and_exact_stride():
    assert seeds.block_seed(1) == 4_294_967_338
    assert seeds.block_seed(20) == 85_899_345_962
    assert all(seeds.block_seed(i + 1) - seeds.block_seed(i) == 4_294_967_296
               for i in range(1, 20))


@pytest.mark.parametrize("block", [0, -1, 21, True, False, 1.0, "1", None])
def test_block_identifier_must_be_an_exact_declared_integer(block):
    with pytest.raises(ValueError):
        seeds.block_seed(block)


def test_small_envelope_contains_geometry_and_each_explicit_seed_site(small):
    # Largest address: seed 7 + SDT 600 + immune 100 + last step 40 + index 15.
    assert seeds.rng_namespace(7, small) == {"minimum": 7, "maximum": 762}
    addresses = {7}
    for arm in (0, 300, 600):
        for index in range(16):
            addresses.add(7 + arm + index)
            for step in range(3):
                addresses.add(7 + arm + 16 + index + step * 16)
                addresses.add(7 + arm + 100 + index + step * 20)
    namespace = seeds.rng_namespace(7, small)
    assert min(addresses) == namespace["minimum"]
    assert max(addresses) == namespace["maximum"]


def test_small_circle_reuse_is_actual_initialization_and_potential_later_calls(small):
    report = seeds.audit(small)
    # The circle contains rows/columns 1..3: nine cells, with six consecutive
    # flat-index pairs. A square-envelope shortcut would incorrectly give 15.
    assert report["tumor_cells"] == 9
    old = report["historical_adjacent"]
    assert old["sdt_initialization"]["all_grid"] == {
        "calls_per_run": 16, "shared_seed_addresses": 15,
    }
    assert old["sdt_initialization"]["tumor_only"] == {
        "calls_per_run": 9, "shared_seed_addresses": 6,
    }
    potential = old["potential_matched_tumor_call_sites"]
    assert potential["biochemistry"]["matched_call_sites"] == 6 * 3 * 3
    assert potential["immune"]["matched_call_sites"] == 6 * 2 * 3
    assert potential["total_including_sdt_tumor_initialization"] == 96
    assert old["envelope_overlap"] == {
        "minimum": 43, "maximum": 797, "integer_addresses": 755,
    }


def test_frozen_configuration_counts_and_all_twenty_namespace_envelopes(config):
    original = copy.deepcopy(config)
    report = seeds.audit(config)
    assert config == original
    assert report["grid_cells"] == 250_000
    assert report["tumor_cells"] == 159_033
    assert report["roots"] == [42 + i * (1 << 32) for i in range(1, 21)]
    assert report["disjoint"] is report["no_u64_wrap"] is True
    assert len(report["namespaces"]) == 20
    assert report["namespaces"][-1]["maximum"] == 87_177_595_961
    for index, left in enumerate(report["namespaces"]):
        assert left["maximum"] - left["minimum"] == 1_278_249_999
        assert all(left["maximum"] < right["minimum"]
                   for right in report["namespaces"][index + 1:])
    old = report["historical_adjacent"]
    assert old["sdt_initialization"]["all_grid"]["shared_seed_addresses"] == 249_999
    assert old["sdt_initialization"]["tumor_only"]["shared_seed_addresses"] == 158_582
    assert old["potential_matched_tumor_call_sites"] == {
        "per_step_per_arm": 158_582,
        "biochemistry": {"steps": 180, "arms": 3, "matched_call_sites": 85_634_280},
        "immune": {"steps": 120, "arms": 3, "matched_call_sites": 57_089_520},
        "total_including_sdt_tumor_initialization": 142_882_382,
    }


def test_namespace_allows_exact_u64_boundary_and_rejects_wrapping(small):
    last_seed = seeds.U64_MAX - 755
    assert seeds.rng_namespace(last_seed, small)["maximum"] == seeds.U64_MAX
    with pytest.raises(ValueError, match="wrap"):
        seeds.rng_namespace(last_seed + 1, small)


@pytest.mark.parametrize("seed", [-1, True, 1.0, "42", 1 << 64])
def test_namespace_rejects_non_u64_seed_arguments(small, seed):
    with pytest.raises(ValueError):
        seeds.rng_namespace(seed, small)


def test_drift_that_overlaps_declared_blocks_is_rejected_before_geometry(config):
    config["n_steps"] = 3_000
    with pytest.raises(ValueError, match="overlap"):
        seeds.audit(config)


def test_no_immune_steps_does_not_create_potential_calls(small):
    small["immune_start_step"] = 10
    assert seeds.audit(small)["historical_adjacent"]["potential_matched_tumor_call_sites"]["immune"] == {
        "steps": 0, "arms": 3, "matched_call_sites": 0,
    }


@pytest.mark.parametrize("section,key,value", [
    (None, "grid_rows", True), (None, "grid_cols", 0), (None, "n_steps", 1.5),
    ("rng", "scheme", "splitmix64"), ("rng", "init_offset", -1),
    ("rng", "immune_step_stride", 2_000_000.0),
])
def test_namespace_rejects_unsupported_or_noninteger_address_contract(small, section, key, value):
    (small if section is None else small[section])[key] = value
    with pytest.raises(ValueError):
        seeds.rng_namespace(42, small)


@pytest.mark.parametrize("key,value", [("cell_size_um", 0), ("cell_size_um", True),
                                       ("tumor_radius_um", float("nan"))])
def test_nonphysical_geometry_is_not_counted_as_a_missing_tumor(small, key, value):
    small[key] = value
    with pytest.raises(ValueError):
        seeds.audit(small)
