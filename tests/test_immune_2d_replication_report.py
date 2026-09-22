"""Frozen whole-block inference and offline replication-archive regression tests.

All new observations here are constructed from the small accounting fixture;
these tests never invoke Rust or create production simulation observations.
"""

import copy
import gzip
import json
import math
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import immune_2d_replication_report as report
from test_immune_2d_measurement_report import dataset, fixture, population, source_fixture


def block_data(block=1):
    data, _, contract, _ = dataset()
    contract = copy.deepcopy(contract)
    root = 42 + block * 2**32
    data["replicate_block"] = block
    data["config"]["seed"] = root
    for arm_index, row in enumerate(data["conditions"]):
        row["seed"] = root + arm_index * 10_000_000
    return data, contract


def empty_row(row):
    obs = row["measurements"]
    for key in ("ferroptotic_events", "eligible_cells", "immune_kill_events"):
        obs[key] = []
    for step in obs["steps"]:
        step.update({key: 0 for key in step if key != "step"})
        step.update(population([]))
    obs["terminal"] = {key: 0 for key in obs["terminal"]}
    row["result"].update(total_dead=0, ferroptosis_kills=0, immune_kills=0)
    row["final_damp"] = {"total": 0, "peak": 0}


def synthetic_blocks():
    """Distinct, strongly paired arm levels with an asymmetric difference list."""
    differences = [-5, -3, -1, 0, 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50]
    blocks = []
    for block, difference in enumerate(differences, 1):
        summaries = {}
        for arm, value in zip(report.ARMS, (block**3, 1000 * block, 1000 * block + difference)):
            summaries[arm] = {
                report.ENDPOINT: value,
                "damp_ge_kd_opportunities": int(block % 2 == 0),
                "damp_ge_9kd_opportunities": int(block % 4 == 0),
            }
        blocks.append({"block": block, "root_seed": 42 + block * 2**32,
                       "arms": summaries, "sdt_rsl3_immune_kill_ratio": None if block == 1 else 2.})
    return blocks


def write_synthetic_archive(path):
    raw, baseline, sources = source_fixture()
    canonical = gzip.compress(raw, mtime=0)
    sources.update({
        report.PLAN: json.dumps(report.FIXED_PLAN).encode(),
        report.PROTOCOL: b"synthetic prospectively frozen protocol\n",
        report.CANONICAL_MANIFEST: json.dumps({
            "artifacts": {"observations.json.gz": report.sha(canonical)}}).encode(),
        "scripts/immune_seed_blocks.py": b"synthetic seed audit source\n",
        "scripts/immune_2d_replication_report.py": b"synthetic replication source\n",
    })
    for name in report.source_paths():
        sources.setdefault(name, b"synthetic frozen source\n")
    artifacts = {"baseline-summary.json": baseline,
                 "canonical-observations.json.gz": canonical,
                 "sources.tar.gz": report.pack_files(sources),
                 "capture-logs.tar.gz": report.pack_files({
                     **{"baseline.log": b"", "canonical.log": b""},
                     **{f"block-{block:02d}.log": b"" for block in range(1, 21)}})}
    runs = []
    for block in range(1, 21):
        name = f"block-{block:02d}.json.gz"
        data, _ = block_data(block)
        # A whole empty block exercises null conditional summaries in the report.
        if block == 20:
            for row in data["conditions"]:
                empty_row(row)
        artifacts[name] = gzip.compress(json.dumps(data).encode(), mtime=0)
        runs.append({"block": block, "seed": 42 + block * 2**32,
                     "artifact": name, "elapsed_seconds": .125})
    manifest = {
        "schema_version": 1, "study": "immune-2d-activation-replication", "complete": True,
        "source_commit": "3" * 40, "captured_at_utc": "2026-01-01T00:00:00+00:00",
        "platform": "synthetic", "rustc": "rustc 1.96.0", "execution": "serial",
        "python": {"implementation": "CPython", "version": "3.14.0"},
        "binary_sha256": "4" * 64, "runs": runs,
        "baseline_reference": report.measurement.baseline_reference(
            sources[report.measurement.GOLDEN].decode()),
        "sources": {name: report.sha(blob) for name, blob in sources.items()},
        "artifacts": {name: report.sha(blob) for name, blob in artifacts.items()},
    }
    path.mkdir()
    for name, blob in artifacts.items():
        (path / name).write_bytes(blob)
    (path / "manifest.json").write_text(json.dumps(manifest))
    return path


@pytest.fixture
def archive(tmp_path):
    return write_synthetic_archive(tmp_path / "replication archive")


@pytest.fixture(autouse=True)
def forbid_production_capture(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("tests must not start a production simulation")
    monkeypatch.setattr(report, "capture", forbidden)


def cli_paths(monkeypatch, tmp_path, archive, *args):
    monkeypatch.setattr(report, "ARCHIVE", archive)
    monkeypatch.setattr(report.measurement, "ARCHIVE", tmp_path / "canonical")
    monkeypatch.setattr(report, "OUT_JSON", tmp_path / "published.json")
    monkeypatch.setattr(report, "OUT_MD", tmp_path / "published.md")
    monkeypatch.setattr(sys, "argv", ["replication-report", *args])


def replace_artifact(archive, manifest, name, blob):
    (archive / name).write_bytes(blob)
    manifest["artifacts"][name] = report.sha(blob)


def test_checked_in_plan_is_the_frozen_twenty_block_design():
    plan = json.loads((report.ROOT / report.PLAN).read_text())
    report.validate_plan(plan)
    assert plan["block_ids"] == list(range(1, 21))
    assert plan["block_seed_stride"] == 2**32
    assert plan["bootstrap"] == {
        "unit": "whole_three_arm_block", "statistic": "arithmetic_mean_of_paired_differences",
        "resamples": 10_000, "seed": 20260922, "quantiles": [.025, .975], "interpolation": "linear"}


@pytest.mark.parametrize("mutation", ["blocks", "arms", "contrast", "endpoint", "stride", "extra",
                                     "schema_type", "resamples", "rng_seed", "quantiles", "unit"])
def test_frozen_plan_changes_are_rejected(mutation):
    plan = copy.deepcopy(report.FIXED_PLAN)
    if mutation == "blocks":
        plan["block_ids"][-1] = 21
    elif mutation == "arms":
        plan["arms"].reverse()
    elif mutation == "contrast":
        plan["contrast"].reverse()
    elif mutation == "endpoint":
        plan["primary_endpoint"] = "mean_activation"
    elif mutation == "stride":
        plan["block_seed_stride"] = 1
    elif mutation == "extra":
        plan["optional_stopping"] = True
    elif mutation == "schema_type":
        plan["schema_version"] = True
    elif mutation == "resamples":
        plan["bootstrap"]["resamples"] = 1000
    elif mutation == "rng_seed":
        plan["bootstrap"]["seed"] += 1
    elif mutation == "quantiles":
        plan["bootstrap"]["quantiles"] = [.05, .95]
    elif mutation == "unit":
        plan["bootstrap"]["unit"] = "individual_cell"
    with pytest.raises(ValueError, match="frozen replication plan changed"):
        report.validate_plan(plan)


@pytest.mark.parametrize(("p", "expected"), [(0, 0), (.25, 7.5), (.5, 25), (.975, 95.5), (1, 100)])
def test_quantile_uses_linear_interpolation_of_ordered_values(p, expected):
    assert report.quantile([40, 0, 100, 10], p) == pytest.approx(expected)
    assert report.quantile([7], p) == 7


@pytest.mark.parametrize(("values", "p"), [([], .5), ([1], -.01), ([1], 1.01)])
def test_invalid_quantiles_are_rejected(values, p):
    with pytest.raises(ValueError):
        report.quantile(values, p)


def test_paired_bootstrap_matches_independent_reference_and_ignores_common_arm_shifts():
    blocks = synthetic_blocks()
    actual = report.summarize(blocks, copy.deepcopy(report.FIXED_PLAN))
    primary = actual["primary"]
    assert primary["paired_differences"] == [-5, -3, -1, 0, 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50]
    assert primary["mean_paired_difference"] == 8.35
    # Independently computed with a NumPy matrix of frozen randrange indices,
    # row means, and np.quantile(..., method="linear"); these are not zero-width CIs.
    assert primary["percentile_interval"] == pytest.approx([3.55, 14.25])
    assert primary["interpretation"] == "positive"
    shifted = copy.deepcopy(blocks)
    for block in shifted:
        offset = block["block"] ** 4
        for arm in ("RSL3", "SDT"):
            block["arms"][arm][report.ENDPOINT] += offset
        block["arms"]["Control"][report.ENDPOINT] = 10**8
    assert report.summarize(shifted, report.FIXED_PLAN)["primary"] == primary
    assert actual["sdt_rsl3_immune_kill_ratio"] == {
        "n_defined": 19, "undefined_blocks": [1], "median": 2., "minimum": 2., "maximum": 2.}
    assert actual["threshold_positive_blocks"]["Control"]["damp_ge_9kd_opportunities"] == [4, 8, 12, 16, 20]


@pytest.mark.parametrize(("difference", "interpretation"), [(-2, "negative"), (0, "unresolved"), (2, "positive")])
def test_degenerate_bootstrap_interpretation(difference, interpretation):
    blocks = synthetic_blocks()
    for block in blocks:
        block["arms"]["SDT"][report.ENDPOINT] = block["arms"]["RSL3"][report.ENDPOINT] + difference
    primary = report.summarize(blocks, report.FIXED_PLAN)["primary"]
    assert primary["percentile_interval"] == [difference, difference]
    assert primary["interpretation"] == interpretation


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "order", "seed", "nonfinite"])
def test_summaries_require_all_twenty_correctly_identified_blocks(mutation):
    blocks = synthetic_blocks()
    if mutation == "missing":
        blocks.pop()
    elif mutation == "duplicate":
        blocks[-1] = copy.deepcopy(blocks[0])
    elif mutation == "order":
        blocks.reverse()
    elif mutation == "seed":
        blocks[0]["root_seed"] = 42
    elif mutation == "nonfinite":
        blocks[0]["arms"]["SDT"][report.ENDPOINT] = math.inf
    with pytest.raises(ValueError):
        report.summarize(blocks, report.FIXED_PLAN)


@pytest.mark.parametrize("block", [1, 20])
def test_block_endpoint_uses_initial_tumor_census_and_all_eligible_opportunities(block):
    data, contract = block_data(block)
    before = copy.deepcopy(data)
    actual = report.validate_block(data, block, contract)
    expected_activation = .9 + .5 + .02 / 50.02 + 2 / 3
    _, cfg = fixture()
    assert cfg["grid_rows"] * cfg["grid_cols"] == 100
    assert actual["root_seed"] == 42 + block * 2**32
    for arm in report.ARMS:
        assert actual["arms"][arm]["activation_sum"] == pytest.approx(expected_activation)
        assert actual["arms"][arm][report.ENDPOINT] == pytest.approx(expected_activation / 49)
        assert actual["arms"][arm]["mean_activation"] == pytest.approx(expected_activation / 4)
    assert data == before


@pytest.mark.parametrize("mutation", ["block", "block_bool", "schema_bool", "dimension", "simulator",
                                     "config", "root_seed", "root_float", "arm_seed", "arm_seed_float",
                                     "order", "duplicate_arm", "treatment", "scenario", "ph", "census",
                                     "outside_circle", "activation", "terminal", "death_count"])
def test_block_identity_configuration_seeds_and_accounting_are_enforced(mutation):
    data, contract = block_data()
    row = data["conditions"][1]
    if mutation == "block":
        data["replicate_block"] = 2
    elif mutation == "block_bool":
        data["replicate_block"] = True
    elif mutation == "schema_bool":
        data["schema_version"] = True
    elif mutation == "dimension":
        data["dimension"] = 3
    elif mutation == "simulator":
        data["simulator"] = "sim-tme-3d"
    elif mutation == "config":
        data["config"]["params"]["sdt_ros"] += 1
    elif mutation == "root_seed":
        data["config"]["seed"] = 42
    elif mutation == "root_float":
        data["config"]["seed"] = float(data["config"]["seed"])
    elif mutation == "arm_seed":
        row["seed"] += 1
    elif mutation == "arm_seed_float":
        row["seed"] = float(row["seed"])
    elif mutation == "order":
        data["conditions"].reverse()
    elif mutation == "duplicate_arm":
        data["conditions"][-1] = copy.deepcopy(data["conditions"][0])
    elif mutation == "treatment":
        row["result"]["treatment"] = "SDT"
    elif mutation == "scenario":
        row["result"]["o2_condition"] = "normoxic"
    elif mutation == "ph":
        row["result"]["ph_mode"] = "off"
    elif mutation == "census":
        row["result"]["total_tumor"] = 100
    elif mutation == "outside_circle":
        row["measurements"]["ferroptotic_events"][1]["cell_index"] = 0
    elif mutation == "activation":
        row["measurements"]["eligible_cells"][1]["activation_sum"] += .001
    elif mutation == "terminal":
        row["measurements"]["terminal"]["terminal_damp"] += 1
    elif mutation == "death_count":
        row["result"]["total_dead"] += 1
    with pytest.raises(ValueError):
        report.validate_block(data, 1, contract)


def test_zero_exposure_retains_zero_primary_and_null_conditional_outcomes():
    data, contract = block_data()
    for row in data["conditions"]:
        empty_row(row)
    block = report.validate_block(data, 1, contract)
    assert block["sdt_rsl3_immune_kill_ratio"] is None
    for arm in report.ARMS:
        actual = block["arms"][arm]
        assert actual["activation_sum"] == actual[report.ENDPOINT] == 0
        for key in ("mean_activation", "max_eligible_damp", "damp_ge_kd_opportunities_fraction",
                    "kills_per_unique_eligible_cell", "kills_per_eligible_cell_step"):
            assert actual[key] is None
    assert report.distribution([(1, None), (2, None)]) == {
        "n_defined": 0, "undefined_blocks": [1, 2], "median": None, "minimum": None, "maximum": None}


def test_complete_archive_reconstructs_offline_without_changing_any_input(archive, monkeypatch, tmp_path):
    before = {path.name: path.read_bytes() for path in archive.iterdir()}
    # The loader must use the archived contract and sources, not the live checkout.
    monkeypatch.setattr(report, "ROOT", tmp_path / "missing-checkout")
    monkeypatch.setattr(report.measurement, "ROOT", tmp_path / "missing-checkout")
    actual = report.load_archive(archive)
    assert [block["block"] for block in actual["blocks"]] == list(range(1, 21))
    assert sum(len(block["arms"]) for block in actual["blocks"]) == 60
    assert actual["seed_audit"]["disjoint"] is True
    assert actual["seed_audit"]["tumor_cells"] == 49
    assert actual["primary"]["mean_paired_difference"] == 0
    assert actual["primary"]["percentile_interval"] == [0, 0]
    assert actual["sdt_rsl3_immune_kill_ratio"]["undefined_blocks"] == [20]
    assert actual["arm_distributions"]["RSL3"]["mean_activation"]["n_defined"] == 19
    assert actual["arm_distributions"]["RSL3"][report.ENDPOINT]["n_defined"] == 20
    assert {path.name: path.read_bytes() for path in archive.iterdir()} == before


@pytest.mark.parametrize("mutation", ["artifact", "source_hash", "baseline_parity", "canonical_parity",
                                     "baseline_provenance", "missing", "extra", "duplicate_run", "missing_run",
                                     "run_seed", "run_artifact", "negative_duration", "incomplete", "parallel",
                                     "block_config", "block_ledger", "frozen_plan",
                                     "missing_protocol_source", "missing_rust_source"])
def test_invalid_archives_preserve_both_published_outputs(archive, tmp_path, monkeypatch, mutation):
    manifest_path = archive / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if mutation == "artifact":
        path = archive / "block-20.json.gz"
        path.write_bytes(path.read_bytes() + b"tampered")
    elif mutation == "source_hash":
        manifest["sources"][report.PLAN] = "0" * 64
    elif mutation == "baseline_parity":
        name = "baseline-summary.json"
        replace_artifact(archive, manifest, name, (archive / name).read_bytes() + b"\n")
    elif mutation == "canonical_parity":
        name = "canonical-observations.json.gz"
        # Valid, hash-consistent JSON still cannot replace the frozen canonical bytes.
        raw = gzip.decompress((archive / name).read_bytes())
        replace_artifact(archive, manifest, name, gzip.compress(raw, mtime=1))
    elif mutation == "baseline_provenance":
        manifest["baseline_reference"]["source_commit"] = "0" * 40
    elif mutation == "missing":
        (archive / "block-20.json.gz").unlink()
    elif mutation == "extra":
        (archive / "block-21.json.gz").write_bytes(b"extra")
    elif mutation == "duplicate_run":
        manifest["runs"][-1] = copy.deepcopy(manifest["runs"][0])
    elif mutation == "missing_run":
        manifest["runs"].pop()
    elif mutation == "run_seed":
        manifest["runs"][-1]["seed"] += 1
    elif mutation == "run_artifact":
        manifest["runs"][-1]["artifact"] = "block-01.json.gz"
    elif mutation == "negative_duration":
        manifest["runs"][-1]["elapsed_seconds"] = -1
    elif mutation == "incomplete":
        manifest["complete"] = False
    elif mutation == "parallel":
        manifest["execution"] = "parallel"
    elif mutation in ("block_config", "block_ledger"):
        name = "block-20.json.gz"
        data = json.loads(gzip.decompress((archive / name).read_bytes()))
        if mutation == "block_config":
            data["config"]["params"]["sdt_ros"] += 1
        else:
            data["conditions"][0]["result"]["total_dead"] += 1
        replace_artifact(archive, manifest, name, gzip.compress(json.dumps(data).encode(), mtime=0))
    elif mutation in ("frozen_plan", "missing_protocol_source", "missing_rust_source"):
        sources = report.unpack_sources((archive / "sources.tar.gz").read_bytes(), manifest["sources"])
        if mutation == "frozen_plan":
            plan = json.loads(sources[report.PLAN])
            plan["bootstrap"]["resamples"] = 99
            sources[report.PLAN] = json.dumps(plan).encode()
        elif mutation == "missing_protocol_source":
            del sources[report.PROTOCOL]
        else:
            del sources["simulations/sim-tme/src/immune_measurements.rs"]
        manifest["sources"] = {name: report.sha(blob) for name, blob in sources.items()}
        replace_artifact(archive, manifest, "sources.tar.gz", report.pack_files(sources))
    manifest_path.write_text(json.dumps(manifest))
    cli_paths(monkeypatch, tmp_path, archive)
    report.OUT_JSON.write_bytes(b"previous JSON\n")
    report.OUT_MD.write_bytes(b"previous Markdown\n")
    before = {path.name: path.read_bytes() for path in archive.iterdir()}
    with pytest.raises(ValueError):
        report.main()
    assert report.OUT_JSON.read_bytes() == b"previous JSON\n"
    assert report.OUT_MD.read_bytes() == b"previous Markdown\n"
    assert {path.name: path.read_bytes() for path in archive.iterdir()} == before


@pytest.mark.parametrize("alias", ["child", "symlink", "hardlink", "canonical_child",
                                  "canonical_hardlink", "outputs_same", "outputs_hardlink"])
def test_path_aliases_fail_before_capture_or_output_writes(archive, tmp_path, monkeypatch, alias):
    cli_paths(monkeypatch, tmp_path, archive, "--capture")
    canonical = report.measurement.ARCHIVE
    canonical.mkdir()
    (canonical / "observations.json.gz").write_bytes(b"immutable canonical")
    if alias == "child":
        monkeypatch.setattr(report, "OUT_MD", archive / "report.md")
    elif alias == "symlink":
        report.OUT_MD.symlink_to(archive / "manifest.json")
    elif alias == "hardlink":
        report.OUT_MD.hardlink_to(archive / "block-01.json.gz")
    elif alias == "canonical_child":
        monkeypatch.setattr(report, "OUT_JSON", canonical / "report.json")
    elif alias == "canonical_hardlink":
        report.OUT_JSON.hardlink_to(canonical / "observations.json.gz")
    elif alias == "outputs_same":
        monkeypatch.setattr(report, "OUT_MD", report.OUT_JSON)
    elif alias == "outputs_hardlink":
        report.OUT_JSON.write_bytes(b"previous output")
        report.OUT_MD.hardlink_to(report.OUT_JSON)
    before = {path.name: path.read_bytes() for path in archive.iterdir()}
    with pytest.raises(ValueError):
        report.main()
    assert {path.name: path.read_bytes() for path in archive.iterdir()} == before
    assert (canonical / "observations.json.gz").read_bytes() == b"immutable canonical"


def test_successful_cli_rebuilds_stale_outputs_and_render_only_preserves_json(archive, tmp_path, monkeypatch):
    cli_paths(monkeypatch, tmp_path, archive)
    report.OUT_JSON.write_text('{"stale": true}\n')
    report.OUT_MD.write_text("stale Markdown\n")
    report.main()
    actual = json.loads(report.OUT_JSON.read_text())
    assert actual == report.load_archive(archive)
    assert report.OUT_MD.read_text() == report.render(actual)
    json_bytes = report.OUT_JSON.read_bytes()
    report.OUT_MD.write_text("stale again\n")
    monkeypatch.setattr(sys, "argv", ["replication-report", "--render-only"])
    monkeypatch.setattr(report, "load_archive", lambda *args: pytest.fail("render-only read the archive"))
    report.main()
    assert report.OUT_JSON.read_bytes() == json_bytes
    assert report.OUT_MD.read_text() == report.render(actual)


def test_render_failure_preserves_both_reports(archive, tmp_path, monkeypatch):
    cli_paths(monkeypatch, tmp_path, archive)
    report.OUT_JSON.write_bytes(b"keep JSON\n")
    report.OUT_MD.write_bytes(b"keep Markdown\n")
    def failed_render(*args):
        raise ValueError("synthetic rendering failure")
    monkeypatch.setattr(report, "render", failed_render)
    with pytest.raises(ValueError, match="synthetic rendering failure"):
        report.main()
    assert report.OUT_JSON.read_bytes() == b"keep JSON\n"
    assert report.OUT_MD.read_bytes() == b"keep Markdown\n"


def test_render_orders_blocks_and_mapping_fields_deterministically(archive):
    original = report.load_archive(archive)
    reordered = copy.deepcopy(original)
    reordered["blocks"].reverse()
    for block in reordered["blocks"]:
        block["arms"] = dict(reversed(list(block["arms"].items())))
    reordered["arm_distributions"] = {
        arm: dict(reversed(list(fields.items())))
        for arm, fields in reversed(list(reordered["arm_distributions"].items()))}
    assert report.render(reordered) == report.render(original)


def test_capture_and_render_only_are_mutually_exclusive(tmp_path, monkeypatch):
    cli_paths(monkeypatch, tmp_path, tmp_path / "absent-archive", "--capture", "--render-only")
    with pytest.raises(SystemExit) as error:
        report.main()
    assert error.value.code == 2


@pytest.mark.skipif(not report.ARCHIVE.exists(), reason="prospectively frozen gate; production capture not yet present")
def test_committed_archive_reconstructs_exactly_to_published_json_and_markdown():
    reconstructed = report.load_archive(report.ARCHIVE)
    assert reconstructed == json.loads(report.OUT_JSON.read_text())
    assert report.render(reconstructed).encode() == report.OUT_MD.read_bytes()
