"""Adversarial ledger and publication checks without production observations.

The synthetic fixture deliberately has zero receptor responses despite positive
injected field mass. It tests accounting, not physical transport or experiment
outcomes. No test runs the production controlled comparison.
"""

import copy
import gzip
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import immune_2d_controlled_report as report


def frozen_plan():
    return json.loads((report.ROOT / report.PLAN).read_text())


def zero_ledger():
    """Full declared masks/conditions with analytic mass and zero probe ledgers."""
    plan = frozen_plan()
    masks = report.expected_masks(plan)
    reduced = set(masks["available_recipient_indices_960"])
    conditions = []
    for design in plan["conditions"]:
        mass = 0.0
        steps = []
        for step in range(plan["config"]["n_steps"]):
            injection = design["total_damp"] if step == design["release_step"] else 0.0
            mass = (mass + injection) * (1 - plan["config"]["damp_clearance_rate"])
            steps.append({
                "step": step, "injected_damp": injection, "field_mass": mass,
                "available_cells": design["recipient_count"], "eligible_cells": 0,
                "activation_sum": 0.0, "local_damp_sum": 0.0, "max_eligible_damp": None,
                "damp_ge_kd_opportunities": 0, "damp_ge_9kd_opportunities": 0,
                "field_fingerprint": "0" * 16,
            })
        recipients = [{
            "cell_index": index,
            "available": design["recipient_count"] == 3840 or index in reduced,
            "opportunities": 0, "activation_sum": 0.0, "local_damp_sum": 0.0,
            "max_eligible_damp": None, "damp_ge_kd_opportunities": 0,
            "damp_ge_9kd_opportunities": 0,
        } for index in masks["reference_recipient_indices"]]
        conditions.append({**design, "steps": steps, "recipients": recipients,
                           "field_fingerprint": report.trajectory_fingerprint(steps),
                           "summary": {
                               "activation_sum": 0.0,
                               "activation_sum_per_reference_recipient": 0.0,
                               "opportunities": 0, "unique_eligible_cells": 0,
                               "max_eligible_damp": None,
                               "damp_ge_kd_opportunities": 0,
                               "damp_ge_9kd_opportunities": 0,
                               "final_field_mass": mass,
                           }})
    return {"schema_version": 1, "study": plan["study"], "simulator": "sim-tme",
            "dimension": 2, "config": copy.deepcopy(plan["config"]),
            "mask_config": copy.deepcopy(plan["masks"]), "units": copy.deepcopy(report.UNITS),
            "field_fingerprint_definition": report.FINGERPRINT_DEFINITION,
            "masks": masks, "conditions": conditions}


@pytest.fixture(scope="module")
def ledger_template():
    return zero_ledger()


@pytest.fixture
def ledger(ledger_template):
    return copy.deepcopy(ledger_template)


def add_opportunity(data, *, damp=50.0, step=60):
    """Add one coherent synthetic cell-step to both recipient-availability arms."""
    index = data["masks"]["available_recipient_indices_960"][0]
    activation = damp / (damp + 50.0)
    for row in data["conditions"][:2]:
        counts = {"damp_ge_kd_opportunities": int(damp >= 50),
                  "damp_ge_9kd_opportunities": int(damp >= 450)}
        cell = next(c for c in row["recipients"] if c["cell_index"] == index)
        cell.update(opportunities=1, activation_sum=activation, local_damp_sum=damp,
                    max_eligible_damp=damp, **counts)
        row["steps"][step].update(eligible_cells=1, activation_sum=activation,
                                  max_eligible_damp=damp, local_damp_sum=damp, **counts)
        row["summary"].update(activation_sum=activation,
                              activation_sum_per_reference_recipient=activation / 3840,
                              opportunities=1, unique_eligible_cells=1,
                              max_eligible_damp=damp, **counts)


def test_plan_preserves_factorial_and_matched_quantities():
    plan = frozen_plan()
    report.validate_plan(plan)
    assert len(plan["conditions"]) == 18
    assert len(plan["contrasts"]) == 36
    factorial = plan["conditions"][:16]
    assert {(r["source_count"], r["total_damp"], r["release_step"], r["recipient_count"])
            for r in factorial} == {
                (n, m, t, r) for n in (64, 256) for m in (1280., 5120.)
                for t in (30, 60) for r in (3840, 960)}
    by_id = {row["condition_id"]: row for row in plan["conditions"]}
    for contrast in plan["contrasts"]:
        left, right = (by_id[contrast[k]] for k in ("left", "right"))
        if contrast["kind"] == "equal_total_mass":
            assert left["total_damp"] == right["total_damp"]
        elif contrast["kind"] == "equal_per_source_amount":
            assert left["per_source_damp"] == right["per_source_damp"]


@pytest.mark.parametrize("mutation", ["schema_type", "budget", "step", "denominator", "contrast"])
def test_plan_changes_require_a_new_freeze(mutation):
    plan = frozen_plan()
    if mutation == "schema_type":
        plan["schema_version"] = True
    elif mutation == "budget":
        plan["conditions"][0]["total_damp"] += 1
    elif mutation == "step":
        plan["config"]["immune_start_step"] -= 1
    elif mutation == "denominator":
        plan["reference_recipient_count"] = 960
    else:
        plan["contrasts"][0]["left"], plan["contrasts"][0]["right"] = (
            plan["contrasts"][0]["right"], plan["contrasts"][0]["left"])
    with pytest.raises(ValueError):
        report.validate_plan(plan)


def test_masks_are_nested_disjoint_and_fixed_independently_of_conditions():
    masks = report.expected_masks(frozen_plan())
    source_64 = masks["source_indices_64"]
    reserved = masks["reserved_source_indices"]
    recipients = masks["reference_recipient_indices"]
    reduced = masks["available_recipient_indices_960"]
    assert len(source_64) == 64 and len(reserved) == 256
    assert len(recipients) == 3840 and len(reduced) == 960
    assert reserved == masks["source_indices_256"]
    assert set(source_64) < set(reserved)
    assert set(reduced) < set(recipients)
    assert not set(reserved) & set(recipients)
    assert sorted(reserved + recipients) == [r * 500 + c for r in range(218, 282)
                                             for c in range(218, 282)]
    assert source_64 == [(218 + 8 * r + (2 if r % 2 == 0 else 5)) * 500
                         + 218 + 8 * c + (2 if c % 2 == 0 else 5)
                         for r in range(8) for c in range(8)]
    assert all((i // 500 - 218) % 2 == 0 and (i % 500 - 218) % 2 == 0 for i in reduced)


def test_zero_probe_ledger_reconciles_every_condition_and_keeps_undefined_means(ledger):
    before = copy.deepcopy(ledger)
    summaries = report.validate_observations(ledger, frozen_plan())
    assert len(summaries) == 18
    assert all(row["activation_sum_per_reference_recipient"] == 0 for row in summaries)
    assert all(row["conditional_mean_activation"] is None for row in summaries)
    assert all(row["max_eligible_damp"] is None for row in summaries)
    assert ledger == before
    result = report.summarize(summaries, frozen_plan())
    assert len(result["conditions"]) == 18
    assert len(result["contrasts"]) == 36


@pytest.mark.parametrize("damp", [.01, 50.0, 450.0])
def test_one_opportunity_uses_the_fixed_reference_denominator_and_thresholds(ledger, damp):
    add_opportunity(ledger, damp=damp)
    summaries = report.validate_observations(ledger, frozen_plan())
    assert summaries[0]["activation_sum_per_reference_recipient"] == pytest.approx(damp / (damp + 50) / 3840)
    assert summaries[0]["activation_sum_per_reference_recipient"] == summaries[1]["activation_sum_per_reference_recipient"]
    assert summaries[0]["damp_ge_kd_opportunities"] == int(damp >= 50)
    assert summaries[0]["damp_ge_9kd_opportunities"] == int(damp >= 450)


def test_every_declared_contrast_uses_left_minus_right_without_selective_reporting(ledger):
    expected = {}
    for pair in range(8):
        damp = (pair + 1) / 10
        selected = ledger["conditions"][pair * 2:pair * 2 + 2]
        add_opportunity({**ledger, "conditions": selected}, damp=damp)
        expected.update({row["condition_id"]: damp / (damp + 50) / 3840 for row in selected})
    expected.update({row["condition_id"]: 0. for row in ledger["conditions"][-2:]})
    plan = frozen_plan()
    result = report.summarize(report.validate_observations(ledger, plan), plan)
    assert [c["contrast_id"] for c in result["contrasts"]] == [c["contrast_id"] for c in plan["contrasts"]]
    for comparison in result["contrasts"]:
        assert comparison["difference"] == pytest.approx(
            expected[comparison["left"]] - expected[comparison["right"]])


@pytest.mark.parametrize("mutation", ["missing_condition", "duplicate_condition", "condition_order",
                                     "missing_zero_cell", "cell_order", "duplicate_cell",
                                     "missing_step", "step_order", "mask_coordinate"])
def test_exact_complete_ordered_coverage_is_required(ledger, mutation):
    first = ledger["conditions"][0]
    if mutation == "missing_condition":
        ledger["conditions"].pop()
    elif mutation == "duplicate_condition":
        ledger["conditions"][-1] = copy.deepcopy(first)
    elif mutation == "condition_order":
        ledger["conditions"].reverse()
    elif mutation == "missing_zero_cell":
        first["recipients"].pop()
    elif mutation == "cell_order":
        first["recipients"].reverse()
    elif mutation == "duplicate_cell":
        first["recipients"][-1] = copy.deepcopy(first["recipients"][0])
    elif mutation == "missing_step":
        first["steps"].pop()
    elif mutation == "step_order":
        first["steps"].reverse()
    else:
        ledger["masks"]["source_indices_64"][0] += 1
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


@pytest.mark.parametrize(("field", "value"), [
    ("field_mass", math.nan), ("field_mass", math.inf), ("field_mass", -1.0),
    ("activation_sum", math.nan), ("activation_sum", -1.0),
    ("eligible_cells", True), ("eligible_cells", 0.5),
    ("damp_ge_kd_opportunities", -1), ("max_eligible_damp", math.inf),
])
def test_invalid_step_numeric_values_are_rejected(ledger, field, value):
    ledger["conditions"][0]["steps"][60][field] = value
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


@pytest.mark.parametrize("mutation", ["mass", "injection_amount", "injection_time", "final_mass"])
def test_injection_and_clearance_mass_ledger_are_checked(ledger, mutation):
    row = ledger["conditions"][0]
    if mutation == "mass":
        row["steps"][90]["field_mass"] += 1
    elif mutation == "injection_amount":
        row["steps"][30]["injected_damp"] += 1
    elif mutation == "injection_time":
        row["steps"][29]["injected_damp"], row["steps"][30]["injected_damp"] = 1280.0, 0.0
    else:
        row["summary"]["final_field_mass"] += 1
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


def test_activation_before_immune_start_is_rejected_even_when_ledgers_reconcile(ledger):
    add_opportunity(ledger, step=59)
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


@pytest.mark.parametrize("mutation", ["activation_sum", "opportunities", "maximum", "threshold", "denominator"])
def test_per_cell_per_step_and_summary_ledgers_must_agree(ledger, mutation):
    add_opportunity(ledger)
    row = ledger["conditions"][0]
    if mutation == "activation_sum":
        next(c for c in row["recipients"] if c["opportunities"])["activation_sum"] += .1
    elif mutation == "opportunities":
        row["steps"][60]["eligible_cells"] += 1
    elif mutation == "maximum":
        row["summary"]["max_eligible_damp"] += 1
    elif mutation == "threshold":
        row["steps"][60]["damp_ge_9kd_opportunities"] = 1
    else:
        row["summary"]["activation_sum_per_reference_recipient"] *= 4
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


@pytest.mark.parametrize("mutation", ["available", "unavailable_response", "fingerprint", "final_fingerprint"])
def test_recipient_interventions_preserve_masks_and_field_invariance(ledger, mutation):
    row = ledger["conditions"][1]
    if mutation == "available":
        row["recipients"][0]["available"] = not row["recipients"][0]["available"]
    elif mutation == "unavailable_response":
        cell = next(c for c in row["recipients"] if not c["available"])
        cell["opportunities"] = 1
    elif mutation == "fingerprint":
        row["steps"][90]["field_fingerprint"] = "1" * 16
    else:
        row["field_fingerprint"] = "1" * 16
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


def test_retained_recipients_have_identical_response_ledgers(ledger):
    add_opportunity(ledger)
    reduced = ledger["conditions"][1]
    for cell in reduced["recipients"]:
        if cell["opportunities"]:
            cell.update(local_damp_sum=100., max_eligible_damp=100., activation_sum=2 / 3)
    reduced["steps"][60].update(max_eligible_damp=100., activation_sum=2 / 3, local_damp_sum=100.)
    reduced["summary"].update(max_eligible_damp=100., activation_sum=2 / 3,
                               activation_sum_per_reference_recipient=(2 / 3) / 3840)
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


def test_attained_maximum_rejects_coherently_understated_activation(ledger):
    add_opportunity(ledger, damp=1.)
    for row in ledger["conditions"][:2]:
        for cell in row["recipients"]:
            if cell["opportunities"]:
                cell["activation_sum"] /= 2
        row["steps"][60]["activation_sum"] /= 2
        row["summary"]["activation_sum"] /= 2
        row["summary"]["activation_sum_per_reference_recipient"] /= 2
    with pytest.raises(ValueError):
        report.validate_observations(ledger, frozen_plan())


def legacy_fixture():
    """Independent, empty 10-by-10 accounting fixture for historical parity."""
    config = json.loads((report.ROOT / report.CONTRACT).read_text())["config"]
    config.update(grid_rows=10, grid_cols=10, tumor_radius_um=80., n_steps=8,
                  immune_start_step=2)
    rows = []
    for index, arm in enumerate(("Control", "RSL3", "SDT")):
        result = {"treatment": arm, "o2_condition": "gradient_120um", "o2_lambda_um": 120.,
                  "immune_mode": "immune_on", "stromal_mode": "off", "total_tumor": 49,
                  "total_dead": 0, "ferroptosis_kills": 0, "immune_kills": 0}
        steps = [{"step": step, "ferroptotic_deaths": 0, "completed_releases": 0,
                  "released_damp": 0., "eligible_cells": 0, "eligible_local_damp_sum": 0.,
                  "immune_kills": 0, "activation_sum": 0., "max_local_damp": None,
                  "damp_ge_kd_opportunities": 0, "damp_ge_9kd_opportunities": 0}
                 for step in range(8)]
        rows.append({"condition_name": f"immune_{arm}", "seed": 42 + index * 10_000_000,
                     "result": result, "final_damp": {"total": 0., "peak": 0.},
                     "measurements": {"ferroptotic_events": [], "eligible_cells": [],
                                      "immune_kill_events": [], "steps": steps,
                                      "terminal": {key: 0 for key in (
                                          "censored_deaths", "terminal_additions", "terminal_damp",
                                          "damp_before_terminal", "damp_after_terminal")}}})
    data = {"schema_version": 2, "simulator": "sim-tme", "dimension": 2,
            "config": config, "conditions": rows}
    historical = {"schema_version": 1, "conditions": [copy.deepcopy(r["result"]) for r in rows]}
    baseline = copy.deepcopy(historical)
    baseline["conditions"] += [{"treatment": "Control", "o2_condition": f"other_{i}",
                                "immune_mode": "off"} for i in range(30)]
    raw, baseline_raw = json.dumps(data).encode(), json.dumps(baseline).encode()
    golden = ("# Source commit: " + "1" * 40 + "\n# Binary SHA-256: " + "2" * 64
              + "\n# Toolchain: rustc 1.96.0; serial execution.\n"
              + report.sha(baseline_raw) + "  tme_summary.json\n")
    return raw, baseline_raw, {
        report.measurement.GOLDEN: golden.encode(),
        report.measurement.HISTORICAL: json.dumps(historical).encode(),
        report.CONTRACT: json.dumps({"schema_version": 2, "config": config}).encode(),
    }


def archive_inputs(ledger):
    raw, baseline, overrides = legacy_fixture()
    canonical = gzip.compress(raw, mtime=0)
    sources = {p: (report.ROOT / p).read_bytes() for p in report.source_paths()}
    sources.update(overrides)
    sources[report.CANONICAL_MANIFEST] = json.dumps({
        "artifacts": {"observations.json.gz": report.sha(canonical)}}).encode()
    artifacts = {"observations.json.gz": gzip.compress(json.dumps(ledger).encode(), mtime=0),
                 "canonical-observations.json.gz": canonical,
                 "baseline-summary.json": baseline, "sources.tar.gz": report.pack_files(sources),
                 "capture-logs.tar.gz": report.pack_files({f"{name}.log": b"synthetic\n"
                                                         for name in ("baseline", "canonical", "controlled")})}
    manifest = {"schema_version": 1, "study": report.STUDY, "execution": "serial", "complete": True,
                "source_commit": "3" * 40, "binary_sha256": "4" * 64,
                "captured_at_utc": "2026-09-22T00:00:00+00:00", "platform": "synthetic",
                "rustc": "rustc 1.96.0", "python": {"implementation": "CPython", "version": "3.14.0"},
                "runs": [{"name": name, "elapsed_seconds": .125}
                         for name in ("baseline", "canonical", "controlled")],
                "baseline_reference": report.measurement.baseline_reference(
                    sources[report.measurement.GOLDEN].decode()),
                "sources": {name: report.sha(blob) for name, blob in sources.items()},
                "artifacts": {name: report.sha(blob) for name, blob in artifacts.items()}}
    return artifacts, sources, manifest


@pytest.fixture
def archive(tmp_path, ledger_template):
    path = tmp_path / "synthetic archive"
    path.mkdir()
    artifacts, _, manifest = archive_inputs(ledger_template)
    for name, blob in artifacts.items():
        (path / name).write_bytes(blob)
    (path / "manifest.json").write_text(json.dumps(manifest))
    return path


def test_offline_archive_reconstruction_includes_all_zero_rows_and_contrasts(archive, monkeypatch):
    monkeypatch.setattr(report.common, "build_binary", lambda *a: pytest.fail("offline reconstruction must not build"))
    result = report.load_archive(archive)
    assert len(result["conditions"]) == 18 and len(result["contrasts"]) == 36
    assert all(c["difference"] == 0 for c in result["contrasts"])
    markdown = report.render(result)
    assert "undefined" in markdown
    assert all(c["condition_id"] in markdown for c in result["conditions"])


@pytest.mark.parametrize("mutation", ["artifact_hash", "source_hash", "missing_artifact",
                                     "extra_file", "symlink", "incomplete", "run_coverage"])
def test_offline_archive_rejects_corrupted_or_incomplete_evidence(archive, mutation):
    path = archive / "manifest.json"
    manifest = json.loads(path.read_text())
    if mutation == "artifact_hash":
        manifest["artifacts"]["observations.json.gz"] = "0" * 64
    elif mutation == "source_hash":
        manifest["sources"][report.PROTOCOL] = "0" * 64
    elif mutation == "missing_artifact":
        (archive / "observations.json.gz").unlink()
    elif mutation == "extra_file":
        (archive / "unplanned.json").write_text("{}")
    elif mutation == "symlink":
        (archive / "capture-logs.tar.gz").unlink()
        (archive / "capture-logs.tar.gz").symlink_to(archive / "sources.tar.gz")
    elif mutation == "incomplete":
        manifest["complete"] = False
    else:
        manifest["runs"].pop()
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        report.load_archive(archive)


def synthetic_capture(monkeypatch, tmp_path, ledger, *, fail_at=None, bad_canonical=False):
    artifacts, sources, _ = archive_inputs(ledger)
    root = tmp_path / "source"
    for name, blob in sources.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(blob)
    legacy_archive = root / "analysis/immune-2d-measurements"
    (legacy_archive / "observations.json.gz").write_bytes(artifacts["canonical-observations.json.gz"])
    binary = tmp_path / "fake-binary"
    binary.write_bytes(b"synthetic executable; never launched")
    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report.measurement, "ARCHIVE", legacy_archive)
    monkeypatch.setattr(report, "source_paths", lambda: sorted(sources))
    monkeypatch.setattr(report.common, "build_binary", lambda *a: binary)
    monkeypatch.setattr(report.platform, "platform", lambda: "synthetic platform")
    monkeypatch.setattr(report.subprocess, "check_output", lambda args, **kw:
                        "" if args[1:3] == ["status", "--porcelain"] else
                        "a" * 40 if args[1:3] == ["rev-parse", "HEAD"] else "rustc 1.96.0")
    calls = []

    def run(args, *, cwd, stdout, **kwargs):
        name = "baseline" if len(args) == 1 else "canonical" if args[1] == "--immune-measurements" else "controlled"
        calls.append(name)
        stdout.write(f"synthetic {name}\n".encode())
        if name == fail_at:
            raise subprocess.CalledProcessError(2, args)
        output = cwd / "output/tme"
        output.mkdir(parents=True)
        if name == "baseline":
            (output / "tme_summary.json").write_bytes(artifacts["baseline-summary.json"])
        elif name == "canonical":
            raw = gzip.decompress(artifacts["canonical-observations.json.gz"])
            (output / "immune_measurements.json").write_bytes(raw + (b" " if bad_canonical else b""))
        else:
            (output / "immune_controlled.json").write_text(json.dumps(ledger))

    monkeypatch.setattr(report.subprocess, "run", run)
    return root / "analysis/new-controlled-study", calls


def test_capture_validates_parity_before_one_controlled_run_and_publishes(monkeypatch, tmp_path, ledger):
    destination, calls = synthetic_capture(monkeypatch, tmp_path, ledger)
    report.capture(destination)
    assert calls == ["baseline", "canonical", "controlled"]
    assert len(report.load_archive(destination)["conditions"]) == 18
    assert not list(destination.parent.glob(".immune-controlled-capture-*"))


@pytest.mark.parametrize("failure", ["baseline", "canonical", "controlled"])
def test_capture_failure_retains_evidence_without_publication_or_retry(monkeypatch, tmp_path, ledger, failure):
    destination, calls = synthetic_capture(monkeypatch, tmp_path, ledger, fail_at=failure)
    with pytest.raises(subprocess.CalledProcessError):
        report.capture(destination)
    assert not destination.exists()
    assert calls[-1] == failure and calls.count(failure) == 1
    [evidence] = list(destination.parent.glob(".immune-controlled-capture-*"))
    assert json.loads((evidence / "failure.json").read_text())["complete"] is False
    assert (evidence / f"{failure}.log").read_text() == f"synthetic {failure}\n"


def test_changed_canonical_bytes_stop_before_controlled_experiment(monkeypatch, tmp_path, ledger):
    destination, calls = synthetic_capture(monkeypatch, tmp_path, ledger, bad_canonical=True)
    with pytest.raises(ValueError, match="canonical observer byte parity"):
        report.capture(destination)
    assert calls == ["baseline", "canonical"] and not destination.exists()


def test_source_changes_during_capture_keep_observations_but_prevent_publication(monkeypatch, tmp_path, ledger):
    destination, calls = synthetic_capture(monkeypatch, tmp_path, ledger)
    original = report.validate_observations

    def change_source(data, plan):
        result = original(data, plan)
        (report.ROOT / report.PROTOCOL).write_text("changed during capture\n")
        return result

    monkeypatch.setattr(report, "validate_observations", change_source)
    with pytest.raises(ValueError, match="frozen sources changed"):
        report.capture(destination)
    assert calls == ["baseline", "canonical", "controlled"] and not destination.exists()
    [evidence] = list(destination.parent.glob(".immune-controlled-capture-*"))
    assert (evidence / "archive/observations.json.gz").is_file()


@pytest.mark.parametrize("guard", ["dirty", "environment", "exists", "dangling"])
def test_capture_guards_fail_before_any_build_or_execution(monkeypatch, tmp_path, guard):
    destination = tmp_path / "archive"
    monkeypatch.setattr(report.common, "build_binary", lambda *a: pytest.fail("must not build"))
    monkeypatch.setattr(report.subprocess, "check_output", lambda *a, **kw: " M protocol.md\n")
    if guard == "environment":
        monkeypatch.setenv("FERRO_UNKNOWN", "")
    elif guard == "exists":
        destination.mkdir()
    elif guard == "dangling":
        destination.symlink_to(tmp_path / "missing")
    pattern = {"dirty": "commit the protocol", "environment": "FERRO_",
               "exists": "already exists", "dangling": "already exists"}[guard]
    with pytest.raises(ValueError, match=pattern):
        report.capture(destination)


def test_probe_damp_cannot_exceed_available_whole_field_mass(ledger):
    add_opportunity(ledger, damp=600.)
    with pytest.raises(ValueError, match="field mass"):
        report.validate_observations(ledger, frozen_plan())


def test_changing_both_availability_trajectory_hashes_cannot_evade_step_chain(ledger):
    for row in ledger["conditions"][:2]:
        row["field_fingerprint"] = "1" * 16
    with pytest.raises(ValueError, match="trajectory fingerprint chain"):
        report.validate_observations(ledger, frozen_plan())


@pytest.mark.skipif(not report.ARCHIVE.exists(), reason="capture follows implementation freeze")
def test_committed_reports_match_complete_offline_archive_reconstruction():
    reconstructed = report.load_archive(report.ARCHIVE)
    assert report.OUT_JSON.read_text() == json.dumps(
        reconstructed, indent=2, sort_keys=True, allow_nan=False) + "\n"
    assert report.OUT_MD.read_text() == report.render(reconstructed)


def test_source_snapshot_collects_its_tests_without_checkout_dependencies(tmp_path):
    for name in report.source_paths():
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((report.ROOT / name).read_bytes())
    collected = subprocess.run(
        [sys.executable, "-I", "-m", "pytest", "--collect-only", "-q",
         "--rootdir", str(tmp_path), "--confcutdir", str(tmp_path),
         "tests/test_immune_2d_controlled_report.py"],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
    )
    assert collected.returncode == 0, collected.stdout + collected.stderr
