"""Exercise publication and failure evidence using a synthetic executable, never Rust."""

import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import immune_2d_replication_report as report
from test_immune_2d_measurement_report import source_fixture


def synthetic_capture(monkeypatch, tmp_path, *, fail_at=None, corrupt_canonical=False):
    root = tmp_path / "source"
    root.mkdir()
    canonical_raw, baseline, sources = source_fixture()
    canonical = gzip.compress(canonical_raw, mtime=0)
    sources[report.PLAN] = json.dumps(report.FIXED_PLAN).encode()
    sources[report.CANONICAL_MANIFEST] = json.dumps({
        "artifacts": {"observations.json.gz": report.sha(canonical)}}).encode()
    sources[report.PROTOCOL] = b"synthetic protocol"
    for name in report.source_paths():
        sources.setdefault(name, b"synthetic frozen source")
    for name, content in sources.items():
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    old_archive = root / "analysis/immune-2d-measurements"
    (old_archive / "observations.json.gz").write_bytes(canonical)
    binary = tmp_path / "fake-binary"
    binary.write_bytes(b"synthetic executable, not executed")
    monkeypatch.setattr(report, "ROOT", root)
    monkeypatch.setattr(report.measurement, "ARCHIVE", old_archive)
    monkeypatch.setattr(report, "source_paths", lambda: sorted(sources))
    monkeypatch.setattr(report.common, "build_binary", lambda *args: binary)
    monkeypatch.setattr(report.platform, "platform", lambda: "synthetic platform")
    monkeypatch.setattr(report.subprocess, "check_output", lambda args, **kwargs:
                        "" if args[1:3] == ["status", "--porcelain"] else
                        "a" * 40 if args[1:3] == ["rev-parse", "HEAD"] else "rustc 1.96.0")
    calls = []

    def run(args, *, cwd, stdout, **kwargs):
        args = args[1:]
        name = "baseline" if not args else "canonical" if args == ["--immune-measurements"] else int(args[1])
        calls.append(name)
        stdout.write(f"synthetic run {name}\n".encode())
        if name == fail_at:
            raise subprocess.CalledProcessError(2, args)
        out = cwd / "output/tme"
        out.mkdir(parents=True)
        if name == "baseline":
            (out / "tme_summary.json").write_bytes(baseline)
        elif name == "canonical":
            (out / "immune_measurements.json").write_bytes(canonical_raw + (b" " if corrupt_canonical else b""))
        else:
            data = json.loads(canonical_raw)
            data["replicate_block"] = name
            data["config"]["seed"] = report.seeds.block_seed(name)
            for index, row in enumerate(data["conditions"]):
                row["seed"] = data["config"]["seed"] + index * 10000000
            (out / "immune_replicate.json").write_text(json.dumps(data))

    monkeypatch.setattr(report.subprocess, "run", run)
    return root / "analysis/new-study", calls


def test_success_runs_every_block_once_after_parity_and_publishes_verified_archive(monkeypatch, tmp_path):
    destination, calls = synthetic_capture(monkeypatch, tmp_path)
    report.capture(destination)
    assert calls == ["baseline", "canonical", *range(1, 21)]
    reconstructed = report.load_archive(destination)
    assert len(reconstructed["blocks"]) == 20
    assert reconstructed["provenance"]["source_commit"] == "a" * 40
    assert not list(destination.parent.glob(".immune-replication-capture-*"))


@pytest.mark.parametrize("failure", ["baseline", "canonical", 1, 7, 20])
def test_failure_keeps_logs_and_prior_observations_without_publishing_or_retrying(monkeypatch, tmp_path, failure):
    destination, calls = synthetic_capture(monkeypatch, tmp_path, fail_at=failure)
    with pytest.raises(subprocess.CalledProcessError):
        report.capture(destination)
    assert not destination.exists()
    assert calls[-1] == failure and calls.count(failure) == 1
    [evidence] = list(destination.parent.glob(".immune-replication-capture-*"))
    assert json.loads((evidence / "failure.json").read_text())["complete"] is False
    name = f"block-{failure:02d}" if type(failure) is int else failure
    assert (evidence / f"{name}.log").read_text() == f"synthetic run {failure}\n"
    if type(failure) is int:
        assert sorted(p.name for p in (evidence / "archive").glob("block-*.json.gz")) == [
            report.block_name(b) for b in range(1, failure)]


def test_even_semantically_identical_canonical_json_must_match_bytes_before_new_blocks(monkeypatch, tmp_path):
    destination, calls = synthetic_capture(monkeypatch, tmp_path, corrupt_canonical=True)
    with pytest.raises(ValueError, match="canonical observer byte parity"):
        report.capture(destination)
    assert calls == ["baseline", "canonical"] and not destination.exists()


@pytest.mark.parametrize("guard", ["dirty", "environment", "exists", "dangling"])
def test_capture_preconditions_fail_before_build_or_simulation(monkeypatch, tmp_path, guard):
    destination, calls = synthetic_capture(monkeypatch, tmp_path)
    if guard == "dirty":
        monkeypatch.setattr(report.subprocess, "check_output", lambda *args, **kwargs: " M frozen.py\n")
    elif guard == "environment":
        monkeypatch.setenv("FERRO_UNKNOWN", "")
    elif guard == "exists":
        destination.mkdir()
    else:
        destination.symlink_to(tmp_path / "missing")
    monkeypatch.setattr(report.common, "build_binary", lambda *args: pytest.fail("must not build"))
    with pytest.raises(ValueError):
        report.capture(destination)
    assert not calls


def test_source_change_during_capture_keeps_evidence_without_publication(monkeypatch, tmp_path):
    destination, calls = synthetic_capture(monkeypatch, tmp_path)
    original = report.validate_block

    def change_source(data, block, contract):
        result = original(data, block, contract)
        if block == 20:
            (report.ROOT / report.PROTOCOL).write_text("changed during capture")
        return result

    monkeypatch.setattr(report, "validate_block", change_source)
    with pytest.raises(ValueError, match="sources changed"):
        report.capture(destination)
    assert calls[-1] == 20 and not destination.exists()
    [evidence] = list(destination.parent.glob(".immune-replication-capture-*"))
    assert len(list((evidence / "archive").glob("block-*.json.gz"))) == 20
