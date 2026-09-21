"""Capture must use the executable from the actual pinned Cargo build."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import immune_measurement_report as report


def artifact(executable, name="sim-tme-3d", kind="bin", fresh=False):
    return {"reason": "compiler-artifact", "target": {"name": name, "kind": [kind]},
            "executable": str(executable) if executable else None, "fresh": fresh}


@pytest.mark.parametrize("target", [None, "aarch64-apple-darwin"])
@pytest.mark.parametrize("fresh", [False, True])
def test_build_uses_emitted_executable_including_configured_target(tmp_path, monkeypatch, target, fresh):
    cache = tmp_path / "shared-target"
    stale = cache / "release" / "sim-tme-3d"
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b"stale host build")
    actual = cache / target / "release" / "sim-tme-3d" if target else stale
    actual.parent.mkdir(parents=True, exist_ok=True)
    actual.write_bytes(b"requested build")
    if target:
        monkeypatch.setenv("CARGO_BUILD_TARGET", target)
    cargo = ["rustup", "run", report.TOOLCHAIN, "cargo"]
    messages = [artifact(None, name="dep", kind="lib"),
                artifact(tmp_path / "other-program", name="other-program"),
                artifact(actual, fresh=fresh), {"reason": "build-finished", "success": True}]

    def build(command, **kwargs):
        assert command == [*cargo, "build", "--locked", "--release", "-p", "sim-tme-3d",
                           "--message-format=json-render-diagnostics"]
        assert kwargs == {"cwd": tmp_path, "check": True, "stdout": subprocess.PIPE, "text": True}
        return subprocess.CompletedProcess(command, 0, "\n".join(map(json.dumps, messages)) + "\n")

    monkeypatch.setattr(report.subprocess, "run", build)
    assert report.build_binary(cargo, tmp_path) == actual
    if target:
        assert stale.read_bytes() == b"stale host build"


@pytest.mark.parametrize("case", ["missing", "ambiguous", "absent_file", "failed"])
def test_build_cannot_fall_back_to_an_unidentified_binary(tmp_path, monkeypatch, case):
    first, second = tmp_path / "first", tmp_path / "second"
    first.write_bytes(b"first")
    second.write_bytes(b"second")
    messages = []
    if case == "ambiguous":
        messages = [artifact(first), artifact(second)]
    elif case == "absent_file":
        messages = [artifact(tmp_path / "does-not-exist")]

    def build(command, **kwargs):
        if case == "failed":
            raise subprocess.CalledProcessError(1, command)
        return subprocess.CompletedProcess(command, 0, "\n".join(map(json.dumps, messages)))

    monkeypatch.setattr(report.subprocess, "run", build)
    expected = subprocess.CalledProcessError if case == "failed" else ValueError
    with pytest.raises(expected):
        report.build_binary(["cargo"], tmp_path)
