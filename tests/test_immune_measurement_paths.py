"""Protect immutable capture archives and reconstruct custom report locations."""

import re
import shlex
import shutil
import sys
from pathlib import Path
from urllib.parse import unquote

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import immune_measurement_report as report


@pytest.mark.parametrize("location", [
    "archive", "new_child", "manifest.json", *sorted(report.ARTIFACTS),
    "symlink", "hardlink", "parent_symlink", "archive_symlink",
])
@pytest.mark.parametrize("capture", [False, True])
def test_report_cannot_modify_archive_before_any_work(tmp_path, monkeypatch, location, capture):
    archive = tmp_path / "archive"
    archive.mkdir()
    for name in report.ARTIFACTS | {"manifest.json"}:
        (archive / name).write_bytes(name.encode())
    before = {path.name: path.read_bytes() for path in archive.iterdir()}
    output = archive / location
    if location == "archive":
        output = archive
    elif location == "symlink":
        output = tmp_path / "report.md"
        output.symlink_to(archive / "manifest.json")
    elif location == "hardlink":
        output = tmp_path / "report.md"
        output.hardlink_to(archive / "baseline-summary.json")
    elif location == "parent_symlink":
        alias = tmp_path / "alias"
        alias.symlink_to(archive, target_is_directory=True)
        output = alias / "nested" / "report.md"
    elif location == "archive_symlink":
        alias = tmp_path / "alias"
        alias.symlink_to(archive, target_is_directory=True)
        output = archive / "report.md"
        archive = alias

    def unexpected_work(*args):
        pytest.fail("output validation must precede capture and archive loading")

    monkeypatch.setattr(report, "capture", unexpected_work)
    monkeypatch.setattr(report, "load_archive", unexpected_work)
    args = ["report", "--archive", str(archive), "--report", str(output)]
    monkeypatch.setattr(sys, "argv", args + (["--capture"] if capture else []))
    with pytest.raises(ValueError, match="immutable archive"):
        report.main()
    assert {path.name: path.read_bytes() for path in archive.iterdir()} == before


def test_capture_rejects_report_in_future_archive_through_parent_symlink(tmp_path, monkeypatch):
    parent = tmp_path / "real"
    parent.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(parent, target_is_directory=True)
    archive = parent / "new-archive"
    monkeypatch.setattr(sys, "argv", ["report", "--capture", "--archive", str(archive),
                                     "--report", str(alias / "new-archive" / "report.md")])
    monkeypatch.setattr(report, "capture", lambda *args: pytest.fail("capture started"))
    with pytest.raises(ValueError, match="outside the immutable archive"):
        report.main()
    assert not archive.exists()


@pytest.mark.parametrize("path_style", ["ordinary", "parent_traversal", "symlink_parent_traversal"])
def test_custom_report_links_and_command_reconstruct_the_same_archive(tmp_path, monkeypatch,
                                                                     path_style):
    archive = tmp_path / "archive #1 (frozen)"
    shutil.copytree(report.ARCHIVE, archive)
    before = {path.name: report.sha(path.read_bytes()) for path in archive.iterdir()}
    output = tmp_path / "reports with spaces" / "replicate.md"
    monkeypatch.chdir(tmp_path)
    output_arg = output.relative_to(tmp_path)
    archive_arg = Path(archive.name)
    if path_style == "parent_traversal":
        output_arg = Path(archive.name) / ".." / output_arg
    elif path_style == "symlink_parent_traversal":
        # Filesystem traversal follows the symlink before '..'; lexical
        # normalization would instead point the command into 'links'.
        nested = tmp_path / "nested"
        nested.mkdir()
        alias = tmp_path / "links" / "alias"
        alias.parent.mkdir()
        alias.symlink_to(nested, target_is_directory=True)
        output_arg = Path("links/alias") / ".." / output_arg
        archive_arg = Path("links/alias") / ".." / archive_arg
    monkeypatch.setattr(sys, "argv", ["report", "--archive", str(archive_arg),
                                     "--report", str(output_arg)])
    report.main()
    text = output.read_text()
    links = dict(re.findall(r"\[([^\]]+)\]\(([^)]+)\)", text))
    for label, target in {
        "protocol": report.ROOT / "docs/IMMUNE_MEASUREMENT_PROTOCOL.md",
        "Archive and provenance": archive / "manifest.json",
        "Unchanged production summary": archive / "baseline-summary.json",
    }.items():
        assert (output.parent / unquote(links[label])).resolve() == target.resolve()
        if label != "protocol":
            assert "%23" in links[label]
    command = shlex.split(re.search(r"```bash\n([^\n]+)\n```", text).group(1))
    assert command[:2] == ["python3", str(report.ROOT / "scripts/immune_measurement_report.py")]
    assert command[2:] == ["--archive", str(archive.resolve()), "--report", str(output.resolve())]
    monkeypatch.chdir(report.ROOT)
    monkeypatch.setattr(sys, "argv", command[1:])
    report.main()
    assert output.read_text() == text
    assert {path.name: report.sha(path.read_bytes()) for path in archive.iterdir()} == before
    report.load_archive(archive)


def test_default_report_stays_byte_identical():
    manifest, summaries = report.load_archive(report.ARCHIVE)
    assert report.render(manifest, summaries, report.ARCHIVE, report.REPORT) == report.REPORT.read_text()
