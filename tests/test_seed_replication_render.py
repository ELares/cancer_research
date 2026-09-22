"""Historical resampling prose must preserve its data and render offline."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "seed_replication_render", REPO / "scripts/seed_replication.py")
reporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reporter)
RAW = REPO / "analysis/seed-replication.json"
REPORT = REPO / "analysis/seed-replication-report.md"


def test_historical_numeric_artifacts_and_all_table_rows_are_preserved():
    assert hashlib.sha256(RAW.read_bytes()).hexdigest() == (
        "dd3c9fa10e604ce9dd67e43ebe04c2c2e24b8acb252ed9c86dcbacd083f3dfea")
    text = reporter.render(json.loads(RAW.read_text()))
    numeric_rows = [line for line in text.splitlines()
                    if re.match(r"^\| (Control|RSL3|SDT|[012]) \|", line)]
    assert len(numeric_rows) == 157
    assert hashlib.sha256("\n".join(numeric_rows).encode()).hexdigest() == (
        "16028489f8b83c20e6b105425858de16db05bf46a738e6a81faf9e5dfd49a430")
    assert text == REPORT.read_text()


def test_renderer_preserves_inputs_and_names_the_uncertainty_limit():
    data = json.loads(RAW.read_text())
    original = copy.deepcopy(data)
    text = " ".join(reporter.render(data).split())
    assert data == original
    assert "independent-run uncertainty is unvalidated" in text
    assert "Consecutive roots 42–61 reuse additive per-cell RNG seed arguments" in text
    assert "does not invalidate the observed counts" in text
    assert "magnitude or sign of outcome correlation or bias" in text
    assert "does not establish representativeness across independent runs" in text
    assert "historical 95% bootstrap interval" in text
    assert "IMMUNE_2D_REPLICATION_PROTOCOL.md" in text
    assert "95% CI" not in text


def test_render_only_needs_no_binary_simulation_or_resampling(monkeypatch, tmp_path):
    raw = tmp_path / "historical.json"
    output = tmp_path / "historical.md"
    original = RAW.read_bytes()
    raw.write_bytes(original)
    monkeypatch.setattr(reporter, "RAW", raw)
    monkeypatch.setattr(reporter, "OUT", output)
    monkeypatch.setattr(reporter, "BINARY", tmp_path / "missing-simulator")

    def forbidden(*args, **kwargs):
        raise AssertionError("offline rendering must not simulate or resample")

    monkeypatch.setattr(reporter, "run_seed", forbidden)
    monkeypatch.setattr(reporter, "boot_ci", forbidden)
    monkeypatch.setattr(reporter, "ThreadPoolExecutor", forbidden)
    monkeypatch.setattr(sys, "argv", ["seed_replication.py", "--render-only"])
    reporter.main()
    assert raw.read_bytes() == original
    assert output.read_text() == REPORT.read_text()


def test_renderer_displays_saved_bootstrap_endpoints_without_recalculating():
    data = json.loads(RAW.read_text())
    data["rows"][0]["lo"], data["rows"][0]["hi"] = 123.0, 456.0
    assert "| 123-456 |" in reporter.render(data)
