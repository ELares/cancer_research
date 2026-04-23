"""
Validates FIGURES.yaml stays in sync with article/figures/ on disk.

Catches:
- Figure PDFs on disk without a YAML entry (new figure, forgot to update)
- YAML entries without a PDF on disk (deleted figure, forgot to update)
- Manuscript figure numbers 1-19 present and unique
- Generator functions exist in the claimed scripts

Run: pytest tests/test_figure_traceability.py -v
"""

from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURES_YAML = REPO_ROOT / "FIGURES.yaml"
FIGURES_DIR = REPO_ROOT / "article" / "figures"


@pytest.fixture
def figures_data():
    if not FIGURES_YAML.exists():
        pytest.skip("FIGURES.yaml not found")
    with open(FIGURES_YAML) as f:
        data = yaml.safe_load(f)
    return data["figures"]


class TestFigureTraceability:
    def test_yaml_is_valid(self):
        """FIGURES.yaml parses without error."""
        assert FIGURES_YAML.exists(), "FIGURES.yaml not found"
        with open(FIGURES_YAML) as f:
            data = yaml.safe_load(f)
        assert "figures" in data
        assert len(data["figures"]) > 0

    def test_every_yaml_entry_has_pdf(self, figures_data):
        """Every YAML entry has a corresponding PDF on disk."""
        for entry in figures_data:
            fn = entry["filename"]
            pdf = FIGURES_DIR / f"{fn}.pdf"
            assert pdf.exists(), f"YAML lists {fn} but {pdf} not found on disk"

    def test_every_pdf_has_yaml_entry(self, figures_data):
        """Every PDF in article/figures/ has a YAML entry."""
        yaml_filenames = {e["filename"] for e in figures_data}
        for pdf in sorted(FIGURES_DIR.glob("*.pdf")):
            basename = pdf.stem
            assert basename in yaml_filenames, (
                f"{pdf.name} exists on disk but has no FIGURES.yaml entry"
            )

    def test_manuscript_figures_complete(self, figures_data):
        """Manuscript figures 1-19 are all present."""
        manuscript_nums = {
            e["manuscript_figure"]
            for e in figures_data
            if e["manuscript_figure"] is not None
        }
        expected = set(range(1, 20))
        missing = expected - manuscript_nums
        assert not missing, f"Missing manuscript figures: {missing}"

    def test_manuscript_figures_unique(self, figures_data):
        """No duplicate manuscript figure numbers."""
        nums = [
            e["manuscript_figure"]
            for e in figures_data
            if e["manuscript_figure"] is not None
        ]
        assert len(nums) == len(set(nums)), f"Duplicate manuscript figures: {nums}"

    def test_required_fields(self, figures_data):
        """Every entry has all required fields."""
        required = {"manuscript_figure", "filename", "generator", "inputs", "type", "status"}
        for entry in figures_data:
            missing = required - set(entry.keys())
            assert not missing, (
                f"{entry.get('filename', '?')} missing fields: {missing}"
            )

    def test_generator_functions_exist(self, figures_data):
        """Generator functions claimed in YAML exist in the claimed scripts."""
        for entry in figures_data:
            func = entry["generator"].get("function")
            script = entry["generator"].get("script")
            if func is None or script is None:
                continue  # Rust binaries and orphans — skip
            if script.endswith(".rs"):
                continue  # Rust source — skip function check
            script_path = REPO_ROOT / script
            assert script_path.exists(), (
                f"{entry['filename']}: script {script} not found"
            )
            content = script_path.read_text()
            assert f"def {func}" in content, (
                f"{entry['filename']}: function {func} not found in {script}"
            )

    def test_valid_types(self, figures_data):
        """Figure types are from the allowed set."""
        allowed = {"corpus-derived", "simulation", "conceptual", "unknown"}
        for entry in figures_data:
            assert entry["type"] in allowed, (
                f"{entry['filename']} has invalid type: {entry['type']}"
            )

    def test_valid_statuses(self, figures_data):
        """Figure statuses are from the allowed set."""
        allowed = {"manuscript", "supplementary", "orphan"}
        for entry in figures_data:
            assert entry["status"] in allowed, (
                f"{entry['filename']} has invalid status: {entry['status']}"
            )
