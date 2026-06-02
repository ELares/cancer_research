"""
Guards the manuscript (article/drafts/v1.md) against the inventory drift that
PR #284 had to fix after the fact: the ferroptosis-core version string and the
build-list of simulation binaries silently fell out of sync with the source.

This catches the next drift at PR time instead of in a later audit.

Checks:
- The `ferroptosis-core version X.Y.Z` string in Appendix A matches the actual
  version in simulations/ferroptosis-core/Cargo.toml.
- Every simulation binary crate (simulations/sim-*) is named somewhere in the
  manuscript (Appendix B's build list is the intended home).

Run: pytest tests/test_manuscript_inventory.py -v
"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT = REPO_ROOT / "article" / "drafts" / "v1.md"
CORE_CARGO = REPO_ROOT / "simulations" / "ferroptosis-core" / "Cargo.toml"
SIM_DIR = REPO_ROOT / "simulations"


def _manuscript_text() -> str:
    return MANUSCRIPT.read_text(encoding="utf-8")


def _core_version() -> str:
    for line in CORE_CARGO.read_text(encoding="utf-8").splitlines():
        m = re.match(r'\s*version\s*=\s*"([^"]+)"', line)
        if m:
            return m.group(1)
    pytest.fail(f"Could not find a version in {CORE_CARGO}")


def _sim_binary_crates() -> list[str]:
    return sorted(
        p.name
        for p in SIM_DIR.glob("sim-*")
        if p.is_dir() and (p / "Cargo.toml").exists()
    )


def test_manuscript_version_matches_cargo():
    """Appendix A's stated ferroptosis-core version must match Cargo.toml."""
    text = _manuscript_text()
    stated = re.findall(r"ferroptosis-core version (\d+\.\d+\.\d+)", text)
    assert stated, "No 'ferroptosis-core version X.Y.Z' string found in the manuscript"
    actual = _core_version()
    for v in stated:
        assert v == actual, (
            f"Manuscript states ferroptosis-core version {v}, but "
            f"simulations/ferroptosis-core/Cargo.toml is {actual}. "
            f"Update Appendix A (and CLAUDE.md) when the crate version changes."
        )


def test_all_sim_binaries_appear_in_manuscript():
    """Every simulations/sim-* binary crate must be named in the manuscript."""
    text = _manuscript_text()
    missing = []
    for name in _sim_binary_crates():
        # Match the name only when NOT followed by another word char or hyphen,
        # so "sim-combo" does not spuriously match "sim-combo-mech".
        if not re.search(re.escape(name) + r"(?![\w-])", text):
            missing.append(name)
    assert not missing, (
        "Simulation binaries missing from the manuscript build list "
        f"(Appendix B): {missing}. Add them so the reproduction guide stays complete."
    )
