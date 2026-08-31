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
MODEL_CARD = REPO_ROOT / "MODEL_CARD.md"


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


def test_model_card_version_matches_cargo():
    """MODEL_CARD.md's stated ferroptosis-core version must match Cargo.toml.
    Guards the drift that left it at 0.49.0 while the crate moved on (#333 review);
    MODEL_CARD has its own version line that no other test covered."""
    text = MODEL_CARD.read_text(encoding="utf-8")
    stated = re.findall(r"ferroptosis-core (\d+\.\d+\.\d+)", text)
    assert stated, "No 'ferroptosis-core X.Y.Z' string found in MODEL_CARD.md"
    actual = _core_version()
    for v in stated:
        assert v == actual, (
            f"MODEL_CARD.md states ferroptosis-core {v}, but "
            f"simulations/ferroptosis-core/Cargo.toml is {actual}. "
            f"Update MODEL_CARD.md when the crate version changes."
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


def test_manuscript_scopes_the_acsl4_prevalence_claim_to_the_deep_cut():
    """#616: the z<-1 low-ACSL4 fraction is not gene-specific; the z<-2 one is.

    The manuscript originally called the shallow fraction "the population prior
    for how many tumors fall in the refractory-leaning tail". Six control genes
    return the same figure, so it describes the cut-point rather than ACSL4.
    Guarded because a prevalence number that reads as biology is exactly the
    kind of claim that survives a correction elsewhere in the repo.
    """
    md = (Path(__file__).resolve().parent.parent
          / "article" / "drafts" / "v1.md").read_text()
    assert "#616" in md, "the manuscript does not cite the correction"
    assert "is not gene-specific" in md
    # The shallow figure must never again be presented as the prior on its own.
    assert ("nineteen percent of tumors per cancer type) is the population prior"
            not in md), "the retracted framing has returned"
    # And the surviving half must still be stated, or the correction over-reaches.
    assert "deep-cut figure is genuine evidence about ACSL4" in md


def test_every_binary_is_listed_in_the_simulations_index():
    """The count being right did not make the TABLE right.

    `simulations/README.md` carries the per-binary table, and README.md links to
    it as "[12 Rust binaries](simulations/README.md)". sim-scale was absent from
    that table while every count in the repo already said 12, so the front door
    linked to a page documenting eleven. A count guard cannot see this: it
    compares two numbers and never opens the page the number points at.
    """
    root = Path(__file__).resolve().parent.parent
    crates = sorted(p.name for p in (root / "simulations").glob("sim-*")
                    if p.is_dir() and (p / "Cargo.toml").exists())
    assert crates, "no sim-* crates found; the glob is wrong"
    index = (root / "simulations" / "README.md").read_text()
    missing = [c for c in crates if f"`{c}`" not in index]
    assert not missing, (
        "simulations/README.md does not list " + ", ".join(missing)
        + f" -- it documents {len(crates) - len(missing)} of {len(crates)} "
        "binaries while README.md links to it as the index of all of them")


def test_documented_binary_count_matches_the_crates_that_exist():
    """The "N binaries" figure in README.md and CLAUDE.md must be derived.

    Adding sim-scale made it 12 and three inventory statements still said 11 --
    including the CLAUDE.md line the same PR was editing to bump a different
    counter, so a human touched that exact line and did not see the stale number
    beside it. A count written in prose next to a directory that can grow is a
    number waiting to rot, so this recomputes it.
    """
    import re
    root = Path(__file__).resolve().parent.parent
    actual = len([p for p in (root / "simulations").glob("sim-*")
                  if p.is_dir() and (p / "Cargo.toml").exists()])
    assert actual > 0, "no sim-* crates found; the glob is wrong"
    for name in ("README.md", "CLAUDE.md"):
        txt = (root / name).read_text()
        for m in re.finditer(r"(\d+) binaries", txt):
            assert int(m.group(1)) == actual, (
                f"{name} says {m.group(1)} binaries; {actual} sim-* crates exist "
                f"({', '.join(sorted(p.name for p in (root/'simulations').glob('sim-*') if p.is_dir()))})")


def _core_version() -> str:
    m = re.search(r'^version\s*=\s*"([^"]+)"', CORE_CARGO.read_text(), re.M)
    assert m, "no version in ferroptosis-core/Cargo.toml"
    return m.group(1)


def test_claude_md_states_the_live_crate_version():
    """CLAUDE.md asserted 0.67.0 and 0.69.0 in the same paragraph.

    Appendix A and MODEL_CARD.md were already pinned here and were both
    correct; the one file that is loaded into every session's context was not
    covered, drifted two minor versions, and contradicted itself. A version a
    reader is told is current is a claim like any other.
    """
    claude = (REPO_ROOT / "CLAUDE.md").read_text()
    version = _core_version()
    stated = re.findall(r"current crate version ([0-9]+\.[0-9]+\.[0-9]+)", claude)
    assert stated, "CLAUDE.md no longer states a current crate version"
    for v in stated:
        assert v == version, (
            f"CLAUDE.md says the crate is at {v}; Cargo.toml says {version}")


def test_the_crate_readme_documents_every_module():
    """`README.md` calls this table THE module list, so an omission is a
    reader looking up a module and concluding it does not exist.

    It had 31 rows against 37 files -- five modality modules missing, four of
    them for several PRs and one added by the PR that also bumped the module
    count in six other places. Nothing had ever compared the two.
    """
    src = REPO_ROOT / "simulations" / "ferroptosis-core" / "src"
    on_disk = {p.stem for p in src.glob("*.rs") if p.stem != "lib"}
    readme = (REPO_ROOT / "simulations" / "ferroptosis-core" / "README.md").read_text()
    section = readme.split("## Modules", 1)
    assert len(section) == 2, "the crate README has no Modules section"
    body = section[1].split("\n## ", 1)[0]
    documented = set(re.findall(r"^\| `([a-z_0-9]+)`", body, re.M))
    assert on_disk - documented == set(), (
        "modules on disk with no row in the crate README: "
        f"{sorted(on_disk - documented)}")
    assert documented - on_disk == set(), (
        "the crate README documents modules that do not exist: "
        f"{sorted(documented - on_disk)}")
