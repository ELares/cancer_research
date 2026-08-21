"""Numeric claims the docs make about the repo, checked against the repo.

A recurring defect class in this project, found four separate times in one
session: a claim in one file about the CONTENTS of another, with nothing tying
them. Figure captions describing absent figures. Three files quoting three
different manuscript lengths. A manifest missing the files that added it. A
budget table whose rows were kept correct by moving words between them.

Each was fixed where it was found. This is the class rather than an instance:
every countable claim the reader-facing docs make about the repo, verified
against the thing it counts.

DELIBERATELY NARROW. A guard broad enough to catch every number in prose would
flag measurements, historical figures, and retracted claims quoted as history
-- and this repo has already learned from `audit_heading_only_assertions.py`
that a check at roughly half precision trains exemption-adding rather than
fixing. So the set here is explicit: each entry names a claim, where it is
made, and how to count the thing it describes. Adding a number to a doc does
not enrol it; someone has to decide it is worth pinning.

SCOPE LIMIT WORTH STATING: this checks that a count is right, not that the
sentence around it is. `FIGURES.yaml indexes 30 entries` can be true while the
description beside it is nonsense.
"""
import re
import subprocess
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
README = REPO / "README.md"
CLAUDE = REPO / "CLAUDE.md"


def _figures():
    d = yaml.safe_load((REPO / "FIGURES.yaml").read_text())
    return d["figures"] if isinstance(d, dict) and "figures" in d else d


def _count_figure_entries():
    return len(_figures())


def _count_manuscript_figures():
    return sum(1 for f in _figures() if f.get("status") == "manuscript")


def _count_core_modules():
    src = REPO / "simulations/ferroptosis-core/src"
    if not src.is_dir():
        pytest.skip("ferroptosis-core not present")
    return len([p for p in src.glob("*.rs") if p.stem not in ("lib", "main")])


def _count_sim_binaries():
    sims = REPO / "simulations"
    if not sims.is_dir():
        pytest.skip("simulations not present")
    return len({p.parent.parent.name for p in sims.glob("sim-*/src/main.rs")})


def _count_c04_descriptors():
    tsv = REPO / "corpus/atlas/mesh/c04-descriptors.tsv"
    if not tsv.exists():
        pytest.skip("C04 descriptor list is gitignored census data")
    return sum(1 for ln in tsv.read_text().splitlines()
               if ln.strip() and not ln.startswith("#"))


def _count_sites():
    tsv = REPO / "analysis/site-descriptor-map.tsv"
    return len({ln.split("\t")[0] for ln in tsv.read_text().splitlines()
                if ln.strip() and not ln.startswith("#")})


def _count_chains():
    import sys

    sys.path.insert(0, str(REPO / "scripts"))
    import config

    return len(config.DIAGNOSTIC_THERAPY_KEYWORDS)


def _count_chapters():
    md = (REPO / "article/drafts/v1.md").read_text()
    return len(re.findall(r"^## Chapter ", md, re.M))


# (label, file, regex capturing the claimed number, how to count it)
CLAIMS = [
    ("FIGURES.yaml entries", README, r"FIGURES\.yaml\)? indexes (\d+) entries",
     _count_figure_entries),
    ("manuscript figures", README, r"(\d+) figures", _count_manuscript_figures),
    ("manuscript chapters", README, r"(\d+) chapters", _count_chapters),
    ("ferroptosis-core modules", CLAUDE, r"MIT licensed, (\d+) modules",
     _count_core_modules),
    ("simulation binaries", README, r"(\d+) binaries", _count_sim_binaries),
    ("C04 descriptors", CLAUDE, r"C04 \((\d+) descriptors", _count_c04_descriptors),
    ("diagnostic chains", CLAUDE, r"covers (\d+) chains", _count_chains),
]


@pytest.mark.parametrize("label,path,pattern,counter",
                         CLAIMS, ids=[c[0] for c in CLAIMS])
def test_a_documented_count_matches_the_thing_it_counts(label, path, pattern,
                                                        counter):
    text = path.read_text()
    m = re.search(pattern, text)
    assert m, (
        f"{path.name} no longer states a count for {label} (pattern "
        f"{pattern!r}). Either the claim was removed -- fine, drop this entry "
        f"-- or it was reworded, in which case this guard has gone silent "
        f"while still passing.")
    claimed = int(m.group(1))
    actual = counter()
    assert claimed == actual, (
        f"{path.name} says {claimed} {label}; the repo has {actual}")


def test_the_site_count_agrees_across_every_doc_that_states_it():
    """Four sites in CLAUDE.md quote the 18-site list.

    A number repeated in four places is one that drifts in three of them, which
    is exactly how the manuscript's length ended up with three different wrong
    values in three files.
    """
    actual = _count_sites()
    bad = []
    for path in (README, CLAUDE):
        for m in re.finditer(r"(\d+) (?:major )?sites", path.read_text()):
            if int(m.group(1)) != actual:
                bad.append(f"{path.name}: {m.group(1)} against {actual}")
    assert not bad, "site-count claims disagreeing with the map:\n  " + "\n  ".join(bad)


def test_the_guarded_set_is_declared_rather_than_inferred():
    """A guard that enrolled every number in prose would flag measurements,
    historical figures and retracted claims quoted as history.

    This repo already learned from audit_heading_only_assertions.py that a
    check at roughly half precision trains exemption-adding rather than
    fixing, so enrolment is a decision someone makes.
    """
    assert len(CLAIMS) >= 5
    src = Path(__file__).read_text()
    assert "DELIBERATELY NARROW" in src
    # Every entry must carry a counter that actually reads the repo, not a
    # literal -- otherwise the guard compares a doc to a hardcoded number and
    # both go stale together.
    for label, _, _, counter in CLAIMS:
        assert callable(counter), label
        # It must READ something, not return a literal -- otherwise the
        # guard compares a doc against a hardcoded number and the two go
        # stale together, which is the failure this file exists for one level
        # up. Checked structurally: the function's source must reference the
        # repo root or call another counter that does.
        import inspect

        src_fn = inspect.getsource(counter)
        assert "REPO" in src_fn or "_figures()" in src_fn, (
            f"{label}'s counter reads no path, so it is comparing the doc to "
            "a constant rather than to the repo")
