"""Exactly one document may designate the keystone experiment (#619).

Two planning documents once named different primary experiments, each on its own
authority and neither citing evidence. That is the failure this guards: not a
wrong choice, but two documents each quietly believing they made it.

`PREREGISTRATION.md` is the registering document and holds the designation.
`analysis/p1-wetlab-protocol.md` defers to it.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PREREG = REPO_ROOT / "PREREGISTRATION.md"
PROTOCOL = REPO_ROOT / "analysis" / "p1-wetlab-protocol.md"


def test_exactly_one_experiment_is_designated_the_keystone():
    txt = PREREG.read_text()
    headings = re.findall(r"^### (E\d)\..*$", txt, re.M)
    assert headings, "the experiment briefs are gone"
    designated = [h for h in re.findall(r"^### (E\d)\.([^\n]*)$", txt, re.M)
                  if "keystone" in h[1].lower()]
    assert len(designated) == 1, (
        f"expected exactly one keystone experiment, found {designated}")
    assert designated[0][0] == "E4", (
        f"the keystone is {designated[0][0]}; if that is intentional, update this "
        "guard and the basis stated alongside it")


def test_the_basis_is_stated_at_the_designation():
    """#619's acceptance criterion: keystone-because-decisive and
    keystone-because-novel point at different experiments, so which one is
    meant has to be written down."""
    txt = PREREG.read_text()
    assert "Why this is the keystone:" in txt
    assert "most DECISIVE, not the most" in txt
    # And the experiment that lost the label must keep an accurate one rather
    # than being silently demoted.
    assert "the most novel test" in txt


def test_the_protocol_defers_rather_than_designating():
    txt = PROTOCOL.read_text()
    assert "designation is `PREREGISTRATION.md`'s, not this document's" in txt
    # It must not reassert its own competing claim.
    assert "highest-leverage" not in txt.split("used to call")[0], (
        "the protocol is designating a keystone on its own authority again")


def test_the_literature_position_table_is_retained():
    """The choice must stay auditable: the numbers that informed it, and the
    statement that they are not a quality ranking, both remain."""
    txt = PREREG.read_text()
    assert "P1, P3" in txt and "479" in txt
    assert "not a quality ranking" in txt
    assert "A sparse leg is where novelty lives" in txt


def test_the_sequencing_recruits_for_the_keystone_first():
    txt = PREREG.read_text()
    seq = txt[txt.index("## Sequencing"):]
    assert "E4" in seq and "FIRST" in seq, (
        "sequencing does not put the keystone first")
    e4, e1 = seq.index("E4"), seq.index("E1")
    assert e4 < e1, "E1 is still recruited before the keystone"
