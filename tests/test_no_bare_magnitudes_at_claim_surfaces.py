"""No headline magnitude may appear on a public claim surface without its interval.

WHY THIS FILE EXISTS
--------------------
`analysis/identifiability-report.json` records, per headline, whether it is
point-estimable. For every one it is not: eleven free rate constants, six
non-identifiable from the kill rate, and ZERO of the headline outputs conditioned
on data in the regime that produces them.

Issue #506 fixed that in the manuscript's abstract. It never reached the two
surfaces that travel furthest:

  README.md              the public front door of a repository whose own README
                         invites third parties to take the work and spin it off
  v1.md Chapter 10       the conclusion of a citable ~115-page manuscript

Both stated "1.99x Bliss synergy", "collapses from 3.7% to 0.1%" and "40% ... to
1.8% (CNS/BBB)" as bare numbers under the heading "results that, if validated
experimentally, would have translational implications". README.md was also the
only claim surface in the repository that never linked the audit pricing them.

THE FAILURE MODE THIS GUARDS
----------------------------
Not a wrong number -- every figure is correctly computed. The exposure is a
bystander acting on three significant figures whose interval spans [1.0x, 5.2x].
That is the most probable way this project causes harm, and it costs hours to
close.

WHY IT KEYS ON THE JSON AND NOT ON WORDING
------------------------------------------
tests/test_manuscript_calibration_framing.py records that #589's four-file
relabel drifted back into the manuscript anyway. A guard pinned to a phrase
drifts with the phrase. This one reads the verdicts from
identifiability-report.json, so adding a headline or changing a verdict changes
what is required, and rewording the prose cannot satisfy it.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
IDENT = REPO_ROOT / "analysis" / "identifiability-report.json"
README = REPO_ROOT / "README.md"
MANUSCRIPT = REPO_ROOT / "article" / "drafts" / "v1.md"

# The magnitudes those verdicts are about, and the interval evidence that has to
# sit near each. Keyed by the number as it is written, because that is what a
# reader takes away.
HEADLINE_FIGURES = {
    "1.99": ("Bliss synergy", ("5.2", "1.0")),
    "12.1%": ("penetration gradient", ("93%", "ORDERING", "ordering", "300 of 300")),
    "3.7% to 0.1%": ("hypoxia collapse", ("86.6%", "SONALA-001", "O2-dependent")),
}

# How far from the figure the qualifier may sit and still be read with it.
NEARBY = 700


def _verdicts() -> dict:
    d = json.loads(IDENT.read_text())

    def walk(o, out):
        if isinstance(o, dict):
            if "verdict" in o and isinstance(o["verdict"], str):
                out.append(o["verdict"])
            for v in o.values():
                walk(v, out)
        elif isinstance(o, list):
            for v in o:
                walk(v, out)
    out = []
    walk(d, out)
    return out


def test_the_premise_still_holds():
    """If anything ever becomes point-estimable, this guard must be revisited."""
    v = _verdicts()
    assert v, "identifiability-report.json no longer records per-headline verdicts"
    assert all(x in ("directional_only", "direction_robust_magnitude_not") for x in v), (
        f"a headline now carries a verdict this guard does not know about: {set(v)}. "
        "If something became point-estimable, the surfaces below may legitimately "
        "state it as a magnitude.")


def _near(text: str, needle: str, qualifiers: tuple) -> bool:
    """Does at least one qualifier sit within NEARBY chars of the figure?"""
    for m in re.finditer(re.escape(needle), text):
        window = text[max(0, m.start() - NEARBY): m.end() + NEARBY]
        if any(q in window for q in qualifiers):
            return True
    return False


def test_the_readme_states_intervals_beside_its_headline_numbers():
    txt = README.read_text()
    missing = []
    for fig, (name, quals) in HEADLINE_FIGURES.items():
        if fig in txt and not _near(txt, fig, quals):
            missing.append(f"{name} ({fig})")
    assert not missing, (
        "README.md states these as bare magnitudes with no interval nearby: "
        + ", ".join(missing)
        + ". It is the public front door of a repo that invites others to reuse the work.")


def test_the_readme_links_the_audit_that_prices_its_claims():
    """It was the only claim surface in the repo that never did."""
    txt = README.read_text()
    assert "identifiability-report" in txt, (
        "README.md no longer links analysis/identifiability-report.md, the document "
        "that prices every number it quotes")
    i = txt.index("identifiability-report")
    j = txt.find("**Combination synergy")
    assert j == -1 or i < j, (
        "the pricing disclaimer sits BELOW the results list; a reader meets the "
        "numbers first, which is the arrangement #506 set out to fix")


def test_the_manuscript_conclusion_states_intervals():
    """Chapter 10 is what a citing reader quotes."""
    txt = MANUSCRIPT.read_text()
    start = txt.index("**What the simulations show.**")
    end = txt.find("\n## ", start)
    section = txt[start: end if end > 0 else start + 12000]
    # PROXIMITY, not whole-section membership. A section-wide `in` check passed a
    # mutation that gutted the penetration caveat, because an unrelated fix
    # elsewhere in the same section still contained one of the qualifier strings.
    missing = []
    for fig, (name, quals) in HEADLINE_FIGURES.items():
        if fig in section and not _near(section, fig, quals):
            missing.append(f"{name} ({fig})")
    assert not missing, (
        "the Chapter 10 conclusion states these as bare magnitudes: "
        + ", ".join(missing))
    # NOT an `or`. The first version accepted "none is point-estimable" OR "not
    # point-estimable", and the Bliss sentence a few lines above contains the
    # latter -- so deleting the governing sentence entirely still passed. Third
    # time an `or` in one of these guards has been satisfied by a different fix.
    assert "none is point-estimable" in section, (
        "the conclusion no longer tells the reader that NONE of these magnitudes "
        "is point-estimable, which is the finding that governs all three")


def test_the_hypoxia_qualifier_names_the_right_quantity():
    """The obvious paraphrase of this one is WRONG, and was proposed.

    "3.7% to 0.1%" is the RSL3 kill collapse. What collapses to ~0% under a fully
    O2-dependent SDT is the SDT hypoxic-zone ADVANTAGE -- a different quantity. A
    fix that attached the collapse-to-zero caveat to the RSL3 number would put a
    new error on the front door while appearing to add rigour.
    """
    for path in (README, MANUSCRIPT):
        txt = path.read_text()
        for m in re.finditer(r"SDT hypoxic-zone advantage", txt):
            window = txt[m.start(): m.start() + 400]
            assert "0% to 86.6%" in window or "86.6" in window, (
                f"{path.name}: the SDT advantage is described without its bracket")
