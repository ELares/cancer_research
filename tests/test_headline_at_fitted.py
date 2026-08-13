"""Guards for the headlines-at-the-fitted-cascade comparison.

WHAT THIS ANALYSIS IS
---------------------
`analysis/calibration/in-vivo-prior-provenance.md` named one outstanding piece of
work: re-derive every headline at the CTRPv2-fitted cascade and report both, so a
reader can see which directions survive crossing the bistable tipping point. It
called that "achievable now and the cheaper of the two". It had not been done.

WHAT IT FOUND, AND THE TRAP IN REPORTING IT
-------------------------------------------
Substituting the fitted cascade drives the in-vivo and spatial models into a
regime where UNTREATED cells die at ~100%, against the model's own stated
constraint of under 2%. Every headline then degenerates: the Bliss ratio is
exactly 1.0 because both single arms saturate, and all three penetration tissues
kill 100% so their ordering is trivially "preserved".

The trap is that each of those reads, on its own, like a headline failing. It is
not. A ratio between two saturated arms is arithmetic, not biology. So the
analysis must establish ADMISSIBILITY first and phrase every downstream verdict
conditionally, and these guards exist to keep it that way.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JSON = REPO_ROOT / "analysis" / "headline-at-fitted-cascade.json"
DOC = REPO_ROOT / "analysis" / "headline-at-fitted-cascade.md"
MANUSCRIPT = REPO_ROOT / "article" / "drafts" / "v1.md"
SCRIPT = REPO_ROOT / "scripts" / "headline_at_fitted.py"


def results() -> dict:
    return json.loads(JSON.read_text())


def test_the_baseline_constraint_matches_the_manuscript():
    """The 2% threshold must come from the manuscript, not from a constant.

    The analysis calls a parameter set inadmissible for breaking a constraint
    Chapter 5 states in prose. If the chapter ever changes that number, a
    hardcoded copy would keep judging against the old one.
    """
    m = re.search(r"all phenotypes show less than (\d+(?:\.\d+)?)% death under Control",
                  MANUSCRIPT.read_text())
    assert m, ("Chapter 5 no longer states the untreated-death constraint in the "
               "form this analysis relies on; re-check scripts/headline_at_fitted.py")
    stated = float(m.group(1)) / 100.0
    src = SCRIPT.read_text()
    used = float(re.search(r"BASELINE_MAX = ([\d.]+)", src).group(1))
    assert abs(used - stated) < 1e-9, (
        f"the manuscript states {stated:.3f} but the analysis judges against {used:.3f}")


def test_the_fitted_sets_are_recorded_as_inadmissible():
    """The load-bearing result. If this ever flips, the write-up is wrong."""
    r = results()
    assert r["default"]["admissibility"]["admissible"], (
        "the DEFAULT parameter set now violates the untreated-death constraint, "
        "which would be a far bigger finding than this analysis reports")
    for s in ("ctrpv2_point", "posterior_median"):
        a = r[s]["admissibility"]
        assert not a["admissible"], (
            f"{s} is now admissible; the document's central negative finding "
            "no longer holds and its prose must be rewritten")
        assert a["worst_rate"] > 0.5, (
            f"{s} untreated death is {a['worst_rate']:.3f}; the document says "
            "untreated cells die en masse, which needs a rate near 1")


def test_the_admissibility_verdict_is_computed_from_the_rates():
    """The verdict must follow the data, not be asserted beside it."""
    for name, r in results().items():
        controls = {c: v for c, v in r["single_cell"].items()
                    if c.endswith("/Control")}
        assert controls, f"{name} has no Control conditions to judge"
        worst = max(controls.values())
        assert abs(worst - r["admissibility"]["worst_rate"]) < 1e-12, (
            f"{name}: recorded worst untreated rate does not match its own "
            "single-cell table")
        assert r["admissibility"]["admissible"] == (worst <= r["admissibility"]["constraint"])


def test_the_document_does_not_report_degenerate_arithmetic_as_a_failed_direction():
    """The trap this whole analysis is shaped to avoid.

    With inadmissible sets present, every headline verdict must be marked as not
    a verdict. Otherwise "supra-additive in 1 of 3 sets" reads as evidence the
    synergy claim collapsed, when the set producing 1.0 had already been shown
    to be a model whose untreated cells are dead.
    """
    txt = DOC.read_text()
    any_bad = any(not r["admissibility"]["admissible"] for r in results().values())
    if not any_bad:
        return
    assert "Not a verdict on the headline" in txt, (
        "inadmissible parameter sets are present but the headline sections do "
        "not mark their numbers as non-verdicts")
    # the caveat must appear beside EVERY headline, not once at the top
    for section in ("Bliss synergy", "Penetration gradient"):
        i = txt.index(f"## {section}")
        j = txt.find("\n## ", i + 1)
        assert "Not a verdict on the headline" in txt[i:j if j > 0 else len(txt)], (
            f"the '{section}' section reports numbers from an inadmissible "
            "parameter set without saying so")


def test_the_bliss_collapse_is_explained_as_saturation():
    """1.0 exactly is a saturation signature, and must be named as one."""
    r = results()
    if r["ctrpv2_point"]["bliss"] != 1.0:
        return
    txt = DOC.read_text()
    assert "saturate" in txt, (
        "the Bliss ratio is exactly 1.0 at the fitted set — the value a Bliss "
        "ratio takes when both single arms already kill everything — and the "
        "document does not say so, so a reader will read it as 'no synergy'")


def test_the_analysis_does_not_claim_to_have_calibrated_anything():
    """An in-vitro posterior in an in-vivo model is not a calibration."""
    txt = DOC.read_text().lower()
    assert "does not" in txt and "data-conditioned" in txt, (
        "the document no longer states that this does not make any spatial "
        "headline data-conditioned")
    for forbidden in ("now calibrated", "is calibrated to", "validates the"):
        assert forbidden not in txt, f"overclaim present: {forbidden!r}"
