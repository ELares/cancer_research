"""The preregistration may not claim it preceded the calibration work.

WHAT WAS WRONG
--------------
`analysis/osf-registration-package.md` said the predictions were "registered
BEFORE fitting any layer to the external data", and PREREGISTRATION.md said
"Registering before calibration is the point." Neither is true. Four calibration
legs were committed before PREREGISTRATION.md itself:

    kill-switch-calibration.md        2026-06-06   (the CTRPv2 fit, RMSE 0.0504)
    spheroid-structure-validation.md  2026-06-06
    pk-calibration.md                 2026-06-07
    penetration-validation.md         2026-06-07
    PREREGISTRATION.md                2026-06-13

This matters more than an ordinary stale sentence because the package is meant
for an OSF DOI, which is immutable, and the repository is public — so a reviewer
can check `git log` and find the claim false in a minute. It is the single line
someone would use to discredit the whole registration.

The defensible claim, which the package already made in its own section 6, is
DISJOINTNESS: none of those datasets sets any prediction. That is true, it is
what actually protects calibrated-versus-predicted, and it does not depend on an
ordering that did not happen.

WHY THIS GUARD READS GIT
------------------------
A test that hardcoded "2026-06-06" would be asserting a fact about a fact. It
derives the dates from `git log --diff-filter=A`, so if history is ever rewritten
or a leg is re-added, the guard follows rather than lying.
"""

import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

PREREG = REPO_ROOT / "PREREGISTRATION.md"
PACKAGE = REPO_ROOT / "analysis" / "osf-registration-package.md"
PLAN = REPO_ROOT / "analysis" / "contribution-plan-2026.md"
MANUSCRIPT = REPO_ROOT / "article" / "drafts" / "v1.md"

# Legs the package names as the external-data fits it claimed to precede.
CALIBRATION_LEGS = [
    "analysis/calibration/kill-switch-calibration.md",
    "analysis/calibration/spheroid-structure-validation.md",
    "analysis/calibration/pk-calibration.md",
    "analysis/calibration/penetration-validation.md",
]

# A claim of PRECEDENCE. Deliberately does not match the retraction sentences,
# which quote the old wording in order to correct it -- so those are anchored on
# their own distinctive phrasing rather than excluded by a fragile substring.
PRECEDENCE = re.compile(
    r"registered\s+BEFORE\s+fitting"
    r"|Registering\s+before\s+calibration\s+is\s+the\s+point"
    r"|register(ed)?\s+BEFORE\s+(any\s+)?calibration"
    r"|\*?directional,\s*pre-calibration\*?\s+predictions"
    r"|to be registered before any calibration",
    re.I)

RETRACTION_MARKERS = ("which claimed registration came before",
                      "came BEFORE fitting those layers")


def _history_is_shallow() -> bool:
    """CI checks out with fetch-depth: 1, where add-dates do not exist.

    On a shallow clone every file looks "added" at the checkout commit, so
    `git log --diff-filter=A` returns today's date for everything and the
    premise check below reports the opposite of the truth. That is exactly how
    this guard first failed in CI while passing locally. The PROSE checks need
    no history and run everywhere; only the premise is skipped, loudly.
    """
    r = subprocess.run(["git", "rev-parse", "--is-shallow-repository"],
                       cwd=REPO_ROOT, capture_output=True, text=True)
    return r.stdout.strip() == "true"


def _added(path: str) -> str:
    """The date a file first entered the repository, from git."""
    out = subprocess.run(
        ["git", "log", "--diff-filter=A", "--format=%cd", "--date=short", "--", path],
        cwd=REPO_ROOT, capture_output=True, text=True).stdout.split()
    assert out, f"{path} has no add-commit in git history"
    return out[-1]


def _flat(p: Path) -> str:
    """File text with whitespace collapsed and markdown quote markers dropped.

    Content assertions must not depend on where a line happens to wrap. The
    first version of this guard looked for "independently of" and failed on
    prose that reads exactly that way but wraps as "independently\n> of" -- a
    true claim reported as missing, which is the failure mode that trains people
    to weaken guards.
    """
    return re.sub(r"\s+", " ", p.read_text().replace("\n>", " "))


def _claim_lines(p: Path) -> list:
    out = []
    for i, line in enumerate(p.read_text().splitlines(), 1):
        if any(m in line for m in RETRACTION_MARKERS):
            continue
        if PRECEDENCE.search(line):
            out.append(f"{p.relative_to(REPO_ROOT)}:{i}: {line.strip()[:110]}")
    return out


def test_the_precedence_claim_is_false_and_git_still_says_so():
    """The premise. If a leg ever genuinely predated nothing, revisit the prose."""
    if _history_is_shallow():
        pytest.skip("shallow clone: add-dates are unavailable, so the premise "
                    "cannot be checked here. The prose guards still run.")
    prereg = _added("PREREGISTRATION.md")
    earlier = {leg: _added(leg) for leg in CALIBRATION_LEGS
               if _added(leg) < prereg}
    assert earlier, (
        f"No calibration leg now predates PREREGISTRATION.md ({prereg}). The "
        "history changed; the corrected wording in these documents should be "
        "re-checked rather than left asserting a retraction that no longer applies.")


def test_no_document_claims_registration_preceded_calibration():
    offenders = []
    for p in (PREREG, PACKAGE, PLAN, MANUSCRIPT):
        offenders += _claim_lines(p)
    assert not offenders, (
        "these assert the preregistration came before the calibration work, "
        "which git contradicts:\n  " + "\n  ".join(offenders))


def test_the_disjointness_claim_is_stated_instead():
    """Removing the false claim is not enough; the true one has to be there."""
    pkg = _flat(PACKAGE)
    assert "independently of" in pkg, (
        "the package no longer states that the predictions are independent of "
        "the calibration legs, which is the property that actually protects "
        "calibrated-versus-predicted")
    assert "NOT used to set these predictions" in pkg, (
        "section 6's disjointness statement is gone; it is the load-bearing one")


def _ledger_bullets() -> list:
    """The CLASSIFICATION bullets of Part 3, not its surrounding prose.

    A whole-section search is not good enough, proven by mutation: the section
    also contains a sentence explaining that an earlier version "silently
    omitted P3 and P8", so deleting P3 from the actual classification left the
    guard finding "P3" in the retraction that describes the bug. The guard was
    satisfied by its own explanation of what it was guarding against.
    """
    txt = PREREG.read_text()
    part3 = txt[txt.index("## Part 3"):txt.index("## Literature position")]
    # Each bullet ends at the first blank line. Without that cut the LAST
    # bullet absorbs every following paragraph -- including the note saying an
    # earlier version "omitted P3 and P8" -- so deleting P3 from the real
    # classification still left the guard finding it. Proven by mutation twice:
    # the same self-satisfying-prose defect, one level further in.
    out = []
    for b in part3.split("\n- **")[1:]:
        out.append(re.sub(r"\s+", " ", b.split("\n\n")[0]))
    return out


def test_the_ledger_classifies_every_prediction():
    """It listed six of eight, dropping P3 and P8 without saying so."""
    defined = set(re.findall(r"^\*\*(P[1-8])\.", PREREG.read_text(), re.M))
    assert len(defined) == 8, f"expected P1-P8, found {sorted(defined)}"
    bullets = _ledger_bullets()
    assert bullets, "Part 3 has no classification bullets; its shape changed"
    classified = {p for p in defined
                  if any(re.search(rf"\b{p}\b", b) for b in bullets)}
    missing = sorted(defined - classified)
    assert not missing, (
        f"Part 3's classification bullets omit {missing}. They are defined in "
        "Part 1, so a reader cannot tell whether they are prior-predictive, "
        "anchored, or derived from an artifact.")


def test_P8_is_flagged_as_derived_from_a_calibration_artifact():
    """P8 reproduces numbers from a file committed the same day it was written.

    That does not make it unfalsifiable, but it does mean it must never be
    counted as an independent hit if the model later agrees with it.

    Anchored on P8's OWN bullet rather than the section, because the section
    mentions the artifact more than once and a mutation removing it from the
    classification survived a whole-section search.
    """
    p8 = [b for b in _ledger_bullets() if re.search(r"\bP8\b", b)]
    assert p8, "no ledger bullet classifies P8"
    bullet = " ".join(p8)
    assert "spheroid-kill-vs-size" in bullet, (
        "P8's bullet does not name the artifact its numbers come from")
    assert "independent hit" in bullet, (
        "P8's bullet does not say what its provenance means for interpreting it")
    # The date that justifies the flag, from git rather than from memory. The
    # artifact must not POSTDATE the preregistration -- if it ever did, P8 would
    # genuinely have preceded it and this classification would be wrong.
    if _history_is_shallow():
        return
    assert _added("analysis/calibration/spheroid-kill-vs-size.md") <= _added("PREREGISTRATION.md"), (
        "spheroid-kill-vs-size.md now postdates PREREGISTRATION.md, so P8 no "
        "longer reproduces a pre-existing artifact and its ledger entry is stale")
