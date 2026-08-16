"""The README's scope statement must equal the counts, not describe them (#731).

WHY THIS EXISTS
---------------
The README now states that the project's analytical, predictive and simulation
work is almost entirely ferroptosis, with a table of counts. That statement is
the answer to a question a reader would otherwise have to answer by counting
files, so it has to stay true as the repository grows.

A hand-written scope paragraph is exactly the shape this repo has been burned by
before: a document that computes some figures and hand-writes the sentence
beside them reads as though the whole line was measured. So the counts are
derived by `scripts/scope_audit.py` and pinned here against the README.

THE FAILURE THIS GUARDS AGAINST, specifically. Someone adds the first
immunotherapy analysis, or the first non-ferroptosis prediction. That is exactly
the outcome the scope statement exists to invite -- and it would silently make
the README wrong in the direction of understating the project's breadth. This
fails at that moment and says so.

WHAT IT DOES NOT GUARD. The bucketing of an individual analysis is a judgement
applied by a stated rule, and reasonable people will disagree about placements.
The artifact lists every member for that reason. What is pinned here is that the
README agrees with whatever the rule produced, and that the two MECHANICAL
counts -- predictions and engine modules -- are right.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "scope_audit.py"
JSON_OUT = REPO_ROOT / "analysis" / "scope-audit.json"
MD = REPO_ROOT / "analysis" / "scope-audit.md"
README = REPO_ROOT / "README.md"


def _doc():
    return json.loads(JSON_OUT.read_text())


def test_the_readme_states_the_scope_at_all():
    """The whole point: a reader must not have to count files."""
    r = README.read_text()
    assert "What the work is actually about" in r, (
        "the README no longer states what the work is about, so the scope is "
        "back to being discoverable only by counting files")
    assert "preregistered predictions" in r and "simulation engine modules" in r


def test_the_readme_counts_equal_the_derived_counts():
    """Pinned to the artifact, so the paragraph cannot drift from the repo."""
    d, r = _doc(), README.read_text()
    nf = len(d["analyses"]["ferroptosis-or-physical"])
    nt = d["n_therapy_subject"]
    nm = len(d["analyses"]["method"])
    pf, np_ = d["n_predictions_ferroptosis"], d["n_predictions"]
    em = d["engine_modules"]
    row = f"| committed analyses | {nf} | **{nt}** | {nm} |"
    assert row in r, f"the README's analysis row is not the derived one: {row}"
    assert f"**{pf} of {np_}**" in r, (
        f"the README does not state the derived prediction count {pf} of {np_}")
    assert f"**{em} of {em}**" in r, (
        f"the README does not state the derived module count {em}")


def test_the_mechanical_counts_are_right():
    """Predictions and modules need no judgement, so they get checked directly."""
    d = _doc()
    preds = re.findall(r"^\*\*(P\d+)\.", (REPO_ROOT / "PREREGISTRATION.md").read_text(), re.M)
    assert len(preds) == d["n_predictions"], (
        f"the audit parsed {d['n_predictions']} predictions, "
        f"PREREGISTRATION.md states {len(preds)}")
    src = REPO_ROOT / "simulations" / "ferroptosis-core" / "src"
    assert d["engine_modules"] == len(list(src.glob("*.rs"))), (
        "the engine module count is not the number of .rs files")


def test_the_ferroptosis_matcher_catches_the_adjectival_form():
    """`ferroptos` does not match `ferroptotic`, and that misfiled P5.

    The error ran in the flattering direction -- it made the project's
    commitments look broader than they are -- which is the direction an audit of
    one's own narrowness must never be wrong in.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("scope", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert m.FERRO.search("dense ferroptotic kill"), (
        "the matcher misses the adjectival form, so any prediction phrased "
        "'ferroptotic' is counted as non-ferroptosis")
    assert m.FERRO.search("ferroptosis induction")
    assert not m.FERRO.search("checkpoint blockade")


def test_a_non_ferroptosis_prediction_would_break_this():
    """The guard must fail when the project widens, which is the good outcome.

    Constructed rather than asserted: if every prediction is ferroptosis today,
    the count is only meaningful if a non-ferroptosis one would change it.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("scope", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    d = _doc()
    assert d["n_predictions_ferroptosis"] == d["n_predictions"], (
        "not every prediction is ferroptosis any more -- which is progress. "
        "Update the README table and this assertion together.")
    # and prove the classifier can say no
    assert not m.FERRO.search(
        "Anti-PD-1 response depends on tumour mutational burden"), (
        "the classifier calls a checkpoint-immunotherapy prediction "
        "ferroptosis, so 'all 8 of 8' is not evidence of anything")


def test_the_method_bucket_is_checked_before_the_vocabulary():
    """Census analyses use FSP1 as their worked example and are not biology."""
    src = SCRIPT.read_text()
    assert "METHOD_STEM.match(stem)" in src, (
        "the instrument-subject test is gone; a body-text scan files "
        "atlas-citation-audit and friends as ferroptosis biology, which "
        "inflates the ferroptosis count in the direction of this file's own "
        "argument")
    d = _doc()
    for stem in ("atlas-citation-audit", "atlas-coverage"):
        assert stem in d["analyses"]["method"], (
            f"{stem} is not in the method bucket; it measures the instrument, "
            "not ferroptosis biology")


def test_the_buckets_are_published_so_a_placement_can_be_disputed():
    """A total nobody can audit is an assertion with a number attached."""
    d, md = _doc(), MD.read_text()
    assert d["analyses"]["therapy-subject"], (
        "no analysis is classified as another therapy's subject; if that is "
        "genuinely true the README claim of 'exactly one' is wrong")
    for stem in d["analyses"]["therapy-subject"]:
        assert stem in md, f"{stem} is counted but not listed in the report"
    assert "can be disputed" in md
