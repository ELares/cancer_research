"""Guards for the measured #617 co-mention regression (#ATLAS-COMENT-REG).

The finding runs against a change this repo made and justified, so it is exactly
the kind of result that quietly reverts. These pin it to recomputed quantities.
"""

import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from comention_regression import wilson  # noqa: E402

RAW = REPO_ROOT / "analysis" / "comention-regression.json"
DOC = REPO_ROOT / "analysis" / "comention-regression.md"
JUDGED = REPO_ROOT / "analysis" / "comention" / "abstract-visible-judgements.csv"
COMENTION = REPO_ROOT / "scripts" / "atlas_comention.py"
GRAPH = REPO_ROOT / "scripts" / "atlas_graph.py"


def _raw():
    return json.loads(RAW.read_text())


def test_the_judgements_are_committed_and_match_the_reported_precision():
    """The precision must be re-derivable from the judged rows, not asserted."""
    rows = list(csv.DictReader(JUDGED.open()))
    assert len(rows) >= 50, f"only {len(rows)} judged mentions"
    assert {r["verdict"] for r in rows} <= {"TP", "FP"}
    # The reported number uses the identifier-level verdicts; `verdict` is kept
    # as the superseded first pass so the correction stays auditable.
    tp = sum(1 for r in rows if r["verdict_v2"] == "TP")
    d = _raw()["after"]
    assert d["judged_n"] == len(rows) and d["judged_tp"] == tp
    assert abs(d["abstract_precision"] - tp / len(rows)) < 1e-12
    lo, hi = wilson(tp, len(rows))
    assert abs(d["abstract_precision_ci"][0] - lo) < 1e-12


def test_precision_fell_and_the_document_says_so():
    """The finding. If a rebuild recovers it, this must fail and be re-read."""
    d = _raw()
    assert d["net_change"] < 0, "precision no longer fell; re-read the finding"
    assert d["before"]["weighted"] > d["after"]["weighted"]
    txt = DOC.read_text()
    assert "**No." in txt, "the document no longer leads with the negative answer"


def test_the_abstract_visible_stratum_grew():
    """The mechanism of the loss: the stratum got cleaner but much larger."""
    d = _raw()
    before, after = d["before"]["strata"], d["after"]["strata"]
    assert after["abstract"] > before["abstract"] * 2, (
        "the stratum no longer tripled; the explanation in the document is stale")
    # And it did get cleaner, which is why the loss is not obvious.
    assert d["after"]["abstract_precision"] > d["before"]["abstract_precision"]


def test_the_carried_over_precisions_are_declared_not_hidden():
    """Only one stratum was re-measured; the other two are assumptions."""
    d = _raw()
    assert set(d["carried_over"]) == {"agree", "body_only"}
    txt = DOC.read_text()
    assert "carried over" in txt and "stated rather than buried" in txt


def test_the_share_bug_is_fixed_at_the_source():
    """The ratio must use a per-(form, identifier) count, not the cross-sense sum.

    Dividing `alias_support` by `ident_mentions` is not a share and can exceed 1
    (274% for `as`), and it admits ambiguous generic words MORE readily than
    specific names -- the opposite of the filter's intent.
    """
    graph = GRAPH.read_text()
    assert "alias_ident_support" in graph, "the corrected numerator is not recorded"
    com = COMENTION.read_text()
    assert "ident_support.get(a," in com, "the share still uses the cross-sense total"
    # The old expression must not come back.
    assert "share = support.get(a, 0) / max(1, ident_tot" not in com


def test_the_lesson_is_recorded_not_just_the_number():
    """A filter justified by an error distribution it then changes must be
    re-measured. That is the transferable part."""
    txt = DOC.read_text()
    assert "wrong about the population it was applied to" in txt
    assert "unblinded" in txt, "the judgement's own bias is not disclosed"


def test_the_baseline_is_pinned_with_provenance_not_read_from_head():
    """Reading the 'before' run from HEAD compared the audit against itself.

    Once the post-rebuild audit was committed, HEAD held the new numbers, so the
    script silently refused to produce anything. The historical run cannot
    change, so it is pinned to the commit it came from and verified against git.
    """
    import comention_regression as cr
    assert len(cr.PRE_REBUILD_COMMIT) == 40
    assert cr.PRE_REBUILD_AUDIT["mentions"] == 1112
    # The verification must be live, not a comment.
    src = (REPO_ROOT / "scripts" / "comention_regression.py").read_text()
    assert "the comparison baseline has moved" in src
    assert "HEAD:analysis/atlas-comention-audit.json" not in src


def test_the_source_fix_is_reported_as_insufficient():
    """The share bug is real, and fixing it does not recover the precision.

    Claiming a repair the measurement does not support is the failure mode this
    whole document is about.
    """
    txt = DOC.read_text()
    assert "It is not enough" in txt
    assert "0 of 37" in txt, "the inert-filter measurement is missing"


def test_judging_is_at_the_identifier_level_not_the_surface_form():
    """The correction that changed the headline.

    The first pass asked whether the sentence contained the matched string. The
    identifier is what the co-mention pair is recorded against, so the right
    question is whether the sentence discusses what that identifier DENOTES.
    `apoptosis` resolving to Malformations of Cortical Development scored
    correct under the first criterion and is plainly wrong.
    """
    rows = list(csv.DictReader(JUDGED.open()))
    assert "verdict_v2" in rows[0], "identifier-level verdicts are missing"
    assert "authority_name" in rows[0], "the NLM label is not attached, so the "\
        "judgement cannot be checked by a reader"
    v2 = sum(1 for r in rows if r["verdict_v2"] == "TP")
    v1 = sum(1 for r in rows if r["verdict"] == "TP")
    assert v2 < v1, "the stricter criterion should not be more permissive"
    d = _raw()["after"]
    assert abs(d["abstract_precision"] - v2 / len(rows)) < 1e-12, \
        "the reported precision does not use the identifier-level verdicts"
    # Spot-check the case that exposed the flaw.
    apop = [r for r in rows if r["surface_form"].lower() == "apoptosis"]
    assert apop and apop[0]["verdict_v2"] == "FP", \
        "apoptosis -> Malformations of Cortical Development is not a true positive"


def test_the_authority_discriminator_is_measured_and_not_yet_recommended():
    """#628's first acceptance criterion, and its honest scope.

    It must beat the filters it replaces on the judged sample, and it must be
    labelled as selected-on-this-sample rather than validated.
    """
    rows = list(csv.DictReader(JUDGED.open()))
    kept = [r for r in rows if r["is_authority_name"] == "True"]
    assert kept, "the discriminator column is missing or never fires"
    kept_tp = sum(1 for r in kept if r["verdict_v2"] == "TP")
    all_tp = sum(1 for r in rows if r["verdict_v2"] == "TP")
    all_fp = len(rows) - all_tp
    cut_fp = all_fp - (len(kept) - kept_tp)
    # It must remove the great majority of false positives -- the support and
    # share filters removed none, which is the comparison that matters.
    assert cut_fp / all_fp > 0.8, f"only {100*cut_fp/all_fp:.0f}% of FPs removed"
    assert kept_tp / len(kept) > all_tp / len(rows), "precision did not improve"
    txt = DOC.read_text()
    assert "NOT yet a recommendation" in txt
    assert "selected on this sample" in txt
