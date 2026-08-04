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
    # It must NOT claim the stratum got cleaner. The two figures were produced
    # under different criteria (surface-form then identifier-level), so the
    # document says only that it did not measurably improve.
    txt = DOC.read_text()
    assert "did not measurably get cleaner" in txt


def test_every_stratum_is_measured_not_assumed():
    """The two carried-over precisions were the biggest hole in the total.

    They were assumed at 92.5% and 30.8% from #617. Measuring them mattered:
    body-only came in at 20%, optimistic by about a third, and it is the stratum
    carrying the most volume.
    """
    import csv as _csv

    d = _raw()
    m = d["measured_strata"]
    assert set(m) == {"agree", "body_only"}
    for name, path in (("agree", "corroborated-judgements.csv"),
                       ("body_only", "body-only-judgements.csv")):
        f = REPO_ROOT / "analysis" / "comention" / path
        rows = list(_csv.DictReader(f.open()))
        assert len(rows) >= 30, f"{path} has only {len(rows)} judged mentions"
        tp = sum(1 for r in rows if r["verdict"] == "TP")
        assert m[name]["tp"] == tp and m[name]["n"] == len(rows)
        assert abs(m[name]["precision"] - tp / len(rows)) < 1e-12
        # Each judgement must be checkable. csv.DictReader supplies every header
        # key on every row, so `"matched_span" in r` tests the HEADER -- blanking
        # every value passed it. Check the values.
        filled = sum(1 for r in rows if (r.get("matched_span") or "").strip())
        assert filled >= 0.9 * len(rows), (
            f"{name}: only {filled}/{len(rows)} rows record the span that fired")
        assert all(r.get("verdict") in ("TP", "FP") for r in rows)
    assert m["body_only"]["precision"] < 0.308, (
        "body-only no longer measures below its old assumed value; re-read")
    assert "All three rows are now hand-judged" in DOC.read_text()


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
    assert 0 < v2, "verdict_v2 is all false positives; that is not a judgement"
    d = _raw()["after"]
    assert abs(d["abstract_precision"] - v2 / len(rows)) < 1e-12, \
        "the reported precision does not use the identifier-level verdicts"
    # Spot-check the case that exposed the flaw.
    apop = [r for r in rows if r["surface_form"].lower() == "apoptosis"]
    assert apop and apop[0]["verdict_v2"] == "FP", \
        "apoptosis -> Malformations of Cortical Development is not a true positive"


def test_the_discriminator_column_is_derivable_not_hand_entered():
    """Nothing in the pipeline generates `is_authority_name`, so it could be

    edited to anything without a failure. Re-derive it here from the two columns
    it depends on, so a hand-edit cannot slip through -- the flagship
    counterexample (`apoptosis` against Malformations of Cortical Development)
    must never be marked a name match.
    """
    import re as _re

    def bag(s):
        return frozenset(w for w in _re.split(r"[^a-z0-9]+", (s or "").lower()) if w)

    rows = list(csv.DictReader(JUDGED.open()))
    for r in rows:
        expected = bool(r["authority_name"]) and bag(r["surface_form"]) == bag(r["authority_name"])
        assert (r["is_authority_name"] == "True") == expected, (
            f"row {r['n']} ({r['surface_form']}): column says "
            f"{r['is_authority_name']}, the rule gives {expected}")


def test_the_discriminator_is_not_degenerate():
    """A filter that keeps almost nothing scores perfectly and is worthless.

    Removing 100% of false positives by keeping one mention would satisfy a
    naive precision guard, so bound the recall side too.
    """
    rows = list(csv.DictReader(JUDGED.open()))
    kept = [r for r in rows if r["is_authority_name"] == "True"]
    all_tp = sum(1 for r in rows if r["verdict_v2"] == "TP")
    kept_tp = sum(1 for r in kept if r["verdict_v2"] == "TP")
    assert kept_tp >= 3, f"only {kept_tp} true matches survive; the filter is degenerate"
    assert kept_tp / all_tp >= 0.3, (
        f"keeps only {100*kept_tp/all_tp:.0f}% of true matches; report that as a "
        "cost rather than a precision win")
    txt = DOC.read_text()
    assert "It removes every gene" in txt, (
        "the gene blind spot is not disclosed; NCBI Gene ids have no MeSH label "
        "so the rule cuts all of them")


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
    assert "NOT a recommendation" in txt
    assert "selected on this sample" in txt
    # And it must describe what it actually is: an identifier-level check on the
    # canonical corpus name, not a per-mention filter on the span that matched.
    assert "IDENTIFIER-LEVEL check, not a per-mention one" in txt, (
        "the discriminator is described as filtering the matched span, which the "
        "audit sample does not record")


def test_the_regression_survives_a_shared_judging_flaw():
    """The comparison's weakest point, worked out rather than hedged.

    If the carried-over precisions share the surface-form flaw, the `before`
    total is overstated too. Scaling both sides' impure strata by the same
    factor k, the change must stay negative for every k in [0,1] -- otherwise
    the regression could be an artifact of correcting only one side.
    """
    import comention_regression as cr

    d = _raw()
    sb, sa = d["before"]["strata"], d["after"]["strata"]
    abs_after = d["after"]["abstract_precision"]
    agree_p = d["measured_strata"]["agree"]["precision"]
    body_p = d["measured_strata"]["body_only"]["precision"]
    deltas = []
    for i in range(11):
        k = i / 10
        B = (sb["agree"] * agree_p
             + k * (sb["abstract"] * cr.PRIOR_ABSTRACT_PRECISION
                    + sb["body"] * body_p))
        A = (sa["agree"] * agree_p
             + k * (sa["abstract"] * abs_after + sa["body"] * body_p))
        deltas.append(A - B)
    assert all(x < 0 for x in deltas), (
        f"the regression reverses at some k: {[round(x,4) for x in deltas]}")
    # And it should widen as k falls -- the `after` run carries more weight in
    # the impure strata, so correcting both sides cannot rescue it.
    assert deltas[0] < deltas[-1], "the gap no longer widens under correction"
    assert "WIDENS the gap" in DOC.read_text()


def test_the_audit_does_not_assert_both_readings_of_the_body_only_stratum():
    """It once said body-only was "the layer doing the job it exists for" in one
    section and "NOT the layer doing its job" in another, three screens apart.

    The second is right. A document holding both is worse than either, because a
    reader takes whichever they reach first.
    """
    md = (REPO_ROOT / "analysis" / "atlas-comention-audit.md").read_text()
    assert "asserted BOTH" in md, "the self-contradiction is no longer acknowledged"
    # The two measurements of this stratum come from different runs AND
    # different criteria, so the document must not present them as comparable.
    assert "not like-for-like" in md
    assert "PRE-FILTER measurement" in md, (
        "the #617 block is not labelled as the pre-filter run, so its 30.8% reads "
        "as current beside the post-filter 20.0%")


def test_matched_forms_reproduces_the_build_matcher():
    """It must be longest-match with consumption, not n-gram candidate generation.

    The build's `sentence_entities` walks tokens, takes the longest alias at each
    position and skips past it. Returning every n-gram hit instead reports spans
    that never fired -- and those spans go straight in front of the next round of
    hand judging, which is what this field exists to prevent.
    """
    import atlas_comention as ac

    alias = {"breast cancer": "D001943", "cancer": "D009369", "breast": "D001940"}
    s = "Patients with breast cancer were enrolled in the trial."
    assert ac.matched_forms(s, "D001943", alias) == ["breast cancer"]
    # The shorter aliases are consumed by the longer match and must NOT appear.
    assert ac.matched_forms(s, "D009369", alias) == []
    assert ac.matched_forms(s, "D001940", alias) == []

    # It must still be token-based, not substring: `oral` is not in `oropharyngeal`.
    assert ac.matched_forms("oropharyngeal carcinoma", "X", {"oral": "X"}) == []
    assert ac.matched_forms("the oral cavity", "X", {"oral": "X"}) == ["oral"]

    # And it must agree with the build's own matcher on the same input.
    ents = ac.sentence_entities(s, alias)
    fired = {e for e in (ents if isinstance(ents, (set, list)) else [])}
    if fired:
        assert "D001943" in fired and "D009369" not in fired


def test_the_judged_sentences_are_long_enough_to_check_the_span():
    """"Every verdict is checkable" is only true if the span is in the sentence.

    The CSVs truncate, and three rows carried a span that fell past the cut --
    including the row this work claimed to have resolved by recovering it.
    """
    import csv as _csv

    for name in ("abstract-visible-judgements.csv", "body-only-judgements.csv",
                 "corroborated-judgements.csv"):
        path = REPO_ROOT / "analysis" / "comention" / name
        for r in _csv.DictReader(path.open()):
            span = (r.get("matched_span") or "").split("|")[0]
            if not span:
                continue
            assert span.lower() in r["sentence"].lower(), (
                f"{name} row {r['n']}: span {span!r} is not in the committed "
                "sentence, so the verdict cannot be checked by a reader")
