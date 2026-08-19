"""Guards for the qualifier-axis measurement (#722 step 1).

WHAT THIS MEASUREMENT IS FOR
----------------------------
`atlas_baseline.parse_articles` reads `DescriptorName` and never
`QualifierName`, so the census carries one of MeSH's two axes. This measures
what that drops, and it is the gate on whether the expensive re-parse in step 2
is worth doing and for which modalities.

WHY IT IS EASY TO GET WRONG, AND WAS
-------------------------------------
The first attempt at this question compared how often a concept is DECLARED in
title/abstract against how often the descriptor layer catches it, and reported a
2.5x sensitivity gap between ferroptosis and radiotherapy. Two errors:

  ASYMMETRIC RULES. One stem against one exact descriptor label on one side,
  a multi-term regex against a multi-descriptor family on the other. Six
  independent recounts produced six different numerators.

  A CONTROL THAT REFUTED THE EXPLANATION. Angiogenesis has no qualifier form and
  scores descriptor recall of the same low order as radiotherapy's (5.9%
  against 2.7% -- an earlier docstring said LOWER, which those figures
  contradict), so low recall is not
  evidence of a qualifier problem. Most of the spread was topicality.

The design that survives compares the SAME articles parsed the SAME way,
differing in one XML element. That difference cannot be topicality, vocabulary
breadth or era. These guards pin that design, not the numbers it produced.

THE RESULT ALSO INVERTED THE ORIGINAL FRAMING, and the report must keep saying
so: radiotherapy is the WEAKEST case of the four, not the strongest. An issue
written to argue radiotherapy is neglected should not quietly keep leading with
it once its own measurement says surgery is the sharp one.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_ingest_sensitivity.py"
MD = REPO_ROOT / "analysis" / "atlas-ingest-sensitivity.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-ingest-sensitivity.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("sens", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_both_arms_are_defined_for_every_modality():
    """Symmetry by construction, since asymmetry is what broke the last try."""
    m = _mod()
    assert m.MODALITIES, "no modalities defined"
    for name, spec in m.MODALITIES.items():
        assert spec.get("qualifiers"), f"{name} has no qualifier side"
        assert spec.get("descriptors"), (
            f"{name} has no descriptor side, so its 'qualifier-only' figure "
            "would be the qualifier count with nothing subtracted")


def test_membership_is_decided_the_ingest_s_way():
    """Both arms must describe the articles the census actually contains."""
    src = SCRIPT.read_text()
    assert "from atlas_baseline import" in src and "ADJACENT_DESCRIPTORS" in src, (
        "the scan no longer imports the ingest's own cancer-membership set, so "
        "it could be measuring a different article population")
    assert "fetch_c04_descriptors" in src, "the C04 set is not the ingest's"
    assert "QualifierName" in src, (
        "the scan does not read QualifierName, which is the entire subject")


def test_qualifiers_cannot_change_membership():
    """The claim that this is a labelling gap, not a selection gap.

    Selection happens on descriptor UIs before any qualifier is read, so the
    two arms are guaranteed to see one article set. If the code ever selected
    on qualifiers the comparison would silently become a different measurement.
    """
    src = SCRIPT.read_text()
    body = src[src.index("def parse_both_axes"):src.index("def scan(")]
    sel = body[body.index("if not any(u in c04"):]
    assert "quals" not in sel.split("yield")[0], (
        "cancer membership is being decided with qualifiers in scope; "
        "membership must be descriptor-only or the arms are not comparable")
    md = MD.read_text()
    assert "not a selection gap" in md, (
        "the report no longer distinguishes a labelling gap from a selection "
        "gap, which is the distinction that keeps this from reading as 'the "
        "census is missing articles'")


def test_the_marginal_gain_is_derived_and_carries_intervals():
    """A sampled per-article rate reported as a point invites over-reading."""
    d, md = _doc(), MD.read_text()
    n = d["cancer_articles"]
    assert n > 1000, f"only {n} cancer articles sampled; too few to gate on"
    m = _mod()
    for name, s in d["modalities"].items():
        qo = s["qualifier_only"]
        assert qo <= s["qualifier"], (
            f"{name}: qualifier-only exceeds the qualifier count, which is "
            "arithmetically impossible")
        assert s["either"] >= max(s["descriptor"], s["qualifier"]), (
            f"{name}: the union is smaller than one of its arms")
        lo, hi = m.wilson(qo, n)
        assert f"{100*lo:.1f}-{100*hi:.1f}" in md, (
            f"{name}'s interval is not the one its own counts give")


def test_the_report_states_which_modality_is_actually_sharpest():
    """The measurement inverted the framing that motivated it.

    #722 was written arguing radiotherapy is the neglected modality. Its own
    step-1 result says surgery gains most and radiotherapy least. The report has
    to name the real winner, or the analysis becomes decoration on a conclusion
    that preceded it.
    """
    d, md = _doc(), MD.read_text()
    ranked = sorted(d["modalities"].items(),
                    key=lambda kv: -kv[1]["qualifier_only"])
    sharpest, weakest = ranked[0][0], ranked[-1][0]
    assert f"**{sharpest}**" in md, (
        f"the report does not name {sharpest} as the sharpest case")
    # AND the renderer must DERIVE it. The assertion above reads the committed
    # report, which does not move when the generator is edited -- the third time
    # this artifact-versus-code trap has appeared in guards written today. A
    # hardcoded winner is exactly how a measurement becomes decoration on a
    # conclusion that preceded it.
    src = SCRIPT.read_text()
    assert 'key=lambda kv: kv[1].get("qualifier_only", 0)' in src, (
        "the sharpest-case sentence no longer picks the maximum; it is naming "
        "a modality chosen by the author rather than by the data")
    assert sharpest != weakest
    # and the numbers must actually rank that way in the artifact
    n = d["cancer_articles"]
    top = 100 * ranked[0][1]["qualifier_only"] / n
    bot = 100 * ranked[-1][1]["qualifier_only"] / n
    assert top > bot, "the ranking is degenerate"


def test_the_superseded_sensitivity_claim_is_not_revived():
    """The 2.5x ratio failed a control and must not come back unqualified."""
    md = MD.read_text()
    assert "asymmetric rules" in md or "asymmetric" in md, (
        "the report no longer records why the earlier text-vs-descriptor "
        "comparison was withdrawn, so someone will recompute it")
    assert "angiogenesis" in md, (
        "the control that refuted the qualifier explanation is no longer "
        "named; it is the reason this design is per-article rather than "
        "per-concept")


def test_an_empty_scan_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'if d["cancer_articles"] == 0:' in src, (
        "the zero-article check no longer tests the count")
    assert "is not a finding" in src and "raise SystemExit" in src


def test_shards_are_sampled_across_the_range_not_taken_as_a_prefix():
    """Baseline files are chronological; a prefix samples only old literature."""
    src = SCRIPT.read_text()
    assert "rng.sample(" in src, (
        "shards are no longer randomly sampled, so a prefix would measure "
        "1970s MeSH indexing practice and call it the census")
    d = _doc()
    assert len(d["shards"]) == d["n_shards"] and d["n_shards"] >= 4, (
        "too few shards to spread across fifty years of indexing practice")


def test_render_only_works_without_the_raw_xml():
    # Run against a COPY of the tree's artifacts. Invoking the generator here
    # rewrites the committed report, which makes the suite non-idempotent and --
    # measured the hard way -- lets a mutation sweep corrupt a committed file:
    # a mutated renderer ran through this test and its text landed on disk.
    import shutil, tempfile
    with tempfile.TemporaryDirectory() as td:
        shutil.copy2(MD, Path(td) / MD.name)
        try:
            res = subprocess.run([sys.executable, str(SCRIPT), "--render-only"],
                                 cwd=REPO_ROOT, capture_output=True, text=True)
        finally:
            shutil.copy2(Path(td) / MD.name, MD)
    assert res.returncode == 0, (
        f"--render-only failed, so the report cannot be rebuilt without "
        f"re-downloading:\n{res.stdout}\n{res.stderr}")

def test_every_descriptor_proxy_entry_can_actually_fire():
    """`neoplasms/surgery` was a descriptor/QUALIFIER composite in a list
    matched against DescriptorName text, so it could never match. It presented
    the surgery arm as five entries when four were live, and the guard whose
    stated purpose was "otherwise its qualifier-only figure would be the
    qualifier count with nothing subtracted" asserted only non-emptiness --
    reducing the arm to that one dead entry left the suite green.
    """
    m = _mod()
    for name, spec in m.MODALITIES.items():
        live = [d for d in spec["descriptors"] if "/" not in d]
        assert live, f"{name}'s descriptor arm has no entry that can fire"
        dead = sorted(set(spec["descriptors"]) - set(live))
        assert not dead, (
            f"{name} lists {dead}, which contain a qualifier and are matched "
            "against DescriptorName, so they can never fire")
        assert len(spec["qualifiers"]) >= 1


def test_the_cross_modality_ordering_is_withdrawn():
    """The gain is measured against a hand-written proxy for an open tree
    while the qualifier arm is a whole closed axis, so ranking the rows ranks
    proxy completeness as much as it ranks the axis.
    """
    md = MD.read_text()
    assert "ORDERING is not a measurement" in md
    assert "THAT ORDERING IS WITHDRAWN" in md
    assert "ranks that incompleteness" in md
    # the two clauses the ordering rested on must not stand unqualified
    for m_ in re.finditer(r"nowhere the ingest looks", md):
        w = md[max(0, m_.start() - 300):m_.end() + 300]
        assert "WITHDRAWN" in w, (
            "the report still says the qualifier axis carries these articles "
            "and nowhere the ingest looks, which holds only against four "
            "descriptors rather than against the ingest")
    for m_ in re.finditer(r"and it says which those are", md):
        w = md[max(0, m_.start() - 400):m_.end() + 200]
        assert "no materiality threshold" in w


def test_the_control_direction_matches_the_only_figures_the_repo_carries():
    """5.9% against 2.7% is HIGHER, and three sites said LOWER."""
    md, src = MD.read_text(), SCRIPT.read_text()
    tst = Path(__file__).read_text()
    for text, where in ((md, "report"), (src, "generator"), (tst, "guard")):
        for m_ in re.finditer(r"[Ll]ower descriptor recall than radiotherapy",
                              text):
            w = text[max(0, m_.start() - 400):m_.end() + 400]
            assert re.search(r"earlier|contradict|withdraw", w, re.I), (
                f"the {where} says angiogenesis has lower descriptor recall "
                "than radiotherapy, and the only figures the repo carries for "
                "it are 5.9% against 2.7%")
    assert "computed NOWHERE in this repo" in src, (
        "the three figures behind the withdrawn claim are quoted as though "
        "this repo produced them")


def test_the_committed_report_is_what_the_generator_produces():
    """No regenerate-and-diff gate existed, so every renderer edit was
    invisible: a mutation dropping a whole bullet passed the suite.
    """
    m = _mod()
    assert m.render(_doc()) == MD.read_text(), (
        "analysis/atlas-ingest-sensitivity.md is not what the current "
        "renderer produces from the committed JSON -- re-run with "
        "--render-only")

