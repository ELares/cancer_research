"""Guards for `analysis/modality-module-depth.{md,json}`.

This page exists because a sentence in Chapter 6 was a measurement written as
an assertion, was true when written, and stopped being true silently. So the
guards are about the two ways it could go wrong again: the count drifting from
the crate, and the page being read as progress.
"""
import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/modality-module-depth.md"
JSON_ = REPO / "analysis/modality-module-depth.json"
CORE = REPO / "simulations/ferroptosis-core/src"
MANUSCRIPT = REPO / "article/drafts/v1.md"


def _load():
    spec = importlib.util.spec_from_file_location(
        "modality_module_depth", REPO / "scripts/modality_module_depth.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MD_ = _load()

# The classification, as a decision record. A module must be added to one of
# these deliberately; see `test_every_module_is_classified_by_a_decision`.
# `photosensitizer_pk` MOVED, AND THE MOVE IS THE LESSON THIS FILE ALMOST
# TAUGHT ITSELF. It is photodynamic therapy's own file and has been since long
# before any arm below existed, and it sat in ENGINE_PIN -- so its lines were
# counted on the COMPARATOR'S side of a ratio measuring distance from the
# comparator, while `arm_parity.py` credited the same lines to PDT. The guard
# below exists for precisely this defect and its docstring describes it
# exactly; it did not fire, because the error was already inside the pin when
# the pin was written. A pin catches a NEW mistake and freezes an OLD one, and
# nothing here distinguished the two.
DEDICATED_PIN = {"ablation", "adc", "adoptive", "checkpoint", "chemo",
                 "oncolytic", "photosensitizer_pk", "radiation", "sonodynamic"}
SHARED_PIN = {"cell", "drug_transport", "immune", "immune_spatial", "nutrient"}
ENGINE_PIN = {
    "acsl4", "alox", "biochem", "clonal", "contact", "copper", "dose_schedule",
    "grid", "ifngamma", "io", "oxygen", "params", "persister", "ph",
    "phenotype_mufa", "physics", "reaction_diffusion",
    "repair", "senescence", "slab", "spheroid", "stats", "stromal",
    "trigger_wave", "tumor_pk", "vasculature",
}



@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


def test_the_counts_are_recomputed_not_stored(d):
    """A hand-edited number is the whole failure mode this page addresses."""
    live = MD_.assemble(MD_.scan())
    for k in ("dedicated_modules", "dedicated_pub_fns", "dedicated_code_lines",
              "ferroptosis_engine_modules", "engine_pub_fns",
              "engine_code_lines"):
        assert live[k] == d[k], f"{k} drifted: live {live[k]} vs stored {d[k]}"
    # PER-MODULE fields too. Only the six scalars above were re-scanned, so
    # `reach`, `pub_types` and `pub_consts` were stored, RENDERED and never
    # checked: hand-editing them made the page publish "99 other files in the
    # workspace call it" -- in a workspace holding 55 .rs files -- and a
    # 64-line module with 77 public types, with the freshness gate green
    # because it compares the .md to a render of the same edited JSON.
    for side in ("dedicated", "shared"):
        by_live = {m["module"]: m for m in live[side]}
        for m in d[side]:
            ref = by_live.get(m["module"])
            assert ref is not None, f"{m['module']} is no longer in {side}"
            for k in ("pub_fns", "pub_types", "pub_consts", "code_lines",
                      "reach"):
                assert m[k] == ref[k], (
                    f"{side}/{m['module']}.{k} is {m[k]} in the artifact and "
                    f"{ref[k]} live")


def test_every_named_module_exists_and_every_module_is_classified(d):
    """A module missing from all three buckets would be uncounted, and one
    named that does not exist would be counted twice over."""
    for m in d["dedicated"] + d["shared"]:
        assert (CORE / f"{m['module']}.rs").exists(), m["module"]
        assert m["serves"].strip(), m["module"]
    on_disk = {p.stem for p in CORE.glob("*.rs") if p.stem != "lib"}
    named = {m["module"] for m in d["dedicated"] + d["shared"]}
    assert named <= on_disk, f"named but absent: {sorted(named - on_disk)}"
    assert d["total_modules"] == len(on_disk)
    assert (d["ferroptosis_engine_modules"] + len(named)) == len(on_disk), (
        "the buckets do not partition the crate, so some module is counted "
        "twice or not at all")


def test_the_ferroptosis_engine_is_still_the_larger_body_of_work(d):
    """The page's whole point is honesty about the gap.

    If the arms ever DO overtake the engine this fails, and the chapter's
    paragraph has to be rewritten rather than the number quietly updated --
    which is the failure that produced this page.
    """
    assert d["engine_code_lines"] > d["dedicated_code_lines"], (
        "the modality arms now exceed the ferroptosis engine in code; "
        "Chapter 6's paragraph claims the opposite and must be re-derived")
    assert d["engine_pub_fns"] > d["dedicated_pub_fns"]
    assert "still the larger body of work" in MD.read_text()


def test_the_page_refuses_to_be_read_as_progress(d):
    """Lines of code are the weakest possible evidence of depth. The page says
    so, and says what it does not measure."""
    md = MD.read_text()
    for frag in ("weakest possible evidence of depth",
                 "A count going up is not progress on its own",
                 "A module can be large and wrong",
                 "none of this appears in a figure or a claim"):
        assert frag in md, f"the page no longer says: {frag}"


def test_the_chapter_quotes_the_measurement_and_not_the_old_sentence(d):
    """The sentence this page replaced must be GONE, not merely contradicted
    somewhere else -- a retraction that leaves the original standing is how
    this repository has shipped two contradictory claims before."""
    text = MANUSCRIPT.read_text()
    assert "one function and a configuration struct each" in text, (
        "the chapter no longer records what the old sentence said, so a "
        "reader cannot tell the claim was retracted")
    # ...but only as a QUOTED retraction, never as a live claim.
    live = re.search(
        r"\*\*Most of these arms are one function and a configuration struct\.\*\*",
        text)
    assert live is None, (
        "the retracted sentence is still standing as a live claim")
    _chapter_paragraph()  # it must still exist to be quoted from


def _chapter_text() -> str:
    """Chapter 6 alone.

    `_chapter_paragraph` used to scan the whole manuscript, so the paragraph it
    "anchored" to could have been anywhere in 60,000 words.
    """
    text = MANUSCRIPT.read_text()
    start = text.index("## Chapter 6:")
    return text[start:text.index("## Chapter 7:", start)]


def _chapter_paragraph() -> str:
    """The ONE paragraph in Chapter 6 that quotes this measurement.

    Anchored deliberately. The first version of the guard below asked whether
    `str(value)` appeared anywhere in a 60,000-word manuscript, and two
    reviewers independently showed it could not fail: 265 of the integers 0-299
    already occur somewhere in `v1.md`, so rewriting the paragraph to claim the
    arms are twenty-five times LARGER than the engine left every assertion
    green. A guard whose own message says "so it can drift from the artifact
    again" and which cannot detect drift is worse than none.
    """
    text = _chapter_text()
    hits = [ln for ln in text.split("\n")
            if "modality-module-depth.md" in ln and "counts it" in ln]
    assert len(hits) == 1, (
        f"expected exactly one Chapter 6 paragraph citing the depth "
        f"measurement, found {len(hits)}")
    return hits[0]


def test_the_chapter_quotes_every_figure_it_uses_from_the_artifact(d):
    """Each number must appear in the paragraph ATTACHED TO THE RIGHT SIDE.

    Two earlier versions failed here for the same underlying reason -- they
    checked that strings were PRESENT and never what a number was attached to.
    The first asked whether `str(value)` appeared anywhere in a 60,000-word
    manuscript. The second anchored to this paragraph and added a direction
    pin, and a reviewer defeated that too: swapping only the NOUNS, so the
    paragraph read "the arms are between 11.5 and 13.8 times LARGER by line
    count ... the engine is that many times smaller", left the whole suite
    green. Every number was still present, in the same order, and
    `"engine is smaller" not in para` is satisfied by any paraphrase.

    So the assertions below are PHRASES that bind a figure to its owner. A
    paragraph that swaps the sides has to rewrite them, and rewriting them is
    the thing the guard exists to notice.
    """
    para = _chapter_paragraph()
    required = [
        f"{d['dedicated_modules']} modules a modality owns outright",
        f"{d['dedicated_pub_fns']} public functions and "
        f"{d['dedicated_code_lines']} lines of production code",
        f"{d['ferroptosis_engine_modules']} modules, {d['engine_pub_fns']} "
        f"functions and {d['engine_code_lines']:,} lines",
        f"{d['engine_modules_with_shared']}, {d['engine_pub_fns_with_shared']} "
        f"and {d['engine_code_lines_with_shared']:,}",
        f"the arms are between {d['line_ratio_narrow']} and "
        f"{d['line_ratio_wide']} times smaller by line",
        f"{d['fn_ratio_narrow']} to {d['fn_ratio_wide']} times smaller by "
        "public function",
    ]
    for phrase in required:
        assert phrase in para, (
            f"Chapter 6 does not attach the live figures to the side they "
            f"belong to; expected the phrase {phrase!r}. Paragraph: {para[:300]}")
    assert "range rather than a number" in para or "the honest figure is the interval" in para

    # A LIST OF `in` CHECKS CANNOT FAIL ON TEXT THAT WAS ADDED. Two mutants
    # kept every phrase above verbatim and were still misleading: one appended
    # "read the other way round ... the modality arms are the substantial body
    # of work here" to this same paragraph, the other put a contradicting
    # paragraph two lines below it, and both left the suite green. So the
    # paragraph is BOUNDED at both ends, and the rest of the chapter is checked
    # for a second claim about the same quantity.
    assert para.endswith("the engine could not answer."), (
        "the measured paragraph has been extended past its own conclusion; "
        f"it now ends: ...{para[-160:]!r}")
    chapter = _chapter_text()
    others = [ln for ln in chapter.split("\n")
              if ln.strip() and ln != para
              and ("times smaller" in ln or "times larger" in ln.lower()
                   or "modality-module-depth" in ln)]
    assert not others, (
        "a second passage in Chapter 6 makes a claim about the same "
        f"measurement, so the guarded paragraph is not the only one a reader "
        f"sees: {others[:1]}")


def test_the_paragraph_cannot_be_satisfied_by_a_different_number(d):
    """Anti-vacuity, run rather than argued.

    Substitutes a wrong value for each figure in turn and requires the check
    above to reject it. Without this the anchored guard could still be passing
    on a coincidence -- a small integer that happens to sit in the paragraph
    for an unrelated reason.
    """
    para = _chapter_paragraph()
    for key in ("dedicated_code_lines", "engine_code_lines",
                "engine_code_lines_with_shared"):
        wrong = d[key] + 1
        rendered = f"{wrong:,}" if d[key] >= 1000 else str(wrong)
        assert rendered not in para, (
            f"the paragraph contains {rendered} as well as the true {key}, so "
            "the guard on that figure proves nothing")


def test_the_retracted_sentence_is_quoted_and_not_standing(d):
    """A retraction that leaves the original standing has shipped here before."""
    text = MANUSCRIPT.read_text()
    assert "one function and a configuration struct each" in text, (
        "the chapter no longer records what the old sentence said, so a "
        "reader cannot tell the claim was retracted")
    live = re.search(
        r"\*\*Most of these arms are one function and a configuration struct\.\*\*",
        text)
    assert live is None, (
        "the retracted sentence is still standing as a live claim")
    para = _chapter_paragraph()
    assert "An earlier version of this paragraph said" in para, (
        "the retraction is no longer attached to the measurement that "
        "replaced it")


def test_no_engine_module_is_named_after_a_treatment_arm(d):
    """The check the pin could not make for itself.

    `test_every_module_is_classified_by_a_decision_somebody_made` catches a
    NEW module landing unclassified and is powerless against one already
    misfiled when the pin was written -- which is how `photosensitizer_pk`
    spent eight sections on the comparator's side of the gap ratio.

    So this asks a question the pin cannot: does any module the engine claims
    share a name with an arm the parity table measures? A modality's own file
    counted as engine inflates the comparator and the gap at once, in the
    direction least likely to be questioned, because it makes this project's
    self-criticism look better founded than it is.
    """
    engine = set(d["engine_module_names"])
    arm_words = {"photo", "sono", "radiat", "chemo", "immun", "onco", "ablat",
                 "adc", "adoptive", "checkpoint", "conjugate", "cart"}
    suspicious = sorted(m for m in engine
                        if any(w in m for w in arm_words)
                        and m not in SHARED_PIN)
    assert not suspicious, (
        f"{suspicious} are counted as ferroptosis engine while named after a "
        "treatment modality. If one of them really is shared machinery, put "
        "it in SHARED (credited to no arm) rather than leaving it on the "
        "comparator's side of a ratio that measures distance from the "
        "comparator.")


def test_every_module_is_classified_by_a_decision_somebody_made(d):
    """A new modality file must NOT silently become ferroptosis engine.

    `engine` is computed as the complement of the two hand-written buckets, so
    the partition assertion is true by construction and cannot fail. A
    reviewer added `bispecific.rs` -- a modality's own file -- and it was
    counted as ferroptosis engine, making the engine look LARGER because a
    modality arm landed. That is the exact self-flattery this page exists to
    prevent, running backwards.

    Pinning all three sets is the fix: any new module fails this until someone
    classifies it. The literal below is the decision record.
    """
    engine = set(d["engine_module_names"])
    dedicated = {m["module"] for m in d["dedicated"]}
    shared = {m["module"] for m in d["shared"]}
    on_disk = {p.stem for p in CORE.glob("*.rs") if p.stem != "lib"}

    assert dedicated == DEDICATED_PIN
    assert shared == SHARED_PIN
    # THE ENGINE SET MUST BE PINNED TOO, and the first version of this test
    # said it was while reading `engine` out of the JSON the generator wrote --
    # so `engine | dedicated | shared == on_disk` was still true by
    # construction. A reviewer added `bispecific.rs`, a modality's own file, and
    # it was counted AS ferroptosis engine with every test green: the engine
    # grew because a modality arm landed, which is the page's own self-flattery
    # running backwards. Pinning the complement is what actually forces a
    # decision.
    assert engine == ENGINE_PIN, (
        "the ferroptosis-engine bucket changed. If a module was added, put it "
        "in DEDICATED, SHARED or ENGINE_PIN deliberately -- do NOT let the "
        "complement absorb it, because an unclassified modality file inflates "
        f"the engine it is being compared against. Added: "
        f"{sorted(engine - ENGINE_PIN)}; removed: {sorted(ENGINE_PIN - engine)}")
    assert engine | dedicated | shared == on_disk, (
        "a module on disk is in no bucket: "
        f"{sorted(on_disk - (engine | dedicated | shared))}")
    assert not (engine & (dedicated | shared)), "a module is in two buckets"
    assert dedicated and shared and engine, (
        "an empty bucket satisfies every structural check in this file")


def test_the_shared_reach_is_measured_not_asserted(d):
    """The page said "four arms reach through `immune.rs`" beside a table row
    naming six, and neither number came from the code."""
    md = MD.read_text()
    # The page keeps the wrong figure as a QUOTED retraction, exactly as the
    # chapter keeps the retracted sentence -- so the guard must forbid the
    # LIVE form and require the quoted one, not ban the string.
    assert "four arms reach through" not in md, (
        "the hand-written reach count is standing as a live claim again")
    assert 'said "four arms" beside a table row naming six' in md, (
        "the page no longer records that the figure was hand-written")
    deepest = max(d["shared"], key=lambda m: m["reach"])
    assert f"**{deepest['reach']} other files in the workspace call it**" in md
    for m in d["shared"]:
        assert m["reach"] >= 1, (
            f"{m['module']} is in the SHARED bucket and nothing calls it, so "
            "it is not machinery several arms reach through")
