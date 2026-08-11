"""The manuscript must not call self-consistency checks calibration.

`simulations/calibration/targets.yaml` holds eight targets and every one is
`target_type: self-consistency` — they regression-guard the model against its own
previously reported behaviour. None is a measurement the model was fitted to.

#589 established that and relabelled it across CLAUDE.md, README.md,
CALIBRATION_STATUS.md and targets.yaml itself, recording in its commit message
that targets.yaml "now correctly holds 0 calibration + 8 self-consistency". It
never touched the manuscript, which went on saying "Five calibration constraints
anchor the model to reality" — the strongest form of the overclaim this repo's
honesty infrastructure exists to prevent, sitting in the document that travels
furthest from the evidence.

MODEL_CARD.md compounded it by disclosing only "the three 3D" targets as
self-consistency, which affirmatively implies the other five are independent.

These guards are pinned to targets.yaml, not to the wording, so they follow the
data if a target is ever genuinely upgraded. Verified by mutation: flipping one
target to `target_type: calibration` makes the two "do not call it calibration"
guards STAND DOWN, while two others fire deliberately -- the count no longer
matches, and `test_the_comment_is_not_mistaken_for_a_target` refuses to pass in
silence. That is the intent. A real upgrade should force someone to re-read the
chapter, not quietly re-permit the old framing.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TARGETS = REPO_ROOT / "simulations" / "calibration" / "targets.yaml"
MS = REPO_ROOT / "article" / "drafts" / "v1.md"
CARD = REPO_ROOT / "MODEL_CARD.md"


def _target_types():
    """(n_targets, n_self_consistency, n_calibration), comments excluded.

    targets.yaml explains in a COMMENT how to upgrade a target to
    `target_type: calibration`, so a naive grep for that string finds a hit and
    reads as one real calibration target. Comment lines are stripped first.
    """
    body = "\n".join(l for l in TARGETS.read_text().split("\n")
                     if not l.strip().startswith("#"))
    return (len(re.findall(r"^\s*-?\s*id:\s*\S+", body, re.M)),
            body.count("target_type: self-consistency"),
            body.count("target_type: calibration"))


def test_the_comment_is_not_mistaken_for_a_target():
    """The parser this file relies on must exclude the how-to comment."""
    raw = TARGETS.read_text()
    assert "target_type: calibration" in raw, (
        "the upgrade instructions are gone; this guard's premise no longer holds")
    _, _, n_cal = _target_types()
    assert n_cal == 0, (
        f"{n_cal} real calibration target(s) now exist -- if that is genuine, the "
        "manuscript's framing should be revisited rather than this test relaxed")


def test_the_manuscript_does_not_call_them_calibration_constraints():
    n, n_sc, n_cal = _target_types()
    txt = MS.read_text()
    if n_cal:
        return                       # genuinely calibrated now; claim defensible
    assert "calibration constraints anchor the model to reality" not in txt, (
        "the manuscript calls self-consistency checks calibration constraints "
        "that anchor the model to reality; targets.yaml says otherwise")
    assert "self-consistency" in txt, (
        "the manuscript never uses the word for what these targets are")


def test_the_manuscript_states_the_count_the_artifact_holds():
    n, n_sc, _ = _target_types()
    flat = " ".join(MS.read_text().split())
    assert f"{n_sc} of {n} carry" in flat, (
        f"the manuscript does not state the measured split ({n_sc} of {n} "
        "self-consistency); a hand-written count here is what drifted before")


def test_the_model_card_does_not_imply_only_the_3d_targets_are_self_consistency():
    n, _, n_cal = _target_types()
    if n_cal:
        return
    txt = CARD.read_text()
    assert 'the three 3D "validation"\ntargets are **self-consistency**' not in txt, (
        "MODEL_CARD names only the 3D targets as self-consistency, implying the "
        "other five are independent measurements")
    assert "ALL EIGHT" in txt or f"all {n}" in txt.lower(), (
        "MODEL_CARD does not say every target is a self-consistency check")


def test_the_figure_7_caption_describes_the_figure_that_ships():
    """The caption claimed error bars and a panel (b) the image does not have.

    Scoped to `v1.md` on purpose. Figure captions are maintained TWICE and
    independently: the markdown carries an inline `[FIGURE n: ...]` description,
    and `scripts/generate_latex.py` substitutes its own caption from a dict when
    it builds the .tex. The LaTeX caption for this figure is one clause and
    claims nothing, so the PDF was never wrong -- only the markdown was. That
    divergence is itself a drift surface worth knowing about, but it is not what
    this guard is for.
    """
    txt = MS.read_text()
    i = txt.index("[FIGURE 7:")
    caption = txt[i:txt.index("]", i) + 1]
    assert "no error bars" in caption, (
        "the Figure 7 caption no longer says the figure carries no error bars")
    assert "tornado" not in caption.lower() or "previously described" in caption, (
        "the caption describes a sensitivity tornado panel; the shipped figure "
        "is a single panel of bars")
    fig = REPO_ROOT / "article" / "figures" / "fig7_monte_carlo_simulation.png"
    assert fig.exists(), "the figure the caption describes is missing"
