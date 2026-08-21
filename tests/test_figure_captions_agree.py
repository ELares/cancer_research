"""Every figure has TWO captions, and nothing tied them together.

`article/drafts/v1.md` carries an inline `[FIGURE N: ...]` caption, which is
what a reader of the markdown sees. `scripts/generate_latex.py` carries its own
captions dict, which is what a reader of the PDF sees. They are written and
edited separately, so they can describe different images -- and three of them
did:

* FIGURE 3's markdown caption described a 19x19 co-occurrence heatmap while
  the asset was a ferroptosis-engagement bar chart.
* FIGURE 4's described a network graph of mechanism convergence while the asset
  was a molecular-pathway chart.
* FIGURE 5's described a mechanism-by-site matrix while the asset was
  "literature disconnect between communities".

The first two predate this campaign: the figures were reordered at some point
and the markdown captions stayed with the old numbers. The third I introduced
by rewriting a section's caption without repointing its figure -- the same
caption/asset defect this project already fixed once at scale, reappearing
because nothing checked it.

THE CHECK IS DELIBERATELY WEAK AND THAT IS THE POINT. Two prose captions for
one image cannot be compared for meaning, so this asks only whether they share
a single distinctive content word. Zero shared words means the two descriptions
have no subject in common, which is unambiguous drift rather than paraphrase.
Measured across all 24 figures at the time of writing, the three broken ones
scored exactly 0 and the lowest surviving score was 1 -- the bar is a floor,
not a threshold tuned to fit.
"""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANUSCRIPT = REPO / "article/drafts/v1.md"
LATEX = REPO / "scripts/generate_latex.py"

# Words too generic to evidence a shared subject. Kept short on purpose: a long
# stoplist would start deciding which agreements count.
STOP = set(
    "the a an and or of in to by for with at on from is are as its it this that "
    "each per not shown showing show over under across between into within "
    "where which what when than left right top bottom axis panel figure".split()
)


def _content_words(text: str) -> set:
    text = re.sub(r"\\\\?[a-z]+\{?|\$|[^a-z0-9 -]", " ", text.lower())
    return {w for w in text.split() if len(w) > 3 and w not in STOP}


def _markdown_captions() -> dict:
    md = MANUSCRIPT.read_text(encoding="utf-8")
    return {int(m.group(1)): m.group(2)
            for m in re.finditer(r"\[FIGURE (\d+): ([^\]]*)", md)}


def _latex_captions() -> dict:
    src = LATEX.read_text(encoding="utf-8")
    return {int(k): (f, c) for k, f, c in
            re.findall(r"'(\d+)': \('([^']+)', '((?:[^'\\]|\\.)*)'", src)}


def test_the_two_caption_sources_cover_the_same_figures():
    md, tex = _markdown_captions(), _latex_captions()
    only_md = sorted(set(md) - set(tex))
    only_tex = sorted(set(tex) - set(md))
    assert not only_md, (
        f"figures {only_md} are captioned in the manuscript but have no LaTeX "
        "entry, so the PDF will not render them")
    assert not only_tex, (
        f"figures {only_tex} have a LaTeX entry but no manuscript caption, so "
        "the markdown reader never learns they exist")


def test_no_figure_has_two_captions_describing_different_images():
    """The check that would have caught all three.

    Reports every failing figure at once rather than the first, because these
    arrive in groups -- a reordering breaks a run of them.
    """
    md, tex = _markdown_captions(), _latex_captions()
    bad = []
    for n in sorted(md):
        if n not in tex:
            continue
        asset, latex_caption = tex[n]
        shared = _content_words(md[n]) & _content_words(latex_caption)
        if not shared:
            bad.append(
                f"FIGURE {n} ({asset}): markdown says "
                f"{md[n][:60]!r} while the LaTeX caption says "
                f"{latex_caption[:60]!r} -- no content word in common")
    assert not bad, (
        "the manuscript and the LaTeX render describe different images:\n  "
        + "\n  ".join(bad)
        + "\nOne of the two captions is describing a figure that is no longer "
          "there. Fix the caption, or repoint the figure.")


def test_the_bar_is_a_floor_not_a_tuned_threshold():
    """A guard calibrated until it passes is not a guard.

    If most figures sat at 1 or 2 shared words, requiring 1 would be a
    threshold fitted to the data rather than a principled floor. The margin is
    checked so that a future rewrite which drifts every caption toward
    generic language fails here rather than silently weakening the check.
    """
    md, tex = _markdown_captions(), _latex_captions()
    scores = sorted(len(_content_words(md[n]) & _content_words(tex[n][1]))
                    for n in md if n in tex)
    assert scores, "no figures to score"
    median = scores[len(scores) // 2]
    assert median >= 5, (
        f"the median caption pair shares only {median} content words, so "
        "requiring 1 is close to the typical case and the check has become a "
        "tuned threshold rather than a floor. The captions have drifted toward "
        "generic language.")


def test_every_captioned_figure_has_an_asset_on_disk():
    """A caption pointing at a file that is not there fails only at build
    time, and there is no LaTeX toolchain here to discover it."""
    tex = _latex_captions()
    missing = sorted(f"{n} ({asset})" for n, (asset, _) in tex.items()
                     if not (REPO / "article/figures" / f"{asset}.pdf").exists())
    assert not missing, f"figures with no asset: {missing}"
