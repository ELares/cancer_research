"""Manuscript captions that quote a statistic must match what the figure draws.

`fig1_ferroptosis_comparison` computes a chi-squared test and PRINTS it; the
caption in `article/drafts/v1.tex` carried the printed value by hand. The
corpus was then halved -- 10,413 files to 4,830 -- and retagged, and both the
figure and the caption went on stating numbers from the superseded corpus for
five months. The figure's own denominators (`n=617`, `n=563`) match the parent
of the commit that introduced it exactly, so it was stale the day it landed.

The figure is now regenerated and the caption reads chi^2 = 38.8. This derives
that number from the corpus rather than trusting either, which is the only
arrangement in which the two cannot drift apart again: a hand-written number
beside a computed one is this repository's most repeated defect.
"""
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TEX = REPO / "article/drafts/v1.tex"
GEN = REPO / "scripts/generate_figures.py"


def _corpus_derived_stems():
    """Figure stems FIGURES.yaml declares corpus-derived from this generator."""
    import yaml

    figs = yaml.safe_load((REPO / "FIGURES.yaml").read_text())["figures"]
    return sorted(
        f["filename"] for f in figs
        if f.get("type") == "corpus-derived"
        and str(f.get("generator", {}).get("script", "")).endswith(
            "generate_figures.py"))


def _generator():
    sys.path.insert(0, str(REPO / "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import generate_figures
    return generate_figures


def _caption_statistic():
    """The chi-squared value and p the caption states, as floats."""
    text = TEX.read_text()
    m = re.search(
        r"\\caption\{Ferroptosis engagement \(\$\\chi\^2=([\d.]+)\$, "
        r"\$p=([\d.]+)\\times10\^\{(-?\d+)\}\$",
        text)
    assert m, (
        "the fig1 caption no longer states a chi-squared value in the form "
        "this test reads; if the statistic was dropped, drop this test with it")
    return float(m.group(1)), float(m.group(2)) * 10 ** int(m.group(3))


def _figure_annotation():
    """The chi-squared and p the committed fig1 actually draws."""
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            return None
    doc = pymupdf.open(REPO / "article/figures/fig1_ferroptosis_comparison.pdf")
    try:
        text = " ".join(" ".join(page.get_text().split()) for page in doc)
    finally:
        doc.close()
    m = re.search(r"\u03c7\u00b2\s*=\s*([\d.]+),\s*p\s*=\s*([\d.]+e?-?\d*)", text)
    if not m:
        return None
    return float(m.group(1)), m.group(2)


def test_the_fig1_caption_matches_the_corpus():
    g = _generator()
    from scipy import stats as sp

    data = g.classify_ferroptosis(g.load_corpus())
    modalities = ["SDT", "IRE", "HIFU", "TTFields", "Frequency"]
    sdt_ferro = data["SDT"]["ferroptosis"]
    sdt_total = data["SDT"]["total"]
    other_ferro = sum(data[m]["ferroptosis"] for m in modalities if m != "SDT")
    other_total = sum(data[m]["total"] for m in modalities if m != "SDT")
    chi2, p, _, _ = sp.chi2_contingency(
        [[sdt_ferro, sdt_total - sdt_ferro],
         [other_ferro, other_total - other_ferro]])

    stated_chi2, stated_p = _caption_statistic()
    assert round(chi2, 1) == stated_chi2, (
        f"the caption states chi^2 = {stated_chi2} and the corpus gives "
        f"{chi2:.1f}. Regenerate fig1 and update the caption together -- the "
        "previous mismatch survived a corpus halving because nothing checked")
    # TWO-SIDED, to two significant figures. The first version asserted
    # `abs(p - stated_p) < stated_p`, which holds for EVERY stated_p >= p --
    # a caption claiming `p = 1.0`, i.e. no significance at all, passed it.
    # Only an understated p could fail, which is the wrong half: a caption
    # that overstates its own p is the one that misleads a reader.
    assert f"{p:.1e}" == f"{stated_p:.1e}", (
        f"the caption states p = {stated_p:.1e} and the corpus gives "
        f"{p:.1e}. These must agree to the digits both display")
    # AND THE FIGURE MUST DISPLAY THE SAME NUMBERS. The caption sits beside a
    # figure that draws its own annotation, and a reviewer found the two
    # disagreeing in the digit both print (caption 4.8, figure 4.7) in the very
    # commit written to stop a caption drifting from its figure. Read out of
    # the committed PDF, so this compares the artifact a reader sees.
    drawn = _figure_annotation()
    assert drawn is not None, "fig1 no longer annotates its chi-squared test"
    drawn_chi2, drawn_p = drawn
    assert drawn_chi2 == stated_chi2 and drawn_p == f"{stated_p:.1e}", (
        f"the figure draws chi^2 = {drawn_chi2}, p = {drawn_p} and the caption "
        f"states chi^2 = {stated_chi2}, p = {stated_p:.1e}")


def test_the_caption_statistic_is_the_one_the_figure_draws():
    """The test above recomputes the contingency table. If the GENERATOR ever
    builds a different one, the caption would agree with this test and disagree
    with the figure -- two artifacts consistent with each other and both wrong,
    which is the shape of defect this file exists to catch."""
    src = GEN.read_text()
    body = src[src.index("def fig1_ferroptosis_comparison"):]
    body = body[:body.index("\ndef ")]
    flat = " ".join(body.split())
    # THE OTHERS ARM TOO. Pinning only the SDT arm left `other_ferro`,
    # `other_total` and the `modalities` list free -- and those decide the
    # table. Measured: rewriting the others arm to drop one modality makes the
    # figure draw chi^2 = 36.3 while the caption and this test both still say
    # 38.8, which is verbatim the "two artifacts agreeing while both diverge
    # from the figure" failure this test is here to exclude.
    for fragment in (
            'modalities = ["SDT", "IRE", "HIFU", "TTFields", "Frequency"]',
            "sdt_ferro = data[\"SDT\"][\"ferroptosis\"]",
            "sdt_total = data[\"SDT\"][\"total\"]",
            "other_ferro = sum(data[m][\"ferroptosis\"] for m in modalities if m != \"SDT\")",
            "other_total = sum(data[m][\"total\"] for m in modalities if m != \"SDT\")",
            "contingency = [[sdt_ferro, sdt_total - sdt_ferro], [other_ferro, other_total - other_ferro]]",
            "chi2, p_value, _, _ = stats.chi2_contingency(contingency)"):
        assert " ".join(fragment.split()) in flat, (
            "fig1 no longer computes its chi-squared the way the caption test "
            f"reproduces it; missing: {fragment}")


def test_the_regenerated_corpus_figures_are_what_the_generator_draws(tmp_path):
    """The corpus figures with TRACKED inputs, regenerated and compared.

    The set is DISCOVERED from `FIGURES.yaml` -- the entries whose `type` is
    `corpus-derived` and whose generator is this script -- and asserted to be
    exactly what this test compares. Hard-coding it meant `tracked[:0]` passed,
    and a new corpus-derived figure would have been added to the manuscript and
    silently skipped here.

    The generator writes 20 figure stems, 13 of which are committed. Eight of
    those thirteen read `simulations/output/`, which is gitignored -- nothing
    under it is tracked but `.gitkeep` -- so regenerating them would commit
    plots drawn from data nobody else has and CI cannot reproduce:
    fig8_simulation_by_treatment, fig10_invivo_comparison, fig11_mufa_sweep,
    fig17_damp_heatmap, fig24_hypoxia_killcurve, fig25_bliss_synergy,
    fig26_vulnerability_window, fig27_resistance_asymmetry. An earlier version
    of this docstring said "three", naming only the last three of the eight.
    Seven more stems are drawn on every run and have never been committed.
    Both sets are reported on issue #788 rather than swept in here.

    PNG CONTENT IS NOT CHECKED, only the PDFs -- the same hole the census gate
    documents, for the same reason (PNG bytes are not portable across
    platforms, and CI compares against figures authored on macOS). Measured:
    restoring a superseded PNG beside a fresh PDF passes. The PNGs are
    committed and hashed into MANIFEST.sha256, so a swap is caught there as a
    changed file, but not here as a stale figure.
    """
    import os
    import subprocess

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_figure_freshness import _drawing

    out = tmp_path / "figs"
    out.mkdir()
    res = subprocess.run([sys.executable, str(GEN)], cwd=REPO,
                         capture_output=True, text=True,
                         env={**os.environ, "MPLBACKEND": "Agg",
                              "FERRO_FIG_DIR": str(out)})
    assert res.returncode == 0, res.stderr[-800:]
    tracked = _corpus_derived_stems()
    assert tracked == sorted([
        "fig1_ferroptosis_comparison", "fig4_molecular_overlap",
        "fig6_sdt_chain_evidence", "fig12_pathway_targets",
        "fig13_gold_set_eval"]), (
        f"FIGURES.yaml now declares {tracked} as corpus-derived from this "
        "generator. Add the new figure to this comparison deliberately, or "
        "correct FIGURES.yaml -- do not let the set drift silently")
    compared = 0
    for stem in tracked:
        produced = out / f"{stem}.pdf"
        assert produced.exists(), f"the generator no longer draws {stem}"
        committed = REPO / "article/figures" / f"{stem}.pdf"
        a, b = _drawing(produced), _drawing(committed)
        if a is None:
            pytest.skip("no PDF reader available")
        assert a == b, (
            f"article/figures/{stem}.pdf is not what "
            "scripts/generate_figures.py draws. Re-run it.")
        compared += 1
    # A COUNT, because slicing the loop to `tracked[:0]` left this green. The
    # census gate learned the same lesson through six review rounds: a check
    # that never runs reports exactly what a check that passes reports.
    assert compared == len(tracked), (
        f"compared {compared} figures, not {len(tracked)}")


# Strings that MUST match, and strings that must NOT, for the acronym patterns
# these figures count with. Both halves are needed: unbounded, `STING` counted
# "suggesting" and 91% of fig4's largest column was artifact; bounded on the
# right, `DAMP` stopped counting "DAMPs", which is how almost everyone writes
# it, and dropped a real SDT article. A fix that trades one error for the other
# is not an improvement, and nothing in this repo pinned either direction.
_PATTERN_CASES = {
    "ICD/DAMPs": (
        ["damage-associated molecular patterns (damps)", "release of damps",
         "a damp signal", "calreticulin exposure", "hmgb1 release",
         "immunogenic cell death"],
        ["dampen the response", "dampened", "hmgb1a", "damping"],
    ),
    "GSH/GPX4": (
        ["gsh depletion", "gsh/gpx4 axis", "glutathione peroxidase",
         "gpx4 inactivation", "slc7a11 expression", "intracellular gsh)"],
        ["changsha university", "gshan", "gpx4l"],
    ),
    "STING/cGAS": (
        ["sting pathway activation", "cgas-sting signalling", "the sting agonist"],
        ["suggesting that", "existing therapies", "boosting immunity",
         "consisting of", "interestingly", "testing the hypothesis"],
    ),
    "ROS": (
        ["ros generation", "ros-generating nanoparticles", "reactive oxygen species"],
        ["necrosis", "prostate cancer", "across the tumour", "sclerosis"],
    ),
    "Ferroptosis": (["ferroptosis", "ferroptosis inducer"], ["apoptosis", "necroptosis"]),
    "Apoptosis": (["apoptosis", "caspase-3 activation"], ["ferroptosis", "autophagy"]),
    "Autophagy": (["autophagy", "autophagic flux"], ["apoptosis", "phagocytosis"]),
    "ER Stress": (
        ["er stress response", "er stressor thapsigargin", "er stresses",
         "endoplasmic reticulum stress"],
        ["other stresses", "under stress conditions", "her stress levels"],
    ),
}


# fig6's chain, which no test touched until a reviewer reinstated the exact
# defect two commits had just removed -- a bare `glutathione` in the depletion
# bar, overstating the chain's second link by 79% -- and the suite stayed
# green. Every regex defect found in five review rounds lived in one of these
# two dicts.
_CHAIN_CASES = {
    "ROS\ngeneration": (
        ["ros generation", "generating ros", "generate high levels of ros",
         "induce ros locally", "elicits ros-mediated apoptosis",
         "reactive oxygen species", "ros production",
         "boost ultrasound (us)-initiated ros production",
         "reactive singlet oxygen species (ros)",
         "generate reactive oxide species (ros)"],
        ["necrosis", "prostate", "across", "ros accumulation is reduced"],
    ),
    "GSH\ndepletion": (
        ["gsh depletion", "depletes intra-tumoral glutathione",
         "depletion of overexpressed glutathione", "consume the reduced glutathione",
         "scavenges gsh via the mno2", "glutathione consumption"],
        ["glutathione peroxidase 4 expression", "intracellular glutathione levels",
         "glutathione (gsh) levels rapidly scavenge the sdt-induced ros"],
    ),
    "GPX4\ninactivation": (
        ["gpx4 inactivation", "glutathione peroxidase 4"],
        ["gpx4l", "agpx4"],
    ),
    "Lipid\nperoxidation": (["lipid peroxidation", "lipid peroxide"], ["peroxidase"]),
    # "ferroptotic" is NOT required: measured, no corpus article says it
    # without also saying "ferroptosis" (79 either way), so demanding it would
    # pin a pattern change nothing observes. The case says what is true.
    "Ferroptosis": (["ferroptosis", "ferroptosis-inducing"], ["apoptosis"]),
    "DAMP\nrelease": (
        ["damage-associated molecular patterns", "damps", "hmgb1", "hmgb-1",
         "calreticulin"],
        ["dampen", "damage associated with replication stress"],
    ),
    "STING\nactivation": (
        ["sting pathway", "cgas-sting"],
        ["suggesting", "existing", "boosting"],
    ),
    "ICD": (["immunogenic cell death"], ["cell death"]),
}


def _pattern_dict(name):
    """A pattern dict as the generator defines it, executed from its source."""
    src = GEN.read_text()
    marker = f"    {name} = {{"
    assert src.count(marker) == 1, (
        f"expected exactly one `{name}` dict in the generator, found "
        f"{src.count(marker)}")
    body = src[src.index(marker):]
    body = body[:body.index("\n    }") + 6]
    ns = {}
    exec("import re\n" + body.strip(), ns)          # noqa: S102 - our own source
    return ns[name]


def test_every_pattern_has_a_case():
    """A case table that does not cover the dict is a check that never runs.

    `_PATTERN_CASES` covered five of `pathways`' eight labels and none of
    `chain_steps`, so three patterns in one figure and all eight in the other
    could be rewritten with the suite green -- which a reviewer demonstrated
    by reinstating a defect two commits had removed.
    """
    for name, cases in (("pathways", _PATTERN_CASES), ("chain_steps", _CHAIN_CASES)):
        declared = set(_pattern_dict(name))
        assert declared == set(cases), (
            f"{name} declares {sorted(declared - set(cases))} with no case and "
            f"{sorted(set(cases) - declared)} cased but absent")
        # AND EACH ENTRY MUST SAY SOMETHING. Key coverage alone let an entry
        # be emptied to `([], [])` with the suite green, which is a case table
        # that covers the dict and tests none of it -- the vacuity the
        # generator comment claims this table forecloses.
        for label in sorted(cases):
            should, should_not = cases[label]
            assert should and should_not, (
                f"{name}[{label!r}] has {len(should)} must-match and "
                f"{len(should_not)} must-not-match cases; both directions are "
                "needed or the pattern is pinned in one direction only")


@pytest.mark.parametrize("which,label", (
    [("pathways", k) for k in sorted(_PATTERN_CASES)]
    + [("chain_steps", k) for k in sorted(_CHAIN_CASES)]))
def test_the_patterns_match_what_they_claim_to(which, label):
    """The patterns are read from the generator, not restated here."""
    import re as _re

    pattern = _pattern_dict(which)[label]
    should, should_not = (_PATTERN_CASES if which == "pathways"
                          else _CHAIN_CASES)[label]
    for text in should:
        assert _re.search(pattern, text, _re.IGNORECASE), (
            f"/{label}/ no longer matches {text!r}; a right-hand word boundary "
            "on an acronym drops its plural, which is how DAMPs was lost")
    for text in should_not:
        assert not _re.search(pattern, text, _re.IGNORECASE), (
            f"/{label}/ matches {text!r}, which is a substring collision -- "
            "the defect that made 91% of one column of fig4 an artifact")


# The pre-cut corpus, recoverable from git. Extracting all 10,413 files takes
# about a second, which is what makes gating the decomposition possible at all.
PRE_CUT_COMMIT = "f65342df"


def _corpus_at(commit, tmp):
    """The by-pmid corpus as it stood at `commit`, parsed like the generator."""
    import subprocess
    import tarfile
    import io
    import yaml

    out = subprocess.run(
        ["git", "archive", commit, "corpus/by-pmid"],
        cwd=REPO, capture_output=True)
    if out.returncode != 0:
        shallow = subprocess.run(
            ["git", "rev-parse", "--is-shallow-repository"],
            cwd=REPO, capture_output=True, text=True).stdout.strip()
        raise AssertionError(
            f"cannot read the corpus at {commit}: "
            f"{out.stderr[-300:].decode(errors='replace')}"
            + ("\n\nThis clone is SHALLOW. The retraction in Section 8.2 is a "
               "claim about this repository's own history, so checking it needs "
               "that history: use `fetch-depth: 0` (which .github/workflows/"
               "python-test.yml sets for exactly this reason) or "
               "`git fetch --unshallow`." if shallow == "true" else ""))
    with tarfile.open(fileobj=io.BytesIO(out.stdout)) as tar:
        tar.extractall(tmp, filter="data")
    arts = []
    for f in sorted((Path(tmp) / "corpus/by-pmid").glob("*.md")):
        m = re.match(r"^---\n(.*?\n)---\n\n?(.*)", f.read_text(encoding="utf-8"),
                     re.DOTALL)
        if not m:
            continue
        fm = yaml.safe_load(m.group(1)) or {}
        fm["_text"] = (fm.get("title", "") + " " + m.group(2).lower()[:4000]).lower()
        fm["_pmid"] = f.stem
        arts.append(fm)
    return arts


def _sdt(arts):
    return [a for a in arts if "sonodynamic" in (a.get("mechanisms") or [])]


def test_the_chapter_8_decomposition_is_derived_not_asserted(tmp_path):
    """Section 8.2 states fifteen numbers about a retraction. Derive all of them.

    This test exists because the sentence it checks is the THIRD attempt. The
    first hand-wrote a count over a corpus that had been halved underneath it
    and stood for five months. The second explained the change with a mechanism
    that could not have produced it -- refuted from this repository's own git
    history by a reviewer in an afternoon. Both were prose beside correct
    figures, which is exactly the position the fig1 caption was in.

    So the same treatment the caption got: every number recomputed here, from
    the current corpus and from the pre-cut one, and the sentence read back.
    """
    import matplotlib
    import subprocess

    matplotlib.use("Agg")
    g = _generator()

    old = _corpus_at(PRE_CUT_COMMIT, tmp_path)
    new = g.load_corpus()
    so, sn = _sdt(old), _sdt(new)
    # READ FROM THE GENERATOR, not restated here. Hard-coded, this test kept
    # passing when the generator's instrument was narrowed -- so the sentence
    # would have gone on describing an axis the code no longer counted, which
    # is the drift it exists to prevent, one level up.
    axis, icd = _axis_patterns()

    def count(pool, pat):
        return sum(1 for a in pool if re.search(pat, a["_text"], re.I))

    def dual(pool):
        return sum(1 for a in pool
                   if "ferroptosis" in a["_text"] and re.search(icd, a["_text"], re.I))

    facts = {
        "old corpus size": len(old),
        "new corpus size": len(new),
        "old SDT": len(so),
        "new SDT": len(sn),
        "retired instrument, old corpus": count(so, r"glutathione"),
        "widened instrument, old corpus": count(so, axis),
        "widened instrument, new corpus": count(sn, axis),
        "dual old": dual(so),
        "dual new": dual(sn),
    }
    # The generator's loader does not synthesise `_pmid`; the frontmatter
    # carries `pmid`, sometimes quoted. One accessor for both corpora.
    def pmid(a):
        return str(a.get("pmid") or a.get("_pmid") or "").strip()

    assert all(pmid(a) for a in so) and all(pmid(a) for a in sn), (
        "an article carries no pmid, so the deleted/retagged split below "
        "would silently compare empty strings")
    lost = {pmid(a) for a in so} - {pmid(a) for a in sn}
    still_held = {pmid(a) for a in new}
    facts["deleted"] = len(lost - still_held)
    facts["retagged"] = len(lost & still_held)
    by_pmid = {pmid(a): a for a in so}
    facts["retagged and GSH-positive"] = sum(
        1 for p in (lost & still_held) if re.search(axis, by_pmid[p]["_text"], re.I))

    # The identity that makes "deleted" and "retagged" mean what they say.
    assert facts["deleted"] + facts["retagged"] + facts["new SDT"] == facts["old SDT"], (
        f"the SDT set does not decompose: {facts}")

    sentence = _chapter_8_sentence()
    pcts = {
        "retired": 100 * facts["retired instrument, old corpus"] / facts["old SDT"],
        "widened old": 100 * facts["widened instrument, old corpus"] / facts["old SDT"],
        "widened new": 100 * facts["widened instrument, new corpus"] / facts["new SDT"],
    }
    required = [
        (f"{facts['widened instrument, new corpus']} GSH/GPX4 articles", "current count"),
        (f"({pcts['widened new']:.1f}% of SDT-tagged articles)", "current share"),
        (f"{facts['dual new']} dual-pathway", "current dual-pathway"),
        (f"read {facts['retired instrument, old corpus']} "
         f"({pcts['retired']:.1f}%) and {facts['dual old']}", "the retired trio"),
        (f"{facts['retired instrument, old corpus']}/{facts['old SDT']} "
         f"= {pcts['retired']:.1f}%", "the retired division"),
        (f"{facts['widened instrument, old corpus']} of those same {facts['old SDT']}, "
         f"or {pcts['widened old']:.1f}%", "the instrument step"),
        (f"{facts['old corpus size']:,} articles to {facts['new corpus size']:,}", "the cut"),
        (f"{facts['widened instrument, new corpus']} of {facts['new SDT']}, "
         f"or {pcts['widened new']:.1f}%", "the corpus step"),
        (f"{facts['deleted']} of the {facts['old SDT']} were deleted", "the deletions"),
        (f"only {facts['retagged']} were retagged", "the retags"),
        (f"{facts['dual old']} to {facts['dual new']}", "dual-pathway movement"),
        (PRE_CUT_COMMIT, "the recovery commit"),
    ]
    for fragment, what in required:
        assert fragment in sentence, (
            f"Section 8.2 does not state {what} as measured -- expected "
            f"{fragment!r}. Recompute the sentence; do not adjust this test.")

    # THE CLAIMS, not only the numerals. Four mutations of this paragraph left
    # the first version of this guard green: "rose" to "fell", "glutathione"
    # to "GPX4", "March 2026" to "January 2019", and "neither of them
    # GSH-positive" to "both of them". A sentence made of fifteen correct
    # numbers can still say the opposite of what they show.
    positive = facts["retagged and GSH-positive"]
    assert positive == 0, (
        "a retagged article is GSH-positive, so the paragraph's 'neither of "
        "them GSH-positive' is now false")
    assert "neither of them GSH-positive" in sentence, (
        f"{positive} of the retagged articles are GSH-positive, and the "
        "paragraph no longer says so")
    # EVERY DIRECTION WORD, derived. Pinning only "rose" left "the count fell"
    # and "Dual-pathway articles fall" free, and flipping both left the whole
    # suite green with the manuscript stating both movements backwards.
    for moved, rose, fell, what in (
            (pcts["widened new"] - pcts["retired"],
             "the share rose by about", "the share fell by about", "the share"),
            (facts["widened instrument, new corpus"]
             - facts["retired instrument, old corpus"],
             "while the count rose because", "while the count fell because",
             "the count"),
            (facts["dual new"] - facts["dual old"],
             "Dual-pathway articles rise", "Dual-pathway articles fall",
             "the dual-pathway count")):
        want, wrong = (rose, fell) if moved > 0 else (fell, rose)
        assert moved != 0, f"{what} did not move; the paragraph claims it did"
        assert want in sentence, (
            f"{what} moved by {moved:+.4g} and the paragraph does not say so "
            f"-- expected {want!r}")
        assert wrong not in sentence, (
            f"{what} moved by {moved:+.4g} and the paragraph states the "
            f"opposite direction: {wrong!r}")
    assert "counting the word *glutathione* alone" in sentence, (
        "the paragraph no longer names the retired instrument, or names a "
        "different one -- it is `glutathione` alone that reproduces 72")
    # The reconstruction is HEDGED, because it is one. No text window over the
    # pre-cut corpus reproduces the source page's 39 ferroptosis and 21 ICD
    # counts alongside its 72, so the instrument is inferred, not recorded.
    # THE SURVEY/PRIMARY SPLIT, derived from the corpus records rather than
    # asserted beside them. Round 4 replaced a deleted PMID by reading the
    # counting instrument instead of reading the article, and picked a
    # text-mining survey; round 5 wrote the correction as prose, which the
    # next reviewer flipped to "three primary reports and no survey" with the
    # whole suite green.
    surveyish = re.compile(
        r"text[- ]mining|bibliometric|knowledge structure|scientometric|"
        r"\bmapping\b.*\bthemes\b", re.I)
    dual_titles = {}
    for a in sn:
        if "ferroptosis" in a["_text"] and re.search(icd, a["_text"], re.I):
            dual_titles[pmid(a)] = a.get("title") or ""
    surveys = {k: v for k, v in dual_titles.items() if surveyish.search(v)}
    assert len(dual_titles) == facts["dual new"]
    # BOTH HALVES DERIVED. Pinning only "one survey" and the literal phrase
    # left the "two" free: adding a fourth dual-pathway article to the corpus
    # and updating the two numerals gave a paragraph claiming two primary
    # reports of a set holding three, with the guard green.
    primaries = len(dual_titles) - len(surveys)
    words = {0: "no", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}
    expected = (f"{words.get(primaries, primaries)} primary report"
                f"{'' if primaries == 1 else 's'} and "
                f"{words.get(len(surveys), len(surveys))} survey"
                f"{'' if len(surveys) == 1 else 's'}")
    assert expected in sentence, (
        f"the dual-pathway set is {primaries} primary and {len(surveys)} "
        f"survey by title, so the paragraph should read {expected!r}. "
        f"Titles: {sorted(dual_titles.values())}")
    assert "reconstruction rather than a record" in sentence, (
        "the paragraph states the retired instrument as fact; it is inferred "
        "from one number and the companion counts do not corroborate it")
    cut = subprocess.run(
        ["git", "log", "-1", "--format=%cs", "37ab0f55"],
        cwd=REPO, capture_output=True, text=True).stdout.strip()
    assert cut[:4] and cut[5:7], f"cannot date the cut commit: {cut!r}"
    month = {"01": "January", "02": "February", "03": "March", "04": "April",
             "05": "May", "06": "June", "07": "July", "08": "August",
             "09": "September", "10": "October", "11": "November",
             "12": "December"}[cut[5:7]]
    assert f"in {month} {cut[:4]}" in sentence, (
        f"the corpus was cut in {month} {cut[:4]} and the paragraph says "
        "something else")
    # The two rounded steps, computed from UNROUNDED shares. Subtracting the
    # rounded operands gives 2.2 and is what the first version of the sentence
    # published.
    assert f"about {pcts['widened old'] - pcts['retired']:.1f} points" in sentence
    assert f"{pcts['widened new'] - pcts['widened old']:.1f} from the corpus" in sentence


def _axis_patterns():
    """The GSH-axis and ICD patterns `classify_ferroptosis` actually uses."""
    src = GEN.read_text()
    body = src[src.index("def classify_ferroptosis"):]
    body = body[:body.index("\n    return results")]
    gsh = re.search(r'gsh = sum\(.*?re\.search\(\s*r"([^"]+)"', body, re.S)
    icd = re.search(r'icd = sum\(.*?re\.search\(\s*r"([^"]+)"', body, re.S)
    assert gsh and icd, (
        "cannot read the GSH/ICD patterns out of classify_ferroptosis; if it "
        "was restructured, update this reader rather than restating them")
    # The sentence names the axis components in prose, so they must agree.
    for token in ("glutathione", "GSH", "GPX4", "SLC7A11"):
        assert token.lower() in gsh.group(1).lower(), (
            f"the GSH axis no longer counts {token}, which Section 8.2 lists "
            "by name as part of the widened instrument")
    return gsh.group(1), icd.group(1)


def _chapter_8_sentence() -> str:
    """The retraction paragraph in the manuscript source, whitespace-normalised.

    BOUNDED to one paragraph inside Chapter 8. An unbounded `.*?` matched from
    Section 8.2 to a copy of the closing sentence parked seven hundred lines
    later under a different heading, swallowing the whole intervening
    manuscript and passing -- so the guard could be satisfied by text that had
    been moved out of the chapter it describes.
    """
    text = (REPO / "article/drafts/v1.md").read_text()
    start = text.index("## Chapter 8:")
    nxt = text.find("\n## ", start + 1)
    chapter = text[start:nxt if nxt != -1 else len(text)]
    paras = [p for p in chapter.split("\n\n") if "SDT-specific data shows" in p]
    assert len(paras) == 1, (
        f"expected one retraction paragraph inside Chapter 8, found {len(paras)}")
    para = " ".join(paras[0].split())
    # SENTENCE COUNT, so an inserted claim cannot hide between checked ones.
    # A reviewer added "The corpus was never actually cut; this paragraph is
    # fiction." mid-paragraph and every fragment assertion still passed --
    # presence checks cannot see what else is there.
    n_sentences = len([x for x in re.split(r"(?<=[.!?]) ", para) if x.strip()])
    assert n_sentences == 10, (
        f"the retraction paragraph has {n_sentences} sentences, not the 10 this "
        "guard was written against. If a sentence was added, add it to the "
        "checks below deliberately -- an unchecked sentence in a retraction is "
        "how this paragraph came to need a guard.")
    assert para.endswith("broader than SDT alone."), (
        "the retraction paragraph does not end where this guard expects; it "
        "may have been split, and half of it would then go unchecked")
    return para


def test_the_latex_manuscript_is_a_regeneration_of_its_source():
    """`v1.tex` must be exactly what `generate_latex.py` produces from `v1.md`.

    THE RELEASE ARTIFACT IS THE TEX. `release-pdf.yml` builds the published PDF
    from `article/drafts/v1.tex`, while every prose guard in this repository --
    including the retraction guard above -- reads `v1.md`. A reviewer reverted
    one sentence in the tex to the retracted "72 GSH/GPX4 articles (11.7%)",
    regenerated `MANIFEST.sha256` exactly as the manifest error message
    instructs, and the whole suite stayed green with the retracted number in
    the shipped LaTeX.

    That is the same .md/.tex split that let the fig1 caption drift for five
    months, one level up: the caption was fixed in the generated file while its
    generator still held the old value. Checking the two agree closes both.
    """
    import subprocess
    import tempfile

    tex = REPO / "article/drafts/v1.tex"
    committed = tex.read_text()
    with tempfile.TemporaryDirectory() as td:
        backup = Path(td) / "v1.tex"
        backup.write_text(committed)
        res = subprocess.run(
            [sys.executable, str(REPO / "scripts/generate_latex.py")],
            cwd=REPO, capture_output=True, text=True)
        regenerated = tex.read_text()
        tex.write_text(committed)          # never leave the tree modified
    assert res.returncode == 0, res.stderr[-600:]
    assert regenerated == committed, (
        "article/drafts/v1.tex is not what scripts/generate_latex.py draws "
        "from v1.md. Re-run it -- an edit made directly to the tex is an edit "
        "no prose guard in this repository can see, and it is what ships.")


def test_cited_titles_match_the_corpus_record():
    """A footnote's title must be the title the corpus holds.

    `34646381` was cited as "…through a mitochondrial targeting **strategy**"
    at three sites in the manuscript; the record says "…through a mitochondrial
    targeting **liposomal nanosystem**". Nothing tied a footnote to the article
    it names, so a paraphrase introduced during an edit stayed for as long as
    nobody opened the file.

    Scoped to the PMIDs Section 8.2 enumerates, which are the ones this PR
    rewrote; a repository-wide version of this check is a larger job (many
    cited PMIDs are not in the corpus at all) and is not attempted here.
    """
    import yaml

    text = (REPO / "article/drafts/v1.md").read_text()
    m = re.search(r"\[\^refs_group8\]:(.*?)(?=\n\[\^|\n\n)", text, re.S)
    assert m, "the dual-pathway footnote is gone or was renamed"
    footnote = " ".join(m.group(1).split())
    # PER SEGMENT, so a title is checked against the PMID it names. Searching
    # the whole footnote for each title made the titles interchangeable: two
    # could be swapped and every assertion still passed -- which is verbatim
    # what this test's docstring says it closes.
    segments = re.findall(r"(.*?PMID:\s*(\d+))", footnote)
    assert segments, f"no PMIDs in the footnote: {footnote[:120]}"
    for segment, pmid in segments:
        record = REPO / f"corpus/by-pmid/{pmid}.md"
        assert record.exists(), (
            f"the footnote cites PMID {pmid}, which is not in the corpus -- "
            "it was cited as evidence and the article is not held")
        fm = yaml.safe_load(
            re.match(r"^---\n(.*?\n)---", record.read_text(), re.S).group(1))
        title = (fm.get("title") or "").rstrip(".")
        assert title, f"{pmid} has no title in the corpus"
        # Compared on words, because the footnote drops the trailing period
        # and the corpus keeps it; the failure this catches is a paraphrase,
        # not punctuation.
        assert " ".join(title.split()).lower() in segment.lower(), (
            f"the footnote's entry for {pmid} does not carry that article's "
            f"corpus title:\n  corpus:  {title}\n  footnote: {segment.strip()}")


def test_fig6_states_the_window_it_actually_matches_over():
    """The caption's description of the instrument, checked against the code.

    Seven review rounds have found a false sentence in this one caption three
    times: it claimed the chain "thins from left to right" when it never did,
    then that the bars differ in strictness in a way they do not, then that
    matching is over "title or abstract" when `load_corpus` slices the first
    4,000 characters of each record's BODY -- a window that reaches
    `## Full Text` in 187 of 187 SDT records, so five of the eight bars would
    be lower if it were the abstract alone.

    Each replacement was written by the commit that removed the previous one.
    A caption nothing reads is prose beside a measurement, which is the defect
    this whole file exists to close, so this reads it.
    """
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            pytest.skip("no PDF reader available")

    src = GEN.read_text()
    m = re.search(r"fm\[\"_text\"\]|_text.*?body\[:(\d+)\]|body\[:(\d+)\]", src)
    assert m, "cannot find the text window in load_corpus"
    window = int(next(g for g in m.groups() if g))
    doc = pymupdf.open(REPO / "article/figures/fig6_sdt_chain_evidence.pdf")
    try:
        caption = " ".join(" ".join(p.get_text().split()) for p in doc)
    finally:
        doc.close()
    assert f"{window:,} body characters" in caption, (
        f"load_corpus matches over the title plus the first {window:,} "
        "characters of each record's BODY and the figure does not say so. "
        "(The record's own first characters are YAML frontmatter, which is "
        "not matched -- an earlier caption said 'first 4,000 characters of "
        "the record', which named the wrong span.)")
    for wrong in ("title or abstract matches",
                  "thins from left to right. Matching",
                  "require the molecule to be named"):
        assert wrong not in caption, (
            f"the caption carries a retracted claim: {wrong!r}")
