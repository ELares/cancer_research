"""Guards for the descriptor-recall asymmetry (#725 / #MS-CENSUS).

THE FINDING
-----------
`Photochemotherapy` and `Ultrasonic Therapy` are comparably PRECISE and differ
markedly in RECALL, so the PDT:SDT ratio `manuscript-vs-census.md` reports is
substantially a measurement of indexing practice. Under one rule applied to
both arms the census cannot be DISTINGUISHED from the manuscript, which
withdraws the "understated by the manuscript" verdict without establishing
agreement -- a distinction an earlier version of this docstring collapsed by
saying the census "REPRODUCES" it.

WHAT WOULD MAKE THIS ANALYSIS ITSELF WRONG
--------------------------------------------
1. AN ASYMMETRIC TEXT RULE. The entire finding is that one arm was measured
   differently from another. If the replacement rule is itself lopsided --
   a broader stem, an extra alternative, an acronym on one side only -- this
   document commits the defect it exists to correct. Guarded structurally,
   not by eye.

2. A VERDICT THAT CANNOT FLIP. If the symmetric ratio had come out ABOVE the
   manuscript's, the report must say the manuscript understated its case. A
   fixed sentence would make the conclusion unfalsifiable.

3. THE INVERTED CAVEAT RETURNING. Nine sites said the SDT count is an
   OVER-estimate and ratios against it are LOWER bounds -- an earlier version
   of this docstring said five, a hand-written count that was itself wrong,
   and two review passes were needed to find the rest. The breadth is real
   and small; the recall shortfall is larger. The direction is settled by
   measurement now and must stay that way.

4. TREATING THE TEXT RULE AS TRUTH. It is not; it is merely applied
   identically to both arms. The report must keep saying so.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_descriptor_recall.py"
MD = REPO_ROOT / "analysis" / "atlas-descriptor-recall.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-descriptor-recall.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("adr", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_text_rule_is_structurally_symmetric():
    """The finding is an asymmetric comparison; the fix must not be one.

    Checked on SHAPE, because the arms necessarily differ in stem. Each must
    contribute the same number of alternatives, and each must carry the same
    three kinds: a plain stem, a hyphenated stem, and a word-bounded acronym.
    """
    m = _mod()
    shapes = {}
    for k, arm in m.ARMS.items():
        alts = arm["text"].split("|")
        shapes[k] = {
            "n": len(alts),
            "hyphenated": sum("-" in a for a in alts),
            "acronym": sum(a.startswith(r"\b") and a.endswith(r"\b") for a in alts),
            "plain": sum("-" not in a and not a.startswith(r"\b") for a in alts),
        }
    first = next(iter(shapes.values()))
    for k, s in shapes.items():
        assert s == first, (
            f"arm {k}'s text rule has shape {s}, others have {first}; an "
            "asymmetric rule reproduces the exact defect this analysis "
            "corrects")
    assert first["n"] >= 3 and first["plain"] and first["acronym"], (
        f"the shared rule shape {first} is too thin to measure recall with")


def test_the_text_rule_is_symmetric_in_COVERAGE_not_only_in_shape():
    """Shape symmetry is necessary and nowhere near sufficient.

    A reviewer changed the PDT stem `photodynamic` to `photo` -- same
    alternative count, same three shapes, same descriptor-set size -- and the
    shape guard passed while PDT text counts went 182 -> 509, the recall
    asymmetry inverted, and the report's headline flipped back to "understated
    by the manuscript". A twelve-character edit reversed the finding with
    every guard green.

    Coverage is checked two ways, both against the corpus rather than the
    rule: each arm's stem must be SPECIFIC (nearly every article it matches
    must also carry that arm's descriptor or the sibling modality vocabulary),
    and the arms' specificities must be comparable.
    """
    import gzip
    m = _mod()
    # THE OFFLINE CONTRACT. corpus/atlas/records/ is gitignored bulk data, so
    # CI has no shards: this guard recounts against the census and can only
    # run where the census exists. Skipping keeps the contract while the
    # check still fires for anyone holding the data -- which is where the
    # scan-level mutations it exists to catch would be introduced.
    import pytest
    if not any((m.ATLAS / "records").glob("*.jsonl.gz")):
        pytest.skip("census shards not present in this checkout")

    pats = {k: re.compile(v["text"], re.I) for k, v in m.ARMS.items()}
    descs = {k: set(v["descriptors"]) for k, v in m.ARMS.items()}
    hit = {k: 0 for k in m.ARMS}
    corroborated = {k: 0 for k in m.ARMS}
    # A STRIDE, not a prefix. Shards are chronological, so the first N are the
    # oldest literature: taking a prefix gave SDT zero hits, because
    # sonodynamic therapy is recent. This repo already recorded that trap once
    # (#722, where a prefix sampled only the oldest articles) and it was
    # reproduced here within the hour.
    every = sorted((m.ATLAS / "records").glob("*.jsonl.gz"))[::12]
    shards = every
    assert len(shards) > 40, f"only {len(shards)} shards sampled"
    for f in shards:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                mesh = {x.lower() for x in (r.get("mesh") or [])}
                blob = (r.get("title") or "") + " " + (r.get("abstract") or "")
                for k in m.ARMS:
                    if pats[k].search(blob):
                        hit[k] += 1
                        if mesh & descs[k]:
                            corroborated[k] += 1
    spec = {k: (corroborated[k] / hit[k]) if hit[k] else None for k in hit}
    measurable = {k: v for k, v in spec.items() if v is not None and hit[k] >= 30}
    assert len(measurable) == len(m.ARMS), (
        f"an arm matched too few articles to judge its specificity: {hit}")
    lo = min(measurable.values())
    assert lo > 0.25, (
        f"arm specificities {({k: round(v, 3) for k, v in measurable.items()})} "
        f"-- one stem matches far more articles than carry its descriptor, so "
        "it is a broader concept than the arm it stands for and the ratio "
        "would measure the stem rather than the modality")
    ratio = max(measurable.values()) / lo
    assert ratio < 2.0, (
        f"the two arms' stems differ {ratio:.2f}x in how specifically they "
        f"pick out their own descriptor ({({k: round(v, 3) for k, v in measurable.items()})}); "
        "that is an asymmetric rule in coverage even though it is symmetric "
        "in shape, which is the defect this analysis exists to correct")


def test_each_arm_has_exactly_one_descriptor_so_breadth_is_comparable():
    m = _mod()
    sizes = {k: len(a["descriptors"]) for k, a in m.ARMS.items()}
    assert len(set(sizes.values())) == 1, (
        f"the arms carry different numbers of descriptors {sizes}; a wider "
        "descriptor set on one side is the asymmetry in another form")


def test_recall_and_precision_are_arithmetically_consistent():
    d = _doc()
    for k, s in d["arms"].items():
        assert s["both"] <= s["text"], f"{k}: both exceeds text"
        assert s["both"] <= s["descriptor"], f"{k}: both exceeds descriptor"
        assert abs(s["recall"] - s["both"] / s["text"]) < 1e-9
        assert abs(s["precision"] - s["both"] / s["descriptor"]) < 1e-9


def test_the_asymmetry_is_in_recall_not_precision():
    """If precision ever became the lopsided axis, the argument changes."""
    d = _doc()
    recs = [s["recall"] for s in d["arms"].values()]
    precs = [s["precision"] for s in d["arms"].values()]
    r_gap = max(recs) / min(recs)
    p_gap = max(precs) / min(precs)
    assert r_gap > p_gap, (
        f"recall gap {r_gap:.2f}x no longer exceeds the precision gap "
        f"{p_gap:.2f}x; the finding that the descriptors differ in RECALL "
        "rather than breadth needs re-stating")
    # RE-DERIVED from the per-arm counts, not read from `recall_asymmetry`.
    # Asserting the artifact's own field appears in the report is a guard
    # computing its own expectation: hardcoding the field to 1.00 moved both
    # sides together and survived.
    truth = max(recs) / min(recs)
    assert abs(d["recall_asymmetry"] - truth) < 1e-9, (
        f"the artifact reports a recall asymmetry of {d['recall_asymmetry']} "
        f"and its own per-arm counts give {truth:.4f}")
    assert f"factor of **{truth:.2f}**" in MD.read_text(), (
        "the report does not state the asymmetry its own counts imply")


def test_the_verdict_follows_the_measurement_and_can_flip():
    """A conclusion that cannot come out the other way is not a conclusion."""
    m, d = _mod(), _doc()
    md = MD.read_text()
    text_r, ms = d["ratio_by_text"], d["manuscript_ratio"]
    ci = d["ratio_by_text_ci"]
    # Decided by the INTERVAL, not by a point comparison. An earlier version
    # branched on `text_r <= ms * 1.05` and demanded the word "reproduces" --
    # a 5% window on a ratio whose 95% interval spans roughly [2.2, 3.9], so
    # the verdict rode a threshold the data could not resolve.
    covers = ci and ci[0] <= ms <= ci[1]
    assert d["symmetric_ratio_covers_manuscript"] == bool(covers), (
        "the artifact's coverage flag disagrees with its own interval")
    if covers:
        assert "cannot distinguish" in md, (
            "the manuscript's ratio lies inside the symmetric interval and "
            "the report does not say the two are indistinguishable")
        # CLAIMING form, not substring: the report necessarily QUOTES
        # "reproduces the manuscript's ratio" in the sentence withdrawing it.
        # A bare ban trips on the retraction -- the sixth time that trap has
        # fired in this repo.
        for mt in re.finditer(r"reproduces the manuscript", md):
            window = md[max(0, mt.start() - 250):mt.end() + 250]
            assert re.search(r"withdraw|earlier draft|not support", window, re.I), (
                "the report claims reproduction as a live claim; overlapping "
                "intervals establish only that the two cannot be told apart")
    # and it must render differently when the numbers say otherwise
    flipped = {**d, "ratio_by_text": ms * 3}
    out = m.render(flipped)
    assert f"{ms * 3:.2f}:1" in out, (
        "the rendered ratio is not read from the data, so the verdict cannot "
        "follow a different measurement")


def test_the_inverted_over_estimate_caveat_does_not_return():
    """Measured: breadth is 3 records, the recall shortfall is an order more.

    Checked as a CLAIMING form. The report necessarily QUOTES the withdrawn
    wording in order to withdraw it, so banning the substring outright trips
    on the correction -- the trap this repo has hit five times now.
    """
    md = MD.read_text()
    for m in re.finditer(r"OVER-estimate|over-estimate", md):
        window = md[max(0, m.start() - 300):m.end() + 300]
        assert re.search(r"under|invert|withdraw|runs the other way|state that",
                         window, re.I), (
            "the report states the SDT count is an over-estimate without the "
            "correction that recall runs the other way")
    d = _doc()
    sdt = d["arms"]["SDT"]
    breadth = sdt["descriptor"] - sdt["both"]
    shortfall = sdt["text"] - sdt["both"]
    assert shortfall > breadth, (
        f"the recall shortfall ({shortfall}) no longer exceeds the breadth "
        f"({breadth}); the direction of the caveat must be re-derived")


def test_the_text_rule_is_not_presented_as_truth():
    md = MD.read_text()
    assert "applied identically to both" in md.lower(), (
        "the report no longer states that the text rule's value is its "
        "symmetry rather than its accuracy")
    assert "Not that the text rule is ground truth" in md


def test_it_does_not_claim_the_manuscript_is_wrong():
    """The manuscript's own figure is reproduced; only a verdict is withdrawn."""
    md = MD.read_text()
    assert "Not that the manuscript is wrong" in md
    assert re.search(r"understat", md, re.I), (
        "the report no longer names the verdict it withdraws, so a reader "
        "cannot tell what changed")
    # and it must NOT claim the census confirms the manuscript, which is a
    # different and unsupported statement -- the intervals overlap, they do
    # not coincide
    assert "does not establish agreement" in md.lower() or \
           "cannot distinguish" in md.lower(), (
        "the report does not distinguish 'cannot tell them apart' from "
        "'they agree'; the data supports only the former")


def test_an_unmeasurable_arm_refuses_to_render():
    src = SCRIPT.read_text()
    assert "matched nothing on one axis" in src
    assert "raise SystemExit" in src
    assert "is not a finding" in src


def test_the_two_artifacts_describe_the_same_census_build():
    """A free fingerprint nobody was comparing.

    Both documents count the same ferroptosis subset. Setting this analysis's
    `subject_articles` to a fabricated value left the page reading "Over the
    99,999 census articles" beside the sibling's 13,346, with every guard
    green -- the `pairs_before` lesson, unlearned one analysis over.
    """
    d = _doc()
    sib = REPO_ROOT / "analysis" / "manuscript-vs-census.json"
    if not sib.exists():
        return
    other = json.loads(sib.read_text())
    theirs = other.get("ferroptosis_records") or other.get("subject_articles")
    if theirs is None:
        return
    assert d["subject_articles"] == theirs, (
        f"this analysis counted {d['subject_articles']:,} ferroptosis "
        f"articles and manuscript-vs-census counted {theirs:,}; one of the "
        "two artifacts is from a different census build, and this document's "
        "recall figures qualify that one's verdict")
    assert f"{d['subject_articles']:,}" in MD.read_text()


def test_the_recall_table_matches_the_artifact_row_by_row():
    """The central table was entirely unguarded.

    Printing the same recall for both arms, while the sentence below still
    stated a 1.74x gap, shipped a self-contradicting page with the suite
    green.
    """
    d, md = _doc(), MD.read_text()
    for k, s in d["arms"].items():
        row = (f"| {k} | `{s['descriptors'][0].title()}` | {s['text']:,} | "
               f"{s['descriptor']:,} | {s['both']:,} | "
               f"**{100*s['recall']:.1f}%** | {100*s['precision']:.1f}% |")
        assert row in md, (
            f"the rendered row for {k} does not match its artifact counts; "
            f"expected:\n  {row}")
    recs = [s["recall"] for s in d["arms"].values()]
    if abs(max(recs) - min(recs)) < 1e-12:
        assert "factor of" not in md, (
            "the table shows equal recalls while the report still states a "
            "gap between them")


def test_the_manuscript_ratio_constant_matches_the_manuscript():
    """A free-floating duplicate of a figure the sibling pins to v1.md.

    Changing it from 2.93 to 2.50 reverted the headline to "understated" with
    every guard green, because nothing tied this constant to anything.
    """
    m = _mod()
    # The figure comes from the manuscript's Section 8.2 TABLE, whose counts
    # `manuscript_vs_census.MODALITIES` holds and whose own guards pin to
    # v1.md. Deriving it from the census-intersection sentence instead gives
    # 5.5 and is a different claim -- a mistake made writing this guard.
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mvc", REPO_ROOT / "scripts" / "manuscript_vs_census.py")
    mvc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mvc)
    counts = {row[0]: row[2] for row in mvc.MODALITIES}
    if not (counts.get("PDT") and counts.get("SDT")):
        return
    implied = counts["PDT"] / counts["SDT"]
    assert abs(m.MANUSCRIPT_RATIO - implied) < 0.01 * implied, (
        f"MANUSCRIPT_RATIO is {m.MANUSCRIPT_RATIO} but the manuscript's own "
        f"Section 8.2 counts ({counts['PDT']}/{counts['SDT']}) imply "
        f"{implied:.2f}; this constant decides the verdict and was "
        "free-floating")
