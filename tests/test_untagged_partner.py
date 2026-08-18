"""Guards for the untagged-partner finding (#723).

THE CLAIM
---------
A modality the mechanism taxonomy cannot name does not become an untagged
article. It becomes someone else's article: it arrives attached to a modality
that DOES have a tag, and the whole paper is filed under the partner.

Measured on the frozen corpus, for articles TITLED about each modality:
radiotherapy 75 of 88 filed under another modality, chemotherapy 140 of 166,
surgery 17 of 24. Immunotherapy is the dominant absorber in all three.

WHY THIS NEEDS GUARDING RATHER THAN JUST STATING
-------------------------------------------------
1. IT WOULD BE EASY TO MANUFACTURE. Widen the title pattern and the subject
   count rises; widen it enough and it sweeps in papers that merely mention the
   modality, at which point "filed under another modality" is trivially true
   because the paper really is about the partner. The patterns must stay
   specific, and the finding must hold for MORE THAN ONE modality -- a single
   modality could be a quirk of one pattern.

2. THE THRESHOLD IS A JUDGEMENT. Five mentions was fixed before the result was
   seen. The any-mention column has to stay so a reader can apply their own,
   and the threshold has to stay stated.

3. IT MUST NOT SILENTLY INSTALL A TAG. Adding a radiotherapy lane would move
   articles' existing partner tags and change figures the manuscript quotes.
   This analysis measures; adopting is a separate act. A guard checks that the
   production vocabulary is untouched.
"""

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_untagged_partner.py"
MD = REPO_ROOT / "analysis" / "atlas-untagged-partner.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-untagged-partner.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("up", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_finding_holds_for_more_than_one_modality():
    """One modality could be an artifact of one pattern. Three is a pattern."""
    d = _doc()
    mods = d["modalities"]
    assert len(mods) >= 3, f"only {len(mods)} modalities measured"
    absorbed = {m: s["subject_tagged_as_other"] / max(s["subject_titled"], 1)
                for m, s in mods.items() if s["subject_titled"] >= 20}
    assert len(absorbed) >= 3, (
        "fewer than three modalities have enough titled articles to judge")
    for m, frac in absorbed.items():
        assert frac > 0.5, (
            f"{m}: only {100*frac:.0f}% of titled articles are filed under "
            "another modality, so the partner-attribution claim does not hold "
            "for it")


def test_the_arithmetic_partitions():
    """Untagged plus other-tagged must equal the titled set."""
    for m, s in _doc()["modalities"].items():
        assert s["subject_untagged"] + s["subject_tagged_as_other"] == \
            s["subject_titled"], (
            f"{m}: the two buckets do not sum to the titled count")
        assert s["strong_fulltext"] <= s["any_fulltext"], (
            f"{m}: more articles pass the mention threshold than mention it")


def test_the_patterns_stay_specific():
    """A wide pattern manufactures the finding.

    Checked by running the real patterns against strings that must NOT match.
    A pattern loose enough to catch these would inflate the titled set with
    papers that merely touch the modality, and the partner-attribution result
    would then be trivially true.
    """
    m = _mod()
    must_not = [
        "Radiation-induced bystander effects in normal tissue",
        "A chemical screen identifies novel inhibitors",
        "Surgical margins were not assessed in this cohort",
    ]
    rt = m.CANDIDATES["radiotherapy"]
    assert not rt.search(must_not[1]), (
        "the radiotherapy pattern fires on an unrelated chemical-screen title")
    # and each pattern must be more than a bare stem
    for name, pat in m.CANDIDATES.items():
        assert len(pat.pattern) > 30, (
            f"{name}'s pattern is a bare stem; it will over-match")


def test_it_does_not_install_a_tag():
    """Measuring is not adopting; adopting moves numbers the manuscript quotes."""
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    import config
    for banned in ("radiotherapy", "chemotherapy", "surgery"):
        assert banned not in config.MECHANISM_KEYWORDS, (
            f"a `{banned}` lane has been installed in MECHANISM_KEYWORDS. That "
            "is a legitimate decision, but it changes the corpus tagging and "
            "the figures built on it, so it needs its own review -- and this "
            "analysis's numbers describe the world before it.")
    md = MD.read_text()
    assert "does not install a tag" in md, (
        "the report no longer says it is measuring rather than adopting")


def test_the_threshold_is_stated_and_the_alternative_is_shown():
    """A judgement call must be visible and a reader must be able to redo it."""
    d, md = _doc(), MD.read_text()
    assert d["strong_threshold"] >= 2
    assert f"{d['strong_threshold']}-mention threshold" in md, (
        "the report does not state the mention threshold it used")
    assert '"any_fulltext": len(fulltext_counts[name])' in SCRIPT.read_text(), (
        "the any-mention column is no longer the unfiltered count, so a reader "
        "cannot apply a threshold other than the one this script chose")
    for m, s in d["modalities"].items():
        assert s["any_fulltext"] > s["strong_fulltext"], (
            f"{m}: the any-mention column adds nothing, so a reader cannot "
            "apply a different threshold")


def test_the_partner_tags_are_named_not_just_counted():
    """'Filed under another modality' is only useful if you say which."""
    d, md = _doc(), MD.read_text()
    for m, s in d["modalities"].items():
        if s["subject_tagged_as_other"] == 0:
            continue
        assert s["partner_tags"], f"{m}: no partner tags recorded"
        # Pinned in the SOURCE as well as the artifact. This is the fourth
        # time today a guard has read only the committed output, which does not
        # move when the generator is edited -- so the mutation that reverted
        # this to a dict passed every assertion below.
        assert '"partner_tags": [[k, v] for k, v in tags.most_common(8)]' \
            in SCRIPT.read_text(), (
            "partner_tags is no longer built as an ordered list of pairs; a "
            "dict is reordered alphabetically by sort_keys and the 'dominant "
            "partner' becomes whichever tag sorts first")
        assert isinstance(s["partner_tags"], list), (
            f"{m}: partner_tags is a dict, so sort_keys reorders it "
            "alphabetically and the 'dominant partner' is whichever tag sorts "
            "first rather than whichever is most common")
        counts = [c for _t, c in s["partner_tags"]]
        assert counts == sorted(counts, reverse=True), (
            f"{m}: partner tags are not in descending order")
        top = s["partner_tags"][0][0]
        assert f"`{top}`" in md, (
            f"{m}'s dominant partner tag `{top}` is not named in the report")


def test_the_report_does_not_claim_this_is_evidence_about_the_field():
    """The corpus was retrieved by 33 queries, none about these modalities."""
    md = MD.read_text()
    assert ("not evidence about the field" in md
            or "arrived attached to something else" in md), (
        "the report no longer states that corpus presence says nothing about "
        "the field, which is the caveat that stops this being read as a "
        "prevalence measurement")


def test_an_empty_result_refuses_to_render():
    src = SCRIPT.read_text()
    assert "is not a finding" in src and "raise SystemExit" in src
    assert 'subject_titled"] == 0' in src, (
        "the empty check no longer tests the titled count")


def test_the_corpus_base_rate_is_measured_and_published():
    """The control this analysis never ran.

    Carrying another modality's tag is the corpus DEFAULT. Without the
    default, a rate of 85% reads as a finding when the corpus runs at 88.7%.
    """
    d, md = _doc(), MD.read_text()
    n, ct = d["corpus_size"], d.get("corpus_tagged")
    assert ct, "the corpus tagged count is not measured"
    # RE-DERIVED from INDEX.jsonl, not read back from the artifact
    idx = [json.loads(l) for l in (REPO_ROOT / "corpus" / "INDEX.jsonl").open()
           if l.strip()]
    want = sum(1 for r in idx if (r.get("mechanisms") or []))
    want_multi = sum(1 for r in idx if len(r.get("mechanisms") or []) >= 2)
    assert ct == want, f"artifact says {ct} tagged, INDEX.jsonl gives {want}"
    assert d.get("corpus_multi_tagged") == want_multi
    assert f"**the corpus** | {n:,}" in md, (
        "the corpus row is not in the table, so the per-modality rates have "
        "no denominator beside them")
    assert f"({100*ct/n:.1f}%)" in md


def test_the_causal_reading_follows_the_comparison():
    """Every modality sits BELOW the default, so nothing is elevated."""
    d, md = _doc(), MD.read_text()
    n, ct = d["corpus_size"], d["corpus_tagged"]
    base = 100 * ct / n
    rates = {m: 100 * s["subject_tagged_as_other"] / s["subject_titled"]
             for m, s in d["modalities"].items()}
    below = [m for m, r in rates.items() if r < base]
    if len(below) == len(rates):
        assert "The causal reading is withdrawn" in md, (
            f"all of {sorted(rates)} sit below the corpus default of "
            f"{base:.1f}% and the report still presents the rate as a finding")
        assert "sits BELOW the corpus default" in md or \
               "sit BELOW the corpus default" in md
        # the withdrawn sentence must not stand as a live claim
        for m_ in re.finditer(r"becomes \*\*someone else's article\*\*", md):
            w = md[max(0, m_.start() - 400):m_.end() + 400]
            assert re.search(r"withdraw|earlier claim|not supported", w, re.I)
    else:
        assert "The causal reading is withdrawn" not in md


def test_the_multi_tag_rate_is_reported_against_the_corpus():
    """The combination-paper story predicts MORE tags; they carry fewer."""
    d, md = _doc(), MD.read_text()
    n, cm = d["corpus_size"], d.get("corpus_multi_tagged")
    if not cm:
        return
    for m, s in d["modalities"].items():
        mt = s.get("subject_multi_tagged")
        assert mt is not None, f"{m} has no multi-tag count"
        t_ = s["subject_titled"]
        assert f"{mt:,} ({100*mt/t_:.1f}%)" in md, (
            f"{m}'s multi-tag rate is measured and not rendered")
    rates = [100 * s["subject_multi_tagged"] / s["subject_titled"]
             for s in d["modalities"].values()]
    if all(r < 100 * cm / n for r in rates):
        assert "carry FEWER tags" in md, (
            "every modality carries fewer tags than the corpus average, which "
            "contradicts the combination-paper reading, and the report does "
            "not say so")


def test_retired_tags_are_named_from_config_not_inferred():
    """A tag can also vanish because the re-run reads less text.

    Attributing a scope effect to the vocabulary would be the asymmetric
    comparison this repo keeps making, so only tags config.py records as
    narrowed may be listed.
    """
    d, md = _doc(), MD.read_text()
    src = (REPO_ROOT / "scripts" / "config.py").read_text()
    assert "#MECH-PRECISION" in src, (
        "config.py no longer records the keyword retirements this section "
        "depends on")
    for m, s in d["modalities"].items():
        old = {k for k, _v in s.get("partner_tags", [])}
        new = {k for k, _v in s.get("partner_tags_current_vocabulary", [])}
        for tag in sorted(old - new):
            if f"`{tag}`" in md.split("since withdrawn")[-1][:800]:
                assert re.search(rf"#MECH-PRECISION.{{0,600}}\b{re.escape(tag)}\b",
                                 src, re.S), (
                    f"`{tag}` is listed as withdrawn under #MECH-PRECISION and "
                    "config.py does not record it as narrowed; it may simply "
                    "have dropped out because this pass reads less text")
    assert "mix the vocabulary with the text scope" in md, (
        "the report no longer discloses that the re-run reads a narrower text "
        "than the frozen tagging, so its deltas would read as vocabulary")
