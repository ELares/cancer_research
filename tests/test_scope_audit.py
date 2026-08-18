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


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("sa", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_readme_states_the_scope_at_all():
    """The whole point: a reader must not have to count files."""
    r = README.read_text()
    assert "What the work is actually about" in r, (
        "the README no longer states what the work is about, so the scope is "
        "back to being discoverable only by counting files")
    assert "preregistered predictions" in r and "engine modules" in r


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
    # THE README ROW USED TO BE `N of N`, and this guard asserted exactly
    # that -- the same variable on both sides, so it could only ever pass.
    # It now pins the MEASURED pair.
    if isinstance(em, dict):
        assert f"**{em['mention']} of {em['modules']}**" in r, (
            f"the README does not state the measured module count "
            f"{em['mention']} of {em['modules']}")
        assert f"**{em['in_code']} of {em['modules']}**" in r, (
            "the README does not state the stricter in-code count")
    else:
        assert f"**{em} of {em}**" in r


def test_the_mechanical_counts_are_right():
    """Predictions and modules need no judgement, so they get checked directly."""
    d = _doc()
    preds = re.findall(r"^\*\*(P\d+)\.", (REPO_ROOT / "PREREGISTRATION.md").read_text(), re.M)
    assert len(preds) == d["n_predictions"], (
        f"the audit parsed {d['n_predictions']} predictions, "
        f"PREREGISTRATION.md states {len(preds)}")
    src = REPO_ROOT / "simulations" / "ferroptosis-core" / "src"
    em = d["engine_modules"]
    n_mods = len([f for f in src.glob("*.rs") if f.stem != "lib"])
    # `lib.rs` is the crate root, not a module. This guard used to assert the
    # count equalled the FILE count, which is what made the `N of N` row
    # unfalsifiable: it pinned a file count as if it were a content measure.
    got = em["modules"] if isinstance(em, dict) else em
    assert got == n_mods, (
        f"the audit counts {got} engine modules; the tree holds {n_mods} "
        "(excluding the crate root lib.rs)")


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


def test_the_artifact_is_fresh_against_a_live_classification():
    """A REAL regenerate-and-diff gate.

    Every other guard here compares the .md to the .json, and the two go stale
    TOGETHER: the shipped row said 13/1/95 while a live run gave 14/1/103, the
    README carried the stale row on the front door, and all fifteen guards
    passed. Nine analyses landed -- this campaign's own -- and nothing fired.

    Classification is a walk over ~118 small files and costs well under a
    second, so there is no reason for this gate not to exist.
    """
    m = _mod()
    live = m.classify_analyses()
    committed = _doc()["analyses"]
    for bucket in ("ferroptosis-or-physical", "therapy-subject", "method"):
        got, want = sorted(committed.get(bucket, [])), sorted(live[bucket])
        assert got == want, (
            f"`{bucket}` differs between a live classification and the "
            f"committed artifact -- analysis/scope-audit.json is stale. "
            f"Re-run `python scripts/scope_audit.py`.\n"
            f"  only live      = {sorted(set(want) - set(got))[:8]}\n"
            f"  only committed = {sorted(set(got) - set(want))[:8]}")
    n = len(list((REPO_ROOT / "analysis").glob("*.md")))
    tot = sum(len(committed.get(b, [])) for b in
              ("ferroptosis-or-physical", "therapy-subject", "method"))
    assert tot == n, (
        f"the artifact classifies {tot} analyses and analysis/ holds {n}")


def test_the_engine_module_row_measures_content():
    """`N of N` is the same variable twice and cannot come out otherwise."""
    m, d, md = _mod(), _doc(), MD.read_text()
    em = d["engine_modules"]
    assert isinstance(em, dict), (
        "engine_modules is a bare count again; the row it feeds printed "
        "`N of N`, arithmetic that cannot fail, from a file count that never "
        "opened a file")
    src = REPO_ROOT / "simulations" / "ferroptosis-core" / "src"
    mods = [f for f in sorted(src.glob("*.rs")) if f.stem != "lib"]
    assert em["modules"] == len(mods), (
        f"the audit counts {em['modules']} modules, the tree holds {len(mods)}")
    # re-derive the content measure independently
    silent = [f.name for f in mods
              if not m.FERRO.search(f.read_text(errors="ignore"))]
    assert sorted(em["silent"]) == sorted(silent), (
        f"the audit reports {sorted(em['silent'])} as silent; an independent "
        f"scan gives {sorted(silent)}")
    assert em["mention"] == len(mods) - len(silent)
    assert em["in_code"] <= em["mention"], (
        "more modules mention it in code than mention it at all")
    # and the page must not claim universality when it is not universal
    if em["mention"] < em["modules"]:
        assert "every module of its simulation engine, concerns" not in md, (
            f"{len(silent)} modules mention neither ferroptosis nor a "
            "physical-ROS modality, and the page still claims every one does")
        assert "mostly but not entirely" in md
        # and the paragraph NAMING them: hiding it left the reader with a
        # count and no way to check it
        assert f"**{len(silent)} modules mention neither" in md, (
            "the report states the shortfall as a number and no longer names "
            "which modules it is")
        for name in silent:
            assert f"`{name}`" in md, (
                f"`{name}` mentions neither and the report does not name it")


def test_the_headline_sentence_follows_the_table():
    """It was unconditional: a doctored input printed the same claim."""
    m, d = _mod(), _doc()
    p = d["predictions"]
    doctored = {**d,
                "predictions": {k: (i % 2 == 0) for i, k in enumerate(p)},
                "engine_modules": {**d["engine_modules"],
                                   "mention": 1, "silent": ["x.rs"]}}
    out = m.render(doctored)
    assert "Every falsifiable commitment the project makes, and every module" \
        not in out, (
        "with half the predictions non-ferroptosis and one module mentioning "
        "it, the page still prints the universal claim -- the sentence is not "
        "a function of the table beside it")
    assert f"{sum(doctored['predictions'].values())} of {len(p)}" in out


def test_the_therapy_asymmetry_is_disclosed_with_both_numbers():
    """`1` is a filename marker; the naive symmetric fix over-claims."""
    d, md = _doc(), MD.read_text()
    body = d["analyses"].get("_therapy_by_body")
    assert body is not None, (
        "the body-route therapy count is gone, so the asymmetry is asserted "
        "rather than measured")
    n_ther = len(d["analyses"]["therapy-subject"])
    assert f"gives **{len(body)}** analyses instead of {n_ther}" in md, (
        "the report does not state both counts, so a reader cannot see how "
        "much the admission rule is worth")
    assert "neither rule measures subject" in md, (
        "the report presents one of the two as the corrected number")


def test_the_mechanism_denominators_are_derived_and_named():
    """Both this page and the README called corpus shares 'tagged' shares."""
    d, md = _doc(), MD.read_text()
    m = d.get("mechanism_denominators")
    if not m:
        return
    idx = REPO_ROOT / "corpus" / "INDEX.jsonl"
    recs = tagged = tags = 0
    counts = {}
    with idx.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            recs += 1
            ms = r.get("mechanisms") or []
            tagged += bool(ms)
            tags += len(ms)
            for x in ms:
                counts[x] = counts.get(x, 0) + 1
    assert m["corpus_records"] == recs and m["tagged"] == tagged and m["tags"] == tags, (
        "the stored denominators disagree with a live count of INDEX.jsonl")
    assert f"{100*m['share_of_corpus']:.1f}% of corpus articles" in md
    assert f"{100*m['share_of_tagged']:.1f}% of the {m['tagged']:,} tagged" in md
    assert "not shares of the cancer literature" in md
    # and it must carry the retraction: reverting the wording to "shares of
    # TAGGED articles" left every other assertion satisfied
    for mt in re.finditer(r"shares of TAGGED articles", md):
        w = md[max(0, mt.start() - 300):mt.end() + 300]
        assert re.search(r"earlier version|they are not|An earlier", w), (
            "the page calls the mechanism shares 'shares of TAGGED articles' "
            "again; they are shares of the corpus, and the tagged-article "
            "figure is a different number")
    assert "An earlier version of this bullet" in md, (
        "the correction is no longer recorded, so a reader cannot tell the "
        "denominator was wrong")
