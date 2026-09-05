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
    # THIS ASSERTION HAS NOW FIRED, and it fired for the reason it was written:
    # P9 to P13 registered falsifiable predictions for the modality arms, so
    # the denominator widened for the first time. What replaces it is the same
    # property one level out -- the count must still be a MEASUREMENT, so it
    # must be able to move in both directions.
    assert d["n_predictions_ferroptosis"] < d["n_predictions"], (
        "every preregistered prediction concerns ferroptosis again; either the "
        "modality predictions were removed or the classifier stopped seeing "
        "them, and the widened denominator this project worked for is gone")
    assert d["n_predictions_ferroptosis"] >= 8, (
        "the ferroptosis predictions have shrunk, which would make the ratio "
        "improve by deletion rather than by addition")
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
    live, live_body = m.classify_analyses()
    committed = _doc()["analyses"]
    for bucket in ("ferroptosis-or-physical", "therapy-subject", "method"):
        got, want = sorted(committed.get(bucket, [])), sorted(live[bucket])
        assert got == want, (
            f"`{bucket}` differs between a live classification and the "
            f"committed artifact -- analysis/scope-audit.json is stale. "
            f"Re-run `python scripts/scope_audit.py`.\n"
            f"  only live      = {sorted(set(want) - set(got))[:8]}\n"
            f"  only committed = {sorted(set(got) - set(want))[:8]}")
    committed_body = _doc().get("therapy_by_body")
    assert committed_body is not None, "the body-route list is gone"
    assert sorted(committed_body) == sorted(live_body), (
        "`therapy_by_body` differs between a live classification and the "
        "committed artifact -- the number carrying the whole asymmetry "
        "section had no freshness gate, so it could go stale exactly the way "
        "the headline row did.\n"
        f"  only live      = {sorted(set(live_body) - set(committed_body))[:6]}\n"
        f"  only committed = {sorted(set(committed_body) - set(live_body))[:6]}")
    # Counted through the generator's own discovery, not a fresh disk glob.
    # A committed inventory describes what is IN THE REPOSITORY, so an
    # untracked or gitignored page sitting in analysis/ must not be expected in
    # it -- and reimplementing the rule here is how the two drift apart. The
    # concrete case is the expansion crawl's progress dashboard: gitignored
    # deliberately, present on a developer's disk, absent in CI.
    n = len(m._tracked_analyses())
    tot = sum(len(committed.get(b, [])) for b in
              ("ferroptosis-or-physical", "therapy-subject", "method"))
    assert tot == n, (
        f"the artifact classifies {tot} analyses and the repository tracks {n}")


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

    # RE-DERIVED, because artifact-versus-page cannot see a scan that stopped
    # stripping comments -- both sides move together and the page published a
    # figure four higher than the truth with everything green.
    def _strip_comments(s):
        return "\n".join(ln for ln in s.split("\n")
                          if not ln.strip().startswith(("//", "*", "/*")))

    want_code = want_prod = 0
    for f in mods:
        body = f.read_text(errors="ignore")
        if m.FERRO.search(_strip_comments(body)):
            want_code += 1
        i = body.find("#[cfg(test)]")
        if m.FERRO.search(_strip_comments(body if i < 0 else body[:i])):
            want_prod += 1
    assert em["in_code"] == want_code, (
        f"the audit reports {em['in_code']} modules mentioning it in code; an "
        f"independent scan stripping comment lines gives {want_code}")
    if "in_production_code" in em:
        assert em["in_production_code"] == want_prod, (
            f"the audit reports {em['in_production_code']} in production code; "
            f"an independent scan excluding #[cfg(test)] gives {want_prod}")
    # and the page must not claim universality when it is not universal
    if em["mention"] < em["modules"]:
        assert "every module of its simulation engine, concerns" not in md, (
            f"{len(silent)} modules mention neither ferroptosis nor a "
            "physical-ROS modality, and the page still claims every one does")
        # The renderer has three branches, and "mostly but not entirely" only
        # appears in the one where EVERY prediction is ferroptosis. Registering
        # P9-P13 moved the page to the third branch, so the guard was pinning
        # a sentence the page correctly stops printing. What must hold in any
        # branch is that the shortfall is stated as a ratio, not swallowed.
        assert (f"{em['mention']} of {em['modules']} engine modules" in md
                or "mostly but not entirely" in md), (
            "the page no longer reports how many modules mention ferroptosis "
            "against how many exist")
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
    # THE PREDICATE MATTERS, not just the conditionality. Branching on
    # `mention` meant three doc-comment edits could restore the universal
    # claim: a module citing a PMID in a comment does not "concern"
    # ferroptosis in the sense the sentence asserts.
    all_mention = {**d, "engine_modules": {**d["engine_modules"],
                                           "mention": d["engine_modules"]["modules"],
                                           "silent": []}}
    assert "and every module of its simulation engine" not in m.render(all_mention), (
        "with every module MENTIONING it but only some using it in production "
        "code, the page prints the universal claim again; the headline is "
        "derived from the wrong predicate")

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
    body = d.get("therapy_by_body")
    assert body is not None, (
        "the body-route therapy count is gone, so the asymmetry is asserted "
        "rather than measured")
    n_ther = len(d["analyses"]["therapy-subject"])
    assert f"admits **{len(body)}** analyses" in md, (
        "the report does not state the body-route count")
    assert ("neither rule measures subject" in md
            or "not established by either rule" in md), (
        "the report presents one of the two rules as giving the corrected "
        "number; neither measures subject")
    # AND the decomposition, which is what shows the upper figure is not a
    # bound: an earlier version published "between 1 and 29" as a range when
    # several of the 29 sit in the page's own ferroptosis column.
    overlap_f = sorted(set(body) & set(d["analyses"]["ferroptosis-or-physical"]))
    if overlap_f:
        assert f"**{len(overlap_f)} of them are in this page's own FERROPTOSIS column**" in md, (
            f"{len(overlap_f)} body-route matches are already classified as "
            "ferroptosis, so the body-route count cannot bound the therapy "
            "count, and the report does not say so")
        assert "not an upper bound" in md


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


def test_the_readme_denominators_are_pinned_to_the_artifact():
    """The front-door half of the finding was guarded by nothing.

    All four figures were hand-typed into README.md and a reviewer falsified
    every one of them with the suite green. Derived on one page, retyped on the
    other, is the exact shape this finding is about.

    THE DENOMINATOR HAS MOVED. The README no longer quotes shares of the
    4,830-record retrieval; it quotes shares of the census records carrying a
    discriminative MeSH descriptor. So this pins the census figures, and pins
    that the retrieval figure appears ONLY as the contrast it is now used for.
    """
    import json as _j
    land = REPO_ROOT / "analysis" / "atlas-landscape.json"
    r = README.read_text()
    if not land.exists():
        return
    rows = _j.loads(land.read_text())["rows"]
    tot = sum(x["mesh_census"] or 0 for x in rows)
    top = max(rows, key=lambda x: x["mesh_census"] or 0)
    im = next(x for x in rows if x["mechanism"].lower() == "immunotherapy")
    assert f"{tot:,}" in r, (
        f"the README does not carry the census mechanism denominator {tot:,}")
    assert f"{im['mesh_census']:,}" in r and f"{100*im['mesh_census']/tot:.1f}%" in r, (
        f"the README does not carry immunotherapy's derived census share "
        f"{im['mesh_census']:,} = {100*im['mesh_census']/tot:.1f}%")
    # the leading mechanism by volume is a SCOPE ARTIFACT if one descriptor
    # carries most of it; the README must say so rather than rank it silently
    if (top.get("top_share") or 0) >= 0.5:
        assert top["mechanism"] in r and "scope artifact" in r.lower(), (
            f"`{top['mechanism']}` leads the census ranking on "
            f"{100*top['top_share']:.0f}% of one descriptor and the README "
            "presents the ranking without saying so")
    # and the retrieval share may appear only as a contrast, never as the claim
    d = _doc()
    m = d.get("mechanism_denominators")
    if m:
        old_share = f"{100*m['share_of_corpus']:.1f}%"
        for i in range(len(r)):
            j = r.find(old_share, i)
            if j < 0:
                break
            w = r[max(0, j - 400):j + 200]
            assert "keyword retrieval" in w or "contrast" in w or "gap" in w, (
                f"the README quotes {old_share} without marking it as the "
                "retrieval figure it is being contrasted against")
            i = j + 1


def test_the_module_rows_are_each_pinned_to_their_own_figure():
    """`in_code` could be forced to equal `mention` undetected.

    The README guard asserted `f"**{em['in_code']} of {em['modules']}**" in r`,
    which the MENTION row already satisfies when the two are equal -- so a
    mutation making every module count as in-code passed while the page
    published a figure four higher than the truth.
    """
    d, r, md = _doc(), README.read_text(), MD.read_text()
    em = d["engine_modules"]
    if not isinstance(em, dict):
        return
    assert em["in_code"] <= em["mention"] <= em["modules"]
    assert em.get("in_production_code", 0) <= em["in_code"], (
        "more modules mention it in production code than in code at all")
    # each row must appear with its OWN label, so one cannot stand in for another
    for label, key in (("mentioning it anywhere", "mention"),
                       ("mentioning it in code", "in_code"),
                       ("mentioning it in PRODUCTION code", "in_production_code")):
        if key not in em:
            continue
        assert f"{label} | **{em[key]} of {em['modules']}**" in md, (
            f"the `{label}` row does not carry its own measured value "
            f"{em[key]}")
        # AND on the README, which is the front door and was NOT checked here:
        # `r` was read and never used for these rows, so #799 shipped a
        # PRODUCTION-code figure of 16 against a generator saying 15 -- the
        # value came from a sibling artifact measuring a different term set --
        # and this file stayed green.
        assert f"{label} | **{em[key]} of {em['modules']}**" in r, (
            f"README's `{label}` row is not the measured value {em[key]} of "
            f"{em['modules']}; the front door disagrees with the generator")


def test_the_denominator_verdict_is_checked_against_a_live_count():
    """The guard built `counts` from INDEX.jsonl and never compared it.

    A reviewer made the generator report the SECOND-largest mechanism as top
    -- the page published "nanoparticle is 515 of 4,830 = 10.7%" beside a
    README saying 47.6% -- and the suite stayed green, because the local was
    dead. Three of six fields were free.
    """
    d = _doc()
    m = d.get("mechanism_denominators")
    if not m:
        return
    counts = {}
    with (REPO_ROOT / "corpus" / "INDEX.jsonl").open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            for x in (json.loads(line).get("mechanisms") or []):
                counts[x] = counts.get(x, 0) + 1
    top, top_n = max(counts.items(), key=lambda kv: kv[1])
    assert m["top_mechanism"] == top, (
        f"the artifact reports `{m['top_mechanism']}` as the top mechanism; "
        f"a live count of INDEX.jsonl gives `{top}`")
    assert m["top_n"] == top_n
    assert abs(m["share_of_corpus"] - top_n / m["corpus_records"]) < 1e-12
    assert abs(m["share_of_tagged"] - top_n / m["tagged"]) < 1e-12
    assert abs(m["share_of_tags"] - top_n / m["tags"]) < 1e-12
