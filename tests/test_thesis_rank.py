"""Guards for the all-modality thesis ranking (#725).

THE CLAIM
---------
Ranking the MeSH `Ferroptosis` intersection against every modality descriptor in
the committed partition universe puts the thesis's own central mechanism --
sonodynamic therapy -- at rank 22 of 166, on 32 articles. The largest
intersection is `antineoplastic agents` at 954, which is 29.8x the sonodynamic
leg, and `immunotherapy` is 370.

So the unexamined precedent is chemotherapy, not radiotherapy.

WHY THIS NEEDS GUARDING
-----------------------
1. THE UNIVERSE MUST NOT BE HAND-PICKED. The point is that a hand-written list
   cannot surface a leg nobody thought of. If this script ever builds its own
   modality list, it reproduces the defect it exists to correct.

2. THE COMPARISON MUST BE DESCRIPTOR-TO-DESCRIPTOR. The withdrawn version of
   this issue claimed radiotherapy was ~10x the sonodynamic leg, from a wide
   text-stem count against a narrow descriptor count. That reproduced under no
   definition. Both sides must come from the same axis.

3. THE WHOLE RANKING MUST BE PUBLISHED. A single pair quoted from a ranking
   nobody can see is a conclusion chosen first and evidenced afterwards.

4. UNMEASURABLE MODALITIES MUST BE NAMED, NOT ZEROED. A modality with no usable
   descriptor reported as 0 is a different and false claim.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_thesis_rank.py"
MD = REPO_ROOT / "analysis" / "atlas-thesis-rank.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-thesis-rank.json"
PARTITIONS = REPO_ROOT / "analysis" / "modality-partitions.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("atr", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_universe_comes_from_the_committed_partitions():
    """A self-built list reproduces the defect this exists to correct."""
    src = SCRIPT.read_text()
    assert "d = json.loads(PARTITIONS.read_text())" in src, (
        "the universe is no longer LOADED from the partition file; naming the "
        "path in a docstring while building the list inline is the same defect "
        "wearing the right words")
    assert "modality-partitions.json" in src, (
        "the modality universe is no longer the committed partition file, so "
        "this ranks against a list written by the same person making the claim")
    assert "def modality_universe" in src
    d = _doc()
    parts = json.loads(PARTITIONS.read_text())
    expected = {x.strip().lower() for spec in parts.values()
                for side in ("pharmacological", "physical")
                for x in spec.get(side, [])}
    assert d["universe_size"] == len(expected), (
        f"the audit used {d['universe_size']} descriptors, the partition file "
        f"holds {len(expected)}")


def test_the_comparison_is_descriptor_to_descriptor():
    """The withdrawn 10x claim came from mixing a text stem with a descriptor."""
    src = SCRIPT.read_text()
    # everything counted must come from the record's `mesh` field
    assert 'r.get("mesh")' in src, "the scan no longer reads the descriptor axis"
    assert "title" not in src.split("def scan(")[1].split("def render")[0], (
        "the scan reads titles; mixing a text count with a descriptor count is "
        "exactly the asymmetry that made the withdrawn 10x claim irreproducible")


def test_the_whole_ranking_is_published():
    """A pair quoted from an unpublished ranking is a chosen conclusion."""
    d, md = _doc(), MD.read_text()
    assert len(d["intersections"]) > 50, (
        f"only {len(d['intersections'])} intersections recorded; the ranking is "
        "too short to locate a leg within it")
    assert "| rank |" in md, "the ranking table is gone"
    # the leg table must report a rank OUT OF the full ranking, not in isolation
    assert f"rank of {len(d['intersections'])}" in md, (
        "the leg table does not state the size of the ranking its ranks are "
        "drawn from")


def test_the_thesis_legs_are_located_not_asserted():
    """Their rank must come from the computed ranking."""
    d, md = _doc(), MD.read_text()
    counts = dict(d["intersections"])
    rank = {k: i + 1 for i, (k, _v) in enumerate(d["intersections"])}
    sdt = "ultrasonic therapy"
    assert sdt in counts, "the sonodynamic descriptor is not in the universe"
    assert f"| sonodynamic therapy | {sdt} | {counts[sdt]:,} | {rank[sdt]} |" in md, (
        "the sonodynamic row is not the one the ranking supports")
    # the headline multiple must be derived from the same two numbers
    top_name, top_n = d["intersections"][0]
    assert f"{top_n/counts[sdt]:.1f}x" in md, (
        "the headline multiple is not top-count over sonodynamic count")


def test_the_top_intersection_is_not_the_one_the_thesis_names():
    """If it were, this analysis would have found nothing worth reporting."""
    d = _doc()
    top = d["intersections"][0][0]
    named = {"ultrasonic therapy", "photochemotherapy", "photosensitizing agents",
             "high-intensity focused ultrasound ablation"}
    assert top not in named, (
        f"the largest intersection ({top}) is one the thesis already names, so "
        "the finding as written -- that something unexamined is larger -- does "
        "not hold and the report needs rewriting")


def test_unmeasurable_modalities_are_named_not_zeroed():
    src, md = SCRIPT.read_text(), MD.read_text()
    assert "NOT_MEASURABLE" in src and "named rather than shown as zero" in md, (
        "modalities with no usable descriptor are no longer named, so their "
        "absence reads as a measured zero")
    # THIS GUARD USED TO PIN THE WORD "OVER-estimate", i.e. it pinned a claim
    # that turned out to be BACKWARDS, so the only way to fail it was to keep
    # the error. `Ultrasonic Therapy` IS broader, but measured that is 3
    # records of 32 (precision 90.6%) against a recall of 46.0% -- the count
    # is an UNDER-estimate by roughly twofold. Pin the DIRECTION the sibling
    # measurement supports, and re-derive it rather than trusting either
    # document's wording.
    import json as _json
    rec_path = REPO_ROOT / "analysis" / "atlas-descriptor-recall.json"
    assert rec_path.exists(), (
        "the descriptor-recall measurement this caveat depends on is gone; "
        "the direction of the sonodynamic caveat is unsupported without it")
    rec = _json.loads(rec_path.read_text())
    sdt = rec["arms"]["SDT"]
    breadth = sdt["descriptor"] - sdt["both"]      # over-count: wrong subject
    shortfall = sdt["text"] - sdt["both"]          # under-count: missed papers
    if shortfall > breadth:
        assert "UNDER-estimate" in md and "UPPER bound" in md, (
            f"recall misses {shortfall} sonodynamic papers against a breadth "
            f"of {breadth}, so the count is an under-estimate and ratios "
            "against it are upper bounds; the report does not say so")
        for m in re.finditer(r"OVER-estimate", md):
            w = md[max(0, m.start() - 300):m.end() + 300]
            assert re.search(r"UNDER|earlier version|stated the opposite", w), (
                "the report still calls the count an over-estimate without "
                "the correction beside it")
    else:
        assert "OVER-estimate" in md, (
            "breadth now exceeds the recall shortfall and the report does not "
            "say the count is an over-estimate")


def test_attention_is_not_reported_as_endorsement():
    md = MD.read_text()
    assert "ATTENTION, not endorsement" in md, (
        "the report no longer distinguishes a co-indexed literature from "
        "evidence that a combination works")
    assert "ATTENTION, not endorsement" in SCRIPT.read_text(), (
        "the renderer no longer emits that distinction")


def test_an_empty_scan_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'if d["ferroptosis_total"] == 0:' in src
    assert "is not a finding" in src and "raise SystemExit" in src


def test_absent_from_the_universe_renders_differently_from_zero():
    """A scope boundary is not a measured zero.

    `counts.get(desc, 0)` printed 0 under a column headed "count" for
    `drug resistance, neoplasm`, whose real intersection is published by the
    sibling atlas-thesis-position artifact from the same build and called
    "the strong one". Absence and zero rendered identically.
    """
    d, md = _doc(), MD.read_text()
    legs = d.get("leg_intersections") or {}
    assert legs, "the thesis legs are no longer counted independently"
    ranked = {k for k, _v in d["intersections"]}
    for leg, desc in _mod().THESIS_LEGS.items():
        c = legs.get(desc)
        assert c is not None, f"{desc} has no measured count"
        if desc not in ranked:
            assert f"| {leg} | {desc} | {c:,} |" in md, (
                f"{leg} is outside the ranking universe and the report does "
                f"not print its real count of {c:,}")
            assert "not a modality; outside this universe" in md
            assert f"| {leg} | {desc} | 0 |" not in md, (
                "an out-of-universe leg is rendering as a measured zero again")
            # the EXPLANATORY paragraph, not just the table cell: suppressing
            # it left the reader with a rank cell they cannot interpret, and
            # the suite green
            assert "outside the modality universe by construction" in md, (
                f"{leg} sits outside the universe and the report no longer "
                "explains that the boundary is by construction rather than a "
                "measured absence")
            assert leg in md.split("outside the modality universe")[0][-400:], (
                f"the scope-boundary paragraph does not name {leg}")


def test_the_out_of_universe_count_agrees_with_the_sibling_artifact():
    """Two documents, same build, must not disagree about the same number."""
    d = _doc()
    sib = REPO_ROOT / "analysis" / "atlas-thesis-position.json"
    if not sib.exists():
        return
    other = json.loads(sib.read_text())
    legs = d.get("leg_intersections") or {}
    # COMPARE THE FIELD, not a substring of the whole blob. `str(c) in
    # json.dumps(other)` is a scan in which roughly a third of small integers
    # match something, so changing 479 to 478 passed.
    def _find(obj, want_desc):
        """Every numeric value stored against a key naming this descriptor."""
        out = []
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(k, str) and want_desc in k.lower() and \
                        isinstance(v, (int, float)):
                    out.append(v)
                out += _find(v, want_desc)
        elif isinstance(obj, list):
            # rows like ["drug resistance", 479, ...]
            if any(isinstance(x, str) and want_desc in x.lower() for x in obj):
                out += [x for x in obj if isinstance(x, (int, float))]
            for v in obj:
                out += _find(v, want_desc)
        return out

    checked = 0
    for desc, c in legs.items():
        if "drug resistance" not in desc or not c:
            continue
        vals = _find(other, "drug resistance")
        assert vals, (
            "the sibling atlas-thesis-position artifact stores no number "
            "against a drug-resistance key, so the two cannot be reconciled")
        assert c in vals, (
            f"this report counts {c:,} for `{desc}`; the sibling artifact "
            f"stores {sorted(set(vals))} against its drug-resistance keys. "
            "One of the two is from a different build")
        checked += 1
    assert checked, "the cross-artifact check exercised nothing"


def test_the_headline_ratio_carries_its_prevalence_normalisation():
    """A ratio of an umbrella descriptor to a specific one measures breadth.

    Normalising each count by how common its descriptor is census-wide
    REVERSES the direction, so publishing the raw ratio alone is publishing
    the more flattering of two answers.
    """
    d, md = _doc(), MD.read_text()
    tot = d.get("census_descriptor_totals") or {}
    rows = d["intersections"]
    if not rows:
        return
    top_name, top_c = rows[0]
    sdt_desc = _mod().THESIS_LEGS["sonodynamic therapy"]
    sdt_c = dict(rows).get(sdt_desc)
    if not (sdt_c and tot.get(top_name) and tot.get(sdt_desc)):
        return
    base = d["ferroptosis_total"] / d["census"]
    e_top = (top_c / tot[top_name]) / base
    e_sdt = (sdt_c / tot[sdt_desc]) / base
    assert f"**{e_top:.2f}x**" in md and f"**{e_sdt:.2f}x**" in md, (
        "the report does not state both normalised enrichments its own "
        "artifact implies")
    # and the verdict must follow the numbers, not a fixed sentence
    if e_sdt > e_top:
        assert "The direction reverses" in md, (
            f"normalised, the sonodynamic descriptor is {e_sdt:.2f}x its base "
            f"rate against the umbrella term's {e_top:.2f}x, and the report "
            "does not say the direction reverses")
    else:
        assert "The direction reverses" not in md


def test_a_no_descriptor_claim_is_checked_against_the_census():
    """A hand-written 'no descriptor' list beside a census that can answer it.

    "cold atmospheric plasma (no descriptor)" shipped while `Plasma Gases`
    carried 474 census records and 9 ferroptosis intersections -- more than a
    leg this report ranks.
    """
    m, d, md = _mod(), _doc(), MD.read_text()
    cand = d.get("candidate_intersections")
    assert cand is not None, "the no-descriptor claims are no longer checked"
    assert m.CANDIDATE_DESCRIPTORS, "no candidate descriptors are declared"
    worst = min((v for _k, v in d["intersections"]), default=0)

    # RECOUNT the candidates over a stride, so a zeroed counter cannot make
    # the report claim "no descriptor" by construction. Zeroing `cand_inter`
    # sent every candidate to 0, the `n > worst` branch never ran, and the
    # whole refutation vanished with the suite green -- the guard only
    # checked the strong case. A stride is a SUBSET, so the full count must
    # be at least the strided one.
    import gzip
    # OFFLINE CONTRACT: corpus/atlas/records/ is gitignored bulk data, so CI
    # has no shards. The recount can only run where the census exists; it
    # still fires for anyone holding the data, which is where a scan-level
    # change would be made.
    import pytest
    if not any((m.ATLAS / "records").glob("*.jsonl.gz")):
        pytest.skip("census shards not present in this checkout")

    want = {d2.lower() for d2 in m.CANDIDATE_DESCRIPTORS.values()}
    seen = {k: 0 for k in want}
    for f in sorted((m.ATLAS / "records").glob("*.jsonl.gz"))[::12]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                mesh = {x.lower() for x in (r.get("mesh") or [])}
                if "ferroptosis" not in mesh:
                    continue
                for k in mesh & want:
                    seen[k] += 1
    for k, lo in seen.items():
        assert cand.get(k, 0) >= lo, (
            f"the artifact reports {cand.get(k, 0)} ferroptosis intersections "
            f"for `{k}` and a 1-in-12 stride already finds {lo}; the counter "
            "is not counting")
    assert sum(seen.values()) > 0, (
        "no candidate descriptor intersects ferroptosis even on a stride, so "
        "this guard is not exercising the refutation path")

    for modality, desc in m.CANDIDATE_DESCRIPTORS.items():
        n = cand.get(desc, 0)
        assert f"| {modality} | `{desc}` |" in md, (
            f"{modality}'s candidate descriptor is not shown with its counts")
        if n > worst:
            assert "measurable after all" in md, (
                f"`{desc}` carries {n:,} ferroptosis intersections against a "
                f"smallest ranked entry of {worst:,}, and the report still "
                "presents it as having no usable descriptor")
            assert f"{modality} (no descriptor)" not in md, (
                f"{modality} is still listed as having no descriptor while "
                "one measurably exists")


def test_the_normalisation_denominators_are_recounted_not_trusted():
    """`census_descriptor_totals` is the whole denominator of the new headline.

    Nothing recounted it. A reviewer multiplied the sonodynamic census total by
    ten and the page printed "The direction holds under normalisation" -- the
    PR's entire finding reversed, suite green. Recounted over a 1-in-12 STRIDE
    (a subset, so the full count must be at least the strided one) and checked
    for internal impossibility.
    """
    import gzip
    m, d = _mod(), _doc()
    # OFFLINE CONTRACT: corpus/atlas/records/ is gitignored bulk data, so CI
    # has no shards. The recount can only run where the census exists; it
    # still fires for anyone holding the data, which is where a scan-level
    # change would be made.
    import pytest
    if not any((m.ATLAS / "records").glob("*.jsonl.gz")):
        pytest.skip("census shards not present in this checkout")

    tot = d.get("census_descriptor_totals") or {}
    assert tot, "the census denominators are gone"

    # INTERNAL IMPOSSIBILITY: an intersection cannot exceed the census total.
    # Zeroing the census counter published "0 census records, 9 ferroptosis
    # intersections" with everything green.
    inter = dict(d["intersections"])
    inter.update(d.get("leg_intersections") or {})
    inter.update(d.get("candidate_intersections") or {})
    for k, v in inter.items():
        if k in tot:
            assert v <= tot[k], (
                f"`{k}` has {v:,} ferroptosis intersections against a census "
                f"total of {tot[k]:,}, which is impossible")

    watch = {d["intersections"][0][0],
             m.THESIS_LEGS["sonodynamic therapy"]} | \
            {x.lower() for x in m.CANDIDATE_DESCRIPTORS.values()}
    seen = {k: 0 for k in watch}
    for f in sorted((m.ATLAS / "records").glob("*.jsonl.gz"))[::12]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                mesh = {x.lower() for x in (r.get("mesh") or [])}
                for k in mesh & watch:
                    seen[k] += 1
    for k, lo in seen.items():
        assert tot.get(k, 0) >= lo, (
            f"the artifact reports {tot.get(k, 0):,} census records for `{k}` "
            f"and a 1-in-12 stride already finds {lo:,}; the denominator is "
            "not counting what it claims")
        # and not absurdly larger than a 12x extrapolation
        assert tot.get(k, 0) <= max(lo * 40, 500), (
            f"the artifact reports {tot.get(k, 0):,} census records for `{k}` "
            f"against a strided estimate near {12*lo:,}; the denominator has "
            "been inflated")


def test_the_measurable_and_not_measurable_verdicts_are_mutually_exclusive():
    """Loosening `>` to `>=` published both at once, in adjacent paragraphs."""
    md = MD.read_text()
    if "measurable after all" not in md:
        return
    # Bound each region at the paragraph break. A fixed character window
    # spilled from the refuting paragraph into the ties paragraph and flagged
    # a descriptor that appears in only one of them -- the guard was reading
    # its own overrun.
    after = md.split("measurable after all", 1)[1]
    refute = after.split("\n\n", 1)[0]
    named_refuting = [x for x in _mod().CANDIDATE_DESCRIPTORS.values()
                      if x in refute]
    assert named_refuting, (
        "the refuting paragraph names no candidate descriptor, so the claim "
        "cannot be checked against one")
    # The ties paragraph must EXIST when there are ties: deleting it removed
    # the statement that the not-measurable framing still stands for the
    # candidate that does not refute, with the suite green -- the same hole as
    # the scope-boundary paragraph one section over.
    m2, d2 = _mod(), _doc()
    cand2 = d2.get("candidate_intersections") or {}
    rows2 = d2["intersections"]
    cnames = {x.lower() for x in m2.CANDIDATE_DESCRIPTORS.values()}
    worst2 = min([v for k, v in rows2 if k not in cnames], default=0)
    ties = [x for x in m2.CANDIDATE_DESCRIPTORS.values()
            if 0 < cand2.get(x.lower(), 0) <= worst2]
    if ties:
        assert "the not-measurable framing stands" in md, (
            f"{ties} do not exceed the benchmark, and the report no longer "
            "says the not-measurable framing stands for them")
        for x in ties:
            assert x in md
    if "the not-measurable framing stands" in md:
        para = md.split("the not-measurable framing stands", 1)[0]
        stands = para.rsplit("\n\n", 1)[-1]
        for desc in named_refuting:
            assert desc not in stands, (
                f"`{desc}` is reported BOTH as measurable after all and as one "
                "the not-measurable framing stands for, in adjacent "
                "paragraphs")


def test_the_candidate_benchmark_excludes_the_candidates_themselves():
    """Some candidates ARE in the ranked universe, so min(rows) includes them."""
    m, d = _mod(), _doc()
    cand = {x.lower() for x in m.CANDIDATE_DESCRIPTORS.values()}
    rows = d["intersections"]
    inside = [k for k, _v in rows if k in cand]
    if not inside:
        return
    others = [v for k, v in rows if k not in cand]
    worst = min(others, default=0)
    md = MD.read_text()
    assert f"smallest RANKED entry of {worst:,}" in md, (
        f"{inside} are themselves in the ranking, so the benchmark must "
        f"exclude them; it should be {worst:,}")
