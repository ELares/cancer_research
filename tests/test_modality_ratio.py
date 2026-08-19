"""Guards for the modality-class sensitivity measurement (#724).

THE CLAIM
---------
The manuscript's pharmacological:physical ratio is reported as 9.1:1 by its own
method and 17.6:1 on the census. Five partitions, each moving BOTH classes, give
1.32:1 to 3.93:1 -- and that spread is one undeclared variable: how much
operative surgery each partition admits to the physical class. Held out, the
five agree. The direction survives; the order-of-magnitude does not.

WHAT MAKES THIS EASY TO FAKE, AND THEREFORE WORTH GUARDING
-----------------------------------------------------------
1. ONE-SIDED WIDENING. Adding radiotherapy to PHYSICAL while leaving cytotoxic
   chemotherapy out of PHARMACOLOGICAL produces an inversion by construction.

2. THE INPUT CONTAINING THE RESULT. The partitions are a committed input file.
   If it also carried each panel's ratio the measurement could validate itself
   against a number it was handed.

3. THE HOLD-OUT RULE BECOMING THE NEW UNDECLARED VARIABLE. This is the one a
   review found live. Every guard below that touches the surgical rule was
   written after mutations showed the previous set green while the rule
   measured RADIOTHERAPY, while it carried a per-partition exception, and while
   the comparator was narrowed on one side only.

4. A GUARD THAT READS ONLY THE COMMITTED ARTIFACT. The .md and .json go stale
   TOGETHER, so comparing them cannot fail. The freshness gate re-renders and
   diffs, and the semantic guards recompute from an INDEPENDENT word list
   rather than from the generator's own regex.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_modality_ratio.py"
PARTITIONS = REPO_ROOT / "analysis" / "modality-partitions.json"
MD = REPO_ROOT / "analysis" / "atlas-modality-ratio.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-modality-ratio.json"
RECORDS = REPO_ROOT / "corpus" / "atlas" / "records"

# INDEPENDENT of the generator's regex, on purpose. A guard that imports
# `SURGICAL` to decide whether the held-out descriptors are surgical compares
# the rule to itself and passes for any rule at all.
_SURGICAL_WORDS = re.compile(
    r"surg|ectom|otom|ostom|resect|excis|scop|transplant|amputat|dissect|"
    r"anastomos|castration|conization|curettage|enucleation|exenteration|"
    r"reoperation|salvage|orthopedic|vertebroplasty|operative|"
    r"metastasectomy|craniotomy|laparotomy|thoracotomy|sternotomy", re.I)
_RADIOTHERAPY_WORDS = re.compile(
    r"radiotherap|irradiat|brachytherap|proton therapy|x-ray therapy|"
    r"radioisotope|radiation|boron neutron|radioimmunotherapy", re.I)


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("mr", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _parts():
    return json.loads(PARTITIONS.read_text())


# ---------------------------------------------------------------------------
# The input, and the partitions themselves
# ---------------------------------------------------------------------------

def test_the_input_does_not_carry_the_result():
    """A partition file holding its own ratio lets the measurement self-validate.

    Exercised BEHAVIOURALLY. The previous version asserted that the source text
    contained the condition, which `if False and (...)` satisfies while doing
    nothing -- a planted mutation shipped green that way.
    """
    m = _mod()
    for name, spec in _parts().items():
        assert "ratio" not in spec and "reported_ratio" not in spec, (
            f"{name} carries a ratio in the input file")
    good = {"p": {"pharmacological": ["a"], "physical": ["b"]}}
    assert m.load_partitions(good) == good
    for poison in ("ratio", "reported_ratio"):
        doctored = {"p": {"pharmacological": ["a"], "physical": ["b"],
                          poison: 4.2}}
        with pytest.raises(SystemExit):
            m.load_partitions(doctored)
    with pytest.raises(SystemExit):
        m.load_partitions({"p": {"pharmacological": ["a"], "physical": []}})
    with pytest.raises(SystemExit):
        m.load_partitions({})


def test_every_partition_moves_both_classes():
    """One-sided widening manufactures whatever inversion it finds."""
    for name, spec in _parts().items():
        phys = " | ".join(spec["physical"]).lower()
        pharm = " | ".join(spec["pharmacological"]).lower()
        assert "radio" in phys or "irradiat" in phys or "brachytherapy" in phys, (
            f"{name}: the physical class has no radiotherapy-family member, so "
            "it repeats the omission this analysis exists to correct")
        assert ("chemotherapy" in pharm or "antineoplastic" in pharm
                or "drug therapy" in pharm), (
            f"{name}: the pharmacological class has no cytotoxic-chemotherapy "
            "member while the physical class gained radiotherapy -- that is "
            "the one-sided widening that manufactures an inversion")


def test_every_member_is_a_descriptor_not_a_heading():
    """A leaked section header inflates a member count and a held-out set.

    `minimally-inclusive` carried `MECHANICAL REMOVAL, i.e. SURGERY (62):`,
    which matched the surgical rule and was counted as a physical member.
    """
    for name, spec in _parts().items():
        for cls in ("pharmacological", "physical"):
            for x in spec[cls]:
                assert not x.endswith(":"), f"{name}/{cls}: {x!r} is a heading"
                assert not re.search(r"\(\d+\)\s*:?$", x), (
                    f"{name}/{cls}: {x!r} carries a count, so it is a heading")
                # MeSH descriptors may start with a digit (`3-Iodobenzylguanidine`)
                # but never carry a run of shouted words.
                assert not re.search(r"\b[A-Z]{3,}\b\s+\b[A-Z]{3,}\b", x), (
                    f"{name}/{cls}: {x!r} reads as a section heading")
                assert x == x.strip(), f"{name}/{cls}: {x!r} is not trimmed"


def test_several_partitions_and_they_genuinely_differ():
    """The published spread exists; near-identical partitions do not give one."""
    d = _doc()["partitions"]
    assert len(d) >= 4, f"only {len(d)} partitions"
    ratios = sorted(v["ratio"] for v in d.values())
    assert ratios[-1] / max(ratios[0], 1e-9) > 1.5, (
        f"all partitions land within {ratios[-1]/ratios[0]:.2f}x of each other, "
        "so this does not demonstrate sensitivity to the class boundary")
    parts = _parts()
    sizes = {k: (len(v["pharmacological"]), len(v["physical"]))
             for k, v in parts.items()}
    assert len(set(sizes.values())) >= 3, (
        f"partition sizes are near-identical {sizes}")


# ---------------------------------------------------------------------------
# The hold-out rule: the variable this page introduced
# ---------------------------------------------------------------------------

def test_the_holdout_rule_cannot_depend_on_which_partition_it_is_applied_to():
    """A per-partition exception reproduces, one level up, the defect measured.

    The previous guard asserted `src.count("SURGICAL.search") == 1`. A planted
    mutation added `and not (k == "regulatory-class" and "transplant" in x)`
    inside the same comprehension -- still one call, still green, and a
    per-partition rule.

    Tested as a PROPERTY instead: a rule must be determined by the descriptor,
    so applying it to the union of every member set and intersecting must give
    the same answer as applying it to each set alone.
    """
    m = _mod()
    parts = _parts()
    sets = {k: {x.lower() for x in v["physical"]} for k, v in parts.items()}
    union = set().union(*sets.values())
    rules = dict(m.holdout_rules())
    rules.update(m.CONTROL_RULES)
    assert len(rules) >= 5, "the rule panel has shrunk to a single judgement"
    for name, rule in rules.items():
        on_union = rule(union)
        for k, s in sets.items():
            assert rule(s) == (on_union & s), (
                f"rule {name!r} gives a different answer for {k} than for the "
                "same descriptors seen in the union, so it is a per-partition "
                "judgement rather than one rule")
        # and permuting the input must not change it either
        shuffled = {x for x in sorted(union, reverse=True)}
        assert rule(shuffled) == on_union, (
            f"rule {name!r} depends on iteration order")
        # NOR may it change its mind with the SIZE of the set it is handed.
        # A planted `len(members) > 130` exception survived the union check
        # above, because the descriptors it excluded live only in the two
        # large partitions -- so the check passed vacuously on the other three.
        for pad in (0, 40, 200):
            padded = union | {f"zz filler descriptor {i}" for i in range(pad)}
            assert rule(padded) & union == on_union, (
                f"rule {name!r} gives a different verdict on the same "
                f"descriptors once {pad} unrelated members are added, so it "
                "is size-dependent -- a per-partition judgement in disguise")
        # The strongest form, and the one a `len(members) > 130` exception
        # slips past every weaker check: the verdict on a descriptor must be
        # the same when it is the ONLY member. Padding the union cannot see
        # this, because the union is larger than every partition.
        for x in sorted(union):
            assert (x in rule({x})) == (x in on_union), (
                f"rule {name!r} treats {x!r} differently alone than in "
                "company, so it is not a function of the descriptor")


def test_the_holdout_rules_actually_select_surgery():
    """A rule that measured radiotherapy would leave every number intact.

    A planted mutation repointed the regex at the radiotherapy family and the
    page still read "whether operative surgery is a physical modality ... the
    spread collapses", with every guard green. Judged here against a word list
    written in this file, never imported from the generator.
    """
    d = _doc()
    ho = d.get("holdouts") or {}
    assert ho, "the hold-out rules are gone from the artifact"
    for name, leg in ho.items():
        picked = sorted({x for r in leg["partitions"].values()
                         for x in r["descriptors"]})
        assert picked, f"rule {name!r} selects nothing anywhere"
        surgical = [x for x in picked if _SURGICAL_WORDS.search(x)]
        radio = [x for x in picked if _RADIOTHERAPY_WORDS.search(x)
                 and not _SURGICAL_WORDS.search(x)]
        assert len(surgical) / len(picked) >= 0.9, (
            f"rule {name!r} holds out {len(picked)-len(surgical)} of "
            f"{len(picked)} descriptors that are not operative surgery by an "
            f"independent reading: {sorted(set(picked)-set(surgical))[:8]}")
        assert not radio, (
            f"rule {name!r} holds out radiotherapy descriptors {radio[:8]}; "
            "the page attributes the collapse to surgery")


# Descriptors that are unambiguously operative removal. Written here, checked
# against the partitions, and REQUIRED of every rule: the 90%-surgical floor
# below leaves ~10 free slots per rule, which a mutation filled with real
# energy descriptors while every guard stayed green.
_MUST_HOLD_OUT = {
    "surgical procedures, operative", "mastectomy", "prostatectomy",
    "hepatectomy", "pneumonectomy", "gastrectomy", "colectomy",
    "nephrectomy", "hysterectomy", "laparotomy", "thoracotomy",
    "lymph node excision", "metastasectomy",
}
_ENERGY_WORDS = re.compile(
    r"cryotherap|laser|photo|puva|ultraviolet|hyperthermi|diatherm|"
    r"radiofrequency|electropor|ablat|coagulation|irradiat|radiotherap|"
    r"brachytherap|proton therapy|light therapy", re.I)


def test_every_surgical_descriptor_a_rule_skips_is_a_NAMED_exception():
    """The page's own method, turned on the page: name every exclusion.

    Three mutations passed the previous guards -- one dropping the nine
    largest surgical descriptors NOT on the core list (headline 1.18x ->
    1.53x), one stopping the primary rule excluding the energy leak it
    describes, and one absorbing seven modalities that escape both hand lists.
    A 90%-surgical floor cannot see any of them, and REMOVING surgical terms
    even raises the surgical fraction. So the requirement is containment
    against the generator's own declared exception sets, not a percentage.
    """
    m, d = _mod(), _doc()
    parts = _parts()
    # the generator's OWN declared exception sets, imported not restated
    named = m.ENERGY_CAUGHT_BY_STEM | m.OPERATIVE_MISSED_BY_STEM
    for name, leg in d["holdouts"].items():
        for k, r in leg["partitions"].items():
            members = {x.lower() for x in parts[k]["physical"]}
            picked = set(r["descriptors"])
            skipped = {x for x in members
                       if _SURGICAL_WORDS.search(x) and x not in picked}
            unnamed = {x for x in skipped
                       if x not in named and "transplant" not in x}
            assert not unnamed, (
                f"rule {name!r} skips {sorted(unnamed)} in {k} -- surgical by "
                "an independent reading and not on any exception list the "
                "generator declares, so the page holds out less than it says")

    # THE EXCEPTION SETS THEMSELVES must pass the independent reading, or the
    # guard above is defeated by widening them: adding eight real operative
    # descriptors to `ENERGY_CAUGHT_BY_STEM` moved the headline collapse
    # 1.18x -> 1.54x with every guard green, while the page printed them under
    # "which are energy modalities rather than operative removal".
    for x in sorted(m.ENERGY_CAUGHT_BY_STEM):
        assert _ENERGY_WORDS.search(x) or re.search(
            r"cryosurg|radiosurg|electrosurg|ultrasonic surgical", x, re.I), (
            f"{x!r} is declared an energy modality the stem list wrongly "
            "catches, and reads as operative removal instead")
    for x in sorted(m.OPERATIVE_MISSED_BY_STEM):
        assert _SURGICAL_WORDS.search(x), (
            f"{x!r} is declared operative removal the stem list misses, and "
            "does not read as surgery by an independent word list")
        assert not _ENERGY_WORDS.search(x), (
            f"{x!r} is declared operative and reads as an energy modality")

    # and the PRIMARY rule's two declared corrections, pinned directly: no
    # guard anywhere used to reference what it excludes, only its spread
    rules = m.holdout_rules()
    for k, spec in parts.items():
        members = {x.lower() for x in spec["physical"]}
        got = rules[m.PRIMARY_RULE](members)
        assert not (m.ENERGY_CAUGHT_BY_STEM & got), (
            f"the primary rule holds out {sorted(m.ENERGY_CAUGHT_BY_STEM & got)} "
            f"in {k}, while the page says it excludes them and calls "
            "`Radiosurgery` radiotherapy mass")
        assert (m.OPERATIVE_MISSED_BY_STEM & members) <= got, (
            f"the primary rule does not hold out "
            f"{sorted((m.OPERATIVE_MISSED_BY_STEM & members) - got)} in {k}, "
            "which the page says it includes")


def test_the_zero_surgery_paragraph_is_present_while_a_zero_partition_exists():
    """It vanished silently under a mutation, taking the objection it answers."""
    d, md = _doc(), MD.read_text()
    zero = [k for k, c in d["partitions"].items()
            if c["surgical_share_of_physical"] == 0]
    wz = d.get("without_zero_surgery_partitions") or {}
    assert sorted(zero) == sorted(wz.get("excluded", [])), (
        f"partitions with no surgical mass are {sorted(zero)} but the "
        f"leave-one-out block excludes {sorted(wz.get('excluded', []))}")
    if zero:
        assert wz.get("held_out_spread"), "the leave-one-out was not computed"
        assert "admits no operative surgery at all" in md, (
            "a partition the primary rule cannot move is the top of both "
            "columns and the page does not say so")
        assert f"**{wz['held_out_spread']:.2f}x**" in md


def test_every_rule_holds_out_the_core_and_no_energy_modality():
    """A rule that drops nine real surgical terms, or absorbs nine energy ones,
    moved the headline spread by 60% and 12% respectively and passed every
    guard. Both are pinned against lists written in this file.
    """
    d = _doc()
    parts = _parts()
    for name, leg in d["holdouts"].items():
        for k, r in leg["partitions"].items():
            members = {x.lower() for x in parts[k]["physical"]}
            picked = set(r["descriptors"])
            core = _MUST_HOLD_OUT & members
            missing = core - picked
            assert not missing, (
                f"rule {name!r} does not hold out {sorted(missing)} in {k}, "
                "which are operative removal by any reading")
            energy = {x for x in picked
                      if _ENERGY_WORDS.search(x) and not _SURGICAL_WORDS.search(x)}
            assert not energy, (
                f"rule {name!r} holds out {sorted(energy)} in {k}; those "
                "deliver energy rather than removing tissue, and the page "
                "attributes the collapse to surgery")


def test_the_rules_disagree_with_each_other():
    """Five rules that pick identical sets are one rule wearing five names."""
    d = _doc()["holdouts"]
    sets = {}
    for name, leg in d.items():
        sets[name] = frozenset(
            (k, x) for k, r in leg["partitions"].items() for x in r["descriptors"])
    assert len(set(sets.values())) >= 4, (
        "the hold-out rules pick the same descriptors, so the page's claim "
        "that five disagreeing rules agree on the outcome is vacuous")


def test_the_leak_lists_are_real_and_rendered():
    """The rule's failures are enumerated, not asserted."""
    m, md = _mod(), MD.read_text()
    parts = _parts()
    universe = {x.lower() for v in parts.values() for x in v["physical"]}
    stem = {x for x in universe if m.SURGICAL.search(x)}
    assert m.ENERGY_CAUGHT_BY_STEM <= stem, (
        "the named energy leak is no longer caught by the stem list, so the "
        "page describes a leak that does not exist")
    assert m.OPERATIVE_MISSED_BY_STEM <= (universe - stem), (
        "descriptors listed as missed by the stem list are in fact caught")
    for x in sorted(m.ENERGY_CAUGHT_BY_STEM | m.OPERATIVE_MISSED_BY_STEM):
        assert x in md.lower(), f"{x!r} is measured as a leak and not rendered"
    inert = [k for k, v in (_doc().get("stem_alternative_hits") or {}).items()
             if not v]
    for alt in inert:
        assert re.search(rf"`{re.escape(alt)}`", md), (
            f"stem {alt!r} matches nothing and the page does not say so")
        assert not any(alt in x for x in
                       {y for v in parts.values() for y in
                        map(str.lower, v["physical"])}), (
            f"stem {alt!r} is reported inert but matches a member")


# ---------------------------------------------------------------------------
# The arithmetic, and the prose that describes it
# ---------------------------------------------------------------------------

def test_the_arithmetic_is_internally_consistent():
    d = _doc()
    assert d["census"] > 4_000_000
    for k, c in d["partitions"].items():
        assert c["both"] <= min(c["pharm"], c["phys"]), (
            f"{k}: more articles in both classes than in one of them")
        assert c["pharm"] <= d["census"] and c["phys"] <= d["census"]
        assert abs(c["pharm"] / c["phys"] - c["ratio"]) < 1e-6
        assert c["phys_nosurg"] <= c["phys"]
        if c["phys_nosurg"]:
            assert abs(c["pharm"] / c["phys_nosurg"]
                       - c["ratio_surgery_held_out"]) < 1e-6
        assert abs((c["phys"] - c["phys_nosurg"]) / c["phys"]
                   - c["surgical_share_of_physical"]) < 1e-9
    for name, leg in d["holdouts"].items():
        rr = [r["ratio"] for r in leg["partitions"].values() if r["ratio"]]
        if len(rr) > 1:
            assert abs(max(rr) / min(rr) - leg["spread"]) < 1e-9, (
                f"{name}: the reported spread is not max/min of its own ratios")
        for k, r in leg["partitions"].items():
            assert r["n_descriptors"] == len(r["descriptors"])
            assert r["phys_remaining"] + r["phys_held_out"] == \
                d["partitions"][k]["phys"]


def test_assemble_is_arithmetically_right_on_a_synthetic_census():
    """Exercises the derivation without the census, which CI does not have.

    Hand-built incidence: a physical class of one surgical and one radiotherapy
    descriptor, where holding surgery out must leave exactly the radiotherapy
    articles.
    """
    from collections import Counter
    m = _mod()
    parts = {"t": {"pharmacological": ["Antineoplastic Agents"],
                   "physical": ["Mastectomy", "Radiotherapy"]}}
    sets = {"t": {"ph": {"antineoplastic agents"},
                  "py": {"mastectomy", "radiotherapy"}}}
    combos = {"t": Counter({frozenset({"mastectomy"}): 30,
                            frozenset({"radiotherapy"}): 10,
                            frozenset({"mastectomy", "radiotherapy"}): 5})}
    counts = {"t": {"pharm": 90, "phys": 45, "both": 7}}
    d = m.assemble(parts, sets, counts, combos, 1000)
    p = d["partitions"]["t"]
    assert p["ratio"] == 2.0
    # only `mastectomy` is held out, so 10 + 5 articles survive
    assert p["phys_nosurg"] == 15, p
    assert abs(p["ratio_surgery_held_out"] - 6.0) < 1e-9
    assert abs(p["surgical_share_of_physical"] - 30 / 45) < 1e-9
    rt = d["controls"]["radiotherapy family"]["partitions"]["t"]
    # removing `radiotherapy` leaves the 30 mastectomy-only articles plus the 5
    # that carry both, so 35 survive and only 10 are held out
    assert rt["phys_remaining"] == 35, "the radiotherapy control mis-selects"
    assert rt["descriptors"] == ["radiotherapy"]


def test_the_spread_verb_agrees_with_the_numbers():
    """A mutation repointed the rule and the page still said "collapses".

    The renderer emitted the headline unconditionally while its assertions sat
    inside `if held < pub:`, so when the finding failed the guard skipped
    itself. The verb is derived here by an independent rule.
    """
    d, md = _doc(), MD.read_text()
    ps = d["partitions"]
    ratios = [c["ratio"] for c in ps.values()]
    hos = [c["ratio_surgery_held_out"] for c in ps.values()
           if c.get("ratio_surgery_held_out")]
    assert len(hos) == len(ps), "a partition has no held-out ratio"
    pub, held = max(ratios) / min(ratios), max(hos) / min(hos)
    want = ("collapses" if held < pub * 0.9
            else "widens" if held > pub * 1.1 else "barely moves")
    banned = {"collapses", "widens", "barely moves"} - {want}
    sent = [s for s in md.split("\n") if "the spread " + want in s
            or f"the spread {want}" in s]
    assert sent, (
        f"the measured spread goes {pub:.2f}x -> {held:.2f}x, i.e. it {want}, "
        f"and no sentence in the report says so")
    line = sent[0]
    for b in banned:
        assert f"the spread {b}" not in line, (
            f"the same sentence also claims the spread {b}")
    assert f"**{pub:.2f}x**" in md and f"**{held:.2f}x**" in md, (
        "the report does not quote its own before/after spread")
    for k, c in ps.items():
        assert f"{100*c['surgical_share_of_physical']:.1f}%" in md, (
            f"{k}'s surgical share is measured and not rendered")
        assert f"{c['ratio_surgery_held_out']:.2f}:1" in md


def test_the_classifier_keeps_the_two_classes_apart():
    """A mutation had `scan` count an article pharmacological when it matched
    the PHYSICAL class, inflating every ratio, and nothing fired: the only
    guard over `scan` is census-dependent and SKIPS in CI. Tested here on
    synthetic records instead, so it runs everywhere.
    """
    m = _mod()
    s = {"ph": {"antineoplastic agents"}, "py": {"mastectomy", "radiotherapy"}}
    cases = [
        ({"antineoplastic agents"}, True, set()),
        ({"mastectomy"}, False, {"mastectomy"}),
        ({"antineoplastic agents", "radiotherapy"}, True, {"radiotherapy"}),
        ({"appendicitis"}, False, set()),
        ({"mastectomy", "radiotherapy"}, False, {"mastectomy", "radiotherapy"}),
    ]
    for mesh, want_pharm, want_phys in cases:
        got_pharm, got_phys = m.classify(mesh, s)
        assert got_pharm is want_pharm, (
            f"{sorted(mesh)}: pharmacological={got_pharm}, expected "
            f"{want_pharm} -- the classifier is reading the wrong class")
        assert got_phys == want_phys, f"{sorted(mesh)}: physical={got_phys}"


def test_the_spread_verb_is_a_function_of_the_data_not_a_constant():
    """Hardcoding `collapses` is invisible while the spread really collapses.

    A guard that only checks agreement with the CURRENT data cannot see a
    constant, so the renderer is handed a document where holding surgery out
    WIDENS the spread and required to say so.
    """
    m, d = _mod(), _doc()
    doc = json.loads(json.dumps(d))
    ps = doc["partitions"]
    order = sorted(ps, key=lambda k: ps[k]["ratio"])
    # give the LOWEST-ratio partition a huge held-out ratio, so holding surgery
    # out now spreads the five further apart than they started
    ps[order[0]]["ratio_surgery_held_out"] = \
        max(c["ratio_surgery_held_out"] for c in ps.values()) * 6
    ranked = sorted(ps.items(), key=lambda kv: -kv[1]["ratio"])
    lo, hi = ranked[-1][1]["ratio"], ranked[0][1]["ratio"]
    out = "\n".join(m._spread_narrative(doc, ranked, lo, hi))
    assert "the spread widens" in out, (
        "handed data where the hold-out WIDENS the spread, the renderer does "
        f"not say so; it emitted: {out[:300]!r}")
    assert "the spread collapses" not in out, (
        "the renderer claims the spread collapses on data where it widens -- "
        "the verb is a constant, not a reading of the numbers")


def test_the_report_keeps_the_direction_and_refuses_a_replacement_number():
    d, md = _doc(), MD.read_text()
    ratios = [v["ratio"] for v in d["partitions"].values()]
    assert all(r > 1.0 for r in ratios), (
        "a partition puts physical above pharmacological; the report claims "
        "the direction survives everywhere and would need rewriting")
    assert "direction survives and the magnitude does not" in md
    assert "direction survives and the magnitude does not" in SCRIPT.read_text()
    lo, hi = min(ratios), max(ratios)
    assert f"**{lo:.2f}:1 and {hi:.2f}:1**" in md, (
        "the report does not state its own measured range")


def test_the_deliverable_framing_stays_withdrawn():
    md = MD.read_text()
    for m in re.finditer(r"the deliverable", md):
        w = md[max(0, m.start() - 300):m.end() + 300]
        assert re.search(r"earlier version|withdraw", w, re.I), (
            "the report still calls the spread the deliverable")


def test_both_column_is_reported_not_resolved():
    d, md = _doc(), MD.read_text()
    assert all(v["both"] > 0 for v in d["partitions"].values())
    assert "both" in md.lower() and "rather than resolved" in md


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------

def test_the_controls_that_can_discriminate_are_run_and_reported():
    """The per-partition mass match had no power; the total match and the
    allocation permutation do, and both must be present and quoted.
    """
    d, md = _doc(), MD.read_text()
    perm = d.get("permutation") or {}
    tm = perm.get("total_matched") or {}
    alloc = perm.get("allocation") or {}
    m = _mod()
    surg = d["holdouts"][m.PRIMARY_RULE]["spread"]
    assert tm.get("n_draws", 0) >= 1000, (
        f"the total-mass-matched control ran {tm.get('n_draws')} draws")
    assert alloc.get(m.PRIMARY_RULE), "the allocation permutation is gone"
    for name, a in alloc.items():
        assert a["n_feasible_assignments"] >= 2, (
            f"{name}: the allocation permutation has nothing to compare")
        assert 0 <= a["n_at_or_below_observed"] <= a["n_feasible_assignments"]
    pa = alloc[m.PRIMARY_RULE]
    assert f"**{pa['n_at_or_below_observed']} of " in md or \
        f"{pa['n_at_or_below_observed']} of {pa['n_feasible_assignments']}" in md, (
        "the allocation result is computed and not rendered")
    assert f"{tm['n_at_or_below_surgical']:,} of {tm['n_draws']:,}" in md
    # ARITHMETIC INVARIANTS. Thirteen mutations that made the control look
    # better than it is passed the first version of this guard: a flipped
    # tail printed "1,000 of 1,000 reach" directly under "spread 1.62x to
    # 11.61x", a shrunken target printed a mass the table below contradicted,
    # an unshuffled pool printed 1,000 identical draws.
    assert (tm["n_at_or_below_surgical"] > 0) == (tm["min"] <= surg), (
        f"{tm['n_at_or_below_surgical']} draws are reported at or below "
        f"{surg:.4f} while the minimum drawn spread is {tm['min']:.4f}")
    assert (tm["min"] <= tm["median"] <= tm["max"])
    assert tm["min"] < tm["max"], (
        "every draw returned the same spread, so the pool is not being "
        "shuffled and there is one draw wearing a thousand names")
    assert perm["target_total_mass"] == sum(
        (d["removed_mass"].get("surgical") or {}).values()), (
        "the draws are matched to a different mass than the surgical rule "
        "removes, so 'the total surgery removes' is not that total")
    assert abs(perm["surgical_spread"] - surg) < 1e-9, (
        "the count is taken against a different spread than the one the page "
        "prints beside it")
    assert tm["n_draws"] + tm["degenerate_draws"] == m.N_PERMUTATIONS
    for name, a in alloc.items():
        # the observed assignment is one of the permutations, so the count of
        # assignments at or below it can never be zero
        assert a["n_at_or_below_observed"] >= 1, (
            f"{name}: the observed allocation is not counted among the "
            "reassignments it is compared against")
        assert a["n_feasible_assignments"] <= 120
        # RECOMPUTED HERE. Dropping the feasibility filter divides by negative
        # denominators, inflating the count with assignments that hand a
        # partition more articles than its whole physical class.
        import itertools
        masses = [a["masses"][k] for k in sorted(a["masses"])]
        phys = [d["partitions"][k]["phys"] for k in sorted(a["masses"])]
        want = sum(1 for perm in itertools.permutations(masses)
                   if all(m < p for m, p in zip(perm, phys)))
        assert a["n_feasible_assignments"] == want, (
            f"{name}: {a['n_feasible_assignments']} assignments counted "
            f"feasible, {want} actually are -- an infeasible one has no ratio")

    # the withdrawn design must not come back, and its verdict must not either
    for gone in ("mass-matched random, median", "does NOT rest on the "
                 "permutation", "the mass argument cannot explain"):
        assert gone not in md, f"the withdrawn reading {gone!r} is still here"
    assert "could not work" in md and "per-partition" in md.lower()


def test_the_identity_is_called_algebra_and_the_real_invariant_is_checked():
    """The first version computed one expression twice and called the zero
    difference a verification over 25 cells. Substituting a constant for
    `_remaining` left it reporting `holds: True`.
    """
    d = _doc()
    mi = d.get("mass_identity") or {}
    assert mi.get("identity_is_algebra_not_measurement") is True
    assert mi.get("incidence_table_accounts_for_every_physical_article") is True
    per = mi.get("per_partition") or {}
    assert per, "the invariant that CAN fail is not recorded"
    for k, r in per.items():
        assert r["in_incidence_table"] == d["partitions"][k]["phys"], (
            f"{k}: the incidence table holds {r['in_incidence_table']} of the "
            f"{d['partitions'][k]['phys']} physical articles, so the hold-out "
            "arithmetic runs on a different denominator than the main table")
    assert "max error" not in MD.read_text(), (
        "the page still presents the algebraic identity as a measurement")


def test_the_even_split_column_is_over_all_five_or_refuses():
    """Surgery's even share exceeds one partition's entire physical class, and
    the first version silently dropped that partition -- reporting 80% of "the
    identical total" over four of five under a header saying otherwise.
    """
    d, md = _doc(), MD.read_text()
    alloc = (d.get("permutation") or {}).get("allocation") or {}
    for name, a in alloc.items():
        # RECOMPUTED, not read: the first version filtered short partitions
        # out inside the comprehension and reported `short: []` regardless.
        mean = sum(a["masses"].values()) / len(a["masses"])
        assert abs(a["uniform_mass_mean"] - mean) < 1e-6
        short = sorted(k for k in a["masses"]
                       if d["partitions"][k]["phys"] - mean <= 0)
        assert a.get("uniform_mass_short", []) == short, (
            f"{name}: partitions that cannot give up an even share are "
            f"{short}, reported {a.get('uniform_mass_short')}")
        assert a.get("uniform_mass_feasible") is (not short), (
            f"{name}: an even split is reported feasible while {short} hold "
            "fewer physical articles than an even share")
        if a.get("uniform_mass_feasible"):
            assert a["uniform_mass_spread"], f"{name}: feasible and no value"
            continue
        assert a["uniform_mass_spread"] is None
        assert a["uniform_mass_short"], f"{name}: infeasible for no partition"
        for k in a["uniform_mass_short"]:
            assert d["partitions"][k]["phys"] <= a["uniform_mass_mean"], (
                f"{k} is reported unable to give up an even share of "
                f"{a['uniform_mass_mean']:,.0f} and holds "
                f"{d['partitions'][k]['phys']:,} physical articles")
        assert "cannot be spread evenly at all" in md, (
            f"`{name}`'s even split is not achievable and the page does not "
            "say so")


def test_the_named_family_widening_is_not_claimed_as_descriptor_evidence():
    """Their masses are near-uniform, so subtracting them from unequal
    denominators must amplify the highest ratios. An earlier draft called that
    "the part the mass argument cannot explain".
    """
    d, md = _doc(), MD.read_text()
    m = _mod()
    alloc = (d.get("permutation") or {}).get("allocation") or {}
    for name, a in alloc.items():
        if name == m.PRIMARY_RULE:
            continue
        u, obs = a.get("uniform_mass_spread"), a.get("observed_spread")
        if not (u and obs):
            continue
        if a.get("uniform_mass_feasible") is False:
            assert u is None, (
                f"`{name}`'s even split is infeasible and a number is "
                "reported anyway, over fewer partitions than the column says")
            continue
        if u > obs * 0.5:
            assert f"{u:.2f}x of its {obs:.2f}x" in md, (
                f"spreading `{name}`'s mass EVENLY already gives {u:.2f}x of "
                f"its {obs:.2f}x, so most of its widening is arithmetic, and "
                "the page does not say so")
    assert "Subtracting a roughly constant amount from unequal denominators" in md


def test_the_named_controls_report_the_mass_they_remove():
    """They are NOT mass-matched, and the argument depends on knowing that."""
    d, md = _doc(), MD.read_text()
    rm = d.get("removed_mass") or {}
    surg = sum((rm.get("surgical") or {}).values())
    assert surg > 0
    assert rm.get("controls"), "the controls' removed mass is not recorded"
    for name, per in rm["controls"].items():
        mass = sum(per.values())
        assert f"{mass:,}" in md, (
            f"control `{name}` removes {mass:,} physical articles and the "
            "report does not say so, while comparing its spread to surgery's")
        # each control's per-partition mass must be its OWN, not surgery's
        assert per != rm["surgical"] or name == "surgical"
    assert f"{surg:,}" in md, (
        f"surgery removes {surg:,} physical articles and the report does not "
        "state it beside the controls it is compared against")


def test_a_named_family_control_does_not_reproduce_the_collapse():
    """If radiotherapy held out did this too, surgery would not be the story."""
    d = _doc()
    surg = d["holdouts"][_mod().PRIMARY_RULE]["spread"]
    pub = d["published_spread"]
    spreads = {k: v["spread"] for k, v in d["controls"].items() if v.get("spread")}
    assert spreads, "no named-family control produced a spread"
    rm = (d.get("removed_mass") or {}).get("controls") or {}
    for k, v in spreads.items():
        assert v > surg * 1.5, (
            f"holding out `{k}` gives {v:.2f}x against surgery's {surg:.2f}x, "
            f"so the page's attribution to surgery is not supported; the "
            f"unremoved spread is {pub:.2f}x")
        # the load-bearing half: each removes LESS mass than surgery and still
        # WIDENS past the unremoved baseline, which less mass cannot explain
        assert v > pub, (
            f"`{k}` brings the five closer together ({v:.2f}x against "
            f"{pub:.2f}x unremoved), so it is a rival explanation rather than "
            "a control")
        assert sum(rm.get(k, {}).values()) < sum(
            (d["removed_mass"].get("surgical") or {}).values()), (
            f"`{k}` removes at least as much mass as surgery, so 'less mass "
            "and still wider' is not the argument available")


# ---------------------------------------------------------------------------
# The 17.6:1 comparator
# ---------------------------------------------------------------------------

def test_the_comparator_restriction_is_applied_to_both_classes():
    """Narrowing only the numerator is the error this page polices.

    The previous guard checked `precise_phys`, a LABEL list written separately
    from `den_p`, the number actually divided. A mutation changed the divided
    denominator alone and passed. Both the labels and the arithmetic are
    checked here, against the landscape artifact.
    """
    d = _doc()
    lc = d.get("landscape_composition") or {}
    if not lc.get("precise_ratio"):
        pytest.skip("landscape artifact absent")
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "al", REPO_ROOT / "scripts" / "atlas_landscape.py")
    al = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(al)
    assert set(lc["precise_pharm"]) == set(al.PHARMACOLOGICAL & al.PRECISE)
    assert set(lc["precise_phys"]) == set(al.PHYSICAL & al.PRECISE), (
        "the restriction is not applied to the physical class")
    assert lc["dropped_phys"], "nothing is dropped from the denominator"

    rows = json.loads((REPO_ROOT / "analysis" / "atlas-landscape.json").read_text())
    rows = (rows if isinstance(rows, list)
            else rows.get("mechanisms") or rows.get("rows") or [])
    assert rows, "the landscape artifact has no mechanism rows"
    cen = {r["mechanism"].lower(): r.get("mesh_census") or 0
           for r in rows if r.get("mechanism")}
    assert lc["precise_numerator"] == sum(cen[x] for x in lc["precise_pharm"]), (
        "the numerator divided is not the sum over the labelled numerator set")
    assert lc["precise_denominator"] == sum(cen[x] for x in lc["precise_phys"]), (
        "THE DENOMINATOR DIVIDED IS NOT THE SUM OVER THE LABELLED DENOMINATOR "
        "SET -- the restriction is one-sided in the arithmetic while the "
        "labels claim it is symmetric")
    assert abs(lc["precise_ratio"]
               - lc["precise_numerator"] / lc["precise_denominator"]) < 1e-9
    assert abs(lc["ratio"] - lc["numerator"] / lc["denominator"]) < 1e-9


def test_the_comparator_placement_sentence_is_derived():
    """"within the range" was hand-written and false: 6.43 against 3.83-4.70."""
    d, md = _doc(), MD.read_text()
    lc = d.get("landscape_composition") or {}
    if not lc.get("precise_ratio"):
        pytest.skip("landscape artifact absent")
    hos = sorted(c["ratio_surgery_held_out"] for c in d["partitions"].values()
                 if c.get("ratio_surgery_held_out"))
    pr = lc["precise_ratio"]
    inside = hos[0] <= pr <= hos[-1]
    if inside:
        assert "is inside the" in md and "above the" not in md.split(
            "So the comparator falls")[1][:400]
    else:
        assert "narrows the gap without closing it" in md, (
            f"the restricted comparator is {pr:.2f}:1 and the held-out range "
            f"is {hos[0]:.2f}-{hos[-1]:.2f}:1, i.e. OUTSIDE it, and the report "
            "claims otherwise")
        assert f"{hos[0]:.2f}-{hos[-1]:.2f}:1" in md


def test_the_precise_set_is_not_presented_as_the_criterion_applied_evenly():
    """Two excluded mechanisms satisfy PRECISE's own stated criterion."""
    d, md = _doc(), MD.read_text()
    lc = d.get("landscape_composition") or {}
    if not lc.get("precise_ratio"):
        pytest.skip("landscape artifact absent")
    # Derived HERE, from the landscape artifact's own top descriptors, so
    # `restored = []` in the generator fails instead of silently deleting the
    # section. `Immune Checkpoint Inhibitors` and `Poly(ADP-ribose) Polymerase
    # Inhibitors` name drug classes; `DNA Methylation` and `Glycolysis` name
    # processes.
    tops = lc.get("top_descriptors") or {}
    independently = sorted(
        x for x in lc["dropped_pharm"]
        if re.search(r"inhibitor|therap|antibod|conjugat|vaccine|agent",
                     tops.get(x) or "", re.I))
    assert independently, (
        "no excluded pharmacological mechanism has a therapy-naming dominant "
        "descriptor any more; the page's argument that PRECISE is not its own "
        "criterion applied evenly needs redoing rather than deleting")
    # THE SAME READING MUST BE APPLIED TO THE DENOMINATOR. Testing only the
    # numerator's dropped mechanisms is the one-sided narrowing this page
    # exists to police, structurally, even while the number is unaffected.
    independently_py = sorted(
        x for x in lc["dropped_phys"]
        if re.search(r"inhibitor|therap|antibod|conjugat|vaccine|agent",
                     tops.get(x) or "", re.I))
    assert set(lc.get("criterion_restored_phys") or []) == set(independently_py), (
        f"the generator restores {sorted(lc.get('criterion_restored_phys') or [])} "
        f"on the denominator side but the same reading gives {independently_py}")
    assert abs(lc["criterion_restored_ratio"]
               - lc["criterion_restored_numerator"]
               / lc["criterion_restored_denominator"]) < 1e-9
    rest = lc.get("criterion_restored_pharm") or []
    assert set(rest) == set(independently), (
        f"the generator restores {sorted(rest)} but the landscape artifact's "
        f"own top descriptors say {independently} satisfy PRECISE's stated "
        "criterion")
    assert lc["criterion_restored_ratio"], "the restored ratio was not computed"
    assert f"**{lc['criterion_restored_ratio']:.2f}:1**" in md, (
        "restoring the criterion-satisfying mechanisms changes the restricted "
        "ratio and the report does not say so")
    for x in rest:
        assert f"`{x}`" in md


def test_the_two_quoted_ratios_are_pinned_to_their_sources():
    """9.1 and 17.6 are typed constants rendered into prose beside numbers the
    page derives. Pinned to the manuscript and to the landscape artifact, so
    an edit there fails here rather than silently changing what is compared.
    """
    m = _mod()
    v1 = (REPO_ROOT / "article" / "drafts" / "v1.md")
    if v1.exists():
        assert f"{m.MANUSCRIPT_RATIO}:1" in v1.read_text(), (
            f"the manuscript no longer states {m.MANUSCRIPT_RATIO}:1, so the "
            "constant this page divides by is stale")
    lc = _doc().get("landscape_composition") or {}
    if lc.get("ratio"):
        assert abs(lc["ratio"] - m.LANDSCAPE_CENSUS_RATIO) < 0.1, (
            f"the landscape artifact gives {lc['ratio']:.2f}:1 while this "
            f"page types {m.LANDSCAPE_CENSUS_RATIO}:1 beside it")


def test_the_criterion_restoration_is_symmetric_on_synthetic_data():
    """Today the denominator side restores nothing, so a mutation dropping it
    is INERT and no artifact-reading guard can see it. Exercised on data where
    the physical side DOES have a therapy-naming dropped mechanism.
    """
    m = _mod()
    cen = {"drugA": 100, "drugB": 50, "physA": 10, "physB": 5}
    top = {"drugA": "Widget Inhibitors", "drugB": "Glycolysis",
           "physA": "Ultrasonic Therapy", "physB": "Electrochemotherapy"}
    ph, py, pre = {"drugA", "drugB"}, {"physA", "physB"}, {"drugA", "physA"}
    r = m.restrict(cen, top, ph, py, pre)
    # `drugB`/Glycolysis is a process -> stays dropped.
    # `physB`/Electrochemotherapy names a therapy -> the SAME reading restores it.
    assert r["criterion_restored_pharm"] == []
    assert r["criterion_restored_phys"] == ["physB"], (
        "the denominator's therapy-naming dropped mechanism is not restored, "
        "so the restoration is applied to the numerator only -- the exact "
        "one-sidedness this page exists to police")
    assert r["criterion_restored_numerator"] == 100
    assert r["criterion_restored_denominator"] == 15, (
        f"denominator {r['criterion_restored_denominator']}, expected 10+5")
    assert abs(r["criterion_restored_ratio"] - 100 / 15) < 1e-9
    assert abs(r["precise_ratio"] - 100 / 10) < 1e-9
    # and a therapy-naming NUMERATOR drop is restored too, symmetrically
    r2 = m.restrict(cen, {**top, "drugB": "Widget Antibodies"}, ph, py, pre)
    assert r2["criterion_restored_pharm"] == ["drugB"]
    assert r2["criterion_restored_numerator"] == 150


def test_the_comparator_consequence_for_the_understatement_claim_is_stated():
    """6.43:1 is below the manuscript's 9.1:1, which inverts a claim in 4 sites."""
    d, md = _doc(), MD.read_text()
    m = _mod()
    lc = d.get("landscape_composition") or {}
    if not lc.get("precise_ratio"):
        pytest.skip("landscape artifact absent")
    if lc["precise_ratio"] < m.MANUSCRIPT_RATIO:
        assert "BELOW the manuscript's own" in md, (
            f"the restricted comparator is {lc['precise_ratio']:.2f}:1, below "
            f"the manuscript's {m.MANUSCRIPT_RATIO}:1, so the repo's "
            "'understates its own case' framing is inverted under this "
            "restriction and the page must say so")
        for site in ("article/drafts/v1.md", "analysis/census-findings.md"):
            assert site in md, f"the affected site {site} is not named"


# ---------------------------------------------------------------------------
# Provenance and the qualifier caveat
# ---------------------------------------------------------------------------

def test_the_unsupported_provenance_is_withdrawn():
    """No artifact records a reviewer, a count, or a stated principle."""
    md, src = MD.read_text(), SCRIPT.read_text()
    for k, v in _parts().items():
        assert set(v) <= {"pharmacological", "physical"}, (
            f"{k} now carries {sorted(set(v) - {'pharmacological','physical'})}; "
            "if a principle or reviewer record has been added, the withdrawn "
            "clause can be restored")
    for phrase in ("independent reviewer", "reproduce its own count",
                   "stated principle", "single stated"):
        for text, where in ((md, "report"), (src, "docstring")):
            for m in re.finditer(re.escape(phrase), text, re.I):
                w = text[max(0, m.start() - 500):m.end() + 500]
                assert re.search(r"withdraw|earlier version|NO ARTIFACT|not "
                                 r"restated", w, re.I), (
                    f"the {where} claims the partitions were {phrase!r} and no "
                    "artifact in the repo supports it")


def test_the_qualifier_caveat_reports_every_row_and_settles_nothing():
    """The page once quoted the pair that suited it and dropped the sharpest row."""
    d, md = _doc(), MD.read_text()
    q = (d.get("qualifier_recalls") or {}).get("modalities") or {}
    src = REPO_ROOT / "analysis" / "atlas-ingest-sensitivity.json"
    if not src.exists():
        pytest.skip("#722 artifact absent")
    raw = json.loads(src.read_text())["modalities"]
    assert set(q) == set(raw), (
        f"the page reports {sorted(q)} of the {sorted(raw)} rows #722 "
        "measured; selecting rows is how the invalid inference was made")
    for m_, c in q.items():
        assert abs(c["recall"] - raw[m_]["descriptor"] / raw[m_]["either"]) < 1e-9, (
            f"{m_}: the recall is not descriptor/either from the raw counts")
        assert f"{c['recall']:.3f}" in md, f"{m_}'s recall is not rendered"
    sharpest = min(q.items(), key=lambda kv: kv[1]["recall"])[0]
    assert sharpest == "surgery", (
        "the sharpest row is no longer surgery; the page's argument that the "
        "caveat may run DOWN for a surgery-heavy physical class needs redoing")
    assert "the wrong pair for these classes" in md, (
        "the page infers a direction from the drug-therapy/radiotherapy pair "
        "without the surgery row that points the other way")
    assert "direction is not established here" in md
    assert "THAT INFERENCE WAS INVALID" in md


# ---------------------------------------------------------------------------
# Freshness, refusals, and the scan contract
# ---------------------------------------------------------------------------

def test_an_empty_match_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'v["phys"] == 0' in src, "the empty check no longer tests the count"
    assert "is not a finding" in src and "raise SystemExit" in src


def test_the_committed_report_is_what_the_generator_produces():
    """Regenerate-and-diff. Every other guard compares the .md to the .json,
    and the two go stale TOGETHER, so none of them can fail on a renderer edit.
    """
    m = _mod()
    assert m.render(_doc()) == MD.read_text(), (
        "analysis/atlas-modality-ratio.md is not what the current renderer "
        "produces from the committed JSON -- re-run "
        "`python scripts/atlas_modality_ratio.py --render-only`")


def test_render_only_works_without_the_census():
    import shutil
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        shutil.copy2(MD, Path(td) / MD.name)
        try:
            res = subprocess.run([sys.executable, str(SCRIPT), "--render-only"],
                                 cwd=REPO_ROOT, capture_output=True, text=True)
        finally:
            shutil.copy2(Path(td) / MD.name, MD)
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"


def test_the_scan_reproduces_the_committed_holdouts_on_a_sample():
    """SCAN-CONTRACT. `--render-only` cannot see a change inside `scan`, so a
    mutation there is invisible to every artifact-reading guard. Runs the real
    scan over a few shards -- STRIDED, since shards are chronological -- and
    checks that the descriptor sets each rule selects are the same ones the
    committed artifact reports.
    """
    if not RECORDS.exists():
        pytest.skip("census not present (gitignored); CI reads artifacts only")
    import gzip
    from collections import Counter
    m = _mod()
    parts = m.load_partitions()
    sets = {k: {"ph": {x.lower() for x in v["pharmacological"]},
                "py": {x.lower() for x in v["physical"]}}
            for k, v in parts.items()}
    counts = {k: {"pharm": 0, "phys": 0, "both": 0} for k in sets}
    combos = {k: Counter() for k in sets}
    n = 0
    shards = sorted(RECORDS.glob("*.jsonl.gz"))[::400]
    assert shards, "no shards"
    for f in shards:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                mesh = {x.lower() for x in (r.get("mesh") or [])}
                if not mesh:
                    continue
                for k, s in sets.items():
                    a, hit = m.classify(mesh, s)
                    counts[k]["pharm"] += a
                    if hit:
                        counts[k]["phys"] += 1
                        counts[k]["both"] += a
                        combos[k][frozenset(hit)] += 1
    got = m.assemble(parts, sets, counts, combos, n)
    want = _doc()
    for name, leg in want["holdouts"].items():
        for k, r in leg["partitions"].items():
            assert got["holdouts"][name]["partitions"][k]["descriptors"] == \
                r["descriptors"], (
                f"a live scan selects different descriptors for {name}/{k} "
                "than the committed artifact reports")
    for k in want["partitions"]:
        assert got["partitions"][k]["n_phys_descriptors"] == \
            want["partitions"][k]["n_phys_descriptors"]
    # the sample must actually exercise the physical class
    assert sum(c["phys"] for c in counts.values()) > 100


def test_the_surgery_recall_is_shown_beside_its_family_control():
    """This page leaned on surgery's 0.101 recall as the reason the qualifier
    bias might run DOWN. #722 now recomputes the descriptor arm as the whole
    MeSH family each modality names, and that figure was a property of a
    four-entry list rather than of indexing.
    """
    d, md = _doc(), MD.read_text()
    q = (d.get("qualifier_recalls") or {}).get("modalities") or {}
    src = REPO_ROOT / "analysis" / "atlas-ingest-sensitivity.json"
    if not src.exists():
        pytest.skip("#722 artifact absent")
    j = json.loads(src.read_text())
    tree = j.get("modalities_tree_arm")
    if not tree:
        pytest.skip("#722 has not yet published the family arm")
    for m, c in q.items():
        assert "recall_tree_arm" in c, (
            f"{m}: the family-arm recall is available and not carried, so the "
            "page quotes a proxy figure with no control beside it")
        t = tree[m]
        assert abs(c["recall_tree_arm"] - t["descriptor"] / t["either"]) < 1e-9
        assert f"{c['recall_tree_arm']:.3f}" in md
    sg = q.get("surgery")
    if sg and sg["recall_tree_arm"] > sg["recall"] * 2:
        assert "A PROPERTY OF A FOUR-ENTRY LIST" in md, (
            f"surgery's recall moves {sg['recall']:.3f} -> "
            f"{sg['recall_tree_arm']:.3f} under the control and the page still "
            "quotes the proxy figure as its reason")
        # the conclusion must NOT have moved -- it never rested on the figure
        assert "direction is not established here" in md


def test_the_direction_flip_is_derived_from_both_arms():
    """The page refuses to apply a qualifier correction. The strongest reason
    is not "unmeasured" but that the two descriptor arms give OPPOSITE
    directions, so the refusal must be computed rather than asserted -- and it
    must stop being asserted if they ever agree.
    """
    d, md = _doc(), MD.read_text()
    q = (d.get("qualifier_recalls") or {}).get("modalities") or {}
    sg, dt = q.get("surgery"), q.get("drug therapy")
    if not (sg and dt and sg.get("recall_tree_arm") and dt.get("recall_tree_arm")):
        pytest.skip("#722 has not published the family arm")
    fp = sg["recall"] / dt["recall"]
    ft = sg["recall_tree_arm"] / dt["recall_tree_arm"]
    if (fp - 1) * (ft - 1) < 0:
        assert "OPPOSITE DIRECTIONS" in md, (
            f"the proxy arm multiplies the ratio by {fp:.2f} and the family "
            f"arm by {ft:.2f} -- opposite sides of 1 -- and the report does "
            "not say so")
        assert f"**x{fp:.2f}**" in md and f"**x{ft:.2f}**" in md
        assert "no correction is applied" in md
    else:
        assert "OPPOSITE DIRECTIONS" not in md, (
            f"both arms point the same way ({fp:.2f}, {ft:.2f}) and the page "
            "still claims they conflict")
    # either way the verdict must stay refused: a flip is a reason not to
    # correct, never a licence to correct in the newer direction
    assert "direction is not established here" in md
