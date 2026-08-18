#!/usr/bin/env python3
"""Does the pharmacological:physical claim survive a symmetric class definition? (#724)

THE CLAIM
---------
The manuscript's central corpus claim is that pharmacological cancer research
vastly outweighs physical-modality research. `atlas_landscape.py` computes it
from two curated sets of this project's own mechanism tags and reports 9.1:1 by
the manuscript's method and 17.6:1 on the census, reading the census figure as
STRENGTHENING the case.

WHY THE CLASSES CANNOT SETTLE IT AS WRITTEN
--------------------------------------------
    PHYSICAL = {hifu, sonodynamic, electrochemical-therapy}
    PHARMACOLOGICAL = {immunotherapy, car-t, antibody-drug-conjugate,
                       bispecific-antibody, synthetic-lethality, epigenetic,
                       metabolic-targeting}

Both omit their largest real-world member. PHYSICAL has no radiotherapy -- the
physical modality delivered to roughly half of all cancer patients -- and
PHARMACOLOGICAL has no cytotoxic chemotherapy. Both are missing for the same
reason: neither has a mechanism tag in this project (see #723), not because
either was weighed and excluded.

An earlier proposal was to add radiotherapy to PHYSICAL. That is a ONE-SIDED
WIDENING and it manufactures whatever inversion it finds. A recomputation has to
move both sides by one principle, or neither.

WHAT THIS DOES
--------------
Five partitions were constructed to move both classes together, using MeSH
descriptors -- where radiotherapy and chemotherapy ARE measurable even though
this project has no tags for them. Two claims once made here about their
provenance are WITHDRAWN and are not restated: that each was "checked by a
separate reviewer" and "required to reproduce its own count" (no artifact in
this repo records either), and that each was built "by a single STATED
principle" (`modality-partitions.json` carries two member lists per partition
and no statement of principle, so the reader has the memberships and not the
reasoning).

The partitions are the INPUT (`analysis/modality-partitions.json`), committed
with their member lists so any placement can be disputed. This script is the
measurement. The panels' own reported ratios are deliberately NOT carried in
that file: an input that contains the result it is supposed to produce is how a
measurement ends up validating itself.

WHAT THE ANSWER IS FOR
----------------------
Not to replace 17.6:1 with a better number. The five partitions span 1.32:1 to
3.93:1, and the FIRST version of this page called that spread the deliverable.
It is not: the five differ almost entirely on one membership question nobody
declared -- whether operative surgery is a physical modality -- and holding
surgery out collapses them onto each other.

THE HOLD-OUT RULE IS ITSELF A JUDGEMENT, WHICH IS THE SAME KIND OF THING THE
PAGE JUST FOUND. Deciding what counts as operative surgery is exactly as
disputable as deciding whether it is physical, so naming one regex and calling
the question settled would move the undeclared variable up a level rather than
remove it. Five rules are therefore run, each disagreeing with the others about
a named, listed set of descriptors, and the collapse is reported as a range
across all five rather than as one number. Two things are measured instead of
argued:

  * what each rule LEAKS. The stem list catches `Radiosurgery`, `Cryosurgery`,
    `Electrosurgery` and `Ultrasonic Surgical Procedures` (energy modalities,
    not operative removal) and nine transplantation descriptors (infused or
    grafted tissue), and misses `Reoperation`, `Castration`, `Curettage`,
    `Eye Enucleation`, `Pelvic Exenteration` and others.

  * whether surgery is SPECIFIC or merely LARGE. Removing any comparable mass
    of physical descriptors shrinks the spread somewhat, so a mass-matched
    permutation control is run and reported beside the surgical rules, and it
    reaches the same place, so what the collapse establishes is that the
    AMOUNT of surgical mass differs across the partitions -- not that the
    surgical descriptors are special as descriptors. The controls that DO
    speak to specificity are named non-surgical families, and they widen the
    spread instead.

A MEASURED CAVEAT WHOSE DIRECTION IS NOT ESTABLISHED. The ingest reads MeSH
DescriptorName and never QualifierName (#722), so both classes are understated.
Which class is understated MORE is not settled here, and the page reports every
row of that artifact rather than the pair that suits its conclusion.

Usage:
    python scripts/atlas_modality_ratio.py
    python scripts/atlas_modality_ratio.py --render-only
"""

import argparse
import gzip
import json
import random
import re
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = PROJECT_ROOT / "corpus" / "atlas"
PARTITIONS = PROJECT_ROOT / "analysis" / "modality-partitions.json"
INGEST = PROJECT_ROOT / "analysis" / "atlas-ingest-sensitivity.json"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-modality-ratio.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-modality-ratio.json"

# What the manuscript and the existing landscape analysis report, for contrast.
MANUSCRIPT_RATIO = 9.1
LANDSCAPE_CENSUS_RATIO = 17.6

# Draws for the mass-matched permutation control, and its seed. Fixed here so
# the control is reproducible and cannot be re-rolled until it agrees.
N_PERMUTATIONS = 20
PERMUTATION_SEED = 20260818


def load_partitions(d: dict | None = None) -> dict:
    """Load the committed partitions, refusing an input that carries its result.

    Takes an optional dict so the refusal can be exercised on a doctored input
    rather than asserted about the source text: a guard that greps for the
    condition passes while `if False and (...)` disables it.
    """
    if d is None:
        d = json.loads(PARTITIONS.read_text())
    if not d:
        raise SystemExit(f"no partitions in {PARTITIONS}")
    for name, spec in d.items():
        if "ratio" in spec or "reported_ratio" in spec:
            raise SystemExit(
                f"{name} carries a ratio in the INPUT file. The input must not "
                "contain the result; that is how a measurement validates "
                "itself.")
        if not spec.get("pharmacological") or not spec.get("physical"):
            raise SystemExit(f"{name} is missing a class")
    return d


# ---------------------------------------------------------------------------
# THE ONE UNDECLARED VARIABLE, AND THE JUDGEMENT INSIDE HOLDING IT OUT.
#
# The five partitions were presented as five independent readings whose SPREAD
# was the deliverable. They are not independent: they differ almost entirely in
# whether operative surgery counts as a "physical modality".
#
# Holding surgery out needs a definition of surgery, and there is no
# uncontroversial one. The stems below are a hand-written list -- the same kind
# of artefact as the partitions themselves -- so what it catches and misses is
# enumerated rather than trusted, and four alternative rules disagree with it on
# named sets.
# ---------------------------------------------------------------------------
STEM_ALTERNATIVES = [
    "surg", "ectomy", "otomy", "ostomy", "resect", "excis", "laparoscop",
    "endoscop", "transplant", "amputat", "dissect", "anastomos", "graft",
    "implantation", "reconstructi", "debulk", "lymphadenectom",
]
SURGICAL = re.compile("|".join(STEM_ALTERNATIVES), re.I)

# Descriptors the stem list catches that deliver ENERGY through a beam or a
# probe rather than removing tissue operatively. `Radiosurgery` is stereotactic
# radiotherapy, so holding it out removes radiotherapy mass from the physical
# class -- the opposite of what this page argues was wrongly excluded.
ENERGY_CAUGHT_BY_STEM = {
    "radiosurgery", "cryosurgery", "electrosurgery",
    "ultrasonic surgical procedures",
}
# Cellular and organ transplantation: infused or grafted, not an operative
# removal, and two of the five partitions file several of these as
# PHARMACOLOGICAL rather than physical.
TRANSPLANT = re.compile(r"transplant", re.I)
# Operative procedures no stem reaches. Listed rather than patched into the
# regex, so the disagreement between rules stays visible.
OPERATIVE_MISSED_BY_STEM = {
    "castration", "conization", "curettage", "dilatation and curettage",
    "eye enucleation", "limb salvage", "orthopedic procedures",
    "pelvic exenteration", "reoperation", "vertebroplasty",
}

# CONTROLS. Two named non-surgical families of comparable construction, so
# "removing a family collapses the spread" can be tested against families that
# are not surgery.
RADIOTHERAPY_FAMILY = re.compile(
    r"radiotherap|irradiat|brachytherap|proton therapy|x-ray therapy|"
    r"radioisotope|radiation|boron neutron|radioimmunotherapy", re.I)
ENERGY_ABLATION_FAMILY = re.compile(
    r"ablat|hyperthermi|laser|photo|puva|cryotherap|electropor|"
    r"radiofrequency|diatherm|ultraviolet|light therapy|coagulation", re.I)


def _stem(members: set) -> set:
    return {x for x in members if SURGICAL.search(x)}


def holdout_rules() -> dict:
    """name -> f(physical member set) -> the subset that rule holds out.

    Every rule is a function of the member set alone, so it cannot consult the
    partition it came from. A per-partition exception is the defect this page
    measures, one level up.
    """
    def no_energy(s):
        return _stem(s) - ENERGY_CAUGHT_BY_STEM

    def no_transplant(s):
        return {x for x in _stem(s) if not TRANSPLANT.search(x)}

    def plus_missed(s):
        return _stem(s) | (s & OPERATIVE_MISSED_BY_STEM)

    def operative_only(s):
        return ((no_energy(s) & no_transplant(s)) | (s & OPERATIVE_MISSED_BY_STEM))

    return {
        "operative-only": operative_only,
        "stem-list": _stem,
        "stem-list, energy modalities restored": no_energy,
        "stem-list, transplantation restored": no_transplant,
        "stem-list plus missed operative terms": plus_missed,
    }


# The rule the main table's held-out column uses. The strictest of the five:
# it excludes everything any of the others excludes and includes the operative
# terms the stem list misses. It is not privileged -- the spread is reported
# across all five and they differ by less than the collapse itself.
PRIMARY_RULE = "operative-only"

CONTROL_RULES = {
    "radiotherapy family": lambda s: {x for x in s if RADIOTHERAPY_FAMILY.search(x)},
    "energy/ablation family": lambda s: {x for x in s
                                         if ENERGY_ABLATION_FAMILY.search(x)},
}


def _spread(vals) -> float | None:
    vals = [v for v in vals if v]
    if len(vals) < 2:
        return None
    return max(vals) / min(vals)


def _spearman(xs, ys) -> float | None:
    n = len(xs)
    if n < 3:
        return None

    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else None


def scan(parts: dict) -> dict:
    """One pass over the census, storing the physical class's INCIDENCE STRUCTURE.

    For each partition the scan records, per distinct set of matched physical
    descriptors, how many articles matched exactly that set. Every hold-out
    question -- how many articles still match after removing descriptors X --
    is then answerable offline from that table, so alternative rules and a
    permutation control cost no extra passes and, more importantly, are
    computed from the SAME scan rather than from a second one that could differ.
    """
    sets = {k: {"ph": {x.lower() for x in v["pharmacological"]},
                "py": {x.lower() for x in v["physical"]}}
            for k, v in parts.items()}
    counts = {k: {"pharm": 0, "phys": 0, "both": 0} for k in sets}
    combos = {k: Counter() for k in sets}
    n = 0
    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                mesh = {m.lower() for m in (r.get("mesh") or [])}
                if not mesh:
                    continue
                for k, s in sets.items():
                    a = bool(mesh & s["ph"])
                    counts[k]["pharm"] += a
                    hit = mesh & s["py"]
                    if hit:
                        counts[k]["phys"] += 1
                        counts[k]["both"] += a
                        combos[k][frozenset(hit)] += 1
    return assemble(parts, sets, counts, combos, n)


def _remaining(combo: Counter, removed: set) -> int:
    """Articles still matching the physical class once `removed` is held out."""
    return sum(v for fs, v in combo.items() if fs - removed)


def assemble(parts, sets, counts, combos, n) -> dict:
    rules = holdout_rules()
    out = {
        "census": n,
        "partitions": {},
        "holdouts": {},
        "controls": {},
        "permutation": {},
        "stem_alternative_hits": {},
        "landscape_composition": landscape_composition(),
        "qualifier_recalls": qualifier_recalls(),
    }

    universe = set().union(*(s["py"] for s in sets.values()))
    for alt in STEM_ALTERNATIVES:
        out["stem_alternative_hits"][alt] = sorted(
            x for x in universe if re.search(alt, x, re.I))

    primary = {k: rules[PRIMARY_RULE](s["py"]) for k, s in sets.items()}
    for k, c in counts.items():
        held = _remaining(combos[k], primary[k])
        out["partitions"][k] = {
            **c,
            "ratio": (c["pharm"] / c["phys"]) if c["phys"] else 0.0,
            "n_pharm_descriptors": len(parts[k]["pharmacological"]),
            "n_phys_descriptors": len(parts[k]["physical"]),
            "n_surgical_descriptors": len(primary[k]),
            "phys_nosurg": held,
            "surgical_share_of_physical":
                (c["phys"] - held) / c["phys"] if c["phys"] else None,
            "ratio_surgery_held_out": (c["pharm"] / held) if held else None,
        }

    def leg(rule):
        rows = {}
        for k, s in sets.items():
            rem = rule(s["py"])
            held = _remaining(combos[k], rem)
            rows[k] = {
                "n_descriptors": len(rem),
                "descriptors": sorted(rem),
                "phys_held_out": counts[k]["phys"] - held,
                "phys_remaining": held,
                "share_of_physical":
                    (counts[k]["phys"] - held) / counts[k]["phys"]
                    if counts[k]["phys"] else None,
                "ratio": (counts[k]["pharm"] / held) if held else None,
            }
        return rows

    for name, rule in rules.items():
        rows = leg(rule)
        out["holdouts"][name] = {
            "partitions": rows,
            "spread": _spread([r["ratio"] for r in rows.values()]),
            "spearman_share_vs_published_ratio": _spearman(
                [rows[k]["share_of_physical"] or 0.0 for k in sets],
                [out["partitions"][k]["ratio"] for k in sets]),
        }
    for name, rule in CONTROL_RULES.items():
        rows = leg(rule)
        out["controls"][name] = {
            "partitions": rows,
            "spread": _spread([r["ratio"] for r in rows.values()]),
        }

    out["permutation"] = permutation_control(sets, counts, combos, primary)
    out["published_spread"] = _spread(
        [v["ratio"] for v in out["partitions"].values()])
    return out


def permutation_control(sets, counts, combos, primary) -> dict:
    """Would removing ANY comparable mass of physical descriptors do this?

    The surgical rules remove a large slice of the physical class, and a large
    removal shrinks a spread on its own. Each draw removes RANDOM physical
    descriptors from each partition until it has held out at least as many
    articles as the surgical rule did there, so the draws are mass-matched
    partition by partition, and the spread across the five is recomputed. The
    surgical collapse is only about surgery to the extent that it beats these.
    """
    rng = random.Random(PERMUTATION_SEED)
    target = {k: counts[k]["phys"] - _remaining(combos[k], primary[k])
              for k in sets}
    draws = []
    for _ in range(N_PERMUTATIONS):
        ratios, sizes = [], []
        for k, s in sets.items():
            pool = sorted(s["py"])
            rng.shuffle(pool)
            chosen, i = set(), 0
            while i < len(pool):
                if counts[k]["phys"] - _remaining(combos[k], chosen) >= target[k]:
                    break
                chosen.add(pool[i])
                i += 1
            held = _remaining(combos[k], chosen)
            ratios.append((counts[k]["pharm"] / held) if held else None)
            sizes.append(len(chosen))
        sp = _spread(ratios)
        if sp:
            draws.append({"spread": sp, "n_descriptors": sizes})
    spreads = sorted(d["spread"] for d in draws)
    surgical = _spread([counts[k]["pharm"] / _remaining(combos[k], primary[k])
                        for k in sets if _remaining(combos[k], primary[k])])
    return {
        "n_draws": len(draws),
        "seed": PERMUTATION_SEED,
        "spreads": spreads,
        "median": spreads[len(spreads) // 2] if spreads else None,
        "min": spreads[0] if spreads else None,
        "max": spreads[-1] if spreads else None,
        "surgical_spread": surgical,
        "n_draws_at_or_below_surgical":
            sum(1 for s in spreads if surgical and s <= surgical),
        "articles_held_out_target": target,
    }


def qualifier_recalls() -> dict:
    """Every row of #722, not the pair that suits the conclusion.

    Computed from the raw counts rather than the rounded percentages the first
    version of this page quoted: 3.3/5.7 and 12.4/22.6 are two significant
    figures and were being divided as if they were the measurement.
    """
    if not INGEST.exists():
        return {}
    d = json.loads(INGEST.read_text())
    rows = {}
    for m, c in (d.get("modalities") or {}).items():
        desc, either = c.get("descriptor"), c.get("either")
        if desc is None or not either:
            continue
        rows[m] = {"descriptor": desc, "either": either,
                   "recall": desc / either}
    return {"modalities": rows, "cancer_articles": d.get("cancer_articles"),
            "n_shards": d.get("n_shards")}


def landscape_composition() -> dict:
    """What the 17.6:1 comparator is made of, and the symmetric restriction.

    Its numerator is dominated by two mechanisms `atlas_landscape.py`'s own
    text calls a SCOPE ARTIFACT rather than a therapy. Dropping those from the
    numerator ALONE would be the one-sided narrowing this page exists to
    police -- the denominator's `electrochemical-therapy` is equally
    non-precise. So the restriction is applied to BOTH classes, using that
    script's own PRECISE set rather than a judgement made here.

    PRECISE was not written for this question, and two of its exclusions do not
    follow from its own stated criterion, so the ratio under a criterion-faithful
    restoration is computed alongside rather than instead.
    """
    src = PROJECT_ROOT / "analysis" / "atlas-landscape.json"
    if not src.exists():
        return {}
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "al", PROJECT_ROOT / "scripts" / "atlas_landscape.py")
    al = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(al)
    d = json.loads(src.read_text())
    rows = d if isinstance(d, list) else d.get("mechanisms") or d.get("rows") or []
    cen, top = {}, {}
    for r in rows:
        if isinstance(r, dict) and r.get("mechanism"):
            cen[r["mechanism"].lower()] = r.get("mesh_census") or 0
            top[r["mechanism"].lower()] = r.get("top_descriptor")
    if not cen:
        return {}

    def tot(names):
        return sum(cen.get(x, 0) for x in names)

    ph, py, pre = al.PHARMACOLOGICAL, al.PHYSICAL, al.PRECISE
    num, den = tot(ph), tot(py)
    pre_ph, pre_py = sorted(ph & pre), sorted(py & pre)
    num_p, den_p = tot(pre_ph), tot(pre_py)
    # PRECISE's own criterion is "the dominant descriptor names a therapy or
    # modality, rather than a process, a material or a technique". These two
    # pharmacological mechanisms satisfy it and are excluded anyway, so the
    # restriction is not the criterion applied evenly.
    restored = sorted(x for x in (ph - pre) if _names_a_therapy(top.get(x)))
    num_r = num_p + tot(restored)
    biggest = sorted(((x, cen.get(x, 0)) for x in ph), key=lambda kv: -kv[1])[:2]
    return {
        "numerator": num, "denominator": den,
        "ratio": num / den if den else None,
        "top_two_numerator": biggest,
        "top_two_share": sum(v for _k, v in biggest) / num if num else None,
        "precise_numerator": num_p, "precise_denominator": den_p,
        "precise_ratio": num_p / den_p if den_p else None,
        "precise_pharm": pre_ph, "precise_phys": pre_py,
        "dropped_pharm": sorted(ph - pre), "dropped_phys": sorted(py - pre),
        "top_descriptors": {x: top.get(x) for x in sorted(ph | py)},
        "criterion_restored_pharm": restored,
        "criterion_restored_numerator": num_r,
        "criterion_restored_ratio": num_r / den_p if den_p else None,
    }


# Descriptor names that ARE a therapy or a modality by PRECISE's own criterion:
# an agent class or a named procedure, rather than a process, a material or a
# laboratory technique.
_THERAPY_NAME = re.compile(r"inhibitor|therap|antibod|conjugat|vaccine|"
                           r"agents?\b|blockade|transplant|surgery|ablation",
                           re.I)


def _names_a_therapy(descriptor) -> bool:
    return bool(descriptor and _THERAPY_NAME.search(descriptor))


def render(d: dict) -> str:
    ps = d["partitions"]
    ranked = sorted(ps.items(), key=lambda kv: -kv[1]["ratio"])
    lo, hi = ranked[-1][1]["ratio"], ranked[0][1]["ratio"]
    L = ["# Does the pharmacological:physical claim survive a symmetric class definition?", ""]
    L += ["*Generated by `scripts/atlas_modality_ratio.py` over "
          f"{d['census']:,} census articles. Partitions are the committed input "
          "`analysis/modality-partitions.json`; every member is listed there so "
          "a placement can be disputed.*", ""]

    L += ["## The measured spread", ""]
    L += [f"The held-out column applies the `{PRIMARY_RULE}` rule below, "
          "identically to all five.", ""]
    L += ["| partition | pharmacological | physical | both | **ratio** | "
          "surgical share of physical | ratio, surgery held out |",
          "|---|--:|--:|--:|--:|--:|--:|"]
    for k, c in ranked:
        ss = c.get("surgical_share_of_physical")
        ho = c.get("ratio_surgery_held_out")
        L.append(f"| {k} | {c['pharm']:,} | {c['phys']:,} | {c['both']:,} | "
                 f"**{c['ratio']:.2f}:1** | "
                 f"{100*ss:.1f}% | {ho:.2f}:1 |" if ss is not None and ho
                 else f"| {k} | {c['pharm']:,} | {c['phys']:,} | "
                      f"{c['both']:,} | **{c['ratio']:.2f}:1** | | |")
    L += [""]

    L += _spread_narrative(d, ranked, lo, hi)
    L += _holdout_rule_section(d)
    L += _control_section(d)

    L += [f"| for comparison | | | | |",
          f"|---|--:|--:|--:|--:|",
          f"| manuscript's own method | | | | {MANUSCRIPT_RATIO}:1 |",
          f"| `atlas_landscape.py` on the census | | | | "
          f"{LANDSCAPE_CENSUS_RATIO}:1 |", ""]

    L += ["## What this says", ""]
    L += [f"Every partition gives a ratio between **{lo:.2f}:1 and {hi:.2f}:1**. "
          f"The reported census figure is {LANDSCAPE_CENSUS_RATIO}:1 -- between "
          f"{LANDSCAPE_CENSUS_RATIO/hi:.1f}x and "
          f"{LANDSCAPE_CENSUS_RATIO/lo:.1f}x larger than anything measured "
          f"here.", ""]
    L += ["**The direction survives and the magnitude does not.** Pharmacological "
          "exceeds physical under all five partitions, so the claim's sign is "
          "robust. But the figure that makes it rhetorically powerful -- an "
          "order of magnitude -- is a property of a PHYSICAL class containing "
          "three mechanism tags and excluding radiotherapy, brachytherapy, "
          "hyperthermia, ablation and phototherapy.", ""]
    L += ["That exclusion is not a judgement anyone made. Those modalities have "
          "no mechanism tag in this project, so they could not enter a "
          "tag-based class. The census figure inherits the taxonomy's field of "
          "view rather than measuring the literature (see "
          "`analysis/atlas-taxonomy-reach.md`).", ""]

    L += _comparator_section(d, lo, hi)
    L += _qualifier_section(d)
    L += _what_would_make_this_wrong(d)
    return "\n".join(L) + "\n"


def _what_would_make_this_wrong(d) -> list:
    perm = d.get("permutation") or {}
    ctl = d.get("controls") or {}
    named = {k: v.get("spread") for k, v in ctl.items() if v.get("spread")}
    surg = perm.get("surgical_spread")
    n_below = perm.get("n_draws_at_or_below_surgical")
    # DERIVED. The first draft of this bullet asserted that mass-matched random
    # removals do NOT reproduce the collapse. They do, and the bullet said
    # otherwise because it was written before the control was run.
    if named and surg:
        rests = (f"five rules disagreeing about named descriptor sets all land "
                 f"in the same place, and that "
                 + ", ".join(f"`{k}`" for k in named)
                 + f" -- non-surgical families built the same way -- do not "
                 f"(their spreads are "
                 + ", ".join(f"{v:.2f}x" for v in named.values())
                 + f" against {surg:.2f}x). It does NOT rest on the "
                 f"permutation: {n_below} of {perm['n_draws']} mass-matched "
                 f"random removals reach the same place, so the descriptor "
                 f"identities are not what carries it -- the AMOUNT of "
                 f"surgical mass each partition admits is.")
    else:
        rests = ("five rules disagreeing about named descriptor sets all land "
                 "in the same place.")
    L = ["## What would make this wrong", ""]
    L += ["* If the partitions here are not defensible. They are a committed "
          "member list and every member is disputable; dispute a placement and "
          "re-run. An earlier version of this bullet also claimed each was "
          "\"checked by an independent reviewer for one-sided widening, and "
          "required to reproduce its own count\", and the page's own summary "
          "claimed each was built by a single \"stated\" principle. NO ARTIFACT "
          "IN THIS REPO SUPPORTS ANY OF THAT -- `modality-partitions.json` "
          "carries only the two member lists, there is no reviewer record, no "
          "statement of principle, and the only per-partition counts are this "
          "script's own output, which makes \"reproduce its own count\" "
          "circular. All three clauses are withdrawn.",
          "* If the hold-out rule is doing the work. It is a hand-written list "
          "like the partitions, it leaks and it misses, and the section above "
          "enumerates both rather than asking to be trusted. What the finding "
          "rests on is that " + rests,
          "* Classes are not mutually exclusive, and the `both` column is "
          "reported rather than resolved. Combined-modality treatment is real, "
          "and forcing an article into one class would be the arbitrary step.",
          ""]
    return L


def _spread_narrative(d, ranked, lo, hi) -> list:
    """Derived, so the claim moves with the data.

    An earlier version emitted the headline unconditionally and asserted its
    verbs inside an `if collapse:` block, so a rule that WIDENED the spread
    produced a page still saying it collapsed.
    """
    hos = [c["ratio_surgery_held_out"] for _k, c in ranked
           if c.get("ratio_surgery_held_out")]
    shares = [c.get("surgical_share_of_physical") for _k, c in ranked
              if c.get("surgical_share_of_physical") is not None]
    if len(hos) != len(ranked) or len(shares) != len(ranked):
        return []
    pub, held = hi / lo, max(hos) / min(hos)
    perm = d.get("permutation") or {}
    med = perm.get("median")
    rho = (d["holdouts"].get(PRIMARY_RULE) or {}).get(
        "spearman_share_vs_published_ratio")
    verb = ("collapses" if held < pub * 0.9
            else "widens" if held > pub * 1.1 else "barely moves")
    L = [f"**The spread is one undeclared variable: whether operative surgery "
         f"is a physical modality.** Its share of the physical class runs "
         f"{100*min(shares):.1f}% to {100*max(shares):.1f}% across the five"
         + (f", ranking inversely with the published ratio "
            f"(Spearman {rho:+.2f})" if rho is not None else "") +
         f". Held out under the `{PRIMARY_RULE}` rule applied identically to "
         f"all five, the spread {verb} from **{pub:.2f}x** ({lo:.2f}-{hi:.2f}) "
         f"to **{held:.2f}x** ({min(hos):.2f}-{max(hos):.2f}) -- every "
         f"partition landing near {sum(hos)/len(hos):.1f}:1.", ""]
    if med:
        share = ""
        if med <= held:
            share = (f", at or past where the surgical rule lands, so ALL of "
                     f"the distance from {pub:.2f}x is what removing that "
                     f"much physical mass does regardless of which descriptors "
                     f"carry it")
        elif med < pub:
            frac = (pub - med) / (pub - held) if pub > held else 0.0
            share = (f" -- {100*min(frac, 1.0):.0f}% of the collapse, so most "
                     f"of the distance from {pub:.2f}x is what removing that "
                     f"much physical mass does regardless of which descriptors "
                     f"carry it")
        L += [f"**What that does and does not establish.** Removing a "
              f"mass-matched RANDOM set of physical descriptors takes the "
              f"spread to a median of {med:.2f}x{share}. The control is "
              f"matched to the surgical mass partition by partition, so it "
              f"cannot show surgery is not the driver -- the mass it removes "
              f"IS the surgical mass. What it shows is that the descriptor "
              f"IDENTITIES do not matter once the amount is fixed. The "
              f"section below runs the controls that can speak to "
              f"specificity: named non-surgical families of comparable "
              f"construction.", ""]
    L += ["So the five are not five independent readings whose "
          "disagreement is the finding. They agree about "
          "pharmacological-versus-physical and disagree about one "
          "membership question nobody declared. An earlier version of this "
          "page called that spread \"the deliverable\".", ""]
    return L


def _holdout_rule_section(d) -> list:
    ho = d.get("holdouts") or {}
    if not ho:
        return []
    L = ["## The hold-out rule is itself a judgement", ""]
    L += ["Deciding what counts as operative surgery is as disputable as "
          "deciding whether surgery is physical, so naming one regex would "
          "move the undeclared variable up a level rather than remove it. "
          "Five rules are run. Each disagrees with the others about a named, "
          "listed set of descriptors (`holdouts` in the JSON carries every "
          "member of every set).", ""]
    L += ["| hold-out rule | descriptors held out | spread across the five | range |",
          "|---|--:|--:|---|"]
    for name, leg in sorted(ho.items()):
        rows = leg["partitions"]
        ns = sorted(r["n_descriptors"] for r in rows.values())
        rr = sorted(r["ratio"] for r in rows.values() if r["ratio"])
        sp = leg.get("spread")
        L.append(f"| `{name}` | {ns[0]}-{ns[-1]} | "
                 + (f"**{sp:.2f}x**" if sp else "-") + " | "
                 + (f"{rr[0]:.2f}-{rr[-1]:.2f}:1" if rr else "-") + " |")
    L += [""]
    sps = [leg["spread"] for leg in ho.values() if leg.get("spread")]
    if sps:
        L += [f"The five rules disagree about which descriptors are surgery and "
              f"agree about what happens when they are removed: every one lands "
              f"between **{min(sps):.2f}x and {max(sps):.2f}x**, a disagreement "
              f"smaller than the collapse itself.", ""]

    # what the stem list leaks and misses, derived
    hits = d.get("stem_alternative_hits") or {}
    inert = sorted(k for k, v in hits.items() if not v)
    if inert:
        L += [f"**Inert alternatives.** {len(inert)} of "
              f"{len(STEM_ALTERNATIVES)} stems in the list match nothing "
              f"anywhere in the partition universe (`"
              + "`, `".join(inert) + "`), so the rule is narrower than its "
              "length suggests.", ""]
    L += ["**What it leaks in.** `surg` reaches `"
          + "`, `".join(sorted(ENERGY_CAUGHT_BY_STEM))
          + "`, which are energy modalities rather than operative removal -- "
            "`Radiosurgery` is stereotactic radiotherapy, so holding it out "
            "removes radiotherapy mass from the physical class. `transplant` "
            "reaches cellular and organ transplantation, which is infused or "
            "grafted rather than resected.", ""]
    L += ["**What it leaks out.** No stem reaches `"
          + "`, `".join(sorted(OPERATIVE_MISSED_BY_STEM)) + "`.", ""]
    L += ["The `operative-only` rule excludes the first group and includes the "
          "second; `stem-list` does the opposite; the three intermediates take "
          "one correction each. That is the disagreement the table above "
          "prices.", ""]
    return L


def _control_section(d) -> list:
    perm = d.get("permutation") or {}
    ctl = d.get("controls") or {}
    if not perm.get("n_draws"):
        return []
    L = ["## Is surgery specific, or merely large?", ""]
    surg = perm.get("surgical_spread")
    pub = d.get("published_spread") or 0.0
    L += [f"Two different controls, answering two different questions. The "
          f"first is a permutation: each of {perm['n_draws']} draws (seed "
          f"{perm['seed']}) removes RANDOM physical descriptors from each "
          f"partition until it has held out at least as many articles as the "
          f"surgical rule did THERE, so it is mass-matched partition by "
          f"partition. The second holds out two named non-surgical families "
          f"built the same way the surgical rule is.", ""]
    L += [f"| removal | spread across the five |", "|---|--:|"]
    L += [f"| nothing removed | {pub:.2f}x |"]
    L += [f"| mass-matched random, median of {perm['n_draws']} | "
          f"{perm['median']:.2f}x |",
          f"| mass-matched random, range | {perm['min']:.2f}-{perm['max']:.2f}x |"]
    for name, leg in sorted(ctl.items()):
        sp = leg.get("spread")
        L.append(f"| `{name}` | " + (f"{sp:.2f}x" if sp else "-") + " |")
    if surg:
        L.append(f"| **`{PRIMARY_RULE}`** | **{surg:.2f}x** |")
    L += [""]
    n_below = perm.get("n_draws_at_or_below_surgical", 0)
    if surg:
        many = n_below >= 0.25 * perm["n_draws"]
        L += [f"**{n_below} of {perm['n_draws']}** random draws reach the "
              f"surgical rule's {surg:.2f}x. "
              + (f"So once the AMOUNT of physical mass removed is fixed, which "
                 f"descriptors carry it makes little difference, and the "
                 f"permutation gives no evidence that the surgical "
                 f"descriptors are special as descriptors. This control was "
                 f"matched to the surgical mass, so it could not have: what "
                 f"varies across the five partitions is how much surgery each "
                 f"admits, and holding that amount out is the operation "
                 f"either way."
                 if many else
                 f"So the surgical rule goes beyond what an equally large "
                 f"arbitrary removal achieves. The control is still matched "
                 f"to the surgical mass, so it bounds the descriptor "
                 f"identities rather than the family.") + "", ""]
    named = {k: v.get("spread") for k, v in ctl.items() if v.get("spread")}
    if named:
        worse = {k: v for k, v in named.items() if surg and v > surg}
        L += [f"**The named-family control is the one that speaks to "
              f"specificity, and it does.** Holding out "
              + ", ".join(f"`{k}` gives {v:.2f}x" for k, v in named.items())
              + (f", against the surgical rule's {surg:.2f}x" if surg else "")
              + ". "
              + (f"Every one leaves the five further apart than surgery does"
                 + (f", and {len([v for v in worse.values() if v > pub])} of "
                    f"{len(named)} leave them further apart than removing "
                    f"nothing at all" if any(v > pub for v in worse.values())
                    else "")
                 + ". The partitions do not disagree about radiotherapy or "
                   "about ablation; they disagree about surgery."
                 if len(worse) == len(named) else
                 "At least one non-surgical family does as well, so the "
                 "attribution to surgery is weaker than the spread table "
                 "suggests and should be read as one candidate among "
                 "several."), ""]
    return L


def _comparator_section(d, lo, hi) -> list:
    lc = d.get("landscape_composition") or {}
    if not lc.get("precise_ratio"):
        return []
    hos = sorted(c["ratio_surgery_held_out"] for c in d["partitions"].values()
                 if c.get("ratio_surgery_held_out"))
    L = [f"## What the {LANDSCAPE_CENSUS_RATIO}:1 comparator is made of", ""]
    top = ", ".join(f"`{k}` {v:,}" for k, v in lc["top_two_numerator"])
    L += [f"Its numerator is {lc['numerator']:,} and its denominator "
          f"{lc['denominator']:,}. Two mechanisms supply "
          f"**{100*lc['top_two_share']:.0f}%** of the numerator: {top}. "
          f"`atlas_landscape.py`'s own text calls the largest of them a "
          f"SCOPE ARTIFACT rather than a therapy -- a descriptor carried "
          f"by any paper that MEASURES the process.", ""]
    L += [f"Dropping those from the numerator alone would be the one-sided "
          f"narrowing this page exists to police. Applying that script's "
          f"own `PRECISE` set to BOTH classes -- which also drops "
          f"{', '.join(f'`{x}`' for x in lc['dropped_phys'])} from the "
          f"denominator -- gives **{lc['precise_ratio']:.2f}:1** on "
          f"{lc['precise_numerator']:,} against {lc['precise_denominator']:,}.", ""]

    pr = lc["precise_ratio"]
    if hos:
        inside = hos[0] <= pr <= hos[-1]
        if inside:
            where = (f"which is inside the {hos[0]:.2f}-{hos[-1]:.2f}:1 range "
                     f"the partitions give once surgery is held fixed")
        else:
            side = "above" if pr > hos[-1] else "below"
            edge = hos[-1] if pr > hos[-1] else hos[0]
            where = (f"which is still {pr/edge if side == 'above' else edge/pr:.2f}x "
                     f"{side} the {hos[0]:.2f}-{hos[-1]:.2f}:1 range the "
                     f"partitions give once surgery is held fixed, so the "
                     f"symmetric restriction narrows the gap without closing it")
        L += [f"So the comparator falls from {lc['ratio']:.1f}:1 to "
              f"{pr:.2f}:1 under that restriction, {where}.", ""]

    # PRECISE was not written for this question.
    rest = lc.get("criterion_restored_pharm") or []
    if rest and lc.get("criterion_restored_ratio"):
        names = ", ".join(
            f"`{x}` (`{lc['top_descriptors'].get(x)}`)" for x in rest)
        L += [f"**`PRECISE` is symmetric in rule and was written for a "
              f"different question.** Its stated criterion is that the "
              f"dominant descriptor names a therapy or modality rather than a "
              f"process, a material or a technique, and it was defined for "
              f"`atlas_landscape.py`'s MATURITY comparison, not for this "
              f"volume ratio. Two excluded pharmacological mechanisms satisfy "
              f"the criterion as written: {names}. Restoring them gives "
              f"**{lc['criterion_restored_ratio']:.2f}:1** instead of "
              f"{pr:.2f}:1, so the restricted figure is sensitive to a "
              f"membership judgement made elsewhere for another purpose. Both "
              f"are reported; neither is adopted.", ""]
    L += [f"After the restriction the classes are also no longer what their "
          f"names suggest: pharmacological is `"
          + "`, `".join(lc["precise_pharm"]) + "` and physical is `"
          + "`, `".join(lc["precise_phys"]) + "`.", ""]

    if pr < MANUSCRIPT_RATIO:
        L += [f"**A consequence this page does not resolve.** {pr:.2f}:1 is "
              f"BELOW the manuscript's own {MANUSCRIPT_RATIO}:1, so under a "
              f"symmetric restriction of `atlas_landscape.py`'s own sets the "
              f"census does not understate the manuscript's case -- it "
              f"overstates it. The claim that the manuscript understates "
              f"itself by about 2x is carried in `article/drafts/v1.md`, "
              f"`analysis/atlas-landscape.md`, `analysis/census-findings.md` "
              f"and `CLAUDE.md`. Whether to restrict the comparator at all is "
              f"a choice about what the ratio is FOR, and changing the "
              f"manuscript's framing is an owner decision, so this page states "
              f"the arithmetic and leaves those four sites alone.", ""]
    return L


def _qualifier_section(d) -> list:
    q = d.get("qualifier_recalls") or {}
    rows = q.get("modalities") or {}
    if not rows:
        return []
    L = ["## The qualifier-axis caveat, and why its direction is not established",
         ""]
    L += [f"The ingest reads MeSH DescriptorName and never QualifierName "
          f"(#722), so both classes are understated. An earlier version of "
          f"this page argued from two of that artifact's rows that the bias "
          f"runs AGAINST the finding. THAT INFERENCE WAS INVALID: it compared "
          f"percentage POINTS, and a ratio responds to relative rather than "
          f"additive understatement. Every row, from the raw counts rather "
          f"than the rounded percentages:", ""]
    L += ["| modality | descriptor axis | either axis | descriptor recall |",
          "|---|--:|--:|--:|"]
    for m, c in sorted(rows.items(), key=lambda kv: kv[1]["recall"]):
        L.append(f"| {m} | {c['descriptor']:,} | {c['either']:,} | "
                 f"{c['recall']:.3f} |")
    L += [""]
    dt = rows.get("drug therapy")
    rt = rows.get("radiotherapy")
    sg = rows.get("surgery")
    if dt and rt:
        L += [f"Correcting a drug-therapy numerator and a radiotherapy "
              f"denominator each by its own recall multiplies the ratio by "
              f"{rt['recall']/dt['recall']:.2f} -- UP, not down.", ""]
    if sg and dt:
        L += [f"**But that is the wrong pair for these classes.** #722's own "
              f"headline is that the sharpest case is SURGERY, at recall "
              f"{sg['recall']:.3f} against drug therapy's {dt['recall']:.3f} "
              f"-- and surgery is precisely what dominates the physical class "
              f"here, at up to "
              f"{100*max(c['surgical_share_of_physical'] for c in d['partitions'].values() if c.get('surgical_share_of_physical') is not None):.0f}% "
              f"of it. Correcting by THAT row moves the ratio sharply DOWN.", ""]
    sizes = [c["n_phys_descriptors"] for c in d["partitions"].values()]
    sizes += [c["n_pharm_descriptors"] for c in d["partitions"].values()]
    L += [f"**Neither correction is licensed.** #722 measures proxy sets of "
          f"four to seven descriptors while the classes here hold "
          f"{min(sizes)} to {max(sizes)}, and a broader descriptor set recalls "
          f"more by construction, so those recalls are floors for narrow "
          f"proxies rather than measurements of these classes. Borrowing "
          f"either pair would repeat the category error being retracted. The "
          f"bias is real; its direction is not established here, and the "
          f"measurement that would settle it -- these classes scored on both "
          f"MeSH axes -- is not computable from the committed records, which "
          f"carry the descriptor axis only.", ""]
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan(load_partitions())
        if not d["partitions"] or all(v["phys"] == 0 for v in d["partitions"].values()):
            raise SystemExit(
                "no physical-class articles matched, which is not a finding -- "
                "it is what a descriptor-case mismatch looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
        # Render from the ROUND-TRIPPED artifact, never from the in-memory
        # dict. A full run once wrote an .md the freshness gate could not
        # reproduce from the .json beside it, because JSON normalises key
        # order and tuples and the renderer had seen neither.
        d = json.loads(OUT_JSON.read_text())
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for k, c in sorted(d["partitions"].items(), key=lambda kv: -kv[1]["ratio"]):
        print(f"  {k:24s} {c['ratio']:>6.2f}:1")


if __name__ == "__main__":
    main()
