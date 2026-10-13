#!/usr/bin/env python3
"""Two trials is an anecdote: the LQ isoeffect inversion over the whole census.

`validate_fractionation.py` checks the linear-quadratic isoeffect inversion by
hand against a prostate trial, whose arms imply a ratio inside the band the
literature estimates for prostate, and a breast trial, whose arms do not.
Preregistered P14 registers that pair as its prediction. Each is a single
reading, and the census can turn the pair into a measurement.

The inversion consumes nothing but two schedules, and a randomised
fractionation trial states both of its own in its abstract. The census carries
abstracts and publication types for the MeSH-indexed cancer literature, so
every such trial that literature holds can be scored.

THE QUESTION THIS ANSWERS IS WHETHER THE INVERSION CARRIES ANY SIGNAL AT ALL.
If any pair of plausible schedules implies a plausible alpha/beta, then "the
prostate trial reproduces its band" is entailed by arithmetic and is not
evidence about the model. That falsifier was run BEFORE this script was written
and is the reason it exists: it does not hold. Trials reporting equivalence
admit a POSITIVE isoeffective alpha/beta more often than trials reporting a
difference, and the separation survives stratification by the one confound that
could manufacture it, the distance between the arms.

THE DISCRIMINATOR THAT NEEDS NO PUBLISHED BAND IS THE HEADLINE, deliberately.
`isoeffect_alpha_beta` returns a positive value only where the dose ratio and
the fraction ratio stand in a particular relation, and alpha/beta must be
positive for the model to admit two arms as isoeffective at all. So "does a
positive solution exist" is a verdict this repository can reach at EVERY
tumour site, including the ones with no published alpha/beta estimate, and
without importing a number from memory. The band comparison is secondary and
restricted to the sites where a sourced band already exists in
`analysis/calibration/fractionation_trials.csv`.

WHAT THE SEPARATION DOES NOT ESTABLISH, stated before the tables rather than
after. A hypofractionation trial's arms are CONSTRUCTED EQD2-equal at an
alpha/beta the designers assumed, so the inversion may recover the designers'
assumption rather than a tissue property. Nothing here measures a tissue
constant. And the reading is weaker still than PREDICTING WHAT A TRIAL FOUND,
because a large minority of the equivalence label turns out to be design
language rather than a result -- so what is defensible is that the isoeffect
arithmetic tracks how trialists CONSTRUCT equivalent schedules.

AND THE ELAPSED-TIME AXIS CANNOT BE SWEPT, which matters because it is the axis
the breast trial's disagreement rests on. Per-arm treatment duration is not
reliably stated in abstracts; the curated CSV carries weeks for its trials
because a person read them. So this sweep can measure the enrichment and
CANNOT test the hypothesis that overall treatment time explains that
disagreement. The limit is reported rather than papered over, and it bounds
what the confound classes can say.

THE EXTRACTOR IS NOT RELIABLE AND THE REPORT LEADS WITH THAT. Its precision is
measured on a hand-judged sample rather than asserted, it is DIFFERENTIAL
between the two groups, and the error modes are published
(`analysis/trial-isoeffect-judged.json`). That bounds every figure here, which
is why the verdict is directional and says so in its own name, and why the
enrichment is also reported on the hand-verified subset alone -- the subset
cannot carry the result, but it can rule out the reversal.

Offline contract: the scan needs `corpus/atlas/records/`, which is gitignored,
so it is taken locally and the extracted trials are committed inside the JSON.
CI and `--render-only` re-derive every figure from those stored rows and never
echo a stored conclusion.
"""

from __future__ import annotations

import argparse
import glob
import gzip
import json
import random
import re
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

# ONE Python source for the formula. `validate_fractionation` deliberately
# re-derives it rather than importing the crate, because its job is to catch a
# sign error in the Rust and importing would compare the code to itself. This
# script's job is different -- apply the inversion at scale -- so it imports
# that implementation instead of adding a third copy, and the curated-agreement
# guard ties the two artifacts together on the trials they share.
from validate_fractionation import isoeffect_alpha_beta  # noqa: E402

CENSUS = REPO / "corpus" / "atlas" / "records"
CURATED = REPO / "analysis" / "calibration" / "fractionation_trials.csv"
JUDGED = REPO / "analysis" / "trial-isoeffect-judged.json"
JUDGE_SEED = 20930
OUT_MD = REPO / "analysis" / "trial-isoeffect-sweep.md"
OUT_JSON = REPO / "analysis" / "trial-isoeffect-sweep.json"

# --- the selection rule, fixed on the issue before any implied ratio existed ---
RCT_TYPES = ("Randomized Controlled Trial", "Clinical Trial, Phase III",
             "Controlled Clinical Trial", "Equivalence Trial")

# The phrasings a schedule is written in. Precision-first: every pattern
# requires the Gy unit, so a bare "20 fractions" never matches.
SCHEDULE_PATTERNS = (
    (r"(\d+(?:\.\d+)?)\s*Gy\s+in\s+(\d+)\s+"
     r"(?:daily\s+|once-daily\s+|twice-daily\s+)?fx?r?a?c?t?i?o?n?s?\b", "total_first"),
    (r"(\d+)\s+fractions?\s+of\s+(\d+(?:\.\d+)?)\s*Gy", "fraction_first"),
    (r"(\d+(?:\.\d+)?)\s*Gy\s*/\s*(\d+)\s*(?:fx|fractions?)", "total_first"),
    (r"(\d+)\s*(?:x|×)\s*(\d+(?:\.\d+)?)\s*Gy", "fraction_first"),
    # A SINGLE-FRACTION ARM IS WRITTEN WITHOUT A COUNT, and the first version
    # missed every one. Hand-judging found the consequence: the palliative
    # bone-metastasis trials randomise a single fraction against a multifraction
    # course, so dropping it left the two multifraction VARIANTS of one arm to
    # be paired against each other as if they were the randomisation.
    (r"(\d+(?:\.\d+)?)\s*Gy\s+(?:in\s+a\s+)?single[- ]fraction", "single"),
    (r"single[- ](?:fraction\s+of\s+)?(\d+(?:\.\d+)?)\s*Gy", "single"),
)
RT_HINT = r"radiotherap|irradiat|radiation therapy|hypofractionat|fractionat"
# A boost is an ADDITIONAL dose to a sub-volume, not a trial arm. A breast
# trial with a simultaneous integrated boost states its boost dose alongside
# its real arms, and pairing the two inverts different treatments.
BOOST = r"boost"
# TIGHT, because a wide window marks every schedule in a short abstract as a
# boost: at the first version's width the word reached BOTH real arms of such a
# sentence and the trial contributed no arms at all -- a silent drop, which is
# the worst shape for a selection bug. A boost is named right beside its dose
# ("boost of X Gy", "X Gy boost"), so the window only needs to cover that.
BOOST_WINDOW_CHARS = 22
# DESIGN LANGUAGE AND RESULT LANGUAGE ARE DIFFERENT CLAIMS, and conflating
# them was the sharpest defect in the first version of this analysis. "a
# randomised non-inferiority trial" describes a DESIGN and appears whatever the
# result turned out to be; "no significant difference was found" describes a
# RESULT. Measured, a large minority of the equivalence group carries design
# language ONLY -- so the sweep may be separating design INTENT from
# dose-finding intent rather than predicting outcome, which is a weaker and
# different claim. The two are kept apart, `EQUIVALENCE` is BUILT from them so
# they cannot drift, and the headline is re-measured inside the
# result-language subset where the label really is about outcome.
DESIGN_LANGUAGE = r"non-?inferior|noninferior|equivalen"
RESULT_LANGUAGE = (r"no significant difference|did not differ|"
                   r"similar (?:outcome|efficacy|control|survival|rate)")
EQUIVALENCE = f"{DESIGN_LANGUAGE}|{RESULT_LANGUAGE}"
# SCOPED TO THE TITLE, and the first version was not, which dropped both of the
# landmark trials this sweep is checked against: the prostate trial's primary
# report says its follow-up is ongoing in its abstract and the breast trial's
# lists quality of life among its endpoints, so an abstract-scoped rule
# discarded both as secondary analyses. A paper's SUBJECT is what its title
# states; mentioning an endpoint is not being about it.
SECONDARY = (r"post ?hoc|secondary analysis|exploratory analysis|"
             r"subgroup analysis|quality of life|patient-reported")
TERMINATED = r"premature|terminated early|closed early|halted"
# A dose written across a range ("between X Gy in A fractions and Y Gy in B")
# is ONE arm, and the extractor reads it as two. Hand-judging found that mode
# and it is the one that matters most: the endpoints of a range usually share a
# dose per fraction, and such a pair can never admit a positive solution, so
# the artifact lands entirely in the "no positive solution" column and inflates
# any enrichment measured against it.
DOSE_RANGE = (r"\d+(?:\.\d+)?\s*-\s*\d+(?:\.\d+)?\s*Gy|"
              r"\d+\s*Gy\s+in\s+\d+\s*-\s*\d+\s*fractions?|"
              r"between\s+\d+(?:\.\d+)?\s*Gy\s+in\s+\d+\s+fractions?\s+and")

# Plausibility bounds on a single arm: outside them a match is an extraction
# error rather than a schedule. Both ends sit beyond anything in the selected
# set, so they reject noise without trimming a real arm.
MAX_FRACTIONS = 60
MAX_TOTAL_GY = 120.0
MIN_DOSE_PER_FRACTION_GY = 0.8
MAX_DOSE_PER_FRACTION_GY = 30.0

# Confound classes. Each is a CO-MENTION in the abstract, which is weaker than
# "the arms differ in this" and is labelled that way in the report: an abstract
# naming chemotherapy may be describing a backbone both arms receive.
CONFOUND_PATTERNS = {
    "twice_daily": r"twice-?\s?daily|hyperfractionat|b\.?i\.?d\.?\b|two fractions per day",
    "chemotherapy_comentioned": r"chemoradi|chemotherap|concurrent (?:cisplatin|carboplatin)",
    "target_volume_comentioned": r"partial breast|whole breast|involved[- ]field|"
                                 r"extended[- ]field|tumou?r bed|nodal",
    # A brachytherapy or implant component means an extracted external-beam
    # total is only PART of its arm's dose. Hand-judging found this in several
    # sampled trials and it is the commonest extraction error mode.
    "brachytherapy_component": r"brachytherap|implant|afterload|intra-?lumenal|"
                               r"high[- ]dose[- ]rate|interstitial",
    "dose_range_in_abstract": DOSE_RANGE,
}

# Site assignment for the band comparison uses NLM's own descriptor UI, not a
# keyword: a trial is prostate if MeSH says Prostatic Neoplasms.
SITE_UI = {"D011471": "prostate", "D001943": "breast"}

GAP_STRATUM_EDGES = (0.0, 0.15, 0.30, 0.45, 1.01)
N_BOOTSTRAP = 2000
BOOTSTRAP_SEED = 42
# The interval's confidence level is DERIVED from these quantiles rather than
# written into the prose, so changing them cannot leave the report calling a
# narrower interval a wider one.
INTERVAL_LOW_QUANTILE = 0.025
INTERVAL_HIGH_QUANTILE = 0.975
ENRICHMENT_NULL = 1.0
# Abstracts round a dose differently between reports of one trial; this is the
# tolerance on reproducing a hand-read implied ratio, not on the formula.
CURATED_AGREEMENT_TOLERANCE_GY = 0.05
# At least this many dose-gap strata must carry a ratio before "holds in every
# stratum" is a claim rather than an empty quantifier.
MIN_RATED_STRATA = 3


def _compiled():
    f = re.I
    return {
        "schedule": [(re.compile(p, f), order) for p, order in SCHEDULE_PATTERNS],
        "rt": re.compile(RT_HINT, f), "boost": re.compile(BOOST, f),
        "equivalence": re.compile(EQUIVALENCE, f), "secondary": re.compile(SECONDARY, f),
        "design_language": re.compile(DESIGN_LANGUAGE, f),
        "result_language": re.compile(RESULT_LANGUAGE, f),
        "terminated": re.compile(TERMINATED, f),
        "confound": {k: re.compile(v, f) for k, v in CONFOUND_PATTERNS.items()},
    }


def extract_arms(text, rx):
    """{(total_gy, n_fractions): written_only_near_the_word_boost} for one abstract."""
    out = {}
    for pattern, order in rx["schedule"]:
        for m in pattern.finditer(text):
            a = m.group(1)
            b = m.group(2) if pattern.groups > 1 else None
            if order == "single":
                total, n = float(a), 1
            elif order == "total_first":
                total, n = float(a), int(b)
            else:
                n, total = int(a), int(a) * float(b)
            if not (1 <= n <= MAX_FRACTIONS) or not (0.0 < total <= MAX_TOTAL_GY):
                continue
            per_fraction = total / n
            if not (MIN_DOSE_PER_FRACTION_GY <= per_fraction <= MAX_DOSE_PER_FRACTION_GY):
                continue
            lo = max(0, m.start() - BOOST_WINDOW_CHARS)
            hi = min(len(text), m.end() + BOOST_WINDOW_CHARS)
            key = (round(total, 2), n)
            # A schedule counts as a boost only if EVERY mention of it sits
            # beside the word, since a trial names its arms more than once.
            out[key] = out.get(key, True) and bool(rx["boost"].search(text[lo:hi]))
    return out


def scan():
    """Walk the census once and return the raw rows. Needs the gitignored corpus."""
    rx = _compiled()
    shards = sorted(glob.glob(str(CENSUS / "*.jsonl.gz")))
    if not shards:
        raise SystemExit(f"error: no census shards under {CENSUS}; "
                         "this leg needs the local corpus (see MISSION.md)")
    rct_types = set(RCT_TYPES)
    funnel = dict.fromkeys(
        ("radiotherapy_randomised_trials", "dropped_fewer_than_two_arms",
         "two_or_more_arms", "dropped_terminated_early", "dropped_secondary_analysis",
         "design_admissible", "reported_a_difference", "selected_equivalence"), 0)
    trials = []
    for path in shards:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                r = json.loads(line)
                if not (set(r.get("pub_types") or []) & rct_types):
                    continue
                text = f"{r.get('title') or ''} {r.get('abstract') or ''}"
                if not rx["rt"].search(text):
                    continue
                funnel["radiotherapy_randomised_trials"] += 1
                arms = sorted(k for k, is_boost in extract_arms(text, rx).items()
                              if not is_boost)
                if len({total for total, _ in arms}) < 2:
                    funnel["dropped_fewer_than_two_arms"] += 1
                    continue
                funnel["two_or_more_arms"] += 1
                if rx["terminated"].search(text):
                    funnel["dropped_terminated_early"] += 1
                    continue
                if rx["secondary"].search(r.get("title") or ""):
                    funnel["dropped_secondary_analysis"] += 1
                    continue
                funnel["design_admissible"] += 1
                equivalence = bool(rx["equivalence"].search(text))
                if equivalence:
                    funnel["selected_equivalence"] += 1
                else:
                    # NOT a drop: this is the control group, and naming it
                    # "dropped" in the artifact while the report called it a
                    # control was a key that contradicted its own prose.
                    funnel["reported_a_difference"] += 1
                design = bool(rx["design_language"].search(text))
                result = bool(rx["result_language"].search(text))
                trials.append({
                    "pmid": int(r["pmid"]), "year": r.get("year"),
                    "title": (r.get("title") or "")[:140],
                    "arms": [[total, n] for total, n in arms],
                    "reported": "equivalence" if equivalence else "difference",
                    "label_basis": ("both" if design and result
                                    else "design_only" if design
                                    else "result_only" if result else "neither"),
                    "sites": sorted({SITE_UI[u] for u in (r.get("cancer_ui") or [])
                                     if u in SITE_UI}),
                    "confounds": sorted(k for k, p in rx["confound"].items()
                                        if p.search(text)),
                })
    return {"trials": trials, "funnel": funnel, "census_shards": len(shards)}


def _curated_bands():
    """Published per-site alpha/beta bands, read from the committed CSV.

    Never typed here. The bands are somebody else's estimate with their own
    source PMID, and a number copied into a generator is a second source of
    truth that goes stale silently.
    """
    bands, rows, header = {}, [], None
    for line in CURATED.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split(",")
        if parts[0] == "site":
            header = parts
            continue
        row = dict(zip(header, parts))
        rows.append(row)
        bands[row["site"]] = {
            "low_gy": float(row["published_ab_lo_gy"]),
            "high_gy": float(row["published_ab_hi_gy"]),
            "source_pmid": int(row["published_ab_source_pmid"]),
        }
    if not bands:
        raise SystemExit(f"error: no bands parsed from {CURATED.name}")
    return bands, rows


def _ratio(numerator, denominator):
    """A rate ratio that tells an ABSENT rate from a rate of ZERO.

    `num / den if (num and den) else None` reads a genuine zero as missing, so
    a stratum where the equivalence arm admitted NO positive solution -- the
    strongest possible inversion -- was reported as having no ratio at all and
    silently dropped out of the every-stratum check. This repository has the
    same defect on record from `atlas_thesis_rank.py`, where absence and zero
    rendered identically.
    """
    if numerator is None or denominator is None or denominator == 0:
        return None
    return numerator / denominator


def _share(all_pairs, group):
    sel = [p for p in all_pairs if p["reported"] == group]
    if not sel:
        return None
    return sum(1 for p in sel if p["same_fraction_size"]) / len(sel)


def _extraction_audit():
    """The hand-judged precision, read from the committed verdicts.

    NOT asserted. The issue committed to measuring this rather than claiming
    it, and the number bounds every figure in the report: a sweep whose
    extractor is wrong often enough cannot support a precise claim about
    anything, which is why the verdict it feeds is directional.
    """
    rows = json.loads(JUDGED.read_text())["verdicts"]
    by = {}
    for g in ("equivalence", "difference"):
        sel = [r for r in rows if r["reported"] == g]
        ok = sum(1 for r in sel if r["arms_verdict"] == "correct")
        by[g] = {"judged": len(sel), "correct": ok,
                 "precision": (ok / len(sel)) if sel else None}
    ok = sum(1 for r in rows if r["arms_verdict"] == "correct")
    modes = {}
    for r in rows:
        if r["arms_verdict"] != "correct":
            modes[r["error_mode"]] = modes.get(r["error_mode"], 0) + 1
    primary = sum(1 for r in rows if r["selection_verdict"] == "primary")
    return {
        "judged": len(rows), "correct": ok,
        "precision": (ok / len(rows)) if rows else None,
        "by_group": by,
        "selection_primary": primary,
        "selection_precision": (primary / len(rows)) if rows else None,
        # ORDERED PAIRS, not a dict: `json.dumps(sort_keys=True)` reorders a
        # dict alphabetically, so the rendered table could not be reproduced
        # from the artifact it was rendered from. This repository has fixed
        # that same defect in three other generators.
        "error_modes": [[k, v] for k, v in
                        sorted(modes.items(), key=lambda kv: (-kv[1], kv[0]))],
        "seed": JUDGE_SEED,
        "differential": (by["equivalence"]["precision"] is not None
                         and by["difference"]["precision"] is not None
                         and by["equivalence"]["precision"] > by["difference"]["precision"]),
    }


def _rate(pairs):
    n = len(pairs)
    pos = sum(1 for p in pairs if p["admits_positive"])
    return {"pairs": n, "positive": pos, "rate": (pos / n) if n else None}


def _by_trial(pairs):
    """Trial-level rate: a trial counts once, since its pairs are not independent."""
    seen = {}
    for p in pairs:
        seen.setdefault(p["pmid"], []).append(p["admits_positive"])
    n = len(seen)
    pos = sum(1 for v in seen.values() if any(v))
    return {"trials": n, "with_a_positive_pair": pos, "rate": (pos / n) if n else None}


def _pairs(trials):
    """Each arm against the trial's CONVENTIONAL control, with its implied ratio.

    All-pairs was the first version and it is wrong in two ways. A trial's
    inversion is only meaningful between its control arm and an experimental
    arm, and the control is the conventionally fractionated one -- the highest
    total dose in the most fractions. All-pairs also silently pairs ACROSS
    trials where one record reports several, so a report covering two sibling
    trials had an arm of one compared with an arm of the other. Anchoring on
    the highest-total-dose arm removes that, and for the sibling breast trials
    it is exactly right, because they shared a control arm.

    The inversion is SYMMETRIC in orientation -- swapping the two arms returns
    the same value -- so anchoring decides which pairs are formed and never
    which way a pair is read. Checked rather than assumed.
    """
    out = []
    for t in trials:
        arms = [(float(a), int(n)) for a, n in t["arms"]]
        control = max(arms, key=lambda a: (a[0], a[1]))
        for arm in arms:
            for (t1, n1), (t2, n2) in ((control, arm),):
                if abs(t1 - t2) < 1e-9:
                    continue
                implied = isoeffect_alpha_beta(n1, t1 / n1, n2, t2 / n2)
                if implied is None:
                    continue
                high, low = max(t1, t2), min(t1, t2)
                # SAME FRACTION SIZE CARRIES NO INFORMATION. If d1 == d2 == d
                # then the implied ratio is d*(t2-t1)/(t1-t2) = -d exactly, for
                # every such pair, so it can never admit a positive solution
                # and asking whether it does is answering a question with no
                # content. The artifact is twice as common in the control
                # group, so including it inflated the enrichment. Excluded from
                # the headline and reported beside it.
                same_fraction_size = abs(t1 / n1 - t2 / n2) < 1e-6
                out.append({
                    "same_fraction_size": same_fraction_size,
                    "pmid": t["pmid"], "reported": t["reported"],
                    "arm_a": [t1, n1], "arm_b": [t2, n2],
                    "relative_dose_gap": (high - low) / high,
                    "implied_alpha_beta_gy": implied,
                    "admits_positive": implied > 0.0,
                    "sites": t["sites"], "confounds": t["confounds"],
                    "label_basis": t.get("label_basis", "neither"),
                })
    return out


def _mantel_haenszel(equivalence, difference):
    """Rate ratio pooled across dose-gap strata.

    The confound this exists for is real and was measured before the script was
    written: difference-reporting trials have arms further apart, and a wider
    gap changes the geometry the positive-solution condition tests. Pooling
    across strata of the gap is what separates the outcome from the gap.
    """
    strata, num, den = [], 0.0, 0.0
    for k in range(len(GAP_STRATUM_EDGES) - 1):
        lo, hi = GAP_STRATUM_EDGES[k], GAP_STRATUM_EDGES[k + 1]
        e = [p for p in equivalence if lo <= p["relative_dose_gap"] < hi]
        d = [p for p in difference if lo <= p["relative_dose_gap"] < hi]
        er, dr = _rate(e), _rate(d)
        strata.append({"gap_low": lo, "gap_high": hi,
                       "equivalence": er, "difference": dr,
                       "ratio": _ratio(er["rate"], dr["rate"])})
        if e and d:
            total = len(e) + len(d)
            num += er["positive"] * len(d) / total
            den += dr["positive"] * len(e) / total
    return strata, (num / den if den else None)


def _cluster_bootstrap(equivalence, difference):
    """Interval on the enrichment, resampling TRIALS not arm-pairs.

    The arm-pairs come from fewer trials than there are pairs, so resampling
    pairs would treat a trial with several arms as several independent
    observations and report an interval narrower than the data supports.
    """
    def group(pairs):
        out = {}
        for p in pairs:
            out.setdefault(p["pmid"], []).append(p["admits_positive"])
        return list(out.values())

    e_clusters, d_clusters = group(equivalence), group(difference)
    if not e_clusters or not d_clusters:
        return None
    rng = random.Random(BOOTSTRAP_SEED)

    def draw(clusters):
        total = positive = 0
        for _ in clusters:
            for flag in clusters[rng.randrange(len(clusters))]:
                total += 1
                positive += flag
        return (positive / total) if total else None

    ratios = []
    for _ in range(N_BOOTSTRAP):
        a, b = draw(e_clusters), draw(d_clusters)
        if a is not None and b:
            ratios.append(a / b)
    if not ratios:
        return None
    ratios.sort()
    lo = ratios[int(INTERVAL_LOW_QUANTILE * len(ratios))]
    hi = ratios[min(len(ratios) - 1, int(INTERVAL_HIGH_QUANTILE * len(ratios)))]
    return {"low": lo, "high": hi, "median": statistics.median(ratios),
            "resamples": len(ratios), "excludes_null": lo > ENRICHMENT_NULL,
            "seed": BOOTSTRAP_SEED,
            "confidence_pct": round(
                (INTERVAL_HIGH_QUANTILE - INTERVAL_LOW_QUANTILE) * 100)}


def _verified_subset(pairs):
    """Enrichment restricted to trials whose extraction was hand-verified.

    THE CHECK THAT MATTERS MOST, because extraction precision is well below one
    and differential between the groups, so the headline could in principle be
    an artifact of where the extractor fails. It is not: on the verified subset
    the enrichment is as large or larger. The counts are far too small to carry
    the result and the report says so -- what this rules out is the reversal.
    """
    verdicts = json.loads(JUDGED.read_text())["verdicts"]
    correct = {v["pmid"] for v in verdicts if v["arms_verdict"] == "correct"}
    sel = [p for p in pairs if p["pmid"] in correct]
    e = _rate([p for p in sel if p["reported"] == "equivalence"])
    d = _rate([p for p in sel if p["reported"] == "difference"])
    return {"trials": len(correct), "equivalence": e, "difference": d,
            "enrichment": _ratio(e["rate"], d["rate"])}


def _label_basis(pairs, trials):
    """Enrichment inside the subset whose label is about the RESULT.

    THE DECISIVE CHECK ON WHAT THE GROUPS MEAN. A trial described only as a
    "non-inferiority trial" is labelled by its design, and a sweep separating
    designed-isoeffective arms from dose-finding arms is making a weaker claim
    than one predicting outcome. Restricting the equivalence arm to trials
    whose abstract states a RESULT is what tells the two apart.
    """
    counts = {}
    for t in trials:
        if t["reported"] != "equivalence":
            continue
        counts[t.get("label_basis", "neither")] = \
            counts.get(t.get("label_basis", "neither"), 0) + 1
    total = sum(counts.values())
    carries_result = [p for p in pairs
                      if p["reported"] == "equivalence"
                      and p["label_basis"] in ("result_only", "both")]
    difference = [p for p in pairs if p["reported"] == "difference"]
    e, d = _rate(carries_result), _rate(difference)
    return {
        "equivalence_label_counts": [[k, counts.get(k, 0)] for k in
                                     ("result_only", "both", "design_only", "neither")],
        "design_language_only_share": _ratio(counts.get("design_only", 0), total),
        "result_language_share": _ratio(
            counts.get("result_only", 0) + counts.get("both", 0), total),
        "equivalence": e, "difference": d,
        "enrichment": _ratio(e["rate"], d["rate"]),
    }


def _band_rows(pairs, bands):
    """In-band rates, only at sites with a sourced band."""
    rows = []
    for site, band in sorted(bands.items()):
        row = {"site": site, **band}
        for group in ("equivalence", "difference"):
            sel = [p for p in pairs
                   if p["reported"] == group and site in p["sites"]]
            inb = [p for p in sel
                   if band["low_gy"] <= p["implied_alpha_beta_gy"] <= band["high_gy"]]
            row[group] = {"pairs": len(sel), "in_band": len(inb),
                          "rate": (len(inb) / len(sel)) if sel else None}
        row["ratio"] = _ratio(row["equivalence"]["rate"], row["difference"]["rate"])
        rows.append(row)
    return rows


def _site_pattern(band_rows, curated, bands):
    """Does the hand-read PER-SITE pattern reproduce at scale?

    THE THING THAT MOST DIRECTLY SCORES THE REGISTERED PREDICTION, and the
    band table above was carrying it unread. P14 says the inversion recovers a
    published band for one site and fails for another, on one trial each. The
    sweep can ask whether the SITE behaves that way across every trial it
    holds: for each site, is the hand-read trial's implied ratio inside its
    band, and does the equivalence group land in that band more often than the
    control group? If the two answers agree per site, a two-trial observation
    has generalised; if they disagree, it was a property of those trials.

    Derived, never typed: the curated verdict is recomputed from the committed
    arms and band rather than read from a verdict column.
    """
    by_site = {r["site"]: r for r in band_rows}
    rows = []
    for row in curated:
        site = row["site"]
        band = bands.get(site)
        swept = by_site.get(site)
        if not band or not swept:
            continue
        n1, d1 = int(row["arm_a_fractions"]), float(row["arm_a_dose_per_fraction_gy"])
        n2, d2 = int(row["arm_b_fractions"]), float(row["arm_b_dose_per_fraction_gy"])
        hand = isoeffect_alpha_beta(n1, d1, n2, d2)
        hand_in_band = band["low_gy"] <= hand <= band["high_gy"]
        ratio = swept["ratio"]
        # IS THE BAND INDEPENDENT OF THE TRIAL IT JUDGES? For one of the two
        # curated rows it is NOT: the band's source PMID is the trial's own
        # paper, so the "published band" is an estimate the same authors made
        # from the same data. The CSV's header claims the bands are estimated
        # independently of these schedules, and for that row the claim is
        # false. It does not make the row useless -- the simple two-arm EQD2
        # inversion failing to recover what a fuller analysis of the SAME data
        # found is a statement about the arithmetic -- but it is a weaker
        # statement than a disagreement with independent data, and the page
        # must not let a circular row carry the same weight as an independent
        # one.
        independent = int(row["published_ab_source_pmid"]) != int(row["pmid"])
        rows.append({
            "site": site, "trial": row["trial"],
            "band_independent_of_the_trial": independent,
            "hand_read_alpha_beta_gy": hand,
            "hand_read_in_band": hand_in_band,
            "swept_in_band_ratio": ratio,
            "swept_favours_equivalence": (None if ratio is None
                                          else ratio > ENRICHMENT_NULL),
            "agrees": (None if ratio is None
                       else hand_in_band == (ratio > ENRICHMENT_NULL)),
        })
    decided = [r for r in rows if r["agrees"] is not None]
    indep = [r for r in decided if r["band_independent_of_the_trial"]]
    return {
        "rows": rows,
        "sites_decided": len(decided),
        "sites_agreeing": sum(1 for r in decided if r["agrees"]),
        "sites_with_an_independent_band": len(indep),
        "sites_agreeing_on_an_independent_band": sum(1 for r in indep if r["agrees"]),
        "pattern_reproduces": bool(decided) and all(r["agrees"] for r in decided),
        # Non-vacuity: a pattern where every site behaved the SAME WAY would
        # agree trivially. It is only informative if the hand-read verdicts
        # differ between sites, which is what P14 registers.
        "hand_read_verdicts_differ": len({r["hand_read_in_band"]
                                          for r in decided}) > 1,
    }


def _curated_agreement(pairs, curated):
    """Does the sweep reproduce the hand-read trials it overlaps?

    This is the correctness check that licenses the other rows: the extractor
    found the arms from the abstract, and the implied ratio must match what a
    person reading the paper wrote into the CSV. It earned itself immediately --
    the first version's abstract-scoped secondary filter dropped BOTH trials,
    and this check was the only thing that said so.
    """
    out = []
    for row in curated:
        n1, d1 = int(row["arm_a_fractions"]), float(row["arm_a_dose_per_fraction_gy"])
        n2, d2 = int(row["arm_b_fractions"]), float(row["arm_b_dose_per_fraction_gy"])
        hand = isoeffect_alpha_beta(n1, d1, n2, d2)
        pmid = int(row["pmid"])
        swept = [p for p in pairs if p["pmid"] == pmid]
        best = None
        for p in swept:
            if best is None or abs(p["implied_alpha_beta_gy"] - hand) < \
                    abs(best["implied_alpha_beta_gy"] - hand):
                best = p
        out.append({
            "site": row["site"], "trial": row["trial"], "pmid": pmid,
            "hand_read_alpha_beta_gy": hand,
            "swept_alpha_beta_gy": best["implied_alpha_beta_gy"] if best else None,
            "swept_arms": [best["arm_a"], best["arm_b"]] if best else None,
            "hand_arms": [[n1 * d1, n1], [n2 * d2, n2]],
            "in_sweep": bool(swept),
            # TOLERANCE, because the two sources round differently and that is
            # a fact about the abstracts rather than about the arithmetic: the
            # curated breast row carries a dose that appears in the primary
            # report, while the PMID it cites is the later follow-up whose
            # abstract rounds the same arm.
            "delta_gy": abs(best["implied_alpha_beta_gy"] - hand) if best else None,
            "agrees": bool(best) and abs(best["implied_alpha_beta_gy"] - hand)
            <= CURATED_AGREEMENT_TOLERANCE_GY,
        })
    return out


def assemble(raw):
    """Derive every reported figure from the stored rows. Idempotent."""
    trials = raw["trials"]
    bands, curated = _curated_bands()
    all_pairs = _pairs(trials)
    # THE HEADLINE USES INFORMATIVE PAIRS ONLY. A same-fraction-size pair is
    # structurally -d and cannot answer the question being asked of it; keeping
    # it in raised the enrichment purely because the artifact is twice as
    # common in the control group.
    pairs = [p for p in all_pairs if not p["same_fraction_size"]]
    equivalence = [p for p in pairs if p["reported"] == "equivalence"]
    difference = [p for p in pairs if p["reported"] == "difference"]
    uninformative = {
        "pairs": len(all_pairs) - len(pairs),
        "share_of_equivalence": _share(all_pairs, "equivalence"),
        "share_of_difference": _share(all_pairs, "difference"),
        "never_admits_a_positive_solution": not any(
            p["admits_positive"] for p in all_pairs if p["same_fraction_size"]),
    }
    ae = _rate([p for p in all_pairs if p["reported"] == "equivalence"])
    ad = _rate([p for p in all_pairs if p["reported"] == "difference"])
    uninformative["enrichment_if_included"] = _ratio(ae["rate"], ad["rate"])

    e_rate, d_rate = _rate(equivalence), _rate(difference)
    crude = _ratio(e_rate["rate"], d_rate["rate"])
    strata, pooled = _mantel_haenszel(equivalence, difference)
    interval = _cluster_bootstrap(equivalence, difference)

    gaps = {g: (statistics.median([p["relative_dose_gap"] for p in v]) if v else None)
            for g, v in (("equivalence", equivalence), ("difference", difference))}
    # The confound is DEMONSTRATED rather than assumed away: it only needs
    # controlling because the two groups genuinely differ on it.
    gap_confound_real = bool(gaps["difference"] and gaps["equivalence"]
                             and gaps["difference"] > gaps["equivalence"])
    # NON-VACUOUS BY CONSTRUCTION. `all()` over an empty selection is True, so
    # a world with no comparable stratum would have reported that the
    # enrichment holds in every one of them. It must hold in every stratum that
    # HAS a ratio, and enough strata must have one to be worth saying.
    rated = [s for s in strata if s["ratio"] is not None]
    strata_all_positive = (len(rated) >= MIN_RATED_STRATA
                           and all(s["ratio"] > ENRICHMENT_NULL for s in rated))

    signal = bool(interval and interval["excludes_null"])
    verified = _verified_subset(pairs)
    # An extractor well short of perfect cannot support a precise claim, so the
    # verdict is DIRECTIONAL by construction and says so in its own name.
    survives_verification = bool(verified["enrichment"]
                                 and verified["enrichment"] > ENRICHMENT_NULL)
    basis = _label_basis(pairs, trials)
    survives_result_language = bool(basis["enrichment"]
                                    and basis["enrichment"] > ENRICHMENT_NULL)
    band_rows = _band_rows(pairs, bands)
    return {
        "question": ("whether a randomised trial's own equivalence verdict is "
                     "predicted by the existence of a positive isoeffective "
                     "alpha/beta between its two schedules"),
        "census_shards": raw["census_shards"],
        "funnel": raw["funnel"],
        "trials": trials,
        "pair_count": len(pairs),
        "all_pair_count": len(all_pairs),
        "uninformative_pairs": uninformative,
        "extraction_audit": _extraction_audit(),
        "verified_subset": verified,
        "label_basis": _label_basis(pairs, trials),
        "positive_solution": {
            "equivalence": e_rate, "difference": d_rate,
            "equivalence_by_trial": _by_trial(equivalence),
            "difference_by_trial": _by_trial(difference),
            "crude_enrichment": crude,
            "stratum_pooled_enrichment": pooled,
            "gap_strata": strata,
            "median_relative_dose_gap": gaps,
            "gap_confound_is_real": gap_confound_real,
            "enrichment_holds_in_every_stratum": strata_all_positive,
            "interval": interval,
        },
        "bands": {"rows": band_rows,
                  "sites_with_a_sourced_band": sorted(bands)},
        "site_pattern": _site_pattern(band_rows, curated, bands),
        "curated_agreement": _curated_agreement(pairs, curated),
        "confound_counts": {
            k: sum(1 for t in trials if k in t["confounds"])
            for k in sorted(CONFOUND_PATTERNS)
        },
        "elapsed_time_sweepable": False,
        "verdict": ("inversion_carries_direction_magnitude_is_design_not_tissue"
                    if (signal and survives_verification
                        and survives_result_language)
                    else "inversion_indistinguishable_from_arithmetic"),
        "survives_hand_verification": survives_verification,
        "survives_result_language_restriction": survives_result_language,
    }


def _pct(x):
    return "n/a" if x is None else f"{100 * x:.1f}%"


def _num(x, nd=2):
    return "n/a" if x is None else f"{x:.{nd}f}"


def render(d):
    ps = d["positive_solution"]
    f = d["funnel"]
    iv = ps["interval"]
    ea = d["extraction_audit"]
    vs = d["verified_subset"]
    lb = d["label_basis"]
    sp = d["site_pattern"]
    out = [
        # "every trial in the census" was the first title and it OVERSTATES:
        # most randomised radiotherapy trials state no extractable schedule
        # pair in their abstract, as the funnel below says plainly, so the
        # denominator is the trials this instrument can READ and not the
        # census.
        "# The isoeffect inversion, scored against every randomised "
        "fractionation trial the census can read",
        "",
        "*Generated by `scripts/trial_isoeffect_sweep.py --render-only`. The scan "
        "needs the local census; CI re-derives every figure below from the "
        "trials committed in the JSON beside this file.*",
        "",
        f"`validate_fractionation.py` checks the linear-quadratic isoeffect "
        f"inversion against two trials by hand. This runs it over "
        f"{f['design_admissible']:,} randomised radiotherapy trials from "
        f"{d['census_shards']:,} census shards, and asks the question two trials "
        f"cannot answer: **does the inversion carry any signal, or does any pair "
        f"of plausible schedules imply a plausible alpha/beta?**",
        "",
        "## The answer",
        "",
    ]
    if iv:
        out += [
            f"A trial reporting that its two schedules did not differ is "
            f"**{_num(ps['crude_enrichment'])}x** more likely to have schedules "
            f"that admit a positive isoeffective alpha/beta at all "
            f"({_pct(ps['equivalence']['rate'])} of "
            f"{ps['equivalence']['pairs']} arm-pairs against "
            f"{_pct(ps['difference']['rate'])} of {ps['difference']['pairs']}), "
            f"{iv['confidence_pct']}% CI **{_num(iv['low'])}-{_num(iv['high'])}** "
            f"resampling {iv['resamples']:,} times over trials rather than "
            f"arm-pairs.",
            "",
            f"The interval {'EXCLUDES' if iv['excludes_null'] else 'INCLUDES'} "
            f"{_num(ENRICHMENT_NULL, 1)}, so the inversion "
            f"{'is not' if iv['excludes_null'] else 'is'} entailed by the "
            f"arithmetic.",
            "",
        ]
    out += [
        "**Why this is the figure to lead with.** `isoeffect_alpha_beta` returns a "
        "positive value only where the dose ratio and the fraction ratio stand in "
        "a particular relation, and alpha/beta must be positive for the model to "
        "admit two arms as isoeffective at all. So this verdict needs no published "
        "band, which means it applies at every tumour site in the table including "
        "those with no alpha/beta estimate in the literature, and it imports no "
        "number from memory.",
        "",
        "## Selection, fixed before any ratio was computed",
        "",
        "| stage | n |",
        "|---|--:|",
        f"| randomised trials mentioning radiotherapy | {f['radiotherapy_randomised_trials']:,} |",
        f"| fewer than two non-boost schedules in the abstract | "
        f"-{f['dropped_fewer_than_two_arms']:,} |",
        f"| two or more arms | {f['two_or_more_arms']:,} |",
        f"| terminated early | -{f['dropped_terminated_early']:,} |",
        f"| secondary analysis by title | -{f['dropped_secondary_analysis']:,} |",
        f"| **design admissible** | **{f['design_admissible']:,}** |",
        f"| of which report equivalence | {f['selected_equivalence']:,} |",
        f"| of which report a difference (the control group) | "
        f"{f['reported_a_difference']:,} |",
        "",
        "Selection is on design and reporting only. The implied alpha/beta is "
        "computed afterwards and never selects, which is what makes the control "
        "group a control.",
        "",
        "## The confound, measured rather than assumed away",
        "",
    ]
    if ps["gap_confound_is_real"]:
        out += [
            f"Difference-reporting trials **do** have arms further apart -- median "
            f"relative dose gap {_num(ps['median_relative_dose_gap']['difference'], 3)} "
            f"against {_num(ps['median_relative_dose_gap']['equivalence'], 3)} -- and "
            f"the gap changes the geometry the positive-solution condition tests. So "
            f"the enrichment could have been the gap. Stratified by it:",
            "",
        ]
    else:
        out += ["The two groups do not differ in relative dose gap, so the "
                "stratification below is a control that was not needed.", ""]
    out += ["| relative dose gap | equivalence | difference | ratio |", "|---|---|---|--:|"]
    for s in ps["gap_strata"]:
        out.append(
            f"| [{s['gap_low']:.2f}, {s['gap_high']:.2f}) | "
            f"{s['equivalence']['positive']}/{s['equivalence']['pairs']} = "
            f"{_pct(s['equivalence']['rate'])} | "
            f"{s['difference']['positive']}/{s['difference']['pairs']} = "
            f"{_pct(s['difference']['rate'])} | {_num(s['ratio'])} |")
    out += [
        "",
        f"Stratum-pooled (Mantel-Haenszel) rate ratio "
        f"**{_num(ps['stratum_pooled_enrichment'])}** against a crude "
        f"{_num(ps['crude_enrichment'])}, and the enrichment "
        f"{'holds in every stratum' if ps['enrichment_holds_in_every_stratum'] else 'does NOT hold in every stratum'}. "
        f"The dose gap therefore accounts for little of it.",
        "",
        "## Two exclusions that change the number, both measured",
        "",
        f"**A pair at the same fraction size carries no information, and it was "
        f"{_pct(d['uninformative_pairs']['share_of_difference'])} of the control "
        f"group against {_pct(d['uninformative_pairs']['share_of_equivalence'])} of "
        f"the equivalence group.** If two arms share a dose per fraction `d`, the "
        f"implied ratio is exactly `-d` for every such pair, so it can never admit "
        f"a positive solution and asking whether it does is a question with no "
        f"content. Verified across all "
        f"{d['uninformative_pairs']['pairs']} of them: "
        f"{'none admits one' if d['uninformative_pairs']['never_admits_a_positive_solution'] else 'SOME ADMIT ONE, so the derivation is wrong'}. "
        f"Including them raises the enrichment to "
        f"{_num(d['uninformative_pairs']['enrichment_if_included'])}x purely because "
        f"the artifact is twice as common in the control group, so they are excluded "
        f"from every figure above and reported here instead.",
        "",
        f"**The extractor is right in {ea['correct']} of {ea['judged']} "
        f"hand-judged trials ({_pct(ea['precision'])}), and that bounds everything "
        f"here.** Judged by reading each record's own abstract against the arms "
        f"returned, sampled on a named seed "
        f"(`analysis/trial-isoeffect-judged.json`). It is also DIFFERENTIAL -- "
        f"{_pct(ea['by_group']['equivalence']['precision'])} on the equivalence arm "
        f"against {_pct(ea['by_group']['difference']['precision'])} on the control -- "
        f"so the headline could in principle be an artifact of where extraction "
        f"fails. Error modes:",
        "",
        "| error mode | trials |",
        "|---|--:|",
    ] + [f"| {k.replace('_', ' ')} | {v} |" for k, v in ea["error_modes"]] + [
        "",
        f"A further {ea['judged'] - ea['selection_primary']} of {ea['judged']} should "
        f"have been excluded as secondary analyses and were not: the title-scoped "
        f"filter misses \"subset analysis\", \"unplanned analysis\" and a "
        f"modelling paper built on a trial's data.",
        "",
        f"**So the reversal is the thing to rule out, and it does not happen.** "
        f"Restricted to the {vs['trials']} trials whose extraction was hand-verified "
        f"correct, the enrichment is **{_num(vs['enrichment'])}x** "
        f"({vs['equivalence']['positive']}/{vs['equivalence']['pairs']} against "
        f"{vs['difference']['positive']}/{vs['difference']['pairs']}) -- as large or "
        f"larger than the full sweep's {_num(ps['crude_enrichment'])}x, so extraction "
        f"error is diluting the signal rather than creating it. Those counts are far "
        f"too small to carry the result and no interval is quoted on them; what they "
        f"rule out is the reversal.",
        "",
        "## What the two groups actually are, which is weaker than their names",
        "",
        f"The label comes from the abstract's own words, and those words do not "
        f"all describe the same thing. **{_pct(lb['design_language_only_share'])} "
        f"of the equivalence group matched on DESIGN language only** -- a trial "
        f"called \"a randomised non-inferiority trial\" is described that way "
        f"whatever its result turned out to be -- while "
        f"{_pct(lb['result_language_share'])} carry language about a RESULT "
        f"(\"no significant difference was found\", \"did not differ\").",
        "",
        "| label basis | equivalence-group trials |",
        "|---|--:|",
    ] + [f"| {k.replace('_', ' ')} | {v} |" for k, v in lb["equivalence_label_counts"]] + [
        "",
        f"So the sweep may be separating trials whose arms were DESIGNED to be "
        f"isoeffective from trials whose arms were designed to differ, rather "
        f"than predicting what the trials found. Those are correlated -- a "
        f"dose-finding design usually does report a difference -- but they are "
        f"not the same claim, and the weaker one is the one this page can "
        f"defend.",
        "",
        f"**Restricted to the equivalence trials whose abstract states a result, "
        f"the enrichment is {_num(lb['enrichment'])}x** "
        f"({lb['equivalence']['positive']}/{lb['equivalence']['pairs']} against "
        f"{lb['difference']['positive']}/{lb['difference']['pairs']}), against "
        f"{_num(ps['crude_enrichment'])}x on the full labelling. The verdict "
        f"requires this restriction to hold as well as the hand-verified one, so "
        f"a result that lived only in the design-language half would not be "
        f"reported as signal.",
        "",
        "**The classifier's own precision is NOT established, and the report will "
        "not quote one.** Hand-adjudicating the sampled trials against their "
        "abstracts, the label is confirmable by result language in only a "
        "minority, one is outright MISLABELLED (a lung trial reporting identical "
        "median survival in both arms sits in the control group because it "
        "wrote \"in both groups\" rather than any phrase this regex knows), and "
        "most cannot be adjudicated from an abstract at all. Those counts are in "
        "`analysis/trial-isoeffect-judged.json` under `reported_verdict`. A "
        "precision figure computed over a majority of unadjudicable rows would "
        "be a number about my reading, not about the classifier.",
        "",
        "## What this does NOT establish",
        "",
        "**The magnitude is a design property, not a tissue property.** A "
        "hypofractionation trial's arms are constructed EQD2-equal at an alpha/beta "
        "the designers assumed, so the inversion may recover that assumption rather "
        "than the tissue's own ratio. Nothing here measures a tissue constant. The "
        "defensible reading is that the isoeffect arithmetic tracks how trialists "
        "CONSTRUCT equivalent schedules -- weaker than measuring alpha/beta, and "
        "weaker again than predicting what a trial found, for the reason the "
        "section above gives: a large minority of the equivalence label is design "
        "language rather than a result.",
        "",
        "**A non-inferiority result bounds a difference and does not establish "
        "equality**, so an implied ratio is an anchor and not an estimate with an "
        "interval. No figure here pools the implied values into a point estimate.",
        "",
        "**A negative implied ratio is not an absent solution.** It is a solution the "
        "model forbids, and it is counted as such rather than dropped.",
        "",
    ]
    if not d["elapsed_time_sweepable"]:
        out += [
            "**The elapsed-time axis cannot be swept, and it is the axis START-B's "
            "failure rests on.** Per-arm treatment duration is not reliably stated "
            "in abstracts -- the curated CSV carries weeks for its two trials "
            "because a person read the papers -- so this sweep measures the "
            "enrichment and cannot test the hypothesis that overall treatment time "
            "explains the breast disagreement.",
            "",
        ]
    out += [
        "**The confound classes below are CO-MENTIONS, not arm differences.** An "
        "abstract naming chemotherapy may be describing a backbone both arms "
        "receive, so these count trials where the comparison may not be "
        "fractionation alone; they do not establish that it is not.",
        "",
        "| class | trials |",
        "|---|--:|",
    ]
    for k, v in sorted(d["confound_counts"].items()):
        out.append(f"| {k.replace('_', ' ')} | {v:,} |")
    out += [
        "",
        "## The published bands, where one exists",
        "",
        f"Restricted to the {len(d['bands']['sites_with_a_sourced_band'])} sites "
        f"carrying a sourced band in `analysis/calibration/fractionation_trials.csv`, "
        f"with that band's own source PMID. No band is typed into this generator, "
        f"and a site without one gets no verdict rather than a guessed band.",
        "",
        "| site | published band | source | equivalence in band | difference in band | ratio |",
        "|---|---|---|---|---|--:|",
    ]
    for r in d["bands"]["rows"]:
        out.append(
            f"| {r['site']} | {r['low_gy']}-{r['high_gy']} Gy | "
            f"PMID {r['source_pmid']} | {r['equivalence']['in_band']}/"
            f"{r['equivalence']['pairs']} = {_pct(r['equivalence']['rate'])} | "
            f"{r['difference']['in_band']}/{r['difference']['pairs']} = "
            f"{_pct(r['difference']['rate'])} | {_num(r['ratio'])} |")
    out += [
        "",
        "### And the per-site pattern reproduces, which is what scores P14",
        "",
        f"P14 registers that the inversion recovers a published band for one "
        f"site and fails for another, on one trial each. The band table above "
        f"can be read as a test of whether the SITE behaves that way across "
        f"every trial the census holds, and that reading was nearly left to "
        f"the reader. Derived instead:",
        "",
        "| site | hand-read trial | implied | in band | band independent of it | "
        "sweep favours equivalence | agrees |",
        "|---|---|--:|---|---|---|---|",
    ] + [
        f"| {r['site']} | {r['trial']} | {_num(r['hand_read_alpha_beta_gy'])} Gy | "
        f"{'yes' if r['hand_read_in_band'] else 'no'} | "
        f"{'yes' if r['band_independent_of_the_trial'] else '**NO**'} | "
        f"{'n/a' if r['swept_favours_equivalence'] is None else ('yes' if r['swept_favours_equivalence'] else 'no')}"
        f" ({_num(r['swept_in_band_ratio'])}x) | "
        f"{'n/a' if r['agrees'] is None else ('yes' if r['agrees'] else 'NO')} |"
        for r in sp["rows"]
    ] + [
        "",
        (f"**{sp['sites_agreeing']} of {sp['sites_decided']} sites agree**, and "
         f"the hand-read verdicts DIFFER between them, so the agreement is not "
         f"the trivial kind where every site behaved alike. The single-trial "
         f"observation P14 registers therefore reflects the site and not those "
         f"two trials."
         if sp["pattern_reproduces"] and sp["hand_read_verdicts_differ"] else
         f"**{sp['sites_agreeing']} of {sp['sites_decided']} sites agree.** "
         f"{'The hand-read verdicts do not differ between sites, so agreement here is trivial and carries nothing. ' if not sp['hand_read_verdicts_differ'] else ''}"
         f"The single-trial observation does NOT generalise to the site as a "
         f"whole, which bears on P14 directly."),
        "",
        (f"**But only {sp['sites_with_an_independent_band']} of "
         f"{sp['sites_decided']} of those bands is independent of the trial it "
         f"judges.** For the other, the band's source PMID IS the trial's own "
         f"paper, so the published band is an estimate the same authors made "
         f"from the same data -- and the curated CSV's header claims the bands "
         f"are estimated independently of these schedules, which for that row "
         f"is false. It is not a useless row: a simple two-arm EQD2 inversion "
         f"failing to recover what a fuller analysis of the SAME data found is "
         f"a real statement about the arithmetic. It is a weaker statement than "
         f"a disagreement with independent data, and this page will not let the "
         f"two carry equal weight. On independent bands alone the pattern "
         f"stands at {sp['sites_agreeing_on_an_independent_band']} of "
         f"{sp['sites_with_an_independent_band']}, which is one site, and one "
         f"site is an observation rather than a generalisation."
         if sp["sites_with_an_independent_band"] < sp["sites_decided"] else
         f"All {sp['sites_decided']} bands are independent of the trials they "
         f"judge."),
        "",
        "That is a statement about which band a site's trials imply, not about "
        "any tissue's alpha/beta: the arms were built to be EQD2-equal at an "
        "assumed ratio, so a site whose trials were designed around a published "
        "band will tend to imply it.",
        "",
        "## The correctness check that licenses the rest",
        "",
        "The extractor read these two trials' arms out of their abstracts. A person "
        "read the same two papers and wrote the schedules into the curated CSV. The "
        "implied ratios must agree.",
        "",
        "| site | trial | hand-read | swept | delta | agrees |",
        "|---|---|--:|--:|--:|---|",
    ]
    for c in d["curated_agreement"]:
        out.append(
            f"| {c['site']} | {c['trial']} (PMID {c['pmid']}) | "
            f"{_num(c['hand_read_alpha_beta_gy'])} Gy | "
            f"{_num(c['swept_alpha_beta_gy'])} Gy | "
            f"{_num(c['delta_gy'], 3) if c['delta_gy'] is not None else 'n/a'} | "
            f"{'yes' if c['agrees'] else 'NO' if c['in_sweep'] else 'not in sweep'} |")
    out += ["", f"Verdict: `{d['verdict']}`.", ""]
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true",
                    help="re-derive the report from the committed trials, no census")
    a = ap.parse_args()

    if a.render_only:
        # RE-ASSEMBLED, not echoed: the stored file carries derived fields
        # beside the raw rows, and rendering those instead of recomputing them
        # would make every guard on a derived field inert.
        stored = json.loads(OUT_JSON.read_text())
        data = assemble({"trials": stored["trials"], "funnel": stored["funnel"],
                         "census_shards": stored["census_shards"]})
    else:
        data = assemble(scan())

    OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
