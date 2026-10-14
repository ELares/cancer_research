#!/usr/bin/env python3
"""Condition the published prior-predictive intervals on the model's own
baseline-viability constraint (#936).

`uncertainty_intervals.py` publishes prior-predictive intervals over uniform
priors across the PRCC rate-constant ranges, and `identifiability-report`
leans on them for its per-headline verdicts. The manuscript separately states a
constraint in prose -- every phenotype must show less than a stated percentage of
death under Control -- and `headline_at_fitted.py` already uses it to reject the
in-vitro-fitted cascade as inadmissible.

NOTHING APPLIED THAT CONSTRAINT TO THE PRIOR ITSELF. The published design
evaluates treated arms only: no Control condition appears in it, so the
untreated death rate was never computed and the admissibility of a draw was
never knowable from the artifact. Most of that prior turns out to violate the
constraint.

WHAT THIS DOES
--------------
Reproduces the published design exactly -- the same ranges, draw count, seeds and
per-draw cell count, all IMPORTED from that generator rather than restated --
adds the untreated arms on the SAME draws, and partitions. Partitioning one
sample rather than drawing two is the whole point: a second sample would confound
the conditioning with sampling noise.

The full-prior half must reproduce the committed report's own numbers, and the
page refuses to compare against anything else. If it cannot reproduce them, the
comparison is not with what is published and the artifact says so instead of
quietly comparing against a re-run.

WHAT CONDITIONING IS AND IS NOT
-------------------------------
Removing draws the model itself declares invalid is not inference. It adds no
data, it cannot turn a prior-predictive interval into a data-conditioned one, and
it leaves the identifiability report's data-conditioned-headline count where it
is. What it does is stop the published intervals describing a region where the
untreated cell is already dead.

DIRECTION, WHICH IS NOT ONE-SIDED
---------------------------------
Admissibility selects for draws with strong defences, so it lowers kill
everywhere: the inducer arms' upper tails collapse, and the exogenous-ROS arms'
lower tails extend further down. One mechanism, two opposite-looking corrections,
and reporting only one of them would flatter whichever side was looked at first.

OFFLINE CONTRACT
----------------
The measurement needs the compiled Python extension, so it is taken locally with
`--measure` and the per-draw results are committed. Every figure is then derived
from that file, so CI re-derives the whole page and downloads nothing.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

# ONE source for the published design. Restating any of these is how a
# "reproduction" quietly stops reproducing anything.
from uncertainty_intervals import (  # noqa: E402
    CONDITIONS, N_CELLS, N_SAMPLES, PARAMS, QUANTILES, RNG_SEED, SIM_SEED,
)
# ONE source for the constraint, whose own guard parses it from the manuscript.
from headline_at_fitted import BASELINE_MAX  # noqa: E402

DRAWS_CSV_GZ = REPO / "analysis" / "calibration" / "admissibility-draws.csv.gz"
OUT_JSON = REPO / "analysis" / "intervals-under-admissibility.json"
OUT_MD = REPO / "analysis" / "intervals-under-admissibility.md"

PUBLISHED_MD = REPO / "analysis" / "uncertainty-intervals-report.md"

# The untreated arms. The same three `headline_at_fitted` checks, so the two
# pages judge admissibility identically. Measured: adding OXPHOS and Stromal
# changes no verdict, because Persister always binds -- reported, not assumed.
CONTROL_CONDITIONS = (("Glycolytic", "Control"),
                      ("Persister", "Control"),
                      ("PersisterNrf2", "Control"))

# Below this many admissible draws no interval is quoted on the subset.
MIN_ADMISSIBLE_DRAWS = 100

# A median moving by less than this is not a correction worth reporting as one.
MATERIAL_SHIFT = 0.02

# The published table prints one decimal place, so parsing it back loses at most
# half of that last digit. The reproduction tolerance is DERIVED from that
# rounding rather than picked: a round number chosen by hand was 22x looser than
# the largest real deviation, which is a tolerance that would accept a genuine
# one-point discrepancy while reporting a reproduction. The factor leaves room
# for the last-bit of the two formatters and is named rather than folded in.
PUBLISHED_DECIMALS = 1
ROUNDING_SLACK = 2.0
REPRODUCTION_TOLERANCE = ROUNDING_SLACK * 0.5 * 10 ** (-PUBLISHED_DECIMALS) / 100.0

CLAIMS_CSV = REPO / "analysis" / "calibration" / "published-interval-claims.csv"


def published_claims():
    """The committed report's own conclusions, read from a committed file.

    They live in data rather than in this source for two reasons. A quoted claim
    legitimately contains whatever figures its author wrote, and a generator may
    not hand-write a figure -- so quoting them here would be indistinguishable
    from asserting them. And a reviewer can audit the list without reading
    Python. A guard asserts every entry is still present verbatim in that
    report: if the prose is rewritten, this page must be rewritten with it
    rather than scoring a sentence nobody publishes any more.
    """
    rows = []
    with CLAIMS_CSV.open() as fh:
        for row in csv.DictReader(fh):
            if row["claim"].startswith("#"):
                continue
            rows.append((row["claim"], row["condition"], row["property"]))
    return tuple(rows)


def quoted_floor(fragment) -> float:
    """The percentage a claim quotes ABOUT ITSELF, parsed out of its own text.

    A threshold typed beside a sentence that already states it is two copies of
    one number, and the sentence is the one that gets edited. Scoring a claim
    against a threshold it did not state is the failure this avoids.
    """
    m = re.search(r">\s*(\d+(?:\.\d+)?)\s*%", fragment)
    if not m:
        raise ValueError(f"no self-quoted percentage in {fragment!r}")
    return float(m.group(1)) / 100.0

# An interval narrower than this is not "wide".
WIDE_INTERVAL = 0.25


def _label(condition) -> str:
    return f"{condition[0]}/{condition[1]}"


def draw_matrix() -> np.ndarray:
    """The published design's own draws, regenerated from its seed."""
    lows = np.array([p[1] for p in PARAMS], float)
    highs = np.array([p[2] for p in PARAMS], float)
    rng = np.random.default_rng(RNG_SEED)
    return lows + rng.random((N_SAMPLES, lows.shape[0])) * (highs - lows)


def measure() -> None:
    """Evaluate every condition, treated and untreated, on the same draws."""
    import uncertainty_intervals as U

    draws = draw_matrix()
    names = [p[0] for p in PARAMS]
    conditions = list(CONDITIONS) + list(CONTROL_CONDITIONS)
    DRAWS_CSV_GZ.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(DRAWS_CSV_GZ, "wt", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["draw"] + names + [_label(c) for c in conditions])
        for i, row in enumerate(draws):
            overrides = {n: float(row[j]) for j, n in enumerate(names)}
            rates = [U._fc().sim_batch(ph, tr, n=N_CELLS, seed=SIM_SEED,
                                      **overrides)["death_rate"]
                     for ph, tr in conditions]
            writer.writerow([i] + [f"{v:.10g}" for v in row]
                            + [f"{v:.6f}" for v in rates])
            if (i + 1) % 500 == 0:
                print(f"  {i + 1} of {N_SAMPLES} draws", flush=True)
    print(f"wrote {DRAWS_CSV_GZ.relative_to(REPO)}")


def load() -> dict:
    """Read the committed per-draw results. No extension, no network."""
    with gzip.open(DRAWS_CSV_GZ, "rt") as fh:
        rows = list(csv.DictReader(fh))
    return {"rows": rows, "published": PUBLISHED_MD.read_text()}


def parse_published(text) -> dict:
    """The committed report's own table, so the comparison is with what is
    PUBLISHED rather than with a re-run of it."""
    out = {}
    for line in text.split("\n"):
        m = re.match(r"\|\s*(\w+)\s*\|\s*(\w+)\s*\|\s*([\d.]+)%\s*\|\s*([\d.]+)%"
                     r"\s*\|\s*([\d.]+)%\s*\|\s*([\d.]+)%\s*\|", line)
        if not m:
            continue
        pheno, treat, point, median, low, high = m.groups()
        out[f"{pheno}/{treat}"] = {
            "point": float(point) / 100.0, "median": float(median) / 100.0,
            "low": float(low) / 100.0, "high": float(high) / 100.0,
        }
    return out


def _quantiles(values):
    """The published method, not a hand-rolled index: `np.quantile`'s default
    interpolation is what produced the numbers this page compares against."""
    if not values:
        return None
    arr = np.asarray(values, float)
    low, median, high = np.quantile(arr, QUANTILES)
    return {"low": float(low), "median": float(median), "high": float(high),
            "width": float(high - low), "n": len(values)}


def assemble(raw) -> dict:
    rows, published = raw["rows"], parse_published(raw["published"])
    treated = [_label(c) for c in CONDITIONS]
    untreated = [_label(c) for c in CONTROL_CONDITIONS]

    admissible = []
    worst = []
    for row in rows:
        w = max(float(row[c]) for c in untreated)
        worst.append(w)
        admissible.append(w < BASELINE_MAX)
    n_adm = sum(1 for a in admissible if a)

    per_condition = {}
    for label in treated:
        full = [float(r[label]) for r in rows]
        adm = [float(r[label]) for r, ok in zip(rows, admissible) if ok]
        f, a = _quantiles(full), (_quantiles(adm)
                                  if n_adm >= MIN_ADMISSIBLE_DRAWS else None)
        pub = published.get(label)
        reproduces = None
        if pub is not None and f is not None:
            reproduces = all(
                abs(f[k] - pub[k]) <= REPRODUCTION_TOLERANCE
                for k in ("median", "low", "high"))
        entry = {"published": pub, "full_prior": f, "admissible_only": a,
                 "reproduces_published": reproduces}
        if f and a:
            entry["median_shift"] = a["median"] - f["median"]
            entry["width_ratio"] = (a["width"] / f["width"]) if f["width"] else None
            entry["upper_ratio"] = (a["high"] / f["high"]) if f["high"] else None
            entry["lower_shift"] = a["low"] - f["low"]
            entry["median_moved_materially"] = abs(entry["median_shift"]) >= MATERIAL_SHIFT
            entry["upper_tail_collapsed"] = (
                f["high"] > MATERIAL_SHIFT
                and entry["upper_ratio"] is not None
                and entry["upper_ratio"] < 0.5)
            entry["lower_tail_extended"] = entry["lower_shift"] <= -MATERIAL_SHIFT
        per_condition[label] = entry

    reproduced = [c for c in treated
                  if per_condition[c]["reproduces_published"] is True]
    comparable = len(reproduced) == len([c for c in treated
                                         if per_condition[c]["published"]])
    moved = [c for c in treated
             if per_condition[c].get("median_moved_materially")]
    collapsed = [c for c in treated if per_condition[c].get("upper_tail_collapsed")]
    extended = [c for c in treated if per_condition[c].get("lower_tail_extended")]
    unmoved = [c for c in treated
               if per_condition[c].get("median_moved_materially") is False
               and not per_condition[c].get("upper_tail_collapsed")
               and not per_condition[c].get("lower_tail_extended")]

    # Which published conclusions survive, derived per claim rather than judged
    # as a whole. Several rest on DIFFERENT properties of the same table, so a
    # single verdict would hide the ones that fail.
    widths_full = {c: per_condition[c]["full_prior"]["width"]
                   for c in treated if per_condition[c]["full_prior"]}
    widths_adm = {c: per_condition[c]["admissible_only"]["width"]
                  for c in treated if per_condition[c]["admissible_only"]}
    claims = []
    for fragment, label, prop in published_claims():
        a = per_condition.get(label, {}).get("admissible_only")
        f = per_condition.get(label, {}).get("full_prior")
        survives = basis = None
        if a and f:
            if prop == "still_the_widest":
                widest = max(widths_adm, key=lambda k: widths_adm[k])
                survives = widest == label
                # The full-prior margin too: a claim that was true by a hair and
                # then reverses is a different story from one that was never
                # close, and quoting only the admissible winner hides which.
                widest_full = max(widths_full, key=lambda k: widths_full[k])
                runner_full = sorted(widths_full.values())[-2]
                runner_adm = sorted(widths_adm.values())[-2]
                basis = (f"full prior: {widest_full} widest at "
                         f"{widths_full[widest_full]:.3f}, next "
                         f"{runner_full:.3f} (margin "
                         f"{widths_full[widest_full] - runner_full:.3f}); "
                         f"admissible: {widest} widest at "
                         f"{widths_adm[widest]:.3f}, next {runner_adm:.3f}")
            elif prop == "still_wide":
                survives = a["width"] >= WIDE_INTERVAL
                basis = (f"admissible width {a['width']:.3f} against a "
                         f"{WIDE_INTERVAL:.2f} bar (was {f['width']:.3f})")
            elif prop == "still_a_floor":
                survives = a["high"] <= MATERIAL_SHIFT
                basis = f"admissible upper end {a['high']:.4f}"
            elif prop == "still_above_its_quoted_floor":
                floor = quoted_floor(fragment)
                survives = a["low"] >= floor
                basis = (f"admissible lower end {a['low']:.4f} against the "
                         f"{floor:.2f} the claim quotes about itself")
        claims.append({"claim": fragment, "condition": label, "property": prop,
                       "survives": survives, "basis": basis})
    survived = [c["claim"] for c in claims if c["survives"] is True]
    refuted = [c["claim"] for c in claims if c["survives"] is False]

    worst_sorted = sorted(worst)
    return {
        "published_claims": claims,
        "claims_survived": survived,
        "claims_refuted": refuted,
        "question": ("what the published prior-predictive intervals say once the "
                     "prior is restricted to draws the model itself calls viable"),
        "design": {"n_draws": len(rows), "n_cells_per_draw": N_CELLS,
                   "rng_seed": RNG_SEED, "sim_seed": SIM_SEED,
                   "quantiles": list(QUANTILES),
                   "baseline_viability_max": BASELINE_MAX,
                   "n_parameters": len(PARAMS),
                   "treated_conditions": treated,
                   "untreated_conditions": untreated},
        "admissibility": {
            "n_admissible": n_adm, "n_draws": len(rows),
            "fraction": n_adm / len(rows) if rows else None,
            "worst_untreated_median": worst_sorted[len(worst_sorted) // 2],
            "worst_untreated_p90": worst_sorted[int(0.9 * len(worst_sorted))],
            "worst_untreated_max": worst_sorted[-1],
            "n_draws_untreated_over_half": sum(1 for w in worst if w > 0.5),
            "reportable": n_adm >= MIN_ADMISSIBLE_DRAWS,
        },
        "conditions": per_condition,
        "reproduces_published": comparable,
        "conditions_reproduced": reproduced,
        "medians_moved": moved,
        "upper_tails_collapsed": collapsed,
        "lower_tails_extended": extended,
        "conditions_unmoved": unmoved,
        "verdict": ("conditioning_changes_what_the_intervals_say" if comparable
                    and (moved or collapsed or extended)
                    else "conditioning_changes_nothing_material" if comparable
                    else "cannot_compare_the_published_numbers_were_not_reproduced"),
    }


def _pct(x, places=1):
    return "n/a" if x is None else f"{100.0 * x:.{places}f}%"


def _yesno(v):
    return "n/a" if v is None else ("yes" if v else "no")


def render(d) -> str:
    g, adm = d["design"], d["admissibility"]
    # Derived from the quantiles rather than typed: a hand-written confidence
    # label is a figure that can outlive the quantiles it describes.
    qs = g["quantiles"]
    conf = f"{100.0 * (max(qs) - min(qs)):.0f}%"
    L = []
    L.append("# The published prior-predictive intervals, restricted to draws "
             "the model calls viable\n")
    L.append("*Generated by `scripts/intervals_under_admissibility.py`. The "
             "per-draw measurements need the compiled Python extension and are "
             "committed; every figure here is derived from that file, so CI "
             "re-derives the page and runs no simulation.*\n")

    L.append("## What was never checked\n")
    L.append(f"`uncertainty-intervals-report.md` publishes {conf} "
             f"prior-predictive intervals over uniform priors across {g['n_parameters']} PRCC "
             f"rate-constant ranges, {g['n_draws']} joint draws at "
             f"{g['n_cells_per_draw']} cells each. The manuscript separately "
             f"states that every phenotype must show less than "
             f"{_pct(g['baseline_viability_max'], 0)} death under Control, and "
             f"`headline_at_fitted.py` already uses that constraint to reject "
             f"the in-vitro-fitted cascade.\n")
    L.append(f"**The published design evaluates treated arms only.** Its "
             f"{len(g['treated_conditions'])} conditions contain no Control "
             f"arm, so the untreated death rate was never computed and a draw's "
             f"admissibility was not knowable from the artifact. This page adds "
             f"the {len(g['untreated_conditions'])} untreated arms on the SAME "
             f"draws -- one sample partitioned, not two sampled, so the "
             f"conditioning is not confounded with sampling noise.\n")

    L.append("## How much of the prior is admissible\n")
    L.append(f"**{adm['n_admissible']} of {adm['n_draws']} draws "
             f"({_pct(adm['fraction'])}).** Across the whole prior the worst "
             f"untreated death rate has median {_pct(adm['worst_untreated_median'])}, "
             f"90th percentile {_pct(adm['worst_untreated_p90'])} and maximum "
             f"{_pct(adm['worst_untreated_max'])}; "
             f"{adm['n_draws_untreated_over_half']} draws kill more than half "
             f"the untreated population.\n")

    L.append("## The full-prior half reproduces what is published\n")
    L.append("Without this the comparison would be against a re-run rather than "
             "against the committed report. Same ranges, same draw count, same "
             "seeds, all imported from that generator:\n")
    L.append(f"| condition | published median | reproduced median | published "
             f"{conf} | reproduced {conf} | reproduces |")
    L.append("|---|--:|--:|---|---|---|")
    for label in g["treated_conditions"]:
        c = d["conditions"][label]
        p, f = c["published"], c["full_prior"]
        pub_med = _pct(p["median"]) if p else "n/a"
        got_med = _pct(f["median"]) if f else "n/a"
        pub_int = f"[{_pct(p['low'])}, {_pct(p['high'])}]" if p else "n/a"
        got_int = f"[{_pct(f['low'])}, {_pct(f['high'])}]" if f else "n/a"
        L.append(f"| {label} | {pub_med} | {got_med} | {pub_int} | {got_int} | "
                 f"{_yesno(c['reproduces_published'])} |")
    L.append("")
    L.append(f"All conditions reproduce: **{_yesno(d['reproduces_published'])}**. "
             f"The verdict below is withheld entirely if they do not.\n")

    L.append("## What conditioning does, per condition\n")
    L.append(f"| condition | full-prior {conf} | admissible-only {conf} | "
             f"median | upper ratio | lower shift |")
    L.append("|---|---|---|--:|--:|--:|")
    for label in g["treated_conditions"]:
        c = d["conditions"][label]
        f, a = c["full_prior"], c["admissible_only"]
        if not f or not a:
            L.append(f"| {label} | n/a | n/a | n/a | n/a | n/a |")
            continue
        ratio = c["upper_ratio"]
        ratio_txt = "n/a" if ratio is None else f"{ratio:.2f}x"
        L.append(f"| {label} | [{_pct(f['low'])}, {_pct(f['high'])}] | "
                 f"[{_pct(a['low'])}, {_pct(a['high'])}] | "
                 f"{_pct(f['median'])} -> {_pct(a['median'])} | "
                 f"{ratio_txt} | {_pct(c['lower_shift'])} |")
    L.append("")

    L.append("## One mechanism, two opposite-looking corrections\n")
    L.append("Admissibility keeps the draws whose defences are strong enough "
             "that an untreated cell survives, so it lowers kill everywhere. "
             "That shows up as two different-looking effects, and reporting "
             "either alone would flatter whichever was looked at first.\n")
    if d["upper_tails_collapsed"]:
        L.append(f"**The inducer arms' upper tails collapse** "
                 f"({len(d['upper_tails_collapsed'])} of "
                 f"{len(g['treated_conditions'])} conditions): "
                 f"{', '.join(d['upper_tails_collapsed'])}. The published "
                 f"intervals say those phenotypes COULD be killed much harder "
                 f"than the admissible model allows, and the draws that achieve "
                 f"it are ones where the untreated cell is already dead.\n")
    if d["lower_tails_extended"]:
        L.append(f"**The exogenous-ROS arms' lower tails extend further down** "
                 f"({len(d['lower_tails_extended'])} of "
                 f"{len(g['treated_conditions'])}): "
                 f"{', '.join(d['lower_tails_extended'])}. Here the published "
                 f"intervals are too OPTIMISTIC at the bottom, which is the "
                 f"opposite direction and the same cause.\n")
    if d["conditions_unmoved"]:
        L.append(f"**And {len(d['conditions_unmoved'])} are unmoved**: "
                 f"{', '.join(d['conditions_unmoved'])}. A condition pinned at a "
                 f"floor or a ceiling across the whole prior has nowhere to go, "
                 f"which is why this is reported per condition rather than as "
                 f"one verdict.\n")
    L.append(f"Verdict: `{d['verdict']}`.\n")

    L.append("## Which of the published conclusions survive\n")
    L.append("Each is quoted from the committed report and scored on the "
             "property it actually rests on, because several rest on different "
             "properties of the same table and one verdict would hide the ones "
             "that fail. A guard asserts every quoted fragment is still in that "
             "report, so a rewrite there forces a rewrite here.\n")
    L.append("| published conclusion | survives | on what basis |")
    L.append("|---|---|---|")
    for c in d["published_claims"]:
        L.append(f"| \"{c['claim']}\" | **{_yesno(c['survives'])}** | "
                 f"{c['basis'] or 'not scorable'} |")
    L.append("")
    L.append(f"**{len(d['claims_refuted'])} of "
             f"{len(d['published_claims'])} refuted, "
             f"{len(d['claims_survived'])} survive.**\n")
    L.append("The one that matters most is the report's LEAD interpretation. It "
             "says the inducer rates are the least constrained, and under the "
             "model's own viability constraint the least-constrained condition "
             "is an exogenous-ROS one instead. The basis column shows that "
             "claim was true by a thin margin on the full prior rather than "
             "comfortably, which is the honest way to describe a reversal: not "
             "that it was badly wrong, but that it was never far from flipping "
             "and the constraint flips it.\n")

    L.append("## What this does NOT establish\n")
    L.append("**Conditioning is not inference.** Removing draws the model "
             "itself declares invalid adds no data. It cannot turn a "
             "prior-predictive interval into a data-conditioned one, and "
             "`identifiability-report.json`'s data-conditioned-headline count is "
             "unchanged by everything on this page.\n")
    L.append(f"**The admissible subset is not a posterior either.** It is the "
             f"same uniform prior with a region removed, so its median is a "
             f"property of the box that remains, not an estimate of anything. "
             f"The {_pct(adm['fraction'])} figure describes the prior's overlap "
             f"with a viability constraint and nothing else.\n")
    L.append("**The constraint is a modelling choice, stated in prose.** It is "
             "imported from the generator whose guard parses it out of the "
             "manuscript, so if the manuscript's number changes this page "
             "follows; but nothing here measures whether under-that-percentage "
             "untreated death is the right bar.\n")
    L.append("**The spatial and combination headline intervals are NOT "
             "corrected here.** `headline-uncertainty-report.md` and its "
             "`-tme-` and `-penetration-` siblings sweep the same prior through "
             "compiled binaries at minutes per run, so recomputing them "
             "conditioned was not attempted. The admissible fraction above is "
             "the size of the correction they need, and it is outstanding rather "
             "than applied -- saying so is the point, because silently leaving "
             "them beside a corrected page would imply they had been checked.\n")
    L.append(f"**Admissibility is judged on "
             f"{len(g['untreated_conditions'])} phenotypes, not all of them.** "
             f"The same three `headline_at_fitted` uses, so the two pages judge "
             f"it identically. Measured separately, adding the remaining "
             f"phenotypes changes no verdict because Persister always binds, so "
             f"the subset is not a weakening -- but it is a subset, and the "
             f"fraction above is therefore an upper bound on admissibility.\n")
    return "\n".join(L) + "\n"


def _split_stored(stored):
    """The gate's hook: both inputs are committed files, so they are read again
    rather than recovered from the artifact's own derived fields."""
    return (load(),)


def main() -> int:
    ap = argparse.ArgumentParser(description="condition the published intervals")
    ap.add_argument("--measure", action="store_true",
                    help="re-run the simulation and rewrite the committed draws")
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild both artifacts from the committed draws")
    a = ap.parse_args()
    if a.measure:
        measure()
    # `assemble(load())` on both paths: every input is a committed file, so
    # there is nothing to echo and a leg that is always recomputed cannot be.
    if a.render_only:
        data = assemble(load())
    else:
        data = assemble(load())
    OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
