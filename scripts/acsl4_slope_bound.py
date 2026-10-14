#!/usr/bin/env python3
"""Bound the ACSL4 bridge's slope, and say why the data cannot fix it (#934).

`acsl4::status_from_zscore` maps a within-cohort ACSL4 mRNA z-score onto the
`acsl4_status_boost` the peroxidation term reads, and the module names its own
placeholder: the slope that divides z is "the placeholder that reproduces the
existing constants". It was not merely uncalibrated -- it was UNREACHABLE. The
field is absent from `apply_param_overrides`, so neither the env hook nor the
Python binding could set it, and the only writer was one binary's hard-coded
preset. Driving it is the enabling change; this page is what the drive bought.

THREE THINGS THIS MEASURES, ALL OF WHICH CONSTRAIN HOW THE LAYER CAN BE USED
---------------------------------------------------------------------------
* WHERE THE OBSERVABLE EXISTS AT ALL. The boost raises the peroxidizable
  substrate for treated and untreated cells alike, so pushing it until the
  treated arm responds can push the UNTREATED arm through the model's own
  baseline-viability constraint. Those two requirements are measured against
  each other per phenotype, and they do not both hold everywhere.
* WHAT THE CELL LINES SAY. The committed CTRPv2 join carries ACSL4 expression
  and a sensitivity readout per line, so the measured sensitivity-versus-z
  relation can be binned and compared with the engine's.
* THE BOUND THAT FOLLOWS. The slope enters as `z / slope`, so at the top of the
  bridge's own range the boost is one over the slope. Requiring the layer to do
  something at that point, and the untreated arm to stay admissible there,
  bounds the slope from both sides.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
It does not fit the slope to the cell lines. The engine emits a single-cell
death rate at one implicit exposure; the readout is a population viability
averaged over a dose range. A threshold response averaged over a heterogeneous
population and a range of doses looks shallow, so a shallow measured curve
cannot distinguish a wrong functional form from the right one seen through that
averaging. The data supplies the DIRECTION and the engine's own admissibility
supplies the bound; neither supplies a fitted value, and the page says so.

Nor does it move the identifiability report's data-conditioned-headline count.
The slope feeds an off-by-default stratification layer that no headline output
reads.

OFFLINE CONTRACT
----------------
The engine measurements need the compiled Python extension, so they are taken
locally with `--measure` and committed. The cell-line leg reads a committed CSV
and is therefore recomputed on EVERY render, so it can never be echoed.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

# ONE Python source for the model's own baseline-viability constraint, which
# `headline_at_fitted` already owns and whose guard parses it from the
# manuscript. Re-declaring the number here is how two copies disagree.
from headline_at_fitted import BASELINE_MAX  # noqa: E402

OUT_JSON = REPO / "analysis" / "calibration" / "acsl4-slope-bound.json"
OUT_MD = REPO / "analysis" / "calibration" / "acsl4-slope-bound.md"

PANEL = REPO / "analysis" / "calibration" / "depmap_ctrpv2_panel.csv"

# The parameter the bridge writes, now reachable from the override hook.
BOOST_PARAM = "acsl4_status_boost"
BRIDGE_GENE = "ACSL4"

# The bridge's own range: the integer-z points `acsl4.rs` documents as
# reproducing its shipped status constants.
BRIDGE_Z_LOW = -2.0
BRIDGE_Z_HIGH = 1.0

# Shipped slope: `status_from_zscore(z) = max(0, 1 + z/SHIPPED_SLOPE)`.
SHIPPED_SLOPE = 2.0

PHENOTYPES = ("Glycolytic", "OXPHOS", "Persister", "PersisterNrf2")
TREATMENT = "RSL3"
UNTREATED = "Control"

N_CELLS = 4000
SEED = 42

# A boost grid dense enough to locate both ends of the bound, and extending
# well past the bridge so "unresponsive in range" is distinguishable from
# "unresponsive anywhere".
BOOST_GRID = (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60,
              0.75, 0.90, 1.00, 1.25, 1.50, 1.75, 2.00, 3.00, 4.00)

# A treated kill below this is not an effect the layer could be said to have.
EFFECT_FLOOR = 0.01

# Compounds to read the measured relation from, and the z bin width.
COMPOUNDS = ("ML162", "ML210")
Z_BIN = 0.5
Z_CLAMP = 2.0
MIN_LINES_PER_BIN = 10


def _f(text):
    """A blank field is an absent value, never a zero."""
    if text is None or text == "":
        return None
    return float(text)


def status_from_zscore(z: float, slope: float = SHIPPED_SLOPE) -> float:
    """The bridge, at an arbitrary slope. Reproduces `acsl4.rs` at the default."""
    return max(0.0, 1.0 + z / slope)


def boost_from_zscore(z: float, slope: float = SHIPPED_SLOPE) -> float:
    """The bridge composed with `pufa_boost_from_status`: `status - 1`, floored."""
    return max(-1.0, status_from_zscore(z, slope) - 1.0)


def measure_engine() -> dict:
    """Per phenotype, the untreated and treated kill across the boost grid.

    Needs the compiled extension. Both arms are measured at every point,
    because the question is not how much the treated arm moves but whether the
    untreated arm stays admissible while it does.
    """
    import ferroptosis_core as fc

    out = {}
    for phenotype in PHENOTYPES:
        rows = []
        for boost in BOOST_GRID:
            kwargs = {BOOST_PARAM: boost}
            untreated = fc.sim_batch(phenotype, UNTREATED, N_CELLS, SEED, **kwargs)
            treated = fc.sim_batch(phenotype, TREATMENT, N_CELLS, SEED, **kwargs)
            rows.append({"boost": boost,
                         "untreated_kill": untreated["death_rate"],
                         "treated_kill": treated["death_rate"],
                         # Carried so the effect floor can be shown to clear the
                         # sampling resolution rather than assumed to.
                         "treated_ci_low": treated.get("ci_low"),
                         "treated_ci_high": treated.get("ci_high")})
        out[phenotype] = rows
        print(f"  {phenotype}: {len(rows)} boost points", flush=True)
    return out


def measure_cell_lines() -> dict:
    """The measured sensitivity-versus-z relation, from the committed panel.

    Recomputed on every render rather than stored: it reads a file in the
    repository, so there is nothing to gain by trusting a copy of it, and a leg
    that is always recomputed cannot be echoed.
    """
    rows = list(csv.DictReader(PANEL.open()))
    out = {}
    for compound in COMPOUNDS:
        pairs = []
        for row in rows:
            expression = _f(row.get(f"expr_{BRIDGE_GENE}"))
            viability = _f(row.get(f"viability_auc_{compound}"))
            if expression is None or viability is None:
                continue
            pairs.append((expression, viability))
        if len(pairs) < MIN_LINES_PER_BIN:
            out[compound] = {"n_lines": len(pairs), "bins": [], "usable": False}
            continue
        values = [e for e, _ in pairs]
        mean = statistics.mean(values)
        sd = statistics.pstdev(values)
        binned = {}
        for expression, viability in pairs:
            z = (expression - mean) / sd if sd else 0.0
            key = max(-Z_CLAMP, min(Z_CLAMP, round(z / Z_BIN) * Z_BIN))
            binned.setdefault(key, []).append(viability)
        bins = [{"z": key, "n_lines": len(v),
                 "mean_viability": statistics.mean(v),
                 "mean_sensitivity": 1.0 - statistics.mean(v)}
                for key, v in sorted(binned.items())
                if len(v) >= MIN_LINES_PER_BIN]
        zs = sorted((e - mean) / sd if sd else 0.0 for e, _ in pairs)
        out[compound] = {"n_lines": len(pairs), "expression_mean": mean,
                         "expression_sd": sd, "bins": bins,
                         "z_observed_low": zs[0], "z_observed_high": zs[-1],
                         "usable": len(bins) >= 2}
    return out


def _responsive_and_admissible(rows) -> dict:
    """The two requirements, measured against each other rather than assumed.

    `boost_responsive` is the smallest boost at which the treated arm clears the
    effect floor; `boost_admissible_max` the largest at which the UNTREATED arm
    still satisfies the model's own baseline-viability constraint. A window
    exists only where the first is at or below the second.
    """
    responsive = [r["boost"] for r in rows if r["treated_kill"] >= EFFECT_FLOOR]
    admissible = [r["boost"] for r in rows if r["untreated_kill"] < BASELINE_MAX]
    low = min(responsive) if responsive else None
    high = max(admissible) if admissible else None
    # The grid must have CONTAINED the crossing, or the endpoint is where the
    # grid stopped rather than where the model changed. This is the defect the
    # first run of this analysis shipped: a bound read off the last grid point.
    inadmissible = [r["boost"] for r in rows if r["untreated_kill"] >= BASELINE_MAX]
    bracketed = bool(inadmissible) and high is not None and high < min(inadmissible)
    window = (low is not None and high is not None and low <= high)
    bridge_top_boost = boost_from_zscore(BRIDGE_Z_HIGH)
    return {
        "boost_responsive_from": low,
        "boost_admissible_to": high,
        "admissibility_crossing_bracketed": bracketed,
        "window_exists": window,
        "responsive_inside_bridge": (low is not None and low <= bridge_top_boost),
        "admissible_inside_bridge": (high is not None and high >= bridge_top_boost),
        "usable_for_calibration": bool(
            window and low is not None and low <= bridge_top_boost
            and high is not None and high >= bridge_top_boost),
    }


def _slope_bound(window) -> dict:
    """Turn a boost window into a slope bound.

    The bridge is `1 + z/slope`, so at the top of its own range the boost is
    exactly one over the slope. A boost window therefore inverts to a slope
    window, and the ORDER flips: the smallest responsive boost gives the LARGEST
    admissible slope.
    """
    low, high = window["boost_responsive_from"], window["boost_admissible_to"]
    if not window["usable_for_calibration"] or not low or not high:
        return {"slope_low": None, "slope_high": None,
                "shipped_inside": None, "basis": "no usable window"}
    slope_low, slope_high = BRIDGE_Z_HIGH / high, BRIDGE_Z_HIGH / low
    return {"slope_low": slope_low, "slope_high": slope_high,
            "shipped_slope": SHIPPED_SLOPE,
            "shipped_inside": slope_low <= SHIPPED_SLOPE <= slope_high,
            "width": slope_high - slope_low,
            "basis": "effect floor above, baseline-viability constraint below"}


def _sampling_resolution(engine) -> dict:
    """The widest Wilson half-width at the boost where the effect floor bites.

    A floor below the sampling resolution would call noise an effect. Reported
    rather than assumed, because `n` is a declared constant and someone will
    lower it.
    """
    widths = []
    for rows in engine.values():
        for r in rows:
            low, high = r.get("treated_ci_low"), r.get("treated_ci_high")
            if low is None or high is None:
                continue
            if r["treated_kill"] <= EFFECT_FLOOR * 3:
                widths.append((high - low) / 2.0)
    if not widths:
        return {"usable": False}
    worst = max(widths)
    return {"usable": True, "n_points": len(widths),
            "worst_half_width": worst,
            "effect_floor": EFFECT_FLOOR,
            "floor_clears_resolution": EFFECT_FLOOR > worst}


def _observed_bound(window, lines) -> dict:
    """The bound the DATA tightens, as opposed to the one the bridge's own
    documented range gives.

    The bridge documents z from its floor to its top, but real lines reach
    further out than that, and the admissibility ceiling applies to whatever z a
    line actually has. So requiring EVERY OBSERVED line to stay admissible is a
    stronger requirement than requiring it at the documented top, and it is the
    one the committed join can impose. Reported beside the documented-range
    bound, never instead of it: they answer different questions.
    """
    ceiling = window["boost_admissible_to"]
    highs = [leg.get("z_observed_high") for leg in lines.values()
             if leg.get("z_observed_high") is not None]
    if not window["usable_for_calibration"] or not ceiling or not highs:
        return {"usable": False}
    z_max = max(highs)
    slope_low = z_max / ceiling
    inadmissible = {}
    for compound, leg in lines.items():
        high = leg.get("z_observed_high")
        if high is None:
            continue
        inadmissible[compound] = {
            "z_observed_high": high,
            "boost_at_shipped_slope": high / SHIPPED_SLOPE,
            "admissible_at_shipped_slope": high / SHIPPED_SLOPE <= ceiling,
        }
    return {"usable": True, "z_observed_high": z_max,
            "slope_low": slope_low,
            "shipped_above_it": SHIPPED_SLOPE >= slope_low,
            "tightening_factor": (slope_low / (BRIDGE_Z_HIGH / ceiling)
                                  if ceiling else None),
            "per_compound": inadmissible}


def _measured_shape(leg) -> dict:
    """Does the measured relation rise, and is it a threshold or a slope?

    Reported as a RANGE and a monotonicity count rather than a fitted slope: a
    binned population average is not a dose-response and giving it a slope
    would invite it to be read as one.
    """
    bins = leg.get("bins") or []
    if len(bins) < 2:
        return {"usable": False}
    first, last = bins[0], bins[-1]
    rises = sum(1 for a, b in zip(bins, bins[1:])
                if b["mean_sensitivity"] > a["mean_sensitivity"])
    span = last["mean_sensitivity"] - first["mean_sensitivity"]
    ratio = (last["mean_sensitivity"] / first["mean_sensitivity"]
             if first["mean_sensitivity"] else None)
    return {
        "usable": True,
        "z_low": first["z"], "z_high": last["z"],
        "sensitivity_low": first["mean_sensitivity"],
        "sensitivity_high": last["mean_sensitivity"],
        "absolute_span": span, "fold_change": ratio,
        "steps_rising": rises, "steps_total": len(bins) - 1,
        "monotone": rises == len(bins) - 1,
        "direction_matches_engine": span > 0,
    }


def _engine_shape(rows) -> dict:
    """The engine's treated response across the bridge's own range."""
    low_boost = boost_from_zscore(BRIDGE_Z_LOW)
    high_boost = boost_from_zscore(BRIDGE_Z_HIGH)
    inside = [r for r in rows if low_boost <= r["boost"] <= high_boost]
    if len(inside) < 2:
        return {"usable": False}
    kills = [r["treated_kill"] for r in inside]
    return {"usable": True, "n_points": len(inside),
            "kill_at_bridge_floor": inside[0]["treated_kill"],
            "kill_at_bridge_top": inside[-1]["treated_kill"],
            "absolute_span": kills[-1] - kills[0],
            "flat_inside_bridge": max(kills) - min(kills) < EFFECT_FLOOR}


def assemble(raw) -> dict:
    engine, lines = raw["engine"], raw["cell_lines"]
    windows = {p: _responsive_and_admissible(engine[p]) for p in sorted(engine)}
    usable = [p for p in sorted(windows) if windows[p]["usable_for_calibration"]]
    bounds = {p: _slope_bound(windows[p]) for p in sorted(windows)}
    measured = {c: _measured_shape(lines[c]) for c in sorted(lines)}
    shapes = {p: _engine_shape(engine[p]) for p in sorted(engine)}

    observed = {p: _observed_bound(windows[p], lines) for p in sorted(windows)}
    agreeing = [c for c in sorted(measured)
                if measured[c].get("usable") and measured[c]["direction_matches_engine"]]
    monotone = [c for c in sorted(measured)
                if measured[c].get("usable") and measured[c]["monotone"]]
    bounded = [p for p in usable if bounds[p]["shipped_inside"]
               and observed[p].get("shipped_above_it")]
    return {
        "question": ("what can be said about the ACSL4 bridge's slope now that "
                     "the parameter it writes can be driven at all"),
        "bridge": {"gene": BRIDGE_GENE, "parameter": BOOST_PARAM,
                   "shipped_slope": SHIPPED_SLOPE,
                   "z_low": BRIDGE_Z_LOW, "z_high": BRIDGE_Z_HIGH,
                   "boost_at_z_low": boost_from_zscore(BRIDGE_Z_LOW),
                   "boost_at_z_high": boost_from_zscore(BRIDGE_Z_HIGH)},
        "engine": engine,
        "cell_lines": {c: {k: v for k, v in lines[c].items() if k != "bins"}
                       for c in sorted(lines)},
        "measured_bins": {c: lines[c].get("bins", []) for c in sorted(lines)},
        "windows": windows,
        "slope_bounds": bounds,
        "measured_shape": measured,
        "observed_range_bound": observed,
        "engine_shape": shapes,
        "phenotypes_usable": usable,
        "compounds_agreeing_in_direction": agreeing,
        "compounds_monotone": monotone,
        "sampling_resolution": _sampling_resolution(engine),
        "constraints": {"baseline_viability_max": BASELINE_MAX,
                        "effect_floor": EFFECT_FLOOR,
                        "n_cells": N_CELLS, "seed": SEED},
        "verdict": ("slope_bounded_not_fitted" if bounded
                    else "slope_not_bounded_by_this_route"),
    }


def _num(value, places=4):
    return "n/a" if value is None else f"{value:.{places}f}"


def _pct(fraction, places=2):
    return "n/a" if fraction is None else f"{100.0 * fraction:.{places}f}%"


def _yesno(value):
    return "n/a" if value is None else ("yes" if value else "no")


def render(d) -> str:
    b = d["bridge"]
    cons = d["constraints"]
    L = []
    L.append("# What can be said about the ACSL4 bridge's slope, now that it can "
             "be driven\n")
    L.append("*Generated by `scripts/acsl4_slope_bound.py`. The engine "
             "measurements need the compiled Python extension and are committed; "
             "the cell-line leg reads a committed CSV and is recomputed on every "
             "render, so it can never be echoed.*\n")

    L.append("## The parameter was not uncalibrated, it was unreachable\n")
    L.append(f"`acsl4::status_from_zscore` maps a within-cohort {b['gene']} mRNA "
             f"z-score onto `{b['parameter']}`, which the peroxidation term "
             f"reads, and the module names the divisor as its own placeholder. "
             f"That field was absent from `apply_param_overrides`, so neither "
             f"the environment hook nor the Python binding could set it and the "
             f"only writer was one binary's hard-coded preset. A slope "
             f"calibration was therefore not declined or deferred -- there was "
             f"no path to the parameter. Adding the arm is byte-identical by "
             f"construction, because that function only consults keys it is "
             f"handed.\n")
    L.append(f"**The bound is anchored at one z, and that choice is part of "
             f"it.** At z {b['z_high']:+.0f} -- the top of the bridge's own "
             f"documented range, its ACSL4-high point -- the boost is one over "
             f"the slope, which is what makes the window invertible. At z "
             f"{0:+.0f} the boost is zero for every slope and no bound exists "
             f"at all, so a bound quoted without its anchor would be "
             f"meaningless. The requirement being imposed is that the layer do "
             f"something at its most sensitive point and stay admissible "
             f"there.\n")
    L.append(f"The bridge's own documented range is z {b['z_low']:+.0f} to "
             f"{b['z_high']:+.0f}, which at the shipped slope of "
             f"{b['shipped_slope']:.0f} is a boost of "
             f"{_num(b['boost_at_z_low'], 2)} to {_num(b['boost_at_z_high'], 2)}. "
             f"The slope enters as `z / slope`, so the boost at the top of that "
             f"range is exactly one over the slope -- which is what makes a "
             f"boost window invertible into a slope window.\n")

    L.append("## Where the observable exists at all\n")
    L.append(f"The boost raises the peroxidizable substrate for treated and "
             f"untreated cells alike. So pushing it until {TREATMENT} responds "
             f"can push the UNTREATED arm through the model's own "
             f"baseline-viability constraint of {_pct(cons['baseline_viability_max'], 0)}, "
             f"and those two requirements have to be measured against each "
             f"other rather than assumed compatible:\n")
    L.append("| phenotype | responds from | admissible to | crossing bracketed "
             "| window | usable in bridge range |")
    L.append("|---|--:|--:|---|---|---|")
    for p in sorted(d["windows"]):
        w = d["windows"][p]
        L.append(f"| {p} | {_num(w['boost_responsive_from'], 2)} | "
                 f"{_num(w['boost_admissible_to'], 2)} | "
                 f"{_yesno(w['admissibility_crossing_bracketed'])} | "
                 f"{_yesno(w['window_exists'])} | "
                 f"{_yesno(w['usable_for_calibration'])} |")
    L.append("")
    usable = d["phenotypes_usable"]
    L.append(f"**{len(usable)} of {len(d['windows'])} phenotypes can carry the "
             f"measurement at all**"
             + (f": {', '.join(usable)}." if usable else ".")
             + " That is not a detail of presentation. The phenotype decides "
               "whether the layer has an effect to measure inside the range its "
               "own bridge spans, and in the others the treated arm is flat "
               "there, or the untreated arm is already dying.\n")
    L.append("The `crossing bracketed` column is there because the first run of "
             "this analysis read a bound off the LAST GRID POINT: the untreated "
             "arm had not yet crossed the constraint when the grid stopped, so "
             "the number reported as a measured limit was where I stopped "
             "looking. The grid now has to contain the crossing, and the column "
             "says whether it did.\n")

    L.append("## The bound that follows\n")
    L.append("| phenotype | slope at least | slope at most | shipped inside | basis |")
    L.append("|---|--:|--:|---|---|")
    for p in sorted(d["slope_bounds"]):
        s = d["slope_bounds"][p]
        L.append(f"| {p} | {_num(s['slope_low'], 2)} | {_num(s['slope_high'], 2)} "
                 f"| {_yesno(s.get('shipped_inside'))} | {s['basis']} |")
    L.append("")
    L.append("**And the data tightens the lower end, which is more than "
             "supplying a direction.** The bridge documents its range as z "
             f"{b['z_low']:+.0f} to {b['z_high']:+.0f}, but real lines reach "
             "further out, and the admissibility ceiling applies to whatever z "
             "a line actually has. Requiring EVERY OBSERVED line to stay "
             "admissible is therefore a stronger requirement than requiring it "
             "at the documented top:\n")
    L.append("| phenotype | observed top z | slope at least | shipped above it "
             "| tighter than the documented bound by |")
    L.append("|---|--:|--:|---|--:|")
    for p in sorted(d["observed_range_bound"]):
        o = d["observed_range_bound"][p]
        if not o.get("usable"):
            L.append(f"| {p} | -- | -- | -- | -- |")
            continue
        L.append(f"| {p} | {o['z_observed_high']:+.2f} | "
                 f"{_num(o['slope_low'], 2)} | "
                 f"{_yesno(o['shipped_above_it'])} | "
                 f"{_num(o['tightening_factor'], 2)}x |")
    L.append("")
    for p in sorted(d["observed_range_bound"]):
        o = d["observed_range_bound"][p]
        if not o.get("usable"):
            continue
        rows = o["per_compound"]
        bad = [c for c in sorted(rows) if not rows[c]["admissible_at_shipped_slope"]]
        L.append(f"At the shipped slope the most ACSL4-high line maps to a boost "
                 f"of "
                 + ", ".join(f"{_num(rows[c]['boost_at_shipped_slope'], 3)} ({c})"
                             for c in sorted(rows))
                 + f", against a ceiling of "
                   f"{_num(d['windows'][p]['boost_admissible_to'], 2)} -- so "
                 + ("**no observed line is inadmissible there**"
                    if not bad else
                    f"**{len(bad)} compound's top line is NOT admissible there** "
                    f"({', '.join(bad)})")
                 + ". That is a statement about the shipped value, not an "
                   "argument for it.\n")
        break
    L.append(f"Verdict: `{d['verdict']}`.\n")

    L.append("## What the cell lines say, and why it does not fit the slope\n")
    L.append(f"The committed join carries {b['gene']} expression and a "
             f"sensitivity readout per line. Binned by within-cohort z:\n")
    L.append("| compound | lines | z from | z to | sensitivity from | "
             "sensitivity to | fold | monotone |")
    L.append("|---|--:|--:|--:|--:|--:|--:|---|")
    for c in sorted(d["measured_shape"]):
        m = d["measured_shape"][c]
        if not m.get("usable"):
            L.append(f"| {c} | {d['cell_lines'][c]['n_lines']} | -- | -- | -- | "
                     f"-- | -- | too few bins |")
            continue
        L.append(f"| {c} | {d['cell_lines'][c]['n_lines']} | {m['z_low']:+.1f} | "
                 f"{m['z_high']:+.1f} | {_num(m['sensitivity_low'])} | "
                 f"{_num(m['sensitivity_high'])} | {_num(m['fold_change'], 2)}x | "
                 f"{m['steps_rising']} of {m['steps_total']} |")
    L.append("")
    L.append(f"**The direction agrees with the engine's for "
             f"{len(d['compounds_agreeing_in_direction'])} of "
             f"{len(d['measured_shape'])} compounds** -- higher {b['gene']}, more "
             f"sensitive -- which is the sign `acsl4.rs` fixes in source. The "
             f"SHAPE does not: the measured rise is shallow and spread across "
             f"the whole range, while the engine's treated response inside the "
             f"bridge is flat and then climbs.\n")
    L.append("| phenotype | kill at bridge floor | kill at bridge top | span | "
             "flat inside bridge |")
    L.append("|---|--:|--:|--:|---|")
    for p in sorted(d["engine_shape"]):
        e = d["engine_shape"][p]
        if not e.get("usable"):
            L.append(f"| {p} | -- | -- | -- | too few points |")
            continue
        L.append(f"| {p} | {_num(e['kill_at_bridge_floor'])} | "
                 f"{_num(e['kill_at_bridge_top'])} | {_num(e['absolute_span'])} | "
                 f"{_yesno(e['flat_inside_bridge'])} |")
    L.append("")
    L.append("**And that mismatch is NOT evidence the form is wrong**, which is "
             "why no slope is fitted here. The engine emits a single-cell death "
             "rate at one implicit exposure; the readout is a population "
             "viability averaged over a dose range. A threshold response "
             "averaged over a heterogeneous population and a range of doses "
             "looks shallow, so a shallow measured curve cannot distinguish a "
             "wrong functional form from the right one seen through that "
             "averaging. Separating them needs a per-line dose correspondence "
             "this repository does not have.\n")

    L.append("## What this does NOT establish\n")
    L.append("**The slope is bounded, not fitted, and the bound is the "
             "engine's own.** Its upper end comes from requiring the layer to "
             "have an effect and its lower end from the untreated arm staying "
             "admissible. Neither end is a measurement of ACSL4 biology. The "
             "cell lines supply the DIRECTION and nothing else.\n")
    L.append("**It does not move the data-conditioned-headline count.** The "
             "slope feeds an off-by-default stratification layer, and no "
             "headline output reads it, so `identifiability-report.json` is "
             "unchanged -- which has to be stated rather than left for a reader "
             "to discover.\n")
    res = d["sampling_resolution"]
    L.append(f"**The effect floor is a declared threshold, not a measured "
             f"one.** A treated kill below {_pct(cons['effect_floor'], 0)} is "
             f"called no effect here. That choice sets the upper end of every "
             f"slope bound in the table, and a different floor would move it. "
             f"What it does clear is the SAMPLING resolution: across the points "
             f"near the floor the widest interval half-width is "
             f"{_num(res.get('worst_half_width'))} at n = {cons['n_cells']}, so "
             f"the floor is "
             f"{'above' if res.get('floor_clears_resolution') else 'NOT above'} "
             f"the noise it would otherwise be reading.\n")
    L.append("**One phenotype is not a consensus.** The bound rests on the "
             "phenotypes where the observable exists, and where it does not the "
             "page reports no bound rather than a wider one.\n")
    L.append(f"**And the measured relation is weak.** The join found that "
             f"{b['gene']} does not clear its extremeness bar against the "
             f"all-gene distribution, so this is a shallow monotone trend, not a "
             f"strong predictor, and a slope fitted to it would have had a wide "
             f"admissible range even without the averaging problem.\n")
    return "\n".join(L) + "\n"


def _split_stored(stored):
    """The gate's hook. The engine leg comes back from the artifact; the
    cell-line leg is RE-MEASURED, because it reads a committed file."""
    return ({"engine": stored["engine"], "cell_lines": measure_cell_lines()},)


def main() -> int:
    ap = argparse.ArgumentParser(description="bound the ACSL4 bridge's slope")
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild from the committed engine measurements")
    a = ap.parse_args()
    if a.render_only:
        # RE-ASSEMBLED, not echoed: the stored file carries derived fields next
        # to the raw measurements, and re-rendering those instead of recomputing
        # them would make every guard on a derived field inert. The cell-line
        # leg is recomputed on this path too.
        stored = json.loads(OUT_JSON.read_text())
        data = assemble({"engine": stored["engine"],
                         "cell_lines": measure_cell_lines()})
    else:
        print(f"measuring {len(PHENOTYPES)} phenotypes x {len(BOOST_GRID)} "
              f"boosts x 2 arms at n={N_CELLS}", flush=True)
        data = assemble({"engine": measure_engine(),
                         "cell_lines": measure_cell_lines()})
    OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
