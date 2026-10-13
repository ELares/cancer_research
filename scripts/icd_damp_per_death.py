#!/usr/bin/env python3
"""The per-dead-cell DAMP advantage the engine documents, measured (#928).

`ferroptosis_core::immune::calculate_damp_release` carried a "Key biological
insight" above it, with two citations: SDT/PDT kill via a runaway LP cascade so
"LP at death is MUCH higher than threshold", RSL3 kills via slow GPX4 inhibition
so its "LP at death is NEAR threshold", "therefore SDT/PDT-killed cells release
MORE DAMPs per cell". Two committed documents QUANTIFIED it -- `sim-icd`'s
README told readers to expect ~1.5-2x under a heading reading "Reproducing
manuscript claims", and the v2 plan registered 2-5x as the expected result.

MEASURED, IN EVERY REGIME THIS REPOSITORY CAN RUN, THE RATIO IS UNDER 1.15x.
Three of them, none requiring new simulation code:

  1. `sim-icd` itself, whose own stderr prints the comparison. Persister
     1.1332x, OXPHOS 1.0748x.
  2. `sim-modality-panel --tme-sweep`, 32 readable conditions, 1.0674-1.1475x.
  3. The committed 2D spatial archive, per cell over every ferroptotic death in
     the baseline immune run: 1.0469x on release LP.

THE DIRECTION HOLDS AND ONLY THE MAGNITUDE FAILS, which is why this page
reports both and withdraws neither more than it has to: SDT sits above RSL3 in
every readable comparison.

AND THE MECHANISM IS CONTRADICTED AT THE MOMENT IT NAMES. The old text is about
LP *at death*, and the archive measures that per cell: 10.2657 under SDT against
10.2474 under RSL3, against a `death_threshold` of 10.0. Both are within a few
percent of threshold and differ by 0.18%, because death IS the threshold
crossing -- no arm can be far from it there. The gap opens only in the
POST-DEATH quantity the consumers actually pass (release LP 19.1230 vs
18.2658). So the old contrast read RSL3 at death against SDT after it.

AND IT LANDS ON A REGISTERED PREDICTION, which is where this bites hardest.
`PREREGISTRATION.md`'s P5 registers a falsification threshold on exactly this
quantity -- "DAMP release and DC maturation per dead cell agree within 1.5x
between dense (SDT or 3D) and sparse (RSL3) kill" -- and both registered
readouts sit well inside it (Persister 1.1332x DAMP / 1.0950x DC maturation;
OXPHOS 1.0748x / 1.0523x). So THE MODEL'S OWN PER-CELL OUTPUT SATISFIES P5'S
FALSIFICATION CONDITION. That does not falsify the biology, which is what
experiment E3 tests; it establishes that the model does not predict what P5's
model-output line says it predicts. That line states a per-cell ratio of "about
4:1 in 3D (down from the 2D 104:1 ceiling)", and the manuscript had already
recorded those as total immune-kill ratios without the registration being
corrected with it. On the per-cell axis the model predicts the null, so an E3
returning agreement within 1.5x would CONFIRM the model and falsify P5 as
registered. THE THRESHOLD IS NOT AMENDED and must not be: it is parsed from the
registration so this page tracks the registered band rather than carrying its
own, and the registration gets a dated measurement note that leaves the
sentence verbatim. Amending a registered threshold after measuring against it
is the thing preregistration exists to prevent.

That was already known HALF-WAY: `sim-modality-panel` records in its own source
that measuring at the crossing gave "~10.2" for both arms and that the
documented difference "vanished ... it was being measured before it happens".
The finding was never propagated back to the library comment it refers to,
which is the defect class this repository calls one-artifact-describes-another.

NOT CLAIMED. `damp_per_lp` is 1.0 and uncalibrated, so every ratio here is a
property of the model, not a measurement of immunogenicity. Nothing is
recalibrated. The manuscript is NOT the problem and is not edited: Section 8.2
already states that the linear LP-to-DAMP assumption "is not a measurement of
per-cell immunogenicity", Figure 16's caption already says the 104:1 "compares
total immune-kill counts, not DAMP potency or immunogenicity per dead cell",
and that section already withdrew a twofold HMGB1/ATP expectation of its own.
No manuscript figure or number rests on the 1.5-5x band.

OFFLINE CONTRACT: the spatial leg reads a committed archive and is recomputed on
every render. The `sim-icd` and panel legs need compiled binaries, so they are
measured locally with `--measure` and committed; CI reads only the JSON, and
`--render-only` re-assembles the prose from it idempotently.

Usage:
  icd_damp_per_death.py [--icd PATH] [--panel PATH]
  icd_damp_per_death.py --render-only
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from statistics import fmean

REPO = Path(__file__).resolve().parent.parent

OUT_JSON = REPO / "analysis" / "calibration" / "icd-damp-per-death.json"
OUT_MD = REPO / "analysis" / "calibration" / "icd-damp-per-death.md"

# The committed per-cell record of every ferroptotic death in the 2D baseline
# immune run. This is the ONLY place the repository stores LP at the moment of
# death alongside LP at release, which is what makes the two-time-point
# comparison possible at all.
SPATIAL_ARCHIVE = REPO / "analysis" / "immune-2d-measurements" / "observations.json.gz"

# `death_threshold` is PARSED rather than typed: "near threshold" is a claim
# about a number that lives in the engine, and a hardcoded 10.0 here would stop
# tracking it the moment anyone retuned the switch.
PARAMS_RS = REPO / "simulations" / "ferroptosis-core" / "src" / "params.rs"
DEATH_THRESHOLD_RE = re.compile(r"^\s*death_threshold:\s*([0-9.]+),", re.MULTILINE)

# The two documented predictions, each with the exact text it is quoted from.
# The numbers live here and the QUOTES are pinned against their source files by
# `tests/test_icd_damp_per_death.py`, so an edit to either document fails the
# guard rather than silently leaving this page comparing against a band nobody
# claims any more. Deleting the predictions outright was the alternative and is
# worse: a registered expectation that failed is evidence, and erasing it would
# erase the falsification along with the error.
PREDICTIONS = [
    {
        "name": "sim-icd README",
        "source": "simulations/sim-icd/README.md",
        "quote": "SDT produces ~1.5-2x more DAMP per dead cell than RSL3",
        "low": 1.5,
        "high": 2.0,
        "context": "stated as the binary's key finding, under "
                   '"Reproducing manuscript claims"',
    },
    {
        "name": "simulation-suite v2 plan",
        "source": "plans/simulation-suite-v2-plan.md",
        "quote": "SDT/PDT produce 2-5× more DAMP per dead cell than RSL3",
        "low": 2.0,
        "high": 5.0,
        "context": 'registered under "Expected Result" for Simulation 3',
    },
]

# P5 in `PREREGISTRATION.md` registers a falsification threshold on exactly
# the quantity this page measures, and it names TWO readouts. The number is
# PARSED from the registration rather than typed: a hardcoded 1.5 here would
# stop tracking the registered band, and the band is the whole point.
PREREG = REPO / "PREREGISTRATION.md"
P5_QUOTE = ("DAMP release and DC maturation per dead cell agree within 1.5x "
            "between dense (SDT or 3D) and sparse (RSL3) kill")
P5_BAND_RE = re.compile(
    r"DAMP release and DC maturation per dead cell agree within ([0-9.]+)x")

# A cell is "near" the threshold if it is within this fraction of it. Fixed
# before the measurement: the old text's own word is "NEAR", and the question is
# whether ONE arm is near it and the other far, not whether both are.
NEAR_THRESHOLD_TOLERANCE = 0.10

ARMS = ("SDT", "PDT", "RSL3", "AntibodyDrugConjugate")


def _default_binary(name: str) -> Path:
    """Locate a release binary across the layouts this repo actually uses."""
    env = os.getenv("CARGO_TARGET_DIR")
    roots = [Path(env)] if env else []
    roots += [REPO / "simulations" / "target", Path.home() / "ironauth-target"]
    for r in roots:
        cand = r / "release" / name
        if cand.exists():
            return cand
    return roots[-1] / "release" / name


def death_threshold() -> float:
    """The engine's own death threshold, read from `params.rs`."""
    m = DEATH_THRESHOLD_RE.search(PARAMS_RS.read_text(encoding="utf-8"))
    if m is None:
        raise RuntimeError(f"no death_threshold default found in {PARAMS_RS}")
    return float(m.group(1))


def measure_spatial() -> dict:
    """Per-arm LP at death and at release, from the committed 2D archive.

    Offline and recomputed on every render, including `--render-only`: it reads
    a committed input rather than a binary, so echoing a stored copy would make
    the guards on it inert for no benefit.
    """
    obs = json.loads(gzip.decompress(SPATIAL_ARCHIVE.read_bytes()))
    arms: dict[str, dict] = {}
    for cond in obs["conditions"]:
        name = cond["condition_name"]
        if not name.startswith("immune_"):
            continue
        events = cond["measurements"]["ferroptotic_events"]
        death_lps = [e["death_lp"] for e in events if e.get("death_lp") is not None]
        completed = [e for e in events if e.get("release_lp") is not None]
        release_lps = [e["release_lp"] for e in completed]
        release_damps = [e["release_damp"] for e in completed]
        arms[name[len("immune_"):]] = {
            "deaths": len(events),
            "completed_releases": len(completed),
            "mean_death_lp": fmean(death_lps) if death_lps else None,
            "mean_release_lp": fmean(release_lps) if release_lps else None,
            "total_release_damp": sum(release_damps) if release_damps else 0.0,
        }
    return {"arms": arms, "archive": str(SPATIAL_ARCHIVE.relative_to(REPO))}


def measure_icd(binary: Path) -> dict:
    """Run `sim-icd` and keep its per-row per-death DAMP."""
    with tempfile.TemporaryDirectory(prefix="ferro_icd_") as wd:
        out = Path(wd) / "icd"
        proc = subprocess.run([str(binary), "--output-dir", str(out)],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"sim-icd failed ({proc.returncode}): {proc.stderr[-400:]}")
        rows = json.loads((out / "icd_comparison.json").read_text())
    keep = []
    for r in rows:
        keep.append({
            "phenotype": r["phenotype"],
            "treatment": r["treatment"],
            "n_dead": r["n_dead"],
            "death_rate": r["death_rate"],
            "avg_final_lp": r.get("avg_final_lp"),
            "damp_per_dead_cell": r["immune_no_pd1"]["damp_per_dead_cell"],
            # P5's threshold names TWO readouts, so both are kept. `mature_dcs`
            # divided by the kill count is the DC-maturation-per-dead-cell the
            # prediction registers; `dc_activation_fraction` is the saturating
            # DAMP response it is computed from.
            "dc_activation_fraction": r["immune_no_pd1"]["dc_activation_fraction"],
            "mature_dcs": r["immune_no_pd1"]["mature_dcs"],
            "immune_kills_with_pd1": r["immune_with_pd1"]["immune_kills"],
        })
    return {"rows": keep}


def measure_panel(binary: Path) -> dict:
    """Run the panel's TME sweep and keep its per-death DAMP and denominators."""
    with tempfile.TemporaryDirectory(prefix="ferro_panel_") as wd:
        out = Path(wd) / "panel"
        proc = subprocess.run(
            [str(binary), "--tme-sweep", "--output-dir", str(out)],
            capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"sim-modality-panel failed ({proc.returncode}): {proc.stderr[-400:]}")
        sweep = json.loads((out / "modality_tme_sweep.json").read_text())
    keep = []
    for r in sweep["conditions"]:
        keep.append({
            "phenotype": r["phenotype"],
            "hypoxic": r["hypoxic"],
            "stroma": r["stroma"],
            "acidic": r["acidic"],
            "heterogeneous": r["heterogeneous"],
            "damp_per_death": {a: r["damp_per_death"].get(a) for a in ARMS},
            "deaths": {a: r["deaths"].get(a) for a in ARMS},
        })
    return {"rows": keep}


def assemble(raw: dict) -> dict:
    """Derive every reported field from the measurements.

    Idempotent, so `--render-only` reproduces a fresh run exactly instead of
    echoing stored prose -- re-rendering a stored derived field would make every
    guard on it inert.
    """
    out = {"icd": raw["icd"], "panel": raw["panel"], "spatial": raw["spatial"]}
    thr = death_threshold()

    # --- the three regimes, each reduced to SDT:RSL3 ratios ---
    ratios: list[dict] = []

    sp = raw["spatial"]["arms"]
    sdt, rsl3 = sp.get("SDT"), sp.get("RSL3")
    spatial_release = spatial_at_death = None
    if sdt and rsl3 and sdt["mean_release_lp"] and rsl3["mean_release_lp"]:
        spatial_release = sdt["mean_release_lp"] / rsl3["mean_release_lp"]
        ratios.append({"regime": "2D spatial, release LP", "ratio": spatial_release,
                       "sdt": sdt["mean_release_lp"], "rsl3": rsl3["mean_release_lp"]})
    if sdt and rsl3 and sdt["mean_death_lp"] and rsl3["mean_death_lp"]:
        spatial_at_death = sdt["mean_death_lp"] / rsl3["mean_death_lp"]

    icd_by_key = {(r["phenotype"], r["treatment"]): r for r in raw["icd"]["rows"]}
    phenos = sorted({p for p, _ in icd_by_key})
    for pheno in phenos:
        s = icd_by_key.get((pheno, "SDT"))
        r = icd_by_key.get((pheno, "RSL3"))
        if not s or not r:
            continue
        if not s.get("damp_per_dead_cell") or not r.get("damp_per_dead_cell"):
            continue
        ratios.append({"regime": f"sim-icd, {pheno}",
                       "ratio": s["damp_per_dead_cell"] / r["damp_per_dead_cell"],
                       "sdt": s["damp_per_dead_cell"], "rsl3": r["damp_per_dead_cell"]})

    panel_rows = raw["panel"]["rows"]
    panel_ratios = []
    for row in panel_rows:
        q = row["damp_per_death"]
        if q.get("SDT") and q.get("RSL3"):
            panel_ratios.append(q["SDT"] / q["RSL3"])
    if panel_ratios:
        ratios.append({"regime": f"panel sweep, {len(panel_ratios)} readable rows",
                       "ratio": max(panel_ratios), "sdt": None, "rsl3": None,
                       "span": [min(panel_ratios), max(panel_ratios)]})

    out["ratios"] = ratios
    # COUNTED ONCE EACH. `ratios` is the DISPLAY table, whose panel entry
    # summarises 32 rows, so summing it with `panel_ratios` double-counts the
    # sweep and the first version of this page reported 36 comparisons where
    # there are 35.
    all_vals = [r["ratio"] for r in ratios if "span" not in r] + panel_ratios
    out["ratio_min"] = min(all_vals) if all_vals else None
    out["ratio_max"] = max(all_vals) if all_vals else None
    out["panel_readable_rows"] = len(panel_ratios)
    out["panel_total_rows"] = len(panel_rows)

    # --- the direction, separately from the magnitude ---
    out["readable_comparisons"] = len(all_vals)
    out["direction_holds_everywhere"] = bool(all_vals) and all(v > 1.0 for v in all_vals)

    # --- the predictions ---
    preds = []
    for p in PREDICTIONS:
        reached = out["ratio_max"] is not None and out["ratio_max"] >= p["low"]
        preds.append({**p, "reached": bool(reached),
                      "shortfall_factor": (p["low"] / out["ratio_max"])
                      if out["ratio_max"] else None})
    out["predictions"] = preds
    out["every_prediction_refuted"] = bool(preds) and not any(p["reached"] for p in preds)

    # --- the mechanism: is ONE arm near threshold and the other far? ---
    out["death_threshold"] = thr
    out["at_death_ratio"] = spatial_at_death
    out["release_ratio"] = spatial_release
    near = {}
    for arm in ("SDT", "RSL3"):
        a = sp.get(arm) or {}
        lp = a.get("mean_death_lp")
        near[arm] = None if lp is None else bool(
            abs(lp - thr) / thr <= NEAR_THRESHOLD_TOLERANCE)
    out["at_death_near_threshold"] = near
    out["both_arms_near_threshold_at_death"] = all(v is True for v in near.values())
    out["at_death_contrast_exists"] = (
        len({v for v in near.values() if v is not None}) == 2)

    # --- the two sentinel/ranking defects ---
    out["panel_rows_with_absent_cohort"] = sum(
        1 for row in panel_rows
        for a in ARMS
        if row["damp_per_death"].get(a) is None)
    out["absent_cohorts_all_have_zero_deaths"] = all(
        row["deaths"].get(a) == 0
        for row in panel_rows for a in ARMS
        if row["damp_per_death"].get(a) is None)
    adc = [row for row in panel_rows
           if row["damp_per_death"].get("AntibodyDrugConjugate")
           and row["damp_per_death"].get("SDT")
           and row["damp_per_death"]["AntibodyDrugConjugate"]
           > row["damp_per_death"]["SDT"]]
    out["adc_outranks_sdt_rows"] = len(adc)
    out["adc_outranks_sdt_is_hypoxic"] = bool(adc) and all(r["hypoxic"] for r in adc)

    # --- the saturation trap, measured rather than recalled ---
    sat = [r for r in raw["icd"]["rows"] if r["death_rate"] > 0.99]
    out["icd_saturated_rows"] = len(sat)
    pers_sdt = icd_by_key.get(("Persister", "SDT"))
    pers_rsl3 = icd_by_key.get(("Persister", "RSL3"))
    if pers_sdt and pers_rsl3:
        out["persister_immune_kills"] = {
            "SDT": pers_sdt["immune_kills_with_pd1"],
            "RSL3": pers_rsl3["immune_kills_with_pd1"],
        }
        out["immune_kill_ranking_inverts_under_saturation"] = bool(
            pers_sdt["death_rate"] > 0.99
            and pers_sdt["immune_kills_with_pd1"] < pers_rsl3["immune_kills_with_pd1"])
    else:
        out["persister_immune_kills"] = None
        out["immune_kill_ranking_inverts_under_saturation"] = False

    # --- P5, whose registered falsification threshold is ON this quantity ---
    # NOT an amendment to the registration, and the distinction matters: the
    # threshold stays exactly as registered and the model is measured against
    # it. `PREREGISTRATION.md`'s P5 states the model output as a per-cell ratio
    # "about 4:1 in 3D (down from the 2D 104:1 ceiling)", which the manuscript
    # has already recorded as a mislabel -- those are total immune-kill ratios.
    # What was never measured is the per-cell output itself.
    band_m = P5_BAND_RE.search(PREREG.read_text(encoding="utf-8"))
    if band_m is None:
        raise RuntimeError(
            f"P5's registered falsification band is not in {PREREG}; this page "
            "measures the model against it and must not supply its own number")
    band = float(band_m.group(1))
    legs = []
    for pheno in phenos:
        s_row = icd_by_key.get((pheno, "SDT"))
        r_row = icd_by_key.get((pheno, "RSL3"))
        if not s_row or not r_row or not s_row["n_dead"] or not r_row["n_dead"]:
            continue
        if not s_row.get("damp_per_dead_cell") or not r_row.get("damp_per_dead_cell"):
            continue
        damp = s_row["damp_per_dead_cell"] / r_row["damp_per_dead_cell"]
        dc_s = s_row["mature_dcs"] / s_row["n_dead"]
        dc_r = r_row["mature_dcs"] / r_row["n_dead"]
        # TWO NORMALISATIONS OF THE SAME READOUT, kept apart on purpose.
        # `mature_dcs / kills` and `dc_activation_fraction` must give the same
        # ratio, because the cascade's mature-DC count is proportional to the
        # activation fraction times the kill count. If they disagreed, the
        # quantity would not be per-cell at all, and the page would be
        # measuring P5's readout with the wrong instrument.
        act = (s_row["dc_activation_fraction"] / r_row["dc_activation_fraction"]
               if r_row["dc_activation_fraction"] else None)
        legs.append({"phenotype": pheno, "damp_per_death_ratio": damp,
                     "dc_maturation_per_death_ratio": dc_s / dc_r,
                     "dc_activation_ratio": act})
    out["p5"] = {
        "band": band,
        "quote": P5_QUOTE,
        "legs": legs,
        # The registration falsifies if BOTH readouts agree within the band.
        "both_legs_inside_band": bool(legs) and all(
            lg["damp_per_death_ratio"] < band
            and lg["dc_maturation_per_death_ratio"] < band for lg in legs),
        "dc_normalisations_agree": bool(legs) and all(
            lg["dc_activation_ratio"] is not None
            and abs(lg["dc_maturation_per_death_ratio"] - lg["dc_activation_ratio"])
            < 1e-5 for lg in legs),
        "worst_leg": max(
            (max(lg["damp_per_death_ratio"], lg["dc_maturation_per_death_ratio"])
             for lg in legs), default=None),
    }

    # --- verdict ---
    if out["ratio_max"] is None:
        out["verdict"] = "not_measurable"
    elif out["every_prediction_refuted"] and out["direction_holds_everywhere"]:
        out["verdict"] = "direction_holds_magnitude_refuted"
    elif out["direction_holds_everywhere"]:
        out["verdict"] = "direction_and_magnitude_hold"
    else:
        out["verdict"] = "direction_does_not_hold"
    return out


def _fmt(v, spec=".4f"):
    return "n/a" if v is None else format(v, spec)


def render(d: dict) -> str:
    thr = d["death_threshold"]
    sp = d["spatial"]["arms"]
    L: list[str] = []
    L.append("# The per-dead-cell DAMP advantage, measured\n")
    L.append(
        "Generated by `scripts/icd_damp_per_death.py`. The engine documented a "
        "per-dead-cell DAMP advantage for physical ROS death over GPX4 "
        "inhibition, and two committed documents put a number on it. This page "
        "measures it in every regime the repository can run.\n")

    L.append("## The documented predictions\n")
    L.append("| source | prediction | reached? |")
    L.append("|---|---|---|")
    for p in d["predictions"]:
        L.append(f"| `{p['source']}` ({p['context']}) | {p['low']:g}-{p['high']:g}x | "
                 f"{'**yes**' if p['reached'] else 'no'} |")
    L.append("")

    L.append("## Measured\n")
    L.append("| regime | SDT | RSL3 | SDT:RSL3 |")
    L.append("|---|--:|--:|--:|")
    for r in d["ratios"]:
        span = r.get("span")
        ratio = (f"{span[0]:.4f} - {span[1]:.4f}" if span else f"**{r['ratio']:.4f}**")
        L.append(f"| {r['regime']} | {_fmt(r.get('sdt'))} | {_fmt(r.get('rsl3'))} | "
                 f"{ratio} |")
    L.append("")
    if d["ratio_max"] is not None:
        worst = min(p["shortfall_factor"] for p in d["predictions"]
                    if p["shortfall_factor"] is not None)
        L.append(
            f"Across {d['readable_comparisons']} readable comparisons the ratio spans "
            f"**{d['ratio_min']:.4f}x to {d['ratio_max']:.4f}x**. The lowest "
            f"documented prediction is {min(p['low'] for p in d['predictions']):g}x, "
            f"so the best measurement falls short of the weakest claim by a factor "
            f"of {worst:.2f}.\n")

    L.append("## The direction, which is not what fails\n")
    if d["direction_holds_everywhere"]:
        L.append(
            f"SDT exceeds RSL3 in **all {d['readable_comparisons']}** readable "
            "comparisons. The documented direction survives; the documented "
            "magnitude does not, and this page withdraws only the second.\n")
    else:
        L.append(
            "**The direction does NOT hold in every readable comparison**, so more "
            "than the magnitude is in question and the rows above should be read "
            "before anything else here.\n")

    L.append("## The mechanism names a moment where the contrast does not exist\n")
    L.append("| arm | mean LP at death | mean LP at release | deaths |")
    L.append("|---|--:|--:|--:|")
    for arm in sorted(sp):
        a = sp[arm]
        L.append(f"| {arm} | {_fmt(a['mean_death_lp'])} | "
                 f"{_fmt(a['mean_release_lp'])} | {a['deaths']:,} |")
    L.append("")
    L.append(
        f"`death_threshold` is {thr:g}, read from `params.rs` rather than quoted. "
        f"At the moment of death the two arms differ by a factor of "
        f"{_fmt(d['at_death_ratio'])}; after the post-death grace window they "
        f"differ by {_fmt(d['release_ratio'])}.\n")
    if d["both_arms_near_threshold_at_death"]:
        L.append(
            f"**Both arms are within {NEAR_THRESHOLD_TOLERANCE:.0%} of the "
            "threshold at death**, so a contrast between one arm 'MUCH higher "
            "than threshold' and another 'NEAR threshold' cannot be read at that "
            "moment: death IS the threshold crossing. The documented contrast "
            "compares one arm at death with the other after it.\n")
    elif d["at_death_contrast_exists"]:
        L.append(
            "One arm is near the threshold at death and the other is not, so the "
            "documented contrast IS readable at that moment and only its "
            "magnitude is at issue.\n")

    L.append("## Two defects the measurement surfaced\n")
    L.append(
        f"**An absent cohort is not a potency of zero.** "
        f"{d['panel_rows_with_absent_cohort']} of "
        f"{d['panel_total_rows'] * len(ARMS)} arm-rows in the panel sweep have no "
        "cohort to average over")
    if d["absent_cohorts_all_have_zero_deaths"]:
        L.append(
            "  -- every one of them at a kill count of exactly zero, confirmed "
            "against the `deaths` denominator now emitted beside the value. They "
            "are reported as `null`; they used to be reported as `0.0`, which "
            "under a column headed DAMP-per-death reads as a measurement.\n")
    else:
        L.append(
            "  -- but NOT all of them at a kill count of zero, so `null` is "
            "hiding something other than an absent cohort and that must be "
            "explained before the table is read.\n")
    if d["adc_outranks_sdt_rows"]:
        L.append(
            f"**The ROS cascade is not top of the per-death ranking.** "
            f"`AntibodyDrugConjugate` releases more per death than SDT in "
            f"{d['adc_outranks_sdt_rows']} of {d['panel_total_rows']} sweep rows"
            + (", all of them hypoxic" if d["adc_outranks_sdt_is_hypoxic"] else "")
            + " -- a cytotoxic payload out-releasing the cascade wherever oxygen "
              "limits it. So 'quality of death' does not order the arms by death "
              "mode the way the documented sentence implies.\n")

    if d["immune_kill_ranking_inverts_under_saturation"]:
        k = d["persister_immune_kills"]
        L.append("## A companion claim that inverts for a structural reason\n")
        L.append(
            f"`sim-icd`'s README also stated that SDT+anti-PD-1 immune kills far "
            f"exceed RSL3's. Measured: SDT {k['SDT']:,.0f} against RSL3 "
            f"{k['RSL3']:,.0f}. `immune_cascade` caps its kills at the SURVIVING "
            f"population, and {d['icd_saturated_rows']} of the binary's rows kill "
            "over 99% of cells biochemically, so the arm that kills everything has "
            "nothing left to earn immune kills on. That comparison ranks survivors, "
            "not immunogenicity.\n")

    p5 = d["p5"]
    if p5["legs"]:
        L.append("## The registered prediction this lands on\n")
        L.append(
            f"`PREREGISTRATION.md`'s **P5** registers a falsification threshold "
            f"on exactly this quantity: *\"{p5['quote']}\"*. Both readouts, "
            "measured:\n")
        L.append("| phenotype | DAMP per dead cell | DC maturation per dead cell |")
        L.append("|---|--:|--:|")
        for lg in p5["legs"]:
            L.append(f"| {lg['phenotype']} | {lg['damp_per_death_ratio']:.4f}x | "
                     f"{lg['dc_maturation_per_death_ratio']:.4f}x |")
        L.append("")
        if p5["both_legs_inside_band"]:
            L.append(
                f"**Both readouts agree within {p5['band']:g}x in every readable "
                f"phenotype, the largest being {p5['worst_leg']:.4f}x -- so the "
                "model's own per-cell output satisfies P5's registered "
                "falsification condition.** Read that precisely: P5 is a "
                "prediction about BIOLOGY, to be tested by experiment E3, and "
                "nothing here falsifies the biology. What it establishes is that "
                "the model does not predict what P5 says it predicts. P5's "
                "\"Quantitative model output\" line states a per-cell ratio of "
                "\"about 4:1 in 3D (down from the 2D 104:1 ceiling)\" and the "
                "manuscript already records those as total immune-kill ratios "
                "rather than per-cell readouts; measured, the per-cell output is "
                "a small POSITIVE difference that falls inside the band P5 "
                "registers as falsification -- not a null, but not "
                "distinguishable from one at P5's own resolution.\n")
            L.append(
                "**The threshold is not amended and must not be.** It is parsed "
                "from the registration here rather than restated, so this page "
                "tracks the registered band instead of carrying its own. "
                "Amending a registered threshold after measuring against it is "
                "the thing preregistration exists to prevent; the manuscript "
                "declined to amend it for the same reason, and so does this.\n")
            L.append(
                "**Which leg carries this, stated rather than assumed.** P5's "
                "contrast is DENSE versus SPARSE kill, which is a spatial "
                "property, and `sim-icd` has no geometry at all -- it is a "
                "population of independent cells, so its SDT-versus-RSL3 "
                "contrast is a death-MODE contrast that the registration "
                "identifies with the dense and sparse arms by name but not by "
                "mechanism. The leg that genuinely measures dense against "
                f"sparse is the 2D spatial archive, where SDT kills "
                f"{sp['SDT']['deaths']:,} cells and RSL3 kills "
                f"{sp['RSL3']['deaths']:,} in the SAME geometry, and it gives "
                f"{d['release_ratio']:.4f}x. So the conclusion rests on the leg "
                "whose mechanism matches the registration, and the other two "
                "agree with it.\n")
            L.append(
                "The DC leg is measured through two normalisations that must "
                "agree: `mature_dcs` divided by the kill count, and "
                "`dc_activation_fraction`. "
                + ("They give the SAME ratio to five decimals, because the "
                   "cascade's mature-DC count is proportional to the activation "
                   "fraction times the kill count -- which is the check that "
                   "the quantity is per-cell at all."
                   if p5["dc_normalisations_agree"] else
                   "**They DISAGREE**, so the DC figure above is not a per-cell "
                   "readout and P5's DC leg is not being measured here; nothing "
                   "should be concluded from it until that is explained.")
                + "\n")
        else:
            L.append(
                f"At least one readout exceeds the registered {p5['band']:g}x "
                f"band (largest {p5['worst_leg']:.4f}x), so the model's per-cell "
                "output is distinguishable from the null P5 would be falsified "
                "by, and P5's per-cell leg remains a model prediction.\n")

    L.append("## Scope\n")
    L.append(
        f"- **{d['panel_readable_rows']} of the {d['readable_comparisons']} "
        "comparisons are panel-sweep CONDITIONS, not independent replicates.** "
        "They vary phenotype, depth, hypoxia, stroma and pH, but they all run "
        "the same single-cell biochemistry through the same "
        "`calculate_damp_release`, as does `sim-icd`. The committed 2D spatial "
        "archive is the one row on a different code path. So the count measures "
        "how little the ratio moves across the conditions the engine can pose, "
        "not how many separate times the claim was tested.\n"
        "- `damp_per_lp` is 1.0 and uncalibrated, so every ratio here is a "
        "property of the model and not a measurement of immunogenicity.\n"
        "- Nothing is recalibrated and no parameter is changed. What changed is "
        "what the engine's comments and the binary's README say about the "
        "measurement.\n"
        "- The manuscript already declines the per-cell reading -- Section 8.2 "
        "calls the linear LP-to-DAMP assumption \"not a measurement of per-cell "
        "immunogenicity\" and Figure 16's caption says the 104:1 \"compares total "
        "immune-kill counts, not DAMP potency or immunogenicity per dead cell\" -- "
        "so no manuscript number rests on the refuted band, and the manuscript is "
        "not edited.\n"
        "- The totals DO differ by roughly three orders of magnitude between these "
        "two arms, and that is the kill count, which is the quantity the "
        "manuscript reports.\n")
    L.append(f"Verdict: `{d['verdict']}`.\n")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the prose from the committed JSON, no binaries")
    ap.add_argument("--icd", type=Path, default=None)
    ap.add_argument("--panel", type=Path, default=None)
    a = ap.parse_args()

    # `measure_spatial` is called on BOTH paths, deliberately: it reads a
    # committed archive rather than a binary, so there is nothing to gain by
    # trusting a stored copy of it and a leg that is always recomputed can
    # never be echoed.
    if a.render_only:
        # RE-ASSEMBLED, not echoed: the stored file carries derived fields next
        # to the raw measurements, and re-rendering those instead of recomputing
        # them would make every guard on a derived field inert.
        stored = json.loads(OUT_JSON.read_text())
        data = assemble({"icd": stored["icd"], "panel": stored["panel"],
                         "spatial": measure_spatial()})
    else:
        icd = a.icd or _default_binary("sim-icd")
        panel = a.panel or _default_binary("sim-modality-panel")
        for b, crate in ((icd, "sim-icd"), (panel, "sim-modality-panel")):
            if not b.exists():
                print(f"error: {crate} not built at {b}; "
                      f"cargo build --release -p {crate}", file=sys.stderr)
                return 2
        data = assemble({"icd": measure_icd(icd), "panel": measure_panel(panel),
                         "spatial": measure_spatial()})

    OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
