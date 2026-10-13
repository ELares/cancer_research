#!/usr/bin/env python3
"""The penetration headline is measured where the penetration length cancels (#926).

`headline_uncertainty.py --headline penetration` publishes the manuscript's
per-tissue figure -- 40% in 2D culture dropping to 12.1% / 2.6% / 1.8% at the
vessel wall -- and reports that "the penetration-gradient DIRECTION is therefore
parameter-robust" because the well >= poorly >= CNS ordering held in 300 of 300
prior draws. `identifiability_report.py` carries that forward as the
`penetration_gap` verdict `direction_robust_magnitude_not`, rationale "the
penetration-gradient ordering is parameter-robust". The manuscript repeats the
sentence at Chapter 5.

TWO THINGS ARE WRONG WITH THAT LABEL, AND BOTH ARE MEASURED HERE.

(1) THE OBSERVABLE CONTAINS NO PENETRATION LENGTH. `drug_transport` computes
    `C(r) = vessel_wall_conc * vascular_permeability * exp(-r / lambda)`, and the
    headline reads bin 0, where `r = 0` exactly, so `exp(0) = 1` for EVERY
    lambda. The quantity named a penetration gradient is a `vascular_permeability`
    gradient. `sim-tissue-pk` already ships the control: it runs the same three
    tissues for two drugs whose penetration lengths differ two-fold (RSL3-like
    100.0 um, the doxorubicin reference 50.0 um), and the vessel-wall kills come
    back BIT-IDENTICAL -- 0.121 / 0.026 / 0.018 under both. What DOES move with
    lambda is the kill depth (14.7 -> 3.7 um) and the overall kill fraction
    (0.0658 -> 0.0461), neither of which is the headline.

    This matters because of which parameter carries the evidence. #335 validated
    `lambda` against measured data (Primeau 2005, Tannock 2002) and
    CALIBRATION_STATUS records that validation; `vascular_permeability` is an
    uncalibrated preset. The headline attributes a validated-sounding property to
    the knob that has no anchor, and does not use the knob that has one.

(2) THE 300/300 IS ENTAILED, NOT MEASURED. `cell_seed` in sim-tissue-pk is
    `seed + bin_idx * N_CELLS * 2 + i * 2` -- it does not depend on the tissue --
    so all three tissues run the IDENTICAL cells with the IDENTICAL RNG streams,
    and at bin 0 they differ in nothing but `c0`. The three headline numbers are
    therefore ONE dose-response evaluated at three points, which this file
    verifies by reproducing each of them EXACTLY from the real binding
    (`sim_batch('Persister','RSL3', n=1000, seed=42)`). A paired comparison of one
    monotone function at a fixed decreasing input sequence cannot invert, so the
    ordering is robust to the rate constants by construction and the draw count
    is not evidence about them. The committed report states the mechanism in the
    same breath as the conclusion -- "a shared monotone dose-response, so the
    across-tissue ORDERING is preserved" -- and then reads the consequence as
    parameter-robustness.

AND WHERE THE PENETRATION LENGTH DOES ACT, THE ORDERING IS NOT ROBUST. The three
presets disagree about which tissue is hardest to treat, because access and
distance run opposite ways for one pair: the inter-vessel HALF-distances are
well 60 um, CNS 75 um, poorly 125 um, so by diffusion distance the CNS tumour is
better perfused than the poorly-vascularized one, the reverse of the published
ordering. At each tissue's own deepest cell the two effects compete, and the
poorly/CNS concentration ordering INVERTS below a crossover this file derives
from the stored presets and then checks against the engine: the two shipped
lambdas bracket it, the 50 um arm inverted and the 100 um arm not. The crossover
lands inside Tannock 2002's measured 25-75 um band, while the 100 um RSL3-like
value that keeps the published ordering is the one CALIBRATION_STATUS flags as
UNVALIDATED ("no ferroptosis inducer has published spatial penetration data").

WHAT IS NOT CLAIMED, because overstating this would be the easy error:
  - The ordering claim is not FALSE. It is robust to the 11 rate constants, and
    measured here it holds in every draw. The defect is that it is robust by
    construction and that it describes vascular access, not penetration.
  - The deep-end inversion is in the CONCENTRATION. The kills at those depths are
    indistinguishable (0.009 in both arms at lambda 50) because both sit far below
    the bistable threshold, so the inversion does NOT flip a kill ordering, and
    this file reports the concentrations rather than dressing them as kills.
  - The MANUSCRIPT'S OWN SECTION 6.1 TREATMENT IS CORRECT and is left alone: it
    attributes the numbers to "tissue-specific vascular permeability", names the
    vessel wall as "the highest-concentration position", and says the transport
    parameters are estimated. Only the Chapter 5 sentence inherited from the
    generator is wrong. A sibling issue once accused this manuscript of an
    over-claim it had already caveated; the check came first this time.
  - Nothing is recalibrated. `vascular_permeability` has no published anchor in
    this repo and inventing one would be worse than naming it uncalibrated.

OFFLINE CONTRACT: the measurement needs the compiled `sim-tissue-pk` and the
compiled Python extension, so it is taken locally and the derived artifact
committed; CI reads only the JSON, and `--render-only` rebuilds the prose
idempotently from it.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from headline_sensitivity import (  # noqa: E402
    PENETRATION_TISSUES,
    extract_tissue_pk_observables,
)

OUT_JSON = REPO / "analysis" / "calibration" / "penetration-access-vs-depth.json"
OUT_MD = REPO / "analysis" / "calibration" / "penetration-access-vs-depth.md"

PRCC = REPO / "analysis" / "prcc-results.json"

# The dose-response grid. `f(c)` is bin 0 of sim-tissue-pk: the SAME 1000 cells
# at the SAME seeds, with GPX4 inhibition scaled by `c`. 0.0 is included because
# it is the no-drug control and the CNS headline sits a few cells above it.
DOSE_GRID = [i / 20 for i in range(21)]

# Draws used to re-measure the ordering claim. The committed report used 300; the
# point here is not the count but that the outcome is entailed, so a smaller
# ensemble that reproduces 100% makes the same point far faster.
N_DRAWS = 200
DRAW_SEED = 332

# Tannock 2002 (PMID 11895922) measured 25-75 um for doxorubicin in vivo -- the
# band #335 validated the reference lambda against. Used ONLY to say whether the
# derived crossover falls inside it.
TANNOCK_BAND_UM = (25.0, 75.0)

PHENOTYPE = "Persister"
BIN0_SEED = 42
BIN0_CELLS = 1000


def default_binary() -> Path:
    env = os.environ.get("FERRO_SIM_TISSUE_PK")
    if env:
        return Path(env)
    out = subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--no-deps"],
        cwd=REPO / "simulations", capture_output=True, text=True,
    )
    if out.returncode == 0:
        target = Path(json.loads(out.stdout)["target_directory"])
        return target / "release" / "sim-tissue-pk"
    return REPO / "simulations" / "target" / "release" / "sim-tissue-pk"


def run_binary(binary: Path) -> dict:
    """One unperturbed sim-tissue-pk run. Returns the per-drug/per-tissue wall
    kills, the lambdas, and each tissue's deepest-bin radius and concentration.

    The RSL3 arm is read through the repo's own `extract_tissue_pk_observables`
    so the headline observable cannot drift away from the one the interval used;
    the second drug is read from the same summary, since that extractor filters
    to RSL3 by design."""
    with tempfile.TemporaryDirectory(prefix="pen_access_") as wd:
        env = dict(os.environ)
        env.pop("FERRO_PARAM_OVERRIDES", None)
        proc = subprocess.run([str(binary)], cwd=wd, env=env,
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"sim-tissue-pk failed: {proc.stderr[-400:]}")
        out = Path(wd) / "output" / "tissue-pk"
        summaries = json.loads((out / "tissue_pk_summary.json").read_text())
        rows = list(csv.DictReader((out / "tissue_pk_results.csv").open()))

    # The headline arm, through the shared extractor.
    rsl3_wall = extract_tissue_pk_observables(summaries)

    wall, lam, deepest, half_dist, c0 = {}, {}, {}, {}, {}
    for s in summaries:
        drug, tissue = s["drug"], s["tissue"]
        wall.setdefault(drug, {})[tissue] = s["vessel_wall_death_rate"]
        lam[drug] = s["penetration_length_um"]
        half_dist[tissue] = s["max_distance_um"]
        c0[tissue] = s["vessel_wall_concentration"]
    for drug in wall:
        deepest[drug] = {}
        for tissue in wall[drug]:
            sel = [r for r in rows if r["drug"] == drug and r["tissue"] == tissue]
            far = max(sel, key=lambda r: float(r["distance_um"]))
            deepest[drug][tissue] = {
                "distance_um": float(far["distance_um"]),
                "concentration": float(far["concentration"]),
                "death_rate": float(far["death_rate"]),
            }
    return {
        "wall_death_rate": wall,
        "penetration_length_um": lam,
        "deepest_bin": deepest,
        "half_distance_um": half_dist,
        "vessel_wall_concentration": c0,
        "rsl3_wall_via_shared_extractor": rsl3_wall,
    }


def measure_dose_response() -> dict:
    """`f(c)` on a grid, and the no-drug control, from the real binding."""
    import ferroptosis_core as fc

    base = fc.default_params()["rsl3_gpx4_inhib"]
    grid = {}
    for c in DOSE_GRID:
        r = fc.sim_batch(PHENOTYPE, "RSL3", n=BIN0_CELLS, seed=BIN0_SEED,
                         context="2d", rsl3_gpx4_inhib=base * c)
        grid[f"{c:g}"] = r["death_rate"]
    control = fc.sim_batch(PHENOTYPE, "Control", n=BIN0_CELLS, seed=BIN0_SEED,
                           context="2d")["death_rate"]
    return {"default_rsl3_gpx4_inhib": base, "grid": grid, "control": control}


def measure_prior_ordering() -> dict:
    """Re-measure the ordering claim over the same 11-parameter prior, and ALSO
    ask the question the committed report does not: in how many draws is the CNS
    kill even distinguishable from that draw's own no-drug control?"""
    import random

    import ferroptosis_core as fc

    ranges = json.loads(PRCC.read_text())["metadata"]["parameter_ranges"]
    names = list(ranges)
    concs = {k: c for k, c in zip(
        [k for k, _ in PENETRATION_TISSUES], (0.8, 0.4, 0.15))}
    rng = random.Random(DRAW_SEED)
    n_ord = n_cns = n_poorly = 0
    for _ in range(N_DRAWS):
        draw = {n: rng.uniform(*ranges[n]) for n in names}
        inhib = draw["rsl3_gpx4_inhib"]
        f = {}
        for key, c in concs.items():
            p = dict(draw)
            p["rsl3_gpx4_inhib"] = inhib * c
            f[key] = fc.sim_batch(PHENOTYPE, "RSL3", n=BIN0_CELLS,
                                  seed=BIN0_SEED, context="2d", **p)["death_rate"]
        ctrl = fc.sim_batch(PHENOTYPE, "Control", n=BIN0_CELLS, seed=BIN0_SEED,
                            context="2d", **draw)["death_rate"]
        keys = [k for k, _ in PENETRATION_TISSUES]
        if f[keys[0]] >= f[keys[1]] >= f[keys[2]]:
            n_ord += 1
        if f[keys[2]] > ctrl:
            n_cns += 1
        if f[keys[1]] > ctrl:
            n_poorly += 1
    return {
        "n_draws": N_DRAWS,
        "seed": DRAW_SEED,
        "ordering_preserved": n_ord,
        "cns_above_own_control": n_cns,
        "poorly_above_own_control": n_poorly,
    }


def _tissue_of(keyword: str, names) -> str:
    hit = [n for n in names if keyword.lower() in n.lower()]
    if len(hit) != 1:
        raise RuntimeError(f"expected exactly 1 tissue matching {keyword!r}, got {hit}")
    return hit[0]


def assemble(raw: dict) -> dict:
    """Derive every published claim from the stored measurements. Idempotent:
    re-running it on its own output reproduces the same dict, so `--render-only`
    re-derives rather than echoing (a derived field that is echoed makes every
    guard on it inert)."""
    wall = raw["wall_death_rate"]
    lam = raw["penetration_length_um"]
    deep = raw["deepest_bin"]
    half = raw["half_distance_um"]
    c0 = raw["vessel_wall_concentration"]
    dr = raw["dose_response"]
    prior = raw["prior_ordering"]

    tissues = sorted(c0, key=lambda t: -c0[t])
    drugs = sorted(lam, key=lambda d: -lam[d])

    # (1) lambda cancels at the wall: every tissue's wall kill is identical
    # across drugs whose lambdas differ. Compared EXACTLY -- the claim is
    # cancellation, and a tolerance would accept a small real difference.
    wall_identical = {
        t: len({wall[d][t] for d in drugs}) == 1 for t in tissues
    }
    lam_ratio = max(lam.values()) / min(lam.values())

    # (2) the headline is one dose-response: the binding reproduces each wall
    # kill exactly at the matching concentration.
    grid = dr["grid"]
    reproduced = {}
    for t in tissues:
        key = f"{c0[t]:g}"
        reproduced[t] = {
            "binary": wall[drugs[0]][t],
            "binding": grid.get(key),
            "exact": grid.get(key) == wall[drugs[0]][t],
        }
    vals = [grid[f"{c:g}"] for c in DOSE_GRID]
    monotone = all(b >= a for a, b in zip(vals, vals[1:]))

    ordering_entailed = bool(
        monotone and prior["ordering_preserved"] == prior["n_draws"]
    )

    # (3) how much of each headline number is the no-drug control
    ctrl = dr["control"]
    attributable = {
        t: {
            "headline": wall[drugs[0]][t],
            "control": ctrl,
            "drug_attributable": wall[drugs[0]][t] - ctrl,
            "fold_over_control": (wall[drugs[0]][t] / ctrl) if ctrl else None,
        }
        for t in tissues
    }

    # (4) access and distance disagree for one pair
    access_order = tissues                                  # by c0, descending
    distance_order = sorted(half, key=lambda t: half[t])     # nearest vessel first
    orders_agree = access_order == distance_order

    poorly = _tissue_of("poorly", tissues)
    cns = _tissue_of("CNS", tissues)
    # c_p exp(-r_p/L) == c_c exp(-r_c/L)  =>  L* = (r_p - r_c) / ln(c_p / c_c)
    num = half[poorly] - half[cns]
    den = math.log(c0[poorly] / c0[cns])
    crossover = num / den if den else None
    inside_band = bool(
        crossover is not None and TANNOCK_BAND_UM[0] <= crossover <= TANNOCK_BAND_UM[1]
    )
    # CHECKED against the engine rather than trusted: the two shipped lambdas
    # must fall on opposite sides of the derived crossover, and the measured
    # deepest-cell concentrations must invert on exactly the side that predicts it.
    deep_inverts = {
        d: deep[d][cns]["concentration"] > deep[d][poorly]["concentration"]
        for d in drugs
    }
    crossover_brackets = bool(
        crossover is not None
        and all(deep_inverts[d] == (lam[d] < crossover) for d in drugs)
        and len(set(deep_inverts.values())) == 2
    )

    if ordering_entailed and all(wall_identical.values()):
        verdict = "direction_robust_by_construction_magnitude_not"
    elif not all(wall_identical.values()):
        verdict = "observable_depends_on_penetration_length"
    else:
        verdict = "direction_not_entailed"

    out = dict(raw)
    out["derived"] = {
        "wall_identical_across_drugs": wall_identical,
        "lambda_ratio": lam_ratio,
        "lambda_cancels_at_the_wall": all(wall_identical.values()),
        "wall_reproduced_by_binding": reproduced,
        "headline_is_one_dose_response": all(r["exact"] for r in reproduced.values()),
        "dose_response_monotone": monotone,
        "ordering_entailed_not_measured": ordering_entailed,
        "drug_attributable": attributable,
        "access_order": access_order,
        "distance_order": distance_order,
        "access_and_distance_agree": orders_agree,
        "deep_crossover_lambda_um": crossover,
        "crossover_inside_tannock_band": inside_band,
        "deep_ordering_inverts": deep_inverts,
        "crossover_bracketed_by_shipped_lambdas": crossover_brackets,
        "verdict": verdict,
    }
    return out


def render(d: dict) -> str:
    v = d["derived"]
    lam = d["penetration_length_um"]
    wall = d["wall_death_rate"]
    c0 = d["vessel_wall_concentration"]
    half = d["half_distance_um"]
    deep = d["deepest_bin"]
    prior = d["prior_ordering"]
    drugs = sorted(lam, key=lambda x: -lam[x])
    tissues = v["access_order"]

    L = ["# The penetration headline is measured where lambda cancels (#926)\n"]
    L.append(
        "Generated by `scripts/penetration_access_vs_depth.py`. The "
        "`penetration_gap` headline reads `vessel_wall_death_rate`, which "
        "`drug_transport` evaluates at `r = 0`, where "
        "`C(r) = vessel_wall_conc * vascular_permeability * exp(-r/lambda)` "
        "loses its `lambda` term exactly. This measures what that does to the "
        "published verdict label.\n"
    )

    L.append("## 1. Two drugs, two penetration lengths, one answer\n")
    L.append("| drug | lambda (um) | " + " | ".join(
        f"{t} wall kill" for t in tissues) + " |")
    L.append("|---|---|" + "---|" * len(tissues))
    for dg in drugs:
        L.append(f"| {dg} | {lam[dg]:.1f} | " + " | ".join(
            f"{wall[dg][t]:.6f}" for t in tissues) + " |")
    L.append("")
    if v["lambda_cancels_at_the_wall"]:
        L.append(
            f"The two arms differ {v['lambda_ratio']:.1f}-fold in penetration "
            "length and agree **bit-for-bit** at every tissue. The headline "
            "observable therefore carries no information about penetration; "
            "what orders the three numbers is `vascular_permeability` "
            f"({', '.join(f'{c0[t]:g}' for t in tissues)}).\n"
        )
    else:
        L.append(
            "The wall kills are NOT identical across the two lambdas, so the "
            "observable does depend on penetration length and the premise of "
            "this analysis does not hold.\n"
        )
    L.append(
        "This is the parameter split that matters: #335 validated `lambda` "
        "against measured data (Primeau 2005, Tannock 2002) and "
        "`CALIBRATION_STATUS.md` records it; `vascular_permeability` carries no "
        "anchor in this repository. The headline uses only the unanchored one.\n"
    )

    L.append("## 2. The three numbers are one dose-response at three points\n")
    L.append("| tissue | c0 | binary wall kill | `sim_batch` at c0 | exact |")
    L.append("|---|---|---|---|---|")
    for t in tissues:
        r = v["wall_reproduced_by_binding"][t]
        b = "n/a" if r["binding"] is None else f"{r['binding']:.6f}"
        L.append(f"| {t} | {c0[t]:g} | {r['binary']:.6f} | {b} | "
                 f"{'yes' if r['exact'] else 'NO'} |")
    L.append("")
    L.append(
        "`cell_seed` in `sim-tissue-pk` is "
        "`seed + bin_idx * N_CELLS * 2 + i * 2`, which does not depend on the "
        "tissue, so at bin 0 the three tissues run identical cells with "
        "identical RNG and differ in nothing but `c0`. Reproducing each wall "
        "kill exactly from `sim_batch(seed=42, n=1000)` confirms it.\n"
    )
    L.append(
        f"That dose-response is monotone non-decreasing across all "
        f"{len(DOSE_GRID)} grid points "
        f"(**{'yes' if v['dose_response_monotone'] else 'NO'}**), and the "
        f"well >= poorly >= CNS ordering held in "
        f"{prior['ordering_preserved']}/{prior['n_draws']} prior draws. "
        + ("Those two facts are the same fact: a paired comparison of one "
           "monotone function at a fixed decreasing input sequence cannot "
           "invert, so the ordering is robust to the rate constants **by "
           "construction** and the draw count is not evidence about them.\n"
           if v["ordering_entailed_not_measured"] else
           "The ordering is NOT entailed here, so the draw count does carry "
           "information.\n")
    )

    L.append("## 3. Most of the CNS number is the no-drug control\n")
    L.append("| tissue | headline | no-drug control | drug-attributable | fold over control |")
    L.append("|---|---|---|---|---|")
    for t in tissues:
        a = v["drug_attributable"][t]
        fold = "n/a" if a["fold_over_control"] is None else f"{a['fold_over_control']:.2f}x"
        L.append(f"| {t} | {a['headline']:.3f} | {a['control']:.3f} | "
                 f"{a['drug_attributable']:+.3f} | {fold} |")
    L.append("")
    L.append(
        f"`Control` and `RSL3` at zero inhibition are the same number "
        f"({v['drug_attributable'][tissues[0]]['control']:.3f}), so the "
        "right-hand column is what the drug adds. Read against its own "
        "control, the gradient is far steeper than the headline triple "
        "suggests at the well-vascularized end and far shallower at the CNS "
        f"end, where the drug contributes "
        f"{v['drug_attributable'][tissues[-1]]['drug_attributable']:+.3f}. In "
        f"{prior['n_draws'] - prior['cns_above_own_control']} of "
        f"{prior['n_draws']} prior draws the CNS kill is not even strictly "
        "above that draw's own control.\n"
    )

    L.append("## 4. Access and distance disagree about which tissue is worst\n")
    L.append("| tissue | vascular permeability (access) | inter-vessel half-distance um (distance) |")
    L.append("|---|---|---|")
    for t in tissues:
        L.append(f"| {t} | {c0[t]:g} | {half[t]:.1f} |")
    L.append("")
    L.append(
        "Ranked by access: "
        + " > ".join(v["access_order"])
        + ". Ranked by proximity to a vessel: "
        + " < ".join(v["distance_order"])
        + ". "
        + ("The two orders AGREE, so nothing turns on which one the headline "
           "measures.\n" if v["access_and_distance_agree"] else
           "The two orders DISAGREE, and the published ordering is the access "
           "one -- measured at the one radius where distance cannot act.\n")
    )
    if v["deep_crossover_lambda_um"] is not None:
        L.append("| drug | lambda (um) | " + " | ".join(
            f"deepest-cell conc, {t}" for t in (v['distance_order'][-1], v['access_order'][-1])
        ) + " | poorly/CNS order |")
        L.append("|---|---|---|---|---|")
        poorly = v["distance_order"][-1]
        cns = v["access_order"][-1]
        for dg in drugs:
            inv = v["deep_ordering_inverts"][dg]
            L.append(
                f"| {dg} | {lam[dg]:.1f} | "
                f"{deep[dg][poorly]['concentration']:.4f} | "
                f"{deep[dg][cns]['concentration']:.4f} | "
                f"{'CNS > poorly (INVERTED)' if inv else 'poorly > CNS'} |"
            )
        L.append("")
        L.append(
            f"Solving the two profiles for equality gives a crossover at "
            f"**lambda = {v['deep_crossover_lambda_um']:.1f} um**, derived from "
            "the stored presets and then CHECKED against the engine: the two "
            "shipped lambdas must fall on opposite sides of it and the measured "
            "concentrations must invert on exactly the predicted side "
            f"(**{'confirmed' if v['crossover_bracketed_by_shipped_lambdas'] else 'NOT confirmed'}**). "
            + (f"It falls inside Tannock 2002's measured "
               f"{TANNOCK_BAND_UM[0]:.0f}-{TANNOCK_BAND_UM[1]:.0f} um band, "
               "while the 100 um RSL3-like value that preserves the published "
               "ordering at depth is the one `CALIBRATION_STATUS.md` flags as "
               "unvalidated.\n"
               if v["crossover_inside_tannock_band"] else
               "It falls outside the measured band.\n")
        )
        L.append(
            "**The inversion is in the concentration, not the kill.** At those "
            "depths both arms sit far below the bistable threshold and their "
            f"kills are indistinguishable "
            f"({deep[drugs[-1]][poorly]['death_rate']:.3f} and "
            f"{deep[drugs[-1]][cns]['death_rate']:.3f} at lambda "
            f"{lam[drugs[-1]]:.0f}), so no kill ordering flips. What flips is "
            "the transport quantity the word 'penetration' refers to.\n"
        )

    L.append("## Verdict\n")
    L.append(f"`{v['verdict']}`\n")
    L.append(
        "The ordering claim is not false -- it is robust to the 11 rate "
        "constants, here in every draw. It is robust **by construction**, and "
        "it describes vascular access rather than penetration. "
        "`identifiability_report.py` derives the `penetration_gap` row from "
        "this artifact so the label cannot outrun the measurement again.\n"
    )

    L.append("## Scope\n")
    L.append(
        "- Nothing is recalibrated. `vascular_permeability` has no published "
        "anchor in this repository, and supplying one from memory would be "
        "worse than naming it uncalibrated.\n"
        "- The manuscript's Section 6.1 treatment is CORRECT and untouched: it "
        "attributes the numbers to vascular permeability, names the vessel "
        "wall as the highest-concentration position, and calls the transport "
        "parameters estimated. Only the Chapter 5 sentence inherited from the "
        "uncertainty generator is corrected.\n"
        "- The two drugs differ in transport only; the cell-level pharmacology "
        "is the RSL3/GPX4 pathway in both arms, as `sim-tissue-pk` states.\n"
        "- The deepest-bin comparison is at each tissue's OWN half-distance, "
        "which is the worst-perfused cell it contains, not a common depth. At "
        "a common depth the shared lambda cancels from the comparison too and "
        "the access ordering returns.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description="measure the penetration headline's lambda axis")
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the prose from the committed JSON, no binary")
    a = ap.parse_args()

    if a.render_only:
        # RE-ASSEMBLED, not echoed: the stored file carries derived fields next
        # to the raw measurements, and re-rendering those instead of
        # recomputing them would make every guard on a derived field inert.
        data = assemble(json.loads(OUT_JSON.read_text()))
    else:
        binary = a.binary or default_binary()
        if not binary.exists():
            print(f"error: sim-tissue-pk not built at {binary}; "
                  "(cd simulations && cargo build --release -p sim-tissue-pk)",
                  file=sys.stderr)
            return 2
        raw = run_binary(binary)
        print("  ran sim-tissue-pk for "
              f"{len(raw['penetration_length_um'])} drugs", flush=True)
        raw["dose_response"] = measure_dose_response()
        print(f"  measured f(c) on {len(DOSE_GRID)} points", flush=True)
        raw["prior_ordering"] = measure_prior_ordering()
        print(f"  {N_DRAWS} prior draws: ordering "
              f"{raw['prior_ordering']['ordering_preserved']}/{N_DRAWS}", flush=True)
        data = assemble(raw)
    OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
