#!/usr/bin/env python3
"""How far each physical modality actually reaches, measured against its physics.

WHY THIS EXISTS
---------------
Chapter 6 and Chapter 10 of the manuscript argue that physical ROS-generating
modalities are worth attention because they deliver energy where a systemic
drug cannot. That argument has always been made with three arms -- PDT, SDT and
RSL3 -- and all three are modalities this project's thesis is about. It had no
CONTROL: a physical modality whose reach is not the limiting factor at all.

`Treatment::Radiation` (#726) is that control, and this page is what it buys.
Megavoltage photons lose about 3% per centimetre of soft tissue, so over a
whole tumour the delivered dose is nearly flat. Putting that beside light dying
in millimetres and ultrasound reaching centimetres turns "physical delivery"
from one story into a range, and lets the manuscript say which end of that
range its thesis sits at.

WHAT IS MEASURED, AND WHAT IS ASSUMED
-------------------------------------
Two independent quantities per modality, and they are NOT the same claim:

* **Delivered energy** -- the modality's own attenuation law, evaluated from
  the engine's constants. This is physics and is well characterised.
* **Observed kill** -- what `sim-spatial` actually produced at each depth.
  This inherits every uncalibrated parameter the engine has, and the
  identifiability report prices those: none of these magnitudes is
  point-estimable, the DIRECTIONS are the result.

The two are reported side by side because a modality can deliver energy and
still not kill (RSL3 delivers uniformly and kills 3-4%: its limit is
biochemical, not penetration), and that dissociation is the actual finding.

RADIATION'S KILL IS NOT THE SAME KIND OF NUMBER as the other three. It comes
from the linear-quadratic DNA-damage model, which is the one arm here with an
external published parameterisation behind it; the other three come from the
ferroptosis engine's uncalibrated ROS chemistry. The table says so per row
rather than inviting a reader to compare them as like for like.

Offline: reads `simulations/output/spatial/depth_kill_curves.csv` (gitignored,
regenerate with `cargo run --release -p sim-spatial`) and the crate constants.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CURVES = REPO / "simulations" / "output" / "spatial" / "depth_kill_curves.csv"
PARAMS_RS = REPO / "simulations" / "ferroptosis-core" / "src" / "params.rs"
RADIATION_RS = REPO / "simulations" / "ferroptosis-core" / "src" / "radiation.rs"

OUT_MD = REPO / "analysis" / "depth-reach-comparison.md"
OUT_JSON = REPO / "analysis" / "depth-reach-comparison.json"

BIN_UM = 250.0
# The spheroid's poles hold a handful of cells, so the deepest bin's rate is
# noise at every treatment. Dropped from the summary rows and SAID, rather than
# quietly included -- an earlier figure annotation pointed at it and reported
# 100% kill for two modalities that disagree everywhere else.
DROP_LAST_BINS = 1


def _rust_float(path: Path, pattern: str) -> float:
    m = re.search(pattern, path.read_text())
    if not m:
        raise SystemExit(f"{pattern} not found in {path}")
    return float(m.group(1))


def _constants() -> dict:
    return {
        "pdt_mu_eff_per_mm": _rust_float(PARAMS_RS, r"\bpdt_mu_eff:\s*([0-9.]+)"),
        "sdt_alpha_db_cm_mhz": _rust_float(PARAMS_RS, r"\bsdt_alpha:\s*([0-9.]+)"),
        "sdt_freq_mhz": _rust_float(PARAMS_RS, r"\bsdt_freq_mhz:\s*([0-9.]+)"),
        "radiation_mu_per_cm": _rust_float(
            RADIATION_RS, r"MU_6MV_SOFT_TISSUE_PER_CM: f64 = ([0-9.]+)"),
    }


def delivered(tx: str, z_mm: float, c: dict) -> float:
    """Relative delivered energy at depth, from each modality's own law."""
    if tx == "PDT":
        return math.exp(-c["pdt_mu_eff_per_mm"] * z_mm)
    if tx == "SDT":
        db = c["sdt_alpha_db_cm_mhz"] * c["sdt_freq_mhz"] * (z_mm / 10.0)
        return 10.0 ** (-db / 10.0)
    if tx == "Radiation":
        return math.exp(-c["radiation_mu_per_cm"] * (z_mm / 10.0))
    if tx == "RSL3":
        return 1.0  # systemic: uniform by construction, not by measurement
    raise KeyError(tx)


def scan() -> dict:
    if not CURVES.exists():
        raise SystemExit(
            f"{CURVES} not found -- run `cargo run --release -p sim-spatial` first")
    rows: dict = defaultdict(list)
    with CURVES.open() as f:
        for r in csv.DictReader(f):
            rows[r["treatment"]].append(
                (float(r["depth_um"]), float(r["death_rate"]), int(r["n_cells"])))
    out = {}
    for tx, pts in rows.items():
        agg: dict = defaultdict(lambda: [0.0, 0])
        for d, rate, n in pts:
            if n <= 0:
                continue
            b = int(d // BIN_UM) * BIN_UM + BIN_UM / 2.0
            agg[b][0] += rate * n
            agg[b][1] += n
        out[tx] = [{"depth_mm": b / 1000.0,
                    "kill_pct": agg[b][0] / agg[b][1] * 100.0,
                    "n_cells": agg[b][1]}
                   for b in sorted(agg)]
    return {"constants": _constants(), "binned": out,
            "bin_um": BIN_UM, "dropped_last_bins": DROP_LAST_BINS}


# Which arm's kill comes from which model. Not decoration: radiation's number
# has a published parameterisation behind it and the other three do not, so a
# reader comparing them as like for like would be wrong.
KILL_SOURCE = {
    "PDT": "ferroptosis engine (uncalibrated ROS chemistry)",
    "SDT": "ferroptosis engine (uncalibrated ROS chemistry)",
    "RSL3": "ferroptosis engine (uncalibrated ROS chemistry)",
    "Radiation": "linear-quadratic DNA damage (published parameterisation)",
    "Control": "untreated baseline",
}


def assemble(raw: dict) -> dict:
    c = raw["constants"]
    rows = []
    for tx, series in sorted(raw["binned"].items()):
        usable = series[:-raw["dropped_last_bins"]] if raw["dropped_last_bins"] else series
        if len(usable) < 2:
            continue
        first, last = usable[0], usable[-1]
        tot = sum(b["n_cells"] for b in usable)
        killed = sum(b["kill_pct"] / 100.0 * b["n_cells"] for b in usable)
        rows.append({
            "treatment": tx,
            "surface_kill_pct": first["kill_pct"],
            "deep_kill_pct": last["kill_pct"],
            "deep_mm": last["depth_mm"],
            "overall_kill_pct": killed / tot * 100.0,
            "delivered_at_deep_pct": (delivered(tx, last["depth_mm"], c) * 100.0
                                      if tx in ("PDT", "SDT", "RSL3", "Radiation")
                                      else None),
            "kill_retained_pct": (last["kill_pct"] / first["kill_pct"] * 100.0
                                  if first["kill_pct"] > 0 else None),
            "kill_source": KILL_SOURCE.get(tx, "unknown"),
        })
    rows.sort(key=lambda r: -(r["delivered_at_deep_pct"] or -1))
    return dict(raw, rows=rows)


def render(d: dict) -> str:
    c = d["constants"]
    rows = [r for r in d["rows"] if r["treatment"] != "Control"]
    deep = max(r["deep_mm"] for r in rows)
    L = ["# How far each modality reaches, and whether reach is what limits it",
         "",
         "*Generated by `scripts/depth_reach_comparison.py --render-only`. "
         "Offline; reads the `sim-spatial` depth curves and the crate "
         "constants.*", "",
         "The manuscript's tissue-access argument has always been made with "
         "three arms, and all three are modalities this project's thesis is "
         "about. It had no CONTROL — a physical modality whose reach is not "
         "the limiting factor. `Treatment::Radiation` (#726) is that control, "
         "and this is what it buys.", "",
         f"| modality | delivered at {deep:.1f} mm | surface kill | "
         f"kill at {deep:.1f} mm | kill retained | overall | kill model |",
         "|---|--:|--:|--:|--:|--:|---|"]
    for r in rows:
        dv = (f"{r['delivered_at_deep_pct']:.1f}%"
              if r["delivered_at_deep_pct"] is not None else "—")
        kr = (f"{r['kill_retained_pct']:.0f}%"
              if r["kill_retained_pct"] is not None else "—")
        L.append(f"| {r['treatment']} | {dv} | {r['surface_kill_pct']:.1f}% | "
                 f"{r['deep_kill_pct']:.1f}% | {kr} | "
                 f"{r['overall_kill_pct']:.1f}% | {r['kill_source']} |")
    L += ["",
          f"Attenuation laws, read from the crate: PDT "
          f"`exp(-{c['pdt_mu_eff_per_mm']}/mm · z)`, SDT "
          f"`10^(-{c['sdt_alpha_db_cm_mhz']}·{c['sdt_freq_mhz']}·z_cm/10)`, "
          f"radiation `exp(-{c['radiation_mu_per_cm']}/cm · z)`, RSL3 uniform "
          "by construction. Rows are binned at "
          f"{d['bin_um']:.0f} µm and the deepest "
          f"{d['dropped_last_bins']} bin(s) dropped — the spheroid's poles "
          "hold a handful of cells, so their rate is noise at every "
          "treatment.", ""]

    rad = next((r for r in rows if r["treatment"] == "Radiation"), None)
    pdt = next((r for r in rows if r["treatment"] == "PDT"), None)
    if rad and pdt:
        L += ["## The finding is a dissociation, not a ranking", "",
              f"Radiation delivers {rad['delivered_at_deep_pct']:.0f}% of its "
              f"surface dose at {rad['deep_mm']:.1f} mm and retains "
              f"{rad['kill_retained_pct']:.0f}% of its surface kill. PDT "
              f"delivers {pdt['delivered_at_deep_pct']:.1f}% and retains "
              f"{pdt['kill_retained_pct']:.0f}%. That is the range the phrase "
              "\"physical modality\" spans, and this project's thesis sits at "
              "the shallow end of it.", "",
              "**But delivered energy is not kill, and the table shows both "
              "because they come apart.** RSL3 delivers uniformly — a systemic "
              "drug has no attenuation — and kills a few percent everywhere: "
              "its limit is biochemical, not penetration. So \"reaches deeper\" "
              "and \"works deeper\" are different claims, and an argument that "
              "moves between them without saying so is not measured.", ""]

    L += ["## What this does not say", "",
          "**The four kill columns are not like for like.** Radiation's comes "
          "from the linear-quadratic model with a published parameterisation "
          "behind it; the other three come from the ferroptosis engine's "
          "uncalibrated ROS chemistry. `analysis/identifiability-report.md` "
          "prices those: with eleven free rate constants and none of these "
          "outputs conditioned on data in the regime that produces them, no "
          "magnitude here is point-estimable. The DIRECTIONS are the result.", "",
          "**Radiation's ferroptosis channel is OFF in these runs.** "
          "`ros_per_gy` is uncalibrated — the literature gives a direction and "
          "no gray-to-ROS conversion — so the radiation column is DNA damage "
          "alone. Switching it on would raise that column by an amount nobody "
          "can currently justify, which is why the default is zero.", "",
          "**Single fraction, and no regrowth.** Fractionation is what makes "
          "radiotherapy work clinically, and it needs cell division between "
          "fractions, which this engine does not model (#727). A single 2 Gy "
          "dose is not a course of radiotherapy and this page is not a claim "
          "about one.", "",
          "**The 2D spatial model has no oxygen field**, so every arm here runs "
          "at full O₂ supply. Radiation's oxygen effect — the one quantity in "
          "this modality with eighty years of measurement behind it — is "
          "therefore invisible in this table. It lives in the 3D model, and "
          "measuring it is the next step rather than something this page "
          "reports.", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for r in d["rows"]:
        if r["treatment"] == "Control":
            continue
        print(f"  {r['treatment']:10s} delivered {r['delivered_at_deep_pct']:6.1f}% "
              f"kill {r['surface_kill_pct']:5.1f}% -> {r['deep_kill_pct']:5.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
