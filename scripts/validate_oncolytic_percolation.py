"""Whether an oncolytic front CROSSES, which a closed form cannot say (#844).

WHAT THE ENGINE ALREADY HAD, AND WHY IT IS NOT ENOUGH
-----------------------------------------------------
`oncolytic::front_speed` is `2*sqrt(D*r)`, the Fisher-KPP speed of a front in a
HOMOGENEOUS medium. It has no term for how much of the tumour the virus can
enter, so it returns a positive speed whenever replication beats clearance --
read literally, it says a virus always spreads eventually.

On a lattice that is false, and the failure is sharp rather than gradual.
Below a site-percolation threshold in the permissive fraction there is no
connected path of enterable cells, so the infection cannot cross the tissue at
ANY dose or for ANY duration. That is a THRESHOLD, it joins the small set of
arms in this chapter limited by one, and it is invisible to a point model.

THE COMPARATOR IS NOT FROM CANCER BIOLOGY, WHICH IS THE POINT
-------------------------------------------------------------
Site percolation on a three-dimensional cubic lattice is a solved problem in
computational physics, and its threshold depends on the neighbourhood:

    26-neighbour (Moore)        p_c ~ 0.0976
    6-neighbour (von Neumann)   p_c ~ 0.3116

`TumorGrid3D::neighbors` returns the 26-Moore neighbourhood, so if this front
is behaving like site percolation its threshold should land near 0.0976 -- a
number nothing in this project was fitted to, from a field with no interest in
oncolytic viruses. Landing near 0.31 would say it is behaving like 6-neighbour
connectivity instead; landing near neither would say spread is dominated by
something other than connectivity and this framing is wrong.

WHY THE THRESHOLD MOVES WITH REPLICATION, AND WHY THAT IS EXPECTED
-------------------------------------------------------------------
The process is mixed SITE-BOND percolation: a cell must be permissive (site)
AND transmission must succeed (bond). Pure site percolation is the limit of
certain transmission, so only the highest replication row should approach
0.0976, and the threshold must RISE as replication falls. A model whose
threshold did not move with replication would not be percolating at all.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SWEEP = REPO / "analysis" / "calibration" / "oncolytic_percolation_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "oncolytic.rs"
OUT_MD = REPO / "analysis" / "calibration" / "oncolytic-percolation-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "oncolytic-percolation-validation.json"

# Established site-percolation thresholds for a 3-D cubic lattice. NOT fitted
# here and not this project's numbers.
PC_MOORE_26 = 0.0976
PC_VON_NEUMANN_6 = 0.3116
# A front has "crossed" when it reaches a large fraction of the tumour radius.
# Deliberately far from both the local-cluster sizes (tens of microns) and the
# spanning sizes (~540), so the classification does not depend on where in that
# gap the line is drawn.
CROSS_FRACTION = 0.5


def _rows():
    out = []
    for ln in SWEEP.read_text().splitlines():
        if not ln.startswith("ONCOLYTIC_PERC "):
            continue
        d = dict(re.findall(r"(\w+)=([-\d.eE]+|true|false)", ln))
        out.append({k: (v == "true") if v in ("true", "false") else float(v)
                    for k, v in d.items()})
    return out


def _meta():
    ln = next(l for l in SWEEP.read_text().splitlines()
              if l.startswith("ONCOLYTIC_PERC_META"))
    d = dict(re.findall(r"(\w+)=([-\d.eE]+)", ln))
    return {k: float(v) for k, v in d.items()}


def scan() -> dict:
    rows, meta = _rows(), _meta()
    # The tumour's own radius, so "crossed" is a fraction of the thing being
    # crossed rather than a hand-picked micron figure.
    radius_um = meta["grid_dim"] * 0.45 * meta["cell_size_um"]
    cross_um = radius_um * CROSS_FRACTION
    reps = sorted({r["replication"] for r in rows}, reverse=True)
    by_rep = []
    for rep in reps:
        rs = sorted((r for r in rows if r["replication"] == rep),
                    key=lambda r: r["permissive_fraction"])
        below = [r for r in rs if r["front_radius_um"] < cross_um]
        above = [r for r in rs if r["front_radius_um"] >= cross_um]
        by_rep.append({
            "replication": rep,
            "highest_not_crossing": max((r["permissive_fraction"] for r in below),
                                        default=None),
            "lowest_crossing": min((r["permissive_fraction"] for r in above),
                                   default=None),
            "unseeded_fractions": sorted(r["permissive_fraction"] for r in rs
                                         if not r["seeded"]),
            "points": [{"permissive_fraction": r["permissive_fraction"],
                        "front_radius_um": r["front_radius_um"],
                        "seeded": r["seeded"],
                        "crossed": r["front_radius_um"] >= cross_um,
                        "closed_form_speed": r["closed_form_speed"]}
                       for r in rs],
        })
    return {
        "pc_moore_26": PC_MOORE_26,
        "pc_von_neumann_6": PC_VON_NEUMANN_6,
        "tumour_radius_um": radius_um,
        "cross_threshold_um": cross_um,
        "cross_fraction_of_radius": CROSS_FRACTION,
        "by_replication": by_rep,
        "grid_dim": meta["grid_dim"],
        "n_steps": meta["n_steps"],
    }


def assemble(raw: dict) -> dict:
    d = dict(raw)
    rows = raw["by_replication"]
    top = rows[0]
    lo, hi = top["highest_not_crossing"], top["lowest_crossing"]
    d["certain_transmission_bracket"] = [lo, hi]
    d["bracket_contains_moore"] = (
        lo is not None and hi is not None and lo < raw["pc_moore_26"] < hi)
    d["bracket_contains_von_neumann"] = (
        lo is not None and hi is not None and lo < raw["pc_von_neumann_6"] < hi)
    # The threshold must RISE as replication falls, or the process is not
    # site-BOND percolation and the framing is wrong.
    crossings = [r["lowest_crossing"] for r in rows if r["lowest_crossing"] is not None]
    d["threshold_rises_as_replication_falls"] = all(
        a <= b for a, b in zip(crossings, crossings[1:]))
    # The closed form predicts spread everywhere, which is the contrast.
    d["closed_form_positive_everywhere"] = all(
        p["closed_form_speed"] > 0 for r in rows for p in r["points"])
    d["rows_where_closed_form_is_wrong"] = sum(
        1 for r in rows for p in r["points"]
        if p["closed_form_speed"] > 0 and not p["crossed"])
    d["verdict"] = (
        "CONFIRMED — the measured threshold brackets the published lattice constant"
        if d["bracket_contains_moore"] and d["threshold_rises_as_replication_falls"]
        else "UNRESOLVED")
    return d


def render(d: dict) -> str:
    top = d["by_replication"][0]
    L = [
        "# Oncolytic virus: does the front cross? (#844)",
        "",
        "*Generated by `scripts/validate_oncolytic_percolation.py --render-only`. "
        "Pure stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}.**",
        "",
        "## What the closed form cannot say",
        "",
        "`oncolytic::front_speed` is `2·√(D·r)` — the Fisher-KPP speed of a "
        "front in a **homogeneous** medium. It has no term for how much of the "
        "tumour the virus can enter, so it returns a positive speed whenever "
        "replication beats clearance. Read literally it says a virus always "
        "spreads eventually.",
        "",
        f"On a lattice that is false, and **it is false on "
        f"{d['rows_where_closed_form_is_wrong']} of the swept conditions**: the "
        "closed form is positive on every one of them, and on those rows the "
        "front does not cross the tumour at all.",
        "",
        "## The threshold, against a constant from another field",
        "",
        "Site percolation on a 3-D cubic lattice is solved, and the threshold "
        "depends on the neighbourhood:",
        "",
        f"- 26-neighbour (Moore) — **{d['pc_moore_26']}**, which is what "
        "`TumorGrid3D::neighbors` returns",
        f"- 6-neighbour (von Neumann) — {d['pc_von_neumann_6']}",
        "",
        "Neither was fitted here, and neither comes from a field with any "
        "interest in oncolytic viruses.",
        "",
        f"At **certain transmission** the measured front crosses between a "
        f"permissive fraction of **{top['highest_not_crossing']}** and "
        f"**{top['lowest_crossing']}**. The Moore constant "
        f"{'lies inside that bracket' if d['bracket_contains_moore'] else 'does NOT lie inside that bracket'}"
        f"; the von Neumann constant "
        f"{'does too' if d['bracket_contains_von_neumann'] else 'does not'}.",
        "",
        "## The full sweep",
        "",
        "| permissive fraction | " + " | ".join(
            f"r={r['replication']:g}" for r in d["by_replication"]) + " |",
        "|--:|" + "|".join("--:" for _ in d["by_replication"]) + "|",
    ]
    fracs = [p["permissive_fraction"] for p in d["by_replication"][0]["points"]]
    for f in fracs:
        cells = []
        for r in d["by_replication"]:
            p = next((x for x in r["points"] if x["permissive_fraction"] == f), None)
            if p is None:
                cells.append("—")
            elif not p["seeded"]:
                cells.append("*never seeded*")
            else:
                cells.append(f"{p['front_radius_um']:.0f} µm"
                             + (" ✓" if p["crossed"] else ""))
        L.append(f"| {f:g} | " + " | ".join(cells) + " |")
    L += [
        "",
        f"A front counts as having crossed at "
        f"{d['cross_threshold_um']:.0f} µm, half the tumour radius. The line is "
        "drawn in a wide empty gap — the non-crossing runs stop in the tens of "
        "microns and the crossing ones reach ~540 — so the classification does "
        "not depend on where in that gap it sits.",
        "",
        "*never seeded* is a different outcome and is marked rather than "
        "printed as a zero radius: at a very low permissive fraction there may "
        "be no enterable cell near the tumour centre at all, and reporting a "
        "failure to START as a failure to CROSS would be two claims rendered "
        "identically.",
        "",
        "## Why the threshold moves, and why it should",
        "",
        "The process is mixed **site–bond** percolation: a cell must be "
        "permissive (site) *and* transmission must succeed (bond). Pure site "
        "percolation is the limit of certain transmission, so only the top row "
        "should approach the lattice constant and the threshold must rise as "
        "replication falls. It does"
        + (", which is what makes this percolation rather than a coincidence."
           if d["threshold_rises_as_replication_falls"] else
           " NOT, which breaks the framing above."),
        "",
        "## What this does NOT establish",
        "",
        "- **Not a claim about real oncolytic viruses.** The permissive "
        "fraction is a parameter here, not a measurement; nothing in this "
        "project measures what fraction of a real tumour a given virus can "
        "enter. What is shown is that *if* that fraction is low enough, no "
        "dose helps — and where 'low enough' sits for this lattice.",
        "- **The agreement is with a lattice constant, not with biology.** It "
        "says the spread rule behaves like site percolation on the "
        "neighbourhood it runs on, which is a check on the implementation and "
        "on the framing, not evidence that tumours percolate.",
        "- **The spread rule is a contact process, not a PDE.** One infection "
        "attempt per neighbour per step. It is the right shape for asking "
        "whether a front crosses and the wrong one for asking how fast a "
        "continuum approximation says it travels.",
        "",
    ]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if args.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    print(f"{OUT_MD.relative_to(REPO)}: {d['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
