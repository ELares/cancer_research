#!/usr/bin/env python3
"""Validate the engine's OER form against published radiobiology (#726).

`analysis/oxygen-form-check.md` established that the engine's linear O2
dependence is not the oxygen enhancement ratio and that the two disagree most
below 10 mmHg -- the regime carrying manuscript section 7.1 and prediction P4,
the leg this project calls its weakest. `ferroptosis-core` now implements the
measured form; this checks the implementation against the published anchors and
against the Rust source, so neither can drift from the other silently.

WHAT IS BEING VALIDATED, precisely. The FUNCTIONAL FORM, not a magnitude. The
Alper-Howard-Flanders relation `m(p) = (3p + 3)/(p + 3)` is a shape with two
published anchors that fix it completely: `m(0) = 1` (anoxia) and an asymptote
of 3, with half the total rise reached at `p = 3` mmHg. Those are the checks
below. The Type I/Type II fraction the form scales is a separate quantity and
stays uncalibrated; nothing here claims otherwise.

WHY A DRIFT GUARD RATHER THAN A FIT. There is nothing to fit: the form is not
free, it is taken from the literature whole. What CAN go wrong is the Rust and
this file disagreeing after an edit to either, which is the failure mode
`validate_trigger_wave.py` guards for the trigger-wave constants. So the Rust
source is parsed and compared, and a mismatch fails.

OFFLINE: pure stdlib, no compiled extension, no network. Reads the Rust source
as text.
"""
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OXY = REPO / "simulations/ferroptosis-core/src/oxygen.rs"
OUT_MD = REPO / "analysis/calibration/oer-form-validation.md"
OUT_JSON = REPO / "analysis/calibration/oer-form-validation.json"

# Published anchors for the Alper-Howard-Flanders OER. These are properties of
# the relation, not fitted values: m(0)=1 by construction, the asymptote is the
# maximum enhancement, and K is the pO2 at which half the rise is reached.
ANCHOR_ANOXIC = 1.0
ANCHOR_ASYMPTOTE = 3.0
ANCHOR_HALF_MAX_MMHG = 3.0
TOL = 1e-9


def oer(p_mmhg: float) -> float:
    """m(p) = (3p + 3)/(p + 3)."""
    p = max(0.0, p_mmhg)
    return (3.0 * p + 3.0) / (p + 3.0)


def linear(o2_supply: float, dependence: float = 1.0) -> float:
    """The form the engine used before this: (1 - d) + d*s."""
    d = min(1.0, max(0.0, dependence))
    s = min(1.0, max(0.0, o2_supply))
    return min(1.0, max(0.0, 1.0 - d + d * s))


def oer_factor(o2_supply: float, dependence: float, p_full: float) -> float:
    """The engine's `oer_exo_factor`, reimplemented from the same definition."""
    d = min(1.0, max(0.0, dependence))
    s = min(1.0, max(0.0, o2_supply))
    rel = min(1.0, max(0.0, oer(s * p_full) / oer(p_full)))
    return min(1.0, max(0.0, 1.0 - d + d * rel))


def rust_constants() -> dict:
    """Read the constants and the formula out of the Rust source.

    A validation that reimplements the formula and never looks at the code is
    checking itself. This parses what actually ships.
    """
    src = OXY.read_text(encoding="utf-8")
    m = re.search(r"pub const OER_REFERENCE_PO2_MMHG: f64 = ([\d.]+);", src)
    if not m:
        raise SystemExit("OER_REFERENCE_PO2_MMHG not found in oxygen.rs")
    body = re.search(
        r"pub fn oxygen_enhancement_ratio\(p_mmhg: f64\) -> f64 \{(.*?)\n\}",
        src, re.S)
    if not body:
        raise SystemExit("oxygen_enhancement_ratio not found in oxygen.rs")
    formula = " ".join(body.group(1).split())
    return {"reference_po2_mmhg": float(m.group(1)), "formula": formula,
            "has_oer_exo_factor": "pub fn oer_exo_factor(" in src,
            "has_relative_efficacy": "pub fn oer_relative_efficacy(" in src}


def build() -> dict:
    rc = rust_constants()
    p_full = rc["reference_po2_mmhg"]
    checks = [
        ("anoxic value m(0)", oer(0.0), ANCHOR_ANOXIC),
        ("half-maximal at K mmHg", oer(ANCHOR_HALF_MAX_MMHG),
         (ANCHOR_ANOXIC + ANCHOR_ASYMPTOTE) / 2.0),
        ("asymptote", oer(1.0e9), ANCHOR_ASYMPTOTE),
    ]
    rows = [{"pO2_mmhg": p, "o2_supply": round(p / p_full, 4),
             "linear": round(linear(p / p_full), 4),
             "oer": round(oer_factor(p / p_full, 1.0, p_full), 4)}
            for p in (0, 0.5, 1, 2, 3, 5, 7.5, 10, 20, p_full)]
    for r in rows:
        r["gap"] = round(r["oer"] - r["linear"], 4)
    worst = max(rows, key=lambda r: r["gap"])
    # The formula must be the published one, read from the source.
    expected = "let p = p_mmhg.max(0.0); (3.0 * p + 3.0) / (p + 3.0)"
    return {
        "rust": rc,
        "anchors": [{"name": n, "computed": round(v, 12),
                     "published": e, "ok": abs(v - e) < 1e-6}
                    for n, v, e in checks],
        "formula_matches_published": rc["formula"] == expected,
        "expected_formula": expected,
        "rows": rows,
        "worst_gap_pO2": worst["pO2_mmhg"],
        "worst_gap": worst["gap"],
        "anoxic_linear": rows[0]["linear"],
        "anoxic_oer": rows[0]["oer"],
        "endpoints_agree": abs(rows[-1]["oer"] - rows[-1]["linear"]) < TOL,
    }


def render(d: dict) -> str:
    L = ["# The engine's oxygen-effect form against published radiobiology\n"]
    L.append(
        f"Generated by `scripts/validate_oer_form.py`. Validates the FORM the "
        f"engine uses for the O2-dependent arm of the exogenous-ROS yield, not "
        f"the fraction it scales. Reference pO2 read from the Rust source: "
        f"**{d['rust']['reference_po2_mmhg']:.0f} mmHg**.\n")
    L.append("## Published anchors\n")
    L.append("| anchor | computed | published | |")
    L.append("|---|--:|--:|---|")
    for a in d["anchors"]:
        L.append(f"| {a['name']} | {a['computed']:.6f} | {a['published']:.6f} "
                 f"| {'OK' if a['ok'] else '**MISMATCH**'} |")
    L.append("")
    L.append(
        "These fix the relation completely -- it has no free parameter, which "
        "is why this is a drift guard rather than a fit. There is nothing to "
        "tune: the shape is taken from the literature whole.\n")
    L.append("## What the form change does\n")
    L.append("| pO2 (mmHg) | o2_supply | linear | OER | gap |")
    L.append("|--:|--:|--:|--:|--:|")
    for r in d["rows"]:
        L.append(f"| {r['pO2_mmhg']} | {r['o2_supply']} | {r['linear']} | "
                 f"{r['oer']} | {r['gap']} |")
    L.append("")
    L.append(
        f"The forms agree exactly at full supply (normalisation is at "
        f"`p_full`, so an A/B isolates SHAPE rather than scale) and part "
        f"company hardest at **{d['worst_gap_pO2']} mmHg** "
        f"(gap {d['worst_gap']}). At anoxia the linear form gives "
        f"{d['anoxic_linear']} and the measured hyperbola "
        f"{d['anoxic_oer']}: the straight line says a fully anoxic cell "
        f"generates NO exogenous ROS, and the measurement says it retains "
        f"about a third.\n")
    L.append("## Scope\n")
    L.append(
        "This validates a functional form against published anchors. It does "
        "NOT calibrate the Type I/Type II fraction the form scales, which is a "
        "sonosensitizer-specific quantity and stays uncalibrated -- see the "
        "row above it in `CALIBRATION_STATUS.md`. It also says nothing about "
        "whether the exogenous-ROS mechanism is right; it says that IF the "
        "yield depends on oxygen, this is the shape eighty years of "
        "measurement give for that dependence.\n")
    return "\n".join(L)


def main() -> int:
    d = build()
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    bad = [a["name"] for a in d["anchors"] if not a["ok"]]
    print(f"wrote {OUT_MD}")
    print(f"  anchors: {len(d['anchors']) - len(bad)}/{len(d['anchors'])} ok; "
          f"formula matches published: {d['formula_matches_published']}; "
          f"worst gap {d['worst_gap']} at {d['worst_gap_pO2']} mmHg")
    return 1 if (bad or not d["formula_matches_published"]) else 0


if __name__ == "__main__":
    raise SystemExit(main())
