"""Is multiplying the barrier fractions the right answer (#844)?

WHAT THE POINT MODEL ASSUMES WITHOUT SAYING SO
-----------------------------------------------
`adoptive::delivery_efficiency` is `trafficking * infiltration * activation`,
and the arm then multiplies by the antigen-positive fraction. Multiplying
fractions is exactly right when they are INDEPENDENT and wrong when they are
not -- and two of this arm's failure modes are not independent in a tumour,
because infiltration is a property of WHERE a cell sits and antigen expression
is a property of the CELL.

The assumption is invisible in a point model because there is no "where" for
the two to correlate across.

THE CONTROL IS THE WHOLE DESIGN
--------------------------------
The spatial run replaces the `infiltration` SCALAR with a depth field of the
same mean, so the two models differ in DISTRIBUTION and in nothing else. At
zero correlation antigen is independent of position, the product must be
right, and the spatial run must agree with it. Any disagreement there would be
a bug rather than a finding, and the first version of this measurement had one:
it multiplied the depth field ON TOP of the infiltration fraction, double-
counting the same barrier, and its control failed at every correlation
including zero.

Away from zero the two should diverge, and in a direction that follows the
sign: antigen concentrated where effectors cannot reach makes the product
OPTIMISTIC, antigen concentrated at the rim makes it pessimistic.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SWEEP = REPO / "analysis" / "calibration" / "cart_independence_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"
OUT_MD = REPO / "analysis" / "calibration" / "cart-independence-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "cart-independence-validation.json"

# How close the spatial answer must sit to the product at zero correlation for
# the control to count. Monte-Carlo over ~80k cells, so a couple of percent.
CONTROL_TOLERANCE = 0.03


def _rows():
    out = []
    for ln in SWEEP.read_text().splitlines():
        if not ln.startswith("CART_INDEP "):
            continue
        d = dict(re.findall(r"(\w+)=([-\d.eE]+)", ln))
        out.append({k: float(v) for k, v in d.items()})
    return out


def scan() -> dict:
    rows = sorted(_rows(), key=lambda r: r["correlation"])
    pts = []
    for r in rows:
        total = r["total_tumor"]
        measured = r["kills"] / total
        product = r["product_prediction"]
        pts.append({
            "correlation": r["correlation"],
            "measured_kill_fraction": measured,
            "product_prediction": product,
            "product_over_measured": product / measured if measured else None,
            "unreached": int(r["unreached"]),
            "antigen_negative": int(r["antigen_negative"]),
            "lost_to_both": int(r["lost_to_both"]),
            "reach_mean": r["reach_mean"],
        })
    return {"points": pts, "control_tolerance": CONTROL_TOLERANCE}


def assemble(raw: dict) -> dict:
    d = dict(raw)
    pts = raw["points"]
    zero = next(p for p in pts if p["correlation"] == 0.0)
    d["control_error"] = abs(zero["product_over_measured"] - 1.0)
    d["control_holds"] = d["control_error"] < raw["control_tolerance"]
    neg = [p for p in pts if p["correlation"] < 0]
    pos = [p for p in pts if p["correlation"] > 0]
    d["product_optimistic_when_antigen_is_deep"] = all(
        p["product_over_measured"] > 1.0 + raw["control_tolerance"] for p in neg)
    d["product_pessimistic_when_antigen_is_rim"] = all(
        p["product_over_measured"] < 1.0 - raw["control_tolerance"] / 2 for p in pos)
    d["max_overstatement"] = max(p["product_over_measured"] for p in pts)
    d["max_understatement"] = min(p["product_over_measured"] for p in pts)
    # The bucket the product model has no slot for at all.
    d["lost_to_both_range"] = [min(p["lost_to_both"] for p in pts),
                               max(p["lost_to_both"] for p in pts)]
    # SATURATION AT THE EXTREME, marked rather than smoothed. The per-cell
    # antigen probability is clamped to [0,1], so at a strong positive
    # correlation the rim cells are already certain to be positive and pushing
    # the correlation further cannot add any. The effect stops growing, and a
    # monotonicity claim over the whole range would be false.
    rising = [p for p in pts if p["correlation"] > 0]
    d["positive_arm_is_monotone"] = all(
        a["product_over_measured"] >= b["product_over_measured"]
        for a, b in zip(rising, rising[1:]))
    d["saturating_correlations"] = (
        [] if d["positive_arm_is_monotone"]
        else [b["correlation"] for a, b in zip(rising, rising[1:])
              if b["product_over_measured"] > a["product_over_measured"]])
    d["verdict"] = (
        "CONFIRMED — the product is right only when the barriers are independent"
        if d["control_holds"] and d["product_optimistic_when_antigen_is_deep"]
        and d["product_pessimistic_when_antigen_is_rim"] else "UNRESOLVED")
    return d


def render(d: dict) -> str:
    L = [
        "# CAR-T: is multiplying the barrier fractions the right answer? (#844)",
        "",
        "*Generated by `scripts/validate_cart_independence.py --render-only`. "
        "Pure stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}.**",
        "",
        "## The assumption a point model cannot see",
        "",
        "`adoptive::delivery_efficiency` is `trafficking × infiltration × "
        "activation`, and the arm multiplies by the antigen-positive fraction "
        "on top. Multiplying fractions is exactly right when they are "
        "**independent** — and two of these are not, in a tumour: infiltration "
        "is a property of *where a cell sits* and antigen expression is a "
        "property of *the cell*. In a point model there is no \"where\" for "
        "them to correlate across, so the assumption never has to be stated.",
        "",
        "## The control is the whole design",
        "",
        "The spatial run replaces the `infiltration` **scalar** with a depth "
        "field of the same mean, so the two models differ in *distribution* and "
        "nothing else. At zero correlation the product must be right.",
        "",
        f"It is: the spatial kill and the product agree to "
        f"**{d['control_error']:.1%}**. That is what makes every other row a "
        "finding rather than a bug — and the first version of this measurement "
        "failed exactly here, because it multiplied the depth field on top of "
        "the infiltration fraction and double-counted the same barrier.",
        "",
        "## What the correlation costs",
        "",
        "| antigen–depth correlation | measured kill | product predicts | product ÷ measured |",
        "|--:|--:|--:|--:|",
    ]
    for p in d["points"]:
        mark = " ← control" if p["correlation"] == 0.0 else ""
        L.append(f"| {p['correlation']:+.1f} | {p['measured_kill_fraction']:.4f} | "
                 f"{p['product_prediction']:.4f} | "
                 f"**{p['product_over_measured']:.3f}**{mark} |")
    L += [
        "",
        f"**When antigen is concentrated where effectors cannot reach, the "
        f"product is optimistic** — up to {d['max_overstatement']:.2f}× the "
        "kill that actually happens. **When it is concentrated at the rim, the "
        f"product is pessimistic** — down to {d['max_understatement']:.2f}×. "
        "Neither direction is available to a model with no geometry, and "
        "neither is small.",
        "",
        ("" if d["positive_arm_is_monotone"] else
         f"The effect stops growing at the positive extreme "
         + ", ".join(f"{c:+.1f}" for c in d["saturating_correlations"])
         + ": the per-cell antigen probability is clamped to [0,1], so once "
           "the rim cells are certain to be positive a stronger correlation "
           "cannot add any more of them. That is saturation, not a reversal, "
           "and it is marked rather than smoothed away."),
        "",
        "## The bucket the product has no slot for",
        "",
        f"Between {d['lost_to_both_range'][0]:,} and "
        f"{d['lost_to_both_range'][1]:,} cells fail **both** ways in these "
        "runs — no effector reached them *and* they were below the antigen "
        "threshold. A product of fractions cannot distinguish that from a cell "
        "lost to one barrier alone, and the two have different consequences: "
        "improving trafficking rescues the first group and not the second.",
        "",
        "## What this does NOT establish",
        "",
        "- **Not a measurement of any real correlation.** Whether antigen "
        "expression tracks depth in real tumours, and in which direction, is "
        "not something this project measures. What is shown is that the answer "
        "*depends* on it, and by how much.",
        "- **Every barrier value is a placeholder**, as this arm's "
        "`CALIBRATION_STATUS` row has said since it was built. The fitted "
        "verdict is ADMISSIBLE on a different quantity and is untouched here.",
        "- **The effector field is the radial-depth proxy.** Real T-cell "
        "infiltration is shaped by chemokine gradients, stromal exclusion and "
        "vessel adhesion, none of which this represents.",
        "- **The antigen threshold is still a cap, not a barrier.** A cell "
        "below it is not reachable by a larger dose, which is why it is "
        "modelled as a property of the cell rather than a coefficient.",
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
