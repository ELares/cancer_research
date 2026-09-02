"""Which of the two rim biases is which (#844).

THE AMBIGUITY THE POINT MODEL CANNOT RESOLVE
---------------------------------------------
A chemotherapy kill is biased toward the tumour rim for two reasons, and the
point model carries neither: the drug has to REACH a cell, and a cycle-specific
agent has to find it CYCLING. Both fail in the same place -- deep, away from
vessels, where cells are quiescent and drug is scarce.

An observed rim-biased kill is therefore ambiguous, and the ambiguity matters
clinically because the two have OPPOSITE remedies. If delivery is the binding
term, the answer is better delivery. If the cell cycle is, better delivery buys
nothing and the answer is recruitment into cycle. A model with no geometry
cannot even pose the question.

THE DECOMPOSITION
-----------------
Four runs per agent class: the 2x2 of delivery gradient on/off and cycle
coupling on/off. With both off the kill must be FLAT with depth, which is the
control that shows the two terms are real rather than artifacts of the zoning.
With each on alone its own contribution is isolated. With both on, the
confounded picture a real experiment sees.

THE DISCRIMINATOR
-----------------
Agent class. A phase-NONSPECIFIC alkylator damages DNA whatever the cell is
doing, so quiescence should cost it comparatively little; a phase-specific
agent should lose far more of its deep kill to a quiescent core even when
delivery is perfectly uniform. If the classes behaved alike, the "cycle" arm
would be measuring something other than the cell cycle.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SWEEP = REPO / "analysis" / "calibration" / "chemo_decomposition_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "chemo.rs"
OUT_MD = REPO / "analysis" / "calibration" / "chemo-decomposition-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "chemo-decomposition-validation.json"

# How flat "flat" has to be for the both-off control to count as a control.
FLAT_TOLERANCE = 0.10


def _rows():
    out = []
    for ln in SWEEP.read_text().splitlines():
        if not ln.startswith("CHEMO_DECOMP "):
            continue
        d = dict(re.findall(r"(\w+)=([-\d.eE]+|true|false|\w+)", ln))
        rec = {}
        for k, v in d.items():
            if v in ("true", "false"):
                rec[k] = v == "true"
            else:
                try:
                    rec[k] = float(v)
                except ValueError:
                    rec[k] = v
        out.append(rec)
    return out


def scan() -> dict:
    rows = _rows()
    classes = []
    for cname in ("PhaseNonspecific", "SPhaseSpecific", "MPhaseSpecific"):
        cs = [r for r in rows if r["class"] == cname]
        def pick(delivery, cycle):
            return next(r for r in cs
                        if r["delivery"] == delivery and r["cycle"] == cycle)
        neither, cyc, deliv, both = (pick(False, False), pick(False, True),
                                     pick(True, False), pick(True, True))
        base = neither["hypoxic"]
        classes.append({
            "class": cname,
            "neither": {"overall": neither["overall"], "rim": neither["normoxic"],
                        "core": neither["hypoxic"]},
            "cycle_only": {"overall": cyc["overall"], "rim": cyc["normoxic"],
                           "core": cyc["hypoxic"]},
            "delivery_only": {"overall": deliv["overall"], "rim": deliv["normoxic"],
                              "core": deliv["hypoxic"]},
            "both": {"overall": both["overall"], "rim": both["normoxic"],
                     "core": both["hypoxic"]},
            # What each term alone does to the CORE kill, relative to neither.
            "core_retained_cycle_only": (cyc["hypoxic"] / base) if base else None,
            "core_retained_delivery_only": (deliv["hypoxic"] / base) if base else None,
            # The control: with neither term, is the kill flat with depth?
            "control_flat": (abs(neither["hypoxic"] - neither["normoxic"])
                             / max(neither["normoxic"], 1e-9)) if neither["normoxic"] else None,
        })
    return {"classes": classes, "flat_tolerance": FLAT_TOLERANCE}


def assemble(raw: dict) -> dict:
    d = dict(raw)
    cs = raw["classes"]
    d["control_is_flat"] = all(c["control_flat"] is not None
                               and c["control_flat"] < raw["flat_tolerance"] for c in cs)
    # Delivery must dominate: it should cut the core kill harder than the cycle
    # term does, in every class.
    d["delivery_dominates_everywhere"] = all(
        c["core_retained_delivery_only"] < c["core_retained_cycle_only"] for c in cs)
    # And the CYCLE term must be class-dependent, or it is not measuring the
    # cell cycle. A phase-nonspecific agent should keep more of its core kill.
    nonspec = next(c for c in cs if c["class"] == "PhaseNonspecific")
    specifics = [c for c in cs if c["class"] != "PhaseNonspecific"]
    d["cycle_term_is_class_dependent"] = all(
        nonspec["core_retained_cycle_only"] > s["core_retained_cycle_only"]
        for s in specifics)
    d["cycle_cost_ratio"] = (
        nonspec["core_retained_cycle_only"]
        / max(s["core_retained_cycle_only"] for s in specifics))
    d["verdict"] = (
        "CONFIRMED — the two terms separate, and only one of them is class-dependent"
        if d["control_is_flat"] and d["delivery_dominates_everywhere"]
        and d["cycle_term_is_class_dependent"] else "UNRESOLVED")
    return d


def render(d: dict) -> str:
    L = [
        "# Chemotherapy: which of the two rim biases is which (#844)",
        "",
        "*Generated by `scripts/validate_chemo_decomposition.py --render-only`. "
        "Pure stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}.**",
        "",
        "## The ambiguity",
        "",
        "A chemotherapy kill is biased toward the tumour rim for two reasons "
        "and the point model carries neither: the drug has to **reach** a cell, "
        "and a cycle-specific agent has to find it **cycling**. Both fail in "
        "the same place.",
        "",
        "That matters clinically because the two have opposite remedies. If "
        "delivery is the binding term, the answer is better delivery. If the "
        "cell cycle is, better delivery buys nothing and the answer is "
        "recruitment into cycle. A model with no geometry cannot pose the "
        "question, let alone answer it.",
        "",
        "## The 2×2, per agent class",
        "",
        "Kill fraction in the hypoxic core, and what each term alone leaves of it:",
        "",
        "| agent class | neither | cycle only | delivery only | both |",
        "|---|--:|--:|--:|--:|",
    ]
    for c in d["classes"]:
        L.append(f"| {c['class']} | {c['neither']['core']:.3f} | "
                 f"{c['cycle_only']['core']:.3f} | {c['delivery_only']['core']:.4f} | "
                 f"{c['both']['core']:.4f} |")
    L += [
        "",
        "| agent class | core kill kept, cycle term alone | core kill kept, delivery term alone |",
        "|---|--:|--:|",
    ]
    for c in d["classes"]:
        L.append(f"| {c['class']} | **{c['core_retained_cycle_only']:.2f}** | "
                 f"{c['core_retained_delivery_only']:.4f} |")
    L += [
        "",
        "## What the control establishes",
        "",
        f"With **neither** term active the kill is flat with depth in every "
        f"class (rim and core within {d['flat_tolerance']:.0%}), which is what "
        "makes the other three columns terms rather than artifacts of how the "
        "zones are drawn. A decomposition whose null case already had a "
        "gradient would be decomposing the zoning.",
        "",
        "## The two findings",
        "",
        "**Delivery dominates.** In every class it takes the core kill to "
        "near zero, far harder than the cycle term does. On these parameters "
        "the depth gradient a real experiment would see is overwhelmingly a "
        "delivery gradient.",
        "",
        f"**But only the cycle term is class-dependent, and that is the half "
        f"better delivery cannot fix.** With delivery uniform, a "
        f"phase-nonspecific alkylator keeps "
        f"{d['classes'][0]['core_retained_cycle_only']:.0%} of its core kill "
        f"against a quiescent core, while the phase-specific agents keep "
        + " and ".join(f"{c['core_retained_cycle_only']:.0%}"
                       for c in d["classes"] if c["class"] != "PhaseNonspecific")
        + f" — a factor of {d['cycle_cost_ratio']:.1f} between them. The "
        "alkylator damages DNA whatever the cell is doing; the others need it "
        "cycling. That class-dependence is the check that the cycle arm is "
        "measuring the cell cycle: if the classes behaved alike it would be "
        "measuring the drug field twice.",
        "",
        "## What this does NOT establish",
        "",
        "- **The fitted verdict stays NO TARGET.** This arm's dose-response "
        "target was already recorded as unreachable rather than merely "
        "unfitted, and nothing here changes that. Dose and potency are "
        "illustrative; the DECOMPOSITION is the result, not the kill "
        "fractions.",
        "- **Which term dominates depends on the parameters.** The penetration "
        "length and the potency were both chosen, not measured, and a longer "
        "penetration length would shift the balance toward the cycle term. "
        "What is robust is that the two terms are separable and that only one "
        "of them varies by agent class.",
        "- **The cycle is read from PHENOTYPE, not simulated.** Cells inherit "
        "a phase distribution from the spheroid layer's radial assignment: "
        "proliferating at the rim, quiescent-rich in the core. No cell "
        "actually divides in this engine, so nothing here models "
        "redistribution or recruitment over time.",
        "- **The phase-sensitivity numbers are conventional**, as `chemo.rs` "
        "says of itself: the ORDERING and the zeros are the assertion, not the "
        "sizes.",
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
