#!/usr/bin/env python3
"""When suppressing immunity helps an oncolytic virus, and when it does not.

THE TENSION, AND WHY NAMING IT IS NOT ENOUGH
--------------------------------------------
An oncolytic virus has two mechanisms that want opposite things from the immune
system. Lysis wants the virus to spread, and an immune response clears it. The
durable arm is the anti-TUMOUR immunity that lysis primes, and that needs an
immune system to exist.

Every review of the field states that tension. Stating it settles nothing,
because it does not say which side wins -- and the clinically consequential
question is exactly that: whether suppressing immunity to let the virus work is
a good idea.

WHAT THE MODEL ADDS
-------------------
A CONDITION. The primed response has to repay the lysis it costs, and whether
it does depends on how efficiently lysis primes. Below a crossover this page
computes, full suppression is optimal and the intuition is right. Above it, the
same move throws away the only durable arm.

The crossover is the finding. It is not a number anyone should carry away --
the priming efficiency it is expressed in is a placeholder, like everything
else here -- but its EXISTENCE is a statement about the structure of the
problem rather than about parameters, and a model that produced an interior
optimum at every setting would be asserting the conclusion instead.

THE TRIAL, AND WHAT IT CANNOT ANCHOR
------------------------------------
OPTiM (Andtbacka 2015, PMID 26014293) reported a durable response rate of 16.3%
for talimogene laherparepvec against 2.1% for GM-CSF. That is a real contrast
and it CANNOT anchor this model. The ratio trick that worked for the checkpoint
arm needs two strata of one treatment; these are two different agents, so no
mapping constant is shared and the comparison measures two mechanisms rather
than one mechanism at two settings. It is recorded as DIRECTION-anchored, which
is what it is.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "oncolytic.rs"
OUT_MD = REPO / "analysis" / "calibration" / "oncolytic-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "oncolytic-validation.json"

CFG = {"initial_infected": 0.01, "replication_rate": 0.9, "clearance_rate": 0.2,
       "interferon_competence": 0.3, "lysis_rate": 0.15}
STEPS = 60
CLEARANCE_SENSITIVITY = 0.5
MAX_COMPETENCE = 16.0
INTERIOR_TOLERANCE = 0.05

OPTIM = {"trial": "OPTiM", "pmid": "26014293", "agent": "talimogene laherparepvec",
         "durable_response_pct": 16.3, "comparator": "GM-CSF",
         "comparator_pct": 2.1}


def rust_constants() -> dict:
    src = RUST.read_text()
    m = re.search(r"pub const ENTRY_RECEPTOR_THRESHOLD: f64 = ([0-9.]+);", src)
    if not m:
        raise SystemExit("ERROR: ENTRY_RECEPTOR_THRESHOLD not found")
    return {"ENTRY_RECEPTOR_THRESHOLD": float(m.group(1))}


def simulate_spread(cfg, steps):
    """Re-implemented from the crate's own stepper: infected grows by
    replication scaled by interferon competence, loses to clearance and to
    lysis, and lysed accumulates."""
    infected = cfg["initial_infected"]
    lysed = 0.0
    repl = cfg["replication_rate"] * (1.0 - cfg["interferon_competence"])
    for _ in range(steps):
        new = infected * repl * max(0.0, 1.0 - infected - lysed)
        died = infected * cfg["lysis_rate"]
        cleared = infected * cfg["clearance_rate"]
        infected = max(0.0, min(1.0, infected + new - died - cleared))
        lysed = min(1.0, lysed + died)
    return infected, lysed


def outcome(competence, priming_efficiency, permissive=1.0):
    survival = 1.0 / (1.0 + CLEARANCE_SENSITIVITY * max(competence, 0.0))
    _infected, lysed_frac = simulate_spread(CFG, STEPS)
    lysed = min(1.0, lysed_frac * survival * permissive)
    primed = min(1.0, priming_efficiency * lysed * competence / (1.0 + competence))
    return min(1.0, lysed + primed * (1.0 - lysed))


def optimum(priming_efficiency):
    best = (0.0, float("-inf"))
    n = 200
    for i in range(n + 1):
        c = MAX_COMPETENCE * i / n
        v = outcome(c, priming_efficiency)
        if v > best[1]:
            best = (c, v)
    return best


def scan() -> dict:
    consts = rust_constants()
    efficiencies = [round(0.2 * i, 2) for i in range(1, 51)]
    rows = []
    crossover = None
    for e in efficiencies:
        c, v = optimum(e)
        # SATURATION IS A THIRD STATE, and lumping it in with "the optimum is
        # at full suppression" was wrong. Once the outcome reaches 1 the scan
        # returns the LOWEST competence that gets there, which looks like a
        # boundary optimum and means something entirely different: the
        # trade-off has stopped operating because everything dies either way.
        # Reported separately, the two regimes are ordered; mixed together they
        # interleave and the crossover looks unstable.
        saturated = v >= 1.0 - 1e-9
        interior = (not saturated) and c > MAX_COMPETENCE * INTERIOR_TOLERANCE
        rows.append({"priming_efficiency": e, "optimal_competence": round(c, 3),
                     "outcome": round(v, 4), "interior": interior,
                     "saturated": saturated})
        if interior and crossover is None:
            crossover = e

    curve = [{"competence": round(c, 2),
              "weak_priming": round(outcome(c, 0.5), 4),
              "strong_priming": round(outcome(c, 4.0), 4)}
             for c in [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0]]

    return {
        "constants": consts,
        "config": CFG,
        "steps": STEPS,
        "clearance_sensitivity": CLEARANCE_SENSITIVITY,
        "interior_tolerance_fraction": INTERIOR_TOLERANCE,
        "max_competence": MAX_COMPETENCE,
        "crossover_priming_efficiency": crossover,
        "optimum_by_efficiency": rows,
        "outcome_curves": curve,
        "trial": OPTIM,
    }


def assemble(raw: dict) -> dict:
    rows = [r for r in raw["optimum_by_efficiency"] if not r["saturated"]]
    raw["n_saturated"] = sum(1 for r in raw["optimum_by_efficiency"]
                             if r["saturated"])
    raw["verdict"] = {
        "interior_optimum_is_conditional":
            "YES" if any(r["interior"] for r in rows)
                     and any(not r["interior"] for r in rows) else "NO",
        "crossover": raw["crossover_priming_efficiency"],
        "trial_anchor": "DIRECTION only -- two agents, so no mapping cancels",
    }
    return raw


def render(d: dict) -> str:
    t = d["trial"]
    L = ["# When suppressing immunity helps an oncolytic virus, and when it does not",
         "",
         "*Generated by `scripts/validate_oncolytic.py --render-only`. Pure "
         "stdlib; runs offline in CI. The spread stepper is re-implemented "
         "here rather than imported, so the two implementations can disagree.*",
         "",
         "## The tension, and why naming it settles nothing", "",
         "Lysis wants the virus to spread and an immune response clears it. The "
         "durable arm is the anti-tumour immunity that lysis primes, and that "
         "needs an immune system to exist. Every review of the field states "
         "that. What it does not say is which side wins -- and the "
         "consequential question is exactly that.", "",
         "## The condition", "",
         f"Scanning the efficiency with which lysis primes anti-tumour "
         f"immunity, the optimum moves off full suppression at "
         f"**{d['crossover_priming_efficiency']}**:", "",
         "| priming efficiency | optimal immune competence | outcome | interior? |",
         "|--:|--:|--:|---|"]
    for r in d["optimum_by_efficiency"]:
        if r["priming_efficiency"] in (0.2, 1.0, 2.0, 2.6, 3.0, 4.0, 6.0, 10.0):
            state = ("saturated" if r["saturated"]
                     else "yes" if r["interior"] else "no")
            L.append(f"| {r['priming_efficiency']:.1f} | "
                     f"{r['optimal_competence']:.2f} | {r['outcome']:.3f} | "
                     f"{state} |")
    L += ["",
          "**Below the crossover the model recommends suppression and above it "
          "the same move throws away the durable arm.** That is the finding: "
          "not that there is an optimum, but that whether there is one depends "
          "on a quantity nobody has measured -- so the clinical intuition is "
          "conditionally right, and the condition is nameable.", "",
          "## What the two regimes look like", "",
          "| immune competence | weak priming (0.5) | strong priming (4.0) |",
          "|--:|--:|--:|"]
    for row in d["outcome_curves"]:
        L.append(f"| {row['competence']:.2f} | {row['weak_priming']:.3f} | "
                 f"{row['strong_priming']:.3f} |")
    L += ["",
          "## The trial, and what it cannot anchor", "",
          f"{t['trial']} (PMID {t['pmid']}) reported a durable response rate of "
          f"{t['durable_response_pct']}% for {t['agent']} against "
          f"{t['comparator_pct']}% for {t['comparator']}. That contrast is "
          "real and it cannot anchor this model.", "",
          "The ratio construction that worked for the checkpoint arm needs two "
          "STRATA of one treatment, where the response-to-kill mapping is "
          "shared and cancels. These are two different agents, so nothing "
          "cancels: a ratio between them measures two mechanisms rather than "
          "one mechanism at two settings. The anchor is recorded as "
          "DIRECTION-only.", "",
          "## What this does not establish", "",
          "- **The crossover is not a measurement.** Priming efficiency is a "
          "placeholder in units this model invented; the crossover's EXISTENCE "
          "is the claim, not its position.",
          "- **`interior` is a definition.** An optimum counts as interior when "
          f"it exceeds {d['interior_tolerance_fraction']:.0%} of the scanned "
          "competence range, and the reported crossover moves if that "
          "tolerance does.",
          f"- **The saturating regime is not meaningful, and it is marked.** "
          f"At high priming efficiency the outcome reaches 1 and the scan "
          f"returns the lowest competence that gets there -- which LOOKS like "
          f"a boundary optimum and means the opposite. "
          f"{d['n_saturated']} of {len(d['optimum_by_efficiency'])} rows are "
          "in that regime and are excluded from the crossover, because mixed "
          "in they make the two regimes interleave and the crossover look "
          "unstable.",
          "- **Nothing here is fitted**, including the spread parameters, and "
          "no absolute outcome should be read as a response rate.", ""]
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
    print(f"  interior optimum conditional: "
          f"{d['verdict']['interior_optimum_is_conditional']}, crossover at "
          f"{d['crossover_priming_efficiency']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
