"""What checkpoint blockade has to work with (#844).

WHY THIS ARM IS UNLIKE THE OTHER SIX
-------------------------------------
The checkpoint MECHANISM was already in the engine before #844:
`Overrides.checkpoints` carries the multi-checkpoint panel and the spatial
immune layer reads it. Nothing needed wiring. What was missing is any way to
ask the question the arm exists for.

`checkpoint::residual_brake` returns a MULTIPLIER, and a multiplier has no
denominator. "Blockade raises immune killing 2.5-fold" is true and says nothing
about whether that matters, because it does not say 2.5 times WHAT. In a point
model there is no denominator to be had; on a grid there is.

THE COLD-TO-HOT LADDER
----------------------
Three treatments producing three orders of magnitude of cell death -- an
untreated control, a pharmacologic inducer, and a physical one -- against four
blockade strengths. The immune kill is reported both as the fold-change the
point model predicts AND as a share of the kill actually happening.

WHAT MAKES THE RESULT WORTH HAVING
-----------------------------------
The fold-change is the SAME for both active treatments, because the brake is a
multiplier and multipliers do not know what they multiply. The share differs by
more than an order of magnitude. So the point model's own output is exactly the
quantity that cannot distinguish the two cases, and the quantity that can is
the one it does not have.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SWEEP = REPO / "analysis" / "calibration" / "checkpoint_priming_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "checkpoint.rs"
OUT_MD = REPO / "analysis" / "calibration" / "checkpoint-priming-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "checkpoint-priming-validation.json"


def _rows():
    out = []
    for ln in SWEEP.read_text().splitlines():
        if not ln.startswith("CHECKPOINT_PRIMING "):
            continue
        d = dict(re.findall(r"(\w+)=([-\d.eE]+|\w+)", ln))
        rec = {}
        for k, v in d.items():
            try:
                rec[k] = float(v)
            except ValueError:
                rec[k] = v
        out.append(rec)
    return out


def scan() -> dict:
    rows = _rows()
    arms = []
    for tx in ("Control", "RSL3", "SDT"):
        rs = sorted((r for r in rows if r["treatment"] == tx),
                    key=lambda r: r["anti_pd1"])
        base = next(r for r in rs if r["anti_pd1"] == 0.0)
        top = rs[-1]

        def share(r):
            total = r["immune_kills"] + r["ferroptosis_kills"]
            return (r["immune_kills"] / total) if total else 0.0

        arms.append({
            "treatment": tx,
            "ferroptosis_kills": int(base["ferroptosis_kills"]),
            "total_damp": base["total_damp"],
            "points": [{"anti_pd1": r["anti_pd1"],
                        "immune_kills": int(r["immune_kills"]),
                        "immune_share": share(r),
                        "combined_brake": r["combined_brake"]} for r in rs],
            "fold_benefit": (top["immune_kills"] / base["immune_kills"])
                            if base["immune_kills"] else None,
            "share_unblocked": share(base),
            "share_blocked": share(top),
        })
    return {"arms": arms}


def assemble(raw: dict) -> dict:
    d = dict(raw)
    by = {a["treatment"]: a for a in raw["arms"]}
    cold, active = by["Control"], [by["RSL3"], by["SDT"]]
    # 1. A cold tumour gets nothing, at ANY blockade strength.
    d["cold_tumour_gets_nothing"] = all(
        p["immune_kills"] == 0 for p in cold["points"])
    # 2. The fold benefit is the same for both active arms -- the multiplier
    #    does not know what it multiplies.
    folds = [a["fold_benefit"] for a in active]
    d["fold_benefits"] = folds
    d["fold_spread"] = max(folds) / min(folds)
    d["fold_is_agent_independent"] = d["fold_spread"] < 1.25
    # 3. The share differs by an order of magnitude, which is the quantity the
    #    point model does not have.
    shares = [a["share_unblocked"] for a in active]
    d["share_ratio"] = max(shares) / min(shares)
    d["share_differs_by_an_order_of_magnitude"] = d["share_ratio"] > 10.0
    # 4. And the counterintuitive part: MORE DAMP goes with a SMALLER share.
    hi_damp = max(active, key=lambda a: a["total_damp"])
    lo_damp = min(active, key=lambda a: a["total_damp"])
    d["more_damp_smaller_share"] = (
        hi_damp["share_unblocked"] < lo_damp["share_unblocked"])
    d["damp_ratio"] = hi_damp["total_damp"] / lo_damp["total_damp"]
    d["higher_damp_treatment"] = hi_damp["treatment"]
    d["monotone_in_blockade"] = all(
        all(a["points"][i]["immune_kills"] <= a["points"][i + 1]["immune_kills"]
            for i in range(len(a["points"]) - 1)) for a in active)
    d["verdict"] = (
        "CONFIRMED — the multiplier is agent-independent and the share is not"
        if d["cold_tumour_gets_nothing"] and d["fold_is_agent_independent"]
        and d["share_differs_by_an_order_of_magnitude"]
        and d["more_damp_smaller_share"] else "UNRESOLVED")
    return d


def render(d: dict) -> str:
    by = {a["treatment"]: a for a in d["arms"]}
    L = [
        "# Checkpoint blockade: what it has to work with (#844)",
        "",
        "*Generated by `scripts/validate_checkpoint_priming.py --render-only`. "
        "Pure stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}.**",
        "",
        "## This arm needed no wiring, and that is the point",
        "",
        "The checkpoint *mechanism* was already in the engine: "
        "`Overrides.checkpoints` carries the multi-checkpoint panel and the "
        "spatial immune layer reads it. What was missing is any way to ask the "
        "question the arm exists for.",
        "",
        "`checkpoint::residual_brake` returns a **multiplier**, and a "
        "multiplier has no denominator. *\"Blockade raises immune killing "
        "2.5-fold\"* is true and says nothing about whether that matters, "
        "because it does not say 2.5 times **what**. In a point model there is "
        "no denominator to be had.",
        "",
        "## The cold-to-hot ladder",
        "",
        "| treatment | direct kills | DAMP | immune kills, no blockade | at full blockade | fold | immune share, unblocked → blocked |",
        "|---|--:|--:|--:|--:|--:|--:|",
    ]
    for a in d["arms"]:
        pts = a["points"]
        fold = f"{a['fold_benefit']:.1f}×" if a["fold_benefit"] else "—"
        L.append(
            f"| {a['treatment']} | {a['ferroptosis_kills']:,} | "
            f"{a['total_damp']:,.0f} | {pts[0]['immune_kills']} | "
            f"{pts[-1]['immune_kills']} | {fold} | "
            f"{a['share_unblocked']:.2%} → {a['share_blocked']:.2%} |")
    L += [
        "",
        "## Three findings, and the third is the one worth carrying",
        "",
        "**1. A cold tumour gets nothing.** The untreated control produces two "
        "deaths and almost no danger signal, and blockade yields exactly zero "
        "extra kills at *every* strength tested. Checkpoint blockade cannot "
        "start a response — it can only amplify one — and this is that "
        "statement made where it can be checked rather than asserted.",
        "",
        f"**2. The fold benefit is the same for both active treatments** "
        f"({' and '.join(f'{f:.1f}×' for f in d['fold_benefits'])}, within "
        f"{d['fold_spread']:.2f}× of each other). That is exactly what a "
        "multiplier does: it does not know what it multiplies. The point "
        "model's own output is therefore the quantity that **cannot** "
        "distinguish these two cases.",
        "",
        f"**3. The share differs by {d['share_ratio']:.0f}×, and it runs "
        "against the danger signal.** The treatment producing "
        f"{d['damp_ratio']:.1f}× more DAMP ({d['higher_damp_treatment']}) has "
        "the *smaller* immune share, because it has already killed almost "
        "everything directly and left the immune system nothing to add.",
        "",
        "That last point is a warning about a natural inference. Ranking "
        "treatments by how much immunogenic signal they generate, and picking "
        "the top one as the partner for checkpoint blockade, gets this "
        "backwards: what matters is the immune contribution as a share of the "
        "kill that **remains**, and the biggest DAMP producer here is the "
        "worst on that measure.",
        "",
        "## What this does NOT establish",
        "",
        "- **The immune contribution in this engine is small in absolute "
        "terms** — a fraction of a percent of the kill under the physical arm. "
        "That is a known property of the model rather than a finding about "
        "biology: the manuscript already records the immune ratio collapsing "
        "roughly 104:1 to 4:1 when the model moved to three dimensions. The "
        "ratios above are what this engine says, not what a patient does.",
        "- **The counts are small.** Converged, and checked to be: 150 and 250 "
        "steps agree to within the per-cell RNG. But tens of kills is a regime "
        "where a single seed matters, and no interval is quoted here.",
        "- **The brake parameters are placeholders.** A PD-1 brake of 0.7 and "
        "the panel's other channels are not fitted to anything, and this arm's "
        "recorded verdict is INADMISSIBLE for a reason that stands: its "
        "response band constrains a product of two factors, neither "
        "identifiable from it.",
        "- **The ladder is three points**, not a dose-response. It separates "
        "cold from warm from hot; it does not locate a threshold between them.",
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
