"""Where a bystander payload is worth most, which a scalar cannot say (#844).

WHAT THE POINT MODEL CANNOT EXPRESS
------------------------------------
`adc::bystander_kill_fraction` takes `neighbours_in_reach` as a bare scalar.
That is not a small simplification: the payload has nowhere to go and no
distance to travel, so the arm cannot state the reason bystander payloads exist
at all. An antibody is large and penetrates a tumour badly; the payload it
releases is small and does not. The payload's job is to reach cells the
conjugate never did.

The prediction that follows, and which needs a "where": the bystander effect's
MARGINAL value should be largest exactly where the antibody penetrates worst,
and should fall away as penetration improves -- a conjugate that reaches
everything has little left for the payload to add.

HOW IT IS MEASURED
------------------
Every condition is run twice, identical but for the linker: cleavable (payload
escapes) and non-cleavable (it does not). The difference IS the bystander's
contribution, measured rather than parameterised, and the pair is matched on
everything else including the RNG stream.

TWO COLUMNS, BECAUSE THEY DISAGREE AND ONE OF THEM FLATTERS
------------------------------------------------------------
The FOLD advantage -- total kill with the payload over total kill without --
falls as penetration improves, which is the prediction. The ABSOLUTE bystander
count RISES over the same range, because a better-penetrating conjugate makes
more dying cells and therefore more payload sources. Reporting only the ratio
would publish the more flattering of two true answers, which this repository
has done before and had to retract.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SWEEP = REPO / "analysis" / "calibration" / "adc_bystander_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "adc.rs"
OUT_MD = REPO / "analysis" / "calibration" / "adc-bystander-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "adc-bystander-validation.json"


def _rows():
    out = []
    for ln in SWEEP.read_text().splitlines():
        if not ln.startswith("ADC_BYSTANDER "):
            continue
        d = dict(re.findall(r"(\w+)=([-\d.eE]+|true|false)", ln))
        out.append({k: (v == "true") if v in ("true", "false") else float(v)
                    for k, v in d.items()})
    return out


def scan() -> dict:
    rows = _rows()
    lambdas = sorted({r["antibody_lambda_um"] for r in rows})
    reaches = sorted({int(r["payload_reach_cells"]) for r in rows})
    by_reach = []
    for reach in reaches:
        pts = []
        for lam in lambdas:
            on = next(r for r in rows if r["antibody_lambda_um"] == lam
                      and int(r["payload_reach_cells"]) == reach and r["cleavable"])
            off = next(r for r in rows if r["antibody_lambda_um"] == lam
                       and int(r["payload_reach_cells"]) == reach and not r["cleavable"])
            pts.append({
                "antibody_lambda_um": lam,
                "direct_kills": int(off["total_dead"]),
                "with_payload": int(on["total_dead"]),
                "bystander_kills": int(on["bystander_kills"]),
                "fold_advantage": on["total_dead"] / off["total_dead"]
                                  if off["total_dead"] else None,
                # The scalar the point model would have returned, which does not
                # move with penetration at all because it has no term for it.
                "scalar_bystander_fraction": on["scalar_bystander_fraction"],
            })
        by_reach.append({"payload_reach_cells": reach, "points": pts})
    return {"by_reach": by_reach, "lambdas": lambdas}


def assemble(raw: dict) -> dict:
    d = dict(raw)
    checks = []
    for r in raw["by_reach"]:
        folds = [p["fold_advantage"] for p in r["points"]]
        absol = [p["bystander_kills"] for p in r["points"]]
        scalars = {p["scalar_bystander_fraction"] for p in r["points"]}
        checks.append({
            "payload_reach_cells": r["payload_reach_cells"],
            "fold_falls_with_penetration": all(a >= b for a, b in zip(folds, folds[1:])),
            "absolute_rises_with_penetration": all(a <= b for a, b in zip(absol, absol[1:])),
            "fold_range": [min(folds), max(folds)],
            "scalar_is_constant": len(scalars) == 1,
        })
    d["checks"] = checks
    d["fold_falls_everywhere"] = all(c["fold_falls_with_penetration"] for c in checks)
    d["absolute_rises_everywhere"] = all(c["absolute_rises_with_penetration"] for c in checks)
    d["scalar_never_moves"] = all(c["scalar_is_constant"] for c in checks)
    # THE PREDICTION HAS A CONDITION, and the sweep found it rather than the
    # author choosing it. At a payload reach of one cell the advantage is flat
    # in penetration: the payload lands only on immediate neighbours, which are
    # cells the conjugate very largely reached anyway, so there is no
    # penetration-dependent gap for it to fill. Reporting the reach-2 row alone
    # would have been picking the arm that agreed.
    d["reaches_where_fold_falls"] = sorted(
        c["payload_reach_cells"] for c in checks if c["fold_falls_with_penetration"])
    d["reaches_where_fold_is_flat"] = sorted(
        c["payload_reach_cells"] for c in checks
        if not c["fold_falls_with_penetration"])
    flat_spans = [c["fold_range"] for c in checks
                  if not c["fold_falls_with_penetration"]]
    d["flat_arm_span_ratio"] = (max(hi / lo for lo, hi in flat_spans)
                                if flat_spans else None)
    d["verdict"] = (
        "PARTIAL — the marginal value falls as penetration improves, but only "
        "where the payload reaches beyond one cell"
        if d["reaches_where_fold_falls"] and d["reaches_where_fold_is_flat"]
        and d["scalar_never_moves"]
        else "CONFIRMED — the marginal value falls as penetration improves"
        if d["fold_falls_everywhere"] and d["scalar_never_moves"]
        else "UNRESOLVED")
    return d


def render(d: dict) -> str:
    L = [
        "# Antibody-drug conjugate: where a bystander payload is worth most (#844)",
        "",
        "*Generated by `scripts/validate_adc_bystander.py --render-only`. Pure "
        "stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}.**",
        "",
        "## What a scalar cannot say",
        "",
        "`adc::bystander_kill_fraction` takes `neighbours_in_reach` as a bare "
        "number. The payload has nowhere to go and no distance to travel, so "
        "the arm cannot state the reason bystander payloads exist: an antibody "
        "is large and penetrates a tumour badly, the payload it releases is "
        "small and does not, and the payload's job is to reach cells the "
        "conjugate never did.",
        "",
        f"The scalar's value is the same at every penetration length in this "
        f"sweep ({'confirmed' if d['scalar_never_moves'] else 'NOT constant'}) — "
        "it has no term for penetration, so it cannot vary with it.",
        "",
        "## Measured by difference, not by parameter",
        "",
        "Every condition runs twice, identical but for the linker: cleavable "
        "(payload escapes) and non-cleavable (it does not). The difference **is** "
        "the bystander's contribution.",
        "",
    ]
    for r in d["by_reach"]:
        L += [
            f"### Payload reach {r['payload_reach_cells']} cell(s)",
            "",
            "| antibody λ | kill without payload | with payload | fold advantage | bystander kills |",
            "|--:|--:|--:|--:|--:|",
        ]
        for p in r["points"]:
            L.append(f"| {p['antibody_lambda_um']:.0f} µm | {p['direct_kills']:,} | "
                     f"{p['with_payload']:,} | **{p['fold_advantage']:.2f}×** | "
                     f"{p['bystander_kills']:,} |")
        L.append("")
    demo = next(r for r in d["by_reach"]
                if r["payload_reach_cells"] in d["reaches_where_fold_falls"])
    first = demo["points"]
    L += [
        "## The finding, and the condition on it",
        "",
        f"**The fold advantage falls as the antibody penetrates better** — at a "
        f"payload reach of {demo['payload_reach_cells']} cells, "
        f"{first[0]['fold_advantage']:.2f}× at λ = {first[0]['antibody_lambda_um']:.0f} µm "
        f"down to {first[-1]['fold_advantage']:.2f}× at "
        f"{first[-1]['antibody_lambda_um']:.0f} µm, monotonically. That is the "
        "prediction: a bystander payload multiplies the conjugate's reach most "
        "when the conjugate reaches least, and a conjugate that already "
        "reaches everything has little left for the payload to add.",
        "",
        f"**But not at every reach.** It holds at "
        + ", ".join(f"{v} cells" for v in d["reaches_where_fold_falls"])
        + f" and is FLAT at "
        + ", ".join(f"{v} cell" for v in d["reaches_where_fold_is_flat"])
        + f" — a span of only {d['flat_arm_span_ratio']:.2f}× across a "
        "sixteenfold range of antibody penetration. The reason is the "
        "prediction's own logic running out: a payload that reaches only its "
        "immediate neighbours lands almost entirely on cells the conjugate "
        "reached anyway, so there is no penetration-dependent gap left for it "
        "to fill. **The sweep found that condition; it was not chosen.** "
        "Reporting the agreeing row alone would have been picking the arm that "
        "agreed.",
        "",
        "**And the absolute bystander count moves the other way where it "
        "moves at all** — it RISES with penetration at "
        + ", ".join(f"{c['payload_reach_cells']} cell(s)" for c in d["checks"]
                    if c["absolute_rises_with_penetration"])
        + (", and does not at "
           + ", ".join(f"{c['payload_reach_cells']} cell(s)" for c in d["checks"]
                       if not c["absolute_rises_with_penetration"])
           if not d["absolute_rises_everywhere"] else "")
        + ". A better-penetrating conjugate kills more cells and so makes more "
        "payload sources, which pushes the count up while the ratio falls. "
        "Both are true of the same runs and they point opposite ways; "
        "reporting only the ratio would publish the more flattering of two "
        "answers, which this repository has done before and had to retract, so "
        "both columns are in the tables above.",
        "",
        "## What this does NOT establish",
        "",
        "- **No calibration.** This arm's fitted verdict is NO TARGET and stays "
        "there: no published bystander dose-response exists to fit. The kill "
        "magnitudes here are illustrative and uncalibrated; what is measured is "
        "how the bystander's contribution MOVES with penetration.",
        "- **The payload reach is a parameter.** Nothing here measures how far "
        "a released payload actually diffuses before it is cleared or bound.",
        "- **The antibody field is the radial-depth proxy**, the same geometry "
        "the oxygen field uses with a shorter length. A real conjugate's "
        "distribution is shaped by binding-site consumption as well as "
        "diffusion, which this does not represent.",
        "- **Payload release is one hop, deliberately.** A bystander-killed "
        "cell releases nothing, because the payload came from the conjugate and "
        "a cell that never took up conjugate has none to give. Letting it "
        "release would be a chain reaction, not a bystander effect.",
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
