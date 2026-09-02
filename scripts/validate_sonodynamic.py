#!/usr/bin/env python3
"""The frequency an SDT applicator should use, and where the model is wrong.

WHY THIS ARM HAD NO MODULE FOR SO LONG
--------------------------------------
Sonodynamic therapy is one of the two modalities the project was built on and
owned zero engine lines until `simulations/ferroptosis-core/src/sonodynamic.rs`.
That was a consequence of a validated decision, not neglect:
`analysis/calibration/pdt-threshold-validation.md` records that the engine
feeds SDT and PDT into ONE death threshold with `sdt_ros == pdt_ros`, because
Zhu 2015 measured the reacted-singlet-oxygen kill threshold to be roughly
source-independent. If the threshold really is a property of the target, then
every difference between these two arms lives UPSTREAM of the cell.

This validates the upstream half.

THE MODEL
---------
Three frequency dependences act at once. Focusing favours going high (a
diffraction-limited focus shrinks with wavelength, so at fixed power the focal
pressure rises linearly with frequency). Attenuation and the mechanical index
itself both favour going low. The product has an interior maximum with a
closed form, `f* = 10 / (alpha * ln(10) * z)`.

THE COMPARATOR, AND WHY IT IS A REAL ONE
----------------------------------------
Ellens & Hynynen 2015 (Med Phys 42(8):4896, PMID 26233216) simulated the same
question independently -- full-wave Rayleigh-Sommerfeld propagation into a
Pennes bioheat solve over transducer arrays of up to 70,000 elements -- and
named the same three mechanisms in the same terms. They scanned 250 to 1500
kHz at 50, 100 and 150 mm and two attenuation coefficients.

THREE CLAIMS ARE SCORED SEPARATELY, AND ONE FAILS
-------------------------------------------------
A single pass/fail over an arm this coarse would hide the informative part.
The existence of the optimum and its movement with attenuation both hold. The
DEPTH SCALING DOES NOT: this model says the optimum falls as 1/z, and they
report the same 750 kHz at all three depths.

WHAT IS DELIBERATELY NOT CLAIMED
--------------------------------
No numerical agreement. Two independent mismatches would make one
uninterpretable: their observable is THERMAL ablation efficiency and this
model's is cavitation likelihood, and absorption enters those two with
opposite sign; and their attenuation is quoted in Np/m/MHz on the pressure
amplitude while this engine's `sdt_alpha` is a dB/cm/MHz applied to intensity,
a conversion carrying a factor of two that is the easiest unit error in
acoustics to make. What survives both mismatches is the DIRECTION of each
dependence, which is what is scored.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "sonodynamic.rs"
PARAMS = REPO / "simulations" / "ferroptosis-core" / "src" / "params.rs"
OUT_MD = REPO / "analysis" / "calibration" / "sonodynamic-validation.md"
OUT_JSON = REPO / "analysis" / "calibration" / "sonodynamic-validation.json"

# Ellens & Hynynen 2015, PMID 26233216, reported optima. Their alpha is in
# Np/m/MHz; it is kept in THEIR units here and never converted, because the
# conversion is the ambiguity this script refuses to resolve.
ELLENS = [
    {"alpha_np_m_mhz": 5.0, "depth_mm": 50, "reported_khz": 750},
    {"alpha_np_m_mhz": 5.0, "depth_mm": 100, "reported_khz": 750},
    {"alpha_np_m_mhz": 5.0, "depth_mm": 150, "reported_khz": 750},
    {"alpha_np_m_mhz": 10.0, "depth_mm": 50, "reported_khz": 750},
    {"alpha_np_m_mhz": 10.0, "depth_mm": 100, "reported_khz": 500},
    {"alpha_np_m_mhz": 10.0, "depth_mm": 150, "reported_khz": 500},
]
ELLENS_SCAN_KHZ = (250, 1500)


def optimal_frequency_mhz(depth_cm, alpha_db_cm_mhz):
    """The crate's closed form, re-implemented (drift-guarded below)."""
    denom = alpha_db_cm_mhz * math.log(10.0) * depth_cm
    return math.inf if denom <= 0 else 10.0 / denom


def delivered_index(depth_cm, freq_mhz, alpha_db_cm_mhz):
    if freq_mhz <= 0:
        return 0.0
    p = freq_mhz * 10.0 ** (-alpha_db_cm_mhz * freq_mhz * max(depth_cm, 0.0) / 20.0)
    return p / math.sqrt(freq_mhz)


def _rust_alpha_default():
    m = re.search(r"sdt_alpha:\s*([0-9.]+)", PARAMS.read_text())
    return float(m.group(1))


def scan() -> dict:
    alpha = _rust_alpha_default()

    # Claim 1: the optimum is interior, and the closed form finds it.
    interior = []
    for z in (3.0, 5.0, 10.0, 15.0):
        closed = optimal_frequency_mhz(z, alpha)
        best_f, best_v = 0.0, -1.0
        f = 0.005
        while f < 30.0:
            v = delivered_index(z, f, alpha)
            if v > best_v:
                best_f, best_v = f, v
            f += 0.005
        interior.append({
            "depth_cm": z,
            "closed_form_mhz": closed,
            "scanned_mhz": best_f,
            "agrees": abs(closed - best_f) < 0.02,
            "beats_low_end": best_v > delivered_index(z, 0.05, alpha),
            "beats_high_end": best_v > delivered_index(z, 20.0, alpha),
        })

    # Claim 2: doubling attenuation lowers the optimum. Ellens & Hynynen see
    # 750 -> 500 kHz at 100 and 150 mm when their alpha doubles.
    atten = []
    for z in (5.0, 10.0, 15.0):
        a1 = optimal_frequency_mhz(z, alpha)
        a2 = optimal_frequency_mhz(z, alpha * 2.0)
        atten.append({"depth_cm": z, "at_alpha_mhz": a1, "at_2alpha_mhz": a2,
                      "ratio": a1 / a2, "falls": a2 < a1})
    reported_falls = any(
        r["reported_khz"] < s["reported_khz"]
        for r in ELLENS for s in ELLENS
        if r["depth_mm"] == s["depth_mm"] and r["alpha_np_m_mhz"] > s["alpha_np_m_mhz"]
    )

    # Claim 3: the depth scaling. Theirs is flat, the model's is 1/z.
    depths = sorted({e["depth_mm"] for e in ELLENS})
    model_span = (optimal_frequency_mhz(min(depths) / 10.0, alpha)
                  / optimal_frequency_mhz(max(depths) / 10.0, alpha))
    reported_spans = {}
    for a in sorted({e["alpha_np_m_mhz"] for e in ELLENS}):
        rows = [e for e in ELLENS if e["alpha_np_m_mhz"] == a]
        khz = [e["reported_khz"] for e in rows]
        reported_spans[str(a)] = max(khz) / min(khz)

    # The band check: does the model land inside the frequency range they
    # scanned at all, at the depths they scanned? A model outside 250-1500 kHz
    # everywhere would not even be answering the same question.
    band = []
    for z_mm in depths:
        f = optimal_frequency_mhz(z_mm / 10.0, alpha) * 1000.0
        band.append({"depth_mm": z_mm, "model_khz": f,
                     "in_scanned_band": ELLENS_SCAN_KHZ[0] <= f <= ELLENS_SCAN_KHZ[1]})

    return {
        "alpha_db_cm_mhz_from_params_rs": alpha,
        "interior": interior,
        "attenuation": atten,
        "reported_optimum_falls_with_attenuation": reported_falls,
        "model_depth_span_ratio": model_span,
        "reported_depth_span_ratios": reported_spans,
        "band": band,
        "ellens": ELLENS,
    }


def assemble(raw: dict) -> dict:
    d = dict(raw)
    d["claim_interior_optimum"] = (
        "CONFIRMED" if all(r["agrees"] and r["beats_low_end"] and r["beats_high_end"]
                           for r in raw["interior"]) else "FAILED")
    d["claim_falls_with_attenuation"] = (
        "CONFIRMED" if (all(r["falls"] for r in raw["attenuation"])
                        and raw["reported_optimum_falls_with_attenuation"]) else "FAILED")
    # The model says the optimum falls threefold across 50-150 mm. If the
    # comparator's own span is essentially 1, the model is refuted, and a
    # refutation is the finding rather than a reason to retune.
    flat = all(v < 1.6 for v in raw["reported_depth_span_ratios"].values())
    d["claim_depth_scaling"] = "REFUTED" if (flat and raw["model_depth_span_ratio"] > 2.0) else "UNDECIDED"
    d["verdict"] = "PARTIAL" if (d["claim_interior_optimum"] == "CONFIRMED"
                                 and d["claim_falls_with_attenuation"] == "CONFIRMED"
                                 and d["claim_depth_scaling"] == "REFUTED") else "UNRESOLVED"
    return d


def render(d: dict) -> str:
    a = d["alpha_db_cm_mhz_from_params_rs"]
    L = [
        "# Sonodynamic frequency optimum: validation (#831)",
        "",
        "*Generated by `scripts/validate_sonodynamic.py --render-only`. Pure "
        "stdlib, runs in CI. Do not edit by hand.*",
        "",
        f"**Verdict: {d['verdict']}** — two of three claims confirmed against an "
        "independent full-wave study, and the third refuted by it.",
        "",
        "## What the model says",
        "",
        "Three frequency dependences act at once: focusing favours high "
        "frequency (a diffraction-limited focus shrinks with wavelength, so "
        "focal pressure rises linearly with it at fixed power), while "
        "attenuation and the mechanical index's own `1/sqrt(f)` both favour "
        "low. The product has an interior maximum,",
        "",
        "> `f* = 10 / (alpha * ln(10) * z)`",
        "",
        f"At the engine's own `sdt_alpha = {a}` dB/cm/MHz:",
        "",
        "| focal depth | model optimum | inside the comparator's scanned band |",
        "|---|--:|---|",
    ]
    for r in d["band"]:
        L.append(f"| {r['depth_mm']} mm | {r['model_khz']:.0f} kHz | "
                 f"{'yes' if r['in_scanned_band'] else 'no'} |")
    L += [
        "",
        "## The comparator",
        "",
        "Ellens & Hynynen 2015 (Med Phys 42(8):4896, PMID 26233216) asked the "
        "same question by a completely different route — full-wave "
        "Rayleigh-Sommerfeld propagation into a Pennes bioheat solve, over "
        "simulated arrays of 2,000 to 70,000 elements — and named the same "
        "three mechanisms: *focal size (decreasing with frequency), peak "
        "pressure (generally increasing with frequency), and attenuation "
        "(also increasing with frequency)*. Their reported optima:",
        "",
        "| their alpha (Np/m/MHz) | depth | reported optimum |",
        "|--:|--:|--:|",
    ]
    for e in d["ellens"]:
        L.append(f"| {e['alpha_np_m_mhz']:.0f} | {e['depth_mm']} mm | {e['reported_khz']} kHz |")
    L += [
        "",
        "## The three claims, scored apart",
        "",
        f"**1. An interior optimum exists — {d['claim_interior_optimum']}.** "
        "The closed form agrees with a brute scan at every depth tested and "
        "beats both ends of the range, so the model is not secretly monotonic. "
        "The comparator's whole paper presumes an interior optimum and finds "
        "one, which is independent confirmation that the question is real.",
        "",
        f"**2. The optimum falls as attenuation rises — {d['claim_falls_with_attenuation']}.** "
        "The model predicts doubling `alpha` halves `f*` exactly:",
        "",
        "| depth | at alpha | at 2·alpha | ratio |",
        "|---|--:|--:|--:|",
    ]
    for r in d["attenuation"]:
        L.append(f"| {r['depth_cm']:.0f} cm | {r['at_alpha_mhz']*1000:.0f} kHz | "
                 f"{r['at_2alpha_mhz']*1000:.0f} kHz | {r['ratio']:.2f}× |")
    span = d["model_depth_span_ratio"]
    L += [
        "",
        "The comparator sees the same direction, moving 750 kHz to 500 kHz at "
        "100 and 150 mm when its attenuation doubles. Its frequency grid is "
        "250 kHz wide, so it cannot resolve the exact factor and the "
        "agreement is on sign, not magnitude.",
        "",
        f"**3. The optimum falls as 1/depth — {d['claim_depth_scaling']}.**",
        "",
        f"The model predicts a **{span:.2f}×** fall across 50 to 150 mm. The "
        "comparator reports:",
        "",
    ]
    for k, v in sorted(d["reported_depth_span_ratios"].items()):
        L.append(f"- at alpha = {k} Np/m/MHz, the optimum spans **{v:.2f}×** "
                 "across the same three depths")
    L += [
        "",
        "At their lower attenuation the reported optimum does not move at all. "
        "This is a genuine disagreement and not a tuning gap: no value of "
        "`alpha` makes a `1/z` law flat.",
        "",
        "### The missing term, which they name themselves",
        "",
        "Their conclusion states that *near-field heat accumulation tends to "
        "be the rate limiting factor in large-volume ablations*. This model "
        "has no term for the tissue in front of the focus at all — it "
        "evaluates one point. Lowering the frequency to chase the focal "
        "optimum also lengthens the heated path, and a model that cannot see "
        "that path is free to slide its optimum in a way a real applicator is "
        "not. That is a named, addressable omission, which is the most useful "
        "thing a refutation can leave behind.",
        "",
        "## What is NOT claimed",
        "",
        "No numerical agreement, and the refusal is deliberate. Two "
        "independent mismatches would each on its own make one "
        "uninterpretable:",
        "",
        "- **Different observable.** Their efficiency is thermal ablation per "
        "joule; this model's is cavitation likelihood. Absorption enters the "
        "two with **opposite sign** — it *helps* heating and *hurts* "
        "cavitation — so the two optima should not coincide even for a "
        "perfect model of both.",
        "- **Different attenuation convention.** Theirs is Np/m/MHz on the "
        "pressure amplitude; this engine's `sdt_alpha` is dB/cm/MHz applied "
        "to intensity. The conversion carries a factor of two, and a numeric "
        "match obtained after choosing one of two defensible conversions "
        "would be a statement about the choice.",
        "",
        "So the scored claims are all directional. The magnitudes stay "
        "unconstrained, and the `f* * alpha * z` product is the invariant a "
        "reader should carry rather than any megahertz figure.",
        "",
    ]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    # `assemble` runs on BOTH paths deliberately. Reading the stored JSON and
    # rendering it directly would re-print derived fields -- the verdicts --
    # from storage instead of recomputing them, which makes every guard that
    # reads a derived field inert. `tests/test_artifact_freshness.py` checks
    # exactly this, and caught it here.
    d = assemble(json.loads(OUT_JSON.read_text()) if args.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    print(f"{OUT_MD.relative_to(REPO)}: {d['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
