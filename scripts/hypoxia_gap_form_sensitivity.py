#!/usr/bin/env python3
"""The O2-FORM axis the §7.1 kill-collapse interval holds fixed (#925).

`headline_uncertainty.py --headline sim-tme` publishes a 95% prior-predictive
interval for the §7.1 kill-collapse gap (SDT minus RSL3 hypoxic-zone kill):
~[0.163, 0.999], median 0.960, default point 0.866. That interval is computed by
sampling the 11 PRCC RATE CONSTANTS through `FERRO_PARAM_OVERRIDES`.

IT NEVER PASSES `--sdt-o2-dependence`. The invocation is `[str(binary)]` with no
arguments, so every draw runs at the default `sdt_o2_dependence = 0.0`, where
SDT keeps its full hypoxic efficacy -- the single O2 setting most favourable to
a positive gap. The form is not in the swept set, so the interval cannot see it.

WHAT THIS MEASURES, AND WHY IT IS NOT A RESTATEMENT. `identifiability-report.md`
already brackets the assumption: it records that the magnitude "collapses from
~87.8% to ~0% under full SDT O2-dependence". So the COLLAPSE was known. What was
never measured is the SIGN at the collapsed end, and the verdict label turns on
exactly that -- the hypoxia row is published as `direction_robust_magnitude_not`,
whose rationale reads "the asymmetry SIGN is parameter-robust, but the magnitude
is ... assumption-dependent". A difference that collapses to zero does not have a
robust sign to report, and measured at full dependence the gap is NEGATIVE.

So this file supplies one number the repo did not have, and the verdict is
derived from it rather than asserted: `identifiability_report.py` reads this
artifact and refuses to call the direction robust when the measured gap changes
sign across the form axis.

SCOPE, stated rather than implied:
  - sim-tme carries only the LINEAR `o2_dependent_exo_factor`. The measured
    Alper-Howard-Flanders hyperbola (#726, `oer_exo_factor`) is NOT expressible
    in this binary, so the form axis sampled here is off-to-linear, which is the
    axis the published interval holds fixed. sim-tme-3d covers the OER form
    (`the_hypoxia_leg_is_decided_by_the_o2_shapes_and_one_pairing_inverts_it`).
  - Rate constants stay at their DEFAULTS here. This varies the form alone, so
    it is the complement of the #332 sweep rather than a replacement: neither
    crosses the other, and a joint sweep is not attempted.
  - `sdt_ros` remains the uncalibrated #336 knob. What is varied is the SHAPE of
    the oxygen response, not the dose.

OFFLINE CONTRACT: the measurement needs the compiled `sim-tme` (~4 min/run), so
it is taken locally and the derived artifact is committed. CI reads only the
JSON. `--render-only` rebuilds the prose from it without the binary.

Usage:
  hypoxia_gap_form_sensitivity.py [--binary PATH] [--render-only]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from headline_sensitivity import extract_tme_observables  # noqa: E402

OUT_JSON = REPO / "analysis" / "calibration" / "hypoxia-gap-form-sensitivity.json"
OUT_MD = REPO / "analysis" / "calibration" / "hypoxia-gap-form-sensitivity.md"

# The form axis. 0.0 is what the #332 interval used; 1.0 is the fully
# oxygen-dependent Type II sonosensitizer the lead clinical agent occupies.
DEPENDENCES = [0.0, 0.25, 0.5, 0.75, 1.0]

# The published interval this complements, quoted for comparison only.
PUBLISHED = {
    "interval_95": [0.163, 0.999],
    "median": 0.960,
    "default_point": 0.866,
    "source": "analysis/headline-uncertainty-tme-report.md (#332)",
}


def default_binary() -> Path:
    env = os.getenv("CARGO_TARGET_DIR")
    roots = [Path(env)] if env else []
    roots += [REPO / "simulations" / "target", Path.home() / "ironauth-target"]
    for r in roots:
        cand = r / "release" / "sim-tme"
        if cand.exists():
            return cand
    return roots[0] / "release" / "sim-tme"


def measure_one(binary: Path, dependence: float) -> float:
    """The §7.1 gap at one O2-dependence, via the repo's OWN extractor.

    `extract_tme_observables` is imported rather than reimplemented: a local
    copy of the arithmetic would let the observable drift away from the one the
    published interval measured, which is the whole basis of the comparison.
    """
    args = [str(binary)]
    if dependence > 0.0:
        args += ["--sdt-o2-dependence", str(dependence)]
    with tempfile.TemporaryDirectory(prefix="ferro_gapform_") as wd:
        proc = subprocess.run(args, cwd=wd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"sim-tme failed at dependence={dependence} "
                f"({proc.returncode}): {proc.stderr[-400:]}"
            )
        summary = Path(wd) / "output" / "tme" / "tme_summary.json"
        conditions = json.loads(summary.read_text())["conditions"]
        return float(extract_tme_observables(conditions)["hypoxia"])


def assemble(raw: dict) -> dict:
    """Derive every reported field from the measurements. Idempotent, so
    `--render-only` reproduces a fresh run exactly rather than echoing stored
    prose (the `test_artifact_freshness` lesson)."""
    pts = {float(k): float(v) for k, v in raw["gap_by_dependence"].items()}
    gaps = [pts[d] for d in sorted(pts)]
    lo, hi = PUBLISHED["interval_95"]
    out = dict(raw)
    out["gap_by_dependence"] = {f"{d:g}": pts[d] for d in sorted(pts)}
    out["published"] = PUBLISHED
    out["gap_at_zero"] = gaps[0]
    out["gap_at_full"] = gaps[-1]
    out["sign_changes_across_form_axis"] = (gaps[0] > 0.0) and (gaps[-1] <= 0.0)
    out["full_dependence_outside_published_interval"] = not (lo <= gaps[-1] <= hi)
    out["monotone_decreasing"] = all(
        b <= a + 1e-12 for a, b in zip(gaps, gaps[1:])
    )
    out["reproduces_published_point"] = abs(gaps[0] - PUBLISHED["default_point"]) < 5e-3
    return out


def render(d: dict) -> str:
    pts = {float(k): v for k, v in d["gap_by_dependence"].items()}
    lo, hi = d["published"]["interval_95"]
    rows = "\n".join(
        f"| {dep:g} | {pts[dep]:+.6f} | "
        f"{'inside' if lo <= pts[dep] <= hi else '**outside**'} |"
        for dep in sorted(pts)
    )
    sign = d["sign_changes_across_form_axis"]
    L = []
    L.append("# The O2-form axis the §7.1 interval holds fixed\n")
    L.append(
        "Generated by `scripts/hypoxia_gap_form_sensitivity.py`. The §7.1 "
        "kill-collapse gap is `SDT - RSL3` hypoxic-zone kill, extracted by the "
        "same `extract_tme_observables` the published interval uses.\n"
    )
    L.append(
        f"The published 95% prior-predictive interval is "
        f"[{lo}, {hi}] (median {d['published']['median']}, default point "
        f"{d['published']['default_point']}), computed by sampling the 11 PRCC "
        "rate constants. **It never passes `--sdt-o2-dependence`**, so every "
        "draw sits at `0.0`, where SDT keeps full hypoxic efficacy.\n"
    )
    L.append("## Measured\n")
    L.append("| `--sdt-o2-dependence` | §7.1 gap | vs published interval |")
    L.append("|---|--:|---|")
    L.append(rows + "\n")
    if d["reproduces_published_point"]:
        L.append(
            f"The zero-dependence row reproduces the published default point "
            f"({d['published']['default_point']}) to within 5e-3, so these runs "
            "are comparable with the interval rather than a different "
            "measurement.\n"
        )
    else:
        L.append(
            "**The zero-dependence row does NOT reproduce the published default "
            "point**, so this table is not comparable with the interval and "
            "nothing should be concluded from it until that is explained.\n"
        )
    L.append("## What it means\n")
    if sign:
        L.append(
            f"**The gap changes sign across the form axis** "
            f"({d['gap_at_zero']:+.6f} at zero dependence, "
            f"{d['gap_at_full']:+.6f} at full), and the full-dependence value "
            f"falls {'outside' if d['full_dependence_outside_published_interval'] else 'inside'} "
            "the published interval. A difference that collapses through zero "
            "does not have a robust sign to report, so the hypoxia row's "
            "`direction_robust` reading is scoped to the rate constants at a "
            "fixed form, not to the form itself.\n"
        )
    else:
        L.append(
            f"The gap does NOT change sign across the form axis "
            f"({d['gap_at_zero']:+.6f} -> {d['gap_at_full']:+.6f}), so the "
            "direction survives this axis and the `direction_robust` reading "
            "stands on it.\n"
        )
    L.append(
        "The magnitude collapse was already recorded in "
        "`identifiability-report.md` (\"~87.8% to ~0%\"). What is added here is "
        "the **sign** at the collapsed end, which is what the verdict label "
        "turns on.\n"
    )
    L.append("## Scope\n")
    L.append(
        "- sim-tme carries only the LINEAR `o2_dependent_exo_factor`; the "
        "measured OER hyperbola (#726) is not expressible in this binary, so "
        "the axis here is off-to-linear. sim-tme-3d covers the OER form.\n"
        "- Rate constants stay at their defaults, so this varies the form "
        "alone and is the complement of the #332 sweep, not a replacement. No "
        "joint sweep is attempted.\n"
        "- `sdt_ros` remains the uncalibrated #336 knob. What varies is the "
        "SHAPE of the oxygen response, not the dose.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the prose from the committed JSON, no binary")
    a = ap.parse_args()

    if a.render_only:
        # RE-ASSEMBLED, not echoed. The stored file carries derived fields
        # alongside the raw measurements, and re-rendering them instead of
        # recomputing would make every guard on a derived field inert.
        data = assemble(json.loads(OUT_JSON.read_text()))
    else:
        binary = a.binary or default_binary()
        if not binary.exists():
            print(f"error: sim-tme not built at {binary}; "
                  f"cargo build --release -p sim-tme", file=sys.stderr)
            return 2
        gaps = {}
        for dep in DEPENDENCES:
            gaps[f"{dep:g}"] = measure_one(binary, dep)
            print(f"  --sdt-o2-dependence {dep:g}: gap {gaps[f'{dep:g}']:+.6f}",
                  flush=True)
        data = assemble({"gap_by_dependence": gaps})
    OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_JSON.relative_to(REPO)} and {OUT_MD.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
