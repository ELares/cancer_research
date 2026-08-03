#!/usr/bin/env python3
"""Seed replication for sim-tme: dispersion for numbers that had none (#SEED-REP).

WHY
---
`sim-tme` ran with a hard-coded `SEED = 42` and no replicate loop anywhere, so
every number this engine contributes to the manuscript is a SINGLE stochastic
draw reported as a point estimate. Several headline quantities rest on small
event counts -- the immune-coupling ratio has a denominator of five simulated
events -- where one draw carries almost no information, and the depth-kill table
already shows the noise floor (RSL3, depth-independent by construction, swings
1.2-4.9%).

This runs the same matrix across N seeds via the `FERRO_SEED` override and
reports medians with percentile bootstrap intervals, so a reader can tell which
differences survive resampling and which are noise.

`FERRO_SEED` unset reproduces seed 42 byte-identically, so the committed
baseline and the #253 regression hash are untouched; replicate runs are written
to a scratch directory and never overwrite `output/tme`.

Usage:
    python scripts/seed_replication.py --seeds 20
    python scripts/seed_replication.py --seeds 5 --jobs 3     # quick look
"""

import argparse
import json
import os
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

BINARY = PROJECT_ROOT / "simulations" / "target" / "release" / "sim-tme"
OUT = PROJECT_ROOT / "analysis" / "seed-replication-report.md"
RAW = PROJECT_ROOT / "analysis" / "seed-replication.json"

# The quantities the manuscript actually quotes from this engine.
METRICS = ["overall_kill_rate", "hypoxic_kill_rate", "normoxic_kill_rate",
           "ferroptosis_kills", "immune_kills"]


def run_seed(seed: int) -> list:
    """Run the matrix at one seed in an isolated cwd; return its conditions."""
    tmp = Path(tempfile.mkdtemp(prefix=f"simtme-seed{seed}-"))
    try:
        env = dict(os.environ)
        env["FERRO_SEED"] = str(seed)
        proc = subprocess.run([str(BINARY)], cwd=tmp, env=env,
                              capture_output=True, text=True, timeout=7200)
        summary = tmp / "output" / "tme" / "tme_summary.json"
        if proc.returncode != 0 or not summary.exists():
            print(f"  seed {seed}: FAILED rc={proc.returncode} {proc.stderr[-200:]}",
                  file=sys.stderr)
            return []
        return json.loads(summary.read_text())["conditions"]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def key(c: dict) -> tuple:
    return (c["treatment"], c["o2_condition"], c.get("immune_mode", ""))


def boot_ci(vals, n_boot=4000, lo=2.5, hi=97.5, rng=None):
    """Percentile bootstrap CI for the median."""
    if len(vals) < 2:
        return (vals[0], vals[0]) if vals else (0.0, 0.0)
    rng = rng or random.Random(12345)
    meds = []
    n = len(vals)
    for _ in range(n_boot):
        meds.append(statistics.median(rng.choices(vals, k=n)))
    meds.sort()
    return meds[int(lo / 100 * n_boot)], meds[min(int(hi / 100 * n_boot), n_boot - 1)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--start", type=int, default=42, help="first seed (42 = the baseline)")
    args = ap.parse_args()

    if not BINARY.exists():
        raise SystemExit(f"{BINARY} not built; run: cargo build --release -p sim-tme")

    seeds = list(range(args.start, args.start + args.seeds))
    print(f"running sim-tme at {len(seeds)} seeds, {args.jobs} at a time ...", flush=True)
    runs = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for seed, conds in zip(seeds, ex.map(run_seed, seeds)):
            if conds:
                runs[seed] = conds
                print(f"  seed {seed} done ({len(conds)} conditions)", flush=True)

    if len(runs) < 2:
        raise SystemExit("need at least 2 successful runs")

    # gather per-condition, per-metric samples
    samples = {}
    for seed, conds in runs.items():
        for c in conds:
            for m in METRICS:
                v = c.get(m)
                if v is None:
                    continue
                samples.setdefault((key(c), m), []).append(float(v))

    baseline = {key(c): c for c in runs.get(args.start, [])}

    rows = []
    for (k, m), vals in sorted(samples.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if len(vals) < 2:
            continue
        med = statistics.median(vals)
        lo, hi = boot_ci(vals)
        base = baseline.get(k, {}).get(m)
        spread = (max(vals) - min(vals))
        rows.append(dict(treatment=k[0], o2=k[1], immune=k[2], metric=m,
                         n=len(vals), seed42=base, median=med, lo=lo, hi=hi,
                         min=min(vals), max=max(vals),
                         rel_spread=(spread / med if med else 0.0)))

    RAW.write_text(json.dumps({"seeds": sorted(runs), "rows": rows}, indent=1), encoding="utf-8")

    volatile = [r for r in rows if r["rel_spread"] > 0.25 and r["median"] > 0]
    L = [
        "# Seed replication for sim-tme (#SEED-REP)", "",
        "Generated by `scripts/seed_replication.py`.", "",
        f"`sim-tme` ran at **{len(runs)} seeds** ({min(runs)}-{max(runs)}). Seed 42 is the",
        "committed baseline every manuscript number comes from; the rest are replicates.",
        "`FERRO_SEED` unset still reproduces seed 42 byte-identically.", "",
        "## Why this exists", "",
        "The engine had a hard-coded seed and no replicate loop, so every number it",
        "contributes was a single draw reported as a point estimate. Where an interval",
        "below spans a comparison, that comparison is not resolved by one run.", "",
        f"## Metrics whose spread exceeds 25% of the median ({len(volatile)})", "",
    ]
    if volatile:
        L += ["| treatment | O2 | immune | metric | seed 42 | median | 95% CI | min-max |",
              "|---|---|---|---|---|---|---|---|"]
        for r in sorted(volatile, key=lambda r: -r["rel_spread"])[:25]:
            b = f"{r['seed42']:.4g}" if r["seed42"] is not None else "-"
            L.append(f"| {r['treatment']} | {r['o2']} | {r['immune']} | `{r['metric']}` | "
                     f"{b} | {r['median']:.4g} | {r['lo']:.4g}-{r['hi']:.4g} | "
                     f"{r['min']:.4g}-{r['max']:.4g} |")
    else:
        L.append("None. Every metric is stable across seeds at this sample size.")

    L += ["", "## Full table", "",
          "| treatment | O2 | immune | metric | n | seed 42 | median | 95% CI |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        b = f"{r['seed42']:.4g}" if r["seed42"] is not None else "-"
        L.append(f"| {r['treatment']} | {r['o2']} | {r['immune']} | `{r['metric']}` | "
                 f"{r['n']} | {b} | {r['median']:.4g} | {r['lo']:.4g}-{r['hi']:.4g} |")

    L += ["", "## How to read this", "",
          "The interval is a percentile bootstrap on the MEDIAN across seeds. It captures",
          "stochastic run-to-run variation ONLY. It is not a parameter-uncertainty",
          "interval -- those are the prior-predictive intervals in",
          "`analysis/headline-uncertainty-report.md`, and they are much wider. A number",
          "needs both before it can be quoted as a magnitude.", ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT} and {RAW}")
    print(f"metrics with >25% relative spread: {len(volatile)}/{len(rows)}")
    for r in sorted(volatile, key=lambda r: -r["rel_spread"])[:8]:
        print(f"  {r['treatment']:<8} {r['o2']:<16} {r['metric']:<20} "
              f"seed42={r['seed42']:.4g} median={r['median']:.4g} "
              f"range={r['min']:.4g}-{r['max']:.4g}")


if __name__ == "__main__":
    main()
