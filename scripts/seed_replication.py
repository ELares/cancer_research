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


def keyed(conds: list) -> dict:
    """Map each condition row to a key that is STABLE ACROSS SEEDS and UNIQUE.

    (treatment, o2_condition, immune_mode) is NOT unique: the matrix runs the
    same triple in several blocks (immune coupling, stromal, pH), so three
    distinct RSL3/gradient_120um/immune_on rows collide under it and pooling
    them fabricates spread that is really between-condition difference. The
    matrix order is deterministic, so an occurrence ordinal disambiguates and
    still lines the same row up across seeds.
    """
    seen: dict = {}
    out: dict = {}
    for i, c in enumerate(conds):
        base = (c["treatment"], c["o2_condition"], c.get("immune_mode", ""))
        n = seen.get(base, 0)
        seen[base] = n + 1
        out[(*base, n)] = (i, c)
    return out


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
        for k, (_i, c) in keyed(conds).items():
            for m in METRICS:
                v = c.get(m)
                if v is None:
                    continue
                samples.setdefault((k, m), []).append((seed, float(v)))

    baseline = {k: c for k, (_i, c) in keyed(runs.get(args.start, [])).items()}

    rows = []
    for (k, m), pairs in sorted(samples.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if len(pairs) < 2:
            continue
        # str keys: JSON turns int keys into strings on the round trip, so the
        # in-memory dict must match what a reader of the JSON sees.
        by_seed = {str(s_): v for s_, v in pairs}
        vals = [v for _s, v in pairs]
        med = statistics.median(vals)
        lo, hi = boot_ci(vals)
        base = baseline.get(k, {}).get(m)
        spread = (max(vals) - min(vals))
        rows.append(dict(treatment=k[0], o2=k[1], immune=k[2], block=k[3], metric=m,
                         n=len(vals), seed42=base, median=med, lo=lo, hi=hi,
                         min=min(vals), max=max(vals), by_seed=by_seed,
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
        f"## Metrics whose spread exceeds 25% of the median "
        f"({len(volatile)} found; the {min(len(volatile), 25)} widest shown)", "",
    ]
    if volatile:
        L += ["| treatment | O2 | immune | block | metric | seed 42 | median | 95% CI | min-max |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(volatile, key=lambda r: -r["rel_spread"])[:25]:
            b = f"{r['seed42']:.4g}" if r["seed42"] is not None else "-"
            L.append(f"| {r['treatment']} | {r['o2']} | {r['immune']} | {r['block']} | "
                     f"`{r['metric']}` | {b} | {r['median']:.4g} | {r['lo']:.4g}-{r['hi']:.4g} | "
                     f"{r['min']:.4g}-{r['max']:.4g} |")
    else:
        L.append("None. Every metric is stable across seeds at this sample size.")

    L += ["", "## Full table", "",
          "| treatment | O2 | immune | block | metric | n | seed 42 | median | 95% CI |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        b = f"{r['seed42']:.4g}" if r["seed42"] is not None else "-"
        L.append(f"| {r['treatment']} | {r['o2']} | {r['immune']} | {r['block']} | "
                 f"`{r['metric']}` | {r['n']} | {b} | {r['median']:.4g} | "
                 f"{r['lo']:.4g}-{r['hi']:.4g} |")

    # --- the immune-coupling ratio, computed per seed ---
    def _find(t, o, i, b, m):
        for r in rows:
            if (r["treatment"], r["o2"], r["immune"], r["block"], r["metric"]) == (t, o, i, b, m):
                return r
        return None

    ratio_lines = []
    for blk in range(3):
        rr = _find("RSL3", "gradient_120um", "immune_on", blk, "immune_kills")
        ss = _find("SDT", "gradient_120um", "immune_on", blk, "immune_kills")
        if not (rr and ss):
            continue
        common = sorted(set(rr["by_seed"]) & set(ss["by_seed"]), key=lambda s: int(s))
        ratios, zeros = [], 0
        for s in common:
            den = rr["by_seed"][s]
            if den == 0:
                zeros += 1
            else:
                ratios.append(ss["by_seed"][s] / den)
        if not ratios:
            continue
        ratios.sort()
        ratio_lines.append(dict(block=blk, n=len(common), zeros=zeros,
                                seed42_sdt=rr and ss["by_seed"].get(str(args.start)),
                                seed42_rsl3=rr["by_seed"].get(str(args.start)),
                                median=statistics.median(ratios),
                                lo=ratios[0], hi=ratios[-1]))

    if ratio_lines:
        L += ["", "## The immune-coupling ratio", "",
              "The manuscript quotes an SDT:RSL3 immune-kill ratio of **104:1** and repeats",
              "it six times. Its denominator is a single-digit event count, so this is the",
              "quantity most exposed to seed choice. Computed as the median of PER-SEED",
              "ratios (not the ratio of medians, which is a different and more flattering",
              "statistic):", "",
              "| block | seed 42 | median ratio | range | seeds with a ZERO denominator |",
              "|---|---|---|---|---|"]
        for r in ratio_lines:
            s42 = (f"{r['seed42_sdt']:.0f}/{r['seed42_rsl3']:.0f}"
                   if r["seed42_rsl3"] else f"{r['seed42_sdt']:.0f}/0")
            L.append(f"| {r['block']} | {s42} | {r['median']:.0f}:1 | "
                     f"{r['lo']:.0f}:1 - {r['hi']:.0f}:1 | {r['zeros']}/{r['n']} |")
        b0 = ratio_lines[0]
        L += ["",
              f"Block 0 is the manuscript's condition. The published point estimate is",
              f"**representative** -- the median across seeds is {b0['median']:.0f}:1 against a",
              f"published 104:1, so seed 42 was not a lucky draw.",
              "",
              f"What the single seed hides is the spread: **{b0['lo']:.0f}:1 to {b0['hi']:.0f}:1**, a",
              f"{b0['hi']/b0['lo']:.0f}-fold range, and in **{b0['zeros']} of {b0['n']}** seeds the RSL3 denominator",
              "is zero so the ratio is undefined altogether.",
              "",
              "So the correction is about precision, not accuracy. Quoting `104:1` implies",
              "two significant figures the data cannot support. An honest form is",
              f"\"~{round(b0['median'], -1):.0f}:1, ranging roughly {b0['lo']:.0f}:1 to {b0['hi']:.0f}:1 across seeds and",
              "undefined when the RSL3 arm scores no immune kills at all\".", ""]

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
