#!/usr/bin/env python3
"""Sweep the rare-event tail of the ferroptosis kill switch across burden scales.

WHY
---
Several conditions in Figure 7 report a death rate of exactly 0% at n = 1e6.
Since the per-cell parameters are drawn from unbounded Gaussians, that 0% is an
upper bound set by the sample size (about 3e-6), not a measurement. This sweeps
n upward and records how the bound moves.

The result is the SHAPE, not any single row. A rate that keeps falling in step
with n is still resolution-limited and the true value is unknown; a rate that
settles is a real estimate; and a condition that produces its first events
part-way up the sweep tells you where its tail actually lives.

WHY THESE SCALES
----------------
They are tumor burdens, chosen so each row has a clinical denominator:

    1e9    ~1 g, the smallest clinically detectable lesion
    1e10   an intermediate burden
    1e11   ~100 g, advanced metastatic disease

1e12 upward is typically lethal and is deliberately excluded. At a survival rate
of 1e-9, a 1e11-cell burden still harbours about a hundred survivors, which is
the resistance-escape question with a real denominator under it.

WHAT THIS IS NOT
----------------
1e11 INDEPENDENT cells is not a 100-gram tumor. These cells never interact: no
diffusion, no crowding, no vasculature, no clonal competition. The honest phrase
is "a population the size of an advanced burden, sampled as independent cells".

INCREMENTAL BY DESIGN
---------------------
Each (condition, scale) result is appended to the output file as soon as it
finishes, and a completed cell is never recomputed on a re-run. A sweep that
takes hours must survive being interrupted, and on a metered instance with an
auto-terminate timer that is not optional.

Usage:
    python scripts/rare_event_sweep.py --scales 1e9
    python scripts/rare_event_sweep.py --scales 1e9,1e10 --out analysis/...jsonl
"""

import argparse
import json
import signal
import subprocess
import sys
import time
from pathlib import Path



# Computed here, deliberately, rather than imported from scripts/config.py:
# that module pulls in `requests` and `python-dotenv`, and this has to run on
# a bare interpreter on a freshly booted instance with nothing pip-installed.
# Importing it would have failed AFTER the Rust build and the self-check --
# the most expensive possible place to discover a missing dependency.
PROJECT_ROOT = Path(__file__).resolve().parent.parent

BIN = PROJECT_ROOT / "simulations" / "target" / "release" / "sim-scale"
OUT = PROJECT_ROOT / "analysis" / "rare-event-sweep.jsonl"

# The conditions worth spending cells on, and why each is here. Chosen from a
# scan at n=1e6: these are the ones whose reported rate is at or below the
# resolution limit, so they are the only ones a larger n can inform.
CONDITIONS = [
    ("Glycolytic", "RSL3",
     "the manuscript's own selectivity constraint -- RSL3 is supposed to spare "
     "this phenotype. 0 deaths at 1e6 AND at 1e8, so the bound is already 3e-8."),
    ("Glycolytic", "Control",
     "baseline viability, the constraint that all phenotypes stay under 2% death "
     "untreated. 0 deaths at 1e6; the true basal rate is unmeasured."),
    ("PersisterNrf2", "Control",
     "exactly 1 death in 1e6 -- sitting on the resolution limit, where a single "
     "event carries a Poisson interval spanning two orders of magnitude."),
]


def done_cells(path: Path) -> set:
    """(phenotype, treatment, n) already recorded, so a re-run resumes."""
    if not path.exists():
        return set()
    out = set()
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        try:
            r = json.loads(line)
        except ValueError:
            continue
        out.add((r["phenotype"], r["treatment"], int(r["n_cells"])))
    return out


def poisson_ci(k: int, n: int):
    """Exact Poisson 95% interval on a rate, which is what small counts need.

    A Wilson interval is fine at 42%, but at 1 event in 1e6 the binomial normal
    approximation is meaningless. Uses the chi-square relationship, computed
    from the gamma quantile via a simple bisection so this stays stdlib-only.
    """
    import math

    def gamma_q(p, shape):
        if shape <= 0:
            return 0.0
        lo, hi = 0.0, max(10.0, shape * 10)
        for _ in range(200):
            mid = (lo + hi) / 2
            # regularised lower incomplete gamma by series/continued fraction
            s, term, k_ = 0.0, 1.0 / shape, 0
            x = mid
            if x < shape + 1:
                term = 1.0 / shape
                s = term
                while k_ < 500:
                    k_ += 1
                    term *= x / (shape + k_)
                    s += term
                    if abs(term) < abs(s) * 1e-15:
                        break
                cdf = s * math.exp(-x + shape * math.log(max(x, 1e-300)) - math.lgamma(shape))
            else:
                cdf = 1.0
                a, b, c, d = 1.0, x + 1 - shape, 1e300, 1.0 / (x + 1 - shape)
                h = d
                for i in range(1, 500):
                    an = -i * (i - shape)
                    b += 2
                    d = an * d + b
                    if abs(d) < 1e-300:
                        d = 1e-300
                    c = b + an / c
                    if abs(c) < 1e-300:
                        c = 1e-300
                    d = 1.0 / d
                    delta = d * c
                    h *= delta
                    if abs(delta - 1) < 1e-15:
                        break
                cdf = 1 - math.exp(-x + shape * math.log(max(x, 1e-300)) - math.lgamma(shape)) * h
            if cdf < p:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    lo = 0.0 if k == 0 else gamma_q(0.025, k) / n
    hi = gamma_q(0.975, k + 1) / n
    return lo, hi


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scales", default="1e9",
                    help="comma-separated cell counts, e.g. 1e9,1e10,1e11")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--params", default=None,
                    help="JSON override map; omitted means Params::default()")
    args = ap.parse_args()

    if not BIN.exists():
        print(f"build the harness first: cargo build --release -p sim-scale\n  ({BIN})",
              file=sys.stderr)
        return 1
    # The gate: refuse to spend hours on an engine that cannot reproduce a
    # known result.
    v = subprocess.run([str(BIN), "--verify"], capture_output=True, text=True)
    if v.returncode != 0:
        print(v.stdout + v.stderr, file=sys.stderr)
        print("self-check failed; not starting the sweep", file=sys.stderr)
        return 1
    print("self-check passed; starting sweep", flush=True)

    out_p = Path(args.out)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    already = done_cells(out_p)
    failures = []
    scales = [int(float(s)) for s in args.scales.split(",")]

    for n in scales:
        for pheno, tx, why in CONDITIONS:
            if (pheno, tx, n) in already:
                print(f"  skip {pheno}/{tx} @ {n:.0e} (already recorded)", flush=True)
                continue
            cmd = [str(BIN), "--cells", str(n), "--phenotype", pheno,
                   "--treatment", tx, "--label", f"sweep-{n:.0e}"]
            if args.params:
                cmd += ["--params", args.params]
            t0 = time.time()
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                # Report the RETURN CODE, not just stderr. A 1e10 run died after
                # 27 minutes and logged "FAILED ... : " with nothing after the
                # colon, because a process killed by a signal writes no stderr at
                # all -- so the one line that could have identified it was the
                # one not printed. Python reports a signal death as a NEGATIVE
                # returncode, which distinguishes "the OS killed it" (SIGKILL
                # under memory pressure, a stray SIGTERM from a parent shell)
                # from "it panicked" in a way stderr alone never can.
                rc = r.returncode
                how = (f"killed by signal {-rc} ({signal.Signals(-rc).name})"
                       if rc < 0 else f"exited {rc}")
                print(f"  FAILED {pheno}/{tx} @ {n:.0e} after "
                      f"{time.time()-t0:.0f}s: {how}", file=sys.stderr)
                for stream, text in (("stderr", r.stderr), ("stdout", r.stdout)):
                    if text and text.strip():
                        print(f"    {stream}: {text.strip()[:500]}", file=sys.stderr)
                    else:
                        print(f"    {stream}: (empty)", file=sys.stderr)
                # Keep going. Each cell is independent and the file is appended
                # per result, so one dead condition must not cost the ones after
                # it -- especially on a metered instance where the remaining
                # budget is wall-clock on a dead-man timer. Rerunning resumes.
                failures.append(f"{pheno}/{tx} @ {n:.0e} ({how})")
                continue
            try:
                rec = json.loads(r.stdout)
            except ValueError as e:
                # Truncated or interleaved stdout: record it as a failure rather
                # than crashing the sweep mid-way through a paid run.
                print(f"  FAILED {pheno}/{tx} @ {n:.0e}: unparseable output "
                      f"({e}); first 200 chars: {r.stdout[:200]!r}", file=sys.stderr)
                failures.append(f"{pheno}/{tx} @ {n:.0e} (unparseable output)")
                continue
            k = rec["n_dead"]
            lo, hi = poisson_ci(k, n)
            rec["poisson_ci_low"] = lo
            rec["poisson_ci_high"] = hi
            rec["rationale"] = why
            with out_p.open("a") as fh:
                fh.write(json.dumps(rec) + "\n")
            print(f"  {pheno:<14}{tx:<8} n={n:.0e}  dead={k:<8} "
                  f"rate={rec['death_rate']:.3e}  "
                  f"poisson95=[{lo:.2e},{hi:.2e}]  {time.time()-t0:.0f}s", flush=True)

    print(f"wrote {out_p}")
    if failures:
        # Exit non-zero so a caller (or a CI step) knows the sweep is INCOMPLETE,
        # while the rows that did land are still on disk and a rerun skips them.
        print(f"{len(failures)} condition(s) did not complete:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
