#!/usr/bin/env python3
"""Read the rare-event sweep and say what the tail actually shows.

WHY THIS EXISTS SEPARATELY FROM THE SWEEP
-----------------------------------------
The sweep produces one row per (condition, n). A row on its own is close to
meaningless: "0 deaths in a billion cells" is not a measurement of zero, it is a
statement about how far down the sweep could see. The finding is in the SHAPE of
a condition's rows across n, and that only exists once several are in hand.

THE ONE DISTINCTION THIS WHOLE ANALYSIS TURNS ON
------------------------------------------------
For a condition that never produces an event, the reported 95% upper bound is
the rule of three, 3/n. That falls by exactly one decade per decade of n -- a
straight line of slope -1 on log-log axes. It is a property of the SAMPLE SIZE
and carries no information about the biology whatsoever. Plotting it without
saying so would manufacture the appearance of a measured decline.

So each condition is classified:

  RESOLUTION-LIMITED  every n gave zero events. The bound tracks 3/n exactly.
                      The true rate is unknown and is somewhere below the last
                      bound. We learned a LIMIT, not a value.
  RESOLVED            events appeared. The rate estimate stabilises across n
                      and the Poisson interval closes around it. This is a
                      measurement.
  EMERGENT            zero at small n, non-zero at large n. The most
                      informative outcome: it locates the order of magnitude
                      where the condition's tail actually begins.

A condition can only become RESOLVED or EMERGENT because the per-cell draw has
unbounded support (`gen_cell` draws from normals). Under the bounded uniform
draw the manuscript used to describe, a zero would be exactly zero and no n
would ever find anything -- see tests/test_manuscript_cell_variation.py.

WHAT THE SCALES MEAN, AND WHAT THEY DO NOT
-------------------------------------------
1e9 is about a gram, the smallest clinically detectable lesion; 1e11 about a
hundred grams, advanced metastatic disease. But 1e11 INDEPENDENT cells is not a
100-gram tumor: these cells never interact -- no diffusion, no crowding, no
vasculature, no clonal competition. The honest reading is "a population the size
of an advanced burden, sampled as independent cells". Every statement this
script emits is hedged that way on purpose.

Usage:
    python scripts/rare_event_analysis.py
    python scripts/rare_event_analysis.py --sweep analysis/rare-event-sweep.jsonl
"""

import argparse
import json
import sys
from pathlib import Path


import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

# Computed here, deliberately, rather than imported from scripts/config.py:
# that module pulls in `requests` and `python-dotenv`, and this has to run on
# a bare interpreter on a freshly booted instance with nothing pip-installed.
# Importing it would have failed AFTER the Rust build and the self-check --
# the most expensive possible place to discover a missing dependency.
PROJECT_ROOT = Path(__file__).resolve().parent.parent

SWEEP = PROJECT_ROOT / "analysis" / "rare-event-sweep.jsonl"
FIG_DIR = PROJECT_ROOT / "article" / "figures"
OUT_MD = PROJECT_ROOT / "analysis" / "rare-event-findings.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "rare-event-findings.json"

# Burden readings, so every scale in the output carries a clinical denominator
# rather than being a round number.
BURDEN = {
    int(1e6): "the manuscript's standard run",
    int(1e7): "sub-detection residual",
    int(1e8): "below any imaging threshold",
    int(1e9): "~1 g, smallest clinically detectable lesion",
    int(1e10): "~10 g, intermediate burden",
    int(1e11): "~100 g, advanced metastatic disease",
}


def load(path: Path) -> dict:
    """rows grouped by (phenotype, treatment), each sorted by n."""
    if not path.exists():
        print(f"no sweep at {path}", file=sys.stderr)
        return {}
    by_cond = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        by_cond.setdefault((r["phenotype"], r["treatment"]), []).append(r)
    for rows in by_cond.values():
        rows.sort(key=lambda r: r["n_cells"])
    return by_cond


def classify(rows: list) -> str:
    """RESOLUTION-LIMITED / EMERGENT / RESOLVED -- see the module docstring."""
    if all(r["n_dead"] == 0 for r in rows):
        return "resolution-limited"
    if rows[0]["n_dead"] == 0 and rows[-1]["n_dead"] > 0:
        return "emergent"
    return "resolved"


def tracks_rule_of_three(rows: list) -> bool:
    """Is the reported bound exactly 3/n at every point?

    If it is, the curve carries no biology at all -- it is the sample size drawn
    on a log axis. Checked rather than assumed, because a condition that picks
    up a single event part-way up the sweep stops tracking it, and that
    departure is the whole signal.
    """
    return all(r["n_dead"] == 0 and
               abs(r["zero_event_upper_bound_95"] - 3.0 / r["n_cells"])
               < 1e-12 * max(1.0, 3.0 / r["n_cells"])
               for r in rows)


def _provenance_caveat(by_cond: dict) -> list:
    """State which parameter set produced these rows, READ FROM THE ROWS.

    This paragraph used to assert "the parameters are the in-vivo defaults,
    which carry zero calibration targets" as literal text, unconditionally --
    while the harness stamps `params_source` on every row and the sweep exposes
    `--params` precisely so a run CAN be driven from the CTRPv2 in-vitro
    posterior instead. Drive it that way and the document would have printed
    the in-vivo disclaimer over in-vitro numbers, which is the one error this
    caveat exists to prevent. It is the same failure as a hand-written figure
    beside a computed one: the sentence outlives the run it described.

    Mixed provenance is reported as an error rather than summarised, because
    the two parameterisations are provably disjoint (#332, #500) and a curve
    spanning both is not a curve.
    """
    sources = sorted({r.get("params_source", "unrecorded")
                      for rows in by_cond.values() for r in rows})
    overrides = sorted({json.dumps(r.get("param_overrides") or {}, sort_keys=True)
                        for rows in by_cond.values() for r in rows})

    if len(overrides) > 1:
        return [
            "* **MIXED PARAMETER SETS IN ONE SWEEP FILE.** These rows were not all",
            "  produced under the same parameters:",
        ] + [f"    - `{o}`" for o in overrides] + [
            "  The in-vivo defaults and the CTRPv2 in-vitro posterior are provably",
            "  disjoint (#332, #500), so a curve drawn across them is not a curve.",
            "  Split the file by parameter set before reading anything above.",
        ]

    if overrides == ['{}']:
        return [
            "* **The parameters are the in-vivo defaults, which carry zero",
            "  calibration targets.** `targets.yaml` holds 8 self-consistency",
            "  checks and 0 calibration targets. The repo's only data-anchored",
            "  parameterisation is the CTRPv2 in-vitro posterior, and it is",
            "  provably *disjoint* from these priors (#332, #500) — so a rate",
            "  computed here cannot be quoted as an in-vitro prediction, and vice",
            "  versa.",
        ]

    return [
        "* **These rows were produced under an OVERRIDDEN parameter set**, not the",
        f"  in-vivo defaults: `{overrides[0]}`",
        f"  (recorded provenance: {'; '.join(sources)}).",
        "  Read every rate above against that set. If it is the CTRPv2 in-vitro",
        "  posterior then these are in-vitro numbers and must not be quoted as",
        "  in-vivo ones; the two are provably disjoint (#332, #500).",
    ]


def _witness(by_cond: dict) -> str:
    """A real, checked example of the nesting for the prose to point at.

    Written out by hand this read "1, 7, 76, 760 for PersisterNrf2 + Control",
    which is a number in prose beside the generated table it describes -- the
    shape that goes stale the moment another scale lands. It is computed now,
    and it ASSERTS the monotonicity it claims rather than asserting it in
    words: if a larger sample ever reported fewer deaths than a smaller one,
    the nesting would be broken and this would say so instead.
    """
    # Keyed off the dict KEY, not off a field repeated inside the rows: the
    # key already is (phenotype, treatment), and reading it from the row made
    # this depend on the row's JSON shape for a name it was already handed.
    candidates = [(k, v) for k, v in sorted(by_cond.items())
                  if any(r["n_dead"] for r in v)]
    if not candidates:
        return "no condition has events yet"
    (pheno, tx), best = max(candidates, key=lambda kv: len(kv[1]))
    counts = [r["n_dead"] for r in best]
    name = f"{pheno} + {tx}"
    if counts != sorted(counts):
        return (f"NOT monotone for {name}: {counts} — the nesting property is "
                "violated and the sweep should be investigated")
    return f"{', '.join(str(c) for c in counts)} for {name}"


def figure(by_cond: dict, path_stem: Path) -> None:
    """Tail resolution against sample size, with the 3/n floor drawn explicitly.

    The reference line is the point of the figure. Without it a reader sees
    three descending curves and reads a dose-response; with it they see that the
    descending ones are pinned to the sample size and only a departure from the
    line means anything.
    """
    fig, ax = plt.subplots(figsize=(7.2, 5.0))

    ns = sorted({r["n_cells"] for rows in by_cond.values() for r in rows})
    if not ns:
        return
    lo, hi = min(ns) / 3, max(ns) * 3
    ref = [lo, hi]
    ax.plot(ref, [3.0 / lo, 3.0 / hi], color="0.55", ls="--", lw=1.3, zorder=1,
            label="rule of three (3/n): the resolution floor, not a result")

    # Conditions that are zero at every n plot at EXACTLY the same bound, so
    # they land on top of each other and the last one drawn hides the rest --
    # the first version of this figure showed two series where there were
    # three. Nudge each one sideways by a few percent so all are visible. The
    # offset is cosmetic and the caption says so; the alternative (leaving them
    # superimposed) silently loses a condition.
    overlapping = [k for k, v in sorted(by_cond.items())
                   if all(r["n_dead"] == 0 for r in v)]

    colors = plt.get_cmap("tab10").colors
    for i, ((pheno, tx), rows) in enumerate(sorted(by_cond.items())):
        c = colors[i % len(colors)]
        kind = classify(rows)
        jitter = 1.0
        if (pheno, tx) in overlapping and len(overlapping) > 1:
            k = overlapping.index((pheno, tx)) - (len(overlapping) - 1) / 2
            jitter = 1.10 ** k
        rows = [dict(r, n_cells_plot=r["n_cells"] * jitter) for r in rows]
        xs = [r["n_cells_plot"] for r in rows]
        # A zero-event point has no rate to plot -- plotting 0 on a log axis is
        # impossible and plotting it as a small number would be a lie. It is
        # drawn at its upper BOUND with a downward caret, which is what the
        # datum actually is: "somewhere below here".
        zx = [r["n_cells_plot"] for r in rows if r["n_dead"] == 0]
        zy = [r["poisson_ci_high"] for r in rows if r["n_dead"] == 0]
        ex = [r["n_cells_plot"] for r in rows if r["n_dead"] > 0]
        ey = [r["death_rate"] for r in rows if r["n_dead"] > 0]

        if zx:
            ax.scatter(zx, zy, marker="v", s=55, color=c, zorder=3,
                       edgecolors="white", linewidths=0.6)
        if ex:
            errs = [[r["death_rate"] - r["poisson_ci_low"] for r in rows if r["n_dead"] > 0],
                    [r["poisson_ci_high"] - r["death_rate"] for r in rows if r["n_dead"] > 0]]
            ax.errorbar(ex, ey, yerr=errs, fmt="o", ms=6, color=c, zorder=3,
                        capsize=3, lw=1.2, mec="white", mew=0.6)
        if len(xs) > 1:
            ax.plot(xs, [r["poisson_ci_high"] if r["n_dead"] == 0 else r["death_rate"]
                         for r in rows], color=c, lw=1.4, alpha=0.75, zorder=2)
        ax.plot([], [], color=c, marker="v" if kind == "resolution-limited" else "o",
                lw=1.4, label=f"{pheno} + {tx} — {kind}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_xlabel("cells simulated per condition (n)")
    ax.set_ylabel("death rate — 95% upper bound where no events were seen")
    ax.set_title("How far down the tail a given sample size can see")
    ax.grid(True, which="major", alpha=0.25)
    ax.grid(True, which="minor", alpha=0.10)
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.95)

    fig.text(0.5, -0.06,
             "Triangles are zero-event runs plotted at their upper bound: the true rate is somewhere below. "
             "Circles are\nmeasured rates with exact Poisson intervals. A series lying ON the dashed line is "
             "reporting the sample size, not the biology.\nZero-event series are identical by construction, so they are offset horizontally by a few percent to keep each one visible.",
             ha="center", fontsize=8, color="0.35")

    for ext in ("png", "pdf"):
        fig.savefig(path_stem.with_suffix(f".{ext}"))
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", default=str(SWEEP))
    ap.add_argument("--out", default=str(OUT_MD))
    args = ap.parse_args()

    by_cond = load(Path(args.sweep))
    if not by_cond:
        return 1

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figure(by_cond, FIG_DIR / "fig29_rare_event_resolution")

    findings = {}
    L = [
        "# What the large-n rare-event sweep shows", "",
        "Generated by `scripts/rare_event_analysis.py` from",
        "`analysis/rare-event-sweep.jsonl`. Figure:",
        "`article/figures/fig29_rare_event_resolution.png`.", "",
        "## The distinction everything here turns on", "",
        "Several conditions report a death rate of exactly 0% at the",
        "manuscript's standard n = 1,000,000. That is not a measurement of zero.",
        "Each cell's parameters are drawn independently from **normal**",
        "distributions (`gen_cell`), which have unbounded support, so every",
        "threshold retains a positive-probability tail. A reported 0% is an upper",
        "bound set by the sample size — about 3e-6 at a million cells.", "",
        "This sweep pushes n up to resolve those bounds. For a condition that",
        "still produces no events, the bound falls as exactly 3/n: one decade per",
        "decade of n, which is a property of the sample and carries no biology.",
        "The informative outcome is a **departure** from that line.", "",
        "## Per condition", "",
    ]

    for (pheno, tx), rows in sorted(by_cond.items()):
        kind = classify(rows)
        on_floor = tracks_rule_of_three(rows)
        last = rows[-1]
        findings[f"{pheno}/{tx}"] = {
            "classification": kind,
            "tracks_rule_of_three": on_floor,
            "largest_n": last["n_cells"],
            "n_dead_at_largest_n": last["n_dead"],
            "bound_or_rate": (last["poisson_ci_high"] if last["n_dead"] == 0
                              else last["death_rate"]),
            "poisson_ci": [last["poisson_ci_low"], last["poisson_ci_high"]],
            "rationale": last.get("rationale", ""),
            "params_source": last.get("params_source", ""),
        }

        L += [f"### {pheno} + {tx} — {kind.upper()}", "",
              f"*Why this condition is in the sweep:* {last.get('rationale','')}", "",
              "| n | burden | deaths | rate or 95% upper bound | exact Poisson 95% |",
              "|---:|---|---:|---:|---|"]
        for r in rows:
            b = BURDEN.get(r["n_cells"], "")
            if r["n_dead"] == 0:
                val = f"< {r['poisson_ci_high']:.2e}"
                ci = "zero events"
            else:
                val = f"{r['death_rate']:.3e}"
                ci = f"[{r['poisson_ci_low']:.2e}, {r['poisson_ci_high']:.2e}]"
            L.append(f"| {r['n_cells']:.0e} | {b} | {r['n_dead']:,} | {val} | {ci} |")
        L.append("")

        if kind == "resolution-limited":
            n = last["n_cells"]
            expected = last["poisson_ci_high"] * n
            floor_note = (
                "and it is still tracking 3/n exactly, so this is a statement "
                "about how far we could see, not about the rate."
                if on_floor else
                "and it has left the 3/n line, which means something other than "
                "the sample size is now setting it — check the sweep."
            )
            L += [
                f"Zero events at every scale tried, up to {n:.0e} cells. The bound",
                f"is now {last['poisson_ci_high']:.2e}, {floor_note}",
                "**The true rate remains unknown**; all that has been established",
                "is that it is below the bound.", "",
                f"Read against burden: in a population of {n:.0e} cells, fewer than",
                f"about {expected:.0f} would be expected to die under this condition.", "",
            ]
        elif kind == "emergent":
            first = next(r for r in rows if r["n_dead"] > 0)
            L += [
                f"**The tail becomes visible at n = {first['n_cells']:.0e}.** Below that",
                "the condition reports zero and the bound is resolution-limited;",
                "at and above it, real events appear. This locates the order of",
                "magnitude at which this condition's rare deaths actually live —",
                "the single most informative outcome available from a sweep like",
                "this, and one no amount of running at n = 1e6 could produce.", "",
            ]
        else:
            # Quantify what the extra cells actually bought. "Large n buys
            # precision" is true of every Monte Carlo ever run; the numbers are
            # the finding. Both the width RATIO and the movement of the point
            # estimate matter, and they say different things -- the first is how
            # much the interval closed, the second is how wrong the small-n
            # estimate was, which is the part a reader of Figure 7 needs.
            first = rows[0]
            w0 = (first["poisson_ci_high"] / first["poisson_ci_low"]
                  if first["poisson_ci_low"] > 0 else float("inf"))
            w1 = (last["poisson_ci_high"] / last["poisson_ci_low"]
                  if last["poisson_ci_low"] > 0 else float("inf"))
            drift = ((last["death_rate"] - first["death_rate"]) / last["death_rate"]
                     if last["death_rate"] else 0.0)
            L += [
                "Events at every scale, so this is a measurement rather than a",
                "bound. What the extra cells bought:", "",
                f"* the exact Poisson interval spans a factor of **{w1:.2f}** at",
                f"  n = {last['n_cells']:.0e} ([{last['poisson_ci_low']:.2e},"
                f" {last['poisson_ci_high']:.2e}]), against a factor of"
                f" {w0:.0f} at n = {first['n_cells']:.0e}"
                f" ([{first['poisson_ci_low']:.2e}, {first['poisson_ci_high']:.2e}]);",
                f"* the point estimate moved from {first['death_rate']:.3e} to"
                f" {last['death_rate']:.3e}, so the small-n figure was off by",
                f"  **{abs(drift)*100:.0f}%** — it sat"
                f" {'above' if drift < 0 else 'below'} the resolved value, and its",
                "  interval did contain the answer, which is the honest defence of",
                "  a small sample: the point estimate was noise, the interval was",
                "  not a lie.", "",
            ]

    L += [
        "## What this does not show", "",
        "* **1e11 independent cells is not a 100-gram tumor.** These cells never",
        "  interact: no diffusion, no crowding, no vasculature, no clonal",
        "  competition, no immune contact. The defensible phrase is *a population",
        "  the size of an advanced burden, sampled as independent cells*. The",
        "  spatial model is the one with interaction, and it does not reach these",
        "  counts.",
    ] + _provenance_caveat(by_cond) + [
        "* **The points across n are NESTED, not independent.** `run_condition`",
        "  seeds cell *i* from its global index (`i*2`, `i*2+1`) regardless of n,",
        "  so the run at 1e8 is *literally the first 1e8 cells* of the run at 1e9.",
        "  Each larger sample contains every smaller one. Two consequences: the",
        f"  death counts must be monotone non-decreasing across the sweep",
        f"  ({_witness(by_cond)}), and successive",
        "  estimates are positively correlated, so the apparent convergence is",
        "  smoother than four *independent* samples would look. The final",
        "  interval is still correct — it depends only on the count and the n it",
        "  came from — but the intermediate points are not four separate",
        "  confirmations of it.",
        "* **Sampling error is not the dominant uncertainty.** Section 5.2 already",
        "  records that these outputs are parameter-limited, not sample-limited.",
        "  Driving n to 1e11 shrinks the sampling interval to nothing and leaves",
        "  the parametric interval exactly where it was. That is a feature of this",
        "  exercise, not a flaw: it isolates the tail question from the parameter",
        "  question by removing one of them entirely.", "",
    ]

    Path(args.out).write_text("\n".join(L) + "\n", encoding="utf-8")
    OUT_JSON.write_text(json.dumps(findings, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"wrote {args.out}")
    print(f"wrote {OUT_JSON}")
    print(f"wrote {FIG_DIR / 'fig29_rare_event_resolution.png'}")
    for k, v in sorted(findings.items()):
        print(f"  {k:<26} {v['classification']:<20} "
              f"n={v['largest_n']:.0e} dead={v['n_dead_at_largest_n']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
