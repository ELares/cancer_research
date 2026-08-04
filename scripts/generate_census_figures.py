#!/usr/bin/env python3
"""Figures for the census layer (#CENSUS-FIG).

WHY THIS EXISTS
---------------
The census produced the finding that bounds the widest set of manuscript claims
-- the frozen corpus is a 213-fold non-uniform sample of the cancer literature --
and it had no figure. Every relative prevalence the manuscript computes over the
4,830-article corpus inherits that unevenness, so a reader deciding how much to
trust a cross-mechanism comparison needs to see its shape, not a sentence about
it.

It is a separate script from `generate_figures.py` because the inputs are
different in kind: this reads only committed atlas JSON, so it runs offline in
seconds, while `generate_figures.py` loads the full corpus and several gitignored
simulation outputs.

FIGURES
-------
fig28  Two panels. (A) per-mechanism capture -- the frozen corpus's share of the
       census for that mechanism -- on a log axis, because a linear one collapses
       everything below 10% into the axis. (B) what that unevenness does to the
       manuscript's central corpus claim.

Usage:
    python scripts/generate_census_figures.py
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_landscape import PHARMACOLOGICAL, PHYSICAL  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

FIG_DIR = PROJECT_ROOT / "article" / "figures"
LANDSCAPE = PROJECT_ROOT / "analysis" / "atlas-landscape.json"

plt.rcParams.update({
    "font.size": 11, "font.family": "serif", "axes.titlesize": 13,
    "axes.labelsize": 12, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def fig28_census_capture():
    """The frozen corpus's capture of the census, per mechanism."""
    d = json.loads(LANDSCAPE.read_text())
    # A mechanism with no census articles has no defined capture. Dropping it is
    # the honest move: 0/0 plotted as 0% would read as "never captured". A
    # mechanism present in the census but absent from the corpus is dropped for
    # the same reason -- capture 0 has no position on a log axis, and it would
    # otherwise divide by zero in the spread below.
    rows, dropped = [], []
    for r in d["rows"]:
        if r["mesh_census"] > 0 and r["mesh_frozen"] > 0:
            rows.append({**r, "capture": r["mesh_frozen"] / r["mesh_census"]})
        else:
            dropped.append(r["mechanism"])
    rows.sort(key=lambda r: r["capture"])

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(14, 6.6), gridspec_kw={"width_ratios": [2, 1.15]})

    names = [r["mechanism"] for r in rows]
    caps = [100 * r["capture"] for r in rows]
    colors = ["#b5651d" if r["mechanism"] in PHYSICAL else "#4a6f8a" for r in rows]

    # A dot-and-stem plot, NOT bars. A bar encodes magnitude by length measured
    # from zero, and a log axis has no zero -- matplotlib picks the left limit,
    # so bar lengths would encode that arbitrary choice rather than the data. At
    # the default limit the longest:shortest bar reads about 21x against a true
    # 213x, flattening the finding by an order of magnitude. Position on a log
    # axis is the honest encoding, and the stems are a reading aid only.
    y = range(len(rows))
    ax.hlines(y, min(caps) * 0.55, caps, color=colors, alpha=0.35, lw=1.6)
    ax.scatter(caps, y, color=colors, s=58, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names)
    ax.set_xscale("log")
    ax.set_xlim(min(caps) * 0.5, max(caps) * 2.6)
    ax.set_xlabel("share of the census this corpus captures (%, log scale)")
    ax.set_title("A. The frozen corpus samples the literature unevenly")
    ax.grid(axis="x", alpha=0.18, which="both")
    ax.set_axisbelow(True)

    for i, c in enumerate(caps):
        ax.text(c * 1.13, i, f"{c:.2f}%", va="center", fontsize=8.5, color="#333")

    spread = caps[-1] / caps[0]
    ax.set_ylim(-2.4, len(rows) - 0.4)
    ax.annotate(
        f"{spread:.0f}-fold spread:  {names[0]} {caps[0]:.2f}%   to   "
        f"{names[-1]} {caps[-1]:.1f}%",
        xy=(0.02, -1.55), xycoords=("axes fraction", "data"), ha="left",
        va="center", fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.45", fc="#f4efe6", ec="#b5651d", lw=1.1))
    ax.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, fc="#b5651d"),
        plt.Rectangle((0, 0), 1, 1, fc="#4a6f8a")],
        labels=["physical modality", "other mechanism"],
        loc="lower right", fontsize=9, framealpha=0.95)

    # Panel B: all THREE arms of the committed design. Showing only the first and
    # last attributes the whole move to corpus selection, and it does not
    # reconcile -- 9.1 x 3.3 is 30, not 17.6. The middle arm is where the missing
    # factor lives, and it runs the OTHER way: holding the corpus fixed and
    # switching keyword labels for MeSH CUTS the ratio, so by the more
    # independent labelling method the manuscript overstates. Recomputed here
    # from the committed rows rather than quoted, so it cannot drift from panel A.
    R = {r["mechanism"]: r for r in rows}
    tot = lambda ks, c: sum(R[k][c] for k in ks if k in R and R[k].get(c))  # noqa: E731
    arms = [
        ("this corpus,\nkeyword labels\n(what the\nmanuscript reports)",
         tot(PHARMACOLOGICAL, "keyword_frozen") / tot(PHYSICAL, "keyword_frozen"), "#4a6f8a"),
        ("same articles,\nMeSH labels\n(method effect\nalone)",
         tot(PHARMACOLOGICAL, "mesh_frozen") / tot(PHYSICAL, "mesh_frozen"), "#9aaebd"),
        ("the census,\nMeSH labels",
         tot(PHARMACOLOGICAL, "mesh_census") / tot(PHYSICAL, "mesh_census"), "#2f4f4f"),
    ]
    cap_phys = 100 * tot(PHYSICAL, "mesh_frozen") / tot(PHYSICAL, "mesh_census")
    cap_pharm = 100 * tot(PHARMACOLOGICAL, "mesh_frozen") / tot(PHARMACOLOGICAL, "mesh_census")

    vals = [a[1] for a in arms]
    bx.bar(range(3), vals, color=[a[2] for a in arms], width=0.6)
    bx.set_xticks(range(3))
    bx.set_xticklabels([a[0] for a in arms], fontsize=8.5)
    bx.set_ylabel("pharmacological : physical articles")
    bx.set_title("B. Method and corpus pull in opposite directions")
    for x, v in enumerate(vals):
        bx.text(x, v + max(vals) * 0.03, f"{v:.1f} : 1", ha="center",
                fontweight="bold", fontsize=11.5)
    bx.set_ylim(0, max(vals) * 1.62)
    bx.annotate(
        f"Switching keyword labels for MeSH on the SAME\n"
        f"articles CUTS the ratio ({vals[1]/vals[0]:.2f}x). Widening to the\n"
        f"census then RAISES it {vals[2]/vals[1]:.1f}x, because the corpus\n"
        f"over-samples physical modalities {cap_phys/cap_pharm:.1f}x\n"
        f"({cap_phys:.2f}% vs {cap_pharm:.2f}% capture).\n"
        f"Net {vals[2]/vals[0]:.1f}x: the claim survives, but the\n"
        f"corpus effect is larger than the net move shows.",
        xy=(0.5, 0.875), xycoords="axes fraction", ha="center", va="center",
        fontsize=8.2,
        bbox=dict(boxstyle="round,pad=0.4", fc="#eef2f5", ec="#4a6f8a", lw=1))

    fig.suptitle(
        "Sampling bias in the frozen corpus, measured against the full cancer census "
        f"({d['census_records']:,} articles)", fontsize=13.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"fig28_census_capture.{ext}")
    plt.close(fig)
    print(f"fig28_census_capture: {len(rows)} mechanisms, {spread:.0f}x spread, "
          f"{vals[0]:.1f}:1 -> {vals[1]:.1f}:1 -> {vals[2]:.1f}:1"
          + (f"; dropped: {', '.join(dropped)}" if dropped else ""))
    return {"spread": spread, "ratios": vals, "dropped": dropped,
            "capture": {r["mechanism"]: r["capture"] for r in rows},
            # Reported so a guard can hold the encoding: a linear axis collapses
            # everything under 10% into the axis and destroys the finding.
            "xscale": ax.get_xscale()}


def main() -> int:
    if not LANDSCAPE.exists():
        print(f"missing {LANDSCAPE}; run scripts/atlas_landscape.py", file=sys.stderr)
        return 1
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig28_census_capture()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
