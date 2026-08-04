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
    # the honest move: 0/0 plotted as 0% would read as "never captured".
    rows = [r for r in d["rows"] if r["mesh_census"] > 0]
    dropped = [r["mechanism"] for r in d["rows"] if r["mesh_census"] == 0]
    for r in rows:
        r["capture"] = r["mesh_frozen"] / r["mesh_census"]
    rows.sort(key=lambda r: r["capture"])

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(13.5, 6.4), gridspec_kw={"width_ratios": [2.15, 1]})

    names = [r["mechanism"] for r in rows]
    caps = [100 * r["capture"] for r in rows]
    colors = ["#b5651d" if r["mechanism"] in PHYSICAL else "#4a6f8a" for r in rows]
    ax.barh(range(len(rows)), caps, color=colors, height=0.72)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(names)
    ax.set_xscale("log")
    ax.set_xlabel("share of the census this corpus captures (%, log scale)")
    ax.set_title("A. The frozen corpus samples the literature unevenly")

    for i, (r, c) in enumerate(zip(rows, caps)):
        ax.text(c * 1.12, i, f"{c:.2f}%", va="center", fontsize=8.5, color="#333")

    # An empty band below the shortest bar, so the caption and legend never sit
    # on top of a bar's own value label.
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

    # Panel B: what the unevenness does to the manuscript's central claim. The
    # ratio is recomputed here from the same committed rows rather than quoted,
    # so it cannot drift away from panel A.
    # PHARMACOLOGICAL is a curated list of drug modalities, imported rather than
    # restated so this panel cannot drift from the committed analysis. It is NOT
    # "everything that is not physical": counting the delivery platforms and
    # genetic tools too would report 12.5:1 instead of 9.1:1 by including things
    # that are neither a drug class nor a physical modality.
    phys = [r for r in rows if r["mechanism"] in PHYSICAL]
    pharm = [r for r in rows if r["mechanism"] in PHARMACOLOGICAL]
    kw_ratio = sum(r["keyword_frozen"] for r in pharm) / sum(
        r["keyword_frozen"] for r in phys)
    census_ratio = sum(r["mesh_census"] for r in pharm) / sum(
        r["mesh_census"] for r in phys)
    cap_phys = 100 * sum(r["mesh_frozen"] for r in phys) / sum(
        r["mesh_census"] for r in phys)
    cap_pharm = 100 * sum(r["mesh_frozen"] for r in pharm) / sum(
        r["mesh_census"] for r in pharm)

    bx.bar([0, 1], [kw_ratio, census_ratio], color=["#4a6f8a", "#2f4f4f"], width=0.56)
    bx.set_xticks([0, 1])
    bx.set_xticklabels(["this corpus,\nkeyword labels\n(what the\nmanuscript reports)",
                        "the census,\nMeSH labels"], fontsize=9.5)
    bx.set_ylabel("pharmacological : physical articles")
    bx.set_title("B. The claim survives, and was understated")
    for x, v in [(0, kw_ratio), (1, census_ratio)]:
        bx.text(x, v + census_ratio * 0.03, f"{v:.1f} : 1", ha="center",
                fontweight="bold", fontsize=12)
    bx.set_ylim(0, census_ratio * 1.62)
    bx.annotate(
        f"The corpus over-samples physical\nmodalities {cap_phys/cap_pharm:.1f}x "
        f"({cap_phys:.2f}% vs {cap_pharm:.2f}%\ncapture), so it understates the\n"
        f"imbalance it set out to measure.",
        xy=(0.5, 0.885), xycoords="axes fraction", ha="center", va="center",
        fontsize=8.8,
        bbox=dict(boxstyle="round,pad=0.4", fc="#eef2f5", ec="#4a6f8a", lw=1))

    fig.suptitle(
        "Sampling bias in the frozen corpus, measured against the full cancer census "
        f"({d['census_records']:,} articles)", fontsize=13.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"fig28_census_capture.{ext}")
    plt.close(fig)
    print(f"fig28_census_capture: {len(rows)} mechanisms, {spread:.0f}x spread, "
          f"{kw_ratio:.1f}:1 -> {census_ratio:.1f}:1"
          + (f"; dropped (no census articles): {', '.join(dropped)}" if dropped else ""))


def main() -> int:
    if not LANDSCAPE.exists():
        print(f"missing {LANDSCAPE}; run scripts/atlas_landscape.py", file=sys.stderr)
        return 1
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig28_census_capture()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
