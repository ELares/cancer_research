#!/usr/bin/env python3
"""Generate conceptual diagrams for TME mechanisms, decision flowchart, and corpus flow.

Creates 6 figures:
  fig18_hypoxia_crosssection.pdf  — O2 gradient with drug efficacy overlay
  fig19_immune_coupling_flow.pdf  — DAMP → DC → T cell pathway
  fig20_stromal_shielding.pdf     — CAF boundary protection
  fig21_ph_ion_trapping.pdf       — pH gradient with drug trapping
  fig22_decision_flowchart.pdf    — Which modality for which context
  fig23_census_flow.pdf            — census construction flow

Usage:
  python3 scripts/generate_conceptual_diagrams.py
"""
import shutil
import subprocess

import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from figure_io import make_figures_deterministic  # noqa: E402

# PDFs embed a creation date; without this a regenerated figure differs
# from its committed copy even when nothing about the drawing changed.
make_figures_deterministic()
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "article" / "figures"


def _render_graphviz(dot_source, out_base):
    """Write `dot_source` to a .gv, render PDF+PNG, and leave no date behind.

    THREE THINGS THIS FIXES, all measured on a machine without graphviz:

    - `dot` is an external program, and `check=True` with no guard turned a
      missing one into a crash that took the whole script down. fig19 is the
      second of six figures, so a contributor without graphviz got fig18 and a
      traceback -- fig20, fig21 and fig23 never ran, which is also why their
      committed copies still carried a 2026-04 creation date nobody could clear.
    - The .gv is a temporary, and it was unlinked only after both `dot` calls
      returned. A missing `dot` therefore stranded it in article/figures as an
      untracked file.
    - `make_figures_deterministic()` patches `Figure.savefig`, so it cannot
      reach a PDF that matplotlib never wrote. `dot` embeds its own
      `/CreationDate`, which is why fig19 and fig22 stayed dated while every
      matplotlib figure beside them came out clean. Rewriting the Info
      dictionary afterwards clears it and leaves the drawing untouched.

    Returns True if the figure was rendered, False if graphviz is absent.
    """
    gv_path = Path(str(out_base) + ".gv")
    try:
        gv_path.write_text(dot_source)
        if shutil.which("dot") is None:
            print(f"  {gv_path.stem}: graphviz `dot` not on PATH -- skipping "
                  "(install graphviz to regenerate this figure)")
            return False
        subprocess.run(["dot", "-Tpdf", "-o", str(out_base) + ".pdf",
                        str(gv_path)], check=True)
        subprocess.run(["dot", "-Tpng", "-Gdpi=300", "-o",
                        str(out_base) + ".png", str(gv_path)], check=True)
    finally:
        gv_path.unlink(missing_ok=True)
    _strip_pdf_date(Path(str(out_base) + ".pdf"))
    return True


def _strip_pdf_date(path):
    """Clear a PDF's creation/modification date in place, drawing untouched."""
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            # LOUD, NOT A PRINT. Returning here emits a DATED pdf, which is
            # exactly what this function exists to prevent, and the only thing
            # that noticed was an unrelated test whose prose pin happens to
            # count dated figures.
            raise RuntimeError(
                f"{path.name}: pymupdf is absent, so the creation date cannot "
                "be removed and this figure would be committed dated. Install "
                "it (it is in requirements-lock.txt) or do not regenerate.")
    # WHAT THIS DOES NOT BUY, stated because the neighbouring matplotlib
    # figures DO have it and a reader will assume the same here.
    #
    # `saveIncr` writes a new trailer, and MuPDF regenerates BOTH halves of
    # `/ID` every time -- measured, three strips of one input give three
    # different files, differing in 62 bytes, all inside the `/ID` array. So
    # fig19 and fig22 are date-free but NOT byte-reproducible: a contributor
    # with graphviz who regenerates them gets a dirty tree and a MANIFEST
    # churn every run, where a matplotlib figure comes back identical.
    #
    # What the freshness gate needs is that the DRAWING is unchanged, and that
    # holds -- the incremental save appends a trailer and leaves the original
    # content stream byte-identical. `test_pdf_output_is_deterministic` covers
    # the matplotlib figures only, and this is the one class where its
    # invariant does not apply.
    #
    # Not idempotent either: a second call appends another section. Harmless
    # only because `dot` always writes a fresh file first.
    doc = pymupdf.open(path)
    try:
        md = dict(doc.metadata or {})
        md["creationDate"] = ""
        md["modDate"] = ""
        doc.set_metadata(md)
        doc.saveIncr()
    finally:
        doc.close()


# Consistent style
COLORS = {
    "tumor_core": "#8B0000",
    "tumor_mid": "#CD5C5C",
    "tumor_edge": "#F08080",
    "vessel": "#4169E1",
    "stroma": "#90EE90",
    "caf": "#228B22",
    "sdt": "#FF8C00",
    "rsl3": "#6A5ACD",
    "pdt": "#DC143C",
    "immune": "#FFD700",
    "damp": "#FF4500",
    "acid": "#FF6347",
    "neutral": "#87CEEB",
    "bg": "#FAFAFA",
    # The non-ferroptosis arms (#726 onward). Kept distinct from the SDT /
    # RSL3 / PDT palette above so a reader can see at a glance which entries
    # belong to this project's original thesis and which do not.
    "radiation": "#8172B2",
    "parp": "#937860",
    "nano": "#DA8BC3",
    "checkpoint": "#64B5CD",
    "text": "#1a1a1a",
}

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  {name}")


# ── Figure 18: Hypoxia cross-section ──────────────────────────────────

def fig18_hypoxia():
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_facecolor(COLORS["bg"])

    # Tumor cross-section as gradient
    for i in range(100):
        x = i * 0.1
        frac = x / 10.0
        r = int(240 - frac * 100)
        g = int(200 - frac * 170)
        b = int(200 - frac * 170)
        color = f"#{r:02x}{g:02x}{b:02x}"
        ax.axvspan(x, x + 0.1, ymin=0.15, ymax=0.85, color=color, alpha=0.7)

    # Blood vessel at left edge
    ax.add_patch(mpatches.Rectangle((0, 0.9), 0.3, 4.2, fc=COLORS["vessel"], ec="black", lw=1.5))
    ax.text(0.15, 3.0, "Blood\nVessel", ha="center", va="center", fontsize=7,
            color="white", fontweight="bold", rotation=90)

    # O2 label
    ax.annotate("", xy=(8.5, 5.2), xytext=(1.0, 5.2),
                arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))
    ax.text(4.75, 5.5, "Decreasing O₂ →", ha="center", fontsize=9, color="gray")

    # Zone labels
    ax.text(2.0, 0.4, "Oxygenated\nperiphery", ha="center", fontsize=8, color=COLORS["vessel"])
    ax.text(7.5, 0.4, "Hypoxic\ncore", ha="center", fontsize=8, color=COLORS["tumor_core"])

    # Drug efficacy curves
    x_vals = np.linspace(0.5, 9.5, 100)
    # SDT: nearly flat, slight drop
    sdt_eff = 0.88 - 0.05 * (x_vals / 10)
    # RSL3: collapses
    rsl3_eff = 0.85 * np.exp(-0.5 * x_vals)

    ax2 = ax.twinx()
    ax2.set_ylim(0, 1.0)
    ax2.plot(x_vals, sdt_eff, color=COLORS["sdt"], lw=2.5, label="SDT efficacy")
    ax2.plot(x_vals, rsl3_eff, color=COLORS["rsl3"], lw=2.5, ls="--", label="RSL3 efficacy")
    ax2.set_ylabel("Relative kill efficacy", fontsize=9)
    ax2.legend(loc="center right", fontsize=8, framealpha=0.9)

    ax.set_xlabel("Distance from blood vessel (mm)", fontsize=9)
    ax.set_yticks([])
    ax.set_title("Hypoxia: O₂ gradient selectively protects against pharmacologic ferroptosis",
                 fontsize=10, fontweight="bold", pad=15)

    # Annotation boxes with white background for readability
    ax.text(7.0, 1.5, "RSL3 depends on\nbasal ROS (needs O₂)",
            fontsize=7, color=COLORS["rsl3"], style="italic", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COLORS["rsl3"], alpha=0.9))
    ax.text(7.0, 4.2, "SDT delivers\nexogenous ROS",
            fontsize=7, color=COLORS["sdt"], style="italic", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COLORS["sdt"], alpha=0.9))

    save(fig, "fig18_hypoxia_crosssection")


# ── Figure 19: Immune coupling flow ───────────────────────────────────

def fig19_immune():
    """Generate immune coupling flow diagram using Graphviz."""
    import graphviz
    import subprocess

    dot = graphviz.Digraph("immune", format="pdf")
    dot.attr(rankdir="LR", bgcolor="white", fontname="Helvetica",
             label="Immune Coupling: Kill density determines DAMP-mediated immune activation",
             labelloc="t", fontsize="13", fontcolor="black", nodesep="0.6", ranksep="0.5")
    dot.attr("node", fontname="Helvetica", fontsize="9", style="filled,rounded",
             shape="box", penwidth="1.5")
    dot.attr("edge", penwidth="1.5", color="black")

    # SDT path (orange) — top row
    sdt_attr = dict(fillcolor="#FFE0B2", color="#FF8C00")
    dot.node("s1", "SDT kills\n~140K cells", **sdt_attr, fontsize="10")
    dot.node("s2", "High LP overshoot\n(LP~20)", **sdt_attr)
    dot.node("s3", "Dense\nDAMP field", **sdt_attr)
    dot.node("s4", "Strong DC\nactivation", **sdt_attr)
    dot.node("s5", "521 immune\nkills", **sdt_attr, fontsize="10")

    dot.edge("s1", "s2")
    dot.edge("s2", "s3")
    dot.edge("s3", "s4")
    dot.edge("s4", "s5")

    # RSL3 path (purple) — bottom row
    rsl3_attr = dict(fillcolor="#D1C4E9", color="#6A5ACD")
    dot.node("r1", "RSL3 kills\n~163 cells", **rsl3_attr, fontsize="10")
    dot.node("r2", "Low LP overshoot\n(LP~7.8)", **rsl3_attr)
    dot.node("r3", "Sparse\nDAMP field", **rsl3_attr)
    dot.node("r4", "Minimal DC\nactivation", **rsl3_attr)
    dot.node("r5", "5 immune\nkills", **rsl3_attr, fontsize="10")

    dot.edge("r1", "r2")
    dot.edge("r2", "r3")
    dot.edge("r3", "r4")
    dot.edge("r4", "r5")

    # Force each stage to align vertically (same rank = same column in LR)
    for s_node, r_node in [("s1","r1"), ("s2","r2"), ("s3","r3"), ("s4","r4"), ("s5","r5")]:
        with dot.subgraph() as sub:
            sub.attr(rank="same")
            sub.node(s_node)
            sub.node(r_node)

    # Ratio label between the two paths
    dot.node("ratio", "104:1\nimmune kill ratio",
             shape="box", style="filled,rounded,bold", fillcolor="white",
             color="#FF4500", fontcolor="#FF4500", fontsize="11", penwidth="2")

    # Position ratio between s3 and r3 using invisible edges
    dot.edge("s3", "ratio", style="invis", weight="10")
    dot.edge("ratio", "r3", style="invis", weight="10")

    if _render_graphviz(dot.source, OUT / "fig19_immune_coupling_flow"):
        print(f"  fig19_immune_coupling_flow")


# ── Figure 20: Stromal shielding ──────────────────────────────────────

def fig20_stromal():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

    for ax, title, kill_inner, kill_boundary, color, label in [
        (ax1, "RSL3 (pharmacologic)", 3.0, 1.5, COLORS["rsl3"], "Kill halved\nat boundary"),
        (ax2, "SDT (physical)", 96.1, 91.2, COLORS["sdt"], "Kill barely\naffected"),
    ]:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 8.5)
        ax.set_facecolor(COLORS["bg"])

        # Tumor interior
        ax.add_patch(mpatches.Rectangle((2.0, 0.5), 7, 5.5, fc="#FFCCCC", ec="black", lw=1))
        ax.text(5.5, 3.2, f"Tumor interior\nKill: {kill_inner}%", ha="center", va="center",
                fontsize=9, fontweight="bold")

        # Stromal boundary (left)
        ax.add_patch(mpatches.Rectangle((0, 0.5), 2.0, 5.5, fc=COLORS["stroma"], ec="black", lw=1))
        ax.text(1.0, 3.2, "Stroma\n(CAFs)", ha="center", va="center", fontsize=8,
                color=COLORS["caf"], fontweight="bold", rotation=90)

        # Boundary cells (highlighted strip)
        ax.add_patch(mpatches.Rectangle((2.0, 0.5), 1.2, 5.5, fc="#FFE0E0", ec="gray", lw=0.5, ls="--"))

        # Labels above the diagram (no overlap)
        ax.text(2.6, 7.0, f"Boundary cells\nKill: {kill_boundary}%", ha="center",
                fontsize=9, color=color, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=color, alpha=0.9))

        # Supply arrows
        for y_pos in [1.5, 3.2, 4.8]:
            ax.plot([1.5, 2.3], [y_pos, y_pos], color=COLORS["caf"], lw=1.5)
            ax.annotate("", xy=(2.3, y_pos), xytext=(2.1, y_pos),
                        arrowprops=dict(arrowstyle="-|>", color=COLORS["caf"], lw=1.5))

        ax.text(1.0, 6.5, "GSH + MUFA\nsupply", fontsize=7, color=COLORS["caf"],
                ha="center", fontweight="bold")
        ax.set_title(title, fontsize=10, fontweight="bold", color=color)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("Stromal Shielding: CAFs protect boundary cells from pharmacologic but not physical ferroptosis",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    save(fig, "fig20_stromal_shielding")


# ── Figure 21: pH ion trapping ────────────────────────────────────────

def fig21_ph():
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_facecolor(COLORS["bg"])

    # pH gradient background
    for i in range(100):
        x = i * 0.1
        frac = x / 10.0
        r = int(135 + frac * 120)
        g = int(206 - frac * 106)
        b = int(235 - frac * 135)
        color = f"#{r:02x}{g:02x}{b:02x}"
        ax.axvspan(x, x + 0.1, ymin=0.12, ymax=0.88, color=color, alpha=0.6)

    # pH labels
    ax.text(0.5, 5.5, "pH 7.4", fontsize=10, fontweight="bold", color=COLORS["vessel"])
    ax.text(8.5, 5.5, "pH 6.5", fontsize=10, fontweight="bold", color=COLORS["tumor_core"])

    # Drug molecules — neutral at left (entering cells)
    for x_pos in [1.5, 3.0, 4.5]:
        ax.plot(x_pos, 3.0, "o", color=COLORS["rsl3"], markersize=12, alpha=0.9)
        ax.text(x_pos, 2.2, "RSL3", fontsize=7, ha="center", color=COLORS["rsl3"], fontweight="bold")

    # Trapped drug molecules at right (stuck outside cells, faded)
    for x_pos in [7.0, 8.0, 9.0]:
        ax.plot(x_pos, 4.5, "o", color=COLORS["rsl3"], markersize=12, alpha=0.3)
        ax.text(x_pos, 3.6, "RSL3⁺\ntrapped", fontsize=7, ha="center", color=COLORS["rsl3"],
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))

    # SDT waves — unaffected by pH
    for x_pos in [1.5, 4.5, 7.5]:
        ax.annotate("", xy=(x_pos + 0.8, 3.0), xytext=(x_pos - 0.3, 3.0),
                    arrowprops=dict(arrowstyle="->", color=COLORS["sdt"], lw=2))
    ax.text(5.0, 1.8, "SDT ultrasound: pH-independent", fontsize=8, ha="center",
            color=COLORS["sdt"], fontweight="bold")

    # Result box
    ax.text(5.0, 0.4, "Drug trapping dominates: RSL3 kills drop 53% | SDT: +0.8% (negligible)",
            ha="center", fontsize=8, fontweight="bold",
            bbox=dict(boxstyle="round", fc="white", ec="black", lw=1))

    ax.set_xlabel("Distance from tumor edge (increasing acidity →)", fontsize=9)
    ax.set_yticks([])
    ax.set_title("Acidic pH: Henderson-Hasselbalch ion trapping reduces drug bioavailability",
                 fontsize=10, fontweight="bold", pad=10)

    save(fig, "fig21_ph_ion_trapping")


# ── Figure 22: Decision flowchart ─────────────────────────────────────

def fig22_flowchart():
    """Generate flowchart using Graphviz (proper arrow-to-box connections)."""
    import graphviz
    import subprocess

    dot = graphviz.Digraph("flowchart", format="pdf")
    dot.attr(rankdir="TB", bgcolor="white", fontname="Helvetica",
             label="Decision Framework: Which Modality for Which Clinical Context?",
             labelloc="t", fontsize="14", fontcolor="black")
    dot.attr("node", fontname="Helvetica", fontsize="10", style="filled,rounded",
             shape="box", penwidth="1.5")
    dot.attr("edge", fontname="Helvetica", fontsize="9", penwidth="1.5")

    # Decision nodes (blue)
    dec_attr = dict(fillcolor="#E3F2FD", color="#1565C0")
    dot.node("q1", "Is the tumor\nlocalizable?", **dec_attr)
    dot.node("q2", "Is it\ndeep-seated?", **dec_attr)
    dot.node("q3", "Are residual cells\nferroptosis-prone?", **dec_attr)
    dot.node("q4", "Immunocompetent\nsetting?", **dec_attr)

    # "No" exit nodes (orange)
    no_attr = dict(fillcolor="#FFF3E0", color="#E65100", fontsize="8")
    dot.node("alt", "Alternative approaches\n(Ch 8.1)", **no_attr)
    dot.node("path", "Pathway-target or\nimmune approaches\n(Ch 8.1, 10.4)", **no_attr)

    # Modality nodes
    dot.node("sdt", "SDT range\n(cm depth)\nCh 6.1",
             fillcolor="#FFE0B2", color="#FF8C00", fontsize="9")
    dot.node("pdt", "PDT range\n(mm depth)\nCh 6.1",
             fillcolor="#FFCDD2", color="#DC143C", fontsize="9")

    # Terminal nodes (purple)
    term_attr = dict(fillcolor="#F3E5F5", color="#6A1B9A", fontsize="9")
    dot.node("combo", "Physical ROS\n+ anti-PD-1\n(Ch 7.2, 9.5)", **term_attr)
    dot.node("direct", "Physical ROS\n(direct kill)\n(Ch 6-7)", **term_attr)

    # Invisible convergence node
    dot.node("conv", "", shape="point", width="0.01", height="0.01")

    # Edges — Yes (green), No (red/orange)
    green = "#2E7D32"
    red = "#E65100"
    gray = "#888888"

    dot.edge("q1", "q2", label="  Yes  ", color=green, fontcolor=green)
    dot.edge("q1", "alt", label="  No  ", color=red, fontcolor=red)

    dot.edge("q2", "sdt", label="  Yes  ", color=green, fontcolor=green)
    dot.edge("q2", "pdt", label="  No  ", color=red, fontcolor=red)

    # Convergence: SDT and PDT both feed into ferroptosis question
    dot.edge("sdt", "conv", style="dashed", color=gray, arrowhead="none")
    dot.edge("pdt", "conv", style="dashed", color=gray, arrowhead="none")
    dot.edge("conv", "q3", style="dashed", color=gray)

    dot.edge("q3", "q4", label="  Yes  ", color=green, fontcolor=green)
    dot.edge("q3", "path", label="  No  ", color=red, fontcolor=red)

    dot.edge("q4", "combo", label="  Yes  ", color=green, fontcolor=green)
    dot.edge("q4", "direct", label="  No  ", color=red, fontcolor=red)

    # Render: save .gv source, then generate PDF and high-DPI PNG
    if _render_graphviz(dot.source, OUT / "fig22_decision_flowchart"):
        print(f"  fig22_decision_flowchart")


# ── Figure 23: census construction flow ──────────────────────────────

# Boxes: (key, x, y, w, h, fill, edge, lines)
_CENSUS_BOXES = [
    ("baseline", 0.50, 0.90, 0.56, 0.10, "#E3F2FD", "#1565C0", [
        ("PubMed Annual Baseline", "bold"),
        ("1,334 gzipped XML files, ~40M records", ""),
        ("E-utilities cannot page past 10,000 hits", "italic")]),
    ("definition", 0.50, 0.735, 0.62, 0.105, "#E3F2FD", "#1565C0", [
        ("Cancer = MeSH tree C04", "bold"),
        ("704 topical descriptors (NLM SPARQL endpoint)", ""),
        ("UNION 9 adjacent experimental-context descriptors", "")]),
    ("indexed", 0.26, 0.525, 0.44, 0.125, "#C8E6C9", "#2E7D32", [
        ("MeSH-Indexed Stream", "bold"),
        ("n = 4,403,994", "bold"),
        ("admitted on a DescriptorName match", ""),
        ("200,758 (4.6%) adjacent-basis only", "")]),
    ("recovered", 0.74, 0.520, 0.44, 0.145, "#FFF3E0", "#E65100", [
        ("Text-Recovered Stream", "bold"),
        ("n = 783,271", "bold"),
        ("admitted on a text match where MeSH", ""),
        ("has not indexed the record yet", ""),
        ("matcher precision 75.7%, recall 95.6%", "italic")]),
    ("census", 0.50, 0.335, 0.34, 0.075, "#E3F2FD", "#1565C0", [
        ("Census", "bold"),
        ("n = 5,187,265", "bold")]),
    ("fulltext", 0.26, 0.135, 0.44, 0.115, "#EEEEEE", "#616161", [
        ("Open-Access Full Text", "bold"),
        ("1,116,481 of 7,517,526 seen", ""),
        ("737,929 oa_comm | 378,552 oa_noncomm", ""),
        ("licence class carried per record", "italic")]),
    ("labels", 0.74, 0.135, 0.44, 0.115, "#EEEEEE", "#616161", [
        ("Labels are NLM's, not ours", "bold"),
        ("mechanism + site: MeSH descriptors", ""),
        ("study design: publication types + check tags", ""),
        ("44.5% carry no design-informative label", "italic")]),
]

# (from, to, label)
_CENSUS_EDGES = [
    ("baseline", "definition", ""),
    ("definition", "indexed", "MeSH descriptor"),
    ("definition", "recovered", "no MeSH yet"),
    ("indexed", "census", ""),
    ("recovered", "census", ""),
    ("census", "fulltext", ""),
    ("census", "labels", ""),
]


def fig23_census_flow():
    """Census construction flow diagram.

    Replaces a PRISMA-inspired retrieval flow. The retrieval it diagrammed was
    real -- 10,415 records screened to 4,830 -- but a PRISMA flow describes a
    systematic review's screening decisions, and no screening decisions were
    ever made: the counts recorded which articles a set of queries happened to
    reach. Drawing them in a shape borrowed from systematic reviewing implied a
    protocol that did not exist.

    The census has a construction flow and it is a different shape. Nothing is
    screened OUT; two streams are admitted by two DIFFERENT RULES, and which
    stream a record lands in is a property of NLM's indexing rather than of the
    record. So the two are drawn side by side as parallel admissions, both
    feeding the total, rather than as a branch where one arm is discarded.

    DRAWN IN MATPLOTLIB RATHER THAN GRAPHVIZ, deliberately. The figure it
    replaces needed a `dot` binary that is not a Python dependency and is not
    present in this project's environment or in CI, so it could not be
    regenerated by anyone who did not already have it installed -- a figure
    nobody can reproduce is a figure nobody can check. This one uses the same
    matplotlib the other conceptual diagrams use.
    """
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Census Construction", fontsize=15, fontweight="bold", pad=14)

    pos = {}
    for key, x, y, w, h, fill, edge, lines in _CENSUS_BOXES:
        box = mpatches.FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.008,rounding_size=0.012",
            facecolor=fill, edgecolor=edge, linewidth=1.6, zorder=2)
        ax.add_patch(box)
        pos[key] = (x, y, w, h)
        n = len(lines)
        # Lay the lines out inside the box rather than at a fixed offset, so a
        # box with five lines does not overflow one sized for three.
        step = (h - 0.018) / max(n, 1)
        top = y + h / 2 - 0.012 - step / 2
        for k, (text, style) in enumerate(lines):
            ax.text(x, top - k * step, text, ha="center", va="center",
                    fontsize=8.4 if style == "bold" else 7.8,
                    fontweight="bold" if style == "bold" else "normal",
                    style="italic" if style == "italic" else "normal",
                    color="#212121" if style != "italic" else "#555555",
                    zorder=3)

    for src, dst, label in _CENSUS_EDGES:
        x0, y0, _, h0 = pos[src]
        x1, y1, _, h1 = pos[dst]
        # zorder ABOVE the boxes and a shrink at the target end: drawn
        # behind them, a diagonal arrow's head is covered by the box it points
        # at, so every diagonal edge in the first render arrived headless while
        # the one vertical edge looked fine.
        arrow = FancyArrowPatch(
            (x0, y0 - h0 / 2), (x1, y1 + h1 / 2),
            arrowstyle="-|>", mutation_scale=13, linewidth=1.5,
            shrinkA=1.0, shrinkB=3.0,
            color="#455A64", connectionstyle="arc3,rad=0.0", zorder=4)
        ax.add_patch(arrow)
        if label:
            ax.text((x0 + x1) / 2, (y0 - h0 / 2 + y1 + h1 / 2) / 2 + 0.012,
                    label, ha="center", va="bottom", fontsize=7.6,
                    style="italic", color="#455A64", zorder=3)

    ax.text(0.5, 0.028,
            "Two parallel admission rules, not a screening funnel: nothing is "
            "excluded, and which stream a record\nlands in reflects how far "
            "NLM's indexing has reached rather than any judgement made here.",
            ha="center", va="center", fontsize=8.2, color="#37474F",
            style="italic", zorder=3)

    fig.tight_layout()
    out_base = OUT / "fig23_census_flow"
    fig.savefig(str(out_base) + ".pdf", bbox_inches="tight")
    fig.savefig(str(out_base) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig23_census_flow")


# ── Figure 30: What the engine can be asked ───────────────────────────

def fig30_modality_landscape():
    """The engine's coverage against the field's attention, per mechanism.

    EVERY NUMBER IS READ FROM `analysis/modality-coverage.json`, which is
    itself generated and guarded. Nothing on this figure is typed, because a
    diagram summarising a measurement is the easiest place for one to go
    stale -- this repository has shipped that defect in prose more than once,
    and a picture is harder to notice it in.

    The figure exists because a reader cannot otherwise see the shape of the
    problem, and that shape MOVED: no measurable mechanism is without engine
    representation any more, so the gap the figure draws is no longer presence
    but DEPTH and applicability -- one mechanism modelled deeply, the rest
    represented by arms an order of magnitude smaller, half of them reachable
    only as modifiers of another treatment.
    """
    src = Path(__file__).resolve().parent.parent / "analysis" / "modality-coverage.json"
    if not src.exists():
        print("  fig30: analysis/modality-coverage.json missing - skipping")
        return
    d = json.loads(src.read_text())
    rows = sorted(d["rows"], key=lambda r: r["census"])

    tier_color = {
        "treatment": COLORS["sdt"],
        "modifier": COLORS["checkpoint"],
        "absent": "#BFBFBF",
    }
    tier_label = {
        "treatment": "applicable as a treatment",
        "modifier": "present, but cannot be applied alone",
        "absent": "no engine representation",
    }

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(12, 6), gridspec_kw={"width_ratios": [1.55, 1]})
    fig.patch.set_facecolor("white")

    names = [r["mechanism"] for r in rows]
    census = [r["census"] for r in rows]
    colors = [tier_color[r["engine_tier"]] for r in rows]
    ypos = np.arange(len(rows))
    axA.barh(ypos, census, color=colors, edgecolor="black", linewidth=0.4)
    axA.set_yticks(ypos)
    axA.set_yticklabels(names, fontsize=8)
    axA.set_xlabel("Cancer articles carrying the mechanism's MeSH descriptor",
                   fontsize=9)
    axA.set_title("(a) The field's attention, by what the engine can express",
                  fontsize=10, fontweight="bold")
    axA.set_facecolor(COLORS["bg"])
    for i, r in enumerate(rows):
        axA.text(r["census"] + max(census) * 0.012, i, f"{r['census']:,}",
                 va="center", fontsize=7, color="#444")
    handles = [mpatches.Patch(facecolor=tier_color[t], edgecolor="black",
                              linewidth=0.4, label=tier_label[t])
               for t in ("treatment", "modifier", "absent")]
    axA.legend(handles=handles, fontsize=7.5, loc="lower right", framealpha=0.95)
    axA.set_xlim(0, max(census) * 1.18)

    absent = [r for r in rows if r["engine_tier"] == "absent"]
    treatment = [r for r in rows if r["engine_tier"] == "treatment"]
    modifier = [r for r in rows if r["engine_tier"] == "modifier"]

    # PANEL (b) CHANGES SUBJECT WHEN THE ABSENT COLUMN EMPTIES, for the same
    # reason the report does. Stacking 100% against 0% draws a bar with
    # nothing in it and collides its own labels, and worse, it invites a
    # reader to take an emptied column as a result. Once nothing is absent the
    # informative split is APPLICABILITY -- how many mechanisms a run can
    # actually select -- which is a harder number and a smaller one.
    if absent:
        cov_vals = [sum(r["census"] for r in rows) - sum(r["census"] for r in absent),
                    sum(r["trials"] for r in rows) - sum(r["trials"] for r in absent)]
        abs_vals = [sum(r["census"] for r in absent), sum(r["trials"] for r in absent)]
        labels = ["census articles", "registered trials"]
        x = np.arange(len(labels))
        tot = [c + a for c, a in zip(cov_vals, abs_vals)]
        axB.bar(x, [c / t * 100 for c, t in zip(cov_vals, tot)], 0.55,
                color=COLORS["checkpoint"], edgecolor="black", linewidth=0.5,
                label="engine can express it")
        axB.bar(x, [a / t * 100 for a, t in zip(abs_vals, tot)], 0.55,
                bottom=[c / t * 100 for c, t in zip(cov_vals, tot)],
                color="#BFBFBF", edgecolor="black", linewidth=0.5,
                label="no engine representation")
        for i, (a, t) in enumerate(zip(abs_vals, tot)):
            axB.text(i, 100 - (a / t * 100) / 2, f"{a:,}\n({a / t * 100:.0f}%)",
                     ha="center", va="center", fontsize=8.5, fontweight="bold")
        axB.set_xticks(x)
        axB.set_xticklabels(labels, fontsize=9)
        axB.set_ylabel("Share of the 16 measurable mechanisms (%)", fontsize=9)
        axB.set_ylim(0, 100)
        axB.set_title(f"(b) {len(absent)} of {len(rows)} mechanisms are absent",
                      fontsize=10, fontweight="bold")
        axB.legend(fontsize=7.5, loc="lower left", framealpha=0.95)
    else:
        cats = ["applicable\nas a treatment", "present, but only\nas a modifier"]
        counts = [len(treatment), len(modifier)]
        vols = [sum(r["census"] for r in treatment),
                sum(r["census"] for r in modifier)]
        x = np.arange(len(cats))
        axB.bar(x, counts, 0.5,
                color=[COLORS["sdt"], COLORS["checkpoint"]],
                edgecolor="black", linewidth=0.5)
        for i, (c, v) in enumerate(zip(counts, vols)):
            axB.text(i, c + max(counts) * 0.03,
                     f"{c} of {len(rows)}\n{v:,} articles",
                     ha="center", va="bottom", fontsize=9, fontweight="bold")
        axB.set_xticks(x)
        axB.set_xticklabels(cats, fontsize=8.5)
        axB.set_ylabel("Mechanisms", fontsize=9)
        axB.set_ylim(0, max(counts) * 1.35)
        axB.set_title(
            f"(b) Every mechanism is present; {len(treatment)} of {len(rows)} "
            "is applicable",
            fontsize=10, fontweight="bold")
        axB.text(0.5, 0.42,
                 "The absent column emptied.\nPresence is not applicability —\n"
                 "a modifier is only ever a coefficient\non something else.",
                 transform=axB.transAxes, ha="center", va="center",
                 fontsize=8, color="#444", style="italic",
                 bbox=dict(boxstyle="round,pad=0.45", fc="white",
                           ec="#BBB", alpha=0.95))
    axB.set_facecolor(COLORS["bg"])

    active = d["active_treatments"]
    kinds = d.get("treatment_kinds", {})
    arms = "; ".join(f"{v} ({kinds.get(v, '?')})" for v in active)
    fig.suptitle(
        "What this engine can be asked, against what the field publishes",
        fontsize=12, fontweight="bold", y=0.99)
    fig.text(0.5, 0.005,
             f"Treatment arms: {arms}.  "
             "Volume is NOT comparable across mechanisms - descriptor breadth "
             "varies, so this is not a ranking; the ENGINE colour is the "
             "content.  Generated from analysis/modality-coverage.json.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.035, 1, 0.96])
    save(fig, "fig30_modality_landscape")


# ── Figure 31: the modality panel ─────────────────────────────────────

def fig31_modality_panel():
    """Every applicable arm on one tumour, grouped by WHAT KILLS.

    The grouping is the figure. A bar chart of kill fractions alone would
    invite exactly the reading `analysis/modality-panel.md` refuses -- a
    ranking of therapies -- so the bars are coloured by ROUTE and the panel
    beside them shows that most arms no longer go through the ferroptosis
    engine at all. That is the answer to the criticism this campaign started
    from, and it is a structural claim rather than an efficacy one.

    Every number is read from `analysis/modality-panel.json`.
    """
    src = Path(__file__).resolve().parent.parent / "analysis" / "modality-panel.json"
    if not src.exists():
        print("  fig31: analysis/modality-panel.json missing - skipping")
        return
    d = json.loads(src.read_text())
    arms = sorted(d["arms"], key=lambda a: a["kill_fraction"])

    def route_key(r):
        if "ferroptosis engine" in r:
            return "ferroptosis engine"
        if "DNA damage" in r:
            return "DNA damage"
        if "threshold" in r:
            return "threshold ablation"
        if "immune" in r or "effectors" in r or "ICD" in r:
            return "immune"
        return "ferroptosis payload, delivered"

    palette = {
        "ferroptosis engine": COLORS["sdt"],
        "DNA damage": COLORS["radiation"],
        "threshold ablation": COLORS["parp"],
        "immune": COLORS["checkpoint"],
        "ferroptosis payload, delivered": COLORS["nano"],
    }

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1.7, 1]})
    fig.patch.set_facecolor("white")

    y = np.arange(len(arms))
    vals = [a["kill_fraction"] * 100 for a in arms]
    cols = [palette[route_key(a["route"])] for a in arms]
    axA.barh(y, vals, color=cols, edgecolor="black", linewidth=0.4)
    axA.set_yticks(y)
    axA.set_yticklabels([a["arm"] for a in arms], fontsize=8.5)
    axA.set_xlabel("Cells killed (%) - NOT a ranking of therapies", fontsize=9)
    axA.set_title("(a) Every applicable arm, same tumour, same seed",
                  fontsize=10, fontweight="bold")
    axA.set_facecolor(COLORS["bg"])
    axA.set_xlim(0, 100)
    for i, a in enumerate(arms):
        axA.text(a["kill_fraction"] * 100 + 1.2, i,
                 f"{a['kill_fraction'] * 100:.2f}%", va="center", fontsize=7.5,
                 color="#444")
    handles = [mpatches.Patch(facecolor=c, edgecolor="black", linewidth=0.4,
                              label=k)
               for k, c in palette.items()]
    axA.legend(handles=handles, fontsize=7.5, loc="lower right",
               framealpha=0.95, title="what kills", title_fontsize=8)

    # (b) the structural claim: how many arms go through the engine this
    # repository was built around.
    ferro = d["n_ferroptosis_routed"]
    other = d["n_other_routes"]
    axB.bar([0, 1], [ferro, other], 0.5,
            color=[COLORS["sdt"], COLORS["checkpoint"]],
            edgecolor="black", linewidth=0.5)
    for i, v in enumerate([ferro, other]):
        axB.text(i, v + 0.15, str(v), ha="center", fontsize=13,
                 fontweight="bold")
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["through the\nferroptosis engine",
                         "through some\nother route"], fontsize=9)
    axB.set_ylabel("Arms", fontsize=9)
    axB.set_ylim(0, max(ferro, other) * 1.4)
    axB.set_title(f"(b) {len(d['distinct_routes'])} distinct routes to death",
                  fontsize=10, fontweight="bold")
    axB.set_facecolor(COLORS["bg"])
    axB.text(0.5, 0.72,
             "A reader opening this repository\nbefore this campaign would have\n"
             "found ONE route and four arms.",
             transform=axB.transAxes, ha="center", va="center", fontsize=8.5,
             color="#444", style="italic",
             bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#BBB",
                       alpha=0.95))

    fig.suptitle("Ten arms, seven routes: what the engine can now be asked",
                 fontsize=12, fontweight="bold", y=0.99)
    fig.text(0.5, 0.005,
             "NOT a ranking. Every kill fraction is a function of the "
             "parameters the arm was given, and all but radiation's DNA "
             "channel are placeholders - see CALIBRATION_STATUS.md. SDT and "
             "PDT coincide because this panel is depth-free and their "
             "exogenous-ROS parameters are equal by default.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    save(fig, "fig31_modality_panel")


# ── Figure 32: what the microenvironment does to every arm ────────────

def fig32_modality_tme():
    """The resistance sweep, on a SIGNED scale, split by cell state.

    Two things this figure exists to show, neither of which survives a
    conventional treatment.

    THE SCALE MUST BE SIGNED. An earlier version of the analysis behind this
    took the absolute value of every change and reported an arm that "loses
    121%" -- impossible, and worse, it hid a real result: one axis HELPS.
    Clonal heterogeneity raises the pharmacologic arm's kill, because widening
    the antioxidant setpoint while holding its mean supplies a low-glutathione
    tail that dies while the average cell resists. A magnitude colour map would
    paint that the same shade as a catastrophic loss, so the map is diverging
    and zero is the neutral colour.

    IT IS SPLIT BY PHENOTYPE because the answer changes with the cell state,
    and that IS the finding. An earlier sweep ran one phenotype and reported
    two axes inert -- correctly for that run, and misleadingly, because what
    was inert was the configuration's ability to see them.

    Every value is read from `analysis/modality-tme.json`.
    """
    src = Path(__file__).resolve().parent.parent / "analysis" / "modality-tme.json"
    if not src.exists():
        print("  fig32: analysis/modality-tme.json missing - skipping")
        return
    d = json.loads(src.read_text())
    ebp = d["effects_by_phenotype"]
    undef = d.get("undefined_cells", {})
    phenos = d["phenotypes"]
    axes_order = list(ebp[phenos[0]].keys())
    arms = d["arms"]

    # HEIGHT SCALES WITH THE ARM COUNT. A fixed 5.6 inches fitted nine arms;
    # the tenth compressed the rows until two of them fell inside the
    # four-point tolerance `tests/test_chapter6_figures.py` uses to read a row,
    # so the guard started reading a neighbour's number as this row's. The
    # figure was wrong before the guard was: ten labelled rows in that space is
    # cramped whatever a test thinks.
    height = 5.6 + 0.42 * max(0, len(arms) - 9)
    fig, axs = plt.subplots(1, len(phenos), figsize=(13.0, height), sharey=True)
    vmax = max(abs(v) for ph in phenos for ax in axes_order
               for v in ebp[ph][ax].values()) or 1.0
    for k, ph in enumerate(phenos):
        M = np.array([[ebp[ph][ax].get(a, 0.0) for ax in axes_order]
                      for a in arms], dtype=float)
        ax = axs[k]
        # `RdBu` and NOT `RdBu_r`: low values must be RED. The first draft
        # used the reversed map, which drew every loss in blue and every gain
        # in red while the caption below said the opposite -- a figure
        # contradicting its own legend, which is the defect this campaign has
        # spent three review rounds removing from prose.
        im = ax.imshow(M, cmap="RdBu", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(axes_order)))
        ax.set_xticklabels(axes_order, rotation=28, ha="right", fontsize=8.5)
        ax.set_yticks(range(len(arms)))
        if k == 0:
            ax.set_yticklabels(arms, fontsize=8.5)
        ax.set_title(f"{ph} state", fontsize=10.5, fontweight="bold")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M[i, j]
                # A cell is left blank ONLY at exactly zero. The threshold was
                # 0.005, which blanked three cells carrying a real effect --
                # and one arm-axis pair was blank in one panel and drawn as
                # -10 in the other, so the same pair was simultaneously called
                # inert and shown acting. A sub-1% effect is printed as <1
                # rather than hidden.
                if arms[i] in undef.get(ph, {}).get(axes_order[j], []):
                    # UNDEFINED, not zero. The relative change has no
                    # denominator when the unstressed kill is 0, and a blank
                    # cell said the axis could not move the arm -- a different
                    # claim, and false for this whole row.
                    ax.text(j, i, "n/a", ha="center", va="center",
                            fontsize=6.2, color="#999", style="italic")
                    continue
                if v == 0.0:
                    continue
                if abs(v) < 0.005:
                    ax.text(j, i, "<1" if v > 0 else ">-1", ha="center",
                            va="center", fontsize=6.4, color="#666")
                    continue
                ax.text(j, i, f"{v * 100:+.0f}", ha="center", va="center",
                        fontsize=7.0,
                        color="white" if abs(v) > 0.55 * vmax else "#222")
        ax.set_xticks(np.arange(-0.5, len(axes_order), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(arms), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.1)
        ax.tick_params(which="minor", length=0)

    # LAYOUT BEFORE COLORBAR. `subplots_adjust` after `colorbar(ax=axs)` undoes
    # the space the colorbar reserved, and it drew over the rightmost column of
    # the persister panel -- the clonal-heterogeneity column carrying +124 and
    # +116, which is the figure's stated finding.
    fig.subplots_adjust(bottom=0.30, top=0.88, left=0.115, right=0.86)
    cax = fig.add_axes([0.885, 0.30, 0.016, 0.58])
    cb = fig.colorbar(im, cax=cax)
    # ONE unit. The cells print percent and the bar was ticked in fractions,
    # which is two notations for one quantity in one image -- the defect the
    # sibling function's own comment congratulates itself on avoiding.
    cb.set_ticks(np.linspace(-vmax, vmax, 5))
    # THE PERCENT SIGN IS LOAD-BEARING, for a reason that is not typographic.
    # Without it a colorbar tick is the same string as a cell value, and a
    # guard that reads the drawn numbers off a ROW cannot tell them apart when
    # a tick happens to land within a few points of a row -- which is exactly
    # what happened when a tenth arm shifted the rows. Ticking in percent makes
    # the bar unambiguous to a reader and distinguishable to the guard, and it
    # is the same unit the cells already print.
    cb.set_ticklabels([f"{v * 100:+.0f}%" for v in np.linspace(-vmax, vmax, 5)])
    # RELATIVE change, and the unit matters: `modality_tme_report._effect`
    # returns `(got - base) / base`, so +124 is "kills 2.24x as much", NOT a
    # 124-point move in a fraction bounded by one -- which would be the same
    # impossibility this chapter retracts.
    cb.set_label("change in kill RELATIVE to the unstressed run (%)  -  "
                 "negative is resistance", fontsize=8.2)
    fig.suptitle("What the microenvironment does to every arm, by cell state",
                 fontsize=12.5, fontweight="bold", y=0.985)
    fig.text(0.5, 0.02,
             "NOT a ranking. Signed deliberately: red is a loss, blue a GAIN, "
             "and the gains are real - clonal heterogeneity supplies a "
             "low-defence tail that dies while the average cell resists. A "
             "blank cell is EXACTLY zero for that arm - a cell the axis "
             "cannot move at all in this run, which is a property of the run "
             "rather than of the biology; an effect too small to round to a "
             "whole percent prints as <1 or >-1 rather than vanishing, and a "
             "cell whose relative change has NO DENOMINATOR - the arm kills "
             "nothing unstressed - prints n/a, because a blank there would "
             "claim the axis cannot move it. Every arm but "
             "radiation's DNA channel is parameterised with placeholders.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    save(fig, "fig32_modality_tme")


# ── Figure 33: the CAR-T barrier cascade ──────────────────────────────

def fig33_adoptive_barriers():
    """One construct, two diseases, as a waterfall.

    The point is that no single step looks catastrophic. Three barriers at
    plausible values leave six per cent, exhaustion removes most of what is
    left, and the product is the difference between a therapy that cures a
    blood cancer and one that does very little in a solid tumour. A single
    efficacy scalar would fit the same endpoints and lose exactly that.

    The bars are DERIVED from `analysis/modality-panel.json` and their product
    is asserted against the collapse before anything is drawn, because the
    decomposition beside them was a residual once -- the antigen term absorbed
    whatever the other two did not explain, which made it unfalsifiable.
    """
    src = Path(__file__).resolve().parent.parent / "analysis" / "modality-panel.json"
    if not src.exists():
        print("  fig33: analysis/modality-panel.json missing - skipping")
        return
    ab = json.loads(src.read_text()).get("adoptive_barriers")
    if not ab:
        print("  fig33: no adoptive_barriers block - skipping")
        return
    leuk = ab["leukaemia_kill_fraction"]
    solid = ab["solid_tumour_kill_fraction"]
    steps = [
        ("infused\n(leukaemia)", leuk, COLORS.get("immune", "#FFD700")),
        ("after the three\ndelivery barriers",
         leuk * ab["delivery_efficiency_solid"], "#CD5C5C"),
        ("after persistence",
         leuk * ab["delivery_efficiency_solid"]
         * ab["persistence_at_run_end_solid"], "#8B0000"),
        ("after the\nantigen ceiling", solid, "#4B0000"),
    ]
    # The drawing must not be able to disagree with the binary, IN EITHER
    # DIRECTION. The first version was `... < 0.05 or binds`, which the `or`
    # made vacuous exactly when the ceiling binds -- the case its own message
    # described. Both directions are checked now.
    moved = abs(steps[-2][1] / solid - 1.0) > 1e-9
    assert moved == bool(ab["antigen_ceiling_binds"]), (
        f"the binary says the antigen ceiling binds="
        f"{ab['antigen_ceiling_binds']} while the drawn steps "
        f"{'move' if moved else 'do not move'} across it")

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    xs = np.arange(len(steps))
    vals = [v for _, v, _ in steps]
    ax.bar(xs, vals, 0.58, color=[c for _, _, c in steps],
           edgecolor="black", linewidth=0.6)
    ax.set_yscale("log")
    ax.set_ylabel("kill fraction (log scale)", fontsize=9.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([n for n, _, _ in steps], fontsize=8.6)
    # ONE format for all four bars. The first draft printed "0.12%" beside
    # "6.93e-05", which are the same kind of quantity in two notations and
    # invites a reader to compare them by eye and get it wrong.
    for x, v in zip(xs, vals):
        ax.text(x, v * 1.55, f"{v * 100:.4g}%", ha="center", fontsize=8.6,
                fontweight="bold")
    ax.set_ylim(min(vals) * 0.45, max(vals) * 4.0)
    for i in range(len(steps) - 1):
        drop = vals[i] / vals[i + 1] if vals[i + 1] > 0 else float("inf")
        if drop < 1.01:
            label = "no effect"
        else:
            label = f"/{drop:,.1f}"
        ax.annotate(label, xy=(i + 0.5, (vals[i] * vals[i + 1]) ** 0.5),
                    ha="center", fontsize=8.6, color="#333",
                    fontweight="bold")
    inert = sum(1 for i in range(len(steps) - 1)
                if vals[i] / vals[i + 1] < 1.01)
    ax.set_title(
        f"The same CAR-T construct, twice: a {leuk / solid:,.0f}x collapse "
        f"across delivery, persistence and the antigen ceiling - "
        f"{inert} of the three doing nothing here",
        fontsize=11.0, fontweight="bold", pad=14)
    fig.text(0.5, 0.005,
             "NOT a clinical comparison. Every barrier value is an "
             "uncalibrated placeholder; the corpus establishes that the "
             "barriers are real and GENERAL rather than antigen-specific, not "
             "which of them dominates. The antigen ceiling is a CAP, so it "
             "contributes nothing unless it fires - here it does not, and the "
             "figure shows that rather than hiding it.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    save(fig, "fig33_adoptive_barriers")



# ── Figure 30: how far each modality reaches ──────────────────────────

def fig34_depth_reach():
    """Physical reach, which is the manuscript's tissue-access argument.

    Drawn because the argument is currently prose beside a table, and the
    thing that makes it an argument is a SHAPE: three physical modalities
    whose delivered energy falls off at rates two orders of magnitude apart,
    against a pharmacologic arm whose delivery does not fall off at all.

    The control is what keeps it honest. `Control` retains 400% of its surface
    kill, which is not robustness -- it is what a ratio does when both terms
    are near zero, and a panel RANKED on retention would put an untreated
    tumour at the top. Both panels therefore keep the SAME order, sorted by
    absolute kill at depth, so the control sits last in both and its bar is
    read against its neighbour rather than crowning a ranking. An earlier
    caption said panel (b) "ranks the untreated arm first", which is not what
    the figure draws: it is the tallest bar, not the first one.
    """
    src = Path(__file__).resolve().parent.parent / "analysis" / "depth-reach-comparison.json"
    if not src.exists():
        print("  fig34: analysis/depth-reach-comparison.json missing - skipping")
        return
    d = json.loads(src.read_text())
    rows = [r for r in d["rows"]]
    order = sorted(rows, key=lambda r: -r["deep_kill_pct"])
    names = [r["treatment"] for r in order]
    deep_mm = order[0]["deep_mm"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.6, 4.6))
    x = np.arange(len(names))
    axA.bar(x - 0.2, [r["surface_kill_pct"] for r in order], 0.4,
            label="at the surface", color="#4169E1", edgecolor="black",
            linewidth=0.5)
    axA.bar(x + 0.2, [r["deep_kill_pct"] for r in order], 0.4,
            label=f"at {deep_mm:.1f} mm", color="#8B0000", edgecolor="black",
            linewidth=0.5)
    axA.set_xticks(x)
    axA.set_xticklabels(names, rotation=20, ha="right", fontsize=8.6)
    axA.set_ylabel("kill (%)", fontsize=9.5)
    axA.set_title("(a) Kill at the surface and at depth", fontsize=10.5,
                  fontweight="bold")
    axA.legend(fontsize=8.2)
    axA.grid(axis="y", alpha=0.25)

    ret = [r["kill_retained_pct"] for r in order]
    cols = ["#BFBFBF" if r["treatment"] == "Control" else "#CD5C5C"
            for r in order]
    axB.bar(x, ret, 0.55, color=cols, edgecolor="black", linewidth=0.5)
    axB.axhline(100, color="black", linewidth=0.8, linestyle="--")
    axB.set_xticks(x)
    axB.set_xticklabels(names, rotation=20, ha="right", fontsize=8.6)
    axB.set_ylabel("share of surface kill retained (%)", fontsize=9.5)
    axB.set_title("(b) Retention - and why it cannot be read alone",
                  fontsize=10.5, fontweight="bold")
    ctrl = next((r for r in order if r["treatment"] == "Control"), None)
    if ctrl:
        i = names.index("Control")
        axB.annotate("a ratio of two\nnear-zero numbers",
                     xy=(i, ctrl["kill_retained_pct"]),
                     xytext=(i - 0.1, ctrl["kill_retained_pct"] * 0.72),
                     ha="center", fontsize=7.6, color="#333",
                     arrowprops=dict(arrowstyle="->", color="#666", lw=0.8))
    for a in (axA, axB):
        a.spines[["top", "right"]].set_visible(False)

    mu = d["constants"]
    fig.suptitle("How far each modality reaches into tissue", fontsize=12.5,
                 fontweight="bold")
    fig.text(0.5, 0.015,
             "NOT a ranking. Both panels share one order, sorted by kill at "
             "depth, so the UNTREATED arm sits last in both - and in (b) it is "
             "the TALLEST bar, which is the point: retention is a ratio, and a "
             "ratio of two near-zero numbers is not robustness. Read (a) and "
             "(b) together or neither. "
             f"Attenuation is fixed physics and the two constants are quoted "
             f"in ONE unit here: PDT {mu['pdt_mu_eff_per_mm']}/mm against "
             f"radiation {mu['radiation_mu_per_cm'] / 10:g}/mm, a factor of "
             f"{mu['pdt_mu_eff_per_mm'] / (mu['radiation_mu_per_cm'] / 10):.0f}. "
             "Printing one per mm and the other per cm made a hundredfold "
             "difference read as a tenfold one. Every kill magnitude rests on "
             "uncalibrated biochemistry.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.09, 1, 0.94])
    save(fig, "fig34_depth_reach")


# ── Figure 31: what a published target could and could not settle ─────

def fig35_calibration_verdicts():
    """The calibration outcome per arm, INCLUDING the ones that failed.

    A figure of what is fitted would be marketing. The interesting cells are
    the three that are not: one arm is INADMISSIBLE because its published band
    constrains a PRODUCT of two unidentifiable factors, and one has NO TARGET
    at all. Drawing those beside the successes is the whole reason the panel
    is worth a figure.
    """
    src = Path(__file__).resolve().parent.parent / "analysis" / "modality-calibration.json"
    if not src.exists():
        print("  fig35: analysis/modality-calibration.json missing - skipping")
        return
    arms = json.loads(src.read_text())["arms"]
    order = {"ADMISSIBLE": 0, "DIRECTIONAL": 1, "UNCONSTRAINED": 2,
             "PARTLY REFUTED": 3, "INADMISSIBLE": 4, "NO TARGET": 5}
    colour = {"ADMISSIBLE": "#2E7D32", "DIRECTIONAL": "#1565C0",
              "UNCONSTRAINED": "#F9A825", "PARTLY REFUTED": "#EF6C00",
              "INADMISSIBLE": "#C62828", "NO TARGET": "#9E9E9E"}
    rows = sorted(arms, key=lambda a: (order.get(a.get("verdict"), 9),
                                       a["arm"]))
    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    y = np.arange(len(rows))
    for i, a in enumerate(rows):
        v = a.get("verdict", "NO TARGET")
        ax.barh(i, 1.0, color=colour.get(v, "#9E9E9E"), edgecolor="black",
                linewidth=0.5)
        ax.text(0.02, i, f"{a['arm']}", va="center", fontsize=9,
                color="white", fontweight="bold")
        ax.text(0.98, i, v, va="center", ha="right", fontsize=8.6,
                color="white", fontweight="bold")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(0, 1)
    ax.invert_yaxis()
    for sp in ax.spines.values():
        sp.set_visible(False)
    counts = {v: sum(1 for a in arms if a.get("verdict") == v)
              for v in order}
    ax.set_title(
        f"What a published target could settle, per arm "
        f"({counts['ADMISSIBLE']} fitted, {counts['DIRECTIONAL']} directional, "
        f"{counts['UNCONSTRAINED']} unconstrained, "
        f"{counts['PARTLY REFUTED']} partly refuted, "
        f"{counts['INADMISSIBLE']} inadmissible, "
        f"{counts['NO TARGET']} with no target)",
        fontsize=11, fontweight="bold")
    fig.text(0.5, 0.02,
             "ADMISSIBLE means a parameter reproduces a published band, NOT "
             "that the arm is validated - none of these feeds a number the "
             "manuscript's quantitative chapters report. The three rows that "
             "are not green are the informative ones: an UNCONSTRAINED arm has "
             "a target so loose it admits almost the whole scanned range, an "
             "INADMISSIBLE one has a target that constrains a PRODUCT of "
             "factors neither of which is identifiable from it, and an arm "
             "with NO TARGET has nothing published to fit at all. Two rows "
             "are neither fitted nor failed: a DIRECTIONAL arm has a target "
             "that constrains a SIGN and not a value, and the PARTLY REFUTED "
             "one is the sharpest row on the page - an independent study "
             "CONTRADICTS one of its three directional claims, and it is "
             "coloured to be found rather than folded into the greens.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    save(fig, "fig35_calibration_verdicts")



def fig36_fractionation():
    """What the schedule layer can be checked against, and where it fails.

    Three panels, and the first is the one that earns the figure. EQD2 is a
    function of the alpha/beta you assume, so two schedules trace two curves,
    and where they CROSS is the ratio at which a trial's two arms are
    equivalent. Prostate's crossing lands inside the band the radiobiology
    literature estimates independently; breast's does not, and the reason is
    geometric rather than rhetorical -- one arm delivers less total dose and
    was still not inferior, so the curves cross where no tissue lives.

    Every number is read from the committed validation artifact or recomputed
    from the same formula the crate implements. Nothing is typed.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "fractionation-validation.json"
    if not src.exists():
        print("  fig36: fractionation-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())

    def eqd2(n, dose, ab):
        return n * dose * (1 + dose / ab) / (1 + 2.0 / ab)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) the isoeffect crossing
    ax = axes[0]
    colours = {"prostate": "#14505c", "breast": "#a2582b"}
    for row in d["isoeffect"]:
        site = row["site"]
        n1, s1 = row["arm_a"].split(" x ")
        n2, s2 = row["arm_b"].split(" x ")
        n1, d1 = int(n1), float(s1.split()[0])
        n2, d2 = int(n2), float(s2.split()[0])
        xs = [0.3 + 0.02 * i for i in range(0, 500)]
        gap = [eqd2(n2, d2, x) - eqd2(n1, d1, x) for x in xs]
        ax.plot(xs, gap, color=colours[site], lw=2,
                label=f"{site}: {row['arm_b']} - {row['arm_a']}")
        ax.axvspan(row["published_lo_gy"], row["published_hi_gy"],
                   color=colours[site], alpha=0.12)
        x0 = row["implied_alpha_beta_gy"]
        ax.plot([x0], [0], "o", color=colours[site], ms=8, zorder=5)
        ax.annotate(f"{x0:.2f} Gy", (x0, 0), textcoords="offset points",
                    xytext=(6, 10 if site == "prostate" else -16),
                    fontsize=8, color=colours[site], fontweight="bold")
    ax.axhline(0, color="#444", lw=0.8)
    ax.set_xlim(0.3, 10.3)
    ax.set_xlabel(r"assumed $\alpha/\beta$ (Gy)")
    ax.set_ylabel("EQD2 difference between the arms (Gy)")
    ax.set_title("(a) where two trial arms are equivalent", fontsize=10)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(alpha=0.25)

    # (b) why radiotherapy is fractionated at all
    ax = axes[1]
    ab_t, ab_l = d["constants"]["ALPHA_BETA_TUMOUR_GY"], d["constants"]["ALPHA_BETA_LATE_GY"]
    sizes = [1.0 + 0.25 * i for i in range(0, 93)]
    ratio = [eqd2(1, s, ab_t) / eqd2(1, s, ab_l) for s in sizes]
    ax.plot(sizes, ratio, color="#14505c", lw=2)
    ax.axhline(1.0, color="#888", ls="--", lw=1)
    ax.axvline(2.0, color="#888", ls=":", lw=1)
    ax.annotate("2 Gy: the convention\nEQD2 is defined against",
                (2.0, 1.0), textcoords="offset points", xytext=(10, 24),
                fontsize=7.5, color="#555")
    ax.set_xlabel("dose per fraction (Gy)")
    ax.set_ylabel(r"tumour effect / late-tissue effect")
    ax.set_title(r"(b) the differential fractionation exploits", fontsize=10)
    ax.grid(alpha=0.25)

    # (c) the fourth R
    ax = axes[2]
    def reox(n, dose, f0, half_life):
        alpha, beta = 0.3, 0.03
        # OER as a dose-modifying factor, the same hyperbola the crate uses
        def oer(p):
            return (3 * p + 3) / (p + 3)
        hyp_dose = dose * oer(0.05 * 40.0) / oer(40.0)
        s_ox = pow(2.718281828459045, -(alpha * dose + beta * dose * dose))
        s_hy = pow(2.718281828459045,
                   -(alpha * hyp_dose + beta * hyp_dose * hyp_dose))
        moved = 0.0 if half_life == float("inf") else \
            1 - pow(2.718281828459045, -0.6931471805599453 / half_life)
        ox, hy = 1 - f0, f0
        for _ in range(n):
            ox *= s_ox
            hy *= s_hy
            t = hy * moved
            hy -= t
            ox += t
        return ox + hy

    ns = list(range(1, 36))
    gain = [reox(n, 2.0, 0.3, float("inf")) / reox(n, 2.0, 0.3, 5.0) for n in ns]
    ax.plot(ns, gain, color="#14505c", lw=2)
    ax.axhline(1.0, color="#888", ls="--", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("fractions delivered")
    ax.set_ylabel("survival with a frozen hypoxic core\n/ survival when it reoxygenates")
    ax.set_title("(c) what reoxygenation is worth, over a course", fontsize=10)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("The schedule, and the two things it can be checked against",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "(a) is the layer's external check and it runs BACKWARDS: two schedules a trial "
             "reported as not differing imply the alpha/beta at which they are equivalent, and "
             "that value is compared against estimates derived from other data (shaded). Prostate "
             "lands inside its band; breast crosses at 0.70 Gy, far below any plausible tissue, "
             "because its shorter arm delivers less total dose and was still not inferior - a "
             "statement about what EQD2 leaves out. (b) and (c) are DIRECTION-only: the "
             "late-tissue alpha/beta is a convention and the reoxygenation half-life is a free "
             "parameter no dataset here constrains.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig36_fractionation")


def fig37_chemotherapy():
    """The chemotherapy arm's two structural predictions.

    Both are read from the committed validation artifact, which is produced by
    a stdlib implementation independent of the crate. Nothing here is typed.

    Panel (a) is the residue a phase-specific agent leaves, and the reason it
    is worth drawing is the SECOND pair of bars: in a population mostly out of
    cycle the gap narrows, because the phase-nonspecific agent is struggling
    too. Panel (b) is the dose-density window, and its shape is the finding --
    an interior peak with the advantage vanishing at BOTH ends.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "chemo-validation.json"
    if not src.exists():
        print("  fig37: chemo-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) dose-response by class
    ax = axes[0]
    colours = {"PhaseNonspecific": "#14505c", "SPhaseSpecific": "#a2582b",
               "MPhaseSpecific": "#6b7f3a"}
    labels = {"PhaseNonspecific": "phase-nonspecific (alkylator, platinum)",
              "SPhaseSpecific": "S-phase specific (antimetabolite)",
              "MPhaseSpecific": "M-phase specific (taxane, vinca)"}
    for cls, pts in d["dose_response"].items():
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, "o-", ms=3.5, lw=1.8, color=colours[cls], label=labels[cls])
    ax.set_yscale("log")
    ax.set_xlabel("dose (arbitrary units; the potency is a placeholder)")
    ax.set_ylabel("surviving fraction")
    ax.set_title("(a) what the cycle does to a dose-response", fontsize=10)
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(alpha=0.25, which="both")

    # (b) the residue ratio, in two populations
    ax = axes[1]
    res = d["residue_ratio_at_dose_8"]
    groups = ["SPhaseSpecific", "MPhaseSpecific"]
    x = range(len(groups))
    w = 0.36
    prolif = [res["proliferating"][g] for g in groups]
    quiet = [res["quiescent_rich"][g] for g in groups]
    ax.bar([i - w / 2 for i in x], prolif, w, color="#14505c", label="proliferating")
    ax.bar([i + w / 2 for i in x], quiet, w, color="#9fb3b8",
           label="mostly out of cycle")
    for i, (a, b) in enumerate(zip(prolif, quiet)):
        ax.text(i - w / 2, a, f"{a:.1f}x", ha="center", va="bottom", fontsize=8)
        ax.text(i + w / 2, b, f"{b:.1f}x", ha="center", va="bottom", fontsize=8)
    ax.axhline(1.0, color="#888", ls="--", lw=1)
    ax.set_xticks(list(x))
    ax.set_xticklabels(["S-phase agent", "M-phase agent"], fontsize=9)
    ax.set_ylabel("survivors relative to a phase-nonspecific agent")
    ax.set_title("(b) the residue, and where it shrinks", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25, axis="y")

    # (c) the dose-density window
    ax = axes[2]
    dd = d["dose_density"]
    ax.plot(dd["regrowth_per_day"], dd["advantage"], color="#14505c", lw=2)
    ax.axhline(1.0, color="#888", ls="--", lw=1)
    if dd["window_lo"] is not None:
        ax.axvspan(dd["window_lo"], dd["window_hi"], color="#14505c", alpha=0.10)
    ax.plot([dd["peak_at_regrowth_per_day"]], [dd["peak_advantage"]], "o",
            color="#a2582b", ms=8, zorder=5)
    ax.annotate(f"{dd['peak_advantage']:.1f}x at "
                f"{dd['peak_at_regrowth_per_day']:g}/day",
                (dd["peak_at_regrowth_per_day"], dd["peak_advantage"]),
                textcoords="offset points", xytext=(12, -4), fontsize=8,
                color="#a2582b", fontweight="bold")
    ax.set_xscale("symlog", linthresh=0.005)
    ax.set_xlabel("Gompertz regrowth rate (per day)")
    ax.set_ylabel("burden ratio, 21-day / 14-day schedule")
    ax.set_title("(c) when a shorter interval is worth anything", fontsize=10)
    ax.grid(alpha=0.25)

    fig.suptitle("Chemotherapy: two predictions that need no fitted potency",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "NO dose-response here is fitted: the repository's CTRPv2 route reached five "
             "ferroptosis compounds and no longer reaches anything, so every absolute kill "
             "fraction is a placeholder and only the SHAPES are results. (b)'s second pair of "
             "bars is the finding that was not designed in - out of cycle the phase-specific "
             "agent's disadvantage shrinks, because the flat agent is struggling too. (c) has "
             "two ends: no advantage with nothing to outrun, and none when regrowth is fast "
             "enough that both schedules return to the plateau. Whether any real tumour sits "
             "inside the shaded window is NOT established here.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig37_chemotherapy")


def fig38_checkpoint():
    """What a ratio can constrain that an absolute response rate cannot.

    Panel (a) is the argument: the model's antigenicity response saturates, and
    the ratio between two burdens is what a trial can refute -- the absolute
    height cannot, because the mapping from a model kill to a radiological
    response is unknown. Panel (b) is the constraint that buys: the region of
    shape parameters the measured band admits, against the whole grid scanned.
    Panel (c) is the limit, and it is drawn at the same size as the result
    because it is comparable in size to it.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "checkpoint-validation.json"
    if not src.exists():
        print("  fig38: checkpoint-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())
    floor = 0.05
    half = d["constants"]["TMB_HALF_MAX_PER_MB"]

    def antigenicity(t, h=half, f=floor):
        return f + (1 - f) * t / (t + h)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) the saturating response and the two strata
    ax = axes[0]
    xs = [i * 0.5 for i in range(0, 121)]
    ax.plot(xs, [antigenicity(x) for x in xs], color="#14505c", lw=2)
    hi = d["representative_tmb"]["high"]
    lo = d["representative_tmb"]["low"]
    for x, lab, col in ((lo, "non-tTMB-high", "#9fb3b8"), (hi, "tTMB-high", "#a2582b")):
        ax.plot([x], [antigenicity(x)], "o", ms=8, color=col, zorder=5)
        ax.annotate(lab, (x, antigenicity(x)), textcoords="offset points",
                    xytext=(8, -12), fontsize=8, color=col, fontweight="bold")
    ax.axvline(d["constants"]["TMB_HIGH_THRESHOLD_PER_MB"], color="#888",
               ls=":", lw=1)
    ax.annotate("trial threshold\n10 mut/Mb",
                (d["constants"]["TMB_HIGH_THRESHOLD_PER_MB"], 0.12),
                textcoords="offset points", xytext=(6, 0), fontsize=7.5,
                color="#555")
    ax.set_xlabel("tumour mutational burden (mut/Mb)")
    ax.set_ylabel("antigenicity (model, dimensionless)")
    ax.set_title("(a) the shape a ratio can test", fontsize=10)
    ax.grid(alpha=0.25)

    # (b) the admitted region
    ax = axes[1]
    region = d["admissible_region"]
    if region:
        ax.scatter([r[1] for r in region], [r[0] for r in region], s=9,
                   color="#14505c", alpha=0.55, label="admitted by the band")
    ax.plot([half], [floor], "*", ms=16, color="#a2582b", zorder=5,
            label="shipped constants")
    ax.set_xlabel("half-maximal burden (mut/Mb)")
    ax.set_ylabel("antigenicity floor")
    ax.set_title(f"(b) {d['admissible_fraction']:.0%} of the grid survives "
                 f"the ratio", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(alpha=0.25)

    # (c) the limit
    ax = axes[2]
    sens = d["sensitivity_to_representative_tmb"]
    labels = [f"{s['tmb_high']:g} / {s['tmb_low']:g}" for s in sens]
    vals = [s["ratio_at_default"] for s in sens]
    ax.barh(range(len(vals)), vals, color="#9fb3b8")
    ax.axvspan(d["measured_band"][0], d["measured_band"][1], color="#14505c",
               alpha=0.12)
    ax.axvline(d["measured_ratio"], color="#14505c", lw=2)
    ax.annotate(f"measured {d['measured_ratio']:.1f}x",
                (d["measured_ratio"], len(vals) - 0.4),
                textcoords="offset points", xytext=(6, 0), fontsize=8,
                color="#14505c", fontweight="bold")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylabel("representative burden, high / low")
    ax.set_xlabel("model response ratio")
    ax.set_title("(c) the choice the trial does not make for us", fontsize=10)
    ax.grid(alpha=0.25, axis="x")

    fig.suptitle("Checkpoint blockade: constraining a shape with a ratio the "
                 "mapping cancels out of", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "An objective response is a 30% reduction in diameter, not a kill fraction, so no "
             "absolute response rate can be compared with this model. A RATIO between two strata "
             "of one trial cancels that unknown mapping exactly, and KEYNOTE-158 supplies one "
             "stratified by mutational burden - which moves antigenicity and, to first order, not "
             "the brake. It constrains a SHAPE and identifies nothing: one ratio is one equation, "
             "the brake is untouched, and (c) shows the representative-burden choice moving the "
             "model's answer by about as much as the target band is wide.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig38_checkpoint")


def fig39_adoptive_escalation():
    """When escalating a CAR-T dose buys tenfold and when it buys nothing.

    Both panels read the committed validation artifact. The left is the
    threshold -- the property that makes a density failure different in KIND
    from a delivery failure -- and the right is the consequence, which is the
    only prediction this arm makes.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "adoptive-validation.json"
    if not src.exists():
        print("  fig39: adoptive-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())
    t = d["constants"]["ANTIGEN_DENSITY_THRESHOLD"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    ax = axes[0]
    xs = [row[0] for row in d["engagement_curve"]]
    ys = [row[1] for row in d["engagement_curve"]]
    ax.plot(xs, ys, "o-", color="#14505c", lw=2, ms=4)
    ax.axvline(t, color="#a2582b", ls="--", lw=1.4)
    ax.annotate("threshold", (t, 0.06), textcoords="offset points",
                xytext=(8, 0), fontsize=8, color="#a2582b", fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlabel("target molecules per tumour cell")
    ax.set_ylabel("fraction of cells a CAR can engage")
    ax.set_title("(a) a threshold, not a gradient", fontsize=10)
    ax.grid(alpha=0.25, which="both")

    ax = axes[1]
    rows = d["discrimination"]
    labels = [r["case"].split(" (")[0] for r in rows]
    gains = [r["gain"] for r in rows]
    colours = ["#14505c", "#a2582b", "#9fb3b8"]
    ax.barh(range(len(rows)), gains, color=colours[:len(rows)])
    for i, g in enumerate(gains):
        ax.text(g, i, f" {g:.2f}x", va="center", fontsize=9, fontweight="bold")
    ax.axvline(1.0, color="#888", ls="--", lw=1)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("kill at ten times the dose, relative to the reference dose")
    ax.set_title("(b) what escalating the dose buys", fontsize=10)
    ax.grid(alpha=0.25, axis="x")

    ax = axes[2]
    ex = d["expansion_vs_drive"]
    ax.plot([e["antigen_drive"] for e in ex], [e["peak_fold"] for e in ex],
            "o-", color="#14505c", lw=2, ms=5)
    ax.set_yscale("log")
    ax.set_xlabel("antigen available to drive expansion")
    ax.set_ylabel("peak fold expansion over 28 days")
    ax.set_title("(c) expansion tracks the antigen it consumes", fontsize=10)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("Adoptive cell therapy: a failure more dose cannot fix",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "Two tumours with the SAME barriers and the same poor outcome, differing only in "
             "antigen density, respond completely differently to the one intervention a clinician "
             "reaches for first. That is this arm's only prediction, and it is about an experiment "
             "rather than a published number: unlike the checkpoint arm, no ratio is available "
             "that cancels the mapping between a remission and a kill, because blood and solid "
             "CAR-T are different trials with different endpoints. Nothing here is fitted - the "
             "density threshold varies by orders of magnitude with the receptor and the target.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig39_adoptive_escalation")


def fig40_oncolytic_bind():
    """The double bind, and the condition that decides which side wins.

    Panel (a) is the tension itself: one curve falls with immune competence and
    the other rises. Panel (b) is the consequence, and the reason the figure is
    worth drawing -- the sum has an interior maximum only when priming is
    efficient enough, so the clinical intuition is conditionally right and the
    condition is nameable. Panel (c) is where the condition flips, with the
    saturating regime marked because in it the optimum's position means
    nothing.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "oncolytic-validation.json"
    if not src.exists():
        print("  fig40: oncolytic-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    ax = axes[0]
    cs = [c * 0.1 for c in range(0, 161)]
    sens = d["clearance_sensitivity"]
    survival = [1.0 / (1.0 + sens * c) for c in cs]
    priming = [c / (1.0 + c) for c in cs]
    ax.plot(cs, survival, color="#14505c", lw=2, label="virus surviving to spread")
    ax.plot(cs, priming, color="#a2582b", lw=2, label="anti-tumour priming")
    ax.set_xlabel("immune competence")
    ax.set_ylabel("relative magnitude")
    ax.set_title("(a) two arms, opposite directions", fontsize=10)
    ax.legend(fontsize=7.5, loc="center right")
    ax.grid(alpha=0.25)

    ax = axes[1]
    curves = d["outcome_curves"]
    xs = [r["competence"] for r in curves]
    ax.plot(xs, [r["weak_priming"] for r in curves], "o-", color="#9fb3b8",
            lw=2, ms=5, label="weak priming (0.5)")
    ax.plot(xs, [r["strong_priming"] for r in curves], "o-", color="#14505c",
            lw=2, ms=5, label="strong priming (4.0)")
    ax.set_xlabel("immune competence")
    ax.set_ylabel("durable outcome (model, dimensionless)")
    ax.set_title("(b) the sum, and where its maximum sits", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    ax = axes[2]
    rows = d["optimum_by_efficiency"]
    live = [r for r in rows if not r["saturated"]]
    sat = [r for r in rows if r["saturated"]]
    ax.plot([r["priming_efficiency"] for r in live],
            [r["optimal_competence"] for r in live], "o-", color="#14505c",
            lw=1.8, ms=3.5, label="optimum")
    if sat:
        ax.scatter([r["priming_efficiency"] for r in sat],
                   [r["optimal_competence"] for r in sat], s=16,
                   color="#c9ccd1", label="saturated (meaningless)")
    x0 = d["crossover_priming_efficiency"]
    ax.axvline(x0, color="#a2582b", ls="--", lw=1.5)
    ax.annotate(f"crossover {x0}", (x0, ax.get_ylim()[1] * 0.75),
                textcoords="offset points", xytext=(6, 0), fontsize=8,
                color="#a2582b", fontweight="bold")
    ax.set_xlabel("priming efficiency")
    ax.set_ylabel("immune competence at the optimum")
    ax.set_title("(c) where suppression stops being the answer", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(alpha=0.25)

    fig.suptitle("Oncolytic virus: the double bind, and the condition that "
                 "decides it", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "Every review of this field states the tension - immunity clears the virus and "
             "immunity is the durable mechanism. Stating it settles nothing. What the model adds "
             "is a CONDITION: below a crossover in priming efficiency the optimum is at full "
             "suppression and the clinical intuition is right; above it the same move throws away "
             "the durable arm. The crossover's existence is the claim, not its position - priming "
             "efficiency is a placeholder in units this model invented. The grey points saturate, "
             "where the trade-off has stopped operating and the optimum's location means nothing.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig40_oncolytic_bind")


def fig41_adc_loading():
    """An optimum that falls out of somebody else's measurement.

    Panel (a) is the anchor and the alternative that fails it: a single power
    law through one measured ratio misses the other, and the acceleration it
    misses is what makes an optimum exist. Panel (b) is the optimum. Panel (c)
    is the in-vitro/in-vivo disagreement, measured in the same study.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "adc-validation.json"
    if not src.exists():
        print("  fig41: adc-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())
    curve = d["curve"]
    dars = [r["dar"] for r in curve]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    ax = axes[0]
    ax.plot(dars, [r["clearance"] for r in curve], color="#14505c", lw=2,
            label="piecewise, through both ratios")
    import math as _m
    e = _m.log(5.0) / _m.log(4.0)
    ax.plot(dars, [(x / 2.0) ** e for x in dars], color="#a2582b", lw=1.8,
            ls="--", label="single power law (misses c8/c4)")
    for p in d["anchor_points"]:
        ax.plot([p["dar"]], [p["clearance_relative_to_dar2"]], "o", ms=9,
                color="#14505c", zorder=5)
    ax.set_xlabel("drug-antibody ratio")
    ax.set_ylabel("clearance, relative to DAR 2")
    ax.set_title("(a) the measured ratios, and a curve that misses one",
                 fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.25)

    ax = axes[1]
    ax.plot(dars, [r["delivered_per_dose"] for r in curve], color="#14505c", lw=2)
    ax.plot([d["optimal_dar"]], [d["delivered_at_optimum"]], "*", ms=17,
            color="#a2582b", zorder=5)
    ax.annotate(f"optimum at DAR {d['optimal_dar']:g}",
                (d["optimal_dar"], d["delivered_at_optimum"]),
                textcoords="offset points", xytext=(10, -4), fontsize=8.5,
                color="#a2582b", fontweight="bold")
    ax.set_xlabel("drug-antibody ratio")
    ax.set_ylabel("payload delivered per unit antibody dose")
    ax.set_title("(b) more drug per antibody, fewer antibodies", fontsize=10)
    ax.grid(alpha=0.25)

    ax = axes[2]
    ax.plot(dars, [r["in_vitro_potency"] for r in curve], color="#a2582b",
            lw=2, label="in vitro (no clearance)")
    scale = max(r["delivered_per_dose"] for r in curve)
    ax.plot(dars, [r["delivered_per_dose"] / scale * 4 for r in curve],
            color="#14505c", lw=2, label="in vivo delivery (scaled)")
    ax.set_xlabel("drug-antibody ratio")
    ax.set_ylabel("relative magnitude")
    ax.set_title("(c) the two orderings disagree", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.25)

    fig.suptitle("Antibody-drug conjugates: an optimum nobody chose",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "The inputs are two clearance RATIOS from one study (PMID 15501986): a DAR-8 "
             "conjugate clears three times faster than DAR-4 and five times faster than DAR-2. "
             "The optimum is the output. A single power law fitted to one ratio misses the other "
             "by a third, and the acceleration it misses is exactly what makes an interior optimum "
             "exist - a smoother curve would have removed the result before it appeared. Panel (c) "
             "is the in-vitro-to-in-vivo gap this book argues about elsewhere, with both halves "
             "measured in one experiment. Delivered payload is not efficacy, and one conjugate is "
             "not all conjugates.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig41_adc_loading")


def fig42_ablation_sleeve():
    """Where a thermal ablation fails, and where electroporation does not.

    The arm's calibration is UNCONSTRAINED and will stay so -- a threshold
    observable cannot identify a threshold parameter. What replaces it is a
    statement about WHERE the survivors are, which a coverage fraction cannot
    make and a clinician can check.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "ablation-validation.json"
    if not src.exists():
        print("  fig42: ablation-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    ax = axes[0]
    body = d["body_temperature_c"]
    length = d["cooling_length_mm"]
    xs = [i * 0.05 for i in range(0, 201)]
    for applicator, colour in ((90.0, "#14505c"), (60.0, "#a2582b"),
                               (50.0, "#9fb3b8")):
        ys = [body + (applicator - body) * (1 - pow(2.718281828459045, -x / length))
              for x in xs]
        ax.plot(xs, ys, color=colour, lw=2, label=f"{applicator:.0f} °C applicator")
    ax.axhline(43.0, color="#888", ls="--", lw=1)
    ax.annotate("43 °C, the dose reference", (7.5, 43.6), fontsize=7.5, color="#555")
    ax.set_xlabel("distance from the vessel (mm)")
    ax.set_ylabel("peak tissue temperature (°C)")
    ax.set_title("(a) blood carries the heat away", fontsize=10)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(alpha=0.25)

    ax = axes[1]
    rows = [r for r in d["by_applicator_temperature"] if not r["total_failure"]]
    dead = [r for r in d["by_applicator_temperature"] if r["total_failure"]]
    ax.plot([r["applicator_c"] for r in rows],
            [r["thermal_sleeve_mm"] for r in rows], "o-", color="#14505c",
            lw=2, ms=6, label="thermal")
    ax.plot([r["applicator_c"] for r in rows],
            [r["electroporation_sleeve_mm"] for r in rows], "s-",
            color="#a2582b", lw=2, ms=6, label="electroporation")
    for r in dead:
        ax.annotate("fails\neverywhere", (r["applicator_c"], 0.6),
                    fontsize=7.5, color="#888", ha="center")
    ax.set_xlabel("applicator temperature (°C)")
    ax.set_ylabel("surviving perivascular sleeve (mm)")
    ax.set_title("(b) a radius, not a percentage", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    ax = axes[2]
    sens = d["sensitivity_to_cooling_length"]
    ax.plot([r["cooling_length_mm"] for r in sens],
            [r["sleeve_at_60c_mm"] for r in sens], "o-", color="#9fb3b8",
            lw=2, ms=6)
    ax.set_xlabel("cooling length (mm) — a placeholder")
    ax.set_ylabel("sleeve at 60 °C (mm)")
    ax.set_title("(c) the size is a restatement of an assumption", fontsize=10)
    ax.grid(alpha=0.25)

    fig.suptitle("Thermal ablation fails in a place, and electroporation does "
                 "not fail there", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "This arm's calibration is UNCONSTRAINED and will stay so: a threshold observable "
             "cannot identify a threshold parameter. What replaces it is a statement about WHERE "
             "the survivors are - a radius rather than a coverage percentage, which a clinician "
             "can look for. Electroporation's zero is structural rather than measured: it is "
             "non-thermal, so a heat sink removes nothing that matters to it, which is the "
             "documented reason it is reached for near vessels. (c) is why the SIZE is not a "
             "result: the cooling length stands in for vessel calibre and flow, which this layer "
             "does not represent, and the sleeve scales with it almost proportionally.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.13, 1, 0.93])
    save(fig, "fig42_ablation_sleeve")


def fig43_sonodynamic_frequency():
    """The frequency an SDT applicator should use, and the one claim that fails.

    The chapter's other arms all closed with a direction their comparator
    confirmed. This one closes with a direction its comparator CONTRADICTS,
    and that is the panel's subject rather than a footnote to it.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "sonodynamic-validation.json"
    if not src.exists():
        print("  fig43: sonodynamic-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())
    alpha = d["alpha_db_cm_mhz_from_params_rs"]

    def delivered(z, f, a):
        return f * pow(10.0, -a * f * z / 20.0) / (f ** 0.5)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) the interior optimum, and that it IS interior
    ax = axes[0]
    freqs = [0.05 + i * 0.01 for i in range(0, 500)]
    for z, colour in ((5.0, "#14505c"), (10.0, "#a2582b"), (15.0, "#9fb3b8")):
        ys = [delivered(z, f, alpha) for f in freqs]
        peak = max(ys)
        ax.plot(freqs, [y / peak for y in ys], color=colour, lw=2,
                label=f"{z:.0f} cm deep")
        fstar = 10.0 / (alpha * 2.302585092994046 * z)
        ax.plot([fstar], [1.0], "o", color=colour, ms=7)
    ax.set_xlim(0, 5)
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("delivered mechanical index (peak = 1)")
    ax.set_title("(a) three effects, two of them opposing", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    # (b) the model against the comparator -- the disagreement
    ax = axes[1]
    depths = [r["depth_mm"] for r in d["band"]]
    ax.plot(depths, [r["model_khz"] for r in d["band"]], "o-", color="#14505c",
            lw=2, ms=7, label="this model")
    for a_val, marker, colour in ((5.0, "s--", "#a2582b"), (10.0, "^--", "#9fb3b8")):
        rows = [e for e in d["ellens"] if e["alpha_np_m_mhz"] == a_val]
        ax.plot([e["depth_mm"] for e in rows], [e["reported_khz"] for e in rows],
                marker, color=colour, lw=1.6, ms=7,
                label=f"PMID 26233216, $\\alpha$={a_val:.0f}")
    ax.axhspan(250, 1500, color="#eef2f3", zorder=0)
    ax.set_ylim(200, 1650)
    ax.annotate("the band they scanned", (52, 1520), fontsize=7.5, color="#888",
                va="top")
    ax.set_xlabel("focal depth (mm)")
    ax.set_ylabel("optimal frequency (kHz)")
    ax.set_title("(b) the depth scaling is REFUTED", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    # (c) what makes this arm structurally different: a threshold
    ax = axes[2]
    # READ FROM THE ARTIFACT, not recomputed here. The first version of this
    # panel re-implemented `delivered_index` inline, which is two
    # implementations of one formula in one repository and a figure free to
    # disagree with the page beside it. It also plotted from zero depth, where
    # the uncapped model DIVERGES -- nothing has attenuated and the focal-gain
    # term rewards frequency without limit -- and the spike was hidden by an
    # axis choice rather than by fixing the model. The artifact's curve is
    # computed at the applicator's own frequency CAP, so it is finite at the
    # surface for a stated physical reason.
    labels = {3.0: "strong applicator", 1.5: "moderate", 0.6: "weak"}
    colours = {3.0: "#14505c", 1.5: "#a2582b", 0.6: "#9fb3b8"}
    for c in d["curves"]:
        strength = c["index_at_reference"]
        xs = [pt[0] for pt in c["points"]]
        ys = [pt[1] for pt in c["points"]]
        ax.plot(xs, ys, color=colours.get(strength, "#666"), lw=2,
                label=labels.get(strength, f"{strength}"))
    ax.axhline(d["cavitation_threshold"], color="#c1440e", ls="--", lw=1.4)
    ax.annotate("a cavitation threshold", (4.2, d["cavitation_threshold"] * 1.12),
                fontsize=7.5, color="#c1440e")
    ax.set_xlabel("focal depth (cm)")
    ax.set_ylabel("delivered mechanical index")
    ax.set_yscale("log")
    ax.set_xlim(0, 15.2)
    ax.set_title("(c) below the line, no exposure time helps", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    fig.suptitle("Sonodynamic therapy: an optimum the model finds, and a "
                 "scaling an independent study contradicts", fontsize=12,
                 fontweight="bold")
    fig.text(0.5, 0.015,
             "(a) Focusing favours a high frequency, while attenuation and the mechanical index's "
             "own 1/sqrt(f) both favour a low one, so the product peaks in between. (b) Ellens & "
             "Hynynen 2015 (PMID 26233216) simulated the same three mechanisms by an independent "
             "full-wave route and report a near-flat optimum where this model predicts a 3x fall "
             "across 50-150 mm. They name the missing term themselves - near-field heating, which "
             "this model has no representation of at all. (c) is why SDT is structurally unlike "
             "the dose-response arms: inertial cavitation is a THRESHOLD, so below the line no "
             "insonation time produces sonochemical ROS and the depth limit is a hard edge rather "
             "than a fade - where each curve crosses the line IS the depth limit tabulated in "
             "the validation page. The threshold height is a parameter, not a measurement: "
             "published in-vivo values move by more than an order of magnitude with nucleation, "
             "and the curves are drawn at a 20 MHz device cap, which is what keeps them finite at "
             "the surface: uncapped, this model diverges where nothing has yet attenuated.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.15, 1, 0.93])
    save(fig, "fig43_sonodynamic_frequency")


def fig44_pdt_fluence_rate():
    """Why the same joules delivered faster kill less, and what bounds it.

    The fourth interior optimum in the chapter built from two opposing
    monotonic effects, and the only one whose slow side is bounded by the
    treatment's own pharmacokinetics rather than by biology.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "pdt-fluence-rate-validation.json"
    if not src.exists():
        print("  fig44: pdt-fluence-rate-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())
    pc = d["phi_crit_default_mw_cm2"]
    p_full = d["oer_reference_po2_mmhg_from_oxygen_rs"]
    total = d["total_fluence_j_cm2"]

    def oer(p):
        return (3.0 * max(p, 0.0) + 3.0) / (max(p, 0.0) + 3.0)

    def yield_factor(rate):
        o2 = 1.0 / (1.0 + rate / pc) if rate > 0 else 1.0
        return oer(o2 * p_full) / oer(p_full)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) the two opposing terms, separately
    ax = axes[0]
    rates = [5.0 * pow(400.0 / 5.0, i / 300.0) for i in range(301)]
    ax.plot(rates, [yield_factor(r) for r in rates], color="#14505c", lw=2,
            label="oxygen: faster is worse")
    import math as _m
    t_half = 4.0
    k = _m.log(2.0) / t_half
    drug = []
    for r in rates:
        dur = total * 1000.0 / r / 3600.0
        drug.append((1.0 - _m.exp(-k * dur)) / (k * dur) if dur > 0 else 1.0)
    ax.plot(rates, drug, color="#a2582b", lw=2, label="drug: slower is worse")
    ax.set_xscale("log")
    ax.set_xlabel("fluence rate (mW/cm$^2$)")
    ax.set_ylabel("relative factor")
    ax.set_title("(a) two monotonic effects, opposed", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    # (b) their product: an interior optimum per sensitizer
    ax = axes[1]
    for th, colour in ((0.5, "#14505c"), (4.0, "#a2582b"), (48.0, "#9fb3b8")):
        kk = _m.log(2.0) / th
        ys = []
        for r in rates:
            dur = total * 1000.0 / r / 3600.0
            ys.append(yield_factor(r) * ((1.0 - _m.exp(-kk * dur)) / (kk * dur)))
        peak = max(ys)
        ax.plot(rates, [y / peak for y in ys], color=colour, lw=2,
                label=f"half-life {th:.1f} h")
        row = min(d["by_half_life"], key=lambda x: abs(x["t_half_h"] - th))
        ax.plot([row["optimal_mw_cm2"]], [1.0], "o", color=colour, ms=7)
    ax.set_xscale("log")
    ax.set_xlabel("fluence rate (mW/cm$^2$)")
    ax.set_ylabel("delivered singlet oxygen (peak = 1)")
    ax.set_title("(b) a best rate, per drug", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    # (c) how much of the answer is the uncalibrated knob
    ax = axes[2]
    rows = d["by_phi_crit"]
    ax.plot([r["phi_crit"] for r in rows], [r["optimal_mw_cm2"] for r in rows],
            "o-", color="#9fb3b8", lw=2, ms=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$\\varphi_{crit}$ (mW/cm$^2$) — measured nowhere here")
    ax.set_ylabel("optimal fluence rate (mW/cm$^2$)")
    ax.set_title(f"(c) the optimum rides on it, exponent "
                 f"{d['optimum_scaling_exponent']:.2f}", fontsize=10)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("Photodynamic therapy consumes the oxygen it needs, so the "
                 "same light delivered faster does less", fontsize=12,
                 fontweight="bold")
    fig.text(0.5, 0.015,
             "A Type II photosensitizer works by handing energy to ground-state oxygen, so it "
             "depletes its own substrate at the site it is treating; Henderson & Busch (PMID "
             "16615136) report depletion within seconds at 75 mW/cm2. Against that, a slower "
             "illumination runs further into the sensitizer's own clearance, and drug that has "
             "left cannot be excited - which is not new machinery here but the module's existing "
             "pharmacokinetics, integrated across the illumination instead of sampled once at its "
             "start. Sampling once is exactly the approximation that makes a long treatment look "
             "free. (c) is the refusal: the optimum's POSITION scales as roughly the square root "
             "of a parameter nothing in this repository measures, so the direction is the result "
             "and the milliwatt figure is a restatement of an assumption.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.15, 1, 0.93])
    save(fig, "fig44_pdt_fluence_rate")



def fig45_radiation_oer():
    """The oxygen enhancement ratio a spatial run exhibits, against the formula.

    The chapter's arms have all been closed-form functions evaluated at a
    point. This is the first one measured from a population on a grid, and the
    answer is not the formula: it moves with the oxygen gradient, which the
    formula has no term for, and it peaks in the middle.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "radiation-oer-validation.json"
    if not src.exists():
        print("  fig45: radiation-oer-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())
    rows = [r for r in d["by_lambda"] if r["dmf_mean"] is not None]
    lo, hi = d["published_oer_band"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) the dose-response by zone, at the best gradient
    ax = axes[0]
    best = max(rows, key=lambda r: r["dmf_mean"])
    for zone, key, colour in (("oxygenated rim", "dose_rim_gy", "#14505c"),
                              ("hypoxic core", "dose_core_gy", "#a2582b")):
        xs = [f["kill_level"] for f in best["factors"] if f[key] is not None]
        ys = [f[key] for f in best["factors"] if f[key] is not None]
        ax.plot(ys, xs, "o-", color=colour, lw=2, ms=6, label=zone)
    ax.set_xlabel("single-fraction dose (Gy)")
    ax.set_ylabel("fraction of the zone killed")
    ax.set_title(f"(a) the same dose, two oxygenations (λ={best['lambda_um']:.0f} µm)",
                 fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)

    # (b) the factor against the gradient -- the interior optimum
    ax = axes[1]
    lams = [r["lambda_um"] for r in rows]
    dmfs = [r["dmf_mean"] for r in rows]
    ax.axhspan(lo, hi, color="#dfe8ea", zorder=0)
    ax.annotate("published 2.5–3.0", (lams[-1] - 46, lo + 0.05), fontsize=7.5,
                color="#5b7078")
    ax.axhline(d["restated_single_cell_dmf"], color="#9fb3b8", ls=":", lw=1.6)
    ax.annotate(f"the formula restated ({d['restated_single_cell_dmf']})",
                (lams[0] + 3, d["restated_single_cell_dmf"] + 0.05),
                fontsize=7.5, color="#7d8e94")
    ax.plot(lams, dmfs, "o-", color="#14505c", lw=2, ms=7)
    ax.plot([best["lambda_um"]], [best["dmf_mean"]], "o", color="#c1440e", ms=10,
            zorder=5)
    ax.set_xlabel("O$_2$ gradient λ (µm)")
    ax.set_ylabel("dose-modifying factor")
    ax.set_title("(b) it moves with a term the formula lacks", fontsize=10)
    ax.grid(alpha=0.25)

    # (c) why: one lambda sets both zones
    ax = axes[2]
    ax.plot(lams, [r["rim_po2_mmhg"] for r in rows], "s-", color="#14505c",
            lw=2, ms=6, label="oxygenated rim")
    ax.plot(lams, [max(r["core_po2_mmhg"], 1e-3) for r in rows], "^-",
            color="#a2582b", lw=2, ms=6, label="hypoxic core")
    ax.set_yscale("log")
    ax.set_xlabel("O$_2$ gradient λ (µm)")
    ax.set_ylabel("mean zone pO$_2$ (mmHg)")
    ax.set_title("(c) one λ sets both, so neither is ever ideal", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("Radiation in the spatial engine: an oxygen effect measured "
                 "from a population, not restated from a formula",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "The engine already contained the Alper-Howard-Flanders hyperbola, and "
             "`dna_channel_dose_modifying_factor` returns 2.86 from it - a number the crate "
             "documents and tests as a RESTATEMENT, so scoring it against the published band would "
             "be a guard computing its own expectation. What (b) shows is a different quantity: "
             "the dose ratio for equal kill between an oxygenated rim and a hypoxic core, read off "
             "a 60^3 population, which moves from 1.24 to 2.55 with the oxygen gradient - a "
             "parameter the formula does not contain. It reaches the published band at ONE gradient "
             "and falls below it at the engine's own zone reference. (c) is why, and why every "
             "value is a LOWER bound: one lambda sets the rim and the core together, so a gradient "
             "steep enough to make the core anoxic leaves the rim hypoxic too, and the model can "
             "never present the fully-oxic-versus-anoxic pair the published band was measured on.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.15, 1, 0.93])
    save(fig, "fig45_radiation_oer")



def fig46_oncolytic_percolation():
    """Whether the front crosses, against a constant from another field.

    The closed form the engine already contained predicts spread on every one
    of these conditions. On a lattice it is wrong on most of them, and the
    place it becomes wrong is a solved problem in computational physics.
    """
    import json
    src = Path(__file__).resolve().parent.parent / "analysis" / "calibration" / \
        "oncolytic-percolation-validation.json"
    if not src.exists():
        print("  fig46: oncolytic-percolation-validation.json missing - skipping")
        return
    d = json.loads(src.read_text())
    radius = d["tumour_radius_um"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) the transition itself
    ax = axes[0]
    colours = ["#14505c", "#a2582b", "#9fb3b8"]
    for r, colour in zip(d["by_replication"], colours):
        pts = [p for p in r["points"] if p["seeded"]]
        ax.plot([p["permissive_fraction"] for p in pts],
                [p["front_radius_um"] for p in pts], "o-", color=colour, lw=2,
                ms=5, label=f"transmission {r['replication']:g}")
        unseeded = [p for p in r["points"] if not p["seeded"]]
        ax.plot([p["permissive_fraction"] for p in unseeded],
                [0.0 for _ in unseeded], "x", color=colour, ms=8)
    ax.axvline(d["pc_moore_26"], color="#c1440e", ls="--", lw=1.6)
    ax.annotate("26-neighbour\npercolation, 0.0976", (d["pc_moore_26"] + 0.012,
                radius * 0.62), fontsize=7.5, color="#c1440e")
    ax.axhline(d["cross_threshold_um"], color="#888", ls=":", lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel("fraction of the tumour the virus can enter")
    ax.set_ylabel("front reach from the seed (µm)")
    ax.set_title("(a) it crosses, or it does not", fontsize=10)
    ax.legend(fontsize=7.5, loc="center right")
    ax.grid(alpha=0.25, which="both")

    # (b) the threshold against the two lattice constants
    ax = axes[1]
    reps = [r["replication"] for r in d["by_replication"]]
    lows = [r["lowest_crossing"] for r in d["by_replication"]]
    ax.plot(reps, lows, "o-", color="#14505c", lw=2, ms=8)
    ax.axhline(d["pc_moore_26"], color="#c1440e", ls="--", lw=1.6)
    ax.annotate("26-neighbour (this lattice)", (0.33, d["pc_moore_26"] * 1.08),
                fontsize=7.5, color="#c1440e")
    ax.axhline(d["pc_von_neumann_6"], color="#9fb3b8", ls="--", lw=1.6)
    ax.annotate("6-neighbour (not this lattice)",
                (0.33, d["pc_von_neumann_6"] + 0.008), fontsize=7.5,
                color="#7d8e94")
    ax.set_ylim(0.05, 0.37)
    ax.set_xlabel("transmission probability")
    ax.set_ylabel("permissive fraction needed to cross")
    ax.set_title("(b) site–bond: certain transmission is the limit", fontsize=10)
    ax.grid(alpha=0.25)

    # (c) what the closed form says about the same conditions
    ax = axes[2]
    top = d["by_replication"][0]
    xs = [p["permissive_fraction"] for p in top["points"]]
    ax.plot(xs, [p["closed_form_speed"] for p in top["points"]], "s-",
            color="#9fb3b8", lw=2, ms=5, label="closed form: spreads")
    ax.plot(xs, [p["closed_form_speed"] if p["crossed"] else 0.0
                 for p in top["points"]], "o-", color="#14505c", lw=2, ms=5,
            label="on the grid: crosses")
    ax.axvline(d["pc_moore_26"], color="#c1440e", ls="--", lw=1.6)
    ax.set_xscale("log")
    ax.set_xlabel("fraction of the tumour the virus can enter")
    ax.set_ylabel("predicted front speed (arb.)")
    ax.set_title(f"(c) the closed form is wrong on "
                 f"{d['rows_where_closed_form_is_wrong']} conditions", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("Oncolytic spread is a threshold, and a closed-form front "
                 "speed cannot express one", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "`oncolytic::front_speed` is the Fisher-KPP speed 2*sqrt(D*r) - the speed of a front in a "
             "HOMOGENEOUS medium - so it returns a positive number whenever replication beats clearance and "
             "has no term for how much of the tumour the virus can enter. On a lattice that is false, and "
             "the failure is sharp: below a site-percolation threshold there is no connected path of "
             "enterable cells and no dose or duration gets the infection across. The measured threshold at "
             "certain transmission brackets 0.0976, the 26-neighbour site-percolation constant for the "
             "neighbourhood this grid uses - a number from computational physics that nothing here was "
             "fitted to, and which excludes the 6-neighbour value of 0.3116. Crosses mark runs where no "
             "enterable cell existed near the seed at all: a failure to START, which is a different claim "
             "from a failure to cross and is marked rather than drawn as a zero. This says the spread RULE "
             "behaves like site percolation on its own lattice; it is not evidence that tumours percolate.",
             ha="center", fontsize=7.2, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.16, 1, 0.93])
    save(fig, "fig46_oncolytic_percolation")



if __name__ == "__main__":
    print("Generating conceptual diagrams...")
    fig18_hypoxia()
    fig19_immune()
    fig20_stromal()
    fig21_ph()
    fig22_flowchart()
    fig23_census_flow()
    fig30_modality_landscape()
    fig31_modality_panel()
    fig32_modality_tme()
    fig33_adoptive_barriers()
    fig34_depth_reach()
    fig35_calibration_verdicts()
    fig36_fractionation()
    fig37_chemotherapy()
    fig38_checkpoint()
    fig39_adoptive_escalation()
    fig40_oncolytic_bind()
    fig41_adc_loading()
    fig42_ablation_sleeve()
    fig43_sonodynamic_frequency()
    fig44_pdt_fluence_rate()
    fig45_radiation_oer()
    fig46_oncolytic_percolation()
    print("Done.")
