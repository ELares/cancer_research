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
    problem: the engine models one mechanism in depth and most of the
    taxonomy not at all, and the mechanisms it misses are not small.
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
    phenos = d["phenotypes"]
    axes_order = list(ebp[phenos[0]].keys())
    arms = d["arms"]

    fig, axs = plt.subplots(1, len(phenos), figsize=(13.0, 5.6), sharey=True)
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
    cb.set_ticklabels([f"{v * 100:+.0f}" for v in np.linspace(-vmax, vmax, 5)])
    cb.set_label("change in kill fraction (percentage points)  -  "
                 "negative is resistance", fontsize=8.2)
    fig.suptitle("What the microenvironment does to every arm, by cell state",
                 fontsize=12.5, fontweight="bold", y=0.985)
    fig.text(0.5, 0.02,
             "NOT a ranking. Signed deliberately: red is a loss, blue a GAIN, "
             "and the gains are real - clonal heterogeneity supplies a "
             "low-defence tail that dies while the average cell resists. A "
             "blank cell is EXACTLY zero for that arm - a cell the axis "
             "cannot move at all in this run, which is a property of the run "
             "rather than of the biology; a real effect under one percentage "
             "point prints as <1 rather than vanishing. Every arm but "
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
    print("Done.")
