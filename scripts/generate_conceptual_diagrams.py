#!/usr/bin/env python3
"""Generate conceptual diagrams for TME mechanisms and decision flowchart.

Creates 5 figures:
  fig18_hypoxia_crosssection.pdf  — O2 gradient with drug efficacy overlay
  fig19_immune_coupling_flow.pdf  — DAMP → DC → T cell pathway
  fig20_stromal_shielding.pdf     — CAF boundary protection
  fig21_ph_ion_trapping.pdf       — pH gradient with drug trapping
  fig22_decision_flowchart.pdf    — Which modality for which context

Usage:
  python3 scripts/generate_conceptual_diagrams.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "article" / "figures"

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

    out_base = str(OUT / "fig19_immune_coupling_flow")
    gv_path = out_base + ".gv"
    with open(gv_path, "w") as f:
        f.write(dot.source)
    subprocess.run(["dot", "-Tpdf", "-o", out_base + ".pdf", gv_path], check=True)
    subprocess.run(["dot", "-Tpng", "-Gdpi=300", "-o", out_base + ".png", gv_path], check=True)
    Path(gv_path).unlink(missing_ok=True)
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
    out_base = str(OUT / "fig22_decision_flowchart")
    gv_path = out_base + ".gv"

    # Write source
    with open(gv_path, "w") as f:
        f.write(dot.source)

    # Generate PDF
    subprocess.run(["dot", "-Tpdf", "-o", out_base + ".pdf", gv_path], check=True)

    # Generate high-DPI PNG
    subprocess.run(["dot", "-Tpng", "-Gdpi=300", "-o", out_base + ".png", gv_path], check=True)

    # Clean up .gv source (reproducible from the script itself)
    Path(gv_path).unlink(missing_ok=True)

    print(f"  fig22_decision_flowchart")


if __name__ == "__main__":
    print("Generating conceptual diagrams...")
    fig18_hypoxia()
    fig19_immune()
    fig20_stromal()
    fig21_ph()
    fig22_flowchart()
    print("Done.")
