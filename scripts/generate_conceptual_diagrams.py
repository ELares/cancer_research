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

    # Annotation
    ax.annotate("RSL3 depends on\nbasal ROS (needs O₂)",
                xy=(6, 1.2), fontsize=7, color=COLORS["rsl3"],
                style="italic", ha="center")
    ax.annotate("SDT delivers\nexogenous ROS",
                xy=(6, 4.5), fontsize=7, color=COLORS["sdt"],
                style="italic", ha="center")

    save(fig, "fig18_hypoxia_crosssection")


# ── Figure 19: Immune coupling flow ───────────────────────────────────

def fig19_immune():
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_facecolor("white")

    box_style = dict(boxstyle="round,pad=0.4", fc="#E8E8E8", ec="black", lw=1.2)
    sdt_style = dict(boxstyle="round,pad=0.4", fc="#FFE0B2", ec=COLORS["sdt"], lw=1.5)
    rsl3_style = dict(boxstyle="round,pad=0.4", fc="#D1C4E9", ec=COLORS["rsl3"], lw=1.5)

    # SDT path (top)
    ax.text(1.0, 3.8, "SDT kills\n~140K cells", ha="center", fontsize=8, fontweight="bold",
            bbox=sdt_style)
    ax.text(3.2, 3.8, "High LP\novershoot\n(LP~20)", ha="center", fontsize=7, bbox=sdt_style)
    ax.text(5.3, 3.8, "Dense\nDAMP field", ha="center", fontsize=7, bbox=sdt_style)
    ax.text(7.2, 3.8, "DC\nactivation", ha="center", fontsize=7, bbox=box_style)
    ax.text(9.0, 3.8, "521 immune\nkills", ha="center", fontsize=8, fontweight="bold",
            bbox=sdt_style)

    # RSL3 path (bottom)
    ax.text(1.0, 1.2, "RSL3 kills\n~163 cells", ha="center", fontsize=8, fontweight="bold",
            bbox=rsl3_style)
    ax.text(3.2, 1.2, "Low LP\novershoot\n(LP~7.8)", ha="center", fontsize=7, bbox=rsl3_style)
    ax.text(5.3, 1.2, "Sparse\nDAMP field", ha="center", fontsize=7, bbox=rsl3_style)
    ax.text(7.2, 1.2, "Minimal DC\nactivation", ha="center", fontsize=7, bbox=box_style)
    ax.text(9.0, 1.2, "5 immune\nkills", ha="center", fontsize=8, fontweight="bold",
            bbox=rsl3_style)

    # Arrows
    for y in [3.8, 1.2]:
        for x1, x2 in [(1.7, 2.5), (3.9, 4.6), (6.0, 6.5), (7.9, 8.3)]:
            ax.annotate("", xy=(x2, y), xytext=(x1, y),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.2))

    # Ratio label
    ax.text(5.0, 2.5, "104:1 immune kill ratio",
            ha="center", fontsize=10, fontweight="bold", color=COLORS["damp"],
            bbox=dict(boxstyle="round", fc="white", ec=COLORS["damp"], lw=1.5))

    ax.set_title("Immune Coupling: Kill density determines DAMP-mediated immune activation",
                 fontsize=10, fontweight="bold")

    save(fig, "fig19_immune_coupling_flow")


# ── Figure 20: Stromal shielding ──────────────────────────────────────

def fig20_stromal():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

    for ax, title, kill_inner, kill_boundary, color, label in [
        (ax1, "RSL3 (pharmacologic)", 3.0, 1.5, COLORS["rsl3"], "Kill halved\nat boundary"),
        (ax2, "SDT (physical)", 96.1, 91.2, COLORS["sdt"], "Kill barely\naffected"),
    ]:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 8)
        ax.set_facecolor(COLORS["bg"])

        # Tumor interior
        ax.add_patch(mpatches.Rectangle((1.5, 1), 7, 6, fc="#FFCCCC", ec="black", lw=1))
        ax.text(5.0, 4.0, f"Tumor interior\nKill: {kill_inner}%", ha="center", va="center",
                fontsize=9, fontweight="bold")

        # Stromal boundary (left)
        ax.add_patch(mpatches.Rectangle((0, 1), 1.5, 6, fc=COLORS["stroma"], ec="black", lw=1))
        ax.text(0.75, 4.0, "Stroma\n(CAFs)", ha="center", va="center", fontsize=7,
                color=COLORS["caf"], fontweight="bold", rotation=90)

        # Boundary cells
        ax.add_patch(mpatches.Rectangle((1.5, 1), 1.0, 6, fc="#FFE0E0", ec="gray", lw=0.5, ls="--"))
        ax.text(2.0, 7.5, f"Boundary cells\nKill: {kill_boundary}%", ha="center",
                fontsize=8, color=color, fontweight="bold")

        # Supply arrows
        for y_pos in [2.5, 4.0, 5.5]:
            ax.annotate("", xy=(1.8, y_pos), xytext=(1.2, y_pos),
                        arrowprops=dict(arrowstyle="->", color=COLORS["caf"], lw=1.5))

        ax.text(0.3, 7.5, "GSH + MUFA\nsupply →", fontsize=7, color=COLORS["caf"])
        ax.set_title(title, fontsize=9, fontweight="bold", color=color)
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

    # Drug molecules — neutral at left, trapped at right
    for x_pos in [1.5, 3.0, 4.5]:
        ax.plot(x_pos, 3.0, "o", color=COLORS["rsl3"], markersize=10, alpha=0.8)
        ax.text(x_pos, 2.3, "RSL3", fontsize=6, ha="center", color=COLORS["rsl3"])

    # Trapped drug molecules at right (outside cells)
    for x_pos in [7.0, 8.0, 9.0]:
        ax.plot(x_pos, 4.5, "o", color=COLORS["rsl3"], markersize=10, alpha=0.4)
        ax.text(x_pos, 3.8, "RSL3⁺\n(trapped)", fontsize=5, ha="center", color=COLORS["rsl3"])

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
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_facecolor("white")

    decision = dict(boxstyle="round,pad=0.5", fc="#E3F2FD", ec="#1565C0", lw=1.5)
    outcome_yes = dict(boxstyle="round,pad=0.4", fc="#E8F5E9", ec="#2E7D32", lw=1.2)
    outcome_no = dict(boxstyle="round,pad=0.4", fc="#FFF3E0", ec="#E65100", lw=1.2)
    terminal = dict(boxstyle="round,pad=0.4", fc="#F3E5F5", ec="#6A1B9A", lw=1.2)

    # Level 1: Is tumor localizable?
    ax.text(6, 9, "Is the tumor\nlocalizable?", ha="center", fontsize=9,
            fontweight="bold", bbox=decision)

    # No → alternative approaches
    ax.text(10.5, 9, "Alternative\napproaches\n(Ch 8.1)", ha="center", fontsize=7, bbox=outcome_no)
    ax.annotate("No", xy=(9.3, 9), xytext=(7.3, 9),
                arrowprops=dict(arrowstyle="->", color="#E65100", lw=1.2),
                fontsize=8, color="#E65100", ha="center")

    # Yes → Is it deep-seated?
    ax.text(4, 7, "Is it\ndeep-seated?", ha="center", fontsize=9,
            fontweight="bold", bbox=decision)
    ax.annotate("Yes", xy=(4.5, 7.8), xytext=(5.5, 8.5),
                arrowprops=dict(arrowstyle="->", color="#2E7D32", lw=1.2),
                fontsize=8, color="#2E7D32")

    # Deep: SDT range
    ax.text(1.5, 5.5, "SDT range\n(cm depth)\nCh 6.1", ha="center", fontsize=7,
            fontweight="bold", bbox=outcome_yes)
    ax.annotate("Yes", xy=(1.8, 6.3), xytext=(3.2, 6.7),
                arrowprops=dict(arrowstyle="->", color="#2E7D32", lw=1.2),
                fontsize=8, color="#2E7D32")

    # Superficial: PDT range
    ax.text(6.5, 5.5, "PDT range\n(mm depth)\nCh 6.1", ha="center", fontsize=7, bbox=outcome_yes)
    ax.annotate("No", xy=(6.2, 6.2), xytext=(5.0, 6.7),
                arrowprops=dict(arrowstyle="->", color="#E65100", lw=1.2),
                fontsize=8, color="#E65100")

    # Level 3: Are persisters ferroptosis-prone?
    ax.text(4, 3.8, "Are residual cells\nferroptosis-prone?", ha="center", fontsize=9,
            fontweight="bold", bbox=decision)
    ax.annotate("", xy=(4, 4.6), xytext=(4, 5.1),
                arrowprops=dict(arrowstyle="->", color="gray", lw=1))

    # No → other approaches
    ax.text(8.5, 3.8, "Pathway-target\napproaches\n(Ch 8.1)", ha="center", fontsize=7,
            bbox=outcome_no)
    ax.annotate("No", xy=(7.2, 3.8), xytext=(5.5, 3.8),
                arrowprops=dict(arrowstyle="->", color="#E65100", lw=1.2),
                fontsize=8, color="#E65100")

    # Yes → Immunocompetent?
    ax.text(2.5, 2, "Immunocompetent\nsetting?", ha="center", fontsize=9,
            fontweight="bold", bbox=decision)
    ax.annotate("Yes", xy=(2.8, 2.8), xytext=(3.5, 3.4),
                arrowprops=dict(arrowstyle="->", color="#2E7D32", lw=1.2),
                fontsize=8, color="#2E7D32")

    # Yes → ICD potential
    ax.text(1, 0.5, "Physical ROS +\nanti-PD-1\n(Ch 7, 9.5)", ha="center", fontsize=7,
            fontweight="bold", bbox=terminal)
    ax.annotate("Yes", xy=(1.2, 1.0), xytext=(1.8, 1.5),
                arrowprops=dict(arrowstyle="->", color="#2E7D32", lw=1.2),
                fontsize=8, color="#2E7D32")

    # No → Direct kill only
    ax.text(5, 0.5, "Physical ROS\n(direct kill only)\nCh 6-7", ha="center", fontsize=7,
            bbox=terminal)
    ax.annotate("No", xy=(4.3, 0.8), xytext=(3.5, 1.5),
                arrowprops=dict(arrowstyle="->", color="#E65100", lw=1.2),
                fontsize=8, color="#E65100")

    ax.set_title("Decision Framework: Which Modality for Which Clinical Context?",
                 fontsize=11, fontweight="bold", pad=10)

    save(fig, "fig22_decision_flowchart")


if __name__ == "__main__":
    print("Generating conceptual diagrams...")
    fig18_hypoxia()
    fig19_immune()
    fig20_stromal()
    fig21_ph()
    fig22_flowchart()
    print("Done.")
