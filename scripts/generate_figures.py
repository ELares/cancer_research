#!/usr/bin/env python3
"""
Generate publication-quality figures and run statistical analyses.

Produces:
  article/figures/fig1_ferroptosis_comparison.pdf
  article/figures/fig2_mechanism_heatmap.pdf
  article/figures/fig3_literature_disconnect.pdf
  article/figures/fig4_molecular_overlap.pdf
  article/figures/fig5_publication_trends.pdf
  article/figures/fig6_sdt_pdt_depth.pdf
  article/figures/fig13_gold_set_eval.pdf
  article/figures/fig14_tissue_mechanism_heatmap.pdf
  article/figures/fig15_designed_combinations.pdf
  article/figures/fig16_weighted_evidence.pdf
  article/figures/fig17_damp_heatmap.pdf
"""

import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import yaml
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PMID_DIR = PROJECT_ROOT / "corpus" / "by-pmid"
INDEX_FILE = PROJECT_ROOT / "corpus" / "INDEX.jsonl"
FIG_DIR = PROJECT_ROOT / "article" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 2D TME (sim-tme) and combination-mechanism (sim-combo-mech) outputs. Both are
# gitignored; regenerate with `cargo run --release -p sim-tme` / `-p sim-combo-mech`.
TME_SUMMARY = PROJECT_ROOT / "simulations" / "output" / "tme" / "tme_summary.json"
COMBO_SUMMARY = PROJECT_ROOT / "simulations" / "output" / "combo-mech" / "combo_summary.json"
WINDOW_JSON = PROJECT_ROOT / "simulations" / "output" / "window" / "vulnerability_window.json"
# Depth-kill curves (sim-spatial, 2D). Gitignored; regenerate with
# `cargo run --release -p sim-spatial`.
SPATIAL_CURVES = PROJECT_ROOT / "simulations" / "output" / "spatial" / "depth_kill_curves.csv"

plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


def load_corpus():
    """Load all articles with frontmatter and abstract text."""
    articles = []
    for f in sorted(PMID_DIR.glob("*.md")):
        try:
            content = f.read_text(encoding="utf-8")
            match = re.match(r"^---\n(.*?\n)---\n\n?(.*)", content, re.DOTALL)
            if not match:
                continue
            fm = yaml.safe_load(match.group(1)) or {}
            body = match.group(2).lower()
            text = (fm.get("title", "") + " " + body[:4000]).lower()
            fm["_text"] = text
            articles.append(fm)
        except Exception:
            pass
    return articles


def classify_ferroptosis(articles):
    """Classify articles by physical modality and ferroptosis/ICD engagement."""
    physical_mechs = {
        "SDT": "sonodynamic",
        "TTFields": "ttfields",
        "HIFU": "hifu",
        "IRE": "electrochemical-therapy",
        "Frequency": "frequency-therapy",
    }

    results = {}
    for label, mech_key in physical_mechs.items():
        mech_articles = [a for a in articles if mech_key in a.get("mechanisms", [])]
        ferro = sum(1 for a in mech_articles if "ferroptosis" in a.get("_text", ""))
        icd = sum(1 for a in mech_articles if re.search(
            r"immunogenic cell death|calreticulin|HMGB1", a.get("_text", ""), re.IGNORECASE
        ))
        gsh = sum(1 for a in mech_articles if re.search(
            r"glutathione|GSH|GPX4|SLC7A11", a.get("_text", ""), re.IGNORECASE
        ))
        both = sum(1 for a in mech_articles if (
            "ferroptosis" in a.get("_text", "") and
            re.search(r"immunogenic cell death|calreticulin|HMGB1", a.get("_text", ""), re.IGNORECASE)
        ))
        results[label] = {
            "total": len(mech_articles),
            "ferroptosis": ferro,
            "icd": icd,
            "gsh": gsh,
            "both": both,
            "ferro_pct": 100 * ferro / len(mech_articles) if mech_articles else 0,
            "icd_pct": 100 * icd / len(mech_articles) if mech_articles else 0,
            "gsh_pct": 100 * gsh / len(mech_articles) if mech_articles else 0,
        }
    return results


# ============================================================
# FIGURE 1: Ferroptosis Engagement Comparison (Bar Chart + Stats)
# ============================================================

def fig1_ferroptosis_comparison(articles):
    print("Generating Figure 1: Ferroptosis engagement comparison...")
    data = classify_ferroptosis(articles)

    modalities = ["SDT", "IRE", "HIFU", "TTFields", "Frequency"]
    ferro_pcts = [data[m]["ferro_pct"] for m in modalities]
    icd_pcts = [data[m]["icd_pct"] for m in modalities]
    gsh_pcts = [data[m]["gsh_pct"] for m in modalities]

    x = np.arange(len(modalities))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width, ferro_pcts, width, label="Ferroptosis", color="#d62728", alpha=0.85)
    bars2 = ax.bar(x, icd_pcts, width, label="ICD markers", color="#1f77b4", alpha=0.85)
    bars3 = ax.bar(x + width, gsh_pcts, width, label="GSH/GPX4 axis", color="#2ca02c", alpha=0.85)

    ax.set_xlabel("Physical Modality")
    ax.set_ylabel("% of Modality Articles Engaging Pathway")
    ax.set_title("Ferroptosis and ICD Pathway Engagement by Physical Modality\n(% of articles within each modality's corpus)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{m}\n(n={data[m]['total']})" for m in modalities])
    ax.legend()
    ax.set_ylim(0, max(max(ferro_pcts), max(icd_pcts), max(gsh_pcts)) * 1.2)

    # Add value labels
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            if height > 0.5:
                ax.annotate(f'{height:.1f}%', xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8)

    # Statistical test: chi-squared for SDT vs others on ferroptosis
    sdt_ferro = data["SDT"]["ferroptosis"]
    sdt_total = data["SDT"]["total"]
    other_ferro = sum(data[m]["ferroptosis"] for m in modalities if m != "SDT")
    other_total = sum(data[m]["total"] for m in modalities if m != "SDT")

    contingency = [[sdt_ferro, sdt_total - sdt_ferro],
                   [other_ferro, other_total - other_ferro]]
    chi2, p_value, _, _ = stats.chi2_contingency(contingency)

    ax.text(0.98, 0.95, f"SDT vs others (ferroptosis):\nχ² = {chi2:.1f}, p < {p_value:.1e}",
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig1_ferroptosis_comparison.pdf")
    fig.savefig(FIG_DIR / "fig1_ferroptosis_comparison.png")
    plt.close()
    print(f"  χ² = {chi2:.1f}, p = {p_value:.2e}")
    print(f"  SDT ferroptosis: {sdt_ferro}/{sdt_total} ({data['SDT']['ferro_pct']:.1f}%)")


# ============================================================
# FIGURE 2: Mechanism-Cancer Heatmap
# ============================================================

def fig2_mechanism_heatmap(articles):
    print("Generating Figure 2: Mechanism-cancer heatmap...")

    from config import MECHANISM_KEYWORDS, CANCER_TYPE_KEYWORDS
    mechanisms = sorted(MECHANISM_KEYWORDS.keys())
    cancers = sorted(CANCER_TYPE_KEYWORDS.keys())

    matrix = np.zeros((len(mechanisms), len(cancers)))
    for a in articles:
        for i, m in enumerate(mechanisms):
            if m in a.get("mechanisms", []):
                for j, c in enumerate(cancers):
                    if c in a.get("cancer_types", []):
                        matrix[i, j] += 1

    # Sort by total
    mech_totals = matrix.sum(axis=1)
    cancer_totals = matrix.sum(axis=0)
    mech_order = np.argsort(-mech_totals)
    cancer_order = np.argsort(-cancer_totals)

    matrix = matrix[mech_order][:, cancer_order]
    mechanisms = [mechanisms[i] for i in mech_order]
    cancers = [cancers[i] for i in cancer_order]

    fig, ax = plt.subplots(figsize=(16, 10))
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="auto")

    ax.set_xticks(range(len(cancers)))
    ax.set_xticklabels([c[:10] for c in cancers], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(mechanisms)))
    ax.set_yticklabels(mechanisms, fontsize=8)
    ax.set_title("Mechanism × Cancer Type Article Count Matrix")

    plt.colorbar(im, ax=ax, label="Article Count", shrink=0.8)

    # Mark zeros
    for i in range(len(mechanisms)):
        for j in range(len(cancers)):
            if matrix[i, j] == 0:
                ax.plot(j, i, 'x', color='gray', markersize=4, alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig2_mechanism_heatmap.pdf")
    fig.savefig(FIG_DIR / "fig2_mechanism_heatmap.png")
    plt.close()
    print(f"  Matrix: {len(mechanisms)}×{len(cancers)}, {int(matrix.sum())} total instances")


# ============================================================
# FIGURE 3: Literature Disconnect Analysis
# ============================================================

def fig3_literature_disconnect(articles):
    print("Generating Figure 3: Literature disconnect (persister vs SDT)...")

    # Classify articles into communities
    persister_articles = set()
    sdt_ferro_articles = set()
    pdt_ferro_articles = set()
    overlap = set()

    for a in articles:
        pmid = a.get("pmid", "")
        text = a.get("_text", "")
        mechs = a.get("mechanisms", [])

        is_persister = bool(re.search(r"persister|drug.tolerant|minimal residual", text))
        is_ferro = "ferroptosis" in text
        is_sdt = "sonodynamic" in mechs
        is_pdt = bool(re.search(r"photodynamic", text))

        if is_persister and is_ferro:
            persister_articles.add(pmid)
        if is_sdt and is_ferro:
            sdt_ferro_articles.add(pmid)
        if is_pdt and is_ferro:
            pdt_ferro_articles.add(pmid)

    overlap_sdt = persister_articles & sdt_ferro_articles
    overlap_pdt = persister_articles & pdt_ferro_articles

    # Venn-like visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Community sizes
    communities = ["Persister-\nFerroptosis", "SDT-\nFerroptosis", "PDT-\nFerroptosis"]
    sizes = [len(persister_articles), len(sdt_ferro_articles), len(pdt_ferro_articles)]
    overlaps = ["-", str(len(overlap_sdt)), str(len(overlap_pdt))]
    colors = ["#ff7f0e", "#1f77b4", "#2ca02c"]

    bars = ax1.bar(communities, sizes, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.set_ylabel("Number of Articles in Corpus")
    ax1.set_title("Research Community Sizes in Our Corpus")

    for bar, size in zip(bars, sizes):
        ax1.annotate(str(size), xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 5), textcoords="offset points", ha='center', fontsize=12, fontweight='bold')

    # Right: Overlap matrix
    labels = ["Persister-Ferro", "SDT-Ferro", "PDT-Ferro"]
    overlap_matrix = np.array([
        [len(persister_articles), len(overlap_sdt), len(overlap_pdt)],
        [len(overlap_sdt), len(sdt_ferro_articles), len(sdt_ferro_articles & pdt_ferro_articles)],
        [len(overlap_pdt), len(sdt_ferro_articles & pdt_ferro_articles), len(pdt_ferro_articles)],
    ])

    im = ax2.imshow(overlap_matrix, cmap="Blues")
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax2.set_yticks(range(3))
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.set_title("Cross-Community Article Overlap")

    for i in range(3):
        for j in range(3):
            ax2.text(j, i, str(overlap_matrix[i, j]),
                    ha="center", va="center", fontsize=14, fontweight='bold',
                    color="white" if overlap_matrix[i, j] > overlap_matrix.max() * 0.5 else "black")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig3_literature_disconnect.pdf")
    fig.savefig(FIG_DIR / "fig3_literature_disconnect.png")
    plt.close()
    print(f"  Persister-ferro: {len(persister_articles)}, SDT-ferro: {len(sdt_ferro_articles)}, Overlap: {len(overlap_sdt)}")


# ============================================================
# FIGURE 4: Molecular Pathway Overlap Across Physical Modalities
# ============================================================

def fig4_molecular_overlap(articles):
    print("Generating Figure 4: Molecular pathway overlap...")

    physical_mechs = {
        "SDT": "sonodynamic",
        "TTFields": "ttfields",
        "IRE": "electrochemical-therapy",
        "HIFU": "hifu",
    }

    pathways = {
        "Ferroptosis": r"ferroptosis",
        "ICD/DAMPs": r"immunogenic cell death|calreticulin|HMGB1|DAMP",
        "GSH/GPX4": r"glutathione|GSH|GPX4",
        "STING/cGAS": r"STING|cGAS|sting pathway",
        "ROS": r"reactive oxygen species|ROS generation",
        "Apoptosis": r"apoptosis|caspase",
        "ER Stress": r"endoplasmic reticulum stress|ER stress|UPR",
        "Autophagy": r"autophagy|autophagic",
    }

    matrix = np.zeros((len(physical_mechs), len(pathways)))

    for a in articles:
        text = a.get("_text", "")
        mechs = set(a.get("mechanisms", []))

        for i, (label, mech_key) in enumerate(physical_mechs.items()):
            if mech_key not in mechs:
                continue
            for j, (pw_name, pw_pattern) in enumerate(pathways.items()):
                if re.search(pw_pattern, text, re.IGNORECASE):
                    matrix[i, j] += 1

    # Normalize by modality total
    totals = matrix.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    pct_matrix = 100 * matrix / np.array([[
        sum(1 for a in articles if mk in a.get("mechanisms", []))
        for mk in physical_mechs.values()
    ]]).T

    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(pct_matrix, cmap="RdYlBu_r", aspect="auto", vmin=0)

    ax.set_xticks(range(len(pathways)))
    ax.set_xticklabels(list(pathways.keys()), rotation=30, ha="right")
    ax.set_yticks(range(len(physical_mechs)))
    ax.set_yticklabels(list(physical_mechs.keys()))
    ax.set_title("Molecular Pathway Engagement by Physical Modality\n(% of articles within each modality)")

    for i in range(len(physical_mechs)):
        for j in range(len(pathways)):
            val = pct_matrix[i, j]
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=9,
                   color="white" if val > 30 else "black")

    plt.colorbar(im, ax=ax, label="% of modality articles", shrink=0.8)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig4_molecular_overlap.pdf")
    fig.savefig(FIG_DIR / "fig4_molecular_overlap.png")
    plt.close()
    print("  Done")


# ============================================================
# FIGURE 5: Publication Trends
# ============================================================

def fig5_publication_trends(articles):
    print("Generating Figure 5: Publication trends...")

    year_mech = defaultdict(Counter)
    for a in articles:
        year = a.get("year")
        if not year or year < 2015:
            continue
        for m in a.get("mechanisms", []):
            year_mech[year][m] += 1

    years = sorted(y for y in year_mech.keys() if 2015 <= y <= 2025)

    highlight_mechs = ["immunotherapy", "nanoparticle", "sonodynamic", "car-t", "ttfields"]
    colors = ["#1f77b4", "#ff7f0e", "#d62728", "#2ca02c", "#9467bd"]

    fig, ax = plt.subplots(figsize=(12, 6))
    for mech, color in zip(highlight_mechs, colors):
        counts = [year_mech[y].get(mech, 0) for y in years]
        ax.plot(years, counts, '-o', label=mech.replace("-", " ").title(),
                color=color, linewidth=2, markersize=4)

    ax.set_xlabel("Year")
    ax.set_ylabel("Number of Articles")
    ax.set_title("Publication Volume by Mechanism (2015-2025)")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig5_publication_trends.pdf")
    fig.savefig(FIG_DIR / "fig5_publication_trends.png")
    plt.close()
    print("  Done")


# ============================================================
# FIGURE 6: SDT Ferroptosis Chain — Quantified Evidence
# ============================================================

def fig6_sdt_chain_evidence(articles):
    print("Generating Figure 6: SDT ferroptosis-ICD chain evidence...")

    sdt_articles = [a for a in articles if "sonodynamic" in a.get("mechanisms", [])]

    chain_steps = {
        "ROS\ngeneration": r"reactive oxygen species|ROS",
        "GSH\ndepletion": r"glutathione|GSH depletion|GSH consumption",
        "GPX4\ninactivation": r"GPX4|glutathione peroxidase 4",
        "Lipid\nperoxidation": r"lipid peroxid",
        "Ferroptosis": r"ferroptosis",
        "DAMP\nrelease": r"calreticulin|HMGB1|DAMP|damage.associated",
        "STING\nactivation": r"STING|cGAS",
        "ICD": r"immunogenic cell death",
    }

    counts = []
    for step_name, pattern in chain_steps.items():
        n = sum(1 for a in sdt_articles if re.search(pattern, a.get("_text", ""), re.IGNORECASE))
        counts.append(n)

    fig, ax = plt.subplots(figsize=(14, 5))

    x = range(len(chain_steps))
    bars = ax.bar(x, counts, color=plt.cm.Reds(np.linspace(0.3, 0.9, len(counts))),
                  edgecolor='black', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(list(chain_steps.keys()), fontsize=9)
    ax.set_ylabel(f"Number of SDT Articles (n={len(sdt_articles)} total)")
    ax.set_title("Evidence Depth Along the SDT → Ferroptosis → ICD Chain\n"
                 "(Number of SDT articles engaging each step)")

    # Add arrows between bars
    for i in range(len(counts) - 1):
        ax.annotate("→", xy=(i + 0.5, max(counts[i], counts[i+1]) + 5),
                    fontsize=16, ha='center', color='gray')

    # Add count labels
    for bar, count in zip(bars, counts):
        ax.annotate(str(count), xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                   xytext=(0, 3), textcoords="offset points", ha='center', fontsize=10, fontweight='bold')

    # Add the chain narrative
    ax.text(0.5, -0.18, "Each bar = number of SDT articles mentioning this step. "
            "The chain is supported at each link but thins from left to right.",
            transform=ax.transAxes, ha='center', fontsize=9, style='italic', color='gray')

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig6_sdt_chain_evidence.pdf")
    fig.savefig(FIG_DIR / "fig6_sdt_chain_evidence.png")
    plt.close()
    print(f"  Chain: {' → '.join(str(c) for c in counts)}")


# ============================================================
# Fig 9: Evidence Tier Composition
# ============================================================

EVIDENCE_ORDER = [
    "phase3-clinical", "phase2-clinical", "phase1-clinical", "clinical-other",
    "preclinical-invivo", "preclinical-invitro", "theoretical",
]
EVIDENCE_LABELS = {
    "phase3-clinical": "Phase III",
    "phase2-clinical": "Phase II",
    "phase1-clinical": "Phase I",
    "clinical-other": "Clinical (non-phase)",
    "preclinical-invivo": "Preclinical in vivo",
    "preclinical-invitro": "Preclinical in vitro",
    "theoretical": "Theoretical",
}
EVIDENCE_COLORS = {
    "phase3-clinical": "#b71c1c",
    "phase2-clinical": "#e65100",
    "phase1-clinical": "#f9a825",
    "clinical-other": "#4fc3f7",
    "preclinical-invivo": "#388e3c",
    "preclinical-invitro": "#81c784",
    "theoretical": "#bdbdbd",
}

TIER_RANK = {lvl: i for i, lvl in enumerate(reversed(EVIDENCE_ORDER))}


def load_index():
    """Load INDEX.jsonl as a list of dicts."""
    entries = []
    with open(INDEX_FILE) as f:
        for line in f:
            entries.append(json.loads(line))
    return entries


def fig9_evidence_tiers(index):
    """Stacked horizontal bar: evidence tier composition per mechanism."""
    print("Figure 9: Evidence tier composition...")

    mech_tiers = defaultdict(lambda: Counter())
    for e in index:
        ev = e.get("evidence_level", "")
        if not ev:
            continue
        for m in e.get("mechanisms", []):
            mech_tiers[m][ev] += 1

    if not mech_tiers:
        print("  No evidence data — skipping")
        return

    def highest_tier(counts):
        for lvl in EVIDENCE_ORDER:
            if counts.get(lvl, 0) > 0:
                return TIER_RANK[lvl]
        return -1

    mechs = sorted(mech_tiers.keys(), key=lambda m: highest_tier(mech_tiers[m]))

    fig, ax = plt.subplots(figsize=(10, max(6, len(mechs) * 0.35)))

    y_pos = np.arange(len(mechs))
    lefts = np.zeros(len(mechs))

    for lvl in EVIDENCE_ORDER:
        widths = [mech_tiers[m].get(lvl, 0) for m in mechs]
        ax.barh(y_pos, widths, left=lefts, height=0.7,
                color=EVIDENCE_COLORS[lvl], label=EVIDENCE_LABELS[lvl])
        lefts += widths

    ax.set_yticks(y_pos)
    ax.set_yticklabels(mechs, fontsize=9)
    ax.set_xlabel("Number of articles")
    ax.set_title("Evidence Tier Composition by Mechanism")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.text(0.5, -0.08,
            "Only articles with detected evidence tags shown. "
            "Untagged and review-like articles excluded.",
            transform=ax.transAxes, ha='center', fontsize=8, style='italic', color='gray')

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig9_evidence_tiers.pdf")
    fig.savefig(FIG_DIR / "fig9_evidence_tiers.png")
    plt.close()
    print(f"  {len(mechs)} mechanisms plotted")


# ============================================================
# Fig 10: In-Vivo vs 2D Ferroptosis Comparison
# ============================================================

INVIVO_JSON = PROJECT_ROOT / "simulations" / "output" / "invivo" / "invivo_comparison.json"


def fig10_invivo_comparison():
    """Grouped bar: 2D vs in-vivo vs SCD1i for RSL3 and exogenous ROS."""
    print("Figure 10: In-vivo ferroptosis comparison...")

    if not INVIVO_JSON.exists():
        print(f"  {INVIVO_JSON} not found — run sim-invivo first. Skipping.")
        return

    data = json.loads(INVIVO_JSON.read_text())

    phenotype_order = ["Glycolytic", "OXPHOS", "Persister (FSP1↓)", "Persister+NRF2"]
    context_order = ["2d", "invivo", "invivo+scd1i"]
    context_labels = {"2d": "2D culture", "invivo": "In vivo", "invivo+scd1i": "In vivo + SCD1i"}
    context_colors = {"2d": "#1565c0", "invivo": "#c62828", "invivo+scd1i": "#f9a825"}

    # Use SDT rows as the exogenous ROS representative (SDT = PDT in this binary)
    treatments = [("RSL3", "RSL3 (GPX4 inhibitor)"), ("SDT", "Exogenous ROS")]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax_idx, (tx_key, tx_label) in enumerate(treatments):
        ax = axes[ax_idx]
        x = np.arange(len(phenotype_order))
        width = 0.25

        for ci, ctx in enumerate(context_order):
            rates = []
            ci_lo = []
            ci_hi = []
            for pheno in phenotype_order:
                rec = next((r for r in data
                           if r["context"] == ctx and r["phenotype"] == pheno and r["treatment"] == tx_key), None)
                if rec:
                    rates.append(rec["death_rate"] * 100)
                    ci_lo.append((rec["death_rate"] - rec["ci_low"]) * 100)
                    ci_hi.append((rec["ci_high"] - rec["death_rate"]) * 100)
                else:
                    rates.append(0)
                    ci_lo.append(0)
                    ci_hi.append(0)

            offset = (ci - 1) * width
            bars = ax.bar(x + offset, rates, width, yerr=[ci_lo, ci_hi],
                         color=context_colors[ctx], label=context_labels[ctx],
                         capsize=3, error_kw={"lw": 0.8})

        ax.set_xlabel("Cell phenotype")
        ax.set_xticks(x)
        ax.set_xticklabels(phenotype_order, fontsize=9, rotation=15, ha='right')
        ax.set_title(tx_label)
        ax.set_ylim(0, 105)

    axes[0].set_ylabel("Death rate (%)")
    axes[0].legend(fontsize=8, loc="upper left")

    fig.suptitle("Effect of SCD1/MUFA Lipid Remodeling on Ferroptosis Sensitivity", fontsize=13, y=1.02)
    fig.text(0.5, -0.04,
             "SDT and PDT are modeled identically (shared exogenous ROS). "
             "In-vivo cells start at MUFA steady state (40% PUFA displacement).",
             ha='center', fontsize=8, style='italic', color='gray')

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig10_invivo_comparison.pdf")
    fig.savefig(FIG_DIR / "fig10_invivo_comparison.png")
    plt.close()
    print("  2 panels (RSL3 + Exo. ROS)")


# ============================================================
# Fig 11: MUFA Sweep Heatmaps
# ============================================================

SWEEP_CSV = PROJECT_ROOT / "simulations" / "output" / "invivo" / "mufa_sweep.csv"


def fig11_mufa_sweep():
    """Side-by-side heatmaps: RSL3 vs exogenous ROS death rate across MUFA parameter space."""
    print("Figure 11: MUFA parameter sweep heatmaps...")

    if not SWEEP_CSV.exists():
        print(f"  {SWEEP_CSV} not found — run sim-invivo first. Skipping.")
        return

    rows = []
    with open(SWEEP_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    # Filter to steady-state only (initial_mufa_protection > 0)
    steady = [r for r in rows if float(r["initial_mufa_protection"]) > 0]

    if not steady:
        print("  No steady-state sweep data found — skipping")
        return

    # SDT rows represent exogenous ROS
    treatments = [("RSL3", "RSL3 (GPX4 inhibitor)"), ("SDT", "Exogenous ROS")]

    rates_sorted = sorted(set(float(r["scd_mufa_rate"]) for r in steady))
    maxes_sorted = sorted(set(float(r["scd_mufa_max"]) for r in steady))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    for ax_idx, (tx_key, tx_label) in enumerate(treatments):
        ax = axes[ax_idx]
        tx_rows = [r for r in steady if r["treatment"] == tx_key]

        grid = np.zeros((len(rates_sorted), len(maxes_sorted)))
        for r in tx_rows:
            ri = rates_sorted.index(float(r["scd_mufa_rate"]))
            mi = maxes_sorted.index(float(r["scd_mufa_max"]))
            grid[ri, mi] = float(r["death_rate"]) * 100

        im = ax.imshow(grid, cmap="RdYlGn_r", vmin=0, vmax=100, aspect="auto",
                       origin="lower")
        ax.set_xticks(range(len(maxes_sorted)))
        ax.set_xticklabels([f"{m:.2f}" for m in maxes_sorted], fontsize=8)
        ax.set_yticks(range(len(rates_sorted)))
        ax.set_yticklabels([f"{r:.3f}" for r in rates_sorted], fontsize=8)
        ax.set_xlabel("scd_mufa_max")
        ax.set_ylabel("scd_mufa_rate")
        ax.set_title(tx_label)

        for i in range(len(rates_sorted)):
            for j in range(len(maxes_sorted)):
                val = grid[i, j]
                color = "white" if val > 50 else "black"
                ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                       fontsize=7, color=color)

    fig.colorbar(im, ax=axes, shrink=0.8, label="Persister death rate (%)")
    fig.suptitle("MUFA Parameter Sensitivity — Persister Cells (steady-state)", fontsize=13, y=1.02)
    fig.text(0.5, -0.04,
             "Cells start at analytically computed MUFA steady state. "
             "Decay rate fixed at 0.005 across all points.",
             ha='center', fontsize=8, style='italic', color='gray')

    fig.savefig(FIG_DIR / "fig11_mufa_sweep.pdf", bbox_inches='tight')
    fig.savefig(FIG_DIR / "fig11_mufa_sweep.png", bbox_inches='tight')
    plt.close()
    print(f"  {len(rates_sorted)}×{len(maxes_sorted)} grid, 2 panels")


# ============================================================
# Fig 12: Pathway Target Prevalence
# ============================================================

def _is_review_or_protocol(entry: dict) -> bool:
    """Title-based approximation of the canonical review/protocol classification.

    INDEX.jsonl lacks pub_types, so this uses the same title markers as
    evidence_utils.is_review_like / is_protocol_like. This matches the
    primary-study-like definition in analyze_corpus.py (tagged + other_untagged).
    """
    title = entry.get("title", "").lower()
    review_markers = ("review", "systematic review", "meta-analysis", "meta analysis",
                      "scoping review", "narrative review", "evidence map")
    protocol_markers = ("protocol", "study protocol", "trial protocol", "protocol for")
    return any(m in title for m in review_markers) or any(m in title for m in protocol_markers)


def fig12_pathway_targets(index):
    """Horizontal bar: pathway target prevalence, total vs primary-study-like."""
    print("Figure 12: Pathway target prevalence...")

    target_total = Counter()
    target_primary = Counter()

    for e in index:
        for pt in e.get("pathway_targets", []):
            target_total[pt] += 1
            if not _is_review_or_protocol(e):
                target_primary[pt] += 1

    if not target_total:
        print("  No pathway target data — skipping")
        return

    targets = sorted(target_total.keys(), key=lambda t: target_total[t])
    totals = [target_total[t] for t in targets]
    primaries = [target_primary.get(t, 0) for t in targets]

    fig, ax = plt.subplots(figsize=(9, max(4, len(targets) * 0.45)))
    y = np.arange(len(targets))
    height = 0.35

    ax.barh(y + height / 2, totals, height, color="#90a4ae", label="Total (incl. reviews)")
    ax.barh(y - height / 2, primaries, height, color="#1565c0", label="Primary-study-like")

    ax.set_yticks(y)
    ax.set_yticklabels(targets, fontsize=9)
    ax.set_xlabel("Number of articles")
    ax.set_title("Pathway Target Prevalence in Corpus")
    ax.legend(fontsize=9, loc="lower right")
    ax.text(0.5, -0.1,
            "Primary-study-like = non-review, non-protocol articles (title-based classification). "
            "cuproptosis-core and fdx1-cuproptosis-axis overlap ~100%.",
            transform=ax.transAxes, ha='center', fontsize=8, style='italic', color='gray')

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig12_pathway_targets.pdf")
    fig.savefig(FIG_DIR / "fig12_pathway_targets.png")
    plt.close()
    print(f"  {len(targets)} targets plotted")


# ============================================================
# Fig 13: Gold-Set Evaluation
# ============================================================

GOLD_SET_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-set-v1.csv"
GOLD_LABELS_FILE = PROJECT_ROOT / "analysis" / "evidence-gold-labels-v1.csv"

GOLD_LABEL_ORDER = [
    "phase3-clinical", "phase2-clinical", "phase1-clinical", "clinical-other",
    "preclinical-invivo", "preclinical-invitro", "theoretical", "none-applicable",
]


def fig13_gold_set_eval():
    """Grouped horizontal bar: precision and recall per evidence label from gold-set evaluation."""
    print("Figure 13: Gold-set evaluation...")

    if not GOLD_SET_FILE.exists() or not GOLD_LABELS_FILE.exists():
        print("  Gold-set CSV files not found — skipping")
        return

    with open(GOLD_SET_FILE, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    with open(GOLD_LABELS_FILE, newline="", encoding="utf-8") as f:
        labels = {r["pmid"]: r for r in csv.DictReader(f)}

    labeled = []
    for row in rows:
        label_row = labels.get(row["pmid"])
        if not label_row:
            continue
        predicted = row["predicted_evidence_level"] or "none-applicable"
        gold = label_row["gold_evidence_level"]
        labeled.append({"predicted": predicted, "gold": gold})

    if not labeled:
        print("  No labeled rows — skipping")
        return

    per_label = {}
    for label in GOLD_LABEL_ORDER:
        tp = sum(1 for r in labeled if r["gold"] == label and r["predicted"] == label)
        fp = sum(1 for r in labeled if r["gold"] != label and r["predicted"] == label)
        fn = sum(1 for r in labeled if r["gold"] == label and r["predicted"] != label)
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        per_label[label] = {"precision": precision, "recall": recall, "tp": tp, "fp": fp, "fn": fn}

    # Binary metrics
    tp_bin = sum(1 for r in labeled if r["gold"] != "none-applicable" and r["predicted"] != "none-applicable")
    fp_bin = sum(1 for r in labeled if r["gold"] == "none-applicable" and r["predicted"] != "none-applicable")
    fn_bin = sum(1 for r in labeled if r["gold"] != "none-applicable" and r["predicted"] == "none-applicable")
    exact = sum(1 for r in labeled if r["predicted"] == r["gold"])
    bin_prec = tp_bin / (tp_bin + fp_bin) if (tp_bin + fp_bin) else 0.0
    bin_rec = tp_bin / (tp_bin + fn_bin) if (tp_bin + fn_bin) else 0.0
    bin_f1 = 2 * bin_prec * bin_rec / (bin_prec + bin_rec) if (bin_prec + bin_rec) else 0.0

    display_labels = [EVIDENCE_LABELS.get(l, l.replace("-", " ").title()) for l in GOLD_LABEL_ORDER]
    prec_vals = [per_label[l]["precision"] for l in GOLD_LABEL_ORDER]
    rec_vals = [per_label[l]["recall"] for l in GOLD_LABEL_ORDER]

    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(GOLD_LABEL_ORDER))
    height = 0.35

    ax.barh(y + height / 2, prec_vals, height, color="#1565c0", label="Precision")
    ax.barh(y - height / 2, rec_vals, height, color="#e65100", label="Recall")

    ax.set_yticks(y)
    ax.set_yticklabels(display_labels, fontsize=10)
    ax.set_xlabel("Score")
    ax.set_xlim(0, 1.15)
    ax.set_title("Evidence Tagger Performance: Gold-Set Evaluation\n(100-article stratified sample)")
    ax.legend(fontsize=10, loc="lower right")

    # Annotation box with overall metrics
    summary = (
        f"Overall (n={len(labeled)})\n"
        f"Exact accuracy: {exact}/{len(labeled)} ({exact/len(labeled):.0%})\n"
        f"Binary precision: {bin_prec:.0%}\n"
        f"Binary recall: {bin_rec:.0%}\n"
        f"Binary F1: {bin_f1:.3f}"
    )
    ax.text(0.98, 0.98, summary, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    for i, (p, r) in enumerate(zip(prec_vals, rec_vals)):
        if p > 0:
            ax.text(p + 0.02, i + height / 2, f"{p:.0%}", va='center', fontsize=8)
        if r > 0:
            ax.text(r + 0.02, i - height / 2, f"{r:.0%}", va='center', fontsize=8)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig13_gold_set_eval.pdf")
    fig.savefig(FIG_DIR / "fig13_gold_set_eval.png")
    plt.close()
    print(f"  {len(labeled)} labeled rows, exact accuracy {exact/len(labeled):.0%}")


# ============================================================
# Fig 14: Tissue × Mechanism Heatmap
# ============================================================

TISSUE_ORDER = ["epithelial", "hematologic", "mesenchymal", "neuroectodermal", "mesothelial"]
TISSUE_DISPLAY = {
    "epithelial": "Epithelial",
    "hematologic": "Hematologic",
    "mesenchymal": "Mesenchymal",
    "neuroectodermal": "Neuroectodermal",
    "mesothelial": "Mesothelial",
}


def fig14_tissue_mechanism_heatmap(index):
    """Heatmap: tissue-of-origin × mechanism article counts."""
    print("Figure 14: Tissue × mechanism heatmap...")

    from config import MECHANISM_KEYWORDS
    mechanisms = sorted(MECHANISM_KEYWORDS.keys())

    matrix = defaultdict(Counter)
    tissue_totals = Counter()
    assigned = 0

    for e in index:
        tissues = e.get("tissue_categories", [])
        if not tissues:
            continue
        assigned += 1
        for tissue in tissues:
            tissue_totals[tissue] += 1
            for mech in e.get("mechanisms", []):
                matrix[tissue][mech] += 1

    data = np.zeros((len(TISSUE_ORDER), len(mechanisms)))
    for i, tissue in enumerate(TISSUE_ORDER):
        for j, mech in enumerate(mechanisms):
            data[i, j] = matrix[tissue].get(mech, 0)

    fig, ax = plt.subplots(figsize=(16, 5))
    im = ax.imshow(data, cmap="YlOrRd", aspect="auto", interpolation="nearest")

    ax.set_xticks(np.arange(len(mechanisms)))
    ax.set_xticklabels([m.replace("-", "\n") for m in mechanisms], fontsize=7, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(TISSUE_ORDER)))
    ax.set_yticklabels(
        [f"{TISSUE_DISPLAY[t]} ({tissue_totals[t]})" for t in TISSUE_ORDER],
        fontsize=10,
    )

    # Cell value annotations
    for i in range(len(TISSUE_ORDER)):
        for j in range(len(mechanisms)):
            val = int(data[i, j])
            if val > 0:
                color = "white" if val > data.max() * 0.6 else "black"
                ax.text(j, i, str(val), ha="center", va="center", fontsize=6, color=color)

    ax.set_title(
        f"Tissue-of-Origin × Mechanism Article Counts\n"
        f"(Coverage: {assigned:,}/{len(index):,} records = {assigned/len(index):.0%})"
    )
    plt.colorbar(im, ax=ax, label="Article count", shrink=0.8, pad=0.02)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig14_tissue_mechanism_heatmap.pdf")
    fig.savefig(FIG_DIR / "fig14_tissue_mechanism_heatmap.png")
    plt.close()
    print(f"  {assigned} tissue-tagged records, {len(TISSUE_ORDER)}×{len(mechanisms)} matrix")


# ============================================================
# Fig 15: Designed-Combination Breakdown
# ============================================================

COMBINATION_CATEGORIES = [
    ("designed-combination-clinical", "Clinical\ndesigned", "#b71c1c"),
    ("designed-combination-preclinical", "Preclinical\ndesigned", "#e65100"),
    ("co-mention-only", "Co-mention\nonly", "#4fc3f7"),
    ("review-or-perspective-multi-lane", "Review /\nperspective", "#bdbdbd"),
]


def fig15_designed_combinations(index):
    """Horizontal stacked bar: multi-mechanism article classification."""
    print("Figure 15: Designed-combination breakdown...")

    counts = Counter()
    for e in index:
        combo = e.get("combination_evidence", "")
        if combo:
            counts[combo] += 1

    total = sum(counts.values())
    if total == 0:
        print("  No combination_evidence data — skipping")
        return

    fig, ax = plt.subplots(figsize=(12, 3.5))

    left = 0
    for key, label, color in COMBINATION_CATEGORIES:
        val = counts.get(key, 0)
        pct = val / total * 100 if total else 0
        ax.barh(0, val, left=left, color=color, edgecolor="white", linewidth=0.5, height=0.6)
        if val > 40:
            ax.text(left + val / 2, 0, f"{label}\n{val} ({pct:.1f}%)",
                    ha="center", va="center", fontsize=9, fontweight="bold",
                    color="white" if key in ("designed-combination-clinical",) else "black")
        left += val

    ax.set_xlim(0, total)
    ax.set_yticks([])
    ax.set_xlabel("Number of articles")
    ax.set_title(
        f"Classification of {total:,} Multi-Mechanism Articles\n"
        f"(articles tagged with 2+ therapeutic mechanisms)"
    )

    # Legend
    patches = [mpatches.Patch(color=c, label=f"{l.replace(chr(10), ' ')} ({counts.get(k, 0)})")
               for k, l, c in COMBINATION_CATEGORIES]
    ax.legend(handles=patches, fontsize=9, loc="upper right", ncol=2)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig15_designed_combinations.pdf")
    fig.savefig(FIG_DIR / "fig15_designed_combinations.png")
    plt.close()
    print(f"  {total} multi-mechanism articles classified")


# ============================================================
# Fig 16: Weighted Evidence Score by Mechanism
# ============================================================

EVIDENCE_TIER_WEIGHTS = {
    "phase3-clinical": 12.0,
    "phase2-clinical": 8.0,
    "phase1-clinical": 5.0,
    "clinical-other": 3.0,
    "preclinical-invivo": 2.0,
    "preclinical-invitro": 1.0,
    "theoretical": 0.5,
}


def _evidence_weight(entry: dict) -> float:
    """Heuristic quality weight — mirrors analyze_corpus.evidence_weight()."""
    level = entry.get("evidence_level", "")
    base = EVIDENCE_TIER_WEIGHTS.get(level)
    if not base:
        return 0.0

    pct = entry.get("icite_percentile") or 0
    try:
        pct = max(0.0, min(float(pct), 100.0))
    except (TypeError, ValueError):
        pct = 0.0
    citation_modifier = 1.0 + (pct / 200.0)

    year = entry.get("year") or 0
    if year:
        year = max(2015, min(int(year), 2026))
        recency_modifier = 0.9 + ((year - 2015) / (2026 - 2015)) * 0.2
    else:
        recency_modifier = 1.0

    return base * citation_modifier * recency_modifier


def fig16_weighted_evidence(index):
    """Horizontal lollipop: weighted evidence score per mechanism."""
    print("Figure 16: Weighted evidence scores...")

    from config import MECHANISM_KEYWORDS
    mechanisms = sorted(MECHANISM_KEYWORDS.keys())

    mech_data = {}
    for mech in mechanisms:
        mech_entries = [e for e in index if mech in e.get("mechanisms", [])]
        tagged = [e for e in mech_entries if e.get("evidence_level")]
        primary_like = [e for e in mech_entries if not _is_review_or_protocol(e)]
        total_weight = sum(_evidence_weight(e) for e in tagged)
        coverage = len(tagged) / len(primary_like) if primary_like else 0.0

        # Highest tier for dot color
        best_tier = ""
        for lvl in EVIDENCE_ORDER:
            if any(e.get("evidence_level") == lvl for e in tagged):
                best_tier = lvl
                break

        mech_data[mech] = {
            "weight": total_weight,
            "tagged": len(tagged),
            "coverage": coverage,
            "best_tier": best_tier,
        }

    sorted_mechs = sorted(mechanisms, key=lambda m: mech_data[m]["weight"])
    weights = [mech_data[m]["weight"] for m in sorted_mechs]
    colors = [EVIDENCE_COLORS.get(mech_data[m]["best_tier"], "#bdbdbd") for m in sorted_mechs]
    coverages = [mech_data[m]["coverage"] for m in sorted_mechs]

    fig, ax = plt.subplots(figsize=(10, 8))
    y = np.arange(len(sorted_mechs))

    # Lollipop stems
    for i, (w, c) in enumerate(zip(weights, colors)):
        ax.plot([0, w], [i, i], color=c, linewidth=1.5, alpha=0.6)
    # Dots
    ax.scatter(weights, y, c=colors, s=80, zorder=5, edgecolors="white", linewidths=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels([m.replace("-", " ") for m in sorted_mechs], fontsize=9)
    ax.set_xlabel("Weighted evidence score")
    ax.set_title("Weighted Evidence Score by Mechanism\n(tier × citation percentile × recency)")

    # Coverage labels on right
    for i, (w, cov) in enumerate(zip(weights, coverages)):
        ax.text(max(weights) * 1.02, i, f"{cov:.0%} cov", fontsize=7, va="center", color="gray")

    # Legend for tier colors
    tier_patches = [mpatches.Patch(color=EVIDENCE_COLORS[lvl], label=EVIDENCE_LABELS[lvl])
                    for lvl in EVIDENCE_ORDER if any(mech_data[m]["best_tier"] == lvl for m in mechanisms)]
    ax.legend(handles=tier_patches, fontsize=8, loc="lower right", title="Highest tier", title_fontsize=9)

    ax.set_xlim(0, max(weights) * 1.15)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig16_weighted_evidence.pdf")
    fig.savefig(FIG_DIR / "fig16_weighted_evidence.png")
    plt.close()
    print(f"  {len(mechanisms)} mechanisms, top score: {max(weights):.1f}")


# ============================================================
# Fig 17: DAMP Heatmap (TME Immune Coupling)
# ============================================================

TME_DIR = PROJECT_ROOT / "simulations" / "output" / "tme"


def fig17_damp_heatmap():
    """3-panel DAMP concentration heatmap: Control / RSL3 / SDT."""
    print("Figure 17: DAMP concentration heatmap...")

    treatment_keys = [
        ("control", "Control"),
        ("rsl3", "RSL3"),
        ("sdt", "SDT"),
    ]

    # Check that files exist
    for tx_key, _ in treatment_keys:
        path = TME_DIR / f"damp_field_{tx_key}.csv"
        if not path.exists():
            print(f"  {path} not found — run sim-tme first. Skipping.")
            return

    # Read immune kill counts from tme_summary.json (avoid hardcoding).
    # Since #224 item 2 the file is `{schema_version, conditions: [...]}`;
    # tolerate the legacy bare-array form for forward compatibility with
    # any pre-refactor outputs still on disk.
    summary_path = TME_DIR / "tme_summary.json"
    imm_kills_map = {}
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        conditions = summary.get("conditions", summary) if isinstance(summary, dict) else summary
        for r in conditions:
            if r.get("immune_mode") == "immune_on":
                imm_kills_map[r["treatment"]] = r.get("immune_kills", 0)

    # Load DAMP fields
    # NOTE: sim-tme encodes each treatment's DAMP field independently as u8
    # (each treatment normalized to its own max). The CSVs are NOT on the
    # same absolute scale. We normalize each panel independently so each
    # uses the full colormap range, and rely on the immune kill annotation
    # for quantitative comparison.
    panels = []
    for tx_key, tx_label in treatment_keys:
        data = np.loadtxt(TME_DIR / f"damp_field_{tx_key}.csv", delimiter=",")
        imm_kills = imm_kills_map.get(tx_label, 0)
        panels.append((tx_label, data, int(imm_kills)))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)

    for ax, (label, data, imm_kills) in zip(axes, panels):
        panel_max = data.max()
        if panel_max == 0:
            panel_max = 1.0
        normed = data / panel_max
        ax.imshow(normed, cmap="inferno", vmin=0, vmax=1, aspect="equal",
                  origin="upper")
        ax.set_title(f"{label}\n({imm_kills} immune kills)", fontsize=11)
        ax.set_xlabel("x (cells)")
        ax.set_xticks([0, 249, 499])
        ax.set_xticklabels(["0", "5", "10 mm"])
        ax.set_yticks([0, 249, 499])
        ax.set_yticklabels(["0", "5", "10 mm"])

    axes[0].set_ylabel("y (cells)")

    fig.suptitle(
        "DAMP Spatial Distribution After Immune Coupling\n"
        "(O$_2$ gradient $\\lambda$=120$\\mu$m, 500×500 grid, per-panel scaling)",
        fontsize=13, y=1.04)

    fig.savefig(FIG_DIR / "fig17_damp_heatmap.pdf")
    fig.savefig(FIG_DIR / "fig17_damp_heatmap.png")
    plt.close()
    print(f"  3 panels, immune kills: {[p[2] for p in panels]}")


# ============================================================
# Main
# ============================================================

def fig24_hypoxia_killcurve():
    """Figure 21: RSL3 vs SDT kill, normoxic vs hypoxic, and robustness across O2 penetration length (2D sim-tme)."""
    print("Figure 21 (fig24): Hypoxia kill-collapse (RSL3 vs SDT)...")
    if not TME_SUMMARY.exists():
        print(f"  {TME_SUMMARY} not found — run `cargo run --release -p sim-tme` first. Skipping.")
        return
    data = json.loads(TME_SUMMARY.read_text())
    conds = data["conditions"] if isinstance(data, dict) and "conditions" in data else data
    base = [c for c in conds if c.get("immune_mode") == "off"]

    def kill(treatment, o2):
        for c in base:
            if c["treatment"] == treatment and c.get("o2_condition") == o2:
                return c["overall_kill_rate"] * 100
        return None

    lambdas = [80, 100, 120, 150]
    rsl3_grad = [kill("RSL3", f"gradient_{l}um") for l in lambdas]
    sdt_grad = [kill("SDT", f"gradient_{l}um") for l in lambdas]
    rsl3_norm, sdt_norm = kill("RSL3", "uniform"), kill("SDT", "uniform")
    if None in rsl3_grad + sdt_grad + [rsl3_norm, sdt_norm]:
        print("  Missing conditions in tme_summary.json — skipping")
        return

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5))

    rsl3_hyp, sdt_hyp = float(np.mean(rsl3_grad)), float(np.mean(sdt_grad))
    groups = ["RSL3\n(GPX4 inhibitor)", "SDT\n(exogenous ROS)"]
    x = np.arange(len(groups))
    w = 0.36
    axA.bar(x - w / 2, [rsl3_norm, sdt_norm], w, label="Normoxic (uniform O$_2$)", color="#4C72B0")
    axA.bar(x + w / 2, [rsl3_hyp, sdt_hyp], w, label="Hypoxic (O$_2$ gradient)", color="#C44E52")
    for xi, (n, h) in enumerate([(rsl3_norm, rsl3_hyp), (sdt_norm, sdt_hyp)]):
        axA.text(xi - w / 2, n + 1.5, f"{n:.1f}%", ha="center", fontsize=9)
        axA.text(xi + w / 2, h + 1.5, f"{h:.1f}%", ha="center", fontsize=9)
    axA.set_xticks(x)
    axA.set_xticklabels(groups)
    axA.set_ylabel("Overall tumor kill (%)")
    axA.set_ylim(0, 105)
    axA.set_title("(a) Kill collapse under hypoxia")
    axA.legend(fontsize=8, loc="upper left")
    axA.annotate(f"{rsl3_norm:.1f}% $\\to$ {rsl3_hyp:.1f}%\n(~{rsl3_norm / max(rsl3_hyp, 0.01):.0f}$\\times$ collapse)",
                 xy=(0, rsl3_hyp), xytext=(0.55, 55), fontsize=8, color="#C44E52", ha="center")

    axB.plot(lambdas, sdt_grad, "o-", color="#C44E52", label="SDT (hypoxic)")
    axB.plot(lambdas, rsl3_grad, "s-", color="#4C72B0", label="RSL3 (hypoxic)")
    axB.axhline(sdt_norm, ls="--", color="#C44E52", alpha=0.4, lw=1, label="SDT normoxic")
    axB.axhline(rsl3_norm, ls="--", color="#4C72B0", alpha=0.4, lw=1, label="RSL3 normoxic")
    axB.set_xlabel("O$_2$ penetration length $\\lambda$ ($\\mu$m)")
    axB.set_ylabel("Overall tumor kill (%)")
    axB.set_ylim(-3, 100)
    axB.set_xticks(lambdas)
    axB.set_title("(b) Robust across penetration length")
    axB.legend(fontsize=8, loc="center right")

    fig.suptitle("Hypoxia collapses pharmacologic ferroptosis but not exogenous ROS (2D model)", fontsize=12, y=1.02)
    fig.text(0.5, -0.06,
             "2D sim-tme, immune off. SDT is modeled as O$_2$-independent — an optimistic upper bound; "
             "SDT's own O$_2$-dependence is contested (see manuscript §7.1).",
             ha="center", fontsize=7.5, style="italic", color="gray")
    fig.savefig(FIG_DIR / "fig24_hypoxia_killcurve.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig24_hypoxia_killcurve.png", bbox_inches="tight")
    plt.close()
    print(f"  RSL3 {rsl3_norm:.1f}%->{rsl3_hyp:.2f}%, SDT {sdt_norm:.1f}%->{sdt_hyp:.1f}%")


def fig25_bliss_synergy():
    """Figure 22: Bliss synergy of dual-pathway depletion — RSL3+FSP1i observed vs expected, plus pairwise scores."""
    print("Figure 22 (fig25): Bliss synergy (dual-pathway depletion)...")
    if not COMBO_SUMMARY.exists():
        print(f"  {COMBO_SUMMARY} not found — run `cargo run --release -p sim-combo-mech` first. Skipping.")
        return
    data = json.loads(COMBO_SUMMARY.read_text())
    combos = data["combinations"] if isinstance(data, dict) and "combinations" in data else data

    def find(a, b):
        for c in combos:
            if {c["drug_a"], c["drug_b"]} == {a, b}:
                return c
        return None

    rf = find("RSL3", "FSP1i")
    if rf is None:
        print("  RSL3+FSP1i combination not found — skipping")
        return

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5))

    labels = ["RSL3\nalone", "FSP1i\nalone", "Bliss\nexpected", "Observed\ncombination"]
    vals = [rf["rate_a"] * 100, rf["rate_b"] * 100, rf["bliss_prediction"] * 100, rf["rate_combo"] * 100]
    colors = ["#4C72B0", "#55A868", "#999999", "#C44E52"]
    bars = axA.bar(labels, vals, color=colors)
    for bar, v in zip(bars, vals):
        axA.text(bar.get_x() + bar.get_width() / 2, v + 1.5, f"{v:.1f}%", ha="center", fontsize=9)
    axA.set_ylabel("Persister kill (%)")
    axA.set_ylim(0, 100)
    axA.set_title("(a) RSL3 + FSP1i: dual-pathway synergy")
    axA.annotate(f"{rf['synergy_score']:.2f}$\\times$ synergy", xy=(3, vals[3]), xytext=(1.9, 93),
                 fontsize=10, fontweight="bold", color="#C44E52",
                 arrowprops=dict(arrowstyle="->", color="#C44E52"))

    meaningful = [(c, c["synergy_score"]) for c in combos if "SDT" not in (c["drug_a"], c["drug_b"])]
    meaningful.sort(key=lambda t: t[1], reverse=True)
    names = [f"{c['drug_a']}+{c['drug_b']}" for c, _ in meaningful]
    scores = [s for _, s in meaningful]
    cols = ["#C44E52" if n == "RSL3+FSP1i" else "#4C72B0" for n in names]
    axB.barh(names[::-1], scores[::-1], color=cols[::-1])
    axB.axvline(1.0, ls="--", color="gray", lw=1, label="additive (1.0$\\times$)")
    for i, s in enumerate(scores[::-1]):
        axB.text(s + 0.02, i, f"{s:.2f}$\\times$", va="center", fontsize=8)
    axB.set_xlabel("Bliss synergy score (observed / expected)")
    axB.set_xlim(0, max(scores) * 1.25)
    axB.set_title("(b) Pairwise synergy (SDT pairs excluded: ceiling)")
    axB.legend(fontsize=8, loc="lower right")

    fig.suptitle("Dual-pathway (GPX4 + FSP1) depletion is synergistic", fontsize=12, y=1.02)
    fig.text(0.5, -0.04,
             "1,000 persister cells/condition, 2D culture params. Drug potencies are estimates; the "
             "directional finding (dual-pathway > single) held across the ±50% sensitivity sweep (§5).",
             ha="center", fontsize=7.5, style="italic", color="gray")
    fig.savefig(FIG_DIR / "fig25_bliss_synergy.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig25_bliss_synergy.png", bbox_inches="tight")
    plt.close()
    print(f"  RSL3+FSP1i: {rf['rate_combo'] * 100:.1f}% obs vs {rf['bliss_prediction'] * 100:.1f}% expected = {rf['synergy_score']:.2f}x")


def fig26_vulnerability_window():
    """Figure 23: ferroptosis-sensitive window — RSL3 closes by ~day 3 (GPX4 recovery), SDT stays open ~4 weeks."""
    print("Figure 23 (fig26): Vulnerability window (RSL3 vs SDT timing)...")
    if not WINDOW_JSON.exists():
        print(f"  {WINDOW_JSON} not found — run `cargo run --release -p sim-window` first. Skipping.")
        return
    data = json.loads(WINDOW_JSON.read_text())
    days = sorted(set(r["timepoint_days"] for r in data))
    present = {(r["treatment"], r["timepoint_days"]) for r in data}
    if not all((tx, d) in present for tx in ("RSL3", "SDT") for d in days):
        print("  incomplete window data (missing treatment/timepoint) — skipping")
        return

    def series(tx, key):
        m = {r["timepoint_days"]: r for r in data if r["treatment"] == tx}
        return [m[d][key] for d in days]

    x = np.arange(len(days))  # even spacing; label with actual day values
    labels = [f"{d:g}" for d in days]
    rsl3 = [v * 100 for v in series("RSL3", "death_rate")]
    sdt = [v * 100 for v in series("SDT", "death_rate")]
    rsl3_lo, rsl3_hi = [v * 100 for v in series("RSL3", "ci_low")], [v * 100 for v in series("RSL3", "ci_high")]
    sdt_lo, sdt_hi = [v * 100 for v in series("SDT", "ci_low")], [v * 100 for v in series("SDT", "ci_high")]
    gpx4_rsl3 = series("RSL3", "mean_gpx4")
    win_end = max(i for i, d in enumerate(days) if d <= 3.0)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5))

    axA.plot(x, sdt, "o-", color="#C44E52", label="SDT (exogenous ROS)")
    axA.fill_between(x, sdt_lo, sdt_hi, color="#C44E52", alpha=0.15)
    axA.plot(x, rsl3, "s-", color="#4C72B0", label="RSL3 (GPX4 inhibitor)")
    axA.fill_between(x, rsl3_lo, rsl3_hi, color="#4C72B0", alpha=0.15)
    axA.axvspan(-0.3, win_end, color="#4C72B0", alpha=0.07)
    axA.text(win_end / 2.0, 55, "RSL3 window\nopen", ha="center", fontsize=8, color="#4C72B0")
    axA.set_xticks(x)
    axA.set_xticklabels(labels)
    axA.set_xlabel("Time post-chemotherapy (days)")
    axA.set_ylabel("Persister kill (%)")
    axA.set_ylim(-3, 107)
    axA.set_title("(a) Treatment window: RSL3 closes, SDT stays open")
    axA.legend(fontsize=8, loc="center right")
    axA.annotate("closes ~day 3", xy=(win_end, rsl3[win_end]), xytext=(win_end + 0.4, 24),
                 fontsize=8, color="#4C72B0", arrowprops=dict(arrowstyle="->", color="#4C72B0"))

    axB.plot(x, rsl3, "s-", color="#4C72B0", label="RSL3 kill (%)")
    axB.set_xticks(x)
    axB.set_xticklabels(labels)
    axB.set_xlabel("Time post-chemotherapy (days)")
    axB.set_ylabel("RSL3 kill (%)", color="#4C72B0")
    axB.set_ylim(-3, 50)
    axB.tick_params(axis="y", labelcolor="#4C72B0")
    axB2 = axB.twinx()
    axB2.plot(x, gpx4_rsl3, "^--", color="#55A868", label="mean GPX4 (recovered fraction)")
    axB2.set_ylabel("Mean GPX4 (recovered fraction)", color="#55A868")
    axB2.tick_params(axis="y", labelcolor="#55A868")
    axB.set_title("(b) Why: GPX4 re-expression closes the window")
    lA, llA = axB.get_legend_handles_labels()
    lB, llB = axB2.get_legend_handles_labels()
    axB.legend(lA + lB, llA + llB, fontsize=8, loc="center right")

    fig.suptitle("The ferroptosis-sensitive window: days for RSL3, weeks for SDT", fontsize=12, y=1.02)
    fig.text(0.5, -0.04,
             "100,000 cells/condition; x-axis shows sampled timepoints (not linear in time). Defense-recovery "
             "half-times (GPX4 3 d, FSP1 7 d, NRF2 5 d, GSH 1 d) are literature-estimated, so window durations "
             "are approximate until experimentally validated.",
             ha="center", fontsize=7.5, style="italic", color="gray")
    fig.savefig(FIG_DIR / "fig26_vulnerability_window.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig26_vulnerability_window.png", bbox_inches="tight")
    plt.close()
    print(f"  RSL3 {rsl3[0]:.1f}%@0d -> {rsl3[win_end]:.1f}%@3d -> {rsl3[-1]:.2f}%@28d; SDT {sdt[-1]:.1f}%@28d")


def fig8_simulation_by_treatment():
    """Manuscript Figure 8 (Tier-1, #285): depth-kill curves. PDT (light, Beer-Lambert)
    collapses with depth; SDT (ultrasound, acoustic) penetrates to centimeters; RSL3
    (systemic drug) is a depth-independent uniform baseline whose limit is biochemical,
    not penetration. Calibrated 2D physics (sim-spatial). Replaces the prior
    externally-post-processed Figure 8 with a tracked, reproducible generator."""
    print("Figure 8 (fig8_simulation_by_treatment): Depth-kill curves (PDT vs SDT vs RSL3)...")
    if not SPATIAL_CURVES.exists():
        print(f"  {SPATIAL_CURVES} not found — run `cargo run --release -p sim-spatial` first. Skipping.")
        return
    rows = defaultdict(list)
    with open(SPATIAL_CURVES) as f:
        for r in csv.DictReader(f):
            rows[r["treatment"]].append((float(r["depth_um"]), float(r["death_rate"]), int(r["n_cells"])))
    if not all(t in rows for t in ("PDT", "SDT", "RSL3")):
        print("  missing PDT/SDT/RSL3 in depth_kill_curves.csv — skipping")
        return

    # Pool the 20-µm rows into coarser depth bins, weighting each row by its
    # tumor-cell count, so the spheroid's sparse poles (few cells/row) do not
    # add noise to the curve. Pooled rate = sum(dead) / sum(total) per bin.
    BIN_UM = 250.0

    def binned(tx):
        agg = defaultdict(lambda: [0.0, 0])  # bin_center_um -> [dead_weighted, n]
        for d, rate, n in rows[tx]:
            if n <= 0:
                continue
            b = int(d // BIN_UM) * BIN_UM + BIN_UM / 2.0
            agg[b][0] += rate * n
            agg[b][1] += n
        xs = sorted(agg)
        return ([x / 1000.0 for x in xs], [agg[x][0] / agg[x][1] * 100.0 for x in xs])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5))

    style = {
        "PDT": ("#C44E52", "o-", "PDT (light, Beer-Lambert)"),
        "SDT": ("#4C72B0", "s-", "SDT (ultrasound, acoustic)"),
        "RSL3": ("#55A868", "^-", "RSL3 (systemic drug, uniform)"),
    }
    series = {}
    for tx in ("SDT", "PDT", "RSL3"):
        x, y = binned(tx)
        series[tx] = (x, y)
        col, ls, lab = style[tx]
        axA.plot(x, y, ls, color=col, label=lab, markersize=4, lw=1.8)
    axA.set_xlabel("Depth from irradiated surface (mm)")
    axA.set_ylabel("Tumor kill (%)")
    axA.set_ylim(-4, 106)
    axA.set_title("(a) Observed kill vs depth (2D sim)")
    axA.legend(fontsize=8, loc="center right")
    # PDT collapse annotation (shallow vs deep observed kill).
    px, py = series["PDT"]
    if py:
        axA.annotate(
            f"PDT: {py[0]:.0f}% at surface,\n{py[-1]:.0f}% at {px[-1]:.0f} mm",
            xy=(px[-1], py[-1]), xytext=(px[-1] * 0.45, 60), fontsize=8, color="#C44E52",
            ha="center", arrowprops=dict(arrowstyle="->", color="#C44E52"))

    # Panel (b): the penetration physics that drives panel (a), computed from the
    # model's OWN equations and default constants (ferroptosis-core/src/physics.rs
    # + params.rs SpatialParams::default):
    #   PDT  I(z) = I0 * exp(-mu_eff * z_mm),  mu_eff = 0.31 /mm  (delta ~ 3.2 mm, 630 nm)
    #   SDT  I(z) = I0 * 10^(-alpha * f * z_cm / 10),  alpha = 0.7 dB/cm/MHz, f = 1 MHz
    #   RSL3 uniform = 1.0 (systemic drug; no depth attenuation)
    #
    # DRIFT GUARD: these constants are hardcoded here (a Python re-implementation
    # of the Rust physics), and this panel assumes the DEFAULT sim-spatial flags
    # (`--dli-h 0`, `Photosensitizer::Uniform(1.0)`) so the PDT drug-yield is 1.0
    # and panel (b) equals the sim's per-cell multiplier that produced panel (a).
    # If params.rs retunes these defaults, or the CSV is regenerated with a
    # non-default `--photosensitizer`/`--dli-h`, the two panels desync. The values
    # are pinned against params.rs by tests/test_depth_kill_physics_constants.py;
    # update BOTH if you change them, and regenerate the CSV with default flags.
    PDT_MU_EFF_PER_MM = 0.31  # params.rs SpatialParams::default().pdt_mu_eff
    SDT_ALPHA_DB_CM_MHZ = 0.7  # params.rs SpatialParams::default().sdt_alpha
    SDT_FREQ_MHZ = 1.0  # params.rs SpatialParams::default().sdt_freq_mhz
    z_mm = np.linspace(0.0, 10.0, 200)
    pdt_I = np.exp(-PDT_MU_EFF_PER_MM * z_mm)
    sdt_I = 10.0 ** (-(SDT_ALPHA_DB_CM_MHZ * SDT_FREQ_MHZ * (z_mm / 10.0)) / 10.0)
    rsl3_I = np.ones_like(z_mm)
    axB.plot(z_mm, sdt_I * 100, "-", color="#4C72B0", lw=2,
             label="SDT acoustic ($\\alpha$=0.7 dB/cm/MHz)")
    axB.plot(z_mm, pdt_I * 100, "-", color="#C44E52", lw=2,
             label="PDT light ($\\mu_{\\mathrm{eff}}$=0.31/mm, $\\delta{\\approx}$3.2 mm)")
    # RSL3 dotted (not solid like the attenuating modalities) to flag that this is
    # DRUG AVAILABILITY, not kill: it is flat at 100% here yet flat near zero in
    # panel (a) — a biochemical, not a penetration, limit.
    axB.plot(z_mm, rsl3_I * 100, ":", color="#55A868", lw=2.5,
             label="RSL3 availability (uniform; kills little, see a)")
    axB.axhline(50, ls=":", color="gray", lw=0.8)
    axB.set_xlabel("Depth from irradiated surface (mm)")
    axB.set_ylabel("Relative energy / drug availability (% of surface)")
    axB.set_ylim(0, 105)
    axB.set_title("(b) Why: penetration physics")
    axB.legend(fontsize=7.5, loc="upper right")

    fig.suptitle(
        "Penetration sets modality reach: light is millimeters, ultrasound is centimeters (2D model)",
        fontsize=12, y=1.02)
    fig.text(0.5, -0.05,
             "2D sim-spatial, 1 cm tissue. Depth profiles follow well-measured physics (Beer-Lambert optics, "
             "acoustic attenuation; high confidence). RSL3 reaches every depth but kills little, a biochemical "
             "limit, not a penetration one. Absolute kill % rests on uncalibrated biochemistry, so read the "
             "profile shape, not the magnitudes. SDT is modeled as O$_2$-independent, an optimistic upper bound "
             "(Section 7.1).",
             ha="center", fontsize=7.5, style="italic", color="gray")
    fig.savefig(FIG_DIR / "fig8_simulation_by_treatment.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig8_simulation_by_treatment.png", bbox_inches="tight")
    plt.close()
    _, sy = series["SDT"]
    if py and sy:
        print(f"  PDT {py[0]:.0f}%->{py[-1]:.0f}% over {px[-1]:.0f} mm; "
              f"SDT {sy[0]:.0f}%->{sy[-1]:.0f}%; depth bins={len(px)}")
    else:
        print("  (no tumor cells in one or more depth curves)")


def fig27_resistance_asymmetry():
    """Manuscript Figure 24 (flagship, #285): the resistance-mechanism asymmetry.
    A 2x2 panel — under each TME resistance mechanism (hypoxia / stromal / pH /
    immune) pharmacologic RSL3 collapses while physical SDT holds. Every panel is
    captioned with its confidence tier so the contested hypoxia leg is not
    entrenched. UNCALIBRATED 2D biochemistry (sim-tme): the cross-modality
    direction is the result, not the magnitudes. Each panel uses the SAME metric
    the corresponding manuscript section (§7.1-7.4) reports, so the figure and the
    prose agree. Reads all values from tme_summary.json."""
    print("Figure 24 (fig27): Resistance-mechanism asymmetry (flagship)...")
    if not TME_SUMMARY.exists():
        print(f"  {TME_SUMMARY} not found — run `cargo run --release -p sim-tme` first. Skipping.")
        return
    data = json.loads(TME_SUMMARY.read_text())
    conds = data["conditions"] if isinstance(data, dict) and "conditions" in data else data

    def find(treatment, **kw):
        for c in conds:
            if c["treatment"] == treatment and all(c.get(k) == v for k, v in kw.items()):
                return c
        return None

    G = "gradient_120um"
    RSL3_C, SDT_C, GHOST = "#4C72B0", "#C44E52", "#C9C9C9"

    # Each panel uses the SAME metric its manuscript section reports, so the
    # figure magnitudes match the prose (and avoid the immune confound):
    #  - (a) Hypoxia (§7.1): OVERALL kill, normoxic(uniform) vs hypoxic(gradient).
    #        Computed at immune_mode=off (the clean O2-only comparison).
    #  - (b) Stromal (§7.3): stromal_adjacent_kill_rate (kill among the CAF-adjacent
    #        boundary cells, where the effect lives). immune_on baseline.
    #  - (c) pH (§7.4): ferroptosis_kills COUNT (immune-pure — a separate counter
    #        from immune_kills). immune_on baseline.
    #  - (d) Immune (§7.2): immune-kill COUNT.
    # Panels (b)-(d) sit on the gradient_120um + immune_on baseline (the only
    # conditions the sim runs those mechanisms under); panel (a) is immune_off.
    hyp = {t: (find(t, o2_condition="uniform", immune_mode="off"),
               find(t, o2_condition=G, immune_mode="off", stromal_mode=None, ph_mode=None))
           for t in ("RSL3", "SDT")}
    strm = {t: (find(t, o2_condition=G, immune_mode="immune_on", stromal_mode="off"),
                find(t, o2_condition=G, immune_mode="immune_on", stromal_mode="stromal_on"))
            for t in ("RSL3", "SDT")}
    phh = {t: (find(t, o2_condition=G, immune_mode="immune_on", stromal_mode="off"),
               find(t, o2_condition=G, immune_mode="immune_on", ph_mode="ph_on"))
           for t in ("RSL3", "SDT")}
    imm = {t: find(t, o2_condition=G, immune_mode="immune_on", stromal_mode="off")
           for t in ("RSL3", "SDT")}

    need = ([hyp[t][i] for t in ("RSL3", "SDT") for i in (0, 1)]
            + [strm[t][i] for t in ("RSL3", "SDT") for i in (0, 1)]
            + [phh[t][i] for t in ("RSL3", "SDT") for i in (0, 1)]
            + [imm[t] for t in ("RSL3", "SDT")])
    if any(x is None for x in need):
        print("  missing conditions in tme_summary.json — skipping")
        return

    fig, ((axH, axS), (axP, axI)) = plt.subplots(2, 2, figsize=(11, 8.5))

    def _legend(ax, base_lab, stress_lab):
        ax.legend(handles=[mpatches.Patch(color=GHOST, label=base_lab),
                           mpatches.Patch(color=RSL3_C, label=f"RSL3, {stress_lab}"),
                           mpatches.Patch(color=SDT_C, label=f"SDT, {stress_lab}")],
                  fontsize=6.5, loc="center left")

    def killbars(ax, pair, metric, title, ylabel, base_lab, stress_lab):
        """Grouped baseline/stressed bars for RSL3 and SDT, kill % on a 0-108 axis."""
        base = {t: pair[t][0][metric] * 100 for t in ("RSL3", "SDT")}
        strs = {t: pair[t][1][metric] * 100 for t in ("RSL3", "SDT")}
        x = np.arange(2)
        w = 0.36
        ax.bar(x - w / 2, [base["RSL3"], base["SDT"]], w, color=GHOST)
        ax.bar(x + w / 2, [strs["RSL3"], strs["SDT"]], w, color=[RSL3_C, SDT_C])
        for xi, t in enumerate(("RSL3", "SDT")):
            ax.text(xi - w / 2, base[t] + 2, f"{base[t]:.1f}", ha="center", fontsize=8, color="#666")
            ax.text(xi + w / 2, strs[t] + 2, f"{strs[t]:.1f}", ha="center", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(["RSL3\n(pharmacologic)", "SDT\n(physical)"])
        ax.set_ylim(0, 108)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10.5)
        _legend(ax, base_lab, stress_lab)
        return base, strs

    def countbars(ax, pair, field, title, ylabel, base_lab, stress_lab):
        """Grouped baseline/stressed COUNT bars for RSL3 and SDT, log y-axis."""
        base = {t: int(pair[t][0][field] or 0) for t in ("RSL3", "SDT")}
        strs = {t: int(pair[t][1][field] or 0) for t in ("RSL3", "SDT")}
        x = np.arange(2)
        w = 0.36
        ax.bar(x - w / 2, [max(base["RSL3"], 1), max(base["SDT"], 1)], w, color=GHOST)
        ax.bar(x + w / 2, [max(strs["RSL3"], 1), max(strs["SDT"], 1)], w, color=[RSL3_C, SDT_C])
        ax.set_yscale("log")
        top = max(base["SDT"], strs["SDT"], 10) * 4
        ax.set_ylim(0.7, top)
        for xi, t in enumerate(("RSL3", "SDT")):
            ax.text(xi - w / 2, max(base[t], 1) * 1.3, f"{base[t]}", ha="center", fontsize=7.5, color="#666")
            ax.text(xi + w / 2, max(strs[t], 1) * 1.3, f"{strs[t]}", ha="center", fontsize=7.5)
        ax.set_xticks(x)
        ax.set_xticklabels(["RSL3\n(pharmacologic)", "SDT\n(physical)"])
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10.5)
        _legend(ax, base_lab, stress_lab)
        return base, strs

    # (a) Hypoxia — overall kill %, matches §7.1 (RSL3 3.7->0.1, SDT 91.9->87.8).
    bH, sH = killbars(axH, hyp, "overall_kill_rate",
                      "(a) Hypoxia  [directional; SDT O$_2$-independence contested]",
                      "Overall tumor kill (%)", "normoxic O$_2$", "hypoxic O$_2$")
    # (b) Stromal — CAF-adjacent kill %, matches §7.3 (RSL3 3.0->1.5, SDT 96.1->91.2).
    bS, sS = killbars(axS, strm, "stromal_adjacent_kill_rate",
                      "(b) Stromal / CAF  [uncalibrated]",
                      "CAF-adjacent cell kill (%)", "no CAF", "CAF")
    # (c) pH — ferroptosis-kill COUNT (immune-pure), matches §7.4 (RSL3 163->77).
    bP, sP = countbars(axP, phh, "ferroptosis_kills",
                       "(c) Acidic pH  [low confidence; RSL3 pKa most uncertain]",
                       "Ferroptosis kills [log]", "neutral pH", "acidic pH")

    # (d) Immune — immune-kill COUNT (one condition per treatment), with ratio.
    ik = {t: max(int(imm[t]["immune_kills"] or 0), 0) for t in ("RSL3", "SDT")}
    bars = axI.bar(["RSL3\n(pharmacologic)", "SDT\n(physical)"],
                   [max(ik["RSL3"], 1), max(ik["SDT"], 1)], width=0.55, color=[RSL3_C, SDT_C])
    axI.set_yscale("log")
    axI.set_ylim(0.7, max(ik["SDT"], 10) * 3)
    axI.set_ylabel("Immune (ICD) kills [log]")
    for b, t in zip(bars, ("RSL3", "SDT")):
        axI.text(b.get_x() + b.get_width() / 2, max(ik[t], 1) * 1.15, f"{ik[t]}", ha="center", fontsize=9)
    ratio = ik["SDT"] / max(ik["RSL3"], 1)
    axI.set_title("(d) Immune / ICD coupling  [directional; 2D ceiling]", fontsize=10.5)
    axI.annotate(f"{ratio:.0f}:1 in 2D\n(~4:1 in 3D)", xy=(1, max(ik["SDT"], 1)),
                 xytext=(0.35, max(ik["SDT"], 1) * 0.9), fontsize=9, fontweight="bold",
                 color="#C44E52", ha="center")

    fig.suptitle("Resistance-mechanism asymmetry: RSL3 (pharmacologic) collapses, SDT (physical) holds (2D model)",
                 fontsize=13, y=0.99)
    fig.text(0.5, -0.03,
             "2D sim-tme. Each panel uses the metric its manuscript section reports, so the figure and prose agree: "
             "hypoxia = overall kill (§7.1); stromal = kill among CAF-adjacent boundary cells (§7.3); pH = ferroptosis "
             "kills, an immune-free counter (§7.4); immune = ICD kill count (§7.2). Panel (a) is computed without the "
             "immune layer (clean O$_2$-only comparison); panels (b)-(d) share the gradient-O$_2$ + immune-on baseline "
             "the sim runs those mechanisms under (the pH 'neutral' bar reuses the stromal-off run, the only available "
             "reference). Confidence tiers differ per panel (titles): the hypoxia leg is the most contested (SDT modeled "
             "O$_2$-independent, an optimistic upper bound, §7.1); the immune 2D ratio over-extrapolates (~4:1 under 3D "
             "volumetric dilution). Magnitudes rest on uncalibrated biochemistry; the cross-modality direction is the "
             "result, not the numbers.",
             ha="center", fontsize=7, style="italic", color="gray", wrap=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(FIG_DIR / "fig27_resistance_asymmetry.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig27_resistance_asymmetry.png", bbox_inches="tight")
    plt.close()
    print(f"  hypoxia RSL3 {bH['RSL3']:.1f}->{sH['RSL3']:.1f}% / SDT {bH['SDT']:.1f}->{sH['SDT']:.1f}%; "
          f"stromal RSL3 {bS['RSL3']:.1f}->{sS['RSL3']:.1f}%; pH ferro RSL3 {bP['RSL3']}->{sP['RSL3']}; "
          f"immune RSL3 {ik['RSL3']} vs SDT {ik['SDT']} ({ratio:.0f}:1)")


def main():
    print("Loading corpus...")
    articles = load_corpus()
    print(f"  Loaded {len(articles)} articles")

    print("Loading index...")
    index = load_index()
    print(f"  Loaded {len(index)} index records\n")

    fig1_ferroptosis_comparison(articles)
    fig2_mechanism_heatmap(articles)
    fig3_literature_disconnect(articles)
    fig4_molecular_overlap(articles)
    fig5_publication_trends(articles)
    fig6_sdt_chain_evidence(articles)

    # fig8 (spatial depth-kill curves) is generated here from sim-spatial's
    # depth_kill_curves.csv (run `cargo run --release -p sim-spatial` first).
    fig8_simulation_by_treatment()

    # Note: fig7 (Monte Carlo simulation) is still generated by the Rust binary
    # (sim-original), not this script. Run it separately:
    #   cargo run --release -p sim-original   -> fig7_monte_carlo_simulation

    fig9_evidence_tiers(index)
    fig10_invivo_comparison()
    fig11_mufa_sweep()
    fig12_pathway_targets(index)

    fig13_gold_set_eval()
    fig14_tissue_mechanism_heatmap(index)
    fig15_designed_combinations(index)
    fig16_weighted_evidence(index)
    fig17_damp_heatmap()

    # Tier-1 quantitative simulation figures (#285): manuscript figures 21, 22.
    fig24_hypoxia_killcurve()
    fig25_bliss_synergy()
    # Tier-2 (#285): manuscript figure 23 — treatment-timing window.
    fig26_vulnerability_window()
    # Flagship (#285): manuscript figure 24 — 2x2 resistance-mechanism asymmetry.
    fig27_resistance_asymmetry()

    print(f"\nAll figures saved to {FIG_DIR}/")
    print("Files:")
    for f in sorted(FIG_DIR.glob("fig*")):
        print(f"  {f.name} ({f.stat().st_size // 1024}KB)")


if __name__ == "__main__":
    main()
