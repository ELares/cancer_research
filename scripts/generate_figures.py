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
# Main
# ============================================================

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

    # Note: fig7 (Monte Carlo simulation) and fig8 (spatial depth-kill curves)
    # are generated by the Rust simulation binaries (sim-original, sim-spatial),
    # not by this Python script. Run them separately:
    #   cargo run --release -p sim-original   -> fig7_monte_carlo_simulation
    #   cargo run --release -p sim-spatial    -> fig8_simulation_by_treatment

    fig9_evidence_tiers(index)
    fig10_invivo_comparison()
    fig11_mufa_sweep()
    fig12_pathway_targets(index)

    fig13_gold_set_eval()
    fig14_tissue_mechanism_heatmap(index)
    fig15_designed_combinations(index)
    fig16_weighted_evidence(index)

    print(f"\nAll figures saved to {FIG_DIR}/")
    print("Files:")
    for f in sorted(FIG_DIR.glob("fig*")):
        print(f"  {f.name} ({f.stat().st_size // 1024}KB)")


if __name__ == "__main__":
    main()
