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
"""

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
# Main
# ============================================================

def main():
    print("Loading corpus...")
    articles = load_corpus()
    print(f"  Loaded {len(articles)} articles\n")

    fig1_ferroptosis_comparison(articles)
    fig2_mechanism_heatmap(articles)
    fig3_literature_disconnect(articles)
    fig4_molecular_overlap(articles)
    fig5_publication_trends(articles)
    fig6_sdt_chain_evidence(articles)

    print(f"\nAll figures saved to {FIG_DIR}/")
    print("Files:")
    for f in sorted(FIG_DIR.glob("fig*")):
        print(f"  {f.name} ({f.stat().st_size // 1024}KB)")


if __name__ == "__main__":
    main()
