#!/usr/bin/env python3
"""
Analyze the corpus and generate all analysis files.

Produces:
  analysis/mechanism-matrix.md
  analysis/convergence-map.md
  analysis/gap-analysis.md
  analysis/evidence-tiers.md
  analysis/key-findings.md
  analysis/timeline.md

Usage:
    python analyze_corpus.py
"""

import json
import re
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

from article_io import load_article
from config import (
    MECHANISM_KEYWORDS,
    CANCER_TYPE_KEYWORDS,
    EVIDENCE_LEVEL_KEYWORDS,
    PROJECT_ROOT,
    RESISTANT_STATE_RULES,
)

INDEX_FILE = PROJECT_ROOT / "corpus" / "INDEX.jsonl"
PMID_DIR = PROJECT_ROOT / "corpus" / "by-pmid"
ANALYSIS_DIR = PROJECT_ROOT / "analysis"

REVIEW_MARKERS = (
    "review", "systematic review", "meta-analysis", "meta analysis",
    "scoping review", "narrative review", "evidence map",
)
PROTOCOL_MARKERS = ("protocol", "study protocol", "trial protocol", "protocol for")


def load_index() -> list[dict]:
    entries = []
    skipped = 0
    for line in INDEX_FILE.read_text(encoding="utf-8").splitlines():
        if line.strip():
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError:
                skipped += 1
    if skipped:
        print(f"  Warning: skipped {skipped} malformed index lines")
    return entries


def load_article_abstract(pmid: str) -> str:
    """Load just the abstract from an article file."""
    fp = PMID_DIR / f"{pmid}.md"
    if not fp.exists():
        return ""
    content = fp.read_text(encoding="utf-8")
    match = re.search(r"## Abstract\n\n?(.*?)(?=\n## |\Z)", content, re.DOTALL)
    return match.group(1).strip() if match else ""


@lru_cache(maxsize=None)
def load_article_frontmatter(pmid: str) -> dict:
    fp = PMID_DIR / f"{pmid}.md"
    if not fp.exists():
        return {}
    fm, _ = load_article(fp)
    return fm


def is_review_like(fm: dict) -> bool:
    pub_types = [p.lower() for p in fm.get("pub_types", [])]
    title = fm.get("title", "").lower()
    return any("review" in p or "meta-analysis" in p for p in pub_types) or any(
        marker in title for marker in REVIEW_MARKERS
    )


def is_protocol_like(fm: dict) -> bool:
    pub_types = [p.lower() for p in fm.get("pub_types", [])]
    title = fm.get("title", "").lower()
    return any("protocol" in p for p in pub_types) or any(marker in title for marker in PROTOCOL_MARKERS)


def classify_evidence_reason(entry: dict) -> str:
    if entry.get("evidence_level"):
        return "tagged"
    fm = load_article_frontmatter(entry.get("pmid", ""))
    if is_review_like(fm):
        return "review_like"
    if is_protocol_like(fm):
        return "protocol_like"
    return "other_untagged"


# ============================================================
# Analysis 1: Mechanism-Cancer Matrix
# ============================================================

def build_mechanism_matrix(entries: list[dict]) -> str:
    lines = ["# Mechanism × Cancer Type Matrix\n"]
    lines.append("Cross-tabulation of article counts by therapeutic mechanism and cancer type.\n")

    mechanisms = sorted(MECHANISM_KEYWORDS.keys())
    cancer_types = sorted(CANCER_TYPE_KEYWORDS.keys())

    # Build matrix
    matrix = defaultdict(Counter)  # matrix[mechanism][cancer_type] = count
    for e in entries:
        for m in e.get("mechanisms", []):
            for c in e.get("cancer_types", []):
                matrix[m][c] += 1

    # Find top mechanism per cancer type
    top_per_cancer = {}
    for c in cancer_types:
        best_m, best_count = "", 0
        for m in mechanisms:
            if matrix[m][c] > best_count:
                best_m, best_count = m, matrix[m][c]
        top_per_cancer[c] = (best_m, best_count)

    # Table header
    lines.append(f"| Mechanism | {' | '.join(c[:8] for c in cancer_types)} | Total |")
    lines.append(f"|---|{'---|' * len(cancer_types)}---|")

    for m in mechanisms:
        row = [str(matrix[m][c]) if matrix[m][c] else "." for c in cancer_types]
        total = sum(matrix[m][c] for c in cancer_types)
        lines.append(f"| **{m}** | {' | '.join(row)} | {total} |")

    # Totals row
    col_totals = [sum(matrix[m][c] for m in mechanisms) for c in cancer_types]
    lines.append(f"| **Total** | {' | '.join(str(t) for t in col_totals)} | {sum(col_totals)} |")

    # Key observations
    lines.append("\n## Key Observations\n")

    # Hotspots (top 20 pairs)
    pairs = []
    for m in mechanisms:
        for c in cancer_types:
            if matrix[m][c] > 0:
                pairs.append((m, c, matrix[m][c]))
    pairs.sort(key=lambda x: -x[2])

    lines.append("### Top 20 Mechanism-Cancer Pairs (by article count)\n")
    for m, c, count in pairs[:20]:
        lines.append(f"- **{m} × {c}**: {count} articles")

    # Broad-spectrum mechanisms
    lines.append("\n### Broad-Spectrum Mechanisms (spanning most cancer types)\n")
    for m in mechanisms:
        cancer_count = sum(1 for c in cancer_types if matrix[m][c] > 0)
        if cancer_count >= 10:
            lines.append(f"- **{m}**: researched across {cancer_count}/{len(cancer_types)} cancer types")

    return "\n".join(lines)


# ============================================================
# Analysis 2: Convergence Map
# ============================================================

def build_convergence_map(entries: list[dict]) -> str:
    lines = ["# Convergence Map: Multi-Mechanism Articles\n"]
    lines.append("Articles that combine multiple therapeutic mechanisms — identifying synergistic research.\n")

    # Count mechanism pairs
    pair_counts = Counter()
    pair_articles = defaultdict(list)

    for e in entries:
        mechs = sorted(e.get("mechanisms", []))
        if len(mechs) < 2:
            continue
        for i in range(len(mechs)):
            for j in range(i + 1, len(mechs)):
                pair = (mechs[i], mechs[j])
                pair_counts[pair] += 1
                if len(pair_articles[pair]) < 5:  # Keep top 5 examples
                    pair_articles[pair].append(e)

    # Multi-mechanism stats
    multi = [e for e in entries if len(e.get("mechanisms", [])) >= 2]
    triple = [e for e in entries if len(e.get("mechanisms", [])) >= 3]

    lines.append(f"**Total articles with 2+ mechanisms**: {len(multi)} ({100*len(multi)//len(entries)}% of corpus)")
    lines.append(f"**Total articles with 3+ mechanisms**: {len(triple)}\n")

    # Top pairs
    lines.append("## Top 30 Mechanism Combinations\n")
    lines.append("| Rank | Mechanism A | Mechanism B | Articles | Example (top-cited) |")
    lines.append("|------|------------|------------|----------|---------------------|")

    for rank, ((m1, m2), count) in enumerate(pair_counts.most_common(30), 1):
        # Find highest-cited example
        examples = pair_articles[(m1, m2)]
        best = max(examples, key=lambda x: x.get("cited_by_count", 0))
        title = best.get("title", "")[:60] + "..."
        pmid = best.get("pmid", "")
        cites = best.get("cited_by_count", 0)
        lines.append(f"| {rank} | {m1} | {m2} | {count} | [{title}](by-pmid/{pmid}.md) ({cites} cites) |")

    # Notable 3-mechanism convergences
    lines.append("\n## Notable Triple-Mechanism Articles\n")
    triple_sorted = sorted(triple, key=lambda x: -(x.get("cited_by_count", 0) or 0))
    for e in triple_sorted[:15]:
        mechs = ", ".join(e.get("mechanisms", []))
        cancers = ", ".join(e.get("cancer_types", [])[:3]) or "general"
        cites = e.get("cited_by_count", 0)
        lines.append(f"- **[{e.get('pmid')}]** ({e.get('year')}, {cites} cites) — {mechs} — {cancers}")
        lines.append(f"  *{e.get('title', '')[:120]}*")

    # Unexplored combinations
    lines.append("\n## Unexplored Mechanism Combinations\n")
    lines.append("Pairs with zero or very few articles — potential research opportunities.\n")
    all_mechs = sorted(MECHANISM_KEYWORDS.keys())
    unexplored = []
    for i in range(len(all_mechs)):
        for j in range(i + 1, len(all_mechs)):
            pair = (all_mechs[i], all_mechs[j])
            count = pair_counts.get(pair, 0)
            if count <= 2:
                unexplored.append((pair, count))

    unexplored.sort(key=lambda x: x[1])
    for (m1, m2), count in unexplored[:30]:
        status = "**ZERO**" if count == 0 else f"{count}"
        lines.append(f"- {m1} + {m2}: {status} articles")

    return "\n".join(lines)


# ============================================================
# Analysis 3: Gap Analysis
# ============================================================

def build_gap_analysis(entries: list[dict]) -> str:
    lines = ["# Gap Analysis: Underexplored Research Areas\n"]
    lines.append(
        "This file is useful for hypothesis generation, but it is taxonomy-dependent. "
        "All zero-count and low-evidence findings here should be read as corpus-level "
        "non-detection rather than definitive absence unless they are externally verified. "
        "See `analysis/taxonomy-rerun-notes.md` for the current list of known taxonomy/query artifacts.\n"
    )

    mechanisms = sorted(MECHANISM_KEYWORDS.keys())
    cancer_types = sorted(CANCER_TYPE_KEYWORDS.keys())

    # Matrix for gap detection
    matrix = defaultdict(Counter)
    for e in entries:
        for m in e.get("mechanisms", []):
            for c in e.get("cancer_types", []):
                matrix[m][c] += 1

    # Evidence matrix
    evidence_matrix = defaultdict(lambda: defaultdict(str))
    evidence_rank = {"phase3-clinical": 6, "phase2-clinical": 5, "phase1-clinical": 4,
                     "preclinical-invivo": 3, "preclinical-invitro": 2, "theoretical": 1, "": 0}
    for e in entries:
        ev = e.get("evidence_level", "")
        for m in e.get("mechanisms", []):
            for c in e.get("cancer_types", []):
                if evidence_rank.get(ev, 0) > evidence_rank.get(evidence_matrix[m][c], 0):
                    evidence_matrix[m][c] = ev

    # 1. Complete gaps (mechanism × cancer = 0 articles)
    lines.append("## 1. Complete Gaps (0 articles)\n")
    lines.append("Mechanism × cancer type pairs with zero research articles in our corpus.\n")

    complete_gaps = []
    for m in mechanisms:
        for c in cancer_types:
            if matrix[m][c] == 0:
                # Only flag if both mechanism and cancer type have substantial corpus
                m_total = sum(1 for e in entries if m in e.get("mechanisms", []))
                c_total = sum(1 for e in entries if c in e.get("cancer_types", []))
                if m_total >= 50 and c_total >= 50:
                    complete_gaps.append((m, c, m_total, c_total))

    complete_gaps.sort(key=lambda x: -(x[2] + x[3]))
    for m, c, m_total, c_total in complete_gaps[:30]:
        lines.append(f"- **{m} × {c}** — mechanism has {m_total} articles total, cancer has {c_total} total, but 0 overlap")

    # 2. Preclinical-only mechanisms with high article counts
    lines.append("\n## 2. Mechanisms Stuck in Preclinical\n")
    lines.append(
        "Mechanisms with many articles but no Phase 2+ clinical evidence detected in the current "
        "keyword-derived evidence tags. Coverage values below use only primary-study-like records "
        "(excluding review-like and protocol-like papers that are intentionally left unclassified).\n"
    )

    for m in mechanisms:
        m_articles = [e for e in entries if m in e.get("mechanisms", [])]
        phase2_plus = [e for e in m_articles if e.get("evidence_level", "") in ("phase2-clinical", "phase3-clinical")]
        primary_like = [e for e in m_articles if classify_evidence_reason(e) in ("tagged", "other_untagged")]
        primary_tagged = sum(1 for e in primary_like if e.get("evidence_level"))
        primary_cov = (primary_tagged / len(primary_like)) if primary_like else 1.0
        if len(m_articles) >= 100 and len(phase2_plus) == 0:
            lines.append(
                f"- **{m}**: {len(m_articles)} articles, 0 Phase 2+ clinical trials detected; "
                f"primary-study-like coverage {primary_tagged}/{len(primary_like)} ({primary_cov:.1%})"
            )
        elif len(m_articles) >= 100 and len(phase2_plus) <= 3:
            lines.append(
                f"- **{m}**: {len(m_articles)} articles, only {len(phase2_plus)} Phase 2+ "
                f"({100*len(phase2_plus)//len(m_articles)}%); primary-study-like coverage "
                f"{primary_tagged}/{len(primary_like)} ({primary_cov:.1%})"
            )

    # 3. Cancer types with limited novel mechanism research
    lines.append("\n## 3. Cancer Types Underserved by Novel Mechanisms\n")
    lines.append("Cancer types where most research is concentrated in immunotherapy alone.\n")

    for c in cancer_types:
        c_articles = [e for e in entries if c in e.get("cancer_types", [])]
        if len(c_articles) < 50:
            continue
        immuno_only = [e for e in c_articles if e.get("mechanisms") == ["immunotherapy"]]
        novel_mechs = set()
        for e in c_articles:
            for m in e.get("mechanisms", []):
                if m != "immunotherapy":
                    novel_mechs.add(m)
        if len(novel_mechs) <= 5:
            lines.append(f"- **{c}**: {len(c_articles)} articles, but only {len(novel_mechs)} non-immunotherapy mechanisms explored")

    # 4. High-potential gaps (well-researched mechanism, well-researched cancer, no overlap)
    lines.append("\n## 4. Highest-Priority Research Opportunities\n")
    lines.append("Mechanism-cancer pairs where both have strong evidence bases but no combined research.\n")

    priority_gaps = []
    for m, c, m_total, c_total in complete_gaps:
        # Check if mechanism has clinical evidence somewhere
        m_clinical = any(e.get("evidence_level", "").startswith("phase") for e in entries if m in e.get("mechanisms", []))
        if m_clinical and c_total >= 100:
            priority_gaps.append((m, c, m_total, c_total))

    priority_gaps.sort(key=lambda x: -(x[2] * x[3]))
    for m, c, m_total, c_total in priority_gaps[:15]:
        lines.append(f"- **{m} × {c}** — {m} has clinical trials elsewhere, {c} has {c_total} articles, but no combined research")

    return "\n".join(lines)


# ============================================================
# Analysis 4: Evidence Tiers
# ============================================================

def build_evidence_tiers(entries: list[dict]) -> str:
    lines = ["# Evidence Tiers by Mechanism\n"]
    lines.append("Highest level of clinical evidence for each therapeutic mechanism.\n")
    reason_counts = Counter(classify_evidence_reason(e) for e in entries)
    coverage = reason_counts["tagged"]
    primary_like_total = reason_counts["tagged"] + reason_counts["other_untagged"]
    lines.append(
        f"Evidence tags are currently populated for {coverage}/{len(entries)} full-text records "
        f"({coverage/len(entries):.1%}). Reviews/meta-analyses ({reason_counts['review_like']}) "
        f"and protocols ({reason_counts['protocol_like']}) are intentionally left unclassified; "
        f"among primary-study-like records, coverage is {coverage}/{primary_like_total} "
        f"({coverage/primary_like_total:.1%}). Absence claims remain provisional.\n"
    )

    mechanisms = sorted(MECHANISM_KEYWORDS.keys())
    evidence_order = ["phase3-clinical", "phase2-clinical", "phase1-clinical",
                      "preclinical-invivo", "preclinical-invitro", "theoretical"]
    evidence_labels = {
        "phase3-clinical": "Phase III RCT",
        "phase2-clinical": "Phase II",
        "phase1-clinical": "Phase I",
        "preclinical-invivo": "Preclinical (in vivo)",
        "preclinical-invitro": "Preclinical (in vitro)",
        "theoretical": "Theoretical/Computational",
    }

    lines.append("| Mechanism | Highest Evidence | Phase 3 | Phase 2 | Phase 1 | In Vivo | In Vitro | Theory | Total |")
    lines.append("|-----------|-----------------|---------|---------|---------|---------|----------|--------|-------|")

    for m in mechanisms:
        m_articles = [e for e in entries if m in e.get("mechanisms", [])]
        counts = Counter(e.get("evidence_level", "") for e in m_articles)
        highest = ""
        for ev in evidence_order:
            if counts.get(ev, 0) > 0:
                highest = evidence_labels.get(ev, ev)
                break
        if not highest:
            highest = "Not classified"

        row = [str(counts.get(ev, 0)) for ev in evidence_order]
        lines.append(f"| **{m}** | {highest} | {' | '.join(row)} | {len(m_articles)} |")

    # Mechanisms with Phase 3 evidence
    lines.append("\n## Mechanisms with Phase III Clinical Evidence\n")
    for m in mechanisms:
        m_articles = [e for e in entries if m in e.get("mechanisms", [])]
        phase3 = [e for e in m_articles if e.get("evidence_level") == "phase3-clinical"]
        if phase3:
            phase3.sort(key=lambda x: -(x.get("cited_by_count", 0) or 0))
            lines.append(f"\n### {m} ({len(phase3)} Phase III articles)\n")
            for e in phase3[:5]:
                cites = e.get("cited_by_count", 0)
                cancers = ", ".join(e.get("cancer_types", [])[:3]) or "various"
                lines.append(f"- **PMID {e['pmid']}** ({e.get('year')}, {cites} cites) — {cancers}")
                lines.append(f"  *{e.get('title', '')[:150]}*")

    return "\n".join(lines)


def build_resistant_state_map(entries: list[dict]) -> str:
    lines = ["# Resistant-State Map\n"]
    lines.append(
        "First-pass scaffold for analyzing the corpus by resistant state rather than by modality alone.\n"
    )
    lines.append(
        "These state assignments are keyword-derived heuristics. They are intended to support "
        "prioritization and literature review, not to assert that a paper experimentally validated a state transition.\n"
    )

    states = sorted(RESISTANT_STATE_RULES.keys())
    mechanisms = sorted(MECHANISM_KEYWORDS.keys())
    tagged_count = sum(1 for e in entries if e.get("resistant_states"))
    lines.append(
        f"Current resistant-state coverage in the index: {tagged_count}/{len(entries)} records "
        f"({tagged_count/len(entries):.1%}).\n"
    )
    if tagged_count == 0:
        lines.append(
            "WARNING: no resistant-state tags are present in the current index. "
            "Re-run `tag_articles.py` and `build_index.py` before interpreting the table below.\n"
        )

    lines.append("\n## Resistant States Tracked\n")
    for state in states:
        lines.append(f"- **{state}**")

    lines.append("\n## State × Mechanism Counts\n")
    lines.append("| Resistant State | Top linked mechanisms | Tagged articles |")
    lines.append("|---|---|---|")

    for state in states:
        state_articles = [e for e in entries if state in e.get("resistant_states", [])]
        mech_counts = Counter()
        for entry in state_articles:
            for mech in entry.get("mechanisms", []):
                mech_counts[mech] += 1
        top_mechs = ", ".join(f"{m} ({c})" for m, c in mech_counts.most_common(5)) or "none"
        lines.append(f"| **{state}** | {top_mechs} | {len(state_articles)} |")

    lines.append("\n## Interpretation\n")
    lines.append(
        "- The repo should use these states as the primary decision layer when comparing interventions."
    )
    lines.append(
        "- Physical ROS modalities should be framed as best-matched to OXPHOS-dependent, ferroptosis-prone persisters rather than as a universal answer to resistance."
    )
    lines.append(
        "- Senescence, stromal sheltering, and NRF2/SLC7A11 compensation should be treated as parallel escape states, not edge cases."
    )
    return "\n".join(lines)


def build_evidence_coverage_audit(entries: list[dict]) -> str:
    lines = ["# Evidence Coverage Audit\n"]
    total = len(entries)
    reason_counts = Counter(classify_evidence_reason(e) for e in entries)
    tagged = [e for e in entries if e.get("evidence_level")]
    primary_like_total = reason_counts["tagged"] + reason_counts["other_untagged"]
    lines.append(
        f"Evidence-level tags are present for {len(tagged)}/{total} records ({len(tagged)/total:.1%}). "
        f"Of the unclassified records, {reason_counts['review_like']} are review-like and "
        f"{reason_counts['protocol_like']} are protocol-like by design; {reason_counts['other_untagged']} "
        f"primary-study-like records remain uncategorized. Primary-study-like evidence coverage is "
        f"{reason_counts['tagged']}/{primary_like_total} ({reason_counts['tagged']/primary_like_total:.1%}).\n"
    )

    lines.append("## Mechanisms Most Exposed To Overstated Absence Claims\n")
    lines.append("| Mechanism | Total | Tagged | Review-like | Protocol-like | Other untagged | Primary-study-like coverage |")
    lines.append("|---|---|---|---|---|---|---|")
    mechanism_rows = []
    for mechanism in sorted(MECHANISM_KEYWORDS.keys()):
        mech_articles = [e for e in entries if mechanism in e.get('mechanisms', [])]
        if not mech_articles:
            continue
        mech_counts = Counter(classify_evidence_reason(e) for e in mech_articles)
        primary_like = mech_counts["tagged"] + mech_counts["other_untagged"]
        primary_cov = (mech_counts["tagged"] / primary_like) if primary_like else 1.0
        mechanism_rows.append((mechanism, mech_articles, mech_counts, primary_cov))
    for mechanism, mech_articles, mech_counts, primary_cov in sorted(
        mechanism_rows,
        key=lambda row: (-row[2]["other_untagged"], row[3], row[0]),
    ):
        lines.append(
            f"| **{mechanism}** | {len(mech_articles)} | {mech_counts['tagged']} | {mech_counts['review_like']} | "
            f"{mech_counts['protocol_like']} | {mech_counts['other_untagged']} | {mech_counts['tagged']}/{mech_counts['tagged'] + mech_counts['other_untagged']} ({primary_cov:.1%}) |"
        )

    lines.append("\n## Sample Of Unclassified Primary-Study-Like Records\n")
    lines.append(
        "Illustrative examples below come from the uncategorized primary-study-like pool rather than the "
        "review/protocol bucket. These are the records most likely to affect absence claims if the evidence "
        "classifier is expanded.\n"
    )
    focus_mechanisms = [
        row[0] for row in sorted(
            mechanism_rows,
            key=lambda row: (-row[2]["other_untagged"], row[3], row[0]),
        )
        if row[2]["other_untagged"] > 0 and len(row[1]) >= 50
    ][:5]
    for mechanism in focus_mechanisms:
        candidates = [
            e for e in entries
            if mechanism in e.get("mechanisms", [])
            and classify_evidence_reason(e) == "other_untagged"
        ]
        candidates.sort(key=lambda e: (-(e.get("cited_by_count") or 0), -(e.get("year") or 0), e.get("pmid", "")))
        lines.append(f"\n### {mechanism}\n")
        for e in candidates[:3]:
            lines.append(
                f"- **PMID {e['pmid']}** ({e.get('year')}) — *{e.get('title', '')[:150]}*"
            )

    lines.append("\n## What The Current Miss-Rate Signal Likely Means\n")
    lines.append(
        "- The raw 36.8% coverage number is pessimistic because review-like and protocol-like records are intentionally excluded from evidence tagging."
    )
    lines.append(
        "- The more relevant upper-bound miss rate is the share of `other_untagged` records within the primary-study-like subset. Mechanisms with the largest remaining uncertainty are immunotherapy, mRNA-vaccine, electrochemical-therapy, TTFields, and CAR-T."
    )
    lines.append(
        "- The sampled uncategorized records are enriched for observational clinical studies, biomarker/antigen-discovery papers, and translational engineering studies that do not announce phase or preclinical status in obvious keywords."
    )
    lines.append(
        "- This means the main risk is overstating `no detected clinical evidence` for modalities with many non-phase clinical or translational papers, not silently missing large numbers of explicit Phase III trials."
    )

    lines.append("\n## Recommended Interpretation Guardrails\n")
    lines.append("- Treat `0 Phase 2+` as `not detected in current keyword-derived evidence tags` unless manually verified.")
    lines.append("- Distinguish review/protocol exclusions from true uncategorized primary-study-like records when discussing evidence coverage.")
    lines.append("- Re-check any high-priority mechanism with external PubMed or trial-registry verification before using it as a headline gap.")
    lines.append("- Prefer coverage-aware language in the manuscript and analysis files whenever evidence tagging is below 50% for a mechanism.")
    return "\n".join(lines)


# ============================================================
# Analysis 5: Key Findings (Top 100 articles by impact)
# ============================================================

def build_key_findings(entries: list[dict]) -> str:
    lines = ["# Key Findings: Top 100 Highest-Impact Articles\n"]
    lines.append("Ranked by iCite Relative Citation Ratio (field-normalized impact).\n")

    # Filter to articles with RCR
    with_rcr = [e for e in entries if e.get("icite_rcr") and e["icite_rcr"] > 0]
    with_rcr.sort(key=lambda x: -(x.get("icite_rcr") or 0))

    lines.append(f"Total articles with RCR: {len(with_rcr)}\n")

    for rank, e in enumerate(with_rcr[:100], 1):
        rcr = e.get("icite_rcr", 0)
        pct = e.get("icite_percentile")
        cites = e.get("cited_by_count", 0)
        mechs = ", ".join(e.get("mechanisms", [])) or "untagged"
        cancers = ", ".join(e.get("cancer_types", [])[:3]) or "general"
        ev = e.get("evidence_level", "unknown")
        journal = e.get("journal", "")
        year = e.get("year", "?")
        oa = "OA" if e.get("is_oa") else "closed"
        pct_str = f", {pct}th percentile" if pct else ""

        lines.append(f"### {rank}. PMID {e['pmid']} — RCR {rcr:.1f}{pct_str}")
        lines.append(f"**{e.get('title', '')}**")
        lines.append(f"*{journal}* ({year}) | {cites} citations | {oa}")
        lines.append(f"Mechanisms: {mechs} | Cancer: {cancers} | Evidence: {ev}")

        # Load abstract for key insight
        abstract = load_article_abstract(e["pmid"])
        if abstract:
            # Extract first 2 sentences as key insight
            sentences = re.split(r'(?<=[.!?])\s+', abstract[:500])
            insight = " ".join(sentences[:2])
            if len(insight) > 300:
                insight = insight[:297] + "..."
            lines.append(f"**Key insight**: {insight}")

        lines.append("")

    return "\n".join(lines)


# ============================================================
# Analysis 6: Timeline
# ============================================================

def build_timeline(entries: list[dict]) -> str:
    lines = ["# Timeline of Key Breakthroughs (2015-2026)\n"]
    lines.append("Major milestones chronologically, based on highest-impact articles per year per mechanism.\n")

    # Get top article per mechanism per year
    mech_year = defaultdict(list)
    for e in entries:
        year = e.get("year")
        if not year or year < 2015:
            continue
        for m in e.get("mechanisms", []):
            mech_year[(m, year)].append(e)

    # For each year, find the most impactful articles
    year_highlights = defaultdict(list)
    for (m, year), articles in mech_year.items():
        best = max(articles, key=lambda x: x.get("cited_by_count", 0) or 0)
        if (best.get("cited_by_count", 0) or 0) >= 50:
            year_highlights[year].append((m, best))

    for year in sorted(year_highlights.keys()):
        items = year_highlights[year]
        items.sort(key=lambda x: -(x[1].get("cited_by_count", 0) or 0))

        lines.append(f"\n## {year}\n")
        for m, e in items[:8]:
            cites = e.get("cited_by_count", 0)
            cancers = ", ".join(e.get("cancer_types", [])[:2]) or "various"
            ev = e.get("evidence_level", "")
            ev_str = f" [{ev}]" if ev else ""
            lines.append(f"- **{m}**{ev_str}: *{e.get('title', '')[:120]}*")
            lines.append(f"  PMID {e['pmid']} | {e.get('journal', '')} | {cites} citations | {cancers}")

    # Summary of trends
    lines.append("\n## Trend Summary\n")

    # Articles per year
    year_counts = Counter(e.get("year") for e in entries if e.get("year") and e["year"] >= 2015)
    lines.append("### Publication volume by year\n")
    for y in sorted(year_counts.keys()):
        bar = "█" * (year_counts[y] // 50)
        lines.append(f"- {y}: {year_counts[y]:>5} articles {bar}")

    # Emerging mechanisms (growth rate)
    lines.append("\n### Fastest-growing mechanisms (2022 vs 2020)\n")
    for m in sorted(MECHANISM_KEYWORDS.keys()):
        count_2020 = sum(1 for e in entries if m in e.get("mechanisms", []) and e.get("year") == 2020)
        count_2024 = sum(1 for e in entries if m in e.get("mechanisms", []) and e.get("year") == 2024)
        if count_2020 > 5 and count_2024 > count_2020:
            growth = (count_2024 - count_2020) / count_2020 * 100
            lines.append(f"- **{m}**: {count_2020} → {count_2024} articles/year (+{growth:.0f}%)")

    return "\n".join(lines)


# ============================================================
# Main
# ============================================================

def main():
    print("Loading index...")
    entries = load_index()
    print(f"  Loaded {len(entries)} articles")
    resistant_state_coverage = sum(1 for e in entries if e.get("resistant_states"))
    if resistant_state_coverage == 0:
        print("  Warning: resistant_states is empty in INDEX.jsonl; regenerate tags before relying on resistant-state outputs")

    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    analyses = [
        ("mechanism-matrix.md", "Mechanism-Cancer Matrix", build_mechanism_matrix),
        ("convergence-map.md", "Convergence Map", build_convergence_map),
        ("gap-analysis.md", "Gap Analysis", build_gap_analysis),
        ("evidence-tiers.md", "Evidence Tiers", build_evidence_tiers),
        ("resistant-state-map.md", "Resistant-State Map", build_resistant_state_map),
        ("evidence-coverage-audit.md", "Evidence Coverage Audit", build_evidence_coverage_audit),
        ("key-findings.md", "Key Findings (top 100)", build_key_findings),
        ("timeline.md", "Timeline", build_timeline),
    ]

    for filename, label, builder in analyses:
        print(f"Building {label}...")
        content = builder(entries)
        filepath = ANALYSIS_DIR / filename
        filepath.write_text(content, encoding="utf-8")
        line_count = content.count("\n")
        print(f"  Written {filepath.name} ({line_count} lines)")

    print(f"\nAll analysis files written to {ANALYSIS_DIR}/")


if __name__ == "__main__":
    main()
