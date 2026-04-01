#!/usr/bin/env python3
"""
Auto-tag articles by mechanism, cancer type, and evidence level.
Reads title + abstract + MeSH terms, matches against keyword dictionaries.

Usage:
    python tag_articles.py          # Tag all articles, rebuild tag index files
    python tag_articles.py --dry-run  # Show what would be tagged without writing
"""

import argparse
import re
import sys
from pathlib import Path

from tqdm import tqdm

from article_io import load_article, save_article
from config import (
    BIOLOGY_PROCESS_KEYWORDS, CANCER_TYPE_KEYWORDS, EVIDENCE_LEVEL_KEYWORDS, MECHANISM_KEYWORDS,
    PATHWAY_TARGET_KEYWORDS, PMID_DIR, RADIOLIGAND_TARGET_KEYWORDS, TAGS_DIR,
    RESISTANT_STATE_RULES,
)
from evidence_utils import is_protocol_like, is_review_like, normalize_text

GENERIC_CANCER_TERMS = (
    "cancer", "neoplasm", "carcinoma", "tumor", "tumour",
    "oncology", "malign", "glioblastoma", "melanoma", "sarcoma",
    "leukemia", "lymphoma", "myeloma",
)

EVIDENCE_PUBTYPE_MARKERS = {
    "phase3-clinical": ("clinical trial, phase iii", "clinical trial, phase 3", "clinical trial, phase iv", "clinical trial, phase 4"),
    "phase2-clinical": ("clinical trial, phase ii", "clinical trial, phase 2"),
    "phase1-clinical": ("clinical trial, phase i", "clinical trial, phase 1", "clinical trial, phase i/ii"),
    "clinical-other": ("clinical trial", "controlled clinical trial", "observational study", "case reports"),
}

MRNA_VACCINE_PLATFORM_TERMS = (
    "mrna vaccine", "mrna vaccines", "mrna cancer vaccine",
    "messenger rna vaccine", "messenger rna vaccines",
    "mrna-based vaccine", "mrna-based vaccines",
    "rna neoantigen vaccine", "rna neoantigen vaccines",
    "personalized rna vaccine", "individualized rna vaccine",
    "individualized mrna vaccine", "individualized mrna vaccines",
)

MRNA_VACCINE_THERAPEUTIC_TERMS = (
    "neoantigen", "tumor antigen", "tumour antigen",
    "cancer vaccine", "cancer vaccines",
    "therapeutic vaccination", "therapeutic vaccine", "therapeutic vaccines",
    "cancer treatment",
    "tumor-specific antigen", "tumour-specific antigen",
)

MRNA_VACCINE_STRONG_THERAPEUTIC_PHRASES = (
    "mrna cancer vaccine",
    "neoantigen vaccine", "neoantigen vaccines",
    "rna neoantigen vaccine", "rna neoantigen vaccines",
    "individualized mrna vaccine", "individualized mrna vaccines",
    "mrna vaccines for cancer treatment",
    "mrna vaccination in breast cancer",
    "mRNA-based precision targeting of neoantigens".lower(),
)

MRNA_VACCINE_EXCLUSION_TERMS = (
    "covid", "sars-cov-2", "coronavirus", "omicron",
    "bnt162", "mrna-1273", "moderna", "pfizer", "booster",
    "infectious disease", "viral infection", "bacterial infection",
    "pseudomonas", "rsv vaccine", "influenza vaccine",
)

RADIOLIGAND_CORE_TERMS = (
    "radioligand therapy", "radiopharmaceutical therapy", "radionuclide therapy",
    "targeted radionuclide therapy", "targeted radioligand therapy",
    "peptide receptor radionuclide therapy", "prrt", "radioiodine therapy",
    "psma radioligand", "radioligand therapeutic", "targeted radionuclide",
)

RADIOLIGAND_SUPPORT_TERMS = (
    "radiopharmaceutical", "radionuclide", "radioligand", "radiolabeled",
    "radiolabelled", "radioisotope", "radioisotopic",
    "theranostic", "theranostics",
)

RADIOLIGAND_THERAPY_TERMS = (
    "therapy", "therapeutic", "treatment", "treated", "dose escalation",
    "phase i", "phase ii", "phase iii", "theranostic agent", "theranostic agents",
)

RADIOLIGAND_ISOTOPE_TERMS = (
    "lutetium-177", "lutetium 177", "lu-177", "177lu",
    "actinium-225", "actinium 225", "ac-225", "225ac",
    "yttrium-90", "yttrium 90", "y-90", "90y",
    "iodine-131", "iodine 131", "i-131", "131i",
    "radium-223", "radium 223", "223ra",
    "terbium-161", "terbium 161", "tb-161",
    "lutathera", "pluvicto", "xofigo",
    "vipivotide tetraxetan", "dotatate", "dotatoc",
)

COMBINATION_LANGUAGE_TERMS = (
    "combination therapy", "combination treatment", "combination regimen",
    "in combination with", "combined with", "combined therapy", "combined treatment",
    "combined regimen", "combined strategy", "combined modality", "combination of",
    "co-treatment", "cotherapy", "together with", "synergizes with", "synergy with",
    "synergistic", "added to", "augment",
)


def has_cancer_context(text: str) -> bool:
    """Require generic cancer context for broad mechanism tags."""
    return any(term in text for term in GENERIC_CANCER_TERMS)


def get_searchable_text(fm: dict, body: str, include_full_text: bool = False) -> str:
    """Combine title, abstract, metadata, and optionally full text into searchable text."""
    parts = [
        fm.get("title", ""),
        " ".join(fm.get("mesh_terms", [])),
        " ".join(fm.get("diseases_annotated", [])),
        " ".join(fm.get("genes", [])),
        " ".join(fm.get("drugs", [])),
    ]

    # Extract abstract from body (between ## Abstract and next ##)
    abstract_match = re.search(r"## Abstract\n\n?(.*?)(?=\n## |\Z)", body, re.DOTALL)
    if abstract_match:
        parts.append(abstract_match.group(1))
    if include_full_text:
        parts.append(body)

    return normalize_text(" ".join(parts))


def match_keywords(text: str, keyword_dict: dict) -> list[str]:
    """Match text against keyword dictionary using word-boundary matching.

    Short keywords (<=4 chars) use strict word-boundary regex to avoid
    false positives like 'all' matching 'overall'.
    """
    matched = []
    for tag, keywords in keyword_dict.items():
        for kw in keywords:
            kw_lower = kw.lower()
            if len(kw_lower) <= 4:
                # Strict word-boundary match for short keywords
                if re.search(r'\b' + re.escape(kw_lower) + r'\b', text):
                    matched.append(tag)
                    break
            else:
                if kw_lower in text:
                    matched.append(tag)
                    break
    return sorted(matched)


def text_matches_keyword(text: str, keyword: str) -> bool:
    """Return True when a keyword is present in text with stricter handling for short terms."""
    kw_lower = keyword.lower()
    if len(kw_lower) <= 4:
        return bool(re.search(r'\b' + re.escape(kw_lower) + r'\b', text))
    return kw_lower in text


def match_resistant_states(text: str) -> list[str]:
    """Match resistant states using composite rules rather than single-keyword OR logic."""
    matched = []
    for state, rule in RESISTANT_STATE_RULES.items():
        all_groups = rule.get("all_of", [])
        if all(any(text_matches_keyword(text, kw) for kw in group) for group in all_groups):
            matched.append(state)
    return sorted(matched)


def match_mrna_vaccine(text: str, title_text: str) -> bool:
    """Match therapeutic cancer-vaccine mRNA studies while excluding supportive/infectious vaccine papers."""
    has_platform = any(term in text for term in MRNA_VACCINE_PLATFORM_TERMS)
    has_therapeutic = any(term in text for term in MRNA_VACCINE_THERAPEUTIC_TERMS)
    has_strong_phrase = any(term in text for term in MRNA_VACCINE_STRONG_THERAPEUTIC_PHRASES)
    has_exclusion = any(term in title_text for term in MRNA_VACCINE_EXCLUSION_TERMS)

    if has_strong_phrase:
        return True
    if not has_platform or not has_therapeutic:
        return False
    return not has_exclusion


def match_radioligand_therapy(text: str, title_text: str) -> bool:
    """Require radionuclide-specific therapy signals instead of generic theranostic language."""
    title_has_core = any(term in title_text for term in RADIOLIGAND_CORE_TERMS)
    title_has_support = any(term in title_text for term in RADIOLIGAND_SUPPORT_TERMS)
    title_has_isotope = any(term in title_text for term in RADIOLIGAND_ISOTOPE_TERMS)
    title_has_target = any(
        text_matches_keyword(title_text, kw)
        for keywords in RADIOLIGAND_TARGET_KEYWORDS.values()
        for kw in keywords
    )
    title_has_therapy = any(term in title_text for term in RADIOLIGAND_THERAPY_TERMS)
    has_core = any(term in text for term in RADIOLIGAND_CORE_TERMS)
    has_isotope = any(term in text for term in RADIOLIGAND_ISOTOPE_TERMS)
    has_target = any(
        text_matches_keyword(text, kw)
        for keywords in RADIOLIGAND_TARGET_KEYWORDS.values()
        for kw in keywords
    )
    has_support = any(term in text for term in RADIOLIGAND_SUPPORT_TERMS)

    if title_has_core or title_has_isotope:
        return True
    if title_has_support and (title_has_target or title_has_isotope):
        return True
    if has_core and has_isotope:
        return True
    return False


def match_radioligand_targets(text: str, mechanisms: list[str]) -> list[str]:
    """Tag explicit radioligand targets only inside the cleaned radioligand lane."""
    if "radioligand-therapy" not in mechanisms:
        return []
    return match_keywords(text, RADIOLIGAND_TARGET_KEYWORDS)


def classify_combination_evidence(fm: dict, title_text: str, abstract_text: str, mechanisms: list[str], evidence: str) -> str:
    """Separate broad co-mentions from deliberate combination studies."""
    if len(mechanisms) < 2:
        return ""
    if is_review_like(fm) or is_protocol_like(fm):
        return "review-or-perspective-multi-lane"

    combo_text = f"{title_text} {abstract_text}"
    has_combo_language = any(term in combo_text for term in COMBINATION_LANGUAGE_TERMS)
    if not has_combo_language:
        return "co-mention-only"

    if evidence in ("phase3-clinical", "phase2-clinical", "phase1-clinical", "clinical-other"):
        return "designed-combination-clinical"
    if evidence in ("preclinical-invivo", "preclinical-invitro"):
        return "designed-combination-preclinical"
    return "co-mention-only"


def match_mechanisms(text: str, title_text: str) -> list[str]:
    """Match mechanisms with a coarse cancer-context gate to reduce off-target tags."""
    if not has_cancer_context(text):
        return []
    matched = set(match_keywords(text, MECHANISM_KEYWORDS))
    if match_mrna_vaccine(text, title_text):
        matched.add("mRNA-vaccine")
    else:
        matched.discard("mRNA-vaccine")
    if match_radioligand_therapy(text, title_text):
        matched.add("radioligand-therapy")
    else:
        matched.discard("radioligand-therapy")
    return sorted(matched)


def match_evidence_level(fm: dict, text: str) -> str:
    """Match text against evidence level keywords, return best match.

    Uses word-boundary matching for short keywords to avoid false positives.
    Priority order: phase3 > phase2 > phase1 > clinical-other > invivo > invitro > theoretical.
    """
    if is_review_like(fm) or is_protocol_like(fm):
        return ""

    pub_types = [normalize_text(p) for p in fm.get("pub_types", [])]
    for level in ["phase3-clinical", "phase2-clinical", "phase1-clinical"]:
        if any(marker in pub_types for marker in EVIDENCE_PUBTYPE_MARKERS[level]):
            return level

    if any(marker in pub_types for marker in EVIDENCE_PUBTYPE_MARKERS["clinical-other"]):
        return "clinical-other"

    for level in ["phase3-clinical", "phase2-clinical", "phase1-clinical", "clinical-other",
                   "preclinical-invivo", "preclinical-invitro", "theoretical"]:
        for kw in EVIDENCE_LEVEL_KEYWORDS[level]:
            kw_lower = kw.lower()
            if len(kw_lower) <= 4:
                if re.search(r'\b' + re.escape(kw_lower) + r'\b', text):
                    return level
            else:
                if kw_lower in text:
                    return level
    return ""


def write_tag_files(tag_type: str, tag_pmids: dict[str, list[str]]) -> None:
    """Write tag index files. Each file contains one PMID per line."""
    tag_dir = TAGS_DIR / tag_type
    tag_dir.mkdir(parents=True, exist_ok=True)

    # Clear existing tag files in this directory
    for f in tag_dir.glob("*.txt"):
        f.unlink()

    for tag, pmids in sorted(tag_pmids.items()):
        if pmids:
            filepath = tag_dir / f"{tag}.txt"
            filepath.write_text("\n".join(sorted(pmids)) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Auto-tag articles and build tag indexes")
    parser.add_argument("--dry-run", action="store_true", help="Show tags without writing")
    args = parser.parse_args()

    files = sorted(PMID_DIR.glob("*.md"))
    if not files:
        print("No articles found in corpus.")
        return

    print(f"Tagging {len(files)} articles...")

    # Accumulators for tag index files
    mechanism_pmids: dict[str, list[str]] = {k: [] for k in MECHANISM_KEYWORDS}
    biology_process_pmids: dict[str, list[str]] = {k: [] for k in BIOLOGY_PROCESS_KEYWORDS}
    cancer_pmids: dict[str, list[str]] = {k: [] for k in CANCER_TYPE_KEYWORDS}
    evidence_pmids: dict[str, list[str]] = {k: [] for k in EVIDENCE_LEVEL_KEYWORDS}
    resistant_state_pmids: dict[str, list[str]] = {k: [] for k in RESISTANT_STATE_RULES}
    radioligand_target_pmids: dict[str, list[str]] = {k: [] for k in RADIOLIGAND_TARGET_KEYWORDS}
    combination_pmids: dict[str, list[str]] = {
        "co-mention-only": [],
        "designed-combination-preclinical": [],
        "designed-combination-clinical": [],
        "review-or-perspective-multi-lane": [],
    }
    journal_pmids: dict[str, list[str]] = {}

    pathway_target_pmids: dict[str, list[str]] = {k: [] for k in PATHWAY_TARGET_KEYWORDS}

    stats = {
        "mechanisms": 0, "biology_processes": 0, "pathway_targets": 0,
        "cancer_types": 0, "evidence": 0, "resistant_states": 0,
        "radioligand_targets": 0, "combination_evidence": 0,
    }

    for filepath in tqdm(files, desc="  Tagging"):
        fm, body = load_article(filepath)
        if not fm:
            continue

        pmid = fm.get("pmid", filepath.stem)
        text = get_searchable_text(fm, body)
        title_text = normalize_text(fm.get("title", ""))
        pathway_text = get_searchable_text(fm, body, include_full_text=True)
        abstract_match = re.search(r"## Abstract\n\n?(.*?)(?=\n## |\Z)", body, re.DOTALL)
        abstract_text = normalize_text(abstract_match.group(1) if abstract_match else "")

        # Match
        mechanisms = match_mechanisms(text, title_text)
        biology_processes = match_keywords(text, BIOLOGY_PROCESS_KEYWORDS)
        pathway_targets = match_keywords(pathway_text, PATHWAY_TARGET_KEYWORDS)
        cancer_types = match_keywords(text, CANCER_TYPE_KEYWORDS)
        evidence = match_evidence_level(fm, text)
        resistant_states = match_resistant_states(text)
        radioligand_targets = match_radioligand_targets(pathway_text, mechanisms)
        combination_evidence = classify_combination_evidence(fm, title_text, abstract_text, mechanisms, evidence)

        # Update frontmatter
        fm["mechanisms"] = mechanisms
        fm["biology_processes"] = biology_processes
        fm["pathway_targets"] = pathway_targets
        fm["cancer_types"] = cancer_types
        fm["evidence_level"] = evidence
        fm["resistant_states"] = resistant_states
        if radioligand_targets:
            fm["radioligand_targets"] = radioligand_targets
        else:
            fm.pop("radioligand_targets", None)
        if combination_evidence:
            fm["combination_evidence"] = combination_evidence
        else:
            fm.pop("combination_evidence", None)

        if not args.dry_run:
            save_article(filepath, fm, body)

        # Accumulate for tag files
        for m in mechanisms:
            mechanism_pmids[m].append(pmid)
        for b in biology_processes:
            biology_process_pmids[b].append(pmid)
        for p in pathway_targets:
            pathway_target_pmids[p].append(pmid)
        for c in cancer_types:
            cancer_pmids[c].append(pmid)
        if evidence:
            evidence_pmids[evidence].append(pmid)
        for r in resistant_states:
            resistant_state_pmids[r].append(pmid)
        for target in radioligand_targets:
            radioligand_target_pmids[target].append(pmid)
        if combination_evidence:
            combination_pmids[combination_evidence].append(pmid)

        # Journal tag
        journal = fm.get("journal", "")
        if journal:
            journal_key = re.sub(r"[^a-z0-9]+", "-", journal.lower()).strip("-")
            if journal_key not in journal_pmids:
                journal_pmids[journal_key] = []
            journal_pmids[journal_key].append(pmid)

        # Stats
        if mechanisms:
            stats["mechanisms"] += 1
        if biology_processes:
            stats["biology_processes"] += 1
        if pathway_targets:
            stats["pathway_targets"] += 1
        if cancer_types:
            stats["cancer_types"] += 1
        if evidence:
            stats["evidence"] += 1
        if resistant_states:
            stats["resistant_states"] += 1
        if radioligand_targets:
            stats["radioligand_targets"] += 1
        if combination_evidence:
            stats["combination_evidence"] += 1

    # Write tag index files
    if not args.dry_run:
        print("\nWriting tag index files...")
        write_tag_files("by-mechanism", mechanism_pmids)
        write_tag_files("by-biology-process", biology_process_pmids)
        write_tag_files("by-pathway-target", pathway_target_pmids)
        write_tag_files("by-cancer-type", cancer_pmids)
        write_tag_files("by-evidence-level", evidence_pmids)
        write_tag_files("by-resistant-state", resistant_state_pmids)
        write_tag_files("by-radioligand-target", radioligand_target_pmids)
        write_tag_files("by-combination-evidence", combination_pmids)
        write_tag_files("by-journal", journal_pmids)

    # Print summary
    print(f"\nTagging complete:")
    print(f"  Articles with mechanism tags: {stats['mechanisms']}/{len(files)}")
    print(f"  Articles with biology-process tags: {stats['biology_processes']}/{len(files)}")
    print(f"  Articles with pathway-target tags: {stats['pathway_targets']}/{len(files)}")
    print(f"  Articles with cancer type tags: {stats['cancer_types']}/{len(files)}")
    print(f"  Articles with evidence level: {stats['evidence']}/{len(files)}")
    print(f"  Articles with resistant-state tags: {stats['resistant_states']}/{len(files)}")
    print(f"  Articles with radioligand target tags: {stats['radioligand_targets']}/{len(files)}")
    print(f"  Articles with combination evidence: {stats['combination_evidence']}/{len(files)}")

    print(f"\nMechanism distribution:")
    for tag, pmids in sorted(mechanism_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nBiology-process distribution:")
    for tag, pmids in sorted(biology_process_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nPathway-target distribution:")
    for tag, pmids in sorted(pathway_target_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nCancer type distribution:")
    for tag, pmids in sorted(cancer_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nEvidence level distribution:")
    for tag, pmids in sorted(evidence_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nResistant-state distribution:")
    for tag, pmids in sorted(resistant_state_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nRadioligand-target distribution:")
    for tag, pmids in sorted(radioligand_target_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nCombination-evidence distribution:")
    for tag, pmids in sorted(combination_pmids.items(), key=lambda x: -len(x[1])):
        if pmids:
            print(f"  {tag}: {len(pmids)}")

    print(f"\nJournals represented: {len([j for j, p in journal_pmids.items() if p])}")

    if not args.dry_run:
        print(f"\nNext step:")
        print(f"  python build_index.py    # Rebuild INDEX.jsonl")


if __name__ == "__main__":
    main()
