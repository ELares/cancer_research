"""
Invariant and property tests for corpus analysis outputs.

These tests verify mathematical properties that must hold regardless
of corpus size or content: monotonicity, positivity, ordering,
schema consistency. They catch broken formulas and configuration
drift without relying on exact floating-point matching.

Run: pytest tests/test_invariants.py -v
"""

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

INDEX_FILE = REPO_ROOT / "corpus" / "INDEX.jsonl"


# ============================================================
# INDEX.jsonl schema invariants
# ============================================================

class TestIndexSchema:
    """INDEX.jsonl should have consistent required fields in every entry."""

    REQUIRED_FIELDS = {
        "pmid", "title", "mechanisms", "cancer_types", "evidence_level",
        "journal", "year", "doi", "is_oa", "oa_status", "author_count",
        "biology_processes", "pathway_targets", "resistant_states",
        "tissue_categories", "date_added", "pmcid", "cited_by_count", "month",
    }

    @pytest.fixture
    def entries(self):
        if not INDEX_FILE.exists():
            pytest.skip("INDEX.jsonl not found")
        return [json.loads(line) for line in INDEX_FILE.read_text().splitlines() if line.strip()]

    def test_all_entries_have_required_fields(self, entries):
        """Every entry must contain all 19 required fields."""
        for i, entry in enumerate(entries):
            missing = self.REQUIRED_FIELDS - set(entry.keys())
            assert not missing, (
                f"Entry {i} (PMID {entry.get('pmid', '?')}) missing fields: {missing}"
            )

    def test_all_entries_have_nonempty_pmid(self, entries):
        """Every entry must have a non-empty pmid."""
        for entry in entries:
            assert entry.get("pmid"), (
                f"Entry missing or empty pmid: {entry.get('title', '?')[:60]}"
            )

    def test_mechanisms_is_list(self, entries):
        """Mechanisms field must be a list (possibly empty)."""
        for entry in entries:
            assert isinstance(entry.get("mechanisms"), list), (
                f"PMID {entry['pmid']}: mechanisms is {type(entry.get('mechanisms'))}, expected list"
            )

    def test_evidence_levels_valid(self, entries):
        """Evidence level must be from the valid set or empty string."""
        valid = {
            "phase3-clinical", "phase2-clinical", "phase1-clinical",
            "clinical-other", "preclinical-invivo", "preclinical-invitro",
            "theoretical", "",
        }
        for entry in entries:
            level = entry.get("evidence_level", "")
            assert level in valid, (
                f"PMID {entry['pmid']}: invalid evidence_level '{level}'"
            )


# ============================================================
# Evidence weight formula invariants
# ============================================================

class TestWeightInvariants:
    """Evidence weighting formula must satisfy mathematical properties."""

    def test_tier_monotonicity(self):
        """Higher evidence tier must always produce higher weight (same citation/year)."""
        from analyze_corpus import evidence_weight

        tiers = [
            "theoretical", "preclinical-invitro", "preclinical-invivo",
            "clinical-other", "phase1-clinical", "phase2-clinical", "phase3-clinical",
        ]
        base_entry = {"icite_percentile": 50, "year": 2023}
        prev = 0.0
        for tier in tiers:
            w = evidence_weight({**base_entry, "evidence_level": tier})
            assert w > prev, f"{tier} weight {w:.2f} should exceed previous {prev:.2f}"
            prev = w

    def test_all_tiers_produce_positive_weight(self):
        """Every valid evidence tier must produce a positive weight."""
        from analyze_corpus import evidence_weight, EVIDENCE_TIER_WEIGHTS

        for tier in EVIDENCE_TIER_WEIGHTS:
            entry = {"evidence_level": tier, "icite_percentile": 0, "year": 2015}
            w = evidence_weight(entry)
            assert w > 0, f"Weight for {tier} should be positive, got {w}"

    def test_weight_upper_bound(self):
        """Max weight must not exceed theoretical maximum (~19.8)."""
        from analyze_corpus import evidence_weight

        # Max: phase3(12) × citation(1.5) × recency(1.1) = 19.8
        entry = {"evidence_level": "phase3-clinical", "icite_percentile": 100, "year": 2026}
        w = evidence_weight(entry)
        assert w <= 20.0, f"Max weight {w:.2f} exceeds theoretical max ~19.8"
        assert w > 19.0, f"Max weight {w:.2f} unexpectedly low (expected ~19.8)"

    def test_no_evidence_returns_zero(self):
        """Missing or invalid evidence level must return zero weight."""
        from analyze_corpus import evidence_weight

        assert evidence_weight({}) == 0.0
        assert evidence_weight({"evidence_level": ""}) == 0.0
        assert evidence_weight({"evidence_level": "garbage"}) == 0.0


# ============================================================
# Corpus-level invariants
# ============================================================

class TestCorpusInvariants:
    """Properties that must hold for any valid corpus."""

    @pytest.fixture
    def entries(self):
        if not INDEX_FILE.exists():
            pytest.skip("INDEX.jsonl not found")
        return [json.loads(line) for line in INDEX_FILE.read_text().splitlines() if line.strip()]

    def test_mechanism_counts_are_non_negative(self, entries):
        """No mechanism should have a negative article count (sanity check on Counter logic)."""
        from collections import Counter

        counts = Counter()
        for e in entries:
            for m in e.get("mechanisms", []):
                counts[m] += 1
        for mech, count in counts.items():
            assert count > 0, f"Mechanism {mech} has non-positive count {count}"

    def test_weighted_scores_consistent_with_tagged_counts(self, entries):
        """Mechanisms with more tagged evidence articles cannot have zero total weight."""
        from analyze_corpus import evidence_weight

        from collections import defaultdict
        mech_tagged = defaultdict(list)
        for e in entries:
            if e.get("evidence_level"):
                for m in e.get("mechanisms", []):
                    mech_tagged[m].append(e)

        for mech, tagged in mech_tagged.items():
            total_weight = sum(evidence_weight(e) for e in tagged)
            assert total_weight > 0, (
                f"Mechanism {mech} has {len(tagged)} tagged articles but zero total weight"
            )
