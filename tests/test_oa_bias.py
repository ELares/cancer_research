"""Guard the #348 open-access-bias analysis against corpus drift.

These pin the load-bearing claims the manuscript §3.3.1 retrieval-bias
subsection now cites: the full-text corpus is overwhelmingly OA, the
abstract-only archive is mostly non-OA, immunotherapy is #1 in both rankings,
and the physical/device modalities are OA-suppressed (they rank far higher
once the non-OA abstract literature is included).
"""

import collections
import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "oa_bias_analysis", REPO / "scripts" / "oa_bias_analysis.py"
)
oa = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(oa)

# Load each corpus once (each reads thousands of files); reuse across tests.
_FT = oa.load_fulltext()
_AB = oa.load_abstracts()


def _counts(records):
    c = collections.Counter()
    for r in records:
        for m in r["mechanisms"]:
            c[m] += 1
    return c


def test_fulltext_is_overwhelmingly_oa():
    ft = _FT
    assert len(ft) > 4000, "full-text corpus should be the ~4,830-record set"
    oa_rate = sum(1 for r in ft if r["is_oa"]) / len(ft)
    assert oa_rate > 0.95, f"full-text corpus is OA by construction; got {oa_rate:.3f}"


def test_abstract_archive_is_mostly_non_oa():
    ab = _AB
    assert len(ab) > 5000, "abstract-only archive should be the ~5,585-record set"
    oa_rate = sum(1 for r in ab if r["is_oa"]) / len(ab)
    # The whole point: the abstract archive carries the non-OA literature.
    assert oa_rate < 0.5, f"abstract archive should be mostly non-OA; got {oa_rate:.3f}"


def test_immunotherapy_is_number_one_in_both_rankings():
    ft = _counts(_FT)
    ab = _counts(_AB)
    assert ft.most_common(1)[0][0] == "immunotherapy"
    assert ab.most_common(1)[0][0] == "immunotherapy"


def test_no_pmid_overlap():
    """The combined-count analysis (Total = full-text + abstract) is valid only
    if the two corpora are disjoint. Guard that invariant against corpus drift
    (e.g. an abstract-only PMID promoted to full text without removing the
    abstract would silently double-count it)."""
    ft_pmids = {r["pmid"] for r in _FT if r["pmid"]}
    ab_pmids = {r["pmid"] for r in _AB if r["pmid"]}
    overlap = ft_pmids & ab_pmids
    assert not overlap, (
        f"{len(overlap)} PMIDs are in BOTH corpora, breaking the combined-count "
        f"disjointness assumption: {sorted(overlap)[:5]}"
    )


def test_physical_modalities_are_oa_suppressed():
    """A representative physical modality ranks far higher among abstracts than
    in the OA-biased full-text corpus."""
    ft_rank = oa.ranks(_counts(_FT))
    ab_rank = oa.ranks(_counts(_AB))
    # bioelectric is the clearest case: it should rank much better among abstracts.
    assert ft_rank["bioelectric"] - ab_rank["bioelectric"] >= 5, (
        f"bioelectric should rank >=5 places higher among abstracts: "
        f"full-text {ft_rank['bioelectric']} vs abstract {ab_rank['bioelectric']}"
    )
