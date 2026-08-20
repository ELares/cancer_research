"""The manuscript's Section 4.5 landmark claim, checked against the census.

Section 4.5 used to list five clinically important papers CONFIRMED MISSING from
the local archive, and warned that each was large enough to distort
mechanism-level claims. It now says the opposite: all five are present, and
VISION and PANOVA-3 are indexed with Phase III trial publication types.

A retraction that flips an absence into a presence is exactly the kind a reader
has to take on trust unless something checks it, and this repo has been wrong in
BOTH directions on this file before -- the corpus once reported a landmark
missing that a later recovery added, and an earlier survey concluded full-text
data was absent when a machine was merely switched off. So the claim is pinned
to the census rather than to the sentence.

OFFLINE CONTRACT: the census is gitignored, so these skip when it is absent
rather than failing. A skip is not a pass and is not silent -- the shape of the
guard is that when the data IS present, the manuscript's sentence must be true
of it.
"""
import gzip
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
MANUSCRIPT = REPO / "article/drafts/v1.md"

# The five the manuscript names, and what it says about each.
LANDMARKS = {
    "34161051": "VISION",
    "40448572": "PANOVA-3",
    "33016924": "mRNA vaccine study",
    "36027916": "mRNA vaccine study",
    "35970920": "mRNA vaccine study",
}
# The two the manuscript states are indexed as Phase III.
PHASE_III_CLAIMED = {"34161051", "40448572"}
PHASE_III = "Clinical Trial, Phase III"


def _census_available() -> bool:
    return RECORDS.is_dir() and any(RECORDS.glob("*.jsonl.gz"))


def _find(pmids: set[str]) -> dict[str, dict]:
    """Scan for a small set of PMIDs. The cheap `in line` prefilter is a
    SUBSTRING test and would match a PMID appearing inside another field, so
    every hit is confirmed by parsing and comparing the pmid field."""
    out: dict[str, dict] = {}
    for f in sorted(RECORDS.glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                if not any(p in line for p in pmids):
                    continue
                r = json.loads(line)
                p = str(r.get("pmid"))
                if p in pmids and p not in out:
                    out[p] = r
        if len(out) == len(pmids):
            break
    return out


@pytest.mark.skipif(not _census_available(), reason="census not present (gitignored)")
def test_every_landmark_the_manuscript_calls_present_is_in_the_census():
    found = _find(set(LANDMARKS))
    missing = sorted(set(LANDMARKS) - set(found))
    assert not missing, (
        "Section 4.5 states all five landmark papers are present in the census, "
        f"and these are not: {[f'{p} ({LANDMARKS[p]})' for p in missing]}. Either "
        "the census build changed or the sentence is wrong; do not weaken the "
        "guard to match the data without establishing which.")


@pytest.mark.skipif(not _census_available(), reason="census not present (gitignored)")
def test_the_two_trials_carry_the_publication_type_the_manuscript_quotes():
    """The manuscript quotes the LABEL, not just presence. A record present but
    indexed without a trial type would leave the section's argument -- that the
    earlier analysis could not see a trial that was there and labelled -- false
    while the presence test still passed."""
    found = _find(PHASE_III_CLAIMED)
    for pmid in sorted(PHASE_III_CLAIMED):
        rec = found.get(pmid)
        assert rec is not None, f"{pmid} absent from the census"
        types = rec.get("pub_types") or []
        assert PHASE_III in types, (
            f"the manuscript says {LANDMARKS[pmid]} ({pmid}) is indexed "
            f"{PHASE_III!r}; NLM records {types}")


def test_the_manuscript_still_makes_the_claim_this_file_checks():
    """A guard for a sentence that has been edited away is not a guard.

    Without this, rewriting Section 4.5 back to "confirmed missing" would leave
    two green tests asserting the opposite, and nothing would say so.
    """
    txt = " ".join(MANUSCRIPT.read_text(encoding="utf-8").split())
    assert "All five are present in the census" in txt, (
        "Section 4.5 no longer claims the five landmarks are present, so this "
        "file is checking a claim the manuscript does not make. Update or "
        "delete it deliberately rather than leaving it green.")
    for pmid in sorted(PHASE_III_CLAIMED):
        assert pmid in txt, (
            f"the manuscript no longer names {pmid} ({LANDMARKS[pmid]}), which "
            "this file pins")
