#!/usr/bin/env python3
"""Regenerate ONLY the diagnostic-therapy audit (`analysis/diagnostic-therapy-audit.md`)
from the frozen corpus, without running the full `analyze_corpus.py`.

This is the safe regeneration path for the diagnostic-therapy chain expansion (#441).
Chain membership is recomputed on demand from the frozen `corpus/by-pmid/*.md` text
using the current `config.DIAGNOSTIC_THERAPY_KEYWORDS`, so adding chains never
mutates `corpus/INDEX.jsonl` (and therefore cannot shift the frozen 19-mechanism
counts the way a full re-tag would, since a full re-tag would also apply unrelated
post-freeze config drift). The audit is the "living analysis"; the index stays frozen.

A built-in safety check asserts that the recompute reproduces the stored 6-chain
`diagnostic_therapy_links` field exactly for the chains that were present at freeze
time, so any divergence (config drift or a text-extraction change) fails loudly
before the audit is written.

Usage:
    python3 scripts/diagnostic_therapy_audit.py            # regenerate the audit
    python3 scripts/diagnostic_therapy_audit.py --check    # validate only, write nothing
"""
import argparse
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from analyze_corpus import (  # noqa: E402
    ANALYSIS_DIR,
    INDEX_FILE,
    build_diagnostic_therapy_audit,
    load_index,
    recompute_diagnostic_therapy_links,
)

# The chains that were present in the frozen corpus's stored field. The recompute,
# restricted to these, must reproduce the stored field exactly.
FROZEN_CHAINS = {
    "psma-imaging-to-radioligand",
    "sstr-imaging-to-prrt",
    "pdl1-ihc-to-checkpoint",
    "tmb-msi-to-immunotherapy",
    "neoantigen-profiling-to-mrna-vaccine",
    "oncolytic-susceptibility-to-virotherapy",
}


def stored_links() -> dict[str, list[str]]:
    out = {}
    with open(INDEX_FILE, encoding="utf-8") as f:
        for line in f:
            e = json.loads(line)
            out[str(e["pmid"])] = sorted(e.get("diagnostic_therapy_links") or [])
    return out


def validate_frozen_reproduction(links_by_pmid: dict[str, list[str]]) -> None:
    """Assert the recompute restricted to the frozen chains reproduces the stored
    field exactly. Fails loudly (SystemExit) on any mismatch."""
    stored = stored_links()
    mismatches = []
    for pmid, stored_chains in stored.items():
        recomputed_frozen = sorted(c for c in links_by_pmid.get(pmid, []) if c in FROZEN_CHAINS)
        if recomputed_frozen != stored_chains:
            mismatches.append((pmid, stored_chains, recomputed_frozen))
    if mismatches:
        print(
            f"SAFETY CHECK FAILED: the recompute does not reproduce the frozen "
            f"diagnostic_therapy_links for {len(mismatches)} record(s). The frozen "
            f"corpus or the diagnostic-therapy config may have drifted; refusing to "
            f"regenerate the audit.",
            file=sys.stderr,
        )
        for pmid, st, rc in mismatches[:10]:
            print(f"  PMID {pmid}: stored={st} recomputed(frozen)={rc}", file=sys.stderr)
        raise SystemExit(1)
    print(f"Safety check passed: recompute reproduces the stored field for all "
          f"{len(stored)} records on the {len(FROZEN_CHAINS)} frozen chains.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="Validate the frozen reproduction only; write nothing."
    )
    args = parser.parse_args()

    links_by_pmid = recompute_diagnostic_therapy_links()
    validate_frozen_reproduction(links_by_pmid)
    if args.check:
        return

    entries = load_index()
    content = build_diagnostic_therapy_audit(entries, links_by_pmid=links_by_pmid)
    out_path = ANALYSIS_DIR / "diagnostic-therapy-audit.md"
    out_path.write_text(content, encoding="utf-8")
    total = sum(1 for e in entries if links_by_pmid.get(str(e.get("pmid"))))
    print(f"Wrote {out_path} ({total} articles with at least one link).")


if __name__ == "__main__":
    main()
