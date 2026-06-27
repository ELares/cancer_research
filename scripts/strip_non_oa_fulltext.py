"""Audit + strip non-redistributable full text from the corpus (#526).

Some `corpus/by-pmid/*.md` records are flagged `is_oa: false` by OpenAlex yet were
stored WITH full text. Having a PMCID is not the same as being in the PMC Open
Access subset (the only PMC tier that permits redistribution); only the OA-subset
service (oa.fcgi) is authoritative. Redistributing copyrighted full text under an
MIT / "free forever" repo is a real legal + credibility risk that PROVENANCE.yaml
already committed to fixing.

Two modes:

  python scripts/strip_non_oa_fulltext.py --audit
      Query the NCBI OA service for every is_oa==false record with a PMCID,
      classify in-subset (keep) vs not-in-subset (strip), and write
      analysis/non-oa-fulltext-audit.{csv,md}. NEEDS NETWORK.

  python scripts/strip_non_oa_fulltext.py --apply
      Read the committed audit CSV and, for each not-in-subset record, remove the
      "## Full Text" section from corpus/by-pmid/<pmid>.md IN PLACE (keeping the
      YAML frontmatter, title, and abstract), replacing it with a removal notice.
      OFFLINE + idempotent.

Strip-in-place (not relocate) is deliberate: the record stays in the full-text
index so the frozen corpus count is unchanged, and the abstract + metadata that
the tagging / chain recompute read are preserved (verified: the diagnostic-therapy
chain recompute still matches the stored field, because chain keywords live in the
abstract, not only the body). Only the copyrighted body is removed.
"""

import argparse
import csv
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INDEX = REPO_ROOT / "corpus" / "INDEX.jsonl"
PMID_DIR = REPO_ROOT / "corpus" / "by-pmid"
AUDIT_CSV = REPO_ROOT / "analysis" / "non-oa-fulltext-audit.csv"
AUDIT_MD = REPO_ROOT / "analysis" / "non-oa-fulltext-audit.md"
OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id="

# Non-greedy + stop at the next "## " heading (or EOF) so only the Full Text
# section is removed, never a trailing section (e.g. a `## Source` footer) — the
# NOTICE promises the rest of the record is retained.
FULLTEXT_HEADING = re.compile(
    r"\n##\s+Full[ -]?Text\b.*?(?=\n##\s|\Z)", re.IGNORECASE | re.DOTALL
)
NOTICE = (
    "\n## Full Text\n\n"
    "> **Full text removed (#526, redistribution audit).** This article is not in "
    "the PMC Open Access subset; its full text is under publisher copyright and is "
    "not redistributable. The abstract and metadata above are retained for the "
    "corpus index; see the DOI for the full article.\n"
)


def non_oa_records():
    recs = []
    for line in INDEX.read_text(encoding="utf-8").splitlines():
        d = json.loads(line)
        if d.get("is_oa") is False and d.get("pmcid"):
            recs.append(d)
    return recs


def in_oa_subset(pmcid: str) -> bool:
    """True iff the PMCID is in the redistributable PMC OA subset (oa.fcgi returns a
    licensed package link); False on any error / not-open-access response."""
    with urllib.request.urlopen(OA_URL + pmcid, timeout=20) as f:
        x = f.read().decode("utf-8", "ignore")
    if "<error" in x or "idIsNotOpenAccess" in x or 'returned-count="0"' in x:
        return False
    return ("href=" in x) and ("license" in x)


def audit():
    recs = non_oa_records()
    rows = []
    for i, d in enumerate(recs):
        pmcid = d["pmcid"]
        try:
            keep = in_oa_subset(pmcid)
            err = ""
        except Exception as e:  # noqa: BLE001 — record + continue
            keep, err = True, str(e)[:60]  # on error, default to KEEP (conservative)
        rows.append({
            "pmid": str(d["pmid"]), "pmcid": pmcid,
            "oa_status": d.get("oa_status", ""), "journal": d.get("journal", ""),
            "in_pmc_oa_subset": keep, "error": err,
        })
        time.sleep(0.34)  # NCBI courtesy (no key)
    rows.sort(key=lambda r: (r["in_pmc_oa_subset"], r["pmid"]))
    AUDIT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    strip = [r for r in rows if not r["in_pmc_oa_subset"]]
    keep = [r for r in rows if r["in_pmc_oa_subset"]]
    lines = [
        "# Non-OA full-text redistribution audit (#526)",
        "",
        "Each `is_oa: false` corpus record with a PMCID, checked against the NCBI PMC "
        "Open Access subset service (`oa.fcgi`). Records NOT in the OA subset have "
        "their `## Full Text` stripped in place (abstract + metadata kept; the "
        "frozen corpus count is unchanged) via "
        "`scripts/strip_non_oa_fulltext.py --apply`.",
        "",
        f"- Non-OA records with a PMCID: **{len(rows)}**",
        f"- In the PMC OA subset (redistributable, KEPT): **{len(keep)}**",
        f"- NOT in the OA subset (copyrighted, STRIPPED): **{len(strip)}**",
        "",
        "## Stripped (not in the OA subset)",
        "",
        "| PMID | PMCID | oa_status | Journal |",
        "|---|---|---|---|",
    ]
    for r in strip:
        lines.append(f"| {r['pmid']} | {r['pmcid']} | {r['oa_status']} | {r['journal']} |")
    AUDIT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Audited {len(rows)} records: {len(keep)} kept (OA subset), "
          f"{len(strip)} to strip. Wrote {AUDIT_CSV.name} + {AUDIT_MD.name}.")


def apply():
    if not AUDIT_CSV.exists():
        sys.exit("No audit CSV — run `--audit` first (needs network).")
    to_strip = [r["pmid"] for r in csv.DictReader(AUDIT_CSV.open(encoding="utf-8"))
                if r["in_pmc_oa_subset"].strip().lower() in ("false", "0", "no")]
    stripped = skipped = 0
    for pmid in to_strip:
        fp = PMID_DIR / f"{pmid}.md"
        if not fp.exists():
            continue
        text = fp.read_text(encoding="utf-8")
        if "Full text removed (#526" in text:
            skipped += 1
            continue
        if not FULLTEXT_HEADING.search(text):
            skipped += 1  # no body to strip (already abstract-only)
            continue
        new = FULLTEXT_HEADING.sub(NOTICE.rstrip("\n"), text).rstrip() + "\n"
        fp.write_text(new, encoding="utf-8")
        stripped += 1
    print(f"Stripped {stripped} record bodies; skipped {skipped} "
          f"(already stripped / no body). Total flagged: {len(to_strip)}.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--audit", action="store_true", help="query oa.fcgi, write the audit (network)")
    ap.add_argument("--apply", action="store_true", help="strip flagged bodies in place (offline)")
    args = ap.parse_args()
    if args.audit:
        audit()
    if args.apply:
        apply()
    if not (args.audit or args.apply):
        ap.print_help()


if __name__ == "__main__":
    main()
