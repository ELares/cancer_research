"""Bind new direction labels to the exact records selected for a census run.

Historical title-only CSVs cannot supply this evidence. New exports carry the
complete candidate text, a hash of each original JSON record, and a cohort hash
covering every record in the analysis intersection, including unclassified
records. These hashes detect mismatched inputs; they do not authenticate labels
or establish adjudicator accuracy.
"""
import csv
import hashlib
import io
import json
from pathlib import Path


CSV_FIELDS = (
    "cohort_sha256", "record_sha256", "pmid", "title", "abstract",
    "regex_label", "adjudicated", "reason",
)
IDENTITY_FIELDS = ("record_sha256", "pmid", "title", "abstract", "regex_label")


def _digest(value) -> str:
    encoded = json.dumps(value, sort_keys=True, ensure_ascii=False,
                         separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class AdjudicationCohort:
    """Collect a small analysis intersection, not the whole census in memory."""

    def __init__(self, analysis: str, labels: tuple[str, ...]):
        self.analysis = analysis
        self.labels = labels
        self._members = {}
        self._pmids = set()
        self._candidates = []

    def add(self, record: dict, regex_label: str) -> None:
        if regex_label not in (*self.labels, "both", "neither"):
            raise ValueError(f"Unknown candidate label: {regex_label}")
        fingerprint = _digest(record)
        pmid = str(record.get("pmid") or "").strip()
        if pmid.isdecimal():
            pmid = str(int(pmid))
        if fingerprint in self._members or (pmid and pmid in self._pmids):
            raise SystemExit(
                "Duplicate record in the selected adjudication cohort; "
                "existing reports were not changed.")
        self._members[fingerprint] = regex_label
        if pmid:
            self._pmids.add(pmid)
        if regex_label in self.labels:
            title, abstract = record.get("title") or "", record.get("abstract") or ""
            if not isinstance(title, str) or not isinstance(abstract, str):
                raise SystemExit("Candidate title and abstract must be text.")
            self._candidates.append({
                "record_sha256": fingerprint, "pmid": pmid,
                "title": title, "abstract": abstract, "regex_label": regex_label,
            })

    def as_dict(self) -> dict:
        return {
            "schema_version": 1,
            "analysis": self.analysis,
            "cohort_sha256": _digest({
                "schema_version": 1, "analysis": self.analysis,
                "records": sorted(self._members.items()),
            }),
            "records": len(self._members),
            "candidates": [dict(row) for row in self._candidates],
        }


def candidate_csv(cohort: dict) -> str:
    """Create an explicit review worksheet; no decisions are preselected."""
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    for candidate in cohort["candidates"]:
        writer.writerow({**candidate, "cohort_sha256": cohort["cohort_sha256"],
                         "adjudicated": "", "reason": ""})
    return output.getvalue()


def read_adjudication(path: Path, cohort: dict, allowed_decisions) -> list[dict]:
    """Require exactly one completed decision for every current candidate.

    Full coverage removes any need to extrapolate an old sample's precision.
    Ambiguous (and, where supported, off-topic) decisions remain explicit and
    are excluded by the analysis from its directional denominator.
    """
    path = Path(path)

    def fail(message):
        raise SystemExit(
            f"Adjudication {path}: {message}; existing reports were not changed.")

    expected = {row["record_sha256"]: row for row in cohort["candidates"]}
    if not expected:
        fail("the selected cohort has no singly classified candidates to adjudicate")
    try:
        with path.open(encoding="utf-8-sig", newline="") as source:
            reader = csv.DictReader(source, strict=True)
            if reader.fieldnames is None or (
                len(reader.fieldnames) != len(CSV_FIELDS)
                or set(reader.fieldnames) != set(CSV_FIELDS)
            ):
                fail("expected the complete candidate-export columns; "
                     "historical title-only CSVs cannot be applied to a new cohort")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        fail(f"cannot read a valid CSV ({exc})")

    seen = set()
    for number, row in enumerate(rows, 2):
        if set(row) != set(CSV_FIELDS) or any(v is None for v in row.values()):
            fail(f"row {number} has missing or extra fields")
        identity = row["record_sha256"]
        if identity in seen:
            fail(f"row {number} duplicates a candidate")
        seen.add(identity)
        if row["cohort_sha256"] != cohort["cohort_sha256"]:
            fail(f"row {number} belongs to a different selected cohort")
        if identity not in expected:
            fail(f"row {number} is not a selected candidate")
        candidate = expected[identity]
        if any(row[key] != candidate[key] for key in IDENTITY_FIELDS):
            fail(f"row {number} changes the exported candidate text or label")
        if row["adjudicated"] not in allowed_decisions:
            fail(f"row {number} has an invalid or unfinished adjudicated label")
        if not row["reason"].strip():
            fail(f"row {number} needs an adjudication reason")
    missing = set(expected) - seen
    if missing:
        fail(f"missing decisions for {len(missing)} selected candidate(s)")
    return rows
