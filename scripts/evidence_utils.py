"""Shared helpers for evidence-tier classification and audit reporting."""

import re

REVIEW_MARKERS = (
    "review", "systematic review", "meta-analysis", "meta analysis",
    "scoping review", "narrative review", "evidence map",
)
PROTOCOL_MARKERS = ("protocol", "study protocol", "trial protocol", "protocol for")


def normalize_text(text: str) -> str:
    """Normalize case and whitespace so keyword matching survives Unicode spacing."""
    return re.sub(r"\s+", " ", text).strip().lower()


def is_review_like(fm: dict) -> bool:
    """Return True for reviews, meta-analyses, evidence maps, and similar summaries."""
    pub_types = [normalize_text(p) for p in fm.get("pub_types", [])]
    title = normalize_text(fm.get("title", ""))
    return any("review" in p or "meta-analysis" in p for p in pub_types) or any(
        marker in title for marker in REVIEW_MARKERS
    )


def is_protocol_like(fm: dict) -> bool:
    """Return True for protocols and planned studies that should not count as completed evidence."""
    pub_types = [normalize_text(p) for p in fm.get("pub_types", [])]
    title = normalize_text(fm.get("title", ""))
    return any("protocol" in p for p in pub_types) or any(marker in title for marker in PROTOCOL_MARKERS)
