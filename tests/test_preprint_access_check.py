"""Guards for the preprint-access measurement.

The finding is a NEGATIVE: bioRxiv's per-article full text is blocked to
automated clients, so ~92,000 identified cancer preprints stay unread. A
negative result rots differently from a positive one -- it stops being true the
moment the block is lifted, and nothing fails when that happens. These pin the
one thing that can be checked offline: the page says what the probe found.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import preprint_access_check as pac  # noqa: E402


def _d(reachable: bool):
    return {"preprints_indexed": 108400, "fulltext_in_epmc": 16432,
            "share_pct": 15.2, "unread": 91968,
            "biorxiv": ({"metadata_api": "ok", "fulltext": "HTTP 200",
                         "bytes": 123456, "reachable": True} if reachable else
                        {"metadata_api": "ok", "fulltext": "HTTP 429",
                         "reachable": False})}


def test_the_page_reports_a_block_as_a_block():
    md = pac.render(_d(False))
    assert "HTTP 429" in md and "standing block" in md
    assert "That block is respected" in md


def test_the_page_stops_claiming_a_block_once_it_is_lifted():
    """THE WAY A NEGATIVE RESULT ROTS. If bioRxiv drops the block, ~92,000
    preprints become collectable and nothing would otherwise say so -- the page
    would keep asserting an obstacle that no longer exists."""
    md = pac.render(_d(True))
    assert "currently reachable" in md
    assert "standing block" not in md
    assert "needs revising" in md


def test_every_count_on_the_page_comes_from_the_measurement():
    """A hand-written figure beside a derived one is this repository's most
    persistent defect. Change the inputs and the page must change with them."""
    d = _d(False)
    d.update(preprints_indexed=200, fulltext_in_epmc=50, share_pct=25.0, unread=150)
    md = pac.render(d)
    assert "108,400" not in md and "91,968" not in md
    assert "200" in md and "150" in md and "25.0%" in md


def test_the_probe_is_not_disguised():
    """A browser User-Agent would very likely get a different answer, and
    sending one would circumvent an access control the operator put there on
    purpose. The measurement is only honest if the client identifies itself."""
    src = Path(pac.__file__).read_text()
    assert "cancer-research-corpus" in src
    for word in ("Mozilla", "Chrome", "Safari", "AppleWebKit", "Gecko"):
        assert word not in src, f"the probe impersonates a browser ({word})"


def test_the_probe_makes_one_request_not_a_crawl():
    """The point is to observe the answer, not to argue with it."""
    src = Path(pac.__file__).read_text()
    fn = src[src.index("def biorxiv_fulltext_status"):src.index("def render(")]
    assert "for " not in fn and "while " not in fn, (
        "the probe loops; a blocked endpoint must be asked once")
