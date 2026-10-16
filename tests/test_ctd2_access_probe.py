"""The chemotherapy row's route claim must match the probe that established it (#958).

`CALIBRATION_STATUS.md`'s #833 row previously ended "A figshare deposit for the
CTD^2 release itself was not located; that is where a reader should look" -- the
one route #932 named and did not probe. It is probed now, and the answer is a
different host, so the row carries a DOI, a size and an md5.

Those are QUOTED EXTERNAL FACTS, so they live in data
(`analysis/calibration/ctd2-access-probe.json`) rather than only in prose -- the
`published-interval-claims.csv` precedent -- and these guards exist to stop the
row drifting from the artifact that establishes it.

Nothing here touches the network: the probe was taken locally and committed, so
CI reads the JSON and downloads nothing.
"""

import ast
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PROBE = REPO / "analysis" / "calibration" / "ctd2-access-probe.json"
STATUS = REPO / "simulations" / "calibration" / "CALIBRATION_STATUS.md"


@pytest.fixture(scope="module")
def probe():
    return json.loads(PROBE.read_text())


@pytest.fixture(scope="module")
def row():
    """The one CALIBRATION_STATUS row this probe qualifies."""
    for line in STATUS.read_text().splitlines():
        if line.startswith("|") and "#833" in line:
            return line
    pytest.fail("the chemotherapy (#833) row is gone from CALIBRATION_STATUS.md")


def test_exactly_one_route_is_open(probe):
    """A probe listing no open route, or several, is not a finding."""
    open_routes = [r for r in probe["routes"] if r["status"] == "open"]
    assert len(open_routes) == 1, (
        f"{len(open_routes)} routes marked open; the finding is that exactly "
        "one is")
    assert open_routes[0]["locator"] == probe["open_route"]["doi"], (
        "the open route and the detailed open_route block must be the same one")


def test_every_route_carries_its_evidence(probe):
    """A status with no evidence is an assertion, which is what this replaces."""
    for r in probe["routes"]:
        assert r.get("evidence"), f"{r['route']} has a status and no evidence"
        assert r.get("established_by"), f"{r['route']} does not say who found it"


def test_the_row_quotes_the_artifact_and_cannot_drift(probe, row):
    """Every external fact in the prose must appear in the data."""
    assert probe["open_route"]["doi"] in row
    assert probe["open_route"]["md5"] in row
    assert f"{probe['open_route']['bytes']:,}" in row, (
        "the row must quote the measured size as the artifact records it")


def test_the_row_still_refuses_what_the_probe_does_not_establish(probe, row):
    """Locating a file fits nothing, and the row has to keep saying so."""
    assert probe["verdict_unchanged"].startswith("NO FITTED POTENCY")
    assert "NO fitted potency" in row or "NO FITTED POTENCY" in row
    # The blocker is downgraded, not removed.
    assert probe["remaining_blocker"]["kind"] == "toolchain"
    assert probe["remaining_blocker"]["not"] == "access"
    assert "TOOLCHAIN" in row.upper()


def test_the_unverified_half_is_declared_in_both_places(probe, row):
    """Reachability was measured; the contents were not, and that must show."""
    assert "file contents" in " ".join(probe["open_route"]["not_measured"])
    assert probe["open_route"]["carries_per_dose_basis"].startswith(
        "documentation, not measurement"), (
        "the per-dose claim rests on the format's documentation, so the "
        "artifact must say so rather than implying it was measured")
    assert re.search(r"not been downloaded|NOT been downloaded", row), (
        "the row must state that the file was not downloaded or parsed")


def test_the_ftp_absence_is_a_checked_absence_not_a_failed_guess(probe, row):
    """The distinction is the whole reason that line is believable."""
    ftp = next(r for r in probe["routes"] if "FTP" in r["route"])
    assert ftp["status"] == "removed"
    assert "absent from its own parent" in ftp["evidence"]
    assert "checked absence" in row, (
        "two guessed filenames failed first; what establishes the tree is gone "
        "is listing the parent, and the row must not blur the two")


def test_the_probe_does_not_claim_a_dependency_decision(probe):
    """Adding R is the owner's call, not a side effect of probing a URL."""
    assert "dependency contract" in probe["remaining_blocker"]["decision_deferred"]
    assert "requirements-lock.txt" in " ".join(probe["remaining_blocker"]["detail"])


def test_nothing_here_reaches_the_network():
    """The probe is committed; CI must not re-run it.

    AST over the IMPORTS, not a substring scan. The first version scanned its
    own source for a list of banned module names and failed on the list
    itself -- a text check over a file containing the pattern matches itself,
    the same shape as a `pgrep` watcher finding its own command line.
    """
    tree = ast.parse(Path(__file__).read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    network = imported & {"urllib", "requests", "socket", "http", "httpx",
                          "ftplib", "subprocess"}
    assert not network, (
        f"this guard must read the committed probe, not re-fetch: {sorted(network)}")
