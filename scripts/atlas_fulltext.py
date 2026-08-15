#!/usr/bin/env python3
"""Atlas: pull open-access full text for the cancer census (#ATLAS).

WHY BULK PACKAGES AND NOT PER-ARTICLE FETCHES
---------------------------------------------
PMC publishes its Open Access Subset as plain-text tarballs partitioned by
PMCID range: 61 commercial-use packages, 61 non-commercial, 36 other. Pulling
the cancer subset article-by-article from the AWS bucket would be ~900,000
HTTPS round trips; streaming 122 tarballs is a few hundred. Measured: 34 MB /
3,028 articles for the first package at ~12 MB/s.

Each package is streamed once, members whose PMCID is in the cancer census are
kept, and the tarball is deleted. Nothing is held in memory beyond one article.

LICENCE
-------
`oa_comm` is the commercial-use subset and `oa_noncomm` is non-commercial;
both are readable for analysis, they differ in REDISTRIBUTION terms. The
licence class is recorded on every record so any downstream release can filter.
`oa_other` is excluded by default -- its terms vary per article.

WHERE IT WRITES
---------------
`FERRO_ATLAS_FULLTEXT`, defaulting to the NAS mount, because the cancer slice
of PMC plain text is tens of GB. Output is sharded gzipped JSONL:

    {"pmcid": "PMC176545", "pmid": "12734009", "licence": "comm", "text": "..."}

Requires the census from `scripts/atlas_baseline.py` (for the PMCID -> PMID map).
Resumable via a manifest. Stdlib only.

Usage:
    python scripts/atlas_fulltext.py --limit 1        # one package, smoke test
    python scripts/atlas_fulltext.py                  # full pull
    python scripts/atlas_fulltext.py --status
"""

import argparse
import gzip
import json
import os
import re
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402

BULK = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/deprecated/oa_bulk"
SUBSETS = {"comm": "oa_comm", "noncomm": "oa_noncomm"}
USER_AGENT = "cancer_research-atlas/1.0 (https://github.com/ELares/cancer_research)"
# Full text is bulky enough to live off the repo disk, typically on external
# storage. Derived from the home directory rather than written out, so no
# machine-specific path or account name is committed to a public repository.
# Override with FERRO_ATLAS_FULLTEXT to put it anywhere else.
DEFAULT_FT = Path.home() / "nas" / "cancer-atlas" / "fulltext"

_PKG_RE = re.compile(r'(oa_(?:comm|noncomm)_txt\.PMC\d+xxxxxx\.baseline\.[\d-]+\.tar\.gz)')
_PMCID_RE = re.compile(r"(PMC\d+)")


def fulltext_root() -> Path:
    return Path(os.getenv("FERRO_ATLAS_FULLTEXT", str(DEFAULT_FT)))


def _get(url: str, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def list_packages() -> list[tuple[str, str]]:
    """[(licence, package-name)] across the redistributable subsets."""
    out = []
    for lic, sub in SUBSETS.items():
        html = _get(f"{BULK}/{sub}/txt/", timeout=180).decode("utf-8", "ignore")
        for name in sorted(set(_PKG_RE.findall(html))):
            out.append((lic, name))
    return out


def load_pmcid_map(root: Path) -> dict:
    """PMCID -> PMID for every cancer article in EITHER census stream.

    Reads `records/` (MeSH-indexed), `records_unindexed/` (text-recovered) and
    `records_updates/` (the daily update stream), in whatever combination is
    present.
    Reading only the first silently excluded all 783,271 recovered articles,
    which is exactly the recent literature the recovery layer exists to reach.

    A RECENCY CEILING follows from this map, and it is structural rather than a
    bug. Full text can only be kept for an article the census already knows, so
    the PMC bulk cannot be matched past whatever the PubMed baseline contains.
    Measured on the 2026-06-17 bulk: both `PMC013xxxxxx` packages returned
    EXACTLY zero cancer articles out of 232,890, while every other package
    yielded 14-18%. The census's PMC identifier space stops at `PMC128xxxx`,
    so nothing in the `PMC13` block can match at all.

    That is a cliff, not the gradual decline MeSH indexing lag produces.

    IT IS NOW CLOSED, AND MEASURED RATHER THAN ESTIMATED. Closing it needed a
    newer PubMed baseline AND this map reading it: the update stream writes to
    its own directory precisely so it cannot mutate the frozen census, so a map
    reading only the first two directories never sees it. An earlier version of
    this docstring said "not a change here", which stopped being true the moment
    `atlas_baseline.py --updates` existed. Both halves were required.

    Re-running the two packages against the grown census:

        oa_comm     149,864 articles ->  9,651 cancer
        oa_noncomm   83,026 articles ->  6,612 cancer
        total       232,890 articles -> 16,263 cancer   (6.98%)

    against exactly 0 before. The full-text corpus went 1,100,218 -> 1,116,481.

    THE ESTIMATE THIS REPLACES WAS ABOUT TWICE TOO HIGH. It read "an estimated
    32,000-41,000 cancer full texts (232,890 at the 13.9-17.7% interquartile
    yield of the reachable packages)", which applies a yield rate measured on
    OTHER packages to this one. What can actually match is the census's own
    identifiers in the block, and holding a PMC id does not put an article in
    the open-access bulk: the census now holds 31,950 distinct PMC013xxxxxx
    identifiers, so 31,950 was the ceiling and 16,263 is the realised half of it.

    A NOTE ON READING THE BLOCK, because the first attempt got it backwards
    twice. Update files are numbered chronologically, so the recent end is the
    only part that carries these identifiers: across the 256-file window the
    oldest 30 files carry ZERO block records and the newest 30 carry 34,644,
    with 99,228 block records over 31,950 distinct ids in total. Sampling the
    oldest files says the block is unreachable, which is the opposite
    conclusion. And "starts with PMC13" is NOT the block: it also matches
    7-digit ids like PMC1349338 from 1988, which live in package PMC001xxxxxx.
    Compare numerically, 13,000,000 <= id < 14,000,000.
    """
    files = []
    for d in ("records", "records_unindexed", "records_updates"):
        files += sorted((root / d).glob("*.jsonl.gz"))
    if not files:
        raise SystemExit(f"no census under {root}; run scripts/atlas_baseline.py first")
    m = {}
    for f in files:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                pmc = (r.get("pmcid") or "").strip()
                if pmc:
                    m[pmc] = r.get("pmid", "")
    return m


def download(url: str, dest: Path) -> Path:
    """Resumable download via HTTP Range."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    pos = dest.stat().st_size if dest.exists() else 0
    headers = {"User-Agent": USER_AGENT}
    if pos:
        headers["Range"] = f"bytes={pos}-"
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=900) as resp:
            mode = "ab" if pos and resp.status == 206 else "wb"
            with open(dest, mode) as fh:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
    except urllib.error.HTTPError as e:
        if e.code != 416:  # 416 = already complete
            raise
    return dest


def process_package(tar_path: Path, licence: str, pmc2pmid: dict, out: Path) -> dict:
    """Stream a package, keep cancer members, write one gzipped JSONL shard."""
    kept = seen = 0
    out.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tf, gzip.open(out, "wt", encoding="utf-8") as w:
        for member in tf:
            if not member.isfile() or not member.name.endswith(".txt"):
                continue
            seen += 1
            m = _PMCID_RE.search(Path(member.name).name)
            if not m:
                continue
            pmcid = m.group(1)
            if pmcid not in pmc2pmid:
                continue
            fh = tf.extractfile(member)
            if fh is None:
                continue
            text = fh.read().decode("utf-8", "ignore")
            w.write(json.dumps({
                "pmcid": pmcid, "pmid": pmc2pmid[pmcid],
                "licence": licence, "text": text,
            }, ensure_ascii=False) + "\n")
            kept += 1
    return {"seen": seen, "kept": kept}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--limit", type=int, default=0, help="process at most N packages")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--keep-raw", action="store_true", help="keep downloaded tarballs")
    ap.add_argument("--only", default="", metavar="SUBSTR",
                    help="process only packages whose name contains SUBSTR")
    ap.add_argument("--redo", action="store_true",
                    help="reprocess selected packages even if already done; use "
                         "with --only after the census has grown, since a "
                         "package's yield depends on the census it was matched "
                         "against and a finished package is otherwise skipped")
    args = ap.parse_args()

    ft = fulltext_root()
    man_path = ft / "manifest.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"packages": {}}

    if args.status:
        done = [k for k, v in man["packages"].items() if v.get("done")]
        kept = sum(v.get("kept", 0) for v in man["packages"].values())
        seen = sum(v.get("seen", 0) for v in man["packages"].values())
        print(f"full-text root : {ft}")
        print(f"packages done  : {len(done)}")
        print(f"articles seen  : {seen:,}")
        print(f"cancer kept    : {kept:,}" + (f"  ({kept/seen:.1%})" if seen else ""))
        return

    root = atlas_root()
    print("loading PMCID -> PMID map from the census ...", flush=True)
    pmc2pmid = load_pmcid_map(root)
    print(f"cancer articles with a PMC id: {len(pmc2pmid):,}")

    ft.mkdir(parents=True, exist_ok=True)
    pkgs = list_packages()
    todo = [(l, n) for l, n in pkgs
            if args.redo or not man["packages"].get(n, {}).get("done")]
    if args.only:
        todo = [(l, n) for l, n in todo if args.only in n]
        if not todo:
            raise SystemExit(f"no package name contains {args.only!r}; "
                             f"available: {', '.join(n for _, n in pkgs)}")
    if args.limit:
        todo = todo[:args.limit]
    print(f"packages available: {len(pkgs)}; to process this run: {len(todo)}")

    raw = ft / "_raw"
    for i, (lic, name) in enumerate(todo, 1):
        t0 = time.time()
        url = f"{BULK}/{SUBSETS[lic]}/txt/{name}"
        tar_path = raw / name
        download(url, tar_path)
        shard = ft / "shards" / (name.replace(".tar.gz", "") + ".jsonl.gz")
        stats = process_package(tar_path, lic, pmc2pmid, shard)
        if not args.keep_raw:
            tar_path.unlink(missing_ok=True)
        stats.update(done=True, licence=lic, seconds=round(time.time() - t0, 1),
                     shard=shard.name)
        man["packages"][name] = stats
        man_path.write_text(json.dumps(man, indent=1, sort_keys=True), encoding="utf-8")
        print(f"  [{i}/{len(todo)}] {name}: {stats['seen']:,} articles -> "
              f"{stats['kept']:,} cancer, {stats['seconds']}s", flush=True)

    kept = sum(v.get("kept", 0) for v in man["packages"].values())
    print(f"\ncancer full texts stored: {kept:,}")


if __name__ == "__main__":
    main()
