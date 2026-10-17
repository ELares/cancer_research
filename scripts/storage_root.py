#!/usr/bin/env python3
"""One refusal, shared by every tool that writes to externally-configured storage.

A MOUNT POINT IS JUST A DIRECTORY. With a network share unmounted, its
mountpoint is an ordinary empty local directory, so

    root.mkdir(parents=True, exist_ok=True)

cheerfully creates the whole path on the BOOT DISK and the caller writes its
corpus into the wrong filesystem, silently, at local-disk speed, with no error
anywhere. The ledgers that record what was written live in the repository
rather than on the share, so they go on reporting success -- and the bytes are
then invisible to every later run, because the run looks for them on the share.

Measured (#960): the expansion crawl's shard root and the open-access
full-text root both default under the same mountpoint, and both created it
with `parents=True`. The share holds 51 GB against 16 GB free locally.

THE RULE IS NARROW ON PURPOSE: a process may create its OWN leaf and nothing
above it. On a mounted share the parent exists; on an unmounted mountpoint it
does not. That is one `stat`, it needs no knowledge of which paths are mounts,
and it is NOT a check for "is this a network filesystem" -- pointing a root at
a local directory is legitimate and every test does it, so a rule that
refused local paths would be enforcing a storage preference rather than
catching a mistake.
"""
from __future__ import annotations

from pathlib import Path


def require_output_root(root: Path, env_var: str = "FERRO_EXPAND_OUT") -> None:
    """Create `root`, refusing to invent the directories above it.

    Raises SystemExit with an actionable message when `root.parent` is absent,
    which is what an unmounted mountpoint looks like from inside the process.
    """
    root = Path(root)
    if not root.parent.is_dir():
        raise SystemExit(
            f"REFUSING to run: {root.parent} does not exist, so {root} cannot "
            f"be created without inventing the path to it.\n"
            f"If that is a mount point, the share is not mounted -- mount it "
            f"and re-run; progress ledgers are local, so work resumes.\n"
            f"If the output really should live somewhere else, create the "
            f"parent yourself or set {env_var}.")
    root.mkdir(exist_ok=True)
