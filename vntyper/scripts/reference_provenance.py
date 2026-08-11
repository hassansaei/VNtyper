"""Per-file install provenance for the reference tree.

`canonical_reference_keys()` in `install_references.py` used to bless any file that
existed under `output_dir` into `config.json`, with no record of whether this or any
earlier `install-references` run had actually verified it. A partial install
(`--references hg19`) seeds its staging directory from the existing tree
(`reference_bundle.staged_install`, `seed_from_existing=True`) precisely so successive
`--references` runs accumulate into one tree - but that same seeding step also carried
forward files nobody involved in *any* run had ever verified: an old, hand-copied or
tampered FASTA sitting in the output directory blesses itself into `config.json` on the
next unrelated install.

This module is the fix: every verified file gets one record - its digest, where it came
from, and when - written into a small JSON ledger inside the reference tree itself
(`install_provenance.json`, at the tree's root). Because it lives inside the tree,
`staged_install`'s `shutil.copytree` seeding step carries earlier runs' records forward
automatically, which is exactly what keeps accumulation working: a `--references hg38`
run after a `--references hg19` run still finds hg19's provenance record without
repeating any verification, and a `--references hg38` run does not need to have ever
seen hg19's bytes to leave its record alone.

`canonical_reference_keys()` then requires a record for a path rather than mere
existence. This module never re-hashes a file to decide that - it only checks a JSON
key - so the ~700 MB genome files cost nothing extra on every `config.json` write. The
records themselves are written once, at install time, by callers who have already paid
for the digest verification that produced the value being recorded (a manifest digest
for the bundle path, a `source_sha256` for the from-source path); see
`install_references.py`'s `_install_bundle_asset` and `install_from_source`.

Imported function-locally everywhere it is used in `install_references.py`, never at
module scope - the same reason `reference_registry` is: `docker/Dockerfile.base` copies
only `install_references.py` and `reference_bundle.py` into the `refs` build stage and
runs the installer without installing the package, so a module-scope import of any other
sibling fails that stage's build. `tests/unit/test_docker_installer_imports.py` enforces
this.
"""

from __future__ import annotations

import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

#: Lives at the root of the installed reference tree, so `staged_install`'s
#: `shutil.copytree` seeding step (see module docstring) carries it forward exactly like
#: any other installed file. `scripts/bundle_release.py`'s `EXCLUSIONS` list keeps it out
#: of every published release asset - it is a ledger for the installer, never a release
#: artifact, and a from-source build run by `.github/workflows/build-reference-bundles.yml`
#: writes one into its working tree same as any other `install-references` run.
PROVENANCE_FILENAME = "install_provenance.json"


def provenance_path(root: Path) -> Path:
    """The ledger's path for a given reference-tree root.

    Args:
        root: Reference-tree root - an installed `output_dir`, or a `staged_install`
            staging directory mid-install.

    Returns:
        Path: `root / "install_provenance.json"`.
    """
    return root / PROVENANCE_FILENAME


def load_provenance(root: Path) -> dict[str, dict[str, Any]]:
    """Read the ledger, tolerating absence or corruption.

    Args:
        root: Reference-tree root.

    Returns:
        dict[str, dict[str, Any]]: Tree-relative POSIX path -> its provenance record.
        Empty when the ledger does not exist yet (first-ever install) or cannot be
        parsed. A corrupt ledger is never treated as fatal: it degrades to "no record
        for anything", which is the fail-closed direction `canonical_reference_keys`
        needs - it must never crash the whole install over a damaged provenance file.
    """
    path = provenance_path(root)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return {}
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        logger.warning(f"{path} is not valid JSON ({error}); treating it as an empty ledger")
        return {}
    if not isinstance(data, dict):
        logger.warning(f"{path} does not contain a JSON object; treating it as an empty ledger")
        return {}
    return data


def build_record(
    *,
    sha256: str,
    source: str,
    asset: str | None = None,
    release_tag: str | None = None,
    source_url: str | None = None,
) -> dict[str, Any]:
    """Build one provenance record.

    Args:
        sha256: The digest that was actually verified against this file's own bytes
            (the bundle path's per-file `release-manifest.json` entry), or against the
            upstream artifact this file was produced from (the from-source path's
            `source_sha256`).
        source: `"bundle"` or `"from-source"`.
        asset: Bundle asset file name. Bundle path only.
        release_tag: Bundle release tag. Bundle path only.
        source_url: Upstream URL the bytes were verified against. From-source path only.

    Returns:
        dict[str, Any]: One ledger entry, timestamped now (UTC, second precision).
    """
    return {
        "source": source,
        "sha256": sha256,
        "asset": asset,
        "release_tag": release_tag,
        "source_url": source_url,
        "installed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def merge(root: Path, new_records: dict[str, dict[str, Any]]) -> None:
    """Merge new records into the ledger and write it back atomically.

    A no-op when `new_records` is empty: it does not even read the existing ledger, so a
    caller that verified nothing this run (an empty `--references` selection resolving
    to no work) leaves an existing ledger untouched rather than rewriting it for nothing.

    Args:
        root: Reference-tree root the ledger lives under.
        new_records: Tree-relative POSIX path -> record, for files this run verified.
            Overwrites any existing record for the same path - the newest verification
            wins, matching install's own merge-not-replace contract for the files
            themselves.
    """
    if not new_records:
        return
    existing = load_provenance(root)
    existing.update(new_records)
    path = provenance_path(root)
    tmp_path = path.with_suffix(".json.tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(existing, handle, indent=2, sort_keys=True)
    os.replace(tmp_path, path)


def has_record(records: dict[str, dict[str, Any]], relative_path: str) -> bool:
    """Whether `relative_path` has a provenance record.

    Args:
        records: The ledger, as returned by :func:`load_provenance`.
        relative_path: Tree-relative POSIX path to check.

    Returns:
        bool: True if a record exists for that exact path.
    """
    return relative_path in records


def relative_posix(path: Path, root: Path) -> str:
    """Express `path` relative to `root` as a POSIX string, the ledger's key form.

    Args:
        path: Absolute or root-relative path to a file inside the tree.
        root: Reference-tree root.

    Returns:
        str: POSIX-separated path relative to `root`.

    Raises:
        ValueError: If `path` does not resolve under `root`.
    """
    return path.resolve().relative_to(root.resolve()).as_posix()
