# vntyper/scripts/verify_seeds.py

"""Check the data repository's committed seeds against the release spec.

The three non-derivable artefacts (`MUC1_motifs_Rev_com.fa`, `code-adVNTR_RUs.fa`,
`vntr_db_advntr.zip`) plus `filter_config.json` are the only reference bytes a bundle
build cannot reproduce from an upstream source, so they are the only ones whose
provenance rests on a commit rather than on a download. The release spec pins each by
SHA-256; this is what turns that pin into a gate.

It runs **first** in the bundle build - before six genomes are downloaded and BWA-indexed
- so a tampered or mis-committed seed costs a minute rather than three hours, and a
release is never cut from bytes the spec does not describe.

This module is deliberately kept out of `docker/Dockerfile.base`'s `refs` stage: it is a
release-build tool, not part of the installer. It may therefore import
`reference_bundle` freely, but nothing the Dockerfile *does* copy may import **this**
module - see the AST guard in `tests/unit/test_docker_installer_imports.py`.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from vntyper.scripts.reference_bundle import verify_sha256

logger = logging.getLogger(__name__)


def load_spec(spec_path: Path) -> dict[str, Any]:
    """Read a committed release spec.

    Args:
        spec_path: Path to the caller repository's `releases/<tag>.json`.

    Returns:
        dict[str, Any]: The parsed document.

    Raises:
        ValueError: If the file is absent, unreadable or not a JSON object. The path is
            named in the message, because the spec lives in a different repository from
            this code and "which file" is the first question asked.
    """
    try:
        document = json.loads(spec_path.read_text(encoding="utf-8"))
    except OSError as error:
        message = f"cannot read release spec {spec_path}: {error}"
        logger.error(message)
        raise ValueError(message) from error
    except json.JSONDecodeError as error:
        message = f"release spec {spec_path} is not valid JSON: {error}"
        logger.error(message)
        raise ValueError(message) from error

    if not isinstance(document, dict):
        message = f"release spec {spec_path} must be a JSON object, got {type(document).__name__}"
        logger.error(message)
        raise ValueError(message)
    return document


def seed_digests(spec: dict[str, Any]) -> dict[str, str]:
    """Extract the `seeds` block as a name -> digest mapping.

    Both spellings the design uses are accepted: `{"name": {"sha256": "..."}}` and the
    flat `{"name": "..."}`.

    Args:
        spec: A parsed release spec.

    Returns:
        dict[str, str]: Seed file name -> expected lowercase hex SHA-256.

    Raises:
        ValueError: If the block is absent, empty, not a mapping, or if any entry
            declares no digest. An empty block would let the build proceed having
            verified nothing, which is worse than failing.
    """
    seeds = spec.get("seeds")
    if not isinstance(seeds, dict) or not seeds:
        message = "release spec declares no non-empty 'seeds' block; refusing to build unverified seeds"
        logger.error(message)
        raise ValueError(message)

    digests: dict[str, str] = {}
    for name, entry in seeds.items():
        digest = entry.get("sha256") if isinstance(entry, dict) else entry
        if not isinstance(digest, str) or not digest:
            message = f"seed '{name}' declares no sha256 in the release spec; refusing to accept it unverified"
            logger.error(message)
            raise ValueError(message)
        digests[name] = digest
    return digests


def resolve_seed(seeds_dir: Path, name: str) -> Path:
    """Turn a spec-supplied seed name into a path inside the seeds directory.

    The spec is committed data from *another* repository, so its keys are treated as
    untrusted path fragments: an entry naming `../something` would otherwise make this
    verify - and the build stage, - a file outside the directory it was pointed at.

    Args:
        seeds_dir: Directory holding the committed seeds.
        name: Seed name as written in the spec.

    Returns:
        Path: The seed's location.

    Raises:
        ValueError: If the name is absolute or escapes `seeds_dir`.
    """
    candidate = Path(name)
    if candidate.is_absolute() or ".." in candidate.parts:
        message = f"seed name '{name}' in the release spec escapes {seeds_dir}"
        logger.error(message)
        raise ValueError(message)
    return seeds_dir / candidate


def verify_seeds(spec: dict[str, Any], seeds_dir: Path) -> list[str]:
    """Verify every seed the spec names, failing on the first that does not match.

    Args:
        spec: A parsed release spec.
        seeds_dir: Directory holding the committed seeds.

    Returns:
        list[str]: The seed names verified, in spec order.

    Raises:
        ValueError: On a malformed `seeds` block, an escaping name, or a digest
            mismatch (raised by `verify_sha256`).
        FileNotFoundError: If a seed the spec names is not in `seeds_dir`. Named
            explicitly, because the usual cause is a staging step that forgot a file.
    """
    verified: list[str] = []
    for name, expected in seed_digests(spec).items():
        seed = resolve_seed(seeds_dir, name)
        if not seed.is_file():
            message = f"seed '{name}' named by the release spec is not present in {seeds_dir}"
            logger.error(message)
            raise FileNotFoundError(message)
        verify_sha256(seed, expected)
        verified.append(name)
    return verified


def main(argv: Sequence[str] | None = None) -> int:
    """Verify the seeds and report the outcome as a process exit code.

    Args:
        argv: Command-line arguments; `sys.argv[1:]` when None.

    Returns:
        int: 0 when every seed matched, 1 otherwise. Usage errors exit with argparse's 2.
    """
    parser = argparse.ArgumentParser(
        description="Verify the committed reference seeds against a release spec.",
    )
    parser.add_argument("--spec", type=Path, required=True, help="Committed release spec (releases/<tag>.json).")
    parser.add_argument("--seeds", type=Path, required=True, help="Directory holding the committed seeds.")
    args = parser.parse_args(argv)

    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    try:
        verified = verify_seeds(load_spec(args.spec), args.seeds)
    except (ValueError, FileNotFoundError, OSError) as error:
        logger.error(f"seed verification failed: {error}")
        return 1

    logger.info(f"Verified {len(verified)} seed(s) against {args.spec}: {', '.join(verified)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
