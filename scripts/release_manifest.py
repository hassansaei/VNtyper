#!/usr/bin/env python3
"""Parse registry manifest descriptors for the release workflows."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

CANONICAL_DIGEST = re.compile(r"sha256:[0-9a-f]{64}")


def parse_manifest_digest(payload: object) -> str:
    """Return the canonical digest from an OCI manifest descriptor.

    Args:
        payload: Decoded JSON emitted for the Buildx ``.Manifest`` descriptor.

    Returns:
        The canonical lowercase SHA-256 digest.

    Raises:
        ValueError: If the payload is not a descriptor with a canonical digest.
    """
    if not isinstance(payload, Mapping):
        raise ValueError("manifest descriptor must contain a canonical sha256:<64 lowercase hex characters> digest")
    digest = payload.get("digest")
    if not isinstance(digest, str) or CANONICAL_DIGEST.fullmatch(digest) is None:
        raise ValueError("manifest digest must be canonical sha256:<64 lowercase hex characters>")
    return digest


def _digest(path: Path) -> int:
    try:
        with path.open(encoding="utf-8") as handle:
            payload: Any = json.load(handle)
        print(parse_manifest_digest(payload))
    except (OSError, json.JSONDecodeError, ValueError) as error:
        print(f"manifest descriptor parsing failed: {error}", file=sys.stderr)
        return 2
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the release manifest command-line parser.

    Returns:
        The configured argument parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    digest = subparsers.add_parser("digest", help="print the descriptor's canonical digest")
    digest.add_argument("path", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the release manifest command.

    Args:
        argv: Optional command arguments excluding the executable name.

    Returns:
        Zero on success and non-zero for invalid input.
    """
    args = build_parser().parse_args(argv)
    if args.command == "digest":
        return _digest(args.path)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
