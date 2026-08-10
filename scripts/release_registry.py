"""Fail-closed classification of GHCR manifest inspection errors."""

import argparse
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

ManifestAbsence = Literal["absent", "ambiguous"]

_NAME_COMPONENT = r"[a-z0-9]+(?:(?:[._]|__|[-]+)[a-z0-9]+)*"
_TAG = r"[A-Za-z0-9_][A-Za-z0-9._-]{0,127}"
_GHCR_REFERENCE = re.compile(rf"ghcr\.io/(?P<repository>{_NAME_COMPONENT}(?:/{_NAME_COMPONENT})+):(?P<tag>{_TAG})")
_ABSENCE_STATUS = r"(?i:(?:404(?:[ \t]+not(?:[ \t]+|-)found)?|manifest(?:[ \t]+|[_-])unknown|not(?:[ \t]+|-)found))"


def classify_manifest_absence(reference: str, error_text: str) -> ManifestAbsence:
    """Classify whether one inspection error proves an exact GHCR tag is absent.

    A verdict is authoritative only when one physical error line binds the exact
    requested tag or its derived GHCR manifest URL directly to a recognized
    absence status. All other output is ambiguous and therefore fail-closed.

    Args:
        reference: Exact tagged GHCR reference that was inspected.
        error_text: Complete standard error emitted by that inspection attempt.

    Returns:
        ``"absent"`` for an authoritative bound absence record, otherwise
        ``"ambiguous"``.

    Raises:
        ValueError: If ``reference`` is not one well-formed tagged GHCR reference.
    """
    match = _GHCR_REFERENCE.fullmatch(reference)
    if match is None:
        msg = f"Requested reference {reference!r} must be one exact tagged GHCR reference."
        raise ValueError(msg)
    if "\x00" in error_text:
        return "ambiguous"

    repository = match.group("repository")
    tag = match.group("tag")
    manifest_url = f"https://ghcr.io/v2/{repository}/manifests/{tag}"
    escaped_targets = "|".join(re.escape(target) for target in (reference, manifest_url))
    bound_absence = re.compile(
        rf"(?<![A-Za-z0-9._/-])(?:{escaped_targets})(?![A-Za-z0-9._-])"
        rf"[ \t]*:[ \t]*{_ABSENCE_STATUS}(?![A-Za-z0-9_-])"
    )
    return "absent" if any(bound_absence.search(line) is not None for line in error_text.splitlines()) else "ambiguous"


def main(argv: Sequence[str] | None = None) -> int:
    """Run the release-registry classifier command.

    Args:
        argv: Optional command arguments excluding the executable name.

    Returns:
        Zero only for authoritative absence, one for ambiguity, and two for
        malformed input or an unreadable error file.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    classify = subparsers.add_parser("classify-absence")
    classify.add_argument("reference")
    classify.add_argument("error_path", type=Path)
    args = parser.parse_args(argv)

    try:
        error_text = args.error_path.read_text(encoding="utf-8")
        verdict = classify_manifest_absence(args.reference, error_text)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"release registry classification failed: {error}", file=sys.stderr)
        return 2
    return 0 if verdict == "absent" else 1


if __name__ == "__main__":
    raise SystemExit(main())
