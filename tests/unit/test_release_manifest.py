"""Pure and executable contracts for OCI manifest descriptor parsing."""

import importlib
import json
import subprocess
import sys
from collections.abc import Callable, Mapping
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]
MANIFEST_SCRIPT = ROOT / "scripts" / "release_manifest.py"
DIGEST = "sha256:" + "b" * 64


def _parser() -> Callable[[object], str]:
    assert MANIFEST_SCRIPT.is_file(), "the fail-closed manifest descriptor parser has not been implemented"
    module = importlib.import_module("scripts.release_manifest")
    return module.parse_manifest_digest


def test_manifest_descriptor_returns_one_canonical_digest() -> None:
    """The registry identity comes only from the descriptor's canonical digest field."""
    parse_manifest_digest = _parser()

    assert (
        parse_manifest_digest({"mediaType": "application/vnd.oci.image.manifest.v1+json", "digest": DIGEST}) == DIGEST
    )


@pytest.mark.parametrize(
    "payload",
    (
        None,
        [],
        {},
        {"digest": None},
        {"digest": ""},
        {"digest": "sha512:" + "b" * 64},
        {"digest": "sha256:" + "b" * 63},
        {"digest": "sha256:" + "b" * 65},
        {"digest": "sha256:" + "B" * 64},
        {"digest": DIGEST + "\nDigest: sha256:" + "c" * 64},
    ),
)
def test_manifest_descriptor_rejects_missing_malformed_and_non_sha256_identity(payload: object) -> None:
    """Missing, ambiguous, and noncanonical identities must all fail closed."""
    parse_manifest_digest = _parser()

    with pytest.raises(ValueError, match="canonical sha256"):
        parse_manifest_digest(payload)


@pytest.mark.parametrize(
    "source",
    (
        "",
        "not-json",
        (
            "Name:      ghcr.io/hassansaei/vntyper:main\n"
            "MediaType: application/vnd.docker.distribution.manifest.v2+json\n"
            f"Digest:    {DIGEST}\n"
        ),
        json.dumps({"digest": None}),
    ),
)
def test_manifest_cli_rejects_legacy_report_and_malformed_descriptor(tmp_path: Path, source: str) -> None:
    """The old three-line Buildx report must never be accepted as a digest."""
    assert MANIFEST_SCRIPT.is_file(), "the fail-closed manifest descriptor parser has not been implemented"
    descriptor = tmp_path / "manifest.json"
    descriptor.write_text(source, encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, str(MANIFEST_SCRIPT), "digest", str(descriptor)],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 2
    assert completed.stdout == ""
    assert "manifest descriptor parsing failed" in completed.stderr


def test_manifest_cli_prints_exact_digest_without_other_buildx_output(tmp_path: Path) -> None:
    """A valid descriptor produces exactly one shell-safe digest line."""
    assert MANIFEST_SCRIPT.is_file(), "the fail-closed manifest descriptor parser has not been implemented"
    descriptor = tmp_path / "manifest.json"
    payload: Mapping[str, object] = {"digest": DIGEST, "size": 3464}
    descriptor.write_text(json.dumps(payload), encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, str(MANIFEST_SCRIPT), "digest", str(descriptor)],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout == DIGEST + "\n"
    assert completed.stderr == ""
