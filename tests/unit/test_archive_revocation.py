"""Revocation-only safeguards for failed public archives."""

from pathlib import Path

import pytest

from vntyper.scripts.archive_safety import revoke_public_archive

pytestmark = pytest.mark.unit


def test_revocation_refuses_to_unlink_the_exact_protected_directory_entry(tmp_path: Path) -> None:
    """A caller cannot use revocation to remove the protected path itself."""
    protected = tmp_path / "patient.zip"
    protected.write_bytes(b"patient bytes")

    with pytest.raises(ValueError, match="revocation targets protected input"):
        revoke_public_archive(tmp_path / "patient", "zip", protected_paths=(protected,))

    assert protected.read_bytes() == b"patient bytes"
