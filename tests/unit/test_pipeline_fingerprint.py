"""The background fingerprint helper hashes regular files and nothing else."""

import hashlib

import pytest

from vntyper.scripts.pipeline import _fingerprint_if_file

pytestmark = pytest.mark.unit


def test_a_regular_file_gets_its_sha256(tmp_path):
    path = tmp_path / "input.bam"
    path.write_bytes(b"BAM\x01" * 1000)

    assert _fingerprint_if_file(str(path)) == hashlib.sha256(b"BAM\x01" * 1000).hexdigest()


@pytest.mark.parametrize("make_path", [lambda tmp: None, lambda tmp: str(tmp / "missing.fa"), lambda tmp: str(tmp)])
def test_no_path_a_missing_file_or_a_directory_gives_none(tmp_path, make_path):
    assert _fingerprint_if_file(make_path(tmp_path)) is None
