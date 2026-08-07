"""On-disk cohort fixtures shared by the cohort discovery, identity and pseudonym tests.

`vntyper cohort` is pointed at directories and ZIP archives holding
`pipeline_summary.json` files, so almost every test of it begins by writing one of those.
These four builders are that beginning. They live here rather than in one test module
because four modules under `tests/unit/` need them - `test_cohort_pseudonym_config.py`,
`test_cohort_identity.py`, `test_cohort_zip_identity.py` and
`test_cohort_deduplication.py`, which were a single 1,210-line file until it was split
along its own seams.

`pseudonym_of` recomputes the pseudonym from `hashlib` rather than from the code under
test, and that is the point of it: an expectation built with the implementation cannot
fail when the implementation is wrong.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path


def summary_json(input_files: dict[str, str] | None = None) -> str:
    """Render a minimal ``pipeline_summary.json``.

    Args:
        input_files: The ``input_files`` mapping the run recorded.

    Returns:
        str: The serialised summary.
    """
    return json.dumps({"version": "2.0.6", "input_files": input_files or {}, "steps": []})


def sample_on_disk(directory: Path, input_files: dict[str, str] | None = None) -> Path:
    """Write one sample directory holding a ``pipeline_summary.json``.

    Args:
        directory: The sample directory to create.
        input_files: The ``input_files`` mapping the run recorded.

    Returns:
        Path: The directory that was created.
    """
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "pipeline_summary.json").write_text(summary_json(input_files), encoding="utf-8")
    return directory


def zip_of(archive: Path, members: dict[str, str]) -> Path:
    """Build a zip archive from a member -> contents mapping.

    Args:
        archive: Where to write the archive.
        members: Archive-relative path -> file contents.

    Returns:
        Path: The archive.
    """
    archive.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, "w") as handle:
        for member, contents in members.items():
            handle.writestr(member, contents)
    return archive


def pseudonym_of(identity: str) -> str:
    """The ``anon_`` pseudonym of one identity, recomputed from ``hashlib`` rather than
    from the code under test.

    Args:
        identity: The sample identity, qualified or not.

    Returns:
        str: ``anon_`` plus the first twelve hex characters of its SHA-256 digest.
    """
    return "anon_" + hashlib.sha256(identity.encode()).hexdigest()[:12]
