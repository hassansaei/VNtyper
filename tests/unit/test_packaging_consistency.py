"""Cross-file invariant: MANIFEST.in and pyproject package-data must agree.

Follows the pattern established by tests/unit/test_version_consistency.py.
"""

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]


def _manifest_included() -> set[str]:
    text = (REPO_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    return {m.group(1) for m in re.finditer(r"^include\s+(\S+\.json)\s*$", text, re.MULTILINE)}


def _package_data_json() -> set[str]:
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    block = text.split("[tool.setuptools.package-data]", 1)[1].split("\n[", 1)[0]
    names = set()
    for m in re.finditer(r'"([^"]+\.json)"', block):
        names.add(m.group(1))
    return names


def test_every_package_data_json_is_in_manifest() -> None:
    """A JSON shipped by package-data but absent from MANIFEST.in may be missing from the sdist."""
    pkg = _package_data_json()
    manifest_basenames = {Path(p).name for p in _manifest_included()}
    missing = sorted(name for name in pkg if Path(name).name not in manifest_basenames)
    assert not missing, (
        f"These files are in pyproject package-data but not MANIFEST.in: {missing}. "
        "Add an `include` line for each, or the sdist may ship without them."
    )
