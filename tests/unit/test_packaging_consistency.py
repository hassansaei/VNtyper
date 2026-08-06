"""Cross-file invariant: MANIFEST.in and pyproject package-data must agree.

Follows the pattern established by tests/unit/test_version_consistency.py.

Everything here compares **full package-relative paths**, never basenames. Comparing
basenames was the original bug in this guard: ``report_config.json`` lives at
``vntyper/scripts/report_config.json``, and a MANIFEST.in line pointing at
``vntyper/templates/report_config.json`` matched by basename while including nothing at
all. The sdist would have shipped without the file and this test would have stayed green.
"""

import re
from pathlib import Path, PurePosixPath

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]


def _manifest_included() -> set[str]:
    """Full package-relative paths of the JSON files MANIFEST.in ``include``s.

    Returns:
        set[str]: Posix-style paths exactly as written in MANIFEST.in.
    """
    text = (REPO_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    return {
        PurePosixPath(m.group(1)).as_posix() for m in re.finditer(r"^include\s+(\S+\.json)\s*$", text, re.MULTILINE)
    }


def _package_data_json() -> set[str]:
    """Full package-relative paths of the JSON files pyproject declares as package-data.

    ``[tool.setuptools.package-data]`` maps a *dotted package name* to patterns relative
    to that package's directory, so ``"vntyper.scripts" = ["kestrel_config.json"]`` means
    ``vntyper/scripts/kestrel_config.json``. Joining the two is what makes the comparison
    against MANIFEST.in meaningful.

    Returns:
        set[str]: Posix-style repo-relative paths.

    Raises:
        AssertionError: If the package-data table is missing from pyproject.toml.
    """
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    parts = text.split("[tool.setuptools.package-data]", 1)
    assert len(parts) == 2, "pyproject.toml has no [tool.setuptools.package-data] table; this guard cannot run."
    block = parts[1].split("\n[", 1)[0]

    paths: set[str] = set()
    for table in re.finditer(r'"([\w.]+)"\s*=\s*\[(.*?)\]', block, re.DOTALL):
        package_dir = PurePosixPath(table.group(1).replace(".", "/"))
        for entry in re.finditer(r'"([^"]+\.json)"', table.group(2)):
            paths.add((package_dir / entry.group(1)).as_posix())
    return paths


def test_every_package_data_json_is_in_manifest() -> None:
    """A JSON shipped by package-data but absent from MANIFEST.in may be missing from the sdist.

    Compared on the full package-relative path. A MANIFEST.in line that names the right
    basename under the wrong directory includes nothing, so matching on basename would
    let exactly that mistake through.

    Raises:
        AssertionError: If a package-data JSON has no matching MANIFEST.in ``include``.
    """
    pkg = _package_data_json()
    assert pkg, "No JSON files parsed out of [tool.setuptools.package-data]; this guard would assert nothing."

    manifest = _manifest_included()
    missing = sorted(pkg - manifest)
    assert not missing, (
        f"These files are in pyproject package-data but not MANIFEST.in: {missing}. "
        f"MANIFEST.in currently includes: {sorted(manifest)}. Add an `include <path>` line "
        "for each missing path - matching the basename under a different directory does not "
        "count, because MANIFEST.in resolves the path literally and silently includes nothing."
    )


def test_every_manifest_include_points_at_a_file_that_exists() -> None:
    """A MANIFEST.in ``include`` of a path that does not exist silently does nothing.

    setuptools only warns ("no previously-included files found"), so a mistyped directory
    never fails a build - it just quietly drops the file from the sdist.

    Raises:
        AssertionError: If MANIFEST.in includes a path that is not in the tree.
    """
    dangling = sorted(path for path in _manifest_included() if not (REPO_ROOT / path).is_file())
    assert not dangling, (
        f"MANIFEST.in includes path(s) that do not exist: {dangling}. setuptools only warns "
        "about these, so the sdist ships without the file. Fix the path."
    )
