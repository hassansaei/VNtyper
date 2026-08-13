"""Cross-file invariant: MANIFEST.in and pyproject package-data must agree.

Follows the pattern established by tests/unit/test_version_consistency.py.

Everything here compares **full package-relative paths**, never basenames. Comparing
basenames was the original bug in this guard: ``report_config.json`` lives at
``vntyper/scripts/report_config.json``, and a MANIFEST.in line pointing at
``vntyper/templates/report_config.json`` matched by basename while including nothing at
all. The sdist would have shipped without the file and this test would have stayed green.

Why this guard grew past ``.json`` (#242)
-----------------------------------------
It used to parse ``.json`` entries only, and that is exactly the gap a release fell
through. ``vntyper/templates/_report_base.html`` -- the token layer both reports
``{% include %}`` -- was in neither MANIFEST.in nor package-data, so a wheel built from
that commit shipped two templates that could not render:

    jinja2.exceptions.TemplateNotFound: '_report_base.html' not found in search path

Measured, by building the wheel and installing it into a clean virtualenv. Nothing in
the repository could see it: an editable checkout renders from the working tree, where
the file is always present, and a ``.json``-only guard had nothing to say about an
``.html``.

So the guard now works in both directions, and over every extension:

* **declared implies shipped** -- every package-data pattern matches at least one real
  file, every file it matches is covered by MANIFEST.in, and every MANIFEST.in
  ``include`` names a path that exists;
* **present implies declared** -- every file under the two directories whose contents
  are pure payload (``vntyper/templates/`` and ``vntyper/assets/``) is named by
  package-data. Adding a template or an asset therefore fails a test rather than a
  release.

The second direction is the one that needs the explicit lists in pyproject.toml. A glob
would ship a new template silently, which is convenient right up to the release where it
ships something nobody meant to publish.
"""

import fnmatch
import re
from pathlib import Path, PurePosixPath

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Directories whose every file must be declared in package-data. These hold nothing but
#: payload the installed package reads at runtime, so "present but undeclared" is always
#: a packaging bug rather than a working file that belongs in the tree and not the wheel.
PAYLOAD_DIRECTORIES = ("vntyper/templates", "vntyper/assets")


def _manifest_lines() -> list[list[str]]:
    """Tokenise MANIFEST.in, ignoring blank lines and comments.

    Returns:
        list[list[str]]: One token list per directive line.
    """
    text = (REPO_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    lines = []
    for raw in text.splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        lines.append(stripped.split())
    return lines


def _manifest_included() -> set[str]:
    """Full package-relative paths named by MANIFEST.in ``include`` directives.

    Returns:
        set[str]: Posix-style paths exactly as written in MANIFEST.in.
    """
    return {
        PurePosixPath(token).as_posix()
        for tokens in _manifest_lines()
        if tokens[0] == "include"
        for token in tokens[1:]
    }


def _manifest_covers(path: str) -> bool:
    """Whether MANIFEST.in ships one repo-relative path.

    Both directives this repository uses are honoured: ``include`` names a path
    literally, and ``recursive-include <dir> <pattern>...`` matches any file at or below
    ``<dir>`` whose name matches one of the patterns.

    Args:
        path: A posix-style repo-relative path.

    Returns:
        bool: True when some directive ships it.
    """
    if path in _manifest_included():
        return True
    candidate = PurePosixPath(path)
    for tokens in _manifest_lines():
        if tokens[0] != "recursive-include" or len(tokens) < 3:
            continue
        directory = PurePosixPath(tokens[1])
        if not candidate.is_relative_to(directory):
            continue
        if any(fnmatch.fnmatch(candidate.name, pattern) for pattern in tokens[2:]):
            return True
    return False


def _package_data_patterns() -> dict[str, list[str]]:
    """The package-data table, as repo-relative patterns keyed by the package that owns them.

    ``[tool.setuptools.package-data]`` maps a *dotted package name* to patterns relative
    to that package's directory, so ``"vntyper.scripts" = ["kestrel_config.json"]`` means
    ``vntyper/scripts/kestrel_config.json``. Joining the two is what makes the comparison
    against MANIFEST.in meaningful.

    Returns:
        dict[str, list[str]]: Dotted package name -> its repo-relative patterns.

    Raises:
        AssertionError: If the package-data table is missing from pyproject.toml.
    """
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    parts = text.split("[tool.setuptools.package-data]", 1)
    assert len(parts) == 2, "pyproject.toml has no [tool.setuptools.package-data] table; this guard cannot run."
    block = parts[1].split("\n[", 1)[0]

    table: dict[str, list[str]] = {}
    for match in re.finditer(r'"([\w.]+)"\s*=\s*\[(.*?)\]', block, re.DOTALL):
        package_dir = PurePosixPath(match.group(1).replace(".", "/"))
        entries = re.findall(r'"([^"]+)"', match.group(2))
        table[match.group(1)] = [(package_dir / entry).as_posix() for entry in entries]
    return table


def _package_data_files() -> set[str]:
    """Every real file the package-data patterns expand to.

    Returns:
        set[str]: Posix-style repo-relative paths of existing files.
    """
    found: set[str] = set()
    for patterns in _package_data_patterns().values():
        for pattern in patterns:
            found.update(path.relative_to(REPO_ROOT).as_posix() for path in REPO_ROOT.glob(pattern) if path.is_file())
    return found


def test_every_package_data_pattern_matches_a_real_file() -> None:
    """A pattern that matches nothing declares nothing, and setuptools does not complain.

    This is the guard that makes the two below non-vacuous: without it, deleting a file
    and leaving its pyproject entry behind would empty the comparison set instead of
    failing.

    Raises:
        AssertionError: If a package-data pattern matches no file in the tree.
    """
    table = _package_data_patterns()
    assert table, "No [tool.setuptools.package-data] entries parsed; this guard would assert nothing."

    empty = sorted(pattern for patterns in table.values() for pattern in patterns if not list(REPO_ROOT.glob(pattern)))
    assert not empty, (
        f"These pyproject package-data patterns match no file in the tree: {empty}. "
        "setuptools ships nothing for them and says nothing about it."
    )


def test_every_package_data_file_is_shipped_by_the_manifest() -> None:
    """A file shipped by package-data but absent from MANIFEST.in may be missing from the sdist.

    Compared on the full package-relative path. A MANIFEST.in line that names the right
    basename under the wrong directory includes nothing, so matching on basename would
    let exactly that mistake through.

    Raises:
        AssertionError: If a package-data file has no matching MANIFEST.in directive.
    """
    files = _package_data_files()
    assert files, "No files expanded out of [tool.setuptools.package-data]; this guard would assert nothing."

    missing = sorted(path for path in files if not _manifest_covers(path))
    assert not missing, (
        f"These files are in pyproject package-data but not shipped by MANIFEST.in: {missing}. "
        f"MANIFEST.in currently includes: {sorted(_manifest_included())}. Add an `include <path>` line "
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


@pytest.mark.parametrize("directory", PAYLOAD_DIRECTORIES)
def test_every_payload_file_is_declared_in_package_data(directory: str) -> None:
    """Adding a template or an asset must fail a test rather than a release.

    This is the direction the ``.json``-only guard could not see, and the direction the
    ``_report_base.html`` omission fell through: the file was in the tree, both reports
    included it, every test rendered it from the working tree, and no distribution
    carried it.

    Args:
        directory: A repo-relative directory whose every file must be declared.

    Raises:
        AssertionError: If a file under ``directory`` is not named by package-data.
    """
    root = REPO_ROOT / directory
    if not root.is_dir():
        pytest.skip(f"{directory} does not exist in this checkout")

    present = {
        path.relative_to(REPO_ROOT).as_posix()
        for path in root.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    }
    assert present, f"{directory} is empty, so this guard would assert nothing."

    undeclared = sorted(present - _package_data_files())
    assert not undeclared, (
        f"These files are in {directory} but are not declared in pyproject's "
        f"[tool.setuptools.package-data]: {undeclared}. A wheel or sdist built from this "
        "commit ships without them, and the failure only appears at runtime in an "
        "installed package. Name each one in package-data and in MANIFEST.in."
    )
