"""MkDocs macros hook — exposes project variables to all markdown pages."""

import re
from pathlib import Path


def define_env(env):
    """Read version from vntyper/version.py and expose as {{ version }}."""
    version_file = Path(env.project_dir) / "vntyper" / "version.py"
    version = "unknown"
    if version_file.exists():
        match = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', version_file.read_text())
        if match:
            version = match.group(1)

    env.variables["version"] = version
