"""Version discovery and command construction for the SHARK read filter (#312).

SHARK 1.2.0 emits no version information from the binary itself (neither
``--version`` nor ``--help`` reports a version string). For conda/mamba environments,
the package version and build string are determined from conda's JSON environment
listing (e.g. ``mamba list -n shark_env shark --json``).
"""

from __future__ import annotations

import json
import logging
import re
import shlex

logger = logging.getLogger(__name__)


def extract_conda_env_from_command(command: str) -> tuple[str, str] | None:
    """Extract the conda/mamba executable and environment name from a command prefix.

    Args:
        command: Command prefix string (e.g. ``"mamba run -n shark_env shark"``).

    Returns:
        A tuple of ``(conda_executable, env_name)`` if the command uses a recognized
        conda run prefix, or ``None`` if it is a bare command or unsupported launcher.
    """
    try:
        tokens = shlex.split(command)
    except ValueError:
        return None

    if len(tokens) < 4:
        return None

    executable = tokens[0]
    launcher_base = executable.split("/")[-1].lower()
    if launcher_base not in {"mamba", "conda", "micromamba"}:
        return None

    if tokens[1] != "run":
        return None

    env_name = None
    i = 2
    while i < len(tokens) - 1:
        if tokens[i] in {"-n", "--name"}:
            env_name = tokens[i + 1]
            break
        i += 1

    if env_name is None:
        return None

    return executable, env_name


def build_shark_version_command(command: str) -> list[str] | None:
    """Construct the command argv to query the SHARK version via conda listing.

    Args:
        command: Configured tool command string for SHARK.

    Returns:
        List of command tokens to execute, or ``None`` if command is not a conda
        run prefix.
    """
    env_info = extract_conda_env_from_command(command)
    if env_info is None:
        return None
    executable, env_name = env_info
    return [executable, "list", "-n", env_name, "shark", "--json"]


def parse_shark_conda_list_json(json_text: str) -> str | None:
    """Parse the version and build string from `conda list --json` output.

    Args:
        json_text: JSON output string from ``conda list ... --json``.

    Returns:
        Version string formatted as ``<version>+<build>`` (e.g. ``"1.2.0+h077b44d_5"``),
        or ``<version>`` if build string is absent, or ``None`` if parsing fails.
    """
    text = json_text.strip()
    if not text:
        return None

    try:
        data = json.loads(text)
    except (json.JSONDecodeError, ValueError):
        return None

    if not isinstance(data, list) or not data:
        return None

    for record in data:
        if not isinstance(record, dict):
            continue
        if record.get("name") == "shark":
            version = record.get("version")
            if not version:
                return None
            build = record.get("build_string")
            if build:
                return f"{version}+{build}"
            return str(version)

    first = data[0]
    if isinstance(first, dict) and "version" in first:
        version = first["version"]
        build = first.get("build_string")
        if build:
            return f"{version}+{build}"
        return str(version)

    return None


def parse_shark_help_output(output: str) -> str | None:
    """Parse version information from SHARK help or stdout text if present.

    Args:
        output: Raw text output from `shark --version` or `shark --help`.

    Returns:
        Extracted version string, or ``None`` if no version was reported.
    """
    text = output.strip()
    if not text:
        return None

    match = re.search(r"\bshark\s+(?:version\s+)?(\d+\.\d+(?:\.\d+)?(?:[_\+\-][a-zA-Z0-9_]+)?)", text, re.IGNORECASE)
    if match:
        return match.group(1)

    return None


def resolve_shark_version_from_output(json_output: str) -> str:
    """Resolve SHARK version string from conda JSON output, defaulting to 'unknown'.

    Args:
        json_output: JSON output from conda list.

    Returns:
        Parsed version string or ``"unknown"``.
    """
    return parse_shark_conda_list_json(json_output) or "unknown"
