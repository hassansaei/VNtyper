"""Check that the files Kestrel needs exist before the pipeline does any work (#338).

Every tool and reference path in ``config.json`` resolves against the process CWD
(AGENTS.md trap 7). A source checkout where ``vntyper install-references`` was never run,
or a job script that ``cd``s somewhere else, therefore has no motif reference. In fast
mode with a CRAM input nothing else under ``reference/`` is read, so nothing noticed until
Kestrel itself, which logs the error and exits ``0`` (see ``kestrel_log_contract``). The
run then failed on every sample with "no usable VCF", a message that reads like a data
problem.

This module holds the decision (which inputs, resolved where, and what to tell the
operator) and does no I/O beyond ``is_file``, so it is testable with ``tmp_path``.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from vntyper.scripts.kestrel_counting import DEFAULT_KANALYZE_PATH

logger = logging.getLogger(__name__)

#: ``reference_data`` keys Kestrel reads. The motif FASTA is its reference; the reverse
#: complement motif file is read by the post-processing that turns its VCF into calls.
REFERENCE_KEYS: tuple[str, ...] = ("muc1_reference_vntr", "muc1_motifs_rev_com")


@dataclass(frozen=True)
class KestrelInput:
    """One file the Kestrel stage reads.

    Attributes:
        key: The config key that declares it, e.g. ``reference_data.muc1_reference_vntr``.
        configured: The value as written in the config.
        resolved: The absolute path it resolves to against the run's CWD.
    """

    key: str
    configured: str
    resolved: Path


def required_kestrel_inputs(config: Mapping[str, Any], project_root: str | Path, *, split: bool) -> list[KestrelInput]:
    """List the files the Kestrel stage will read, resolved the way the stage resolves them.

    Args:
        config: The loaded ``config.json``.
        project_root: The CWD captured at pipeline entry; relative paths resolve here.
        split: Whether k-mer counting runs as a separate KAnalyze step, which adds
            ``kanalyze.jar`` to what must exist.

    Returns:
        list[KestrelInput]: In a stable order: tool JARs, then references. A key the
        config leaves unset or empty is listed with an empty ``configured`` value so the
        caller reports it as missing rather than skipping it.
    """
    tools = config.get("tools", {})
    references = config.get("reference_data", {})
    declared: list[tuple[str, object]] = [("tools.kestrel", tools.get("kestrel"))]
    if split:
        declared.append(("tools.kanalyze", tools.get("kanalyze", DEFAULT_KANALYZE_PATH)))
    declared.extend((f"reference_data.{key}", references.get(key)) for key in REFERENCE_KEYS)

    root = Path(project_root)
    inputs = []
    for key, value in declared:
        configured = value if isinstance(value, str) else ""
        resolved = (root / configured).resolve() if configured else root
        inputs.append(KestrelInput(key=key, configured=configured, resolved=resolved))
    return inputs


def missing_kestrel_inputs(inputs: list[KestrelInput]) -> list[KestrelInput]:
    """Return the inputs that are unset or are not a regular file.

    Args:
        inputs: From :func:`required_kestrel_inputs`.

    Returns:
        list[KestrelInput]: The missing ones, in input order.
    """
    return [item for item in inputs if not item.configured or not item.resolved.is_file()]


def describe_missing_inputs(missing: list[KestrelInput], project_root: str | Path) -> str:
    """Render the operator-facing message for missing Kestrel inputs.

    Args:
        missing: Non-empty output of :func:`missing_kestrel_inputs`.
        project_root: The CWD relative paths were resolved against.

    Returns:
        str: A message naming each file, where it was looked for, and the fix.

    Raises:
        ValueError: If ``missing`` is empty.
    """
    if not missing:
        raise ValueError("describe_missing_inputs needs at least one missing input.")
    lines = []
    for item in missing:
        if item.configured:
            lines.append(f"  - {item.key} = {item.configured!r} -> {item.resolved} (not found)")
        else:
            lines.append(f"  - {item.key} is not set in the config")
    fix = []
    if any(item.key.startswith("reference_data.") for item in missing):
        fix.append(
            "Install the reference bundle into this directory with "
            "`vntyper install-references --output-dir reference`, or run VNtyper from the directory "
            "where it is installed, or set these keys to absolute paths in the file passed to --config-path."
        )
    if any(item.key.startswith("tools.") for item in missing):
        fix.append(
            "The Kestrel JARs ship inside the VNtyper source tree under vntyper/dependencies/kestrel/; "
            "run from the repository root or set tools.kestrel / tools.kanalyze to absolute paths."
        )
    return (
        "Kestrel cannot run because required files are missing. Relative paths in config.json are "
        f"resolved against the current working directory ({Path(project_root)}):\n"
        + "\n".join(lines)
        + "\n"
        + "\n".join(fix)
    )


def check_kestrel_inputs(
    config: Mapping[str, Any],
    project_root: str | Path,
    *,
    runtime_component: Mapping[str, Any] | None = None,
) -> None:
    """Raise before any pipeline work if a file Kestrel needs does not exist.

    Args:
        config: The loaded ``config.json``.
        project_root: The CWD captured at pipeline entry.
        runtime_component: The run's resolved Kestrel runtime sidecar, or None for the
            packaged one. It decides whether counting is split into a KAnalyze step, and
            so whether ``kanalyze.jar`` must exist; the stage resolves it the same way.

    Raises:
        ValueError: Naming every missing file and how to provide it. A validation
            failure, like a missing BWA reference, so the CLI exits 1 without a traceback.
    """
    from vntyper.scripts.pipeline_kestrel import resolve_kestrel_counting_mode

    split = resolve_kestrel_counting_mode(runtime_component, config) == "split"
    missing = missing_kestrel_inputs(required_kestrel_inputs(config, project_root, split=split))
    if missing:
        msg = describe_missing_inputs(missing, project_root)
        logger.error(msg)
        raise ValueError(msg)
