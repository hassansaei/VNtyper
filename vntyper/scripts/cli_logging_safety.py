"""Pre-open ownership checks for the pipeline's application log."""

from __future__ import annotations

import argparse
import os
import stat
from pathlib import Path
from typing import Any

from vntyper.scripts.alignment_contract import index_candidate_names
from vntyper.scripts.alignment_target_io import bwa_index_paths, reference_index_paths
from vntyper.scripts.reference_resolution import configured_reference_candidates, resolve_from_mapping


def _absolute(path: str | Path) -> Path:
    return Path(os.path.abspath(path))


def _same_file(left: Path, right: Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def _is_within(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def _alignment_input_trees(args: argparse.Namespace) -> tuple[Path, ...]:
    trees: list[Path] = []
    for attribute in ("bam", "cram"):
        value = getattr(args, attribute, None)
        if value is None:
            continue
        absolute = _absolute(value)
        trees.extend((absolute.parent, absolute.resolve(strict=False).parent))
    return tuple(dict.fromkeys(trees))


def _selected_bwa_reference(args: argparse.Namespace, config: dict[str, Any]) -> Path | None:
    """The configured BWA reference path this guard must protect.

    Deliberately **not** the same call as `cli_handlers.select_bwa_reference(...,
    required=False)`: that helper's contract is to degrade a configured-but-missing
    file to None (`cli_handlers.py`, `TestThePresentButMissingFileFailsClosed`), which
    is correct for its own callers but was wrong here - it dropped a not-yet-installed
    reference out of the protected-path set entirely, so `--log-file` naming that exact
    path sailed through the guard and `setup_logging` created a regular file there. A
    configured-but-missing reference is exactly the case this guard exists to protect,
    so this resolves the same key (via the same `resolve_from_mapping` walk
    `select_bwa_reference` itself uses) **without** requiring the path to exist.

    Returns None only when nothing is configured for the assembly at all, or when the
    configured value is present-but-null (an explicit "disabled", not a path to guard)
    or not a string (malformed config, left to later validation - see
    `test_malformed_fastq_bwa_reference_is_left_to_pipeline_validation`).
    """
    assembly = args.reference_assembly or config.get("default_values", {}).get("reference_assembly", "hg19")
    # An unknown assembly still raises and MUST propagate: swallowing it fails the guard
    # open, and the guard runs before setup_logging opens the log file in append mode.
    resolved = resolve_from_mapping("bwa", assembly, config.get("reference_data", {}))
    if resolved is None or not isinstance(resolved.value, str) or not resolved.value:
        return None
    return Path(resolved.value)


def _pipeline_operator_paths(args: argparse.Namespace, config: dict[str, Any]) -> tuple[Path, ...]:
    paths: list[Path] = []
    for attribute in ("fastq1", "fastq2", "bed_file"):
        value = getattr(args, attribute, None)
        if value is not None:
            paths.append(Path(value))

    explicit_reference = getattr(args, "reference_fasta", None)
    if explicit_reference is not None:
        reference = Path(explicit_reference)
        paths.extend((reference, *reference_index_paths(reference)))

    for attribute, file_format in (("bam", "bam"), ("cram", "cram")):
        value = getattr(args, attribute, None)
        if value is None:
            continue
        alignment = Path(value)
        paths.append(alignment)
        paths.extend(Path(candidate) for candidate in index_candidate_names(str(alignment), file_format))

    if getattr(args, "cram", None) is not None:
        assembly = args.reference_assembly or config.get("default_values", {}).get("reference_assembly", "hg19")
        try:
            candidates = configured_reference_candidates(config, assembly)
        except (AttributeError, TypeError, ValueError):
            reference_data = config.get("reference_data", {})
            candidates = (
                tuple(
                    (key, value)
                    for key, value in reference_data.items()
                    if isinstance(key, str)
                    and key.startswith(("cram_reference_", "bwa_reference_"))
                    and isinstance(value, str)
                )
                if isinstance(reference_data, dict)
                else ()
            )
        for _source, path in candidates:
            if isinstance(path, str):
                reference = Path(path)
                paths.extend((reference, *reference_index_paths(reference)))

    if getattr(args, "fastq1", None) is not None or getattr(args, "fastq2", None) is not None:
        bwa_reference = _selected_bwa_reference(args, config)
        if bwa_reference is not None:
            paths.append(bwa_reference)
            paths.extend(bwa_index_paths(bwa_reference, config))
    return tuple(dict.fromkeys(paths))


def validate_pipeline_log_destination(
    log_file: str | Path,
    args: argparse.Namespace,
    config: dict[str, Any],
) -> None:
    """Reject a pipeline log without exclusive regular-file ownership.

    This check must run before the log parent is created or ``FileHandler`` opens
    the destination. Missing paths and single-link regular rerun logs are safe.
    Existing symlinks, non-regular entries, multiply linked files, and aliases of
    operator input state are rejected.

    Args:
        log_file: Explicit or default application-log destination.
        args: Parsed pipeline arguments.
        config: Loaded pipeline configuration used to select BWA sidecars.

    Raises:
        ValueError: If the log entry is unsafe or aliases an operator input.
    """
    log_path = Path(log_file)
    if os.path.lexists(log_path):
        metadata = os.lstat(log_path)
        if stat.S_ISLNK(metadata.st_mode):
            raise ValueError(f"Pipeline log file must not be a symlink: {log_path}")
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"Pipeline log file must be a regular file: {log_path}")
        if metadata.st_nlink > 1:
            raise ValueError(f"Pipeline log file has multiple hard links: {log_path}")
    log_variants = (_absolute(log_path), _absolute(log_path).resolve(strict=False))
    for protected in _pipeline_operator_paths(args, config):
        protected_absolute = _absolute(protected)
        protected_variants = (protected_absolute, protected_absolute.resolve(strict=False))
        if any(log_variant in protected_variants for log_variant in log_variants) or _same_file(log_path, protected):
            raise ValueError(f"Pipeline log file aliases operator-owned input: {log_path}")
    input_trees = _alignment_input_trees(args)
    if any(_is_within(log_variant, tree) for log_variant in log_variants for tree in input_trees):
        raise ValueError(f"Pipeline log file is inside an operator-owned input tree: {log_path}")
    # Names alone cannot see one host directory mounted at two container paths, and this
    # check runs before the log directory is created, so a lexical-only miss here leaves
    # an empty directory inside the patient input tree even when the run is refused (#238).
    for log_variant in log_variants:
        for ancestor in (log_variant, *log_variant.parents):
            for tree in input_trees:
                if ancestor != tree and _same_file(ancestor, tree):
                    raise ValueError(
                        f"Pipeline log file is inside an operator-owned input tree: {log_path} lies under "
                        f"{ancestor}, which is the same directory as the input tree {tree}."
                    )
