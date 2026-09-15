"""Pure shell construction for legacy and approved calibrated adVNTR commands."""

from __future__ import annotations

import logging
import shlex

from vntyper.scripts.command_builders import quote_path

logger = logging.getLogger(__name__)


def build_advntr_command(
    executable: str,
    *,
    vid: int,
    alignment: str,
    result: str,
    model: str,
    working_directory: str,
    threads: int,
    additional_commands: str,
    calibrated: bool,
    calibrated_policy_arguments: tuple[str, ...] | None,
) -> str:
    """Build a native command while preserving the legacy fragment contract.

    Args:
        executable: Configured prefix; calibrated mode uses observed tokens.
        vid: Selected model locus.
        alignment: Prepared alignment path.
        result: Native result path.
        model: Run-owned model snapshot path.
        working_directory: Native stage directory.
        threads: Resolved native thread count.
        additional_commands: Legacy operator extension fragment.
        calibrated: Whether the resolved profile has calibrated native settings.
        calibrated_policy_arguments: Explicit arguments from verified native preflight.

    Returns:
        Shell command with every path and calibrated argument quoted as one token.

    Raises:
        ValueError: If calibrated arguments are missing, unapproved, or conflicting.
    """
    if calibrated != (calibrated_policy_arguments is not None):
        raise ValueError("calibrated native settings require their approved explicit policy arguments")
    if calibrated_policy_arguments is not None:
        if (
            not isinstance(calibrated_policy_arguments, tuple)
            or any(not isinstance(token, str) or not token for token in calibrated_policy_arguments)
            or calibrated_policy_arguments[:2] != ("-t", str(threads))
            or calibrated_policy_arguments.count("-t") != 1
            or additional_commands != ""
        ):
            raise ValueError("calibrated native arguments conflict with effective runtime settings")
        executable = shlex.join(shlex.split(executable))
        policy_fragment = shlex.join(calibrated_policy_arguments)
    else:
        additional_fragment = f" {additional_commands}" if additional_commands else ""
        policy_fragment = f"-t {quote_path(threads)}{additional_fragment}"
    return (
        f"{executable} genotype -fs -vid {quote_path(vid)} "
        f"--alignment_file {quote_path(alignment)} -o {quote_path(result)} "
        f"-m {quote_path(model)} --working_directory {quote_path(working_directory)} {policy_fragment}"
    )
