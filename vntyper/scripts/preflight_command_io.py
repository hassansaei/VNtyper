"""Safe captured-command execution for alignment preflight."""

from __future__ import annotations

import logging
import os
import signal
import subprocess
import tempfile
from collections.abc import Iterable
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class CaptureMetadata:
    """Structured facts established by the command runner, never parsed from output."""

    timed_out: bool = False
    timeout_seconds: float | None = None


def _same_file(left: str | Path, right: str | Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def protected_collision(path: Path, protected_paths: Iterable[str | Path]) -> str | None:
    """Return the protected path aliased by a derived path.

    Args:
        path: Derived path whose lexical name and inode are inspected.
        protected_paths: Input paths that the derived path must not alias.

    Returns:
        The aliased protected path, or ``None`` when there is no collision.
    """
    path_absolute = os.path.abspath(path)
    for protected in protected_paths:
        protected_absolute = os.path.abspath(protected)
        if path_absolute == protected_absolute:
            return protected_absolute
        if os.path.lexists(path) and _same_file(path, protected):
            return protected_absolute
    return None


def validate_log_entry(path: Path, protected_paths: Iterable[str | Path]) -> None:
    """Reject a command-log entry that could mutate protected input data.

    Args:
        path: Final command-log path to validate.
        protected_paths: Input paths that the log must not alias.

    Raises:
        ValueError: If the log aliases protected data or has an unsafe type.
    """
    protected = protected_collision(path, protected_paths)
    if protected is not None:
        raise ValueError(f"Log path {path} aliases protected source {protected}")
    if os.path.lexists(path) and path.is_dir() and not path.is_symlink():
        raise ValueError(f"Log path has a wrong type and will not be replaced: {path}")


def capture_command(
    command: str,
    log_file: str,
    cwd: str | None = None,
    *,
    protected_paths: Iterable[str | Path] = (),
    timeout_seconds: float | None = None,
    metadata: CaptureMetadata | None = None,
) -> tuple[bool, str]:
    """Run a shell command while retaining its combined output for parsing.

    Args:
        command: Complete shell command to run under Bash.
        log_file: File receiving the command's combined stdout and stderr.
        cwd: Working directory for the child, or ``None`` to inherit it.
        protected_paths: Paths whose inodes the final log must not alias.
        timeout_seconds: Optional wall-clock deadline. A timeout kills and reaps
            the complete shell process group before this function returns.
        metadata: Optional mutable result receiving runner-authenticated facts.

    Returns:
        A success flag and complete output or safe diagnostic. Non-zero exits and
        OS failures are returned so preflight can try another candidate.
    """
    if metadata is not None:
        metadata.timed_out = False
        metadata.timeout_seconds = None
    logger.debug(f"Running captured command: {command}")
    final_log = Path(log_file)
    try:
        if not final_log.parent.is_dir():
            raise OSError(f"log directory does not exist: {final_log.parent}")
        validate_log_entry(final_log, protected_paths)
    except (OSError, ValueError) as error:
        diagnostic = f"Unable to write command log safely: {error}"
        logger.error(diagnostic)
        return False, diagnostic
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=final_log.parent,
            prefix=f".{final_log.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_log:
            temporary_path = Path(temporary_log.name)
            try:
                if timeout_seconds is None:
                    completed = subprocess.run(
                        command,
                        shell=True,
                        executable="/bin/bash",
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        cwd=cwd,
                        check=False,
                    )
                    output = completed.stdout or ""
                    return_code = completed.returncode
                else:
                    process = subprocess.Popen(
                        command,
                        shell=True,
                        executable="/bin/bash",
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        cwd=cwd,
                        start_new_session=True,
                    )
                    try:
                        output, _ = process.communicate(timeout=timeout_seconds)
                        output = output or ""
                        return_code = process.returncode
                    except subprocess.TimeoutExpired as error:
                        if metadata is not None:
                            metadata.timed_out = True
                            metadata.timeout_seconds = timeout_seconds
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except OSError:
                            with suppress(OSError):
                                process.kill()
                        reaped_output, _ = process.communicate()
                        partial_output = error.output or reaped_output or ""
                        if isinstance(partial_output, bytes):
                            partial_output = partial_output.decode(errors="replace")
                        timeout_text = f"Command timed out after {timeout_seconds:g} seconds."
                        output = f"{partial_output.rstrip()}\n{timeout_text}" if partial_output else timeout_text
                        return_code = 124
            except OSError as error:
                output = f"Unable to run command: {error}"
                return_code = 1
            temporary_log.write(output)
        os.replace(temporary_path, final_log)
        temporary_path = None
    except OSError as error:
        if metadata is not None:
            metadata.timed_out = False
            metadata.timeout_seconds = None
        diagnostic = f"Unable to write command log safely: {error}"
        logger.error(diagnostic)
        return False, diagnostic
    finally:
        if temporary_path is not None:
            with suppress(OSError):
                os.unlink(temporary_path)
    if return_code != 0:
        logger.error(f"Command failed: {command}")
        return False, output
    return True, output
