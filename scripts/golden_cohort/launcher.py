"""The wrapper that proves which tree a gate run actually executed.

This is the part of the gate that could otherwise make the whole exercise worthless, and
the gate page says so. The ``vntyper`` console script resolves the package through
setuptools' editable finder, which is **appended** to ``sys.meta_path`` and points at
whichever worktree the editable install was made from - irrespective of the current
directory. A baseline checkout therefore silently executes candidate code: with the
process CWD set to the baseline worktree, a plain ``import vntyper`` still resolved to the
candidate worktree and reported ``vntyper.scripts.pipeline_guards`` present.

Every run of every side is launched through :func:`launch`, which

1. changes into its own tree and puts that tree first on ``sys.path`` (a path entry, so
   ``PathFinder`` - which precedes the appended editable finder - resolves it),
2. prints its resolved ``vntyper.__file__`` and marker state **as its first line**,
3. ``sys.exit``s before doing any work unless both agree with its side.

The tree also has to be the working directory, not merely importable: ``vntyper/config.json``
names ``reference/...``, ``vntyper/templates`` and ``vntyper/dependencies/kestrel/kestrel.jar``
by **relative** path, so a run whose CWD is elsewhere reads another tree's references,
templates and Kestrel jar even when the code is right.

The marker is a parameter, not a constant: it differs per comparison. Runs 1-3 used the
module ``vntyper.scripts.pipeline_guards``, which exists only after ``078a6c4``.

A branch that only *modifies* files has no such module to name, and naming one both sides
have distinguishes nothing -- so ``--marker`` also accepts ``module:attribute``. That form
is the stronger witness of the two, because it is reached through the same import machinery
as the code under test rather than through a sibling module; the cost is that observing an
attribute requires importing its module, where ``find_spec`` executes nothing.

While the run is under way every ``subprocess.Popen`` construction is appended to a JSONL
file, which is what run 3 added and what let it show that ``2ae28c5``'s ``shlex.quote``
changes no executed command in this cohort.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shlex
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any

#: Printed as the first line of every launched run.
LAUNCH_PREFIX = "GATE-LAUNCH"

#: Printed instead when the resolution does not agree with the side, before anything runs.
ABORT_PREFIX = "GATE-ABORT"

#: Exit code used when the wrapper refuses to start. Distinct from the CLI's own 0 and 1
#: so an abort can never be mistaken for a pipeline result.
EXIT_ABORT = 97

#: Exit code used when the CLI raised something it did not handle itself.
EXIT_UNHANDLED = 98


def _record_commands(log_path: Path) -> Any:
    """Append every command the run executes to ``log_path`` as JSONL.

    ``subprocess.run``, ``check_output`` and ``call`` all construct a ``Popen``, so
    wrapping its ``__init__`` records all of them and nothing twice.

    Args:
        log_path: The JSONL file to append to.

    Returns:
        Any: The original ``subprocess.Popen.__init__`` method, for restoration after
        the wrapped CLI returns.
    """
    original_init = subprocess.Popen.__init__
    log_path.parent.mkdir(parents=True, exist_ok=True)

    def recording_init(self: Any, args: Any, *rest: Any, **kwargs: Any) -> None:
        rendered = shlex.join(str(part) for part in args) if isinstance(args, (list, tuple)) else str(args)
        try:
            with open(log_path, "a", encoding="utf-8") as handle:
                handle.write(json.dumps({"command": rendered, "shell": bool(kwargs.get("shell"))}) + "\n")
        except OSError:  # pragma: no cover - recording must never break the run
            pass
        original_init(self, args, *rest, **kwargs)

    subprocess.Popen.__init__ = recording_init  # type: ignore[method-assign]
    return original_init


def parse_marker(marker: str) -> tuple[str, str | None]:
    """Split a marker into a module name and, if the marker names one, an attribute.

    ``module`` and ``module:attribute`` are both accepted, the second spelled as an entry
    point is. A half-written marker is refused rather than degraded: ``mod:`` checked as a
    bare ``mod`` would answer *present* on both sides of a comparison whose whole purpose
    is that the two sides differ, and the gate would then compare two sides it believed it
    had verified.

    Args:
        marker: The ``--marker`` value.

    Returns:
        tuple[str, str | None]: The module name, and the attribute or None.

    Raises:
        ValueError: If either half of a ``module:attribute`` marker is empty.
    """
    if marker.count(":") > 1:
        raise ValueError(f"marker {marker!r} has more than one ':'; write module:attribute or a bare module")
    module, separator, attribute = marker.partition(":")
    if not separator:
        if not module:
            raise ValueError("marker is empty")
        return module, None
    if not module or not attribute:
        raise ValueError(f"marker {marker!r} has an empty half; write module:attribute or a bare module")
    return module, attribute


def marker_is_present(marker: str) -> bool:
    """Report whether ``marker`` resolves in the current interpreter.

    A bare module is answered with ``find_spec``, which locates without executing -- runs
    1-3 relied on that, because the check happens before ``vntyper.cli.main`` configures
    logging. An attribute cannot be observed without importing its module, which is the
    price of a marker for a branch that adds no module of its own.

    Args:
        marker: A module name or ``module:attribute``.

    Returns:
        bool: Whether it resolves.
    """
    module, attribute = parse_marker(marker)
    if attribute is None:
        return importlib.util.find_spec(module) is not None
    return hasattr(importlib.import_module(module), attribute)


def resolve(tree: Path, marker: str) -> dict[str, Any]:
    """Import ``vntyper`` from ``tree`` and report what resolved.

    Args:
        tree: The source tree this side must run from.
        marker: A module name, or ``module:attribute``, whose presence distinguishes the two
            sides -- e.g. ``vntyper.scripts.pipeline_guards`` or
            ``vntyper.modules.advntr.advntr_genotyping:resolve_advntr_threads``.

    Returns:
        dict[str, Any]: ``vntyper_file``, ``in_tree``, ``marker``, ``marker_present`` and
        ``error``. ``in_tree`` is False rather than raising when the import resolved
        somewhere else, because that is the failure this whole wrapper exists to catch.
    """
    info: dict[str, Any] = {
        "tree": str(tree),
        "vntyper_file": None,
        "in_tree": False,
        "marker": marker,
        "marker_present": None,
        "error": None,
    }
    try:
        import vntyper
    except Exception as exc:  # noqa: BLE001 - any import failure is a refusal to start
        info["error"] = f"{type(exc).__name__}: {exc}"
        return info

    resolved = Path(vntyper.__file__).resolve()
    info["vntyper_file"] = str(resolved)
    info["in_tree"] = tree in resolved.parents

    try:
        module, attribute = parse_marker(marker)
    except ValueError as exc:
        info["marker_present"] = False
        info["error"] = f"marker({marker}): ValueError: {exc}"
        return info

    try:
        info["marker_present"] = marker_is_present(marker)
    except Exception as exc:  # noqa: BLE001 - any resolution failure is a refusal to start
        # A missing *parent* package raises rather than returning None, and importing a
        # module to reach an attribute can fail in any way that module's import can.
        info["marker_present"] = False
        label = f"find_spec({marker})" if attribute is None else f"import_module({module})"
        info["error"] = f"{label}: {type(exc).__name__}: {exc}"
    return info


def _launch_line(side: str, info: dict[str, Any], expect_marker: bool) -> str:
    """Render the one line every run prints before it does anything.

    Args:
        side: ``before`` or ``after``.
        info: The mapping from :func:`resolve`.
        expect_marker: Whether the marker must be present on this side.

    Returns:
        str: The line.
    """
    state = {True: "present", False: "absent", None: "unknown"}[info["marker_present"]]
    return (
        f"{LAUNCH_PREFIX} side={side} tree={info['tree']} vntyper_file={info['vntyper_file']} "
        f"in_tree={info['in_tree']} marker={info['marker']} marker_state={state} "
        f"expected_marker={'present' if expect_marker else 'absent'}"
    )


def launch(
    *,
    tree: Path,
    side: str,
    marker: str,
    expect_marker: bool,
    commands_log: Path | None,
    argv: list[str],
) -> int:
    """Run one ``vntyper`` invocation, but only from the tree it claims to be running.

    Args:
        tree: The source tree this side must run from.
        side: ``before`` or ``after``, for the printed line.
        marker: The marker, as a module name or ``module:attribute``.
        expect_marker: Whether the marker must be present on this side.
        commands_log: Where to append the executed commands, or None to not record them.
        argv: The ``vntyper`` argument list, without the program name.

    Returns:
        int: The process exit code - the CLI's own, :data:`EXIT_ABORT` if the resolution
        disagreed with the side, or :data:`EXIT_UNHANDLED` if the CLI raised.
    """
    original_cwd = Path.cwd()
    original_path = sys.path.copy()
    original_argv = sys.argv.copy()
    original_popen_init: Any | None = None
    tree = tree.resolve()
    try:
        os.chdir(tree)
        sys.path.insert(0, str(tree))

        info = resolve(tree, marker)
        print(_launch_line(side, info, expect_marker), flush=True)

        if info["error"] is not None and info["vntyper_file"] is None:
            print(f"{ABORT_PREFIX} reason=import-failed detail={info['error']}", flush=True)
            return EXIT_ABORT
        if not info["in_tree"]:
            print(
                f"{ABORT_PREFIX} reason=wrong-tree resolved={info['vntyper_file']} expected-under={tree}",
                flush=True,
            )
            return EXIT_ABORT
        if info["error"] is not None:
            # `resolve` records a marker that would not resolve, and until this check
            # existed nothing consumed the record: an unimportable marker module printed
            # `marker_state=absent`, indistinguishable from a genuinely absent attribute.
            # On the side expecting absent that voided the witness while the launch line
            # asserted it held.
            print(f"{ABORT_PREFIX} reason=marker-unresolvable detail={info['error']}", flush=True)
            return EXIT_ABORT
        if info["marker_present"] is not expect_marker:
            print(
                f"{ABORT_PREFIX} reason=marker-mismatch marker={marker} "
                f"got={info['marker_present']} expected={expect_marker}",
                flush=True,
            )
            return EXIT_ABORT

        if commands_log is not None:
            original_popen_init = _record_commands(commands_log)

        from vntyper.cli import main  # noqa: PLC0415 - deliberately after the resolution check

        sys.argv[:] = ["vntyper", *argv]
        try:
            main()
        except SystemExit as exc:
            code = exc.code
            if code is None:
                return 0
            return code if isinstance(code, int) else 1
        except Exception:  # noqa: BLE001 - the harness records failures, it does not raise them
            traceback.print_exc()
            return EXIT_UNHANDLED
        return 0
    finally:
        if original_popen_init is not None:
            subprocess.Popen.__init__ = original_popen_init  # type: ignore[method-assign]
        sys.argv[:] = original_argv
        sys.path[:] = original_path
        os.chdir(original_cwd)
