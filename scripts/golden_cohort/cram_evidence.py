"""Fail-closed A-178-2 decisions over recorded CRAM read-set evidence."""

from __future__ import annotations

import json
import re
import shlex
from pathlib import Path
from typing import Any

MAX_GUARD_COUNT_DIGITS = 20
MAX_COMMAND_LOG_BYTES = 1024 * 1024
GUARD_MESSAGE_PATTERN = (
    rf"idxstats reports ([0-9]{{1,{MAX_GUARD_COUNT_DIGITS}}}) placed-unmapped reads; using stream scan"
)
PLACED_UNMAPPED_GUARD_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(GUARD_MESSAGE_PATTERN),
    re.compile(rf"ValueError: {GUARD_MESSAGE_PATTERN}"),
    re.compile(
        rf"[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}},[0-9]{{3}} - "
        rf"vntyper\.scripts\.idxstats_parsing - ERROR - {GUARD_MESSAGE_PATTERN}"
    ),
)


def observe_unmapped_scan(commands_log: Path, output_bam: Path) -> tuple[str | None, str | None, list[str]]:
    """Read the bounded command record and identify the executed CRAM scan.

    Args:
        commands_log: JSONL written by the launcher's ``Popen`` recorder.
        output_bam: Exact unmapped BAM destination owned by the conversion.

    Returns:
        Observed mode, its exact command, and fail-closed diagnostic problems.
    """
    try:
        with commands_log.open("rb") as handle:
            raw = handle.read(MAX_COMMAND_LOG_BYTES + 1)
    except OSError:
        return None, None, ["A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"]
    if len(raw) > MAX_COMMAND_LOG_BYTES:
        return None, None, [f"A-178-2 command log exceeds {MAX_COMMAND_LOG_BYTES} bytes"]
    try:
        lines = raw.decode("utf-8").splitlines()
    except UnicodeError:
        return None, None, ["A-178-2 command log is not valid UTF-8"]

    target = str(output_bam)
    observed: list[tuple[str, str]] = []
    target_referenced = False
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        try:
            entry = json.loads(line)
        except ValueError:
            return None, None, [f"A-178-2 command log is malformed at line {line_number}"]
        if (
            not isinstance(entry, dict)
            or set(entry) != {"command", "shell"}
            or not isinstance(entry["command"], str)
            or type(entry["shell"]) is not bool
        ):
            return None, None, [f"A-178-2 command log is malformed at line {line_number}"]
        command = entry["command"]
        target_referenced = target_referenced or target in command
        if not entry["shell"]:
            continue
        try:
            tokens = shlex.split(command)
        except ValueError:
            return None, None, [f"A-178-2 command log is malformed at line {line_number}"]
        writes_target = any(tokens[index : index + 2] == ["-o", target] for index in range(len(tokens) - 1))
        filters_unmapped = any(tokens[index : index + 2] == ["-f", "4"] for index in range(len(tokens) - 1))
        if not writes_target or not filters_unmapped:
            continue
        # What distinguishes the two modes is the `'*'` index query, not the shape of
        # the shell command. `'*'` fetches only the *unplaced* unmapped reads through
        # the index; its absence means a whole-file decode. Keying `stream` on the
        # presence of a shell pipe was incidental to how the stream mode happened to be
        # implemented, and stopped classifying it the moment that became one process.
        #
        # Widening `stream` to "anything that is not indexed" is only safe because the
        # `writes_target` and `filters_unmapped` guards above have already excluded
        # every command that is not an extraction into this exact file.
        if "*" in tokens:
            observed.append(("indexed", command))
        else:
            observed.append(("stream", command))

    modes = sorted({mode for mode, _command in observed})
    if not modes:
        if target_referenced:
            return (
                None,
                None,
                ["A-178-2 command log references the unmapped BAM but no recognized extraction mode was observed"],
            )
        return None, None, ["A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"]
    if len(modes) > 1:
        return None, None, [f"A-178-2 observed contradictory CRAM unmapped-extraction modes: {', '.join(modes)}"]
    commands = list(dict.fromkeys(command for _mode, command in observed))
    if len(commands) != 1:
        return None, None, [f"A-178-2 observed multiple distinct {modes[0]} unmapped-extraction commands"]
    return modes[0], commands[0], []


def parse_placed_unmapped_guard_count(stderr: str) -> int | None:
    """Return the unique placed-unmapped count from the stable guard diagnostic.

    Args:
        stderr: Complete pipeline stderr, including propagated logger output.

    Returns:
        The diagnostic count when every exact occurrence agrees, otherwise ``None``.
    """
    counts: set[int] = set()
    for line in stderr.splitlines():
        for pattern in PLACED_UNMAPPED_GUARD_PATTERNS:
            match = pattern.fullmatch(line)
            if match is not None:
                counts.add(int(match.group(1)))
                break
    if len(counts) != 1:
        return None
    return counts.pop()


def validate_cram_evidence(case: dict[str, Any], record: dict[str, Any]) -> list[str]:
    """Validate one candidate CRAM result against its declared evidence policy.

    Args:
        case: Runtime matrix case with an optional ``cram_evidence_expectation``.
        record: Result record containing raw-indexed and produced-unmapped evidence.

    Returns:
        Human-readable problems. An empty list is the only passing decision.
    """
    expectation = case.get("cram_evidence_expectation")
    if expectation is None:
        return []

    problems: list[str] = []
    expected_raw = expectation.get("raw_indexed_read_set")
    actual_raw = record.get("raw_indexed_read_set")
    if actual_raw != expected_raw:
        problems.append(f"A-178-2 raw indexed evidence differs: expected {expected_raw}, got {actual_raw}")

    observation_problems = record.get("scan_observation_problems")
    if not isinstance(observation_problems, list) or not all(
        isinstance(problem, str) for problem in observation_problems
    ):
        observation_problems = []
    scan = record.get("observed_unmapped_scan")
    command = record.get("observed_unmapped_command")
    actual_stream = record.get("unmapped_read_set")
    expected_stream = expectation.get("stream_read_set")
    if expectation.get("indexed_authorized") is True:
        if scan not in {"indexed", "stream"} or not isinstance(command, str) or not command:
            problems.extend(
                observation_problems or ["A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"]
            )
            return problems
        configured_scan = case.get("effective_unmapped_scan")
        if scan != configured_scan:
            problems.append(f"A-178-2 observed {scan} extraction, but declared {configured_scan}")
            return problems
        if actual_stream != expected_stream:
            problems.append(f"A-178-2 stream evidence differs: expected {expected_stream}, got {actual_stream}")
        actual_loss = record.get("raw_indexed_loss")
        if actual_loss is None:
            problems.append("A-178-2 authorized evidence did not record the raw indexed loss")
        elif isinstance(expected_raw, dict) and isinstance(expected_stream, dict):
            raw_count = expected_raw.get("count")
            stream_count = expected_stream.get("count")
            if isinstance(raw_count, int) and isinstance(stream_count, int):
                expected_loss = stream_count - raw_count
                if actual_loss != expected_loss:
                    problems.append(f"A-178-2 raw indexed loss differs: expected {expected_loss}, got {actual_loss}")
        return problems

    configured_scan = case.get("effective_unmapped_scan")
    if configured_scan == "indexed":
        expected_guard_count = expectation.get("placed_unmapped_guard_count")
        actual_guard_count = record.get("placed_unmapped_guard_count")
        nonzero_observation_problems = [
            problem
            for problem in observation_problems
            if problem != "A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"
        ]
        problems.extend(nonzero_observation_problems)
        if record.get("exit_code") != 1:
            problems.append(f"A-178-2 indexed rejection exit differs: expected 1, got {record.get('exit_code')}")
        if scan is not None:
            problems.append(
                f"A-178-2 indexed rejection observed a {scan} unmapped-extraction command instead of failing before work"
            )
        if command is not None:
            problems.append(
                "A-178-2 indexed rejection executed an unmapped-extraction command instead of failing before work"
            )
        if type(expected_guard_count) is not int or expected_guard_count <= 0:
            problems.append(f"A-178-2 indexed guard expectation is missing or malformed: {expected_guard_count}")
        elif type(actual_guard_count) is not int or actual_guard_count != expected_guard_count:
            problems.append(
                f"A-178-2 indexed guard count differs: expected {expected_guard_count}, got {actual_guard_count}"
            )
        if actual_stream is not None:
            problems.append("A-178-2 indexed rejection produced an unmapped BAM instead of failing before work")
        return problems

    if scan != "stream" or not isinstance(command, str) or not command:
        problems.extend(
            observation_problems or ["A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"]
        )
        return problems
    if configured_scan in {"indexed", "stream"} and scan != configured_scan:
        problems.append(f"A-178-2 observed {scan} extraction, but declared {configured_scan}")
        return problems

    if actual_stream != expected_stream:
        problems.append(f"A-178-2 stream evidence differs: expected {expected_stream}, got {actual_stream}")
    actual_loss = record.get("raw_indexed_loss")
    if actual_loss is None:
        problems.append("A-178-2 stream evidence did not record the raw indexed loss")
    elif isinstance(expected_raw, dict) and isinstance(expected_stream, dict):
        raw_count = expected_raw.get("count")
        stream_count = expected_stream.get("count")
        if isinstance(raw_count, int) and isinstance(stream_count, int):
            expected_loss = stream_count - raw_count
            if actual_loss != expected_loss:
                problems.append(f"A-178-2 raw indexed loss differs: expected {expected_loss}, got {actual_loss}")
    return problems
