"""Reference-file and candidate-policy tests for alignment preflight."""

from __future__ import annotations

import os
import shlex
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts import alignment_preflight
from vntyper.scripts.alignment_preflight import (
    HTSLIB_REFERENCE_SOURCE,
    resolve_reference,
    run_preflight,
)
from vntyper.scripts.preflight_error_io import PreflightErrorContext

pytestmark = pytest.mark.unit


def _reference(path: Path) -> Path:
    """Create a minimal readable reference FASTA."""
    path.write_text(">chr1\nACGT\n")
    return path


def test_a_nonexistent_explicit_fasta_is_not_probed_and_the_next_candidate_can_win(tmp_path: Path) -> None:
    """A mocked command success must not authorize a path that does not exist."""
    missing = tmp_path / "missing.fa"
    usable = _reference(tmp_path / "usable.fa")

    with patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "decoded")) as capture:
        reference, source, _, _binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(missing)), ("config_cram_reference", str(usable))),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
        )

    assert reference is not None
    assert Path(reference).read_bytes() == usable.read_bytes()
    assert source == "config_cram_reference"
    assert capture.call_count == 1
    assert _binding is not None
    _binding.close()


def test_an_unreadable_explicit_fasta_is_named_and_skipped_before_the_probe(tmp_path: Path) -> None:
    """An unreadable candidate cannot win merely because a command double says success."""
    unreadable = _reference(tmp_path / "unreadable.fa")
    usable = _reference(tmp_path / "usable.fa")
    real_open = os.open

    def deny_one_path(path: str | Path, flags: int):
        if Path(path) == unreadable:
            raise PermissionError("permission denied")
        return real_open(path, flags)

    with (
        patch("vntyper.scripts.preflight_input_io.os.open", side_effect=deny_one_path),
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "decoded")) as capture,
    ):
        reference, source, _, _binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(unreadable)), ("config_cram_reference", str(usable))),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
        )

    assert reference is not None
    assert Path(reference).read_bytes() == usable.read_bytes()
    assert source == "config_cram_reference"
    assert capture.call_count == 1
    assert _binding is not None
    _binding.close()


def test_missing_and_unreadable_fastas_are_named_in_the_final_diagnostic(tmp_path: Path) -> None:
    """Operator diagnostics distinguish local file failures from decode failures."""
    missing = tmp_path / "missing.fa"
    unreadable = _reference(tmp_path / "unreadable.fa")
    real_open = os.open

    def deny_one_path(path: str | Path, flags: int):
        if Path(path) == unreadable:
            raise PermissionError("permission denied")
        return real_open(path, flags)

    with (
        patch("vntyper.scripts.preflight_input_io.os.open", side_effect=deny_one_path),
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "ambient failed")) as capture,
        pytest.raises(ValueError) as error,
    ):
        resolve_reference(
            "/run/view.cram",
            (("cli", str(missing)), ("config_cram_reference", str(unreadable))),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
        )

    message = str(error.value)
    assert "reference FASTA not found" in message
    assert "reference FASTA unreadable" in message
    assert "permission denied" in message
    assert capture.call_count == 1


def test_a_missing_fai_warns_that_coverage_is_unavailable_but_does_not_reject_a_successful_probe(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The real decode probe remains authoritative when no FASTA index is available."""
    reference = _reference(tmp_path / "reference.fa")

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "decoded")),
        caplog.at_level("WARNING"),
    ):
        resolved, source, uncovered, _binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(reference)),),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "abc",
        )

    assert resolved is not None
    assert Path(resolved).read_bytes() == reference.read_bytes()
    assert (source, uncovered) == ("cli", ())
    assert "coverage unavailable" in caplog.text.lower()
    assert _binding is not None
    _binding.close()


def test_an_unreadable_fai_warns_that_coverage_is_unavailable_but_keeps_the_probe_winner(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A successfully decoded FASTA is not rejected solely because its FAI cannot be read."""
    reference = _reference(tmp_path / "reference.fa")
    fai = Path(f"{reference}.fai")
    fai.write_text("chr1\t4\t6\t4\t5\n")
    real_open = os.open

    def deny_fai(path: str | Path, flags: int):
        if Path(path) == fai:
            raise PermissionError("permission denied")
        return real_open(path, flags)

    with (
        patch("vntyper.scripts.preflight_input_io.os.open", side_effect=deny_fai),
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "decoded")),
        caplog.at_level("WARNING"),
    ):
        resolved, source, uncovered, _binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(reference)),),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "abc",
        )

    assert resolved is not None
    assert Path(resolved).read_bytes() == reference.read_bytes()
    assert (source, uncovered) == ("cli", ())
    assert "coverage unavailable" in caplog.text.lower()
    assert _binding is not None
    _binding.close()


@pytest.mark.parametrize(
    ("order", "message"),
    [
        (("cli", "htslib_resolved"), "must be a list"),
        (["cli", "cli", "htslib_resolved"], "duplicate"),
        (["none", "htslib_resolved"], "unknown reference candidate source: none"),
        (["header_ur", "htslib_resolved"], "unknown reference candidate source: header_ur"),
        (["htslib", "htslib_resolved"], "unknown reference candidate source: htslib"),
        (["htslib_resolved", "cli"], "terminal htslib_resolved"),
        (["cli"], "exactly one terminal htslib_resolved"),
        (["cli", "htslib_resolved", "htslib_resolved"], "exactly one terminal htslib_resolved"),
    ],
)
def test_invalid_reference_candidate_policies_fail_before_creating_the_view(
    tmp_path: Path, order: object, message: str
) -> None:
    """Configuration is validated verbatim rather than silently reordered or deduplicated."""
    alignment = tmp_path / "sample.cram"
    alignment.write_bytes(b"CRAM\x02")
    (tmp_path / "sample.cram.crai").write_bytes(b"CRAI")
    output = tmp_path / "output"
    config = {"cram": {"reference_candidate_order": order}}
    clean_idxstats = "chr1\t4\t1\t0\n*\t0\t0\t2\n"

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, clean_idxstats)) as capture,
        pytest.raises(ValueError, match=message),
    ):
        run_preflight(
            str(alignment),
            str(output),
            "sample",
            "cram",
            config,
            1,
            region="chr1:1-2",
        )

    assert not output.exists()
    capture.assert_not_called()


def test_explicit_candidates_are_probed_in_order_before_one_no_reference_probe(tmp_path: Path) -> None:
    """The policy produces explicit ``-T`` probes followed by exactly one ambient probe."""
    references = tuple(_reference(tmp_path / f"candidate-{position}.fa") for position in range(3))
    commands: list[str] = []

    def fail_until_ambient(command: str, log_file: str, cwd: str | None = None, **kwargs: object) -> tuple[bool, str]:
        del log_file, cwd, kwargs
        commands.append(command)
        return (" -T " not in f" {command} ", "did not decode")

    candidates = tuple(zip(("cli", "config_cram_reference", "config_bwa_reference"), map(str, references), strict=True))
    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=fail_until_ambient):
        reference, source, _, _binding = resolve_reference(
            "/run/view.cram",
            candidates,
            "chr1:1-2",
            None,
            {},
            4,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
        )

    assert len(commands) == 4
    assert all(" -T " in f" {command} " for command in commands[:3])
    assert " -T " not in f" {commands[-1]} "
    assert (reference, source) == (None, HTSLIB_REFERENCE_SOURCE)


def test_the_first_explicit_candidate_that_decodes_wins(tmp_path: Path) -> None:
    """Resolution stops immediately once the real decode probe succeeds."""
    first = _reference(tmp_path / "first.fa")
    second = _reference(tmp_path / "second.fa")
    outcomes = [(False, "wrong digest"), (True, "decoded")]

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=outcomes) as capture:
        reference, source, _, _binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(first)), ("config_cram_reference", str(second))),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            (),
            None,
        )

    assert reference is not None
    assert Path(reference).read_bytes() == second.read_bytes()
    assert source == "config_cram_reference"
    assert capture.call_count == 2
    assert _binding is not None
    _binding.close()


def test_a_total_decode_failure_names_every_candidate_and_reason(tmp_path: Path) -> None:
    """The terminal diagnostic retains explicit, omitted, and ambient attempts."""
    cli_reference = _reference(tmp_path / "cli.fa")
    bwa_reference = _reference(tmp_path / "bwa.fa")
    outcomes = [(False, "checksum mismatch"), (False, "missing chr2"), (False, "UR unavailable")]
    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=outcomes),
        pytest.raises(ValueError) as error,
    ):
        resolve_reference(
            "/run/view.cram",
            (
                ("cli", str(cli_reference)),
                ("config_cram_reference", None),
                ("config_bwa_reference", str(bwa_reference)),
            ),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "012345",
        )

    message = str(error.value)
    for expected in (
        "contig=chr1",
        "M5=012345",
        f"path={cli_reference}",
        "checksum mismatch",
        "not supplied",
        f"path={bwa_reference}",
        "missing chr2",
        HTSLIB_REFERENCE_SOURCE,
        "UR unavailable",
    ):
        assert expected in message


@pytest.mark.parametrize("target_source", ["region", "bed"])
def test_a_reference_error_names_the_actual_target_contig(tmp_path: Path, target_source: str) -> None:
    """Diagnostics use the selected region or BED target rather than the first header contig."""
    bed_file = tmp_path / "target.bed"
    bed_file.write_text("# target\n\nchr2\t100\t200\n")
    region = "chr1:100-200"

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "not found")),
        pytest.raises(ValueError) as error,
    ):
        resolve_reference(
            "/run/view.cram",
            (),
            region,
            bed_file if target_source == "bed" else None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chrM", "chr1", "chr2"),
            "digest",
        )

    expected_contig = "chr2" if target_source == "bed" else "chr1"
    assert f"contig={expected_contig}" in str(error.value)
    assert "contig=chrM" not in str(error.value)


def test_reference_probes_use_the_runs_region_or_bed_target(tmp_path: Path) -> None:
    """Probe commands preserve target selection and BED precedence."""
    reference = _reference(tmp_path / "full.fa")
    bed_file = tmp_path / "actual regions.bed"
    bed_file.write_text("chr9\t90\t99\n")
    commands: list[str] = []

    def succeed(command: str, log_file: str, cwd: str | None = None, **kwargs: object) -> tuple[bool, str]:
        del log_file, cwd, kwargs
        commands.append(command)
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=succeed):
        resolve_reference(
            "/run/view.cram",
            (("cli", str(reference)),),
            "chr9:90-99",
            None,
            {},
            1,
            str(tmp_path),
            "region",
            (),
            None,
        )
        resolve_reference(
            "/run/view.cram",
            (("cli", str(reference)),),
            "chr9:90-99",
            bed_file,
            {},
            1,
            str(tmp_path),
            "bed",
            (),
            None,
        )

    assert "chr9:90-99" in commands[0]
    assert f"-L '{bed_file}'" in commands[1]
    assert "chr9:90-99" not in commands[1]


def test_reference_winner_must_also_decode_the_independent_coverage_region(tmp_path: Path) -> None:
    """The slice proof cannot authorize a reference that the later depth consumer rejects."""
    target_only = _reference(tmp_path / "target-only.fa")
    full = _reference(tmp_path / "full.fa")
    target_only.write_text(">chr1\nAAAA\n", encoding="utf-8")
    full.write_text(">chr1\nCCCC\n", encoding="utf-8")
    commands: list[str] = []
    depth_reference_bytes: list[bytes] = []

    def decode(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        commands.append(command)
        if " depth " in f" {command} ":
            arguments = shlex.split(command)
            candidate = Path(arguments[arguments.index("--reference") + 1])
            candidate_bytes = candidate.read_bytes()
            depth_reference_bytes.append(candidate_bytes)
            if candidate_bytes == b">chr1\nAAAA\n":
                return False, "Unable to fetch reference chr2"
        return True, "decoded"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=decode):
        reference, source, _, _binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(target_only)), ("config_cram_reference", str(full))),
            "chr1:1-2",
            None,
            {},
            2,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "target-checksum",
            coverage_region="chr2:30-40",
        )

    assert reference is not None
    assert Path(reference).read_bytes() == b">chr1\nCCCC\n"
    assert source == "config_cram_reference"
    assert len(commands) == 4
    assert all(" -P " in f" {commands[position]} " for position in (0, 2))
    for position in (1, 3):
        assert " depth -a " in f" {commands[position]} "
        assert "-r chr2:30-40" in commands[position]
    assert depth_reference_bytes == [b">chr1\nAAAA\n", b">chr1\nCCCC\n"]
    assert _binding is not None
    _binding.close()


def test_known_fai_differences_are_returned_as_uncovered_contigs(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Known coverage gaps remain actionable even though the decode probe succeeds."""
    reference = _reference(tmp_path / "chr1.fa")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n")

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "decoded")),
        caplog.at_level("WARNING"),
    ):
        _, _, uncovered, _binding = resolve_reference(
            "/run/view.cram",
            (("config_bwa_reference", str(reference)),),
            "chr1:1-2",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "abc",
        )

    assert uncovered == ("chr2",)
    assert "chr2" in caplog.text
    assert _binding is not None
    _binding.close()


def test_stream_mode_rejects_a_target_only_winner_and_tries_the_next_candidate(tmp_path: Path) -> None:
    """A candidate must decode both the target and the later whole-file stream."""
    target_only = _reference(tmp_path / "target-only.fa")
    full = _reference(tmp_path / "full.fa")
    target_only.write_text(">chr1\nAAAA\n", encoding="utf-8")
    full.write_text(">chr1\nCCCC\n", encoding="utf-8")
    commands: list[str] = []

    def decode(command: str, log_file: str, cwd: str | None = None, **kwargs: object) -> tuple[bool, str]:
        del log_file, cwd, kwargs
        commands.append(command)
        if " -P " in f" {command} ":
            return True, "target decoded"
        arguments = shlex.split(command)
        candidate = Path(arguments[arguments.index("-T") + 1])
        if candidate.read_bytes() == b">chr1\nAAAA\n":
            return False, "Unable to fetch reference chr2"
        return True, "whole file decoded"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=decode):
        reference, source, _, _binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(target_only)), ("config_cram_reference", str(full))),
            "chr1:1-2",
            None,
            {"cram": {"reference_probe_timeout_seconds": 7}},
            1,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "abc",
            unmapped_scan="stream",
        )

    assert reference is not None
    assert Path(reference).read_bytes() == b">chr1\nCCCC\n"
    assert source == "config_cram_reference"
    assert len(commands) == 4
    assert " -P " in f" {commands[0]} "
    assert " -P " not in f" {commands[1]} "
    assert " -P " in f" {commands[2]} "
    assert " -P " not in f" {commands[3]} "
    assert _binding is not None
    _binding.close()


def test_terminal_stream_failure_names_the_later_missing_contig_and_its_m5(tmp_path: Path) -> None:
    """A whole-file failure must not be attributed to the first otherwise-decodable BED target."""
    bed_file = tmp_path / "targets.bed"
    bed_file.write_text("chr1\t1\t100\nchr2\t1\t100\n", encoding="utf-8")
    outcomes = [
        (True, "target decoded"),
        (False, "Reference file given, but ref 'chr2' not present"),
    ]
    failure_context = PreflightErrorContext(tmp_path)

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=outcomes),
        pytest.raises(ValueError) as error,
    ):
        resolve_reference(
            "/run/view.cram",
            (),
            None,
            bed_file,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "first-checksum",
            header_m5s=(("chr1", "first-checksum"), ("chr2", "missing-checksum")),
            unmapped_scan="stream",
            failure_context=failure_context,
        )

    assert "contig=chr2" in str(error.value)
    assert "M5=missing-checksum" in str(error.value)
    assert failure_context.payload is not None
    assert "contig=chr2" in failure_context.payload["message"]
    assert "M5=missing-checksum" in failure_context.payload["message"]


@pytest.mark.parametrize("configured", [0, -1, 120.1, "30", True, float("inf"), float("nan")])
def test_reference_probe_timeout_rejects_values_outside_the_bounded_numeric_contract(configured: object) -> None:
    """A replacement config cannot disable or exceed the 120-second deadline."""
    with pytest.raises(ValueError, match="reference_probe_timeout_seconds"):
        alignment_preflight._reference_probe_timeout_seconds({"cram": {"reference_probe_timeout_seconds": configured}})


def test_reference_probe_timeout_defaults_to_120_and_accepts_a_smaller_deadline() -> None:
    """The shipped ceiling is also the replacement-config default."""
    assert alignment_preflight._reference_probe_timeout_seconds({}) == 120.0
    assert (
        alignment_preflight._reference_probe_timeout_seconds({"cram": {"reference_probe_timeout_seconds": 0.25}})
        == 0.25
    )
