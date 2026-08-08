"""What has to be true for a CRAM upload to be genotyped as a CRAM (#188).

`docker/app/uploads.py` has accepted `.cram` for the alignment slot since the
allowlist was introduced (`ALIGNMENT_EXTENSIONS = ("bam", "cram")`), but the
worker handed every accepted alignment to the CLI's `--bam` flag. The CLI has
had `--cram` since `cli_parser.py`, `cli_handlers.handle_pipeline` threads it
into `run_pipeline`, and `pipeline.py` branches on it -- so only the task layer
was wrong, and every accepted CRAM took the BAM code path.

Correcting the flag on its own is not enough. Three things move together, and
two of them are only *introduced* by the correction:

1. The flag itself. `--cram` for a CRAM, `--bam` for a BAM, chosen without
   regard to case because the upload allowlist matches case-insensitively and
   `SAMPLE.CRAM` is an accepted name.
2. The sample name. `handle_pipeline` derived it from `--bam` or `--fastq1` and
   never from `--cram`. Today a CRAM arrives through `--bam` and gets its file
   stem; the moment the worker starts sending `--cram`, every CRAM run without
   an explicit `--sample-name` would be named the literal string `"sample"` --
   in the report and in the output filenames. The one-line fix introduces that
   regression, so the two must land together.
3. The generated index name. `samtools index` writes `.crai` beside a CRAM and
   `.bai` beside a BAM. The worker's fallback named `.bai` for both, so for a
   CRAM submitted without an index it named a file that is never created:
   cleanup then removes nothing and leaves the real index on the volume every
   job shares.

The command used to be assembled inline inside `run_vntyper_job`, which made
the flag choice unreachable without a Celery worker; it now lives in
`tasks.build_vntyper_command`, and the index fallback in
`tasks.resolve_index_path`, so both are asserted here directly.

The sample-name assertion drives the real `handle_pipeline` with `run_pipeline`
stubbed. It is a `vntyper/` behaviour rather than a web one and would normally
sit in `tests/unit/test_cli_handlers.py`; it is here because #188 is one defect
with three symptoms and reads better in one place than split across two files.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.tasks` is importable here.
"""

import logging
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

pytestmark = pytest.mark.unit

from app.tasks import build_vntyper_command, resolve_index_path  # noqa: E402
from app.uploads import safe_upload_path  # noqa: E402

from vntyper.scripts import cli_handlers  # noqa: E402
from vntyper.scripts.artifact_names import PIPELINE_BASENAME  # noqa: E402
from vntyper.scripts.cli_parser import build_parser  # noqa: E402
from vntyper.scripts.utils import validate_bam_file  # noqa: E402

MINIMAL_CONFIG: dict[str, dict[str, object]] = {
    "default_values": {
        "threads": 4,
        "reference_assembly": "hg19",
        "output_name": PIPELINE_BASENAME,
        "archive_format": "zip",
    },
    "reference_data": {"bwa_reference_hg19": "/refs/hg19.fa", "bwa_reference_hg38": "/refs/hg38.fa"},
}


def _pipeline_kwargs(argv: list[str]) -> dict[str, Any]:
    """Parse ``argv`` and return the kwargs ``handle_pipeline`` passes on.

    Args:
        argv: Arguments after ``vntyper``, starting with ``pipeline``.

    Returns:
        dict[str, Any]: The keyword arguments ``run_pipeline`` was called with.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    with mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as stub:
        cli_handlers.handle_pipeline(
            args,
            config=MINIMAL_CONFIG,
            parser=parser,
            log_level_value=logging.INFO,
            log_file_str=None,
        )
    return dict(stub.call_args.kwargs)


# ---------------------------------------------------------------------------
# 1. The flag the worker hands the CLI.
# ---------------------------------------------------------------------------


def test_a_cram_upload_is_passed_to_the_cli_as_cram_not_bam() -> None:
    """The web service accepts .cram; it must not hand it to the BAM flag.

    Issue #188. The CLI has had a ``--cram`` flag since ``cli_parser.py``,
    ``cli_handlers.handle_pipeline`` threads it into ``run_pipeline``, and
    ``pipeline.py`` branches on it. Only the task layer hardcoded ``--bam``, so
    every accepted CRAM took the BAM code path.
    """
    command = build_vntyper_command(
        alignment_path="/data/sample.cram", output_dir="/out", thread=4, reference_assembly="hg38"
    )

    assert "--cram" in command
    assert "--bam" not in command
    assert command[command.index("--cram") + 1] == "/data/sample.cram"


def test_a_bam_upload_still_uses_the_bam_flag() -> None:
    """The correction must not send a BAM down the CRAM path instead."""
    command = build_vntyper_command(
        alignment_path="/data/sample.bam", output_dir="/out", thread=4, reference_assembly="hg38"
    )

    assert "--bam" in command
    assert "--cram" not in command
    assert command[command.index("--bam") + 1] == "/data/sample.bam"


def test_the_flag_is_chosen_case_insensitively() -> None:
    """``uploads.py`` accepts ``SAMPLE.CRAM``, so the flag choice must match on case.

    The allowlist in ``docker/app/uploads.py`` compiles its pattern with
    ``re.IGNORECASE``, so an upper-case extension is a name the endpoint stores
    and enqueues. A case-sensitive ``endswith(".cram")`` would send exactly
    those submissions back down the BAM path.
    """
    command = build_vntyper_command(
        alignment_path="/data/SAMPLE.CRAM", output_dir="/out", thread=4, reference_assembly="hg38"
    )

    assert "--cram" in command
    assert "--bam" not in command


def test_every_other_flag_is_unchanged_by_the_extraction() -> None:
    """The extraction must be behaviour-preserving apart from the flag choice.

    Pins the whole argument vector rather than a handful of members: the point
    of lifting the command out of ``run_vntyper_job`` is that nothing else about
    it moved, and only an exact comparison can say so.
    """
    command = build_vntyper_command(
        alignment_path="/data/s.bam",
        output_dir="/out",
        thread=8,
        reference_assembly="hg19",
        fast_mode=True,
        keep_intermediates=True,
        archive_results=True,
        advntr_mode=True,
    )

    assert command == [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "vntyper",
        "vntyper",
        "pipeline",
        "--bam",
        "/data/s.bam",
        "-o",
        "/out",
        "--threads",
        "8",
        "--reference-assembly",
        "hg19",
        "--fast-mode",
        "--keep-intermediates",
        "--archive-results",
        "--extra-modules",
        "advntr",
        "--advntr-max-coverage",
        "300",
    ]


def test_no_optional_flag_is_appended_when_none_was_asked_for() -> None:
    """The defaults produce the base command and nothing more."""
    command = build_vntyper_command(
        alignment_path="/data/s.cram", output_dir="/out", thread=2, reference_assembly="hg19"
    )

    assert command == [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "vntyper",
        "vntyper",
        "pipeline",
        "--cram",
        "/data/s.cram",
        "-o",
        "/out",
        "--threads",
        "2",
        "--reference-assembly",
        "hg19",
    ]


@pytest.mark.parametrize(
    ("filename", "expected_flag"),
    [
        ("SAMPLE.CRAM", "--cram"),
        ("SAMPLE.BAM", "--bam"),
        ("sample.Cram", "--cram"),
        ("sample.cram", "--cram"),
        ("sample.bam", "--bam"),
    ],
    ids=["upper-cram", "upper-bam", "mixed-cram", "lower-cram", "lower-bam"],
)
def test_a_name_the_endpoint_accepts_survives_all_three_layers(
    tmp_path: Path, filename: str, expected_flag: str
) -> None:
    """The three layers that see the same filename must agree about it.

    The flag assertions above stop at the argument vector the worker builds:
    they never run that vector's alignment path through the validator the CLI
    reaches for first. That left a real gap open. ``uploads.py`` accepted
    ``SAMPLE.CRAM``, ``tasks.py`` routed it to ``--cram``, and
    ``vntyper.scripts.utils.validate_bam_file`` then rejected it with
    ``ValueError: Invalid alignment file extension`` because it compared the
    suffix case-sensitively -- so the job was accepted, enqueued, and died in
    validation. Every existing test passed throughout.

    This walks one name through all three in order, which is the only shape of
    test that can catch a disagreement between them.

    Args:
        tmp_path: Pytest temporary directory, standing in for the job input dir.
        filename: A client-supplied upload name.
        expected_flag: The CLI flag that name must be routed to.
    """
    stored_path = safe_upload_path(str(tmp_path), filename)
    Path(stored_path).touch()

    command = build_vntyper_command(
        alignment_path=stored_path, output_dir=str(tmp_path), thread=1, reference_assembly="hg38"
    )
    assert expected_flag in command
    handed_to_the_cli = command[command.index(expected_flag) + 1]

    with mock.patch("vntyper.scripts.utils.run_command", return_value=True):
        validate_bam_file(handed_to_the_cli)  # must not raise


# ---------------------------------------------------------------------------
# 2. The conventional index name worker cleanup falls back to.
# ---------------------------------------------------------------------------


def test_a_cram_upload_with_no_index_falls_back_to_the_crai_name() -> None:
    """CRAM cleanup uses the conventional ``.crai`` fallback.

    Pipeline preflight builds missing indexes only in the output tree. This
    fallback keeps cleanup compatible with conventional CRAM index residue.
    """
    assert resolve_index_path("/data/s.cram", None) == "/data/s.cram.crai"


def test_a_bam_upload_with_no_index_still_falls_back_to_the_bai_name() -> None:
    """The BAM fallback is unchanged; without this the test above passes on a rename."""
    assert resolve_index_path("/data/s.bam", None) == "/data/s.bam.bai"


def test_the_crai_fallback_is_also_chosen_case_insensitively() -> None:
    """``SAMPLE.CRAM`` is an accepted upload name and gets a ``.crai`` too."""
    assert resolve_index_path("/data/SAMPLE.CRAM", None) == "/data/SAMPLE.CRAM.crai"


@pytest.mark.parametrize(
    ("alignment", "supplied"),
    [
        ("/data/s.cram", "/data/given.crai"),
        ("/data/s.bam", "/data/given.bai"),
    ],
)
def test_an_index_the_submission_carried_is_used_unchanged(alignment: str, supplied: str) -> None:
    """The fallback only applies when the submission carried no index.

    Args:
        alignment: The stored alignment path.
        supplied: The index path the endpoint handed the worker.
    """
    assert resolve_index_path(alignment, supplied) == supplied


# ---------------------------------------------------------------------------
# 3. The sample name a CRAM run reports under.
# ---------------------------------------------------------------------------


def test_a_cram_run_derives_its_sample_name_from_the_file_stem(tmp_path: Path) -> None:
    """``handle_pipeline`` derived the name from ``--bam`` or ``--fastq1`` only.

    Switching the worker to ``--cram`` without this would name every CRAM run
    ``"sample"`` -- the fallback at the end of the same chain. Today a CRAM
    arrives via ``--bam`` and gets its stem, so the one-line #188 fix would
    INTRODUCE this regression. The two must land together.

    Args:
        tmp_path: Pytest temporary directory, used as the output directory.
    """
    kwargs = _pipeline_kwargs(["pipeline", "-o", str(tmp_path), "--cram", "/data/NA12878.cram"])

    assert kwargs["cram"] == "/data/NA12878.cram"
    assert kwargs["sample_name"] == "NA12878"


def test_a_bam_run_still_derives_its_sample_name_from_the_file_stem(tmp_path: Path) -> None:
    """The CRAM arm must be added beside the BAM one, not in place of it.

    Args:
        tmp_path: Pytest temporary directory, used as the output directory.
    """
    kwargs = _pipeline_kwargs(["pipeline", "-o", str(tmp_path), "--bam", "/data/NA12878.bam"])

    assert kwargs["sample_name"] == "NA12878"


def test_an_explicit_sample_name_still_wins_over_the_cram_stem(tmp_path: Path) -> None:
    """The derivation is a fallback; a name the caller gave is not overwritten.

    Args:
        tmp_path: Pytest temporary directory, used as the output directory.
    """
    kwargs = _pipeline_kwargs(
        ["pipeline", "-o", str(tmp_path), "--cram", "/data/NA12878.cram", "--sample-name", "given-name"]
    )

    assert kwargs["sample_name"] == "given-name"
