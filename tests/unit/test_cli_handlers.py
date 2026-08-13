"""``handle_pipeline``: the seam between parsed arguments and ``run_pipeline``.

The handler is where CLI values become pipeline arguments, and two of them are pinned
here because losing one costs nothing visible:

* ``--output-name`` is refused rather than threaded. ``run_pipeline`` takes no such
  parameter; see :mod:`vntyper.scripts.artifact_names` for why honouring it is not a
  one-line change.
* ``--extra-modules`` is comma-split and validated against the known module names, so
  ``--extra-modules advntr,shark`` selects both.

A value lost between the two would fail quietly: a green exit code and a report that
reads as a negative genotype. These tests drive the real handler with ``run_pipeline``
stubbed and assert on what it was asked to do.
"""

import ast
import inspect
import logging
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts import cli_handlers, pipeline
from vntyper.scripts.artifact_names import PIPELINE_BASENAME
from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.pipeline import run_pipeline

pytestmark = pytest.mark.unit

#: `_resolve_bwa_reference` now fails closed on a configured-but-missing file (Important
#: 1), so a fake, never-written path like the old "/refs/hg19.fa" would make every
#: FASTQ-input test below raise before `run_pipeline` was ever reached. Pointing at this
#: test module's own file is a real, always-present path that carries no other meaning -
#: nothing here asserts on the literal value, only on whether resolution succeeded.
_EXISTING_FILE = str(Path(__file__).resolve())

MINIMAL_CONFIG: dict[str, dict[str, object]] = {
    "default_values": {
        "output_dir": "out",
        "threads": 4,
        "reference_assembly": "hg19",
        "output_name": PIPELINE_BASENAME,
        "archive_format": "zip",
    },
    "reference_data": {"bwa_reference_hg19": _EXISTING_FILE, "bwa_reference_hg38": _EXISTING_FILE},
}


def _run_handler(argv: list[str], config: dict | None = None) -> mock.MagicMock:
    """Parse ``argv`` and run ``handle_pipeline`` with ``run_pipeline`` stubbed.

    Args:
        argv: Arguments after ``vntyper``, starting with ``pipeline``.
        config: Configuration mapping. Defaults to :data:`MINIMAL_CONFIG`.

    Returns:
        unittest.mock.MagicMock: The ``run_pipeline`` stub.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    with mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as stub:
        cli_handlers.handle_pipeline(
            args,
            config=config if config is not None else MINIMAL_CONFIG,
            parser=parser,
            log_level_value=logging.INFO,
            log_file_str=None,
        )
    return stub


# --------------------------------------------------------------------------------------
# --output-name (defect D2)
# --------------------------------------------------------------------------------------


def test_the_default_output_name_runs_the_pipeline(tmp_path: Path) -> None:
    """An ordinary run is unaffected: no flag, no error, pipeline reached.

    Args:
        tmp_path: Pytest temporary directory.
    """
    stub = _run_handler(["pipeline", "-o", str(tmp_path), "--bam", "in.bam"])
    assert stub.call_count == 1


def test_reference_fasta_is_forwarded_to_the_pipeline(tmp_path: Path) -> None:
    """The explicit CRAM reference must survive the parser/handler boundary.

    Args:
        tmp_path: Pytest temporary directory.
    """
    reference = tmp_path / "full reference.fa"

    stub = _run_handler(
        ["pipeline", "-o", str(tmp_path / "out"), "--cram", "in.cram", "--reference-fasta", str(reference)]
    )

    assert stub.call_args.kwargs["reference_fasta"] == reference


def test_the_pipeline_basename_is_accepted_explicitly(tmp_path: Path) -> None:
    """Asking for the basename the pipeline already uses is a no-op, not an error.

    Args:
        tmp_path: Pytest temporary directory.
    """
    stub = _run_handler(
        ["pipeline", "-o", str(tmp_path), "--bam", "in.bam", "--output-name", PIPELINE_BASENAME],
    )
    assert stub.call_count == 1


def test_a_single_fastq_is_forwarded_without_a_paired_end_usage_error(tmp_path: Path) -> None:
    stub = _run_handler(["pipeline", "-o", str(tmp_path), "--fastq1", "single.fastq.gz"])

    assert stub.call_count == 1
    assert stub.call_args.kwargs["fastq1"] == "single.fastq.gz"
    assert stub.call_args.kwargs["fastq2"] is None


def test_a_single_fastq_with_shark_is_a_usage_error_before_pipeline_dispatch(tmp_path: Path) -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "pipeline",
            "-o",
            str(tmp_path),
            "--fastq1",
            "single.fastq.gz",
            "--extra-modules",
            "shark",
        ]
    )

    with mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as stub, pytest.raises(SystemExit) as excinfo:
        cli_handlers.handle_pipeline(
            args,
            config=MINIMAL_CONFIG,
            parser=parser,
            log_level_value=logging.INFO,
            log_file_str=None,
        )

    assert excinfo.value.code == 2
    stub.assert_not_called()


def test_an_output_name_the_pipeline_cannot_honour_is_refused(tmp_path: Path, caplog) -> None:
    """The defect: ``-n mysample`` used to be resolved, dropped and forgotten.

    Args:
        tmp_path: Pytest temporary directory.
        caplog: Pytest log capture.
    """
    with caplog.at_level(logging.ERROR), pytest.raises(ValueError) as excinfo:
        _run_handler(["pipeline", "-o", str(tmp_path), "--bam", "in.bam", "-n", "mysample"])

    assert "mysample" in str(excinfo.value)
    assert PIPELINE_BASENAME in str(excinfo.value)
    assert any(record.levelno >= logging.ERROR and "mysample" in record.message for record in caplog.records)


def test_a_configured_output_name_the_pipeline_cannot_honour_is_refused(tmp_path: Path) -> None:
    """A ``config.json`` that disagrees with the pipeline must fail loudly too.

    ``default_values.output_name`` was ``"processed"`` for the whole life of the flag,
    silently, because the value never reached ``run_pipeline``. Threading it would have
    renamed every artefact on an ordinary run; refusing it says so instead.

    Args:
        tmp_path: Pytest temporary directory.
    """
    config = {**MINIMAL_CONFIG, "default_values": {**MINIMAL_CONFIG["default_values"], "output_name": "processed"}}
    with pytest.raises(ValueError, match="processed"):
        _run_handler(["pipeline", "-o", str(tmp_path), "--bam", "in.bam"], config=config)


def _extra_modules_passed(argv: list[str]) -> list[str]:
    """Run the handler and return the ``extra_modules`` ``run_pipeline`` received.

    Args:
        argv: Arguments after ``vntyper``, starting with ``pipeline``.

    Returns:
        list[str]: The forwarded module list.
    """
    stub = _run_handler(argv)
    assert stub.call_count == 1, "the handler never reached run_pipeline"
    return stub.call_args.kwargs["extra_modules"]


# --------------------------------------------------------------------------------------
# --extra-modules (defect D3)
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("flags", "expected"),
    [
        pytest.param([], [], id="none"),
        pytest.param(["--extra-modules", "advntr"], ["advntr"], id="single"),
        pytest.param(["--extra-modules", "advntr,shark"], ["advntr", "shark"], id="comma_separated"),
        pytest.param(
            ["--extra-modules", "advntr", "--extra-modules", "shark"],
            ["advntr", "shark"],
            id="repeated_flag",
        ),
        pytest.param(["--extra-modules", " advntr , shark "], ["advntr", "shark"], id="padded"),
        pytest.param(["--extra-modules", "advntr,,shark,"], ["advntr", "shark"], id="empty_fields"),
        pytest.param(["--extra-modules", "AdVNTR,SHARK"], ["advntr", "shark"], id="mixed_case"),
        pytest.param(
            ["--extra-modules", "advntr,shark", "--extra-modules", "advntr"],
            ["advntr", "shark"],
            id="duplicates_collapsed",
        ),
        pytest.param(["--extra-modules", "shark,advntr"], ["shark", "advntr"], id="order_preserved"),
    ],
)
def test_extra_modules_are_normalised(flags: list[str], expected: list[str], tmp_path: Path) -> None:
    """``--extra-modules advntr,shark`` used to match neither name.

    Membership was tested by exact string against the raw ``append`` list, so a
    comma-separated value produced a silent Kestrel-only run with exit code 0.

    Args:
        flags: The ``--extra-modules`` arguments under test.
        expected: The list ``run_pipeline`` must receive.
        tmp_path: Pytest temporary directory.
    """
    argv = ["pipeline", "-o", str(tmp_path), "--fastq1", "r1.fq", "--fastq2", "r2.fq", *flags]
    assert _extra_modules_passed(argv) == expected


@pytest.mark.parametrize("name", ["advtnr", "sharkk", "kestrel", "adVNTR2"])
def test_an_unknown_extra_module_is_refused(name: str, tmp_path: Path, caplog) -> None:
    """A typo used to disable the module silently; now it names what it knows.

    Args:
        name: The misspelled module name.
        tmp_path: Pytest temporary directory.
        caplog: Pytest log capture.
    """
    argv = ["pipeline", "-o", str(tmp_path), "--fastq1", "r1.fq", "--fastq2", "r2.fq", "--extra-modules", name]
    with caplog.at_level(logging.ERROR), pytest.raises(ValueError) as excinfo:
        _run_handler(argv)

    assert name in str(excinfo.value)
    for known in cli_handlers.KNOWN_EXTRA_MODULES:
        assert known in str(excinfo.value), "the error must name the modules that do exist"
    assert any(record.levelno >= logging.ERROR and name in record.message for record in caplog.records)


def test_an_unknown_module_beside_a_known_one_is_still_refused(tmp_path: Path) -> None:
    """One good name must not excuse a bad one in the same comma list.

    Args:
        tmp_path: Pytest temporary directory.
    """
    argv = [
        "pipeline",
        "-o",
        str(tmp_path),
        "--fastq1",
        "r1.fq",
        "--fastq2",
        "r2.fq",
        "--extra-modules",
        "advntr,shrak",
    ]
    with pytest.raises(ValueError, match="shrak"):
        _run_handler(argv)


def test_advntr_module_args_survive_the_comma_form(tmp_path: Path) -> None:
    """The adVNTR module arguments hang off the same membership test.

    ``--advntr-max-coverage`` was silently ignored whenever ``--extra-modules`` was
    comma-separated, because ``"advntr" in ["advntr,shark"]`` is False.

    Args:
        tmp_path: Pytest temporary directory.
    """
    argv = [
        "pipeline",
        "-o",
        str(tmp_path),
        "--fastq1",
        "r1.fq",
        "--fastq2",
        "r2.fq",
        "--extra-modules",
        "advntr,shark",
        "--advntr-max-coverage",
        "300",
    ]
    stub = _run_handler(argv)
    assert stub.call_args.kwargs["module_args"]["advntr"]["max_coverage"] == 300


def test_shark_with_bam_input_is_rejected_in_the_comma_form(tmp_path: Path) -> None:
    """SHARK is FASTQ-only (#62), and the comma form used to slip past that guard.

    Args:
        tmp_path: Pytest temporary directory.
    """
    argv = ["pipeline", "-o", str(tmp_path), "--bam", "in.bam", "--extra-modules", "advntr,shark"]
    with pytest.raises(SystemExit) as excinfo:
        _run_handler(argv)
    assert excinfo.value.code == 1


def test_the_known_modules_are_the_ones_the_pipeline_tests_for() -> None:
    """The CLI's whitelist and ``pipeline.py``'s membership tests must agree.

    A module the CLI accepts but the pipeline never looks for is a silent no-op;
    one the pipeline looks for but the CLI rejects is unreachable.
    """
    tree = ast.parse(inspect.getsource(pipeline))
    tested = {
        node.left.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Compare)
        and isinstance(node.left, ast.Constant)
        and isinstance(node.left.value, str)
        and any(isinstance(op, ast.In) for op in node.ops)
        and any(isinstance(cmp, ast.Name) and cmp.id == "extra_modules" for cmp in node.comparators)
    }
    assert tested, "found no `'x' in extra_modules` tests in pipeline.py; this test would be vacuous"
    assert tested == set(cli_handlers.KNOWN_EXTRA_MODULES)


def test_run_pipeline_still_takes_no_output_name_parameter() -> None:
    """Pins the reason the flag is refused rather than threaded.

    If ``run_pipeline`` ever grows an ``output_name`` parameter, this fails - and that
    is the moment to revisit :func:`vntyper.scripts.artifact_names.validate_output_name`,
    but only once ``generate_report.py``, ``cli_report.py`` and ``kestrel_genotyping.py``
    take the basename too.
    """
    parameters = inspect.signature(run_pipeline).parameters
    assert parameters, "run_pipeline has no parameters; this assertion would be vacuous"
    assert "output_name" not in parameters


# --------------------------------------------------------------------------------------
# Where the sample name came from (#242)
# --------------------------------------------------------------------------------------


def test_an_operators_sample_name_is_recorded_as_explicit(tmp_path: Path) -> None:
    """``--sample-name`` is a name, and the run says so.

    Without the flag the summary carries a bare string, and the report cannot tell
    ``--sample-name sample`` from the ``"sample"`` the CLI falls back to - which is
    how a valid explicit name came to be discarded.

    Args:
        tmp_path: Pytest temporary directory.
    """
    stub = _run_handler(["pipeline", "-o", str(tmp_path), "--bam", "in.bam", "-s", "PATIENT_042"])

    assert stub.call_args.kwargs["sample_name"] == "PATIENT_042"
    assert stub.call_args.kwargs["sample_name_is_explicit"] is True


@pytest.mark.parametrize(
    "argv_tail,expected_name",
    [
        (["--bam", "in.bam"], "in"),
        (["--cram", "in.cram"], "in"),
        (["--fastq1", "S1_R1.fastq.gz"], "S1_R1.fastq"),
    ],
)
def test_a_name_derived_from_an_input_path_is_recorded_as_derived(
    argv_tail: list[str], expected_name: str, tmp_path: Path
) -> None:
    """The three reachable derivations, and the value each one actually records.

    The FASTQ case is the one that matters: ``Path("S1_R1.fastq.gz").stem`` is
    ``S1_R1.fastq``, a half-stripped file name. It stays that here on purpose -
    ``run_kestrel`` builds its output filenames from this same string - and the flag
    is what lets the report finish deriving it into ``S1`` without renaming anything.

    Args:
        argv_tail: The input flags for this shape.
        expected_name: The name the handler records for it.
        tmp_path: Pytest temporary directory.
    """
    stub = _run_handler(["pipeline", "-o", str(tmp_path), *argv_tail])

    assert stub.call_args.kwargs["sample_name"] == expected_name
    assert stub.call_args.kwargs["sample_name_is_explicit"] is False
