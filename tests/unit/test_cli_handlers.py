"""``handle_pipeline``: the seam between parsed arguments and ``run_pipeline``.

The handler is where CLI values become pipeline arguments, and two of them were being
lost on the way:

* ``--output-name`` was resolved from ``config.json`` and then **dropped** - the call to
  ``run_pipeline`` never mentions it. See :mod:`vntyper.scripts.artifact_names` for why
  it is now refused rather than threaded.
* ``--extra-modules`` was tested by exact string membership with no comma split and no
  validation, so ``--extra-modules advntr,shark`` matched neither name and produced a
  silent Kestrel-only run.

Both failed the same way: quietly, with a green exit code and a report that reads as a
negative genotype. These tests drive the real handler with ``run_pipeline`` stubbed and
assert on what it was asked to do.
"""

import inspect
import logging
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts import cli_handlers
from vntyper.scripts.artifact_names import PIPELINE_BASENAME
from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.pipeline import run_pipeline

pytestmark = pytest.mark.unit

MINIMAL_CONFIG = {
    "default_values": {
        "output_dir": "out",
        "threads": 4,
        "reference_assembly": "hg19",
        "output_name": PIPELINE_BASENAME,
        "archive_format": "zip",
    },
    "reference_data": {"bwa_reference_hg19": "/refs/hg19.fa", "bwa_reference_hg38": "/refs/hg38.fa"},
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


def test_the_pipeline_basename_is_accepted_explicitly(tmp_path: Path) -> None:
    """Asking for the basename the pipeline already uses is a no-op, not an error.

    Args:
        tmp_path: Pytest temporary directory.
    """
    stub = _run_handler(
        ["pipeline", "-o", str(tmp_path), "--bam", "in.bam", "--output-name", PIPELINE_BASENAME],
    )
    assert stub.call_count == 1


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
