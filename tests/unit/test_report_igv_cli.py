"""``--report-igv`` reaches the generator down **both** routes, and nothing else stops.

Two entry paths, and only one of them is obvious
------------------------------------------------
``vntyper report`` is the one an option like this looks like it belongs to. It is not
the one most reports come from: an ordinary run calls ``handle_pipeline`` ->
``run_pipeline`` -> ``generate_summary_report`` and never touches the ``report``
subcommand at all. An option wired only into ``cli_report.py`` would therefore be
accepted by the parser, documented, tested, and dead for every report the pipeline
produces - the failure mode that makes this worth a file of its own rather than one
more assertion in ``test_cli_report.py``.

So both routes are driven end to end here, through ``cli.main`` and through a spy that
**binds the real signature** of the function it stands in for. A spy taking
``*args, **kwargs`` records a broken call as happily as a correct one; one that runs
``inspect.signature(...).bind(...)`` fails exactly where the real call would.

The compatibility matrix
------------------------
Adding a keyword argument to a call site is how the *other* keyword arguments get
dropped: a rebased edit, a merge, a hand-applied patch. Every option the report handler
is responsible for forwarding is therefore asserted to arrive, in the same test that
asserts the new one does - ``--report-file``, ``--bed-file``, ``--bam-file``,
``--vcf-file``, ``--reference-fasta``, ``--flanking`` and ``--sample-name``, the last
of which an earlier task added to the same call.
"""

from __future__ import annotations

import ast
import inspect
import logging
from pathlib import Path
from typing import Any

import pytest

from vntyper import cli
from vntyper.scripts import cli_handlers, cli_report, generate_report, pipeline, report_assets
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.pipeline import run_pipeline

pytestmark = pytest.mark.unit


@pytest.fixture(autouse=True)
def _no_logging_reconfiguration(monkeypatch):
    """Keep ``setup_logging`` from tearing down pytest's log capture."""
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)


def _binding_spy(signature: inspect.Signature, calls: list[inspect.BoundArguments]):
    """Return a spy that binds ``signature`` and records the result.

    Args:
        signature: The signature the real call site must satisfy.
        calls: Where each bound call is recorded.

    Returns:
        Callable: The spy.
    """

    def _spy(*args: Any, **kwargs: Any) -> None:
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        calls.append(bound)

    return _spy


@pytest.fixture
def report_calls(monkeypatch) -> list[inspect.BoundArguments]:
    """Spy on ``generate_summary_report`` as ``cli_report`` calls it."""
    calls: list[inspect.BoundArguments] = []
    monkeypatch.setattr(
        cli_report, "generate_summary_report", _binding_spy(inspect.signature(generate_summary_report), calls)
    )
    return calls


@pytest.fixture
def pipeline_calls(monkeypatch) -> list[inspect.BoundArguments]:
    """Spy on ``run_pipeline`` as ``cli_handlers`` calls it.

    The pipeline route has two joints, and both can drop an argument independently:
    ``handle_pipeline`` -> ``run_pipeline``, and ``run_pipeline`` ->
    ``generate_summary_report``. This covers the first;
    :func:`test_run_pipeline_forwards_the_mode_to_the_generator` covers the second.
    """
    calls: list[inspect.BoundArguments] = []
    monkeypatch.setattr(cli_handlers, "run_pipeline", _binding_spy(inspect.signature(run_pipeline), calls))
    return calls


# ---------------------------------------------------------------------------
# Route 1: vntyper report
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("mode", report_assets.REPORT_IGV_MODES)
def test_the_report_subcommand_forwards_every_mode(tmp_path: Path, report_calls, mode: str) -> None:
    """Each of the three, through argparse and the handler."""
    cli.main(["report", "-o", str(tmp_path), "--report-igv", mode])

    assert len(report_calls) == 1, "the report handler did not reach generate_summary_report"
    assert report_calls[0].arguments["report_igv"] == mode


def test_the_report_subcommand_defaults_to_embedding(tmp_path: Path, report_calls) -> None:
    """An operator who asks for nothing gets the self-contained archive."""
    cli.main(["report", "-o", str(tmp_path)])

    assert report_calls[0].arguments["report_igv"] == report_assets.DEFAULT_REPORT_IGV
    assert report_calls[0].arguments["report_igv"] == "embedded"


def test_the_report_subcommand_rejects_an_unknown_mode(tmp_path: Path, report_calls) -> None:
    """A typo is a usage error at the parser, not a silent fallback in the generator.

    Exit code 2 is argparse's, and it is the contract this repository documents for
    argument errors (AGENTS.md, "Exit codes").
    """
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["report", "-o", str(tmp_path), "--report-igv", "embeded"])

    assert exit_info.value.code == 2
    assert report_calls == []


def test_every_report_option_still_reaches_the_generator(tmp_path: Path, report_calls) -> None:
    """The compatibility matrix for ``vntyper report``.

    Every option this handler forwards, in one call, so that adding a keyword argument
    to the call site cannot quietly drop another. ``--sample-name`` is in the list
    because an earlier task added it to this same call and it has the same exposure.
    """
    bed = tmp_path / "regions.bed"
    bam = tmp_path / "reads.bam"
    vcf = tmp_path / "calls.vcf"
    fasta = tmp_path / "ref.fa"
    for path in (bed, bam, vcf, fasta):
        path.write_text("", encoding="utf-8")

    cli.main(
        [
            "report",
            "-o",
            str(tmp_path),
            "--report-file",
            "custom_report.html",
            "--bed-file",
            str(bed),
            "--bam-file",
            str(bam),
            "--vcf-file",
            str(vcf),
            "--reference-fasta",
            str(fasta),
            "--flanking",
            "123",
            "--sample-name",
            "NA12878",
            "--report-igv",
            "sidecar",
        ]
    )

    assert len(report_calls) == 1
    arguments = report_calls[0].arguments
    delivered = {
        "report_file": arguments["report_file"],
        "bed_file": Path(arguments["bed_file"]),
        "bam_file": Path(arguments["bam_file"]),
        "vcf_file": Path(arguments["vcf_file"]),
        "fasta_file": Path(arguments["fasta_file"]),
        "flanking": arguments["flanking"],
        "sample_name": arguments["sample_name"],
        "report_igv": arguments["report_igv"],
    }

    assert delivered == {
        "report_file": "custom_report.html",
        "bed_file": bed,
        "bam_file": bam,
        "vcf_file": vcf,
        "fasta_file": fasta,
        "flanking": 123,
        "sample_name": "NA12878",
        "report_igv": "sidecar",
    }


# ---------------------------------------------------------------------------
# Route 2: vntyper pipeline, which is where most reports come from
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("mode", report_assets.REPORT_IGV_MODES)
def test_the_pipeline_subcommand_forwards_every_mode(tmp_path: Path, pipeline_calls, mode: str) -> None:
    """The route an ordinary run takes, which never touches ``vntyper report``."""
    cli.main(["pipeline", "-o", str(tmp_path), "--bam", "in.bam", "--report-igv", mode])

    assert len(pipeline_calls) == 1, "the pipeline handler did not reach run_pipeline"
    assert pipeline_calls[0].arguments["report_igv"] == mode


def test_the_pipeline_subcommand_defaults_to_embedding(tmp_path: Path, pipeline_calls) -> None:
    """And the default is the same one on both routes.

    Two defaults would be worse than one wrong one: the same run, described the same
    way, would produce different artifacts depending on which subcommand rendered it.
    """
    cli.main(["pipeline", "-o", str(tmp_path), "--bam", "in.bam"])

    assert pipeline_calls[0].arguments["report_igv"] == report_assets.DEFAULT_REPORT_IGV


def test_the_pipeline_subcommand_rejects_an_unknown_mode(tmp_path: Path, pipeline_calls) -> None:
    """Same usage error, same exit code, on the other subcommand."""
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["pipeline", "-o", str(tmp_path), "--bam", "in.bam", "--report-igv", "none"])

    assert exit_info.value.code == 2
    assert pipeline_calls == []


def test_every_pipeline_option_the_report_depends_on_still_reaches_run_pipeline(tmp_path: Path, pipeline_calls) -> None:
    """The compatibility matrix for ``vntyper pipeline``'s call to ``run_pipeline``.

    Narrower than the report handler's, because ``run_pipeline`` takes twenty-odd
    arguments and this is not the place to restate the whole pipeline contract. What is
    checked is the set that shares a call site with the new keyword and that a report
    depends on: the input, the assembly, the sample name and its provenance flag, the
    BED, and the log the report embeds.
    """
    bed = tmp_path / "regions.bed"
    bed.write_text("", encoding="utf-8")

    cli.main(
        [
            "--log-file",
            str(tmp_path / "run.log"),
            "pipeline",
            "-o",
            str(tmp_path),
            "--bam",
            "in.bam",
            "--bed-file",
            str(bed),
            "--reference-assembly",
            "hg38",
            "--sample-name",
            "NA12878",
            "--report-igv",
            "off",
        ]
    )

    arguments = pipeline_calls[0].arguments
    delivered = {
        "bam": arguments["bam"],
        "bed_file": Path(arguments["bed_file"]),
        "reference_assembly": arguments["reference_assembly"],
        "sample_name": arguments["sample_name"],
        "sample_name_is_explicit": arguments["sample_name_is_explicit"],
        "log_file": arguments["log_file"],
        "report_igv": arguments["report_igv"],
    }

    assert delivered == {
        "bam": "in.bam",
        "bed_file": bed,
        "reference_assembly": "hg38",
        "sample_name": "NA12878",
        "sample_name_is_explicit": True,
        "log_file": str(tmp_path / "run.log"),
        "report_igv": "off",
    }


def test_run_pipeline_forwards_the_mode_to_the_generator() -> None:
    """The second joint on the pipeline route, read out of the syntax tree.

    ``run_pipeline`` is 700 lines of orchestration over external binaries and cannot be
    driven in the unit tier, so this reads its call site instead - through ``ast``
    rather than by slicing the text, because the call already contains a nested
    ``config.get("default_values", {})`` and a text slice ends at the wrong bracket.

    The vacuity guard is the point: the assertion is worthless unless the call really is
    there, so the call itself is asserted first.
    """
    tree = ast.parse(Path(pipeline.__file__).read_text(encoding="utf-8"))
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "generate_summary_report"
    ]

    assert len(calls) == 1, f"expected one generate_summary_report call in pipeline.py, found {len(calls)}"

    forwarded = {keyword.arg: ast.unparse(keyword.value) for keyword in calls[0].keywords if keyword.arg is not None}

    assert forwarded.get("report_igv") == "report_igv", (
        "run_pipeline accepts --report-igv and does not pass it on, so the option is accepted, documented "
        f"and dead for every report a pipeline run produces. The call forwards: {sorted(forwarded)}"
    )


@pytest.mark.parametrize("function", [run_pipeline, generate_summary_report])
def test_both_signatures_default_to_the_same_mode(function) -> None:
    """A caller that passes nothing gets the same artifact wherever it called from."""
    default = inspect.signature(function).parameters["report_igv"].default

    assert default == report_assets.DEFAULT_REPORT_IGV


#: The minimum each subcommand needs to parse at all, so the check below is about
#: ``--report-igv`` and not about a missing required option.
_MINIMAL_ARGV: dict[str, list[str]] = {
    "pipeline": ["pipeline", "-o", "out", "--bam", "in.bam"],
    "report": ["report", "-o", "out"],
    "cohort": ["cohort", "-o", "out", "-i", "run"],
    "online": ["online"],
}


@pytest.mark.parametrize("mode", report_assets.REPORT_IGV_MODES)
def test_every_offered_mode_is_one_the_generator_accepts(mode: str) -> None:
    """Two spellings of one list is how a mode is accepted by the parser and refused
    by the generator, three frames deeper, in the middle of a run.

    Asserted behaviourally rather than by reading argparse's private ``_actions``:
    what matters is that a mode the parser takes is a mode the asset layer takes, and
    both halves are exercised here rather than compared as strings.
    """
    parser = cli.build_parser()

    for command in ("pipeline", "report"):
        parsed = parser.parse_args([*_MINIMAL_ARGV[command], "--report-igv", mode])
        assert parsed.report_igv == mode

    # And it is a mode the generator will not raise on.
    generate_report.report_assets.igv_payload(mode)


@pytest.mark.parametrize("command", ["cohort", "online"])
def test_the_option_is_declared_only_where_a_report_is_produced(command: str) -> None:
    """``--report-igv`` belongs to the two subcommands that render a per-sample report.

    The cohort report has no alignment browser and ``online`` submits work to a server,
    so an option accepted there would be one that silently does nothing.
    """
    parser = cli.build_parser()

    with pytest.raises(SystemExit) as exit_info:
        parser.parse_args([*_MINIMAL_ARGV[command], "--report-igv", "embedded"])

    assert exit_info.value.code == 2
