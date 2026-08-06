"""
Unit tests for ``vntyper/scripts/output_paths.py``: report destinations stay inside ``--output-dir``.

``--report-file`` and ``--summary-file`` are documented as *names* - "Name of the output report
file" - but were typed as plain strings and joined onto ``--output-dir`` with no check at all.
``--report-file ../../summary.html`` and ``--report-file /etc/cron.d/vntyper`` both wrote outside the
directory the user nominated, over whatever was already there, silently.

The rule enforced here is deliberately narrower than "does not escape": a single path segment,
inside the output directory, and **rejected rather than repaired**. Quietly rewriting
``../../report.html`` to ``report.html`` would write a file the caller did not ask for and say
nothing, which is the same class of surprise in a friendlier costume.
"""

import os
from pathlib import Path

import pandas as pd
import pytest

import vntyper
from vntyper.cli import load_config
from vntyper.scripts import cohort_summary
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.output_paths import contained_output_path

pytestmark = pytest.mark.unit

TEMPLATE_DIR = Path(vntyper.__file__).resolve().parent / "templates"

#: Values that are not a single path segment inside the output directory.
ESCAPES = [
    ("../report.html", "a parent-relative name climbs out of the output directory"),
    ("../../etc/passwd", "two levels up, and onto a real file"),
    ("sub/report.html", "a subdirectory is still not the directory the user named"),
    ("./report.html", "leading `./` is a path, not a name"),
    ("..", "the parent directory itself"),
    ("", "no name at all"),
    ("   ", "whitespace is not a filename"),
]


@pytest.mark.parametrize("name,reason", ESCAPES, ids=[repr(name) for name, _ in ESCAPES])
def test_a_value_that_is_not_a_single_name_is_refused(tmp_path, name, reason) -> None:
    """Every rejected value raises rather than being quietly turned into something acceptable."""
    with pytest.raises(ValueError, match="single file name"):
        contained_output_path(tmp_path, name, "--report-file")

    assert list(tmp_path.iterdir()) == [], reason


def test_an_absolute_path_is_refused(tmp_path) -> None:
    """`--report-file /etc/cron.d/vntyper` ignored `--output-dir` entirely."""
    with pytest.raises(ValueError, match="single file name"):
        contained_output_path(tmp_path, "/etc/cron.d/vntyper", "--report-file")


def test_the_option_name_is_in_the_message(tmp_path) -> None:
    """Two options share this check, so the error has to say which one was wrong."""
    with pytest.raises(ValueError, match="--summary-file"):
        contained_output_path(tmp_path, "../out.html", "--summary-file")


def test_an_ordinary_name_is_joined_onto_the_output_directory(tmp_path) -> None:
    """The everyday case is unchanged: the name is used exactly as given."""
    assert contained_output_path(tmp_path, "summary_report.html", "--report-file") == tmp_path / "summary_report.html"


def test_a_name_with_spaces_and_dots_is_accepted(tmp_path) -> None:
    """Only path *structure* is refused; the check is not a filename allowlist."""
    result = contained_output_path(tmp_path, "cohort 2024.v2.report.html", "--summary-file")

    assert result == tmp_path / "cohort 2024.v2.report.html"


def test_a_symlinked_output_directory_is_still_its_own_directory(tmp_path) -> None:
    """
    Containment is decided after resolution, so a symlinked output dir is not a false refusal.

    `--output-dir` is the caller's own choice; the check is on what the *name* does to it.
    """
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real)

    result = contained_output_path(link, "report.html", "--report-file")

    assert os.path.realpath(result.parent) == os.path.realpath(real)


def test_a_name_that_is_a_symlink_out_of_the_directory_is_refused(tmp_path) -> None:
    """
    The last containment check is on the resolved path, so a planted symlink does not slip past.

    A single segment satisfies the shape test but can still resolve elsewhere, and writing through
    it overwrites the target.
    """
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "victim.html").write_text("original", encoding="utf-8")
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "report.html").symlink_to(outside / "victim.html")

    with pytest.raises(ValueError, match="output directory"):
        contained_output_path(output_dir, "report.html", "--report-file")

    assert (outside / "victim.html").read_text(encoding="utf-8") == "original"


def test_the_output_directory_may_be_a_string(tmp_path) -> None:
    """`--output-dir` is `type=str` on both subcommands, so both forms have to work."""
    assert contained_output_path(str(tmp_path), "report.html", "--report-file") == Path(tmp_path) / "report.html"


# ---------------------------------------------------------------------------
# The two writers that take one
# ---------------------------------------------------------------------------


def test_the_per_sample_report_refuses_to_write_outside_its_output_directory(tmp_path) -> None:
    """`vntyper report --report-file ../escaped.html` used to overwrite a sibling run's report."""
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    victim = tmp_path / "escaped.html"
    victim.write_text("someone else's report", encoding="utf-8")

    with pytest.raises(ValueError, match="--report-file"):
        generate_summary_report(
            output_dir=str(output_dir),
            template_dir=str(TEMPLATE_DIR),
            report_file="../escaped.html",
            log_file=None,
            config=load_config(None),
        )

    assert victim.read_text(encoding="utf-8") == "someone else's report"


def test_the_cohort_report_refuses_to_write_outside_its_output_directory(tmp_path) -> None:
    """`vntyper cohort --summary-file` has the identical hole and gets the identical check."""
    output_dir = tmp_path / "cohort"
    output_dir.mkdir()
    victim = tmp_path / "escaped.html"
    victim.write_text("someone else's report", encoding="utf-8")

    with pytest.raises(ValueError, match="--summary-file"):
        cohort_summary.generate_cohort_summary_report(
            output_dir=str(output_dir),
            kestrel_df=pd.DataFrame(),
            advntr_df=pd.DataFrame(),
            summary_file="../escaped.html",
            config=load_config(None),
        )

    assert victim.read_text(encoding="utf-8") == "someone else's report"
