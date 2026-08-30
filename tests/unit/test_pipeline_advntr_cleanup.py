"""Exact planning and ownership rules for adVNTR preflight cleanup."""

from __future__ import annotations

from pathlib import Path

import pytest

from vntyper.scripts.pipeline_advntr_cleanup import (
    plan_advntr_cleanup,
    validate_pipeline_log_outside_advntr_preflight,
)

pytestmark = pytest.mark.unit


PUBLIC_OUTPUTS = (
    Path("advntr/output_adVNTR_result.tsv"),
    Path("advntr/output_adVNTR.tsv"),
    Path("advntr/output_adVNTR.vcf"),
    Path("advntr/cross_match_results.tsv"),
    Path("advntr/output_advntr.log"),
    Path("pipeline_summary.json"),
    Path("pipeline_summary.csv"),
    Path("pipeline_summary.tsv"),
    Path("summary_report.html"),
)


@pytest.mark.parametrize(
    ("archive_results", "archive_format", "archive_suffix", "shutil_format"),
    [
        (False, "not-consulted", None, None),
        (True, "zip", ".zip", "zip"),
        (True, "tar.gz", ".tar.gz", "gztar"),
    ],
)
def test_cleanup_plan_is_the_exact_authoritative_destination_set(
    tmp_path: Path,
    archive_results: bool,
    archive_format: str,
    archive_suffix: str | None,
    shutil_format: str | None,
) -> None:
    """Cleanup and both log guards consume this one literal plan."""
    output = tmp_path / "out"

    plan = plan_advntr_cleanup(output, archive_results=archive_results, archive_format=archive_format)

    assert plan.public_outputs == tuple(output / relative for relative in PUBLIC_OUTPUTS)
    assert plan.model_snapshot == output / "advntr" / "advntr_model.db"
    if archive_suffix is None:
        assert plan.archive is None
        assert plan.cleanup_destinations == plan.public_outputs
    else:
        assert plan.archive is not None
        assert plan.archive.destination == Path(f"{output}{archive_suffix}")
        assert plan.archive.base_name == str(output)
        assert plan.archive.shutil_format == shutil_format
        assert plan.cleanup_destinations == (*plan.public_outputs, plan.archive.destination)
    assert plan.model_snapshot not in plan.cleanup_destinations
    assert plan.destructive_destinations == (*plan.cleanup_destinations, plan.model_snapshot)


def test_cleanup_plan_rejects_an_unknown_selected_archive_format(tmp_path: Path) -> None:
    """A selected archive cannot silently fall out of cleanup ownership."""
    with pytest.raises(ValueError, match="Unsupported archive format: rar"):
        plan_advntr_cleanup(tmp_path / "out", archive_results=True, archive_format="rar")


def test_cleanup_guard_rejects_a_hardlink_alias_and_accepts_no_file_log(tmp_path: Path) -> None:
    """Filesystem identity catches aliases whose normalized spellings remain different."""
    output = tmp_path / "out"
    plan = plan_advntr_cleanup(output, archive_results=False, archive_format="not-consulted")
    cleanup_destination = output / "pipeline_summary.json"
    cleanup_destination.parent.mkdir()
    cleanup_destination.write_text("active log\n", encoding="utf-8")
    hardlink_alias = tmp_path / "application.log"
    hardlink_alias.hardlink_to(cleanup_destination)

    with pytest.raises(ValueError, match="aliases an adVNTR destructive preflight destination"):
        validate_pipeline_log_outside_advntr_preflight(hardlink_alias, plan)

    validate_pipeline_log_outside_advntr_preflight(None, plan)
