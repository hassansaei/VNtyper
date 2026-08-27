"""Report templates resolve from the installed package, never the process CWD.

The shipped ``paths.template_dir`` value is the historical relative spelling
``vntyper/templates``.  Treating that spelling as an ordinary filesystem path
made both renderers fail outside a source checkout and let an unrelated writable
CWD substitute the report's top-level template.  Explicit operator overrides
remain supported, but a missing entry template must fail instead of silently
falling back to the packaged report.
"""

from __future__ import annotations

import argparse
import copy
from pathlib import Path

import pandas as pd
import pytest

import vntyper
from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.cli import load_config
from vntyper.scripts import cli_report, report_assets
from vntyper.scripts.cohort_summary import generate_cohort_summary_report
from vntyper.scripts.generate_report import generate_summary_report

pytestmark = pytest.mark.unit

PACKAGED_TEMPLATES = Path(vntyper.__file__).resolve().parent / "templates"


def _render_sample(output_dir: Path, template_dir: str | Path | None) -> str:
    generate_summary_report(
        output_dir=str(output_dir),
        template_dir=template_dir,
        report_file="summary_report.html",
        log_file=None,
        config=load_config(None),
    )
    return (output_dir / "summary_report.html").read_text(encoding="utf-8")


def _render_cohort(output_dir: Path, config: dict) -> str:
    generate_cohort_summary_report(
        output_dir=str(output_dir),
        kestrel_df=pd.DataFrame([{"Sample": "s1", "Confidence": "High_Precision", "Flag": "Not flagged"}]),
        advntr_df=pd.DataFrame(),
        summary_file="cohort_summary.html",
        config=config,
    )
    return (output_dir / "cohort_summary.html").read_text(encoding="utf-8")


def test_absent_and_shipped_template_settings_resolve_to_installed_package() -> None:
    """Replacing ``__file__``-relative resolution with a CWD path breaks this."""
    packaged = str(PACKAGED_TEMPLATES)

    assert report_assets.template_search_paths(None) == [packaged]
    assert report_assets.template_search_paths("") == [packaged]
    assert report_assets.template_search_paths("vntyper/templates") == [packaged]
    assert report_assets.template_search_paths(PACKAGED_TEMPLATES) == [packaged]


def test_an_operator_template_is_first_and_packaged_includes_remain_available(tmp_path: Path) -> None:
    """Dropping either the override or package fallback breaks custom templates."""
    custom = tmp_path / "custom"
    custom.mkdir()
    (custom / "report_template.html").write_text(
        '<!doctype html><style>{% include "_report_base.html" %}</style><p>operator report</p>',
        encoding="utf-8",
    )

    assert report_assets.template_search_paths(custom, entry_template="report_template.html") == [
        str(custom),
        str(PACKAGED_TEMPLATES),
    ]

    rendered = _render_sample(tmp_path, custom)
    assert "operator report" in rendered
    assert "--ink:" in rendered


def test_per_sample_report_ignores_a_shadow_template_in_an_unrelated_cwd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The reproduced substitution defect: the legacy literal followed the CWD."""
    shadow = tmp_path / "vntyper" / "templates"
    shadow.mkdir(parents=True)
    (shadow / "report_template.html").write_text("SHADOW SAMPLE REPORT", encoding="utf-8")
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    monkeypatch.chdir(tmp_path)

    rendered = _render_sample(output_dir, "vntyper/templates")

    assert rendered.lstrip().startswith("<!DOCTYPE html>")
    assert "SHADOW SAMPLE REPORT" not in rendered


def test_cohort_report_ignores_a_shadow_template_in_an_unrelated_cwd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The cohort renderer has the same package-relative top-level contract."""
    shadow = tmp_path / "vntyper" / "templates"
    shadow.mkdir(parents=True)
    (shadow / "cohort_summary_template.html").write_text("SHADOW COHORT REPORT", encoding="utf-8")
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    monkeypatch.chdir(tmp_path)

    rendered = _render_cohort(output_dir, load_config(None))

    assert rendered.lstrip().startswith("<!DOCTYPE html>")
    assert "SHADOW COHORT REPORT" not in rendered


def test_per_sample_override_missing_its_entry_template_fails(tmp_path: Path) -> None:
    """A typoed override must not silently render the packaged sample report."""
    override = tmp_path / "partial"
    override.mkdir()

    with pytest.raises(ValueError, match="operator template directory.*report_template.html"):
        _render_sample(tmp_path, override)

    assert not (tmp_path / "summary_report.html").exists()


def test_cohort_override_missing_its_entry_template_fails(tmp_path: Path) -> None:
    """A typoed override must not silently render the packaged cohort report."""
    override = tmp_path / "partial"
    override.mkdir()

    with pytest.raises(ValueError, match="operator template directory.*cohort_summary_template.html"):
        _render_cohort(tmp_path, {"paths": {"template_dir": str(override)}})

    assert not (tmp_path / "cohort_summary.html").exists()


def test_pipeline_propagates_an_absent_template_setting_as_packaged_default(tmp_path: Path) -> None:
    """Restoring the old literal at the pipeline call site makes CWD significant."""
    config = copy.deepcopy(MINIMAL_CONFIG)
    config.pop("paths")

    harness = run_pipeline_under_harness(tmp_path / "run", config=config)

    assert harness.error is None
    assert harness.positional("generate_summary_report")[1] is None


def test_report_handler_propagates_an_absent_template_setting_as_packaged_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The standalone report command and pipeline must choose the same default."""
    captured: dict[str, object] = {}

    def _capture(**kwargs: object) -> None:
        captured.update(kwargs)

    monkeypatch.setattr(cli_report, "generate_summary_report", _capture)
    args = argparse.Namespace(
        output_dir=str(tmp_path),
        input_dir=None,
        report_file="summary_report.html",
        bed_file=None,
        bam_file=None,
        reference_fasta=None,
        flanking=50,
    )

    cli_report.handle_report(
        args,
        config={"default_values": {}},
        parser=argparse.ArgumentParser(),
        log_level_value=20,
        log_file_str=None,
    )

    assert captured["template_dir"] is None
