"""Focused owned-boundary tests for unreadable alignment headers."""

from __future__ import annotations

import json
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts.pipeline_alignment import prepare_input_alignment_preflight

pytestmark = pytest.mark.unit


def test_unreadable_cram_header_cannot_bypass_remote_uri_or_assembly_guards(tmp_path: Path) -> None:
    """CRAM policy cannot proceed when the header needed to enforce it is absent."""
    output = tmp_path / "run-output"
    output.mkdir()

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", return_value=None),
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_header_reference_policy") as uri_policy,
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_declared_assembly") as assembly,
        mock.patch("vntyper.scripts.pipeline_alignment.prepare_alignment_target") as target,
        mock.patch("vntyper.scripts.pipeline_alignment.run_preflight") as preflight,
        pytest.raises(ValueError, match="CRAM header could not be read"),
    ):
        prepare_input_alignment_preflight(
            in_path=tmp_path / "patient-input" / "sample.cram",
            input_type="CRAM",
            output_dir=output,
            config={},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8"))["code"] == (
        "alignment_header_invalid"
    )
    uri_policy.assert_not_called()
    assembly.assert_not_called()
    target.assert_not_called()
    preflight.assert_not_called()
