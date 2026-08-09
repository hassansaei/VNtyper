"""Focused contracts for golden-cohort causal diagnostics."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest import mock

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from golden_cohort import runner  # noqa: E402
from golden_cohort.case_expectations import MIXED_LAYOUT_FAILURE_DIAGNOSTIC  # noqa: E402

from vntyper.scripts import pipeline_read_routing  # noqa: E402


def _case(case_id: str, **extra: object) -> dict[str, object]:
    """Build the minimum pipeline case consumed by the side runner."""
    return {
        "case_id": case_id,
        "kind": "pipeline",
        "group": "base",
        "bam": f"/data/example_{case_id}.bam",
        "assembly": "hg19",
        "fast_mode": True,
        "advntr": False,
        "expect_exit": "zero",
        "required_artifacts": [],
        **extra,
    }


@pytest.mark.parametrize(
    ("selected", "problem"),
    [
        ({"required_artifacts": []}, "expect_exit"),
        ({"expect_exit": "maybe", "required_artifacts": []}, "expect_exit"),
        ({"expect_exit": "zero", "required_artifacts": "summary.json"}, "required_artifacts"),
        (
            {"expect_exit": "nonzero", "required_artifacts": [], "expected_stderr_contains": ""},
            "expected_stderr_contains",
        ),
        (
            {"expect_exit": "nonzero", "required_artifacts": [], "expected_stderr_contains": ["mixed"]},
            "expected_stderr_contains",
        ),
        (
            {"expect_exit": "nonzero", "required_artifacts": [], "expected_stderr_contains": None},
            "expected_stderr_contains",
        ),
        (
            {"expect_exit": "nonzero", "required_artifacts": [], "expected_stderr_contains": "   "},
            "expected_stderr_contains",
        ),
    ],
)
def test_side_expectation_rejects_malformed_outcome_fields(selected: object, problem: str) -> None:
    """A partial side declaration must not inherit a legacy outcome and pass silently."""
    case = {"case_id": "mixed", "side_expectations": {"after": selected}}

    with pytest.raises(ValueError, match=problem):
        runner.materialize_side_expectation(case, "after")


@pytest.mark.parametrize(
    ("stderr", "expected_met"),
    [
        (f"{MIXED_LAYOUT_FAILURE_DIAGNOSTIC} Produced FASTQs: ...", True),
        ("reference could not be decoded", False),
        ("", False),
    ],
)
def test_run_one_uses_actual_stderr_to_prove_the_declared_mixed_layout_failure(
    tmp_path: Path, stderr: str, expected_met: bool
) -> None:
    """Exit one is causal evidence only when the declared diagnostic was emitted."""
    case = _case(
        "mixed",
        expect_exit="nonzero",
        expected_stderr_contains=MIXED_LAYOUT_FAILURE_DIAGNOSTIC,
    )
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n",
        stderr=stderr,
        returncode=1,
    )

    with (
        mock.patch.object(runner.subprocess, "run", return_value=pipeline_result),
        mock.patch.object(runner.time, "monotonic", side_effect=[100.0, 101.0]),
    ):
        record = runner._run_one(
            case=case,
            argv=["pipeline", "--bam", str(case["bam"])],
            tree=tmp_path,
            side="after",
            marker="vntyper.scripts.alignment_contract",
            expect_marker=True,
            output_dir=tmp_path / "out",
            log_dir=tmp_path / "logs",
            timeout=60,
        )

    assert record["expectation_met"] is expected_met
    assert bool(record["expectation_problems"]) is not expected_met


def test_production_mixed_layout_failure_starts_with_the_golden_diagnostic(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The harness's causal oracle must stay tied to the production refusal."""
    produced = (
        "/r1.fastq.gz",
        "/r2.fastq.gz",
        "/other.fastq.gz",
        "/single.fastq.gz",
    )
    counts = dict(zip(produced, (2, 2, 0, 1), strict=True))
    monkeypatch.setattr(
        pipeline_read_routing,
        "count_fastq_records",
        lambda path, *, lines_per_record: counts[str(path)],
    )

    with pytest.raises(ValueError) as raised:
        pipeline_read_routing.route_converted_fastqs(produced, config={})

    assert str(raised.value).startswith(f"{MIXED_LAYOUT_FAILURE_DIAGNOSTIC} Produced FASTQs:")
