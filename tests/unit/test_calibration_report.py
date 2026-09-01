"""Deterministic static calibration evidence reports."""

from pathlib import Path

import pytest

from vntyper.scripts.calibration_report import (
    decode_calibration_report,
    render_calibration_report,
    write_calibration_report,
)

pytestmark = pytest.mark.unit


def _report() -> dict[str, object]:
    return {
        "schema_version": "calibration-report-v1",
        "phase": "held-out",
        "profile_sha256": "a" * 64,
        "protocol_sha256": "b" * 64,
        "evidence_sha256": "c" * 64,
        "objective": "lexicographic-safety-v1",
        "tier_metrics": [
            {"tier": tier, "displayed": 4 - index, "exact": 3 - index, "wrong": 1}
            for index, tier in enumerate(("A", "B", "C"))
        ],
        "abstentions": [
            {"split": "locked-heldout", "reason": "record-tie", "count": 2, "rate": "1/10"},
            {"split": "locked-heldout", "reason": "xd-missingness", "count": 1, "rate": "1/20"},
        ],
        "provenance": {
            "software_versions": ["VNtyper 2.0.26"],
            "reference_versions": ["hg19 sha256:reference"],
            "sample_composition": ["mutated=20", "controls=10"],
            "assays": ["capture-short-read"],
            "depths": ["30x", "60x"],
            "read_lengths": ["150 bp paired-end"],
            "independent_array_size": ["measured independently for 12/30 members"],
            "mutation_classes": ["duplication", "deletion"],
            "manifest_hashes": ["partition=" + "d" * 64, "dataset=" + "e" * 64],
            "access_attempts": ["locked-heldout: one precommitted attempt"],
            "boundary_coverage": ["count margin 0/1/2", "record share 0.49/0.50/0.51"],
            "seeds": ["bootstrap seed=295"],
        },
        "statistics": {
            "intervals": [
                "two-sided 95% exact recovery interval 0.80–0.99",
                "one-sided paired lower bound 0 percentage points",
            ],
            "roc_rows": ["threshold=1/2 TP=10 FP=0 TN=10 FN=0"],
            "pr_rows": ["threshold=1/2 precision=1 recall=1"],
            "joint_surface_rows": ["count_margin=1 share=1/2 xd_veto=disabled"],
        },
        "limitations": [
            "Local append-only custody guards are not proof of independent custody; closure requires a named external custodian.",
            "Small-n strata have wide uncertainty.",
            "Reporting an interval is not a clinical safety claim.",
        ],
    }


def test_static_report_contains_complete_metrics_provenance_statistics_and_limitations() -> None:
    html = render_calibration_report(decode_calibration_report(_report()))

    for tier in ("A", "B", "C"):
        assert f">{tier}<" in html
    for text in (
        "record-tie",
        "xd-missingness",
        "optional minimum k-mer depth",
        "lexicographic-safety-v1",
        "fitted versus validation versus held-out",
        "capture-short-read",
        "150 bp paired-end",
        "measured independently",
        "one-sided paired lower bound",
        "two-sided 95%",
        "count margin 0/1/2",
        "bootstrap seed=295",
        "Small-n",
        "Local append-only custody guards are not proof of independent custody",
        "Reporting an interval is not a clinical safety claim.",
    ):
        assert text in html
    assert "<script" not in html.casefold()
    assert "http://" not in html and "https://" not in html


def test_report_is_deterministic_escaped_and_written_byte_exactly(tmp_path: Path) -> None:
    raw = _report()
    provenance = raw["provenance"]
    assert isinstance(provenance, dict)
    provenance["assays"] = ["<script>alert(1)</script>"]
    report = decode_calibration_report(raw)

    first = render_calibration_report(report)
    second = render_calibration_report(report)
    destination = tmp_path / "report.html"
    write_calibration_report(destination, report)

    assert first == second
    assert destination.read_text(encoding="utf-8") == first
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in first
    assert "<script>alert" not in first


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("phase", "training", "phase"),
        ("profile_sha256", "A" * 64, "SHA-256"),
        ("objective", "f1", "objective"),
        ("limitations", [], "limitations"),
    ],
)
def test_report_contract_rejects_unknown_semantics(field: str, value: object, message: str) -> None:
    raw = _report()
    raw[field] = value

    with pytest.raises(ValueError, match=message):
        decode_calibration_report(raw)
