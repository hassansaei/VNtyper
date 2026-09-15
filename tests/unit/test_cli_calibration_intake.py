"""CLI adapter for strict local calibration intake bundles."""

import json
import os
import stat
from copy import deepcopy
from importlib import import_module
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.unit.test_calibration_intake_contract import synthetic_intake
from vntyper import cli
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

pytestmark = pytest.mark.unit


def test_intake_parser_exposes_the_complete_explicit_adapter_contract() -> None:
    parser = import_module("vntyper.scripts.cli_parser").build_parser()
    args = parser.parse_args(
        [
            "calibrate",
            "intake",
            "--manifest",
            "intake.json",
            "--output",
            "bundle",
            "--preprocessing-priority",
            "raw-v1",
            "processed-v2",
            "--cram-references",
            "references.json",
            "--temporary-directory",
            "temporary",
        ]
    )

    assert args.calibration_operation == "intake"
    assert args.manifest == Path("intake.json")
    assert args.output == Path("bundle")
    assert args.preprocessing_priority == ["raw-v1", "processed-v2"]
    assert args.cram_references == Path("references.json")
    assert args.temporary_directory == Path("temporary")


@pytest.mark.parametrize("missing", ["manifest", "output", "priority"])
def test_intake_parser_reports_missing_required_arguments_as_usage_errors(missing: str) -> None:
    arguments = {
        "manifest": ["--manifest", "intake.json"],
        "output": ["--output", "bundle"],
        "priority": ["--preprocessing-priority", "raw-v1"],
    }
    command = ["calibrate", "intake"]
    for name, values in arguments.items():
        if name != missing:
            command.extend(values)

    with pytest.raises(SystemExit) as failure:
        import_module("vntyper.scripts.cli_parser").build_parser().parse_args(command)
    assert failure.value.code == 2


def test_cram_reference_document_is_closed_strict_and_immutable(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.cli_calibration_intake")
    reference = tmp_path / "synthetic.fa"
    document = tmp_path / "references.json"
    document.write_bytes(
        canonical_json_bytes(
            {
                "schema_version": "calibration-cram-references-v1",
                "references": {"assembly-a": {"path": str(reference.resolve()), "sha256": "a" * 64}},
            }
        )
    )

    decoded = module.load_cram_references(document)

    assert decoded == {"assembly-a": module.PinnedCramReference(reference.resolve(), "a" * 64)}
    with pytest.raises(TypeError):
        decoded["other"] = module.PinnedCramReference(reference.resolve(), "b" * 64)
    assert module.load_cram_references(None) == {}


@pytest.mark.parametrize(
    "invalid",
    [
        "duplicate",
        "nonfinite",
        "root_field",
        "row_field",
        "version",
        "references_type",
        "assembly",
        "path",
        "digest",
    ],
)
def test_cram_reference_document_rejects_ambiguous_or_open_content_before_producer(
    tmp_path: Path, invalid: str
) -> None:
    module = import_module("vntyper.scripts.cli_calibration_intake")
    document = tmp_path / "references.json"
    reference: dict[str, object] = {"path": "/synthetic/reference.fa", "sha256": "a" * 64}
    raw: dict[str, object] = {
        "schema_version": "calibration-cram-references-v1",
        "references": {"assembly-a": reference},
    }
    if invalid == "duplicate":
        document.write_text('{"schema_version":"x","schema_version":"y","references":{}}', encoding="utf-8")
    elif invalid == "nonfinite":
        document.write_text(
            '{"schema_version":"calibration-cram-references-v1","references":{"a":{"path":NaN,"sha256":"x"}}}',
            encoding="utf-8",
        )
    else:
        if invalid == "root_field":
            raw["extra"] = True
        elif invalid == "row_field":
            reference["extra"] = True
        elif invalid == "version":
            raw["schema_version"] = "calibration-cram-references-v2"
        elif invalid == "references_type":
            raw["references"] = []
        elif invalid == "assembly":
            raw["references"] = {" ": reference}
        elif invalid == "path":
            reference["path"] = "relative.fa"
        else:
            reference["sha256"] = "A" * 64
        document.write_text(json.dumps(raw), encoding="utf-8")

    with patch.object(module, "prepare_intake_bundle") as producer, pytest.raises(ValueError):
        module.run_calibration_intake(
            _arguments(tmp_path, cram_references=document),
        )
    producer.assert_not_called()


def test_intake_adapter_passes_explicit_arguments_to_the_single_atomic_producer(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.cli_calibration_intake")
    references = tmp_path / "references.json"
    references.write_bytes(
        canonical_json_bytes(
            {
                "schema_version": "calibration-cram-references-v1",
                "references": {"assembly-a": {"path": "/synthetic/reference.fa", "sha256": "a" * 64}},
            }
        )
    )
    args = _arguments(tmp_path, cram_references=references)

    with patch.object(module, "prepare_intake_bundle") as producer:
        module.run_calibration_intake(args)

    producer.assert_called_once_with(
        args.manifest,
        args.output,
        preprocessing_priority=("raw-v1", "processed-v2"),
        cram_references={"assembly-a": module.PinnedCramReference(Path("/synthetic/reference.fa"), "a" * 64)},
        temporary_parent=args.temporary_directory,
    )


@pytest.mark.parametrize("invalid", ["namespace", "manifest", "priority", "references", "temporary"])
def test_intake_adapter_revalidates_direct_argument_objects_before_production(tmp_path: Path, invalid: str) -> None:
    module = import_module("vntyper.scripts.cli_calibration_intake")
    args = _arguments(tmp_path)
    if invalid == "namespace":
        args = object()
    elif invalid == "manifest":
        args.manifest = "intake.json"
    elif invalid == "priority":
        args.preprocessing_priority = ("raw-v1",)
    elif invalid == "references":
        args.cram_references = "references.json"
    else:
        args.temporary_directory = "temporary"

    with patch.object(module, "prepare_intake_bundle") as producer, pytest.raises(ValueError):
        module.run_calibration_intake(args)
    producer.assert_not_called()


def test_cli_malformed_intake_content_exits_one_without_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "intake.json"
    manifest.write_text('{"schema_version":"calibration-intake-v1","unknown":true}', encoding="utf-8")
    output = tmp_path / "bundle"
    monkeypatch.setattr(cli, "setup_logging", lambda log_level, log_file: None)

    with pytest.raises(SystemExit) as failure:
        cli.main(
            [
                "calibrate",
                "intake",
                "--manifest",
                str(manifest),
                "--output",
                str(output),
                "--preprocessing-priority",
                "raw-v1",
            ]
        )

    assert failure.value.code == 1
    assert not output.exists()


def test_cram_reference_fifo_fails_before_intake_artifact_reads(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.cli_calibration_intake")
    references = tmp_path / "references.fifo"
    os.mkfifo(references)

    with patch.object(module, "prepare_intake_bundle") as producer, pytest.raises(ValueError, match="regular"):
        module.run_calibration_intake(_arguments(tmp_path, cram_references=references))
    producer.assert_not_called()


def test_cli_builds_real_synthetic_fastq_bundle_without_reading_locked_membership(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = tmp_path / "invented R1.fastq"
    second = tmp_path / "invented R2.fastq"
    first.write_text("@invented-read/1\nACGT\n+\nIIII\n", encoding="ascii")
    second.write_text("@invented-read/2\nTGCA\n+\nJJJJ\n", encoding="ascii")
    raw = synthetic_intake()
    raw["artifacts"][0].update(
        path=str(first.resolve()),
        format="FASTQ_PAIR",
        mate_path=str(second.resolve()),
        expected_sha256=None,
    )
    locked_specimen = deepcopy(raw["specimens"][0])
    locked_specimen.update(key="locked-member", individual_key="locked-individual")
    locked_assignment = deepcopy(raw["assignments"][0])
    locked_assignment.update(
        specimen_key="locked-member",
        role="locked-heldout",
        provenance="external-custodian",
        groups={name: [f"{name}:locked"] for name in locked_assignment["groups"]},
    )
    raw["specimens"].append(locked_specimen)
    raw["assignments"].append(locked_assignment)
    manifest = tmp_path / "intake.json"
    manifest.write_bytes(canonical_json_bytes(raw))
    output = tmp_path / "bundle"
    temporary = tmp_path / "temporary"
    temporary.mkdir()
    monkeypatch.setattr(cli, "setup_logging", lambda log_level, log_file: None)

    cli.main(
        [
            "calibrate",
            "intake",
            "--manifest",
            str(manifest),
            "--output",
            str(output),
            "--preprocessing-priority",
            "synthetic-preprocessing-v1",
            "--temporary-directory",
            str(temporary),
        ]
    )

    assert {path.name for path in output.iterdir()} == {
        "normalized.json",
        "dedup_audit.json",
        "partitions.json",
        "provenance.json",
    }
    assert stat.S_IMODE(output.stat().st_mode) == 0o700
    assert all(stat.S_IMODE(path.stat().st_mode) == 0o600 for path in output.iterdir())
    normalized = load_strict_json_object((output / "normalized.json").read_bytes())
    audit = load_strict_json_object((output / "dedup_audit.json").read_bytes())
    assert {row["key"] for row in normalized["specimens"]} == {"sample-001", "locked-member"}
    assert set(audit["fingerprints"]) == {"artifact-001"}
    assert audit["fingerprints"]["artifact-001"]["primary_record_count"] == 2


def _arguments(tmp_path: Path, *, cram_references: Path | None = None):
    module = import_module("argparse")
    return module.Namespace(
        manifest=tmp_path / "intake.json",
        output=tmp_path / "bundle",
        preprocessing_priority=["raw-v1", "processed-v2"],
        cram_references=cram_references,
        temporary_directory=tmp_path,
    )
