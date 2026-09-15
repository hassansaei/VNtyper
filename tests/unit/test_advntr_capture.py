"""Tests for the shell-free adVNTR completed-capture adapter."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from tests.unit.test_advntr_calibration_policy import caller_values, capabilities, capture_policy
from vntyper.modules.advntr import advntr_capture as capture
from vntyper.modules.advntr.advntr_calibration_policy import (
    AdvntrToolPin,
    decode_advntr_capabilities,
    decode_capture_policy,
)
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values

pytestmark = pytest.mark.unit


def _plan(tmp_path: Path, *, mode: str = "legacy", additional_commands: str = "") -> capture.CapturePlan:
    model = tmp_path / "model db.sqlite"
    model.write_bytes(b"synthetic model")
    background = None
    background_sha = None
    if mode == "exact":
        background = tmp_path / "background model.json"
        background.write_bytes(b"synthetic background")
        background_sha = hashlib.sha256(background.read_bytes()).hexdigest()
    raw_capture = capture_policy()
    raw_capture["parameters"]["caller_mode"] = mode  # type: ignore[index]
    caller = decode_caller_policy_values(caller_values(mode=mode))
    projected = capture.capture_policy_for_caller(decode_capture_policy(raw_capture), caller)
    return capture.build_capture_plan(
        ("/env with spaces/python", "-m", "advntr"),
        alignment_path=tmp_path / "input specimen.bam",
        model_path=model,
        result_path=tmp_path / "production result.vcf",
        sink_path=tmp_path / "complete capture.jsonl",
        working_directory=tmp_path / "working directory",
        vntr_ids=(17, 25561),
        capture_policy=projected,
        caller_policy=caller,
        model_sha256=hashlib.sha256(model.read_bytes()).hexdigest(),
        background_path=background,
        background_sha256=background_sha,
        additional_commands=additional_commands,
    )


def _record(plan: capture.CapturePlan, vntr_id: int) -> dict[str, object]:
    capabilities_document = capabilities()
    loaded_background: dict[str, object] | None = None
    if plan.background_sha256 is not None:
        loaded_background = {
            "schema": "advntr.frameshift.background",
            "version": 1,
            "default_probability": 0.1,
            "states": {},
        }
    model_locus = {"opaque": "preserved"}
    return {
        "schema_version": "advntr-frameshift-capture-v2",
        "completion": "completed-vntr",
        "producer": {key: capabilities_document[key] for key in ("package_version", "build_id", "source_revision")},
        "assets": {
            "model_sha256": plan.model_sha256,
            "model_locus_sha256": capture.advntr_canonical_sha256(model_locus),
            "background_sha256": plan.background_sha256,
            "loaded_background_sha256": (
                None if loaded_background is None else capture.advntr_canonical_sha256(loaded_background)
            ),
        },
        "model_locus": model_locus,
        "loaded_background": loaded_background,
        "capture_policy": capture.capture_policy_document(plan.capture_policy),
        "caller_policy": capture.caller_policy_document(plan.caller_policy),
        "locus": {"vntr_id": vntr_id, "read_length": 150, "is_haploid": False, "selected_read_count": 0},
        "unit_geometry": {},
        "reference_order": [],
        "flank_boundaries": {},
        "warnings": ["synthetic audit warning"],
        "occurrences": [],
        "spans": [],
        "evidence_rows": [],
        "candidate_traversal": {},
        "decision_visits": [],
    }


def test_builds_exact_shell_free_capture_argv_with_path_spaces(tmp_path: Path) -> None:
    plan = _plan(tmp_path)

    assert plan.argv[:5] == ("/env with spaces/python", "-m", "advntr", "genotype", "-fs")
    assert plan.argv[plan.argv.index("--alignment_file") + 1] == str(tmp_path / "input specimen.bam")
    assert plan.argv[plan.argv.index("--frameshift-calibration-out") + 1] == str(tmp_path / "complete capture.jsonl")
    assert plan.argv[plan.argv.index("-vid") + 1] == "17,25561"
    assert "--frameshift-capture-version" in plan.argv
    assert "--exact-frameshift-caller" not in plan.argv


def test_capture_plan_refuses_extension_flags_and_wrong_background_mode(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="additional_commands"):
        _plan(tmp_path, additional_commands="--threads 99")

    with pytest.raises(ValueError, match="background"):
        capture.build_capture_plan(
            ("advntr",),
            alignment_path=tmp_path / "a.bam",
            model_path=tmp_path / "model",
            result_path=tmp_path / "result",
            sink_path=tmp_path / "sink",
            working_directory=tmp_path / "work",
            vntr_ids=(17,),
            capture_policy=_plan(tmp_path).capture_policy,
            caller_policy=_plan(tmp_path).caller_policy,
            model_sha256="a" * 64,
            background_path=tmp_path / "background",
            background_sha256="b" * 64,
        )


def test_decodes_complete_sink_and_preserves_each_opaque_record(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    rows = [_record(plan, 17), _record(plan, 25561)]
    raw = b"".join(json.dumps(row, sort_keys=True, separators=(",", ":")).encode() + b"\n" for row in rows)

    decoded = capture.decode_capture_sink(raw, plan, decode_advntr_capabilities(capabilities()))

    assert decoded.sink_sha256 == hashlib.sha256(raw).hexdigest()
    assert tuple(row.vntr_id for row in decoded.records) == (17, 25561)
    assert capture.capture_record_document(decoded.records[0]) == rows[0]
    assert capture.capture_record_document(decoded.records[0])["warnings"] == ["synthetic audit warning"]


def test_sink_requires_canonical_record_bytes_and_loaded_asset_digests(tmp_path: Path) -> None:
    plan = _plan(tmp_path, mode="exact")
    rows = [_record(plan, 17), _record(plan, 25561)]
    noncanonical = b"".join(json.dumps(row).encode() + b"\n" for row in rows)

    with pytest.raises(ValueError, match="canonical"):
        capture.decode_capture_sink(noncanonical, plan, decode_advntr_capabilities(capabilities()))

    rows[0]["assets"]["loaded_background_sha256"] = "0" * 64  # type: ignore[index]
    raw = b"".join(json.dumps(row, sort_keys=True, separators=(",", ":")).encode() + b"\n" for row in rows)
    with pytest.raises(ValueError, match="loaded background"):
        capture.decode_capture_sink(raw, plan, decode_advntr_capabilities(capabilities()))


@pytest.mark.parametrize("failure", ["partial", "duplicate", "wrong_model", "wrong_policy"])
def test_sink_rejects_partial_duplicate_or_mismatched_evidence(tmp_path: Path, failure: str) -> None:
    plan = _plan(tmp_path)
    rows = [_record(plan, 17), _record(plan, 25561)]
    if failure == "duplicate":
        rows[1]["locus"]["vntr_id"] = 17  # type: ignore[index]
    elif failure == "wrong_model":
        rows[0]["assets"]["model_sha256"] = "0" * 64  # type: ignore[index]
    elif failure == "wrong_policy":
        rows[0]["caller_policy"]["minimum_read_support"] = 99  # type: ignore[index]
    raw = b"".join(json.dumps(row, sort_keys=True, separators=(",", ":")).encode() + b"\n" for row in rows)
    if failure == "partial":
        raw = raw.rstrip(b"\n")

    with pytest.raises(ValueError):
        capture.decode_capture_sink(raw, plan, decode_advntr_capabilities(capabilities()))


def test_run_preflights_then_executes_and_requires_both_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan = _plan(tmp_path)
    plan.alignment_path.write_bytes(b"alignment")
    plan.working_directory.mkdir()
    events: list[str] = []

    def probe(*_args: object, **_kwargs: object) -> object:
        events.append("probe")
        return decode_advntr_capabilities(capabilities())

    def runner(argv: tuple[str, ...], **kwargs: object) -> subprocess.CompletedProcess[str]:
        events.append("run")
        assert kwargs == {"cwd": plan.working_directory, "capture_output": True, "text": True, "check": False}
        plan.result_path.write_text("synthetic result\n")
        rows = [_record(plan, 17), _record(plan, 25561)]
        plan.sink_path.write_text(
            "".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in rows)
        )
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(capture, "probe_advntr_capabilities", probe)
    result = capture.run_capture(plan, AdvntrToolPin("2.4.0", "a" * 64, "b" * 40), runner=runner)

    assert events == ["probe", "run"]
    assert len(result.records) == 2


def test_failed_process_does_not_decode_partial_sink(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    plan = _plan(tmp_path)
    plan.alignment_path.write_bytes(b"alignment")
    plan.working_directory.mkdir()
    monkeypatch.setattr(
        capture, "probe_advntr_capabilities", lambda *_args, **_kwargs: decode_advntr_capabilities(capabilities())
    )

    def runner(argv: tuple[str, ...], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        plan.sink_path.write_text("partial private evidence")
        return subprocess.CompletedProcess(argv, 1, "", "private stderr")

    with pytest.raises(RuntimeError, match="capture process failed"):
        capture.run_capture(plan, AdvntrToolPin("2.4.0", "a" * 64, "b" * 40), runner=runner)
