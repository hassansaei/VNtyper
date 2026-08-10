"""Behavioral tests for exact-SHA Docker release evidence."""

import json
from pathlib import Path

import pytest

from scripts.release_evidence import (
    eligible_successful_runs,
    main,
    summarize_evidence_state,
    validate_evidence,
    validate_image,
)

pytestmark = pytest.mark.unit

SHA = "a" * 40
OTHER_SHA = "c" * 40
DIGEST = "sha256:" + "b" * 64
OTHER_DIGEST = "sha256:" + "d" * 64
VERSION = "2.0.10"


def _evidence(**overrides: object) -> dict[str, object]:
    payload: dict[str, object] = {
        "contract_version": 1,
        "sha": SHA,
        "digest": DIGEST,
        "run_id": 41,
        "run_attempt": 2,
        "revision": SHA,
        "version": VERSION,
    }
    payload.update(overrides)
    return payload


def _image(revision: str = SHA, version: str = VERSION) -> dict[str, object]:
    return {
        "config": {
            "Labels": {
                "org.opencontainers.image.revision": revision,
                "org.opencontainers.image.version": version,
            }
        }
    }


def test_evidence_binds_successful_exact_sha_run_and_image_labels() -> None:
    """A schedule run or substituted image must not become release provenance."""
    runs = {
        "workflow_runs": [
            {
                "id": 41,
                "run_attempt": 2,
                "head_sha": SHA,
                "head_branch": "main",
                "event": "push",
                "status": "completed",
                "conclusion": "success",
            },
            {
                "id": 99,
                "run_attempt": 1,
                "head_sha": SHA,
                "head_branch": "main",
                "event": "schedule",
                "status": "completed",
                "conclusion": "success",
            },
        ]
    }

    eligible = eligible_successful_runs(runs, SHA)

    assert tuple(run["id"] for run in eligible) == (41,)
    assert (eligible[0]["id"], eligible[0]["run_attempt"]) == (41, 2)
    evidence = validate_evidence(_evidence(), sha=SHA, version=VERSION, run_id=41, run_attempt=2)
    result = validate_image(
        evidence,
        _image(),
        _image(),
        immutable_digest=DIGEST,
        short_digest=DIGEST,
    )
    assert result == {
        "source_digest": DIGEST,
        "source_revision": SHA,
        "source_version": VERSION,
        "short_tag_collision": False,
    }


def test_eligible_runs_filter_and_sort_all_pages_by_descending_numeric_id() -> None:
    """Newest-only selection would prevent recovery from a missing newest artifact."""
    base = {
        "run_attempt": 1,
        "head_sha": SHA,
        "head_branch": "main",
        "event": "push",
        "status": "completed",
        "conclusion": "success",
    }
    pages = [
        {"workflow_runs": [base | {"id": 40}, base | {"id": 42, "html_url": "https://runs/42"}]},
        {
            "workflow_runs": [
                base | {"id": 41, "run_attempt": 3, "html_url": "https://runs/41"},
                base | {"id": 90, "head_sha": OTHER_SHA},
                base | {"id": 91, "head_branch": "feature"},
                base | {"id": 92, "event": "workflow_dispatch"},
                base | {"id": 93, "status": "in_progress", "conclusion": None},
                base | {"id": 94, "conclusion": "failure"},
                base | {"id": "95"},
                base | {"id": 96, "run_attempt": "2"},
            ]
        },
    ]

    selected = eligible_successful_runs(pages, SHA)

    assert tuple(run["id"] for run in selected) == (42, 41, 40)
    assert selected[1]["run_attempt"] == 3
    assert selected[1]["html_url"] == "https://runs/41"


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("contract_version", 2),
        ("run_id", 40),
        ("run_attempt", 3),
        ("sha", OTHER_SHA),
        ("revision", OTHER_SHA),
        ("version", "2.0.9"),
        ("digest", "sha256:not-a-digest"),
    ],
)
def test_validate_evidence_rejects_each_stale_or_malformed_field(field: str, value: object) -> None:
    """Every contract identity field must independently bind the selected run and candidate."""
    with pytest.raises(ValueError, match=field):
        validate_evidence(_evidence(**{field: value}), sha=SHA, version=VERSION, run_id=41, run_attempt=2)


@pytest.mark.parametrize(
    ("image", "immutable_digest", "message"),
    [
        ({}, DIGEST, "revision"),
        ({"config": {"Labels": {"org.opencontainers.image.version": VERSION}}}, DIGEST, "revision"),
        (_image(OTHER_SHA), DIGEST, "revision"),
        (_image(SHA, "2.0.9"), DIGEST, "version"),
        (_image(), OTHER_DIGEST, "immutable_digest"),
    ],
)
def test_validate_image_rejects_missing_labels_and_immutable_manifest_drift(
    image: dict[str, object],
    immutable_digest: str,
    message: str,
) -> None:
    """Registry configuration and manifest identity must both match evidence."""
    with pytest.raises(ValueError, match=message):
        validate_image(
            _evidence(),
            image,
            _image(),
            immutable_digest=immutable_digest,
            short_digest=DIGEST,
        )


def test_equal_short_prefix_with_other_digest_is_reported_not_substituted() -> None:
    """A proven seven-character collision is diagnostic, not a provenance substitute."""
    other_sha = SHA[:7] + "d" * 33

    result = validate_image(
        _evidence(),
        _image(),
        _image(other_sha, "2.0.9"),
        immutable_digest=DIGEST,
        short_digest=OTHER_DIGEST,
    )

    assert result["source_digest"] == DIGEST
    assert result["source_revision"] == SHA
    assert result["short_tag_collision"] is True


def test_collision_requires_a_canonical_registry_short_digest() -> None:
    """A revision label alone cannot legitimize a malformed registry observation."""
    other_sha = SHA[:7] + "d" * 33

    with pytest.raises(ValueError, match="short_digest"):
        validate_image(
            _evidence(),
            _image(),
            _image(other_sha, "2.0.9"),
            immutable_digest=DIGEST,
            short_digest="not-a-manifest-digest",
        )


@pytest.mark.parametrize(
    "short_image",
    [
        {},
        {"config": {"Labels": {"org.opencontainers.image.version": VERSION}}},
        _image(OTHER_SHA),
        _image(SHA),
        _image(SHA[:7] + "g" * 33),
    ],
)
def test_unproven_short_tag_digest_mismatch_fails_closed(short_image: dict[str, object]) -> None:
    """Only a different valid full SHA sharing the exact prefix proves collision."""
    with pytest.raises(ValueError, match="unexplained short-tag drift"):
        validate_image(
            _evidence(),
            _image(),
            short_image,
            immutable_digest=DIGEST,
            short_digest=OTHER_DIGEST,
        )


@pytest.mark.parametrize(
    ("evidence", "message"),
    [
        ({"digest": "invalid", "revision": SHA, "version": VERSION}, "evidence digest"),
        ({"digest": DIGEST, "revision": "invalid", "version": VERSION}, "evidence revision"),
        ({"digest": DIGEST, "revision": SHA, "version": ""}, "evidence version"),
    ],
)
def test_image_validation_rejects_unvalidated_evidence_identity(
    evidence: dict[str, object],
    message: str,
) -> None:
    """Direct callers cannot bypass the evidence validator with malformed identity."""
    with pytest.raises(ValueError, match=message):
        validate_image(evidence, _image(), _image(), immutable_digest=DIGEST, short_digest=DIGEST)


def test_matching_short_digest_still_requires_matching_labels() -> None:
    """A contradictory config observation must fail even when manifest digests match."""
    with pytest.raises(ValueError, match="unexplained short-tag drift"):
        validate_image(
            _evidence(),
            _image(),
            _image(OTHER_SHA),
            immutable_digest=DIGEST,
            short_digest=DIGEST,
        )


def test_summary_marks_missing_phases_and_preserves_available_identity() -> None:
    """Failure summaries must retain the selected run and evidence even if image inspection failed."""
    summary = summarize_evidence_state(
        sha=SHA,
        step_outcome="failure",
        selected_run={"id": 41, "run_attempt": 2, "html_url": "https://runs/41"},
        validated_evidence=_evidence(),
        image=None,
    )

    assert summary == {
        "sha": SHA,
        "state": "failed",
        "step_outcome": "failure",
        "selected_run": {
            "available": True,
            "id": 41,
            "run_attempt": 2,
            "html_url": "https://runs/41",
        },
        "evidence": {
            "available": True,
            "digest": DIGEST,
            "revision": SHA,
            "version": VERSION,
        },
        "image": {"available": False, "revision": None, "version": None},
    }


def test_summary_makes_all_missing_local_phase_files_explicit() -> None:
    """An early shell failure must not yield an ambiguous empty evidence summary."""
    summary = summarize_evidence_state(
        sha=SHA,
        step_outcome="failure",
        selected_run=None,
        validated_evidence=None,
        image=None,
    )

    assert summary["state"] == "failed"
    assert summary["selected_run"] == {
        "available": False,
        "id": None,
        "run_attempt": None,
        "html_url": None,
    }
    assert summary["evidence"] == {
        "available": False,
        "digest": None,
        "revision": None,
        "version": None,
    }
    assert summary["image"] == {"available": False, "revision": None, "version": None}


def test_summary_reports_verified_only_when_every_phase_is_available() -> None:
    """A successful step with complete local identity is explicitly verified."""
    summary = summarize_evidence_state(
        sha=SHA,
        step_outcome="success",
        selected_run={"id": 41, "run_attempt": 2, "html_url": "https://runs/41"},
        validated_evidence=_evidence(),
        image=_image(),
    )

    assert summary["state"] == "verified"
    assert summary["image"] == {"available": True, "revision": SHA, "version": VERSION}


def test_cli_validates_image_and_writes_exact_five_outputs(tmp_path: Path) -> None:
    """The workflow boundary must receive only the validated immutable source identity."""
    evidence_path = tmp_path / "evidence.json"
    image_path = tmp_path / "image.json"
    short_image_path = tmp_path / "short-image.json"
    output_path = tmp_path / "github-output"
    evidence_path.write_text(json.dumps(_evidence()), encoding="utf-8")
    image_path.write_text(json.dumps(_image()), encoding="utf-8")
    short_image_path.write_text(json.dumps(_image()), encoding="utf-8")

    result = main(
        [
            "validate-image",
            str(evidence_path),
            str(image_path),
            str(short_image_path),
            "--immutable-digest",
            DIGEST,
            "--short-digest",
            DIGEST,
            "--image",
            "ghcr.io/hassansaei/vntyper",
            "--output",
            str(output_path),
        ]
    )

    assert result == 0
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        f"source_ref=ghcr.io/hassansaei/vntyper@{DIGEST}",
        f"source_digest={DIGEST}",
        f"source_revision={SHA}",
        f"source_version={VERSION}",
        "short_tag_collision=false",
    ]


def test_cli_artifact_present_requires_the_exact_unexpired_name(tmp_path: Path) -> None:
    """A similarly named, expired, or unrelated artifact must not satisfy recovery."""
    artifacts = tmp_path / "artifacts.json"
    artifacts.write_text(
        json.dumps(
            [
                {
                    "artifacts": [
                        {"name": f"docker-release-evidence-{SHA}-2-old", "expired": False},
                        {"name": f"docker-release-evidence-{SHA}-2", "expired": True},
                    ]
                }
            ]
        ),
        encoding="utf-8",
    )

    assert main(["artifact-present", str(artifacts), "--name", f"docker-release-evidence-{SHA}-2"]) == 1
    payload = json.loads(artifacts.read_text(encoding="utf-8"))
    payload[0]["artifacts"][1]["expired"] = False
    artifacts.write_text(json.dumps(payload), encoding="utf-8")
    assert main(["artifact-present", str(artifacts), "--name", f"docker-release-evidence-{SHA}-2"]) == 0


def test_cli_serializes_eligible_runs_validated_evidence_and_failure_summary(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Every workflow CLI command must preserve typed JSON across local-file boundaries."""
    runs_path = tmp_path / "runs.json"
    evidence_path = tmp_path / "evidence.json"
    summary_path = tmp_path / "summary.json"
    run = {
        "id": 41,
        "run_attempt": 2,
        "head_sha": SHA,
        "head_branch": "main",
        "event": "push",
        "status": "completed",
        "conclusion": "success",
        "html_url": "https://runs/41",
    }
    runs_path.write_text(json.dumps({"workflow_runs": [run]}), encoding="utf-8")
    evidence_path.write_text(json.dumps(_evidence()), encoding="utf-8")

    assert main(["eligible-runs", str(runs_path), "--sha", SHA]) == 0
    assert json.loads(capsys.readouterr().out) == [run]
    assert (
        main(
            [
                "validate-evidence",
                str(evidence_path),
                "--sha",
                SHA,
                "--version",
                VERSION,
                "--run-id",
                "41",
                "--run-attempt",
                "2",
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out) == _evidence()
    assert (
        main(
            [
                "summarize",
                "--sha",
                SHA,
                "--step-outcome",
                "failure",
                "--selected-run",
                str(tmp_path / "missing-run.json"),
                "--validated-evidence",
                str(tmp_path / "missing-evidence.json"),
                "--image",
                str(tmp_path / "missing-image.json"),
                "--output",
                str(summary_path),
            ]
        )
        == 0
    )
    assert json.loads(summary_path.read_text(encoding="utf-8"))["state"] == "failed"


@pytest.mark.parametrize("payload", [[], [1]])
def test_cli_rejects_non_object_evidence_and_pages(tmp_path: Path, payload: list[object]) -> None:
    """Malformed local JSON shapes must fail before policy evaluation."""
    path = tmp_path / "payload.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="JSON object"):
        main(
            ["validate-evidence", str(path), "--sha", SHA, "--version", VERSION, "--run-id", "41", "--run-attempt", "2"]
        )
    if payload:
        with pytest.raises(ValueError, match="array of objects"):
            main(["eligible-runs", str(path), "--sha", SHA])
