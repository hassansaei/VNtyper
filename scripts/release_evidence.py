"""Pure adapters for exact-SHA Docker release evidence."""

import argparse
import json
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

CANONICAL_DIGEST = re.compile(r"sha256:[0-9a-f]{64}")
FULL_SHA = re.compile(r"[0-9a-f]{40}")
REVISION_LABEL = "org.opencontainers.image.revision"
VERSION_LABEL = "org.opencontainers.image.version"


def _pages(payload: Mapping[str, object] | Sequence[Mapping[str, object]]) -> tuple[Mapping[str, object], ...]:
    if isinstance(payload, Mapping):
        return (payload,)
    return tuple(page for page in payload if isinstance(page, Mapping))


def _records(
    payload: Mapping[str, object] | Sequence[Mapping[str, object]],
    key: str,
) -> tuple[Mapping[str, object], ...]:
    records: list[Mapping[str, object]] = []
    for page in _pages(payload):
        values = page.get(key)
        if isinstance(values, Sequence) and not isinstance(values, (str, bytes)):
            records.extend(value for value in values if isinstance(value, Mapping))
    return tuple(records)


def eligible_successful_runs(
    payload: Mapping[str, object] | Sequence[Mapping[str, object]],
    sha: str,
) -> tuple[Mapping[str, object], ...]:
    """Select completed successful main-push Docker runs for one exact SHA.

    Args:
        payload: One workflow-runs response or the page array emitted by ``gh api --slurp``.
        sha: Exact full candidate commit SHA.

    Returns:
        Every eligible run ordered by descending numeric run ID.
    """
    eligible: list[tuple[int, Mapping[str, object]]] = []
    for run in _records(payload, "workflow_runs"):
        run_id = run.get("id")
        run_attempt = run.get("run_attempt")
        if type(run_id) is not int or type(run_attempt) is not int:
            continue
        if (
            run.get("head_sha") == sha
            and run.get("head_branch") == "main"
            and run.get("event") == "push"
            and run.get("status") == "completed"
            and run.get("conclusion") == "success"
            and isinstance(run.get("html_url"), str)
            and bool(run.get("html_url"))
        ):
            eligible.append((run_id, run))
    eligible.sort(key=lambda item: item[0], reverse=True)
    return tuple(run for _, run in eligible)


def artifact_present(
    payload: Mapping[str, object] | Sequence[Mapping[str, object]],
    name: str,
) -> bool:
    """Return whether an exact, unexpired artifact name is present.

    Args:
        payload: One artifacts response or the page array emitted by ``gh api --slurp``.
        name: Exact attempt-qualified artifact name.

    Returns:
        ``True`` only for an exact active artifact.
    """
    return any(
        artifact.get("name") == name and artifact.get("expired") is False for artifact in _records(payload, "artifacts")
    )


def exact_artifact(
    payload: Mapping[str, object] | Sequence[Mapping[str, object]],
    name: str,
) -> Mapping[str, object] | None:
    """Select the exact active artifact record including its registry URLs.

    Args:
        payload: One artifacts response or the page array emitted by ``gh api --slurp``.
        name: Exact attempt-qualified artifact name.

    Returns:
        A plain copy of the first exact active artifact, otherwise ``None``.
    """
    return next(
        (
            dict(artifact)
            for artifact in _records(payload, "artifacts")
            if artifact.get("name") == name
            and artifact.get("expired") is False
            and type(artifact.get("id")) is int
            and isinstance(artifact.get("url"), str)
            and bool(artifact.get("url"))
            and isinstance(artifact.get("archive_download_url"), str)
            and bool(artifact.get("archive_download_url"))
        ),
        None,
    )


def validate_evidence(
    payload: Mapping[str, object],
    *,
    sha: str,
    version: str,
    run_id: int,
    run_attempt: int,
) -> Mapping[str, object]:
    """Validate contract-v1 evidence against the selected run and candidate.

    Args:
        payload: Downloaded Docker evidence JSON object.
        sha: Exact candidate commit SHA.
        version: Candidate package version.
        run_id: Selected GitHub Actions run ID.
        run_attempt: Selected GitHub Actions attempt.

    Returns:
        A plain copy of the validated evidence mapping.

    Raises:
        ValueError: If any contract or identity field is invalid.
    """
    expected: tuple[tuple[str, object], ...] = (
        ("contract_version", 1),
        ("run_id", run_id),
        ("run_attempt", run_attempt),
        ("sha", sha),
        ("revision", sha),
        ("version", version),
    )
    for field, value in expected:
        if payload.get(field) != value or (isinstance(value, int) and type(payload.get(field)) is not int):
            msg = f"evidence {field} does not match the selected exact-SHA Docker Build run"
            raise ValueError(msg)
    digest = payload.get("digest")
    if not isinstance(digest, str) or CANONICAL_DIGEST.fullmatch(digest) is None:
        msg = "evidence digest must be canonical sha256:<64 lowercase hex characters>"
        raise ValueError(msg)
    return dict(payload)


def _image_labels(image: Mapping[str, object]) -> Mapping[str, object]:
    config = image.get("config")
    if not isinstance(config, Mapping):
        return {}
    labels = config.get("Labels")
    return labels if isinstance(labels, Mapping) else {}


def validate_image(
    evidence: Mapping[str, object],
    image: Mapping[str, object],
    short_image: Mapping[str, object],
    *,
    immutable_digest: str,
    short_digest: str,
) -> Mapping[str, object]:
    """Validate immutable image identity and diagnose short-prefix collisions.

    Args:
        evidence: Previously validated contract-v1 evidence.
        image: OCI image configuration inspected by immutable digest.
        short_image: OCI image configuration inspected through the short-SHA tag.
        immutable_digest: Registry-reported digest for the immutable reference.
        short_digest: Registry-reported digest for the short-SHA tag.

    Returns:
        Validated source identity plus the collision verdict.

    Raises:
        ValueError: If image identity drifts or a short-tag mismatch is not a proven prefix collision.
    """
    digest = evidence.get("digest")
    revision = evidence.get("revision")
    version = evidence.get("version")
    if not isinstance(digest, str) or CANONICAL_DIGEST.fullmatch(digest) is None:
        msg = "evidence digest must be canonical before image validation"
        raise ValueError(msg)
    if not isinstance(revision, str) or FULL_SHA.fullmatch(revision) is None:
        msg = "evidence revision must be a full lowercase hexadecimal SHA"
        raise ValueError(msg)
    if not isinstance(version, str) or not version:
        msg = "evidence version must be a non-empty string"
        raise ValueError(msg)
    if type(evidence.get("contract_version")) is not int or evidence.get("contract_version") != 1:
        msg = "evidence contract_version must be exactly 1 before image validation"
        raise ValueError(msg)
    if immutable_digest != digest:
        msg = "immutable_digest does not match the evidence digest"
        raise ValueError(msg)
    if CANONICAL_DIGEST.fullmatch(short_digest) is None:
        msg = "short_digest must be canonical sha256:<64 lowercase hex characters>"
        raise ValueError(msg)

    labels = _image_labels(image)
    if labels.get(REVISION_LABEL) != revision:
        msg = "immutable image revision label does not match evidence"
        raise ValueError(msg)
    if labels.get(VERSION_LABEL) != version:
        msg = "immutable image version label does not match evidence"
        raise ValueError(msg)

    short_labels = _image_labels(short_image)
    short_revision = short_labels.get(REVISION_LABEL)
    short_version = short_labels.get(VERSION_LABEL)
    collision = False
    if short_digest == digest:
        if short_revision != revision or short_version != version:
            msg = "unexplained short-tag drift"
            raise ValueError(msg)
    else:
        collision = (
            isinstance(short_revision, str)
            and FULL_SHA.fullmatch(short_revision) is not None
            and short_revision != revision
            and short_revision[:7] == revision[:7]
            and isinstance(short_version, str)
            and bool(short_version)
        )
        if not collision:
            msg = "unexplained short-tag drift"
            raise ValueError(msg)

    return {
        "source_digest": digest,
        "source_revision": revision,
        "source_version": version,
        "evidence_contract_version": 1,
        "short_tag_collision": collision,
    }


def _recovery_instruction(
    sha: str,
    selected_run: Mapping[str, object] | None,
    selected_artifact: Mapping[str, object] | None,
    candidate_runs: Sequence[Mapping[str, object]],
) -> str:
    if selected_run is not None:
        run_url = selected_run.get("html_url")
        run_id = selected_run.get("id")
        run_attempt = selected_run.get("run_attempt")
        artifact_name = (
            f"docker-release-evidence-{sha}-{run_attempt}"
            if selected_artifact is None
            else selected_artifact.get("name")
        )
        instruction = (
            f"Rerun existing Docker Build run {run_url} (run {run_id}, attempt {run_attempt}) to regenerate exact "
            f"artifact {artifact_name}"
        )
        if selected_artifact is None:
            return instruction + "."
        return instruction + f"; artifact API URL: {selected_artifact.get('url')}."
    candidate_urls: list[str] = []
    for run in candidate_runs:
        run_url = run.get("html_url")
        if isinstance(run_url, str):
            candidate_urls.append(run_url)
    if candidate_urls:
        return (
            "Rerun one of these existing exact-SHA Docker Build runs to regenerate its exact attempt-qualified "
            f"evidence artifact: {', '.join(candidate_urls)}."
        )
    return (
        "Wait for or rerun an exact-SHA main-push Docker Build run to regenerate its exact attempt-qualified "
        "evidence artifact."
    )


def _valid_run(record: Mapping[str, object] | None) -> bool:
    return (
        record is not None
        and type(record.get("id")) is int
        and type(record.get("run_attempt")) is int
        and isinstance(record.get("html_url"), str)
        and bool(record.get("html_url"))
    )


def _valid_artifact(record: Mapping[str, object] | None) -> bool:
    return (
        record is not None
        and type(record.get("id")) is int
        and isinstance(record.get("name"), str)
        and bool(record.get("name"))
        and isinstance(record.get("url"), str)
        and bool(record.get("url"))
        and isinstance(record.get("archive_download_url"), str)
        and bool(record.get("archive_download_url"))
    )


def _valid_evidence(record: Mapping[str, object] | None) -> bool:
    if record is None:
        return False
    digest = record.get("digest")
    revision = record.get("revision")
    version = record.get("version")
    return (
        type(record.get("contract_version")) is int
        and record.get("contract_version") == 1
        and isinstance(digest, str)
        and CANONICAL_DIGEST.fullmatch(digest) is not None
        and isinstance(revision, str)
        and FULL_SHA.fullmatch(revision) is not None
        and isinstance(version, str)
        and bool(version)
    )


def _valid_image(record: Mapping[str, object] | None) -> bool:
    labels = _image_labels(record or {})
    revision = labels.get(REVISION_LABEL)
    version = labels.get(VERSION_LABEL)
    return (
        isinstance(revision, str)
        and FULL_SHA.fullmatch(revision) is not None
        and isinstance(version, str)
        and bool(version)
    )


def _valid_source(record: Mapping[str, object] | None) -> bool:
    if record is None:
        return False
    source_ref = record.get("source_ref")
    digest = record.get("source_digest")
    revision = record.get("source_revision")
    version = record.get("source_version")
    return (
        isinstance(source_ref, str)
        and bool(source_ref)
        and isinstance(digest, str)
        and CANONICAL_DIGEST.fullmatch(digest) is not None
        and isinstance(revision, str)
        and FULL_SHA.fullmatch(revision) is not None
        and isinstance(version, str)
        and bool(version)
        and record.get("evidence_contract_version") == 1
        and type(record.get("short_tag_collision")) is bool
    )


def summarize_evidence_state(
    *,
    sha: str,
    step_outcome: str,
    selected_run: Mapping[str, object] | None,
    validated_evidence: Mapping[str, object] | None,
    image: Mapping[str, object] | None,
    candidate_runs: Sequence[Mapping[str, object]] = (),
    selected_artifact: Mapping[str, object] | None = None,
    validated_image: Mapping[str, object] | None = None,
) -> Mapping[str, object]:
    """Build an explicit evidence summary after success or partial failure.

    Args:
        sha: Exact candidate SHA.
        step_outcome: GitHub's evidence-step outcome.
        selected_run: Selected run state, when written before failure.
        validated_evidence: Validated evidence state, when available.
        image: Immutable image configuration, when inspection succeeded.
        candidate_runs: Every eligible exact-SHA run considered for evidence recovery.
        selected_artifact: Exact active artifact metadata, when selected.
        validated_image: Validated immutable source identity, when available.

    Returns:
        A JSON-compatible summary preserving every available identity field.
    """
    selected_run = selected_run if _valid_run(selected_run) else None
    selected_artifact = selected_artifact if _valid_artifact(selected_artifact) else None
    validated_evidence = validated_evidence if _valid_evidence(validated_evidence) else None
    image = image if _valid_image(image) else None
    validated_image = validated_image if _valid_source(validated_image) else None
    candidate_runs = tuple(run for run in candidate_runs if _valid_run(run))
    selected_available = selected_run is not None
    artifact_available = selected_artifact is not None
    evidence_available = validated_evidence is not None
    image_available = image is not None
    source_available = validated_image is not None
    labels = _image_labels(image or {})
    candidate_state = [
        {"id": run.get("id"), "run_attempt": run.get("run_attempt"), "html_url": run.get("html_url")}
        for run in candidate_runs
    ]
    return {
        "sha": sha,
        "state": "verified"
        if step_outcome == "success"
        and all((selected_available, artifact_available, evidence_available, image_available, source_available))
        else "failed",
        "step_outcome": step_outcome,
        "selected_run": {
            "available": selected_available,
            "id": None if selected_run is None else selected_run.get("id"),
            "run_attempt": None if selected_run is None else selected_run.get("run_attempt"),
            "html_url": None if selected_run is None else selected_run.get("html_url"),
        },
        "selected_artifact": {
            "available": artifact_available,
            "id": None if selected_artifact is None else selected_artifact.get("id"),
            "name": None if selected_artifact is None else selected_artifact.get("name"),
            "url": None if selected_artifact is None else selected_artifact.get("url"),
            "archive_download_url": (
                None if selected_artifact is None else selected_artifact.get("archive_download_url")
            ),
        },
        "candidate_runs": candidate_state,
        "evidence": {
            "available": evidence_available,
            "contract_version": (None if validated_evidence is None else validated_evidence.get("contract_version")),
            "digest": None if validated_evidence is None else validated_evidence.get("digest"),
            "revision": None if validated_evidence is None else validated_evidence.get("revision"),
            "version": None if validated_evidence is None else validated_evidence.get("version"),
        },
        "image": {
            "available": image_available,
            "revision": labels.get(REVISION_LABEL),
            "version": labels.get(VERSION_LABEL),
        },
        "source_ref": None if validated_image is None else validated_image.get("source_ref"),
        "evidence_contract_version": (
            None if validated_evidence is None else validated_evidence.get("contract_version")
        ),
        "short_tag_collision": (None if validated_image is None else validated_image.get("short_tag_collision")),
        "recovery_instruction": _recovery_instruction(sha, selected_run, selected_artifact, candidate_runs),
    }


def _load_mapping(path: Path) -> Mapping[str, object]:
    payload: object = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        msg = f"{path} must contain a JSON object"
        raise ValueError(msg)
    return payload


def _load_pages(path: Path) -> Mapping[str, object] | Sequence[Mapping[str, object]]:
    payload: object = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, Mapping):
        return payload
    if (
        isinstance(payload, Sequence)
        and not isinstance(payload, (str, bytes))
        and all(isinstance(page, Mapping) for page in payload)
    ):
        return payload  # type: ignore[return-value]
    msg = f"{path} must contain a JSON object or array of objects"
    raise ValueError(msg)


def _optional_mapping(path: Path) -> Mapping[str, object] | None:
    if not path.is_file():
        return None
    try:
        return _load_mapping(path)
    except (OSError, ValueError):
        return None


def _optional_pages(path: Path) -> tuple[Mapping[str, object], ...]:
    if not path.is_file():
        return ()
    try:
        payload = _load_pages(path)
    except (OSError, ValueError):
        return ()
    if isinstance(payload, Mapping):
        return (payload,)
    return tuple(payload)


def _write_json_atomic(path: Path, payload: object) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    eligible = commands.add_parser("eligible-runs")
    eligible.add_argument("payload", type=Path)
    eligible.add_argument("--sha", required=True)

    artifact = commands.add_parser("artifact-present")
    artifact.add_argument("payload", type=Path)
    artifact.add_argument("--name", required=True)

    selected_artifact = commands.add_parser("select-artifact")
    selected_artifact.add_argument("payload", type=Path)
    selected_artifact.add_argument("--name", required=True)

    evidence = commands.add_parser("validate-evidence")
    evidence.add_argument("payload", type=Path)
    evidence.add_argument("--sha", required=True)
    evidence.add_argument("--version", required=True)
    evidence.add_argument("--run-id", required=True, type=int)
    evidence.add_argument("--run-attempt", required=True, type=int)

    image = commands.add_parser("validate-image")
    image.add_argument("evidence", type=Path)
    image.add_argument("image_config", type=Path)
    image.add_argument("short_image_config", type=Path)
    image.add_argument("--immutable-digest", required=True)
    image.add_argument("--short-digest", required=True)
    image.add_argument("--image", required=True)
    image.add_argument("--selected-run", required=True, type=Path)
    image.add_argument("--selected-artifact", required=True, type=Path)
    image.add_argument("--state-output", required=True, type=Path)
    image.add_argument("--output", required=True, type=Path)

    summary = commands.add_parser("summarize")
    summary.add_argument("--sha", required=True)
    summary.add_argument("--step-outcome", required=True)
    summary.add_argument("--selected-run", required=True, type=Path)
    summary.add_argument("--validated-evidence", required=True, type=Path)
    summary.add_argument("--image", required=True, type=Path)
    summary.add_argument("--candidate-runs", type=Path, default=Path("evidence-eligible-runs.json"))
    summary.add_argument("--selected-artifact", type=Path, default=Path("selected-artifact.json"))
    summary.add_argument("--validated-image", type=Path, default=Path("validated-image.json"))
    summary.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the local-file release evidence adapter.

    Args:
        argv: Optional CLI arguments excluding the program name.

    Returns:
        Process exit status.
    """
    args = _parser().parse_args(argv)
    if args.command == "eligible-runs":
        print(json.dumps(eligible_successful_runs(_load_pages(args.payload), args.sha), separators=(",", ":")))
        return 0
    if args.command == "artifact-present":
        return 0 if artifact_present(_load_pages(args.payload), args.name) else 1
    if args.command == "select-artifact":
        selected = exact_artifact(_load_pages(args.payload), args.name)
        if selected is None:
            return 1
        print(json.dumps(selected, separators=(",", ":")))
        return 0
    if args.command == "validate-evidence":
        result = validate_evidence(
            _load_mapping(args.payload),
            sha=args.sha,
            version=args.version,
            run_id=args.run_id,
            run_attempt=args.run_attempt,
        )
        print(json.dumps(result, separators=(",", ":")))
        return 0
    if args.command == "validate-image":
        evidence_payload = _load_mapping(args.evidence)
        result = validate_image(
            evidence_payload,
            _load_mapping(args.image_config),
            _load_mapping(args.short_image_config),
            immutable_digest=args.immutable_digest,
            short_digest=args.short_digest,
        )
        selected_run = _load_mapping(args.selected_run)
        selected_artifact = _load_mapping(args.selected_artifact)
        source_ref = f"{args.image}@{result['source_digest']}"
        validated_image = {"source_ref": source_ref, **result}
        recovery_instruction = _recovery_instruction(
            str(evidence_payload.get("sha")), selected_run, selected_artifact, ()
        )
        lines = (
            f"source_ref={source_ref}",
            f"source_digest={result['source_digest']}",
            f"source_revision={result['source_revision']}",
            f"source_version={result['source_version']}",
            f"evidence_contract_version={result['evidence_contract_version']}",
            f"short_tag_collision={str(result['short_tag_collision']).lower()}",
            f"source_run_id={selected_run.get('id')}",
            f"source_run_attempt={selected_run.get('run_attempt')}",
            f"source_run_url={selected_run.get('html_url')}",
            f"source_artifact_name={selected_artifact.get('name')}",
            f"source_artifact_url={selected_artifact.get('url')}",
            f"source_artifact_download_url={selected_artifact.get('archive_download_url')}",
            f"recovery_instruction={recovery_instruction}",
        )
        _write_json_atomic(args.state_output, validated_image)
        args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return 0
    if args.command == "summarize":
        summary = summarize_evidence_state(
            sha=args.sha,
            step_outcome=args.step_outcome,
            selected_run=_optional_mapping(args.selected_run),
            validated_evidence=_optional_mapping(args.validated_evidence),
            image=_optional_mapping(args.image),
            candidate_runs=_optional_pages(args.candidate_runs),
            selected_artifact=_optional_mapping(args.selected_artifact),
            validated_image=_optional_mapping(args.validated_image),
        )
        _write_json_atomic(args.output, summary)
        return 0
    raise AssertionError("argparse accepted an unknown release evidence command")


if __name__ == "__main__":
    raise SystemExit(main())
