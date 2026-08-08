"""Safe persistence of curated alignment-preflight failures."""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts.preflight_error_io import (
    public_preflight_error_payload,
    public_reference_error_payload,
    write_preflight_error,
)

pytestmark = pytest.mark.unit


def _payload() -> dict:
    """Return one literal payload for writer tests."""
    return {
        "code": "reference_unresolved",
        "message": "Unable to resolve CRAM reference: contig=chr1, M5=digest.",
        "candidates": [["cli", "full.fa", "probe exited non-zero"]],
    }


def test_public_reference_payload_has_exact_keys_and_no_absolute_worker_path() -> None:
    """Raw samtools diagnostics and worker roots must never reach the public message."""
    payload = public_reference_error_payload(
        "chr1",
        "digest",
        (
            ("cli", "/opt/vntyper/private/full.fa", "decode failed for /opt/vntyper/private/input.cram"),
            ("config_cram_reference", None, "not supplied"),
            ("config_bwa_reference", "/refs/chr1.fa", "reference FASTA not found"),
        ),
    )

    assert set(payload) == {"code", "message", "candidates"}
    assert payload["code"] == "reference_unresolved"
    assert "contig=chr1" in payload["message"]
    assert "M5=digest" in payload["message"]
    assert [candidate[0] for candidate in payload["candidates"]] == [
        "cli",
        "config_cram_reference",
        "config_bwa_reference",
    ]
    serialized = json.dumps(payload)
    assert "/opt/vntyper" not in serialized
    assert "/refs" not in serialized
    assert "decode failed for" not in serialized


def test_public_reference_payload_sanitizes_hostile_header_identifiers() -> None:
    """Contig and checksum fields stay useful without transporting path-like header text."""
    payload = public_reference_error_payload(
        "/private/worker/contig",
        "digest\\private\\worker",
        (("htslib-resolved", None, "probe exited non-zero"),),
    )

    assert (
        payload["message"] == "Unable to resolve reference for CRAM input: contig=unknown, M5=unknown. Candidates: "
        "source=htslib-resolved, path=no path, reason=probe exited non-zero"
    )
    assert "/private/" not in json.dumps(payload)
    assert "\\private\\" not in json.dumps(payload)


@pytest.mark.parametrize(
    ("internal_message", "code", "public_message"),
    [
        (
            "Alignment view must stay inside output directory: /private/view.bam",
            "preflight_output_unsafe",
            "Alignment preflight rejected an unsafe alignment view; choose a safe output directory.",
        ),
        (
            "Log path /private/sample.log aliases protected source /private/patient.bam",
            "preflight_output_unsafe",
            "Alignment preflight rejected an unsafe log destination; remove the conflicting entry and retry.",
        ),
        (
            "Generated index /private/sample.bai resolves to the same file as protected source /private/input.bai",
            "preflight_output_unsafe",
            "Alignment preflight rejected an unsafe index destination or ownership record; "
            "remove the conflicting entry and retry.",
        ),
        (
            "cram.local_ref_path must not contain a remote URL when ambient reference resolution is disabled: "
            "https://private.invalid/ref",
            "reference_policy_invalid",
            "CRAM reference policy rejected a remote REF_PATH; configure a local cache or explicitly allow "
            "ambient resolution.",
        ),
    ],
)
def test_public_failure_curator_preserves_the_safe_actionable_failure_category(
    internal_message: str,
    code: str,
    public_message: str,
) -> None:
    """Unsafe destinations remain actionable without copying their private paths."""
    payload = public_preflight_error_payload(ValueError(internal_message))

    assert payload == {"code": code, "message": public_message, "candidates": []}
    assert "/private/" not in json.dumps(payload)


def test_writer_atomically_replaces_a_destination_symlink_without_following_it(tmp_path: Path) -> None:
    """A planted link is replaced as an entry; its external target stays unchanged."""
    protected = tmp_path / "worker-private.txt"
    protected.write_text("do not overwrite", encoding="utf-8")
    destination = tmp_path / "preflight_error.json"
    destination.symlink_to(protected)

    written = write_preflight_error(tmp_path, _payload())

    assert written == destination
    assert not destination.is_symlink()
    assert json.loads(destination.read_text(encoding="utf-8")) == _payload()
    assert protected.read_text(encoding="utf-8") == "do not overwrite"


def test_writer_atomically_replaces_an_existing_single_link_regular_artifact(tmp_path: Path) -> None:
    """A prior owned artifact is replaced completely rather than rewritten in place."""
    destination = tmp_path / "preflight_error.json"
    destination.write_text('{"previous":true}\n', encoding="utf-8")

    write_preflight_error(tmp_path, _payload())

    assert json.loads(destination.read_text(encoding="utf-8")) == _payload()
    assert destination.stat().st_nlink == 1


def test_writer_rejects_a_payload_with_extra_or_missing_contract_keys(tmp_path: Path) -> None:
    """Only the reviewed three-field artifact may cross the worker boundary."""
    with pytest.raises(ValueError, match="exactly code, message, candidates"):
        write_preflight_error(tmp_path, {"code": "reference_unresolved", "message": "unsafe"})


def test_failed_atomic_install_preserves_the_previous_artifact_and_cleans_the_temporary_file(
    tmp_path: Path,
) -> None:
    """An interrupted install cannot expose partial JSON or leave patient-adjacent scratch."""
    destination = tmp_path / "preflight_error.json"
    destination.write_text('{"previous": true}\n', encoding="utf-8")

    with (
        mock.patch("vntyper.scripts.preflight_error_io.os.replace", side_effect=OSError("disk full")),
        pytest.raises(OSError, match="disk full"),
    ):
        write_preflight_error(tmp_path, _payload())

    assert destination.read_text(encoding="utf-8") == '{"previous": true}\n'
    assert set(tmp_path.iterdir()) == {destination}


@pytest.mark.parametrize("entry_kind", ["directory", "fifo", "hardlink"])
def test_writer_rejects_wrong_type_or_multiply_linked_destination_without_mutation(
    tmp_path: Path, entry_kind: str
) -> None:
    """Only an unlinked name, symlink entry, or single-link regular file is replaceable."""
    destination = tmp_path / "preflight_error.json"
    protected = tmp_path / "protected.json"
    protected.write_text("protected bytes", encoding="utf-8")
    if entry_kind == "directory":
        destination.mkdir()
    elif entry_kind == "fifo":
        os.mkfifo(destination)
    else:
        os.link(protected, destination)

    with pytest.raises(ValueError, match="safe final entry"):
        write_preflight_error(tmp_path, _payload())

    assert protected.read_text(encoding="utf-8") == "protected bytes"
    assert not tuple(tmp_path.glob(".preflight_error.json.*.tmp"))


def test_writer_rejects_a_symlinked_output_directory_without_touching_its_target(tmp_path: Path) -> None:
    """The fixed artifact name stays in the literal run directory supplied by the caller."""
    protected = tmp_path / "protected"
    protected.mkdir()
    output = tmp_path / "output"
    output.symlink_to(protected, target_is_directory=True)

    with pytest.raises(ValueError, match="real directory"):
        write_preflight_error(output, _payload())

    assert not (protected / "preflight_error.json").exists()
