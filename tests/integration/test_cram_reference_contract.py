"""A-209 reference-dependent and reference-free CRAM pipeline contracts."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

import pysam
import pytest

from scripts.make_cram_fixtures import build_reference_dependent_fixture

logger = logging.getLogger(__name__)

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_DATA_CONFIG = json.loads((REPO_ROOT / "tests" / "test_data_config.json").read_text(encoding="utf-8"))
PURPOSE_FIXTURES = TEST_DATA_CONFIG["purpose_fixtures"]["cram_reference_contract"]


def _require_purpose_fixture(name: str) -> Path:
    """Resolve one registered fixture or name the ordinary derivation command."""
    path = REPO_ROOT / PURPOSE_FIXTURES[name]
    assert path.is_file(), f"Missing purpose fixture {path}; run `make cram-fixtures`."
    return path


def _local_config(tmp_path: Path) -> Path:
    """Pin reference lookup and coverage to the purpose fixture's compact contig."""
    payload = json.loads((REPO_ROOT / "vntyper" / "config.json").read_text(encoding="utf-8"))
    payload["cram"]["local_ref_path"] = str(tmp_path / "local-ref" / "%2s" / "%2s" / "%s")
    payload["bam_processing"]["assemblies"]["GRCh37"]["vntr_region_coords"] = "1-10000"
    config = tmp_path / "config.json"
    config.write_text(json.dumps(payload), encoding="utf-8")
    return config


def _no_ref_config(tmp_path: Path) -> Path:
    """Remove every explicit hg19 reference while keeping terminal resolution last."""
    config = _local_config(tmp_path)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["reference_data"]["cram_reference_hg19"] = None
    payload["reference_data"]["bwa_reference_hg19"] = None
    payload["cram"]["reference_candidate_order"] = [
        "cli",
        "config_cram_reference",
        "config_bwa_reference",
        "htslib_resolved",
    ]
    config.write_text(json.dumps(payload), encoding="utf-8")
    return config


def _tree_digest(root: Path) -> dict[str, str]:
    """Return the relative listing and bytes of every regular file in a tree."""
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _run_cram(
    cram: Path,
    output_dir: Path,
    config: Path,
    *,
    reference: Path | None = None,
    fast_mode: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run the purpose CRAM over its own compact target interval."""
    command = [
        sys.executable,
        "-m",
        "vntyper.cli",
        "--config-path",
        str(config),
        "-l",
        "DEBUG",
        "pipeline",
        "--cram",
        str(cram),
        "--reference-assembly",
        "hg19",
        "--custom-regions",
        "chr1:1-10000",
        "--threads",
        "2",
        "--keep-intermediates",
        "--output-dir",
        str(output_dir),
    ]
    if fast_mode:
        command.append("--fast-mode")
    if reference is not None:
        command.extend(["--reference-fasta", str(reference)])
    environment = os.environ.copy()
    environment["NO_PROXY"] = "127.0.0.1,localhost"
    environment["no_proxy"] = "127.0.0.1,localhost"
    result = subprocess.run(command, capture_output=True, text=True, check=False, env=environment)
    logger.info("CRAM contract stdout:\n%s", result.stdout)
    logger.info("CRAM contract stderr:\n%s", result.stderr)
    return result


def _reheader_with_reference_uri(
    source: Path,
    destination: Path,
    uri: str,
    *,
    trailing_uri: str | None = None,
) -> None:
    """Rewrite the first SQ UR without decoding the reference-compressed records."""
    header = subprocess.check_output(["samtools", "view", "-H", str(source)], text=True)
    rewritten: list[str] = []
    changed = False
    for line in header.splitlines():
        if line.startswith("@SQ\t") and not changed:
            fields = line.split("\t")
            fields = [f"UR:{uri}" if field.startswith("UR:") else field for field in fields]
            if trailing_uri is not None:
                fields.append(f"UR:{trailing_uri}")
            changed = any(field == f"UR:{uri}" for field in fields)
            line = "\t".join(fields)
        rewritten.append(line)
    assert changed, "The purpose CRAM header has no SQ UR tag to rewrite."
    header_path = destination.with_suffix(".header.sam")
    header_path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    with destination.open("wb") as output_handle:
        subprocess.run(["samtools", "reheader", str(header_path), str(source)], stdout=output_handle, check=True)


def _reheader_without_reference_uri(source: Path, destination: Path) -> None:
    """Remove SQ UR tags so the terminal probe can exercise REF_PATH alone."""
    header = subprocess.check_output(["samtools", "view", "-H", str(source)], text=True)
    rewritten = [
        "\t".join(field for field in line.split("\t") if not field.startswith("UR:"))
        if line.startswith("@SQ\t")
        else line
        for line in header.splitlines()
    ]
    header_path = destination.with_suffix(".header.sam")
    header_path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    with destination.open("wb") as output_handle:
        subprocess.run(["samtools", "reheader", str(header_path), str(source)], stdout=output_handle, check=True)


def _assert_no_reference_flags_in_cram_commands(result: subprocess.CompletedProcess[str]) -> None:
    """Require every logged CRAM consumer to use the terminal no-reference form."""
    markers = (
        "Running captured command: samtools view ",
        "Executing region slicing with command: samtools view ",
        "Executing filtering with command: set -o pipefail; samtools view ",
        "Calculating VNTR coverage with command: samtools depth ",
    )
    for marker in markers:
        commands = [line for line in result.stderr.splitlines() if marker in line]
        assert commands, f"No logged CRAM command matched {marker!r}"
        assert all(" -T " not in command and " --reference " not in command for command in commands), commands


@pytest.mark.parametrize("read_only", [False, True], ids=["writable", "read-only"])
def test_local_header_ur_reference_is_bound_without_mutating_its_operator_tree(
    tmp_path: Path,
    ensure_test_data: None,
    read_only: bool,
) -> None:
    """A header-only unindexed FASTA succeeds without an operator-side FAI."""
    fixture = build_reference_dependent_fixture(tmp_path / "patient-input")
    reference_index = Path(f"{fixture.reference}.fai")
    reference_index.unlink()
    operator_root = fixture.reference.parent
    before = _tree_digest(operator_root)
    if read_only:
        fixture.reference.chmod(0o444)
        operator_root.chmod(0o555)

    try:
        result = _run_cram(
            fixture.cram,
            tmp_path / "run-output",
            _no_ref_config(tmp_path),
        )
    finally:
        if read_only:
            operator_root.chmod(0o755)
            fixture.reference.chmod(0o644)

    assert result.returncode == 0, result.stdout + result.stderr
    assert _tree_digest(operator_root) == before
    assert "Resolved CRAM reference from header_ur" in result.stderr


def test_a209_1_missing_reference_names_the_digest_and_candidates_before_stages(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """The actual header ``UR:`` target is renamed before the failing decode."""
    _require_purpose_fixture("reference_dependent_cram")
    _require_purpose_fixture("reference_fasta")
    fixture = build_reference_dependent_fixture(tmp_path / "purpose")
    with pysam.AlignmentFile(str(fixture.cram), "rc", reference_filename=str(fixture.reference)) as alignment:
        sequence = alignment.header.to_dict()["SQ"][0]
    ur_target = Path(sequence["UR"])
    expected_m5 = hashlib.md5(("A" * 10_000).encode()).hexdigest()
    assert ur_target.samefile(fixture.reference)
    assert sequence["M5"] == expected_m5

    hidden_target = ur_target.with_name(f"{ur_target.name}.a209-missing")
    ur_target.rename(hidden_target)
    output = tmp_path / "missing-output"
    try:
        config = _local_config(tmp_path)
        result = _run_cram(fixture.cram, output, config)
    finally:
        hidden_target.rename(ur_target)

    diagnostic = f"{result.stdout}\n{result.stderr}"
    assert result.returncode != 0
    assert "contig=chr1" in diagnostic
    assert f"M5={expected_m5}" in diagnostic
    for source in (
        "source=cli",
        "source=config_cram_reference",
        "source=config_bwa_reference",
        "source=header_ur",
    ):
        assert source in diagnostic
    assert "source=htslib-resolved" not in diagnostic
    assert not (output / "kestrel").exists()
    assert not (output / "coverage").exists()
    assert not (output / "pipeline_summary.json").exists()


def test_a209_2_explicit_reference_completes_the_reference_dependent_cram(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    cram = _require_purpose_fixture("reference_dependent_cram")
    reference = _require_purpose_fixture("reference_fasta")
    output = tmp_path / "explicit-reference-output"

    result = _run_cram(cram, output, _local_config(tmp_path), reference=reference)

    assert result.returncode == 0
    assert (output / "pipeline_summary.json").is_file()
    assert (output / "kestrel" / "kestrel_result.tsv").is_file()


def test_a209_3_no_ref_cram_completes_without_an_explicit_reference(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    cram = _require_purpose_fixture("no_ref_cram")
    output = tmp_path / "no-ref-output"

    result = _run_cram(cram, output, _no_ref_config(tmp_path), fast_mode=False)

    assert result.returncode == 0
    assert "Resolved CRAM reference through htslib-resolved" in result.stderr
    _assert_no_reference_flags_in_cram_commands(result)
    assert (output / "pipeline_summary.json").is_file()
    assert (output / "kestrel" / "kestrel_result.tsv").is_file()


def test_a178_1_blackhole_reference_probe_exits_within_its_configured_deadline(
    tmp_path: Path,
    ensure_test_data: None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A server that accepts and never responds cannot hold a CRAM run forever."""
    fixture = build_reference_dependent_fixture(tmp_path / "purpose")
    hidden_reference = fixture.reference.with_name("hidden-reference.fa")
    fixture.reference.rename(hidden_reference)

    listener = socket.socket()
    listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listener.bind(("127.0.0.1", 0))
    listener.listen()
    listener.settimeout(0.1)
    port = listener.getsockname()[1]
    accepted = threading.Event()
    stop = threading.Event()

    def blackhole() -> None:
        connection: socket.socket | None = None
        try:
            while not stop.is_set():
                try:
                    connection, _ = listener.accept()
                    break
                except TimeoutError:
                    continue
            if connection is not None:
                accepted.set()
                stop.wait(10)
        finally:
            if connection is not None:
                connection.close()

    server = threading.Thread(target=blackhole, daemon=True)
    server.start()
    ref_path_cram = fixture.cram.with_name("ref-path.cram")
    _reheader_without_reference_uri(fixture.cram, ref_path_cram)
    monkeypatch.setenv("REF_PATH", f"http://127.0.0.1:{port}/%s")
    monkeypatch.setenv("NO_PROXY", "127.0.0.1,localhost")
    config = _no_ref_config(tmp_path)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["cram"]["allow_ambient_reference_resolution"] = True
    payload["cram"]["reference_probe_timeout_seconds"] = 0.25
    config.write_text(json.dumps(payload), encoding="utf-8")

    started = time.monotonic()
    try:
        result = _run_cram(ref_path_cram, tmp_path / "blackhole-output", config)
        accepted.wait(timeout=0.25)
    finally:
        elapsed = time.monotonic() - started
        stop.set()
        listener.close()
        server.join(timeout=2)
        hidden_reference.rename(fixture.reference)

    diagnostic = f"{result.stdout}\n{result.stderr}"
    assert accepted.is_set(), "the CRAM decode never reached the blackhole reference server"
    assert result.returncode != 0
    assert elapsed < 5
    assert "timed out after 0.25 seconds" in diagnostic


def test_default_mode_rejects_a_remote_header_uri_before_network_or_probe_work(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """A default CRAM run refuses a localhost SQ UR before any captured command."""
    input_root = tmp_path / "patient-input"
    fixture = build_reference_dependent_fixture(input_root / "purpose")
    listener = socket.socket()
    listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listener.bind(("127.0.0.1", 0))
    listener.listen()
    listener.settimeout(0.1)
    port = listener.getsockname()[1]
    accepted = threading.Event()
    stop = threading.Event()

    def respond_without_blocking() -> None:
        while not stop.is_set():
            try:
                connection, _ = listener.accept()
            except TimeoutError:
                continue
            except OSError:
                break
            accepted.set()
            with connection:
                try:
                    connection.recv(4096)
                    connection.sendall(b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\n\r\n")
                except OSError:
                    pass

    server = threading.Thread(target=respond_without_blocking, daemon=True)
    server.start()
    remote_uri = f"http://127.0.0.1:{port}/private/reference.fa"
    remote_cram = input_root / "remote-header.cram"
    _reheader_with_reference_uri(
        fixture.cram,
        remote_cram,
        remote_uri,
        trailing_uri="file:///local/reference.fa",
    )
    rewritten_header = subprocess.check_output(["samtools", "view", "-H", str(remote_cram)], text=True)
    assert f"UR:{remote_uri}\tUR:file:///local/reference.fa" in rewritten_header
    output = tmp_path / "run-output" / "remote-header-output"
    try:
        result = _run_cram(remote_cram, output, _local_config(tmp_path))
        accepted.wait(timeout=0.25)
    finally:
        stop.set()
        listener.close()
        server.join(timeout=2)

    diagnostic = f"{result.stdout}\n{result.stderr}"
    assert result.returncode != 0
    assert "contig=chr1" in diagnostic
    assert "scheme=http" in diagnostic
    assert remote_uri not in diagnostic
    assert not accepted.is_set(), "default preflight contacted the remote SQ UR"
    assert "Running captured command:" not in diagnostic
    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8")) == {
        "code": "reference_policy_invalid",
        "message": "Remote CRAM header reference is disabled by policy: contig=chr1, scheme=http. Replace the "
        "@SQ UR with a local path, relative path, or file-scheme reference, or set "
        "cram.allow_ambient_reference_resolution=true to accept network access.",
        "candidates": [],
    }
    assert not (output / "custom_regions.bed").exists()
    preflight_dir = output / "fastq_bam_processing"
    assert preflight_dir.is_dir()
    assert list(preflight_dir.iterdir()) == []
    assert not (output / "kestrel").exists()
    assert not (output / "coverage").exists()
    assert not (output / "pipeline_summary.json").exists()
