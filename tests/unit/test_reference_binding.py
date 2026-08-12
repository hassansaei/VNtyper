"""Lifetime, replacement, and cleanup tests for CRAM reference bindings."""

from __future__ import annotations

import hashlib
import os
import shlex
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.alignment_preflight import resolve_reference, run_preflight
from vntyper.scripts.reference_binding import ReferenceBinding
from vntyper.scripts.reference_uri_policy import LocalHeaderReference

pytestmark = pytest.mark.unit


def _reference(path: Path) -> Path:
    """Create a small reference FASTA without sidecars."""
    path.write_text(">chr1\nACGT\n", encoding="utf-8")
    return path


def test_reference_binding_preserves_a_colliding_run_local_entry(tmp_path: Path) -> None:
    """A raced-in entry is never removed or overwritten while installing a view."""
    reference = _reference(tmp_path / "reference.fa")
    collided_paths: list[Path] = []

    def collide(_source: str, destination: str | Path, *_args: object, **_kwargs: object) -> None:
        collision = Path(destination)
        collision.write_bytes(b"not pipeline-owned")
        collided_paths.append(collision)
        raise FileExistsError("simulated path collision")

    with (
        patch("vntyper.scripts.reference_binding.os.symlink", side_effect=collide),
        pytest.raises(RuntimeError, match="Unable to retain CRAM reference input"),
    ):
        ReferenceBinding(str(reference), tmp_path, "sample", 1)

    assert len(collided_paths) == 1
    assert collided_paths[0].read_bytes() == b"not pipeline-owned"


def test_constructor_cleanup_failure_does_not_replace_a_primary_base_exception(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Cleanup diagnostics cannot replace cancellation or interpreter shutdown."""
    reference = _reference(tmp_path / "reference.fa")

    with (
        patch("vntyper.scripts.reference_binding._InodeView", side_effect=KeyboardInterrupt("cancelled")),
        patch.object(
            ReferenceBinding,
            "_close_after_failed_construction",
            side_effect=RuntimeError("cleanup exploded"),
        ),
        caplog.at_level("ERROR"),
        pytest.raises(KeyboardInterrupt, match="cancelled"),
    ):
        ReferenceBinding(str(reference), tmp_path, "sample", 1)

    assert "incomplete CRAM reference namespace cleanup: cleanup exploded" in caplog.text


def test_binding_collision_is_fatal_before_an_ambient_probe_can_succeed(tmp_path: Path) -> None:
    """Infrastructure failure cannot become a successful unbound reference plan."""
    reference = _reference(tmp_path / "reference.fa")
    collisions: list[Path] = []

    def collide(_source: str, destination: str | Path, *_args: object, **_kwargs: object) -> None:
        collision = Path(destination)
        collision.write_bytes(b"foreign archive payload")
        collisions.append(collision)
        raise FileExistsError("simulated path collision")

    with (
        patch("vntyper.scripts.reference_binding.os.symlink", side_effect=collide),
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "ambient decoded")) as capture,
        pytest.raises(RuntimeError, match="Unable to retain CRAM reference input"),
    ):
        resolve_reference(
            "/run/view.cram",
            (("header_ur", str(reference)),),
            "chr1:1-4",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
        )

    capture.assert_not_called()
    assert collisions == [tmp_path / ".sample_reference_1" / "reference.fa"]
    assert collisions[0].read_bytes() == b"foreign archive payload"
    assert not (tmp_path / "results.zip").exists()

    with pytest.raises(RuntimeError, match="reserved CRAM reference directory already exists"):
        ReferenceBinding(str(reference), tmp_path, "sample", 1)


def test_failed_local_header_ur_decode_cannot_fall_through_to_no_reference_probe(tmp_path: Path) -> None:
    """A no-T retry must not reopen or index a local header reference path."""
    reference = _reference(tmp_path / "wrong-reference.fa")
    outcomes = ((False, "reference mismatch"), (True, "unsafe ambient success"))

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=outcomes) as capture,
        pytest.raises(ValueError, match="Unable to resolve reference"),
    ):
        resolve_reference(
            "/run/view.cram",
            (("header_ur", str(reference)),),
            "chr1:1-4",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
        )

    assert capture.call_count == 1
    assert " -T " in f" {capture.call_args.args[0]} "


def test_ambient_remote_ref_path_cannot_reopen_a_writable_failed_local_header_uri(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The network waiver does not authorize an unbound local-UR fallback."""
    (tmp_path / "operator").mkdir()
    reference = _reference(tmp_path / "operator" / "wrong-reference.fa")
    before = {path.name: path.read_bytes() for path in reference.parent.iterdir()}
    required_digest = hashlib.md5(b"AAAA").hexdigest()
    monkeypatch.setenv("REF_PATH", "https://refget.example/%s")

    with (
        patch(
            "vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "reference mismatch")
        ) as capture,
        pytest.raises(ValueError, match="Unable to resolve reference"),
    ):
        resolve_reference(
            "/run/view.cram",
            (("header_ur", str(reference)),),
            "chr1:1-4",
            None,
            {"cram": {"allow_ambient_reference_resolution": True}},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            required_digest,
            header_m5s=(("chr1", required_digest),),
            header_references=(LocalHeaderReference("chr1", required_digest, str(reference)),),
        )

    assert capture.call_count == 1
    assert {path.name: path.read_bytes() for path in reference.parent.iterdir()} == before


def test_preflight_probes_and_retains_the_opened_reference_and_fai_inodes(tmp_path: Path) -> None:
    """Replacing operator paths during the probe cannot switch reference bytes."""
    alignment = tmp_path / "input.cram"
    alignment.write_bytes(b"CRAM input")
    reference = _reference(tmp_path / "reference.fa")
    fai = Path(f"{reference}.fai")
    fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    gzi = Path(f"{reference}.gzi")
    gzi.write_bytes(b"original GZI")
    replacement = tmp_path / "replacement.fa"
    replacement.write_bytes(b">chr1\nTGCA\n")
    replacement_fai = tmp_path / "replacement.fa.fai"
    replacement_fai.write_text("chr1\t4\t999\t4\t5\n", encoding="utf-8")
    replacement_gzi = tmp_path / "replacement.fa.gzi"
    replacement_gzi.write_bytes(b"replacement GZI")
    output = tmp_path / "run"
    observed: list[tuple[bytes, bytes, bytes]] = []

    def capture(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        if "index" in arguments:
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"fresh CRAI")
            return True, ""
        if "idxstats" in arguments:
            return True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"
        if "-c" in arguments:
            return True, "2\n"
        if "-T" in arguments:
            candidate = Path(arguments[arguments.index("-T") + 1])
            replacement.replace(reference)
            replacement_fai.replace(fai)
            replacement_gzi.replace(gzi)
            observed.append(
                (
                    candidate.read_bytes(),
                    Path(f"{candidate}.fai").read_bytes(),
                    Path(f"{candidate}.gzi").read_bytes(),
                )
            )
        return True, "decoded"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=capture):
        plan = run_preflight(
            str(alignment),
            str(output),
            "input",
            "cram",
            {
                "tools": {"samtools": "samtools"},
                "cram": {"reference_candidate_order": ["cli", "htslib_resolved"]},
            },
            1,
            region="chr1:1-4",
            reference_fasta=str(reference),
            header_contigs=("chr1",),
        )

    assert plan.reference_path is not None
    reference_view = Path(plan.reference_path)
    reference_index_view = Path(f"{plan.reference_path}.fai")
    compressed_index_view = Path(f"{plan.reference_path}.gzi")
    try:
        assert observed == [(b">chr1\nACGT\n", b"chr1\t4\t6\t4\t5\n", b"original GZI")]
        assert reference_view != reference
        assert reference_view.read_bytes() == b">chr1\nACGT\n"
        assert reference_index_view.read_bytes() == b"chr1\t4\t6\t4\t5\n"
        assert compressed_index_view.read_bytes() == b"original GZI"
        assert reference.read_bytes() == b">chr1\nTGCA\n"
        assert fai.read_text(encoding="utf-8") == "chr1\t4\t999\t4\t5\n"
        assert gzi.read_bytes() == b"replacement GZI"
    finally:
        plan.close()

    assert not os.path.lexists(reference_view)
    assert not os.path.lexists(reference_index_view)
    assert not os.path.lexists(compressed_index_view)


def test_generated_fai_is_bound_before_the_coverage_probe_reopens_it(tmp_path: Path) -> None:
    """Deferring FAI binding until all probes finish leaves the depth probe unbound.

    The coverage probe must see the retained inode under its own name: one link, not a
    replaced entry pointing back at this process's descriptor for it (#238).
    """
    reference = _reference(tmp_path / "reference.fa")
    observed_bound_state: list[tuple[int, bytes]] = []
    order: list[str] = []
    unbound_binder = ReferenceBinding.bind_generated_sidecars

    def spy_bind(binding: ReferenceBinding) -> None:
        order.append("bind")
        unbound_binder(binding)

    def capture(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        if "-T" in arguments:
            order.append("target-probe")
            candidate = Path(arguments[arguments.index("-T") + 1])
            Path(f"{candidate}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
        if "--reference" in arguments:
            order.append("coverage-probe")
            candidate = Path(arguments[arguments.index("--reference") + 1])
            index = Path(f"{candidate}.fai")
            observed_bound_state.append((index.stat(follow_symlinks=False).st_nlink, index.read_bytes()))
        return True, "decoded"

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=capture),
        patch.object(ReferenceBinding, "bind_generated_sidecars", spy_bind),
    ):
        resolved, source, uncovered, binding = resolve_reference(
            "/run/view.cram",
            (("cli", str(reference)),),
            "chr1:1-4",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
            coverage_region="chr1:1-4",
        )

    try:
        assert resolved is not None
        assert (source, uncovered) == ("cli", ())
        assert order == ["target-probe", "bind", "coverage-probe"]
        assert observed_bound_state == [(1, b"chr1\t4\t6\t4\t5\n")]
    finally:
        assert binding is not None
        binding.close()


def test_probe_failure_is_not_replaced_by_reference_cleanup_failure(tmp_path: Path) -> None:
    """Cleanup diagnostics cannot hide the command boundary's primary exception."""
    reference = _reference(tmp_path / "reference.fa")
    binding = Mock(spec=ReferenceBinding)
    binding.consumer_path = str(tmp_path / "bound-reference.fa")
    binding.close.side_effect = RuntimeError("cleanup exploded")

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=ValueError("probe exploded")),
        patch("vntyper.scripts.alignment_preflight.ReferenceBinding", return_value=binding),
        pytest.raises(ValueError, match="probe exploded"),
    ):
        resolve_reference(
            "/run/view.cram",
            (("cli", str(reference)),),
            "chr1:1-4",
            None,
            {},
            1,
            str(tmp_path),
            "sample",
            ("chr1",),
            "abc",
        )

    binding.close.assert_called_once_with()


def test_reference_cleanup_refusal_does_not_skip_independent_alignment_cleanup(tmp_path: Path) -> None:
    """One replaced reference view cannot strand the plan's independent alignment view."""
    output = tmp_path / "run"
    output.mkdir()
    alignment = tmp_path / "input.cram"
    alignment.write_bytes(b"alignment")
    alignment_view = output / "input.cram"
    alignment_binding = AlignmentBinding(str(alignment))
    alignment_binding.install_view(alignment_view)
    reference = tmp_path / "reference.fa"
    reference.write_bytes(b">chr1\nACGT\n")
    reference_binding = ReferenceBinding(str(reference), output, "input", 1)
    reference_view = Path(reference_binding.consumer_path)
    plan = AlignmentPlan(
        input_path=str(alignment),
        view_path=str(alignment_view),
        file_format="cram",
        index_path=str(output / "input.cram.crai"),
        reference_path=str(reference_view),
        reference_source="cli",
        uncovered_contigs=(),
        unmapped_scan="indexed",
        binding=alignment_binding,
        reference_binding=reference_binding,
    )
    reference_view.unlink()
    reference_view.write_bytes(b"replacement")

    with pytest.raises(RuntimeError, match="replaced CRAM reference view"):
        plan.close()

    assert reference_view.read_bytes() == b"replacement"
    assert reference_binding.is_open is True
    assert alignment_binding.is_open is False
    assert not os.path.lexists(alignment_view)
    reference_view.unlink()
    reference_binding.close()


def test_one_replaced_sidecar_does_not_skip_other_reference_cleanup(tmp_path: Path) -> None:
    """Reference cleanup attempts every independently owned view before reporting failure."""
    reference = _reference(tmp_path / "reference.fa")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    Path(f"{reference}.gzi").write_bytes(b"original GZI")
    binding = ReferenceBinding(str(reference), tmp_path, "sample", 1)
    reference_view = Path(binding.consumer_path)
    fai_view = Path(f"{binding.consumer_path}.fai")
    gzi_view = Path(f"{binding.consumer_path}.gzi")
    fai_view.unlink()
    fai_view.write_bytes(b"replacement FAI")

    with pytest.raises(RuntimeError, match="replaced CRAM reference view"):
        binding.close()

    assert fai_view.read_bytes() == b"replacement FAI"
    assert not os.path.lexists(reference_view)
    assert not os.path.lexists(gzi_view)
    assert binding.is_open is False
    fai_view.unlink()
    binding.close()


def test_generated_sidecar_keeps_its_name_instead_of_linking_to_its_own_descriptor(
    tmp_path: Path,
) -> None:
    """#238: a sidecar replaced by a link to its own FD only resolves where procfs jumps dentries."""
    reference = _reference(tmp_path / "reference.fa")
    output = tmp_path / "run"
    output.mkdir()
    binding = ReferenceBinding(str(reference), str(output), "sample", 1)
    try:
        generated = Path(f"{binding.consumer_path}.fai")
        generated.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

        binding.bind_generated_sidecars()

        assert not generated.is_symlink()
        assert generated.stat(follow_symlinks=False).st_nlink == 1
        assert generated.read_text(encoding="utf-8") == "chr1\t4\t6\t4\t5\n"
    finally:
        binding.close()
    assert not os.path.lexists(output / ".sample_reference_1")


def test_a_reference_view_replaced_while_being_proven_is_never_removed(tmp_path: Path) -> None:
    """The reachability fallback must not unlink an entry this binding did not install."""
    reference = _reference(tmp_path / "reference.fa")
    output = tmp_path / "run"
    output.mkdir()
    intruder = tmp_path / "intruder"
    intruder.write_bytes(b"not ours")

    def replace_then_report_unreachable(path: str | Path) -> tuple[None, str]:
        os.replace(intruder, path)
        return None, "Too many levels of symbolic links (errno 40)"

    with (
        patch(
            "vntyper.scripts.reference_binding.consumer_reachable_identity",
            side_effect=replace_then_report_unreachable,
        ),
        pytest.raises(RuntimeError, match="Refusing to remove replaced CRAM reference view"),
    ):
        ReferenceBinding(str(reference), str(output), "sample", 1)

    assert (output / ".sample_reference_1" / "reference.fa").read_bytes() == b"not ours"
