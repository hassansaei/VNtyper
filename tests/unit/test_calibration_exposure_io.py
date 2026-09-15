"""External exposure history cannot be reset by role/study/path relabelling."""

import os
from concurrent.futures import ThreadPoolExecutor
from importlib import import_module
from pathlib import Path
from unittest.mock import patch

import pytest

pytestmark = pytest.mark.unit


def _module():
    return import_module("vntyper.scripts.calibration_exposure_io")


def _identities(*, specimen="a", family="b", reads="c"):
    return [
        {"namespace": "family", "sha256": family * 64},
        {"namespace": "named-readset", "sha256": reads * 64},
        {"namespace": "specimen", "sha256": specimen * 64},
    ]


def _claim(path, ledger_id, **changes):
    values = {
        "expected_ledger_id": ledger_id,
        "target": "length",
        "role": "validation",
        "study_sha256": "d" * 64,
        "partition_sha256": "e" * 64,
        "evidence_sha256": "f" * 64,
        "identities": _identities(),
    }
    values.update(changes)
    return _module().record_exposure(path, **values)


def test_initialize_is_exclusive_private_and_claims_bind_immutable_history(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    assert len(ledger_id) == 64
    assert set(ledger_id) <= set("0123456789abcdef")
    assert path.stat().st_mode & 0o777 == 0o600
    receipt = _claim(path, ledger_id)
    assert receipt.sequence == 1
    assert receipt.target == "length" and receipt.role == "validation"
    assert receipt.exposure_ledger_id == ledger_id
    assert receipt.study_sha256 == "d" * 64
    document = _module().exposure_receipt_document(receipt)
    assert "identities" not in document
    assert document["membership_sha256"]
    before = path.read_bytes()
    with pytest.raises(ValueError):
        _module().initialize_exposure_ledger(path)
    assert path.read_bytes() == before
    with pytest.raises(ValueError, match="previously exposed"):
        _claim(path, ledger_id)
    assert path.read_bytes() == before


@pytest.mark.parametrize("prior", ["training", "policy-selection", "development-assessment", "validation"])
@pytest.mark.parametrize("next_role", ["validation", "locked-heldout"])
def test_previous_exposure_cannot_be_relabelled_for_confirmation(tmp_path: Path, prior, next_role):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    _claim(path, ledger_id, role=prior)
    with pytest.raises(ValueError, match="previously exposed"):
        _claim(path, ledger_id, role=next_role, study_sha256="0" * 64, partition_sha256="1" * 64)


@pytest.mark.parametrize("shared", ["specimen", "family", "reads"])
def test_any_stable_identity_collision_blocks_new_confirmation(tmp_path: Path, shared):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    _claim(path, ledger_id, role="development-assessment")
    values = {"specimen": "1", "family": "2", "reads": "3"}
    values[shared] = {"specimen": "a", "family": "b", "reads": "c"}[shared]
    with pytest.raises(ValueError, match="previously exposed"):
        _claim(path, ledger_id, identities=_identities(**values))


def test_history_is_target_specific_and_development_reuse_stays_recorded(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    first = _claim(path, ledger_id)
    other = _claim(path, ledger_id, target="callers")
    repeated = _claim(path, ledger_id, role="development-assessment")
    assert (first.sequence, other.sequence, repeated.sequence) == (1, 2, 3)
    assert len(path.read_bytes().splitlines()) == 4


def test_copying_history_or_changing_input_labels_does_not_reset_exposure(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    _claim(path, ledger_id)
    copied = tmp_path / "renamed.jsonl"
    copied.write_bytes(path.read_bytes())
    copied.chmod(0o600)
    with pytest.raises(ValueError, match="previously exposed"):
        _claim(copied, ledger_id)
    with pytest.raises(ValueError, match="identity"):
        _claim(path, "0" * 64, identities=_identities(specimen="1", family="2", reads="3"))


@pytest.mark.parametrize("corruption", ["torn", "modified", "extra-field", "empty"])
def test_corrupt_or_incomplete_history_is_never_reset(tmp_path: Path, corruption):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    _claim(path, ledger_id)
    original = path.read_bytes()
    if corruption == "torn":
        data = original + b'{"unfinished":'
    elif corruption == "modified":
        data = original.replace(b'"role":"validation"', b'"role":"training"')
    elif corruption == "extra-field":
        data = original.replace(b"{", b'{"unknown":true,', 1)
    else:
        data = b""
    path.write_bytes(data)
    with pytest.raises(ValueError):
        _claim(path, ledger_id, role="development-assessment")
    assert path.read_bytes() == data


def test_concurrent_confirmation_claims_serialize_before_outcome_access(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)

    def attempt(_):
        try:
            return _claim(path, ledger_id).sequence
        except ValueError as error:
            assert "previously exposed" in str(error)
            return None

    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(attempt, (1, 2)))
    assert sorted(result for result in results if result is not None) == [1]
    assert results.count(None) == 1


def test_partial_append_fails_closed_without_truncating_history(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    module = _module()
    ledger_id = module.initialize_exposure_ledger(path)
    real_write = os.write

    def torn_write(descriptor, data):
        real_write(descriptor, data[:20])
        raise OSError("synthetic interrupted append")

    with patch.object(module.os, "write", side_effect=torn_write), pytest.raises((ValueError, OSError)):
        _claim(path, ledger_id)
    damaged = path.read_bytes()
    with pytest.raises(ValueError):
        _claim(path, ledger_id)
    assert path.read_bytes() == damaged


def test_truncating_a_whole_record_cannot_restore_unseen_state(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    _claim(path, ledger_id)
    header_only = path.read_bytes().splitlines(keepends=True)[0]
    path.write_bytes(header_only)
    with pytest.raises(ValueError, match="head|incomplete|corrupt"):
        _claim(path, ledger_id)
    assert path.read_bytes() == header_only


def test_complete_append_without_a_committed_head_is_not_reused(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    module = _module()
    ledger_id = module.initialize_exposure_ledger(path)
    real_seek = os.lseek
    starts = 0

    def fail_head_rewind(descriptor, offset, whence):
        nonlocal starts
        if offset == 0 and whence == os.SEEK_SET:
            starts += 1
            if starts == 2:
                raise OSError("synthetic head update failure")
        return real_seek(descriptor, offset, whence)

    with patch.object(module.os, "lseek", side_effect=fail_head_rewind), pytest.raises(OSError):
        _claim(path, ledger_id)
    incomplete = path.read_bytes()
    with pytest.raises(ValueError):
        _claim(path, ledger_id, role="development-assessment")
    assert path.read_bytes() == incomplete


def test_ledger_must_be_external_and_nonsymlink(tmp_path: Path):
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / ".git").write_text("gitdir: example")
    with pytest.raises(ValueError, match="outside|external"):
        _module().initialize_exposure_ledger(repo / "ledger.jsonl")
    with pytest.raises(ValueError, match="absolute"):
        _module().initialize_exposure_ledger(Path("relative.jsonl"))
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    linked = tmp_path / "linked.jsonl"
    linked.symlink_to(path)
    with pytest.raises(ValueError):
        _claim(linked, ledger_id)
    with pytest.raises(ValueError, match="outside|external"):
        _claim(path, ledger_id, forbidden_roots=(tmp_path,))


@pytest.mark.parametrize(
    "changes",
    [
        {"target": "dominance"},
        {"role": "holdout"},
        {"identities": []},
        {"identities": [{"namespace": "specimen", "sha256": "a" * 64}]},
        {"identities": [{"namespace": "unknown", "sha256": "a" * 64}]},
        {"evidence_sha256": "bad"},
        {"study_sha256": True},
    ],
)
def test_invalid_claims_do_not_append_any_history(tmp_path: Path, changes):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    original = path.read_bytes()
    with pytest.raises(ValueError):
        _claim(path, ledger_id, **changes)
    assert path.read_bytes() == original


def test_nonprivate_ledger_refuses_access(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    path.chmod(0o644)
    before = path.read_bytes()
    with pytest.raises(ValueError, match="private regular"):
        _claim(path, ledger_id)
    assert path.read_bytes() == before


@pytest.mark.parametrize("bad", [1, [], ("directory",)])
def test_invalid_forbidden_roots_fail_before_creation(tmp_path: Path, bad):
    path = tmp_path / "ledger.jsonl"
    with pytest.raises(ValueError, match="Path tuple"):
        _module().initialize_exposure_ledger(path, forbidden_roots=bad)
    assert not path.exists()


def test_zero_write_and_fsync_failures_never_return_receipt(tmp_path: Path):
    path = tmp_path / "ledger.jsonl"
    module = _module()
    ledger_id = module.initialize_exposure_ledger(path)
    original = path.read_bytes()
    with patch.object(module.os, "write", return_value=0), pytest.raises(ValueError, match="did not complete"):
        _claim(path, ledger_id)
    assert path.read_bytes() == original
    with patch.object(module.os, "fsync", side_effect=OSError("synthetic disk failure")), pytest.raises(OSError):
        _claim(path, ledger_id)
    incomplete = path.read_bytes()
    assert incomplete != original
    with pytest.raises(ValueError, match="committed head"):
        _claim(path, ledger_id)
    assert path.read_bytes() == incomplete


@pytest.mark.parametrize("bad", [1, "-0000000000000000001", "1" * 21, "00009007199254740992", "٠" * 20])
def test_invalid_head_sequence_refuses_history(tmp_path: Path, bad):
    from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

    path = tmp_path / "ledger.jsonl"
    ledger_id = _module().initialize_exposure_ledger(path)
    header = load_strict_json_object(path.read_bytes())
    header["head_sequence"] = bad
    path.write_bytes(canonical_json_bytes(header))
    original = path.read_bytes()
    with pytest.raises(ValueError, match="head sequence"):
        _claim(path, ledger_id)
    assert path.read_bytes() == original
