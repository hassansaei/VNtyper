# Milestone 5 — References and Configuration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give VNtyper one source of truth for `(reference kind, assembly) → file path`, and move all reference data into versioned, checksummed, verified release bundles published from `berntpopp/vntyper-data`.

**Architecture:** `reference_registry.py` gains a resolver that maps an assembly *label* to a *physical reference identity* — `(coordinate_system, reference_source)` — and from there to a `reference_data` config key. All five readers and the one writer go through it. A new `reference_bundle.py` handles fetch/verify/stage/activate. `install-references` gains a `--from-source` path that is also what the bundle build workflow runs, so the derivation logic cannot rot.

**Tech Stack:** Python 3.10–3.13, pytest, ruff, mypy, Docker (two-stage base/app split), GitHub Actions reusable workflows, conda (`vntyper`, `envadvntr`, `shark_env`), BWA, samtools.

**Spec:** `docs/superpowers/specs/2026-08-11-milestone-5-references-and-configuration-design.md`

## Global Constraints

- **Base-image content hash.** `conda/**`, `docker/Dockerfile.base`, `docker/requirements-web.txt`, `vntyper/scripts/install_references.py`, `vntyper/scripts/install_references_config.json`, `vntyper/dependencies/advntr/**`, `reference/**`, `.dockerignore` are hashed at three call sites in two workflow files (`docker-base.yml:66`, `docker-build.yml:96`, `docker-build.yml:168`), mirrored in `Makefile` `BASE_INPUTS` and guarded by `tests/unit/test_workflow_consistency.py`. Editing any of them requires a base rebuild.
- **Maintainer branch only.** A fork's `GITHUB_TOKEN` cannot push to GHCR, so this work cannot go green from a fork (AGENTS.md trap 10).
- **Two PRs onto `main`, in order** (owner decision 2026-08-11; spec §4.6):
  1. `feat/milestone-5-reference-builder` — the builder only, all tracked data retained.
     Merge. → base rebuild #1.
  2. Publish `refs-v1` from `berntpopp/vntyper-data`, pinned at PR-1's **merge commit**.
  3. `fix/milestone-5-references-configuration` — the consumers. → base rebuild #2.

  Not stacked: both target `main`, so `ci-tests.yml` fires on each.
- **`install_references.py` and `reference_bundle.py` must import nothing from `vntyper`
  at module scope.** `Dockerfile.base` copies them into `/opt/ir` and runs them without
  installing the package, and anything they import must join the base-image content-hash
  set — where every future edit costs a 25–120 min rebuild. Inline the three-line
  `quote_path`, use `subprocess` directly rather than `utils.run_command`, and import the
  registry **function-locally** inside `canonical_reference_keys` (the `--config-path` path,
  which Docker never takes). A test asserts this.
- **Config is loaded at import time** into module globals in `kestrel_genotyping.py`, `advntr_genotyping.py`, `shark_filtering.py` (AGENTS.md trap 1). Resolvers must read the `config` argument, never the module global; tests patch the argument.
- **Membership, not truthiness**, at every config lookup. A key present with value `null` is an authoritative "disabled" and must not fall through. This matches `reference_resolution.py:50` and is pinned by `tests/unit/test_reference_resolution.py:68-84`.
- **`run_command` uses `shell=True`** with `executable="/bin/bash"`. Every interpolated path goes through `command_builders.quote_path`. Command *prefixes* from `config["tools"]` (e.g. `mamba run -n shark_env shark`) are never quoted.
- **The in-image reference layout must not change.** `/opt/vntyper/reference` with `WORKDIR /opt/vntyper`; `tests/docker/test_image_structure.py:244` pins that `config.json` reference paths stay relative.
- **Verify with `make check-all`.** Fast loop: `make format && make test-unit` from the repo root. `pytest` must come from the `vntyper` conda env — put `~/miniforge3/envs/vntyper/bin` first on `PATH` or `make check-all` reports success against the wrong interpreter.
- **Commit trailers** on every commit:
  ```
  Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR
  ```
  Commits that carry work derived from PR #164 additionally get `Co-Authored-By: ElenaPianfetti <elena.pianfetti@example.invalid>` — replace with the address from `git log` on that PR's head before committing.
- **The eight assembly labels** are exactly `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ncbi`, `hg38_ncbi`, `hg19_ensembl`, `hg38_ensembl` (`reference_registry.py:164-176`). They denote **six** physical files.

---

## Physical identity table

Every task below depends on this mapping. It is the whole design in one table.

| label | coord. system | source | physical id | `bwa_reference_*` keys, in order |
|---|---|---|---|---|
| `hg19` | GRCh37 | ucsc | `hg19` | `_hg19` |
| `hg38` | GRCh38 | ucsc | `hg38` | `_hg38` |
| `GRCh37` | GRCh37 | ncbi | `GRCh37` | `_GRCh37`, `_hg19` |
| `hg19_ncbi` | GRCh37 | ncbi | `GRCh37` | `_hg19_ncbi`, `_GRCh37`, `_hg19` |
| `GRCh38` | GRCh38 | ncbi | `GRCh38` | `_GRCh38`, `_hg38` |
| `hg38_ncbi` | GRCh38 | ncbi | `GRCh38` | `_hg38_ncbi`, `_GRCh38`, `_hg38` |
| `hg19_ensembl` | GRCh37 | ensembl | `hg19_ensembl` | `_hg19_ensembl`, `_hg19` |
| `hg38_ensembl` | GRCh38 | ensembl | `hg38_ensembl` | `_hg38_ensembl`, `_hg38` |

**Three tiers, not two.** The label tier is an existing tested capability:
`tests/unit/test_reference_resolution.py::test_label_specific_reference_keys_override_the_family_fallback_including_null`
pins that a replacement config may set `bwa_reference_hg19_ncbi` to specialise one label or
`cram_reference_hg19_ncbi: null` to disable it. A purely physical scheme routes `hg19_ncbi`
through `*_GRCh37` and silently ignores both — a regression disguised as a simplification.
The **writer** emits only the physical key, because that is the file it installed.
Duplicates are removed and order preserved, so a UCSC label yields one key.

---

## File Structure

**Created**
- `vntyper/scripts/reference_bundle.py` — fetch, verify, safe-extract, stage, activate. Pure of CLI concerns; `install_references.py` orchestrates it.
- `.github/workflows/build-reference-bundles.yml` — reusable (`workflow_call`) bundle builder.
- `tests/unit/test_reference_keys.py` — the resolver's table-driven tests.
- `tests/unit/test_reference_bundle.py` — fetch/verify/extract/rollback.
- `tests/unit/test_install_references_config_writing.py` — the writer (currently 0 % covered).
- `tests/unit/test_reference_key_agreement.py` — end-to-end writer↔reader agreement across 8 labels × 5 readers.

**Modified**
- `vntyper/scripts/reference_registry.py` — `REFERENCE_KINDS`, `physical_reference_id`, `ucsc_family`, `reference_keys`.
- `vntyper/scripts/reference_resolution.py` — `resolve_from_mapping`; `configured_reference_candidates` delegates.
- `vntyper/scripts/cli_handlers.py:261-275` — BWA selection.
- `vntyper/scripts/cli_logging_safety.py:43-51` — the fifth reader (safety boundary).
- `vntyper/scripts/pipeline.py:484-505` — adVNTR selection.
- `vntyper/modules/shark/shark_filtering.py` — layered SHARK selection.
- `vntyper/modules/shark/shark_config.json` — assembly-keyed.
- `vntyper/config.json` — six physical BWA keys, two SHARK keys.
- `vntyper/scripts/install_references.py`, `install_references_config.json` — bundle + `--from-source`.
- `docker/Dockerfile.base`, `docker/requirements-web.txt`, `docker/app/main.py`.
- `.github/workflows/docker-base.yml`, `docker-build.yml`, `scheduled-tests.yml`; `Makefile`; `.dockerignore`; `pyproject.toml`; `SPEC.md`; `AGENTS.md`.

**Deleted**
- `reference/*.fa`, `reference/*.fa.fai`, `reference/vntr_db_advntr.zip`, `reference/md5_checksums.txt`.
- `reference/generate_vntr_reference.py`, `reference/filter_config.json` — moved to `vntyper-data/seeds/`.

---

# PR 1 — the builder  (`feat/milestone-5-reference-builder`)
Tasks 3, 4, 5. **All tracked reference data stays in place**, and no consumer changes.
This PR must merge before anything can be built, because the bundle build runs its
`--from-source` path. → base rebuild #1.

Tasks 1 and 2 (the registry resolver) are *not* here: nothing in PR-1 needs them, and
keeping PR-1 to the builder is the whole point of splitting. They open PR-2.

---

### Task 3: Safe, transactional bundle installation  *(PR-1)*

**Files:**
- Create: `vntyper/scripts/reference_bundle.py`
- Test: `tests/unit/test_reference_bundle.py`

**Interfaces:**
- Consumes: nothing from earlier tasks — and **nothing from `vntyper` at all.** This module
  is copied bare into the Docker `refs` stage (spec §4.5b), so its only imports are
  `hashlib`, `inspect`, `logging`, `shutil`, `tarfile`, `tempfile`, `contextlib`, `pathlib`.
- Produces:
  - `verify_sha256(path: Path, expected: str) -> None` — `logger.error` then `ValueError`
  - `safe_extract(archive: Path, destination: Path) -> None` — `logger.error` then `ValueError`
  - `staged_install(target: Path, *, seed_from_existing: bool = True) -> ContextManager[Path]`

**No custom exception classes** — `AGENTS.md:140`: "Errors: no custom exception classes.
The convention is `logger.error(msg)` followed by a raise." Use `ValueError` for invalid
bundle contents and `RuntimeError` for download or installation failure.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_reference_bundle.py`:

```python
"""Installation must be all-or-nothing, and an archive must not write outside its root.

The installer this replaces skipped existing destinations without revalidating them
(`install_references.download_file`), unpacked archives with an unrestricted
`tar.extractall`, and overwrote config.json in place. A network blip therefore left a
half-populated reference tree that the next run treated as complete.
"""

import io
import tarfile

import pytest

from vntyper.scripts.reference_bundle import safe_extract, staged_install, verify_sha256

pytestmark = pytest.mark.unit


def _tar_with(tmp_path, members):
    """Build a tar whose members are (arcname, payload) pairs.

    Uses addfile rather than add: `TarFile.gettarinfo(arcname="/etc/passwd")` silently
    strips the leading slash, so an `add`-based fixture cannot produce the absolute-path
    member the traversal test is supposed to reject.
    """
    archive = tmp_path / "bundle.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for arcname, payload in members:
            info = tarfile.TarInfo(arcname)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return archive


class TestVerifySha256:
    def test_a_matching_digest_passes_silently(self, tmp_path):
        import hashlib

        target = tmp_path / "f"
        target.write_bytes(b"reference")
        verify_sha256(target, hashlib.sha256(b"reference").hexdigest())

    def test_a_mismatch_names_the_file_and_both_digests(self, tmp_path):
        target = tmp_path / "chr1.hg38.fa"
        target.write_bytes(b"tampered")
        with pytest.raises(ValueError, match="chr1.hg38.fa"):
            verify_sha256(target, "0" * 64)


class TestSafeExtract:
    def test_an_ordinary_archive_extracts(self, tmp_path):
        archive = _tar_with(tmp_path, [("alignment/chr1.hg38.fa", b">chr1\nACGT\n")])
        destination = tmp_path / "out"
        safe_extract(archive, destination)
        assert (destination / "alignment/chr1.hg38.fa").read_bytes() == b">chr1\nACGT\n"

    @pytest.mark.parametrize("arcname", ["../escaped.fa", "/etc/passwd", "a/../../escaped.fa"])
    def test_a_member_escaping_the_root_is_rejected_by_name(self, tmp_path, arcname):
        archive = _tar_with(tmp_path, [(arcname, b"x")])
        with pytest.raises(ValueError, match="escap|outside|absolute"):
            safe_extract(archive, tmp_path / "out")

    def test_a_hard_link_escaping_the_root_is_rejected(self, tmp_path):
        """Hard-link targets are archive-root-relative, unlike symlink targets."""
        archive = tmp_path / "hard.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            info = tarfile.TarInfo("nested/evil")
            info.type = tarfile.LNKTYPE
            info.linkname = "../../outside.fa"
            tar.addfile(info)
        with pytest.raises(ValueError):
            safe_extract(archive, tmp_path / "out")

    def test_a_symlink_pointing_outside_the_root_is_rejected(self, tmp_path):
        archive = tmp_path / "link.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            info = tarfile.TarInfo("evil")
            info.type = tarfile.SYMTYPE
            info.linkname = "/etc/passwd"
            tar.addfile(info)
        with pytest.raises(ValueError):
            safe_extract(archive, tmp_path / "out")


class TestStagedInstall:
    def test_an_existing_installation_is_carried_into_staging(self, tmp_path):
        """Installing hg38 after hg19 must not delete hg19, or the retained README."""
        target = tmp_path / "refs"
        (target / "alignment").mkdir(parents=True)
        (target / "alignment" / "chr1.hg19.fa").write_text("hg19")
        (target / "README.md").write_text("kept")
        with staged_install(target) as staging:
            (staging / "alignment" / "chr1.hg38.fa").write_text("hg38")
        assert (target / "alignment" / "chr1.hg19.fa").read_text() == "hg19"
        assert (target / "alignment" / "chr1.hg38.fa").read_text() == "hg38"
        assert (target / "README.md").read_text() == "kept"

    def test_a_successful_stage_is_activated_atomically(self, tmp_path):
        target = tmp_path / "refs"
        with staged_install(target) as staging:
            (staging / "chr1.fa").write_text(">chr1")
            assert not target.exists(), "target must not appear until activation"
        assert (target / "chr1.fa").read_text() == ">chr1"

    def test_a_failure_leaves_no_partial_tree_behind(self, tmp_path):
        target = tmp_path / "refs"
        with pytest.raises(RuntimeError):
            with staged_install(target) as staging:
                (staging / "half.fa").write_text("x")
                raise RuntimeError("network died mid-extract")
        assert not target.exists()

    def test_a_previous_tree_is_preserved_when_it_cannot_be_restored(self, tmp_path, monkeypatch):
        """Losing the only surviving copy is worse than leaving a stray directory."""
        import pathlib

        target = tmp_path / "refs"
        target.mkdir()
        (target / "chr1.fa").write_text("original")
        real_rename = pathlib.Path.rename

        def fail_restore(self, other):
            if pathlib.Path(other) == target and ".previous." in str(self):
                raise OSError("restore blocked")
            return real_rename(self, other)

        monkeypatch.setattr(pathlib.Path, "rename", fail_restore)
        with pytest.raises(RuntimeError):
            with staged_install(target) as staging:
                (staging / "chr1.fa").write_text("replacement")
                raise RuntimeError("boom")
        preserved = list(tmp_path.glob(".refs.previous.*"))
        assert preserved, "the previous tree must survive under a named path"
        assert (preserved[0] / "chr1.fa").read_text() == "original"

    def test_a_failure_leaves_a_previous_installation_untouched(self, tmp_path):
        target = tmp_path / "refs"
        target.mkdir()
        (target / "chr1.fa").write_text("original")
        with pytest.raises(RuntimeError):
            with staged_install(target) as staging:
                (staging / "chr1.fa").write_text("replacement")
                raise RuntimeError("boom")
        assert (target / "chr1.fa").read_text() == "original"
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
pytest tests/unit/test_reference_bundle.py -v
```

Expected: collection error, `ModuleNotFoundError: vntyper.scripts.reference_bundle`.

- [ ] **Step 3: Write the implementation**

Create `vntyper/scripts/reference_bundle.py`:

```python
"""Fetch, verify and install reference bundles without ever leaving a partial tree.

`install_references.py` orchestrates; this module owns the parts that must not be got
wrong: digest verification against a value committed in this repository, archive
extraction that cannot write outside its root, and an activation that is atomic from
the point of view of the next run.
"""

from __future__ import annotations

import hashlib
import inspect
import logging
import shutil
import tarfile
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger(__name__)

_CHUNK = 1024 * 1024

# `extractall(filter=...)` landed in 3.12 and was backported to 3.10.12 and 3.11.4.
# `requires-python` is ">=3.10", so it cannot be passed unconditionally.
_EXTRACT_KWARGS: dict[str, str] = (
    {"filter": "data"} if "filter" in inspect.signature(tarfile.TarFile.extractall).parameters else {}
)


def sha256_of(path: Path) -> str:
    """Return the hex SHA-256 of a file, read in chunks.

    Args:
        path: File to digest.

    Returns:
        str: Lowercase hex digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str) -> None:
    """Fail closed unless a file matches the digest committed in this repository.

    The release's own `SHA256SUMS` is co-hosted with the assets it describes, so it
    corroborates but cannot be the root of trust. The expected value comes from
    `install_references_config.json`, which is a base-image content-hash input - so
    changing any reference byte necessarily changes the base tag too.

    Args:
        path: File to check.
        expected: Lowercase hex SHA-256.

    Raises:
        ValueError: If the digests differ.
    """
    actual = sha256_of(path)
    if actual != expected:
        message = f"checksum mismatch for {path.name}: expected {expected}, got {actual}"
        logger.error(message)
        raise ValueError(message)
    logger.info(f"  ✓ verified {path.name}")


def _is_within(root: Path, candidate: Path) -> bool:
    """True when `candidate` resolves inside `root`."""
    try:
        candidate.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def safe_extract(archive: Path, destination: Path) -> None:
    """Extract a tar archive, rejecting any member that could write outside the root.

    The explicit member loop below is the guarantee, not `tarfile`'s `filter=`
    argument: `requires-python` is `>=3.10` and `filter=` only exists from 3.10.12,
    3.11.4 and 3.12, so passing it unconditionally would raise `TypeError` on an
    early 3.10. It is applied as defence in depth where available.

    Args:
        archive: `.tar.gz` to unpack.
        destination: Directory to unpack into; created if absent.

    Raises:
        ValueError: On an absolute path, a `..` component, a device or FIFO member,
            or a link resolving outside the destination. Per AGENTS.md the convention is
            `logger.error(message)` then `raise`, with no custom exception class.
    """
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            name = Path(member.name)
            if name.is_absolute():
                message = f"{archive.name}: absolute path in member '{member.name}'"
                logger.error(message)
                raise ValueError(message)
            if ".." in name.parts:
                message = f"{archive.name}: member '{member.name}' escapes the archive root"
                logger.error(message)
                raise ValueError(message)
            if not (member.isfile() or member.isdir() or member.issym() or member.islnk()):
                message = f"{archive.name}: member '{member.name}' is not a regular file or link"
                logger.error(message)
                raise ValueError(message)
            # Symlink targets are relative to the member's own directory; hard-link
            # targets are relative to the archive root. Resolving both the same way
            # misjudges one of them, and on a Python without `filter=` this loop is the
            # only protection there is.
            if member.issym():
                linked = destination / name.parent / member.linkname
            elif member.islnk():
                linked = destination / member.linkname
            else:
                linked = None
            if linked is not None and not _is_within(destination, linked):
                message = f"{archive.name}: link '{member.name}' escapes the archive root"
                logger.error(message)
                raise ValueError(message)
        tar.extractall(path=destination, **_EXTRACT_KWARGS)


@contextmanager
def staged_install(target: Path, *, seed_from_existing: bool = True) -> Iterator[Path]:
    """Build a reference tree beside its destination and activate it only on success.

    The staging directory is a sibling of `target`, so activation is a rename on the
    same filesystem. On any exception the staging directory is removed and any
    pre-existing installation is restored, or - if it cannot be restored - preserved
    under a named path and reported. The previous tree is deleted only after activation
    has succeeded, so no failure path can leave the caller with nothing.

    `seed_from_existing` copies the current tree into the staging directory first, so an
    install of one assembly does not erase another. Without it,
    `install-references --references hg38` after `--references hg19` would silently
    delete hg19, and installing into the tracked `reference/` directory would delete the
    retained `README.md` and `pseudonymize*` files.

    Args:
        target: Final directory.
        seed_from_existing: Copy the existing tree into staging before yielding, so the
            install merges rather than replaces. False only for a deliberate clean install.

    Yields:
        Path: The staging directory to populate.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}.staging.", dir=target.parent))
    previous: Path | None = None
    if seed_from_existing and target.exists():
        shutil.copytree(target, staging, dirs_exist_ok=True, symlinks=True)
    try:
        yield staging
        if target.exists():
            previous = Path(tempfile.mkdtemp(prefix=f".{target.name}.previous.", dir=target.parent))
            previous.rmdir()
            target.rename(previous)
        staging.rename(target)
        staging = None  # type: ignore[assignment]
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        if previous is not None:
            # Restore only into a vacant slot, and NEVER delete `previous` on this path:
            # if something else recreated `target`, or the rename back fails, the previous
            # installation is the only surviving copy. Leaving it on disk under a named
            # path is recoverable; deleting it is not.
            if not target.exists():
                try:
                    previous.rename(target)
                    previous = None
                except OSError as restore_error:
                    logger.error(f"could not restore {target} from {previous}: {restore_error}")
            if previous is not None:
                logger.error(f"previous reference tree preserved at {previous}; restore it by hand")
        raise
    else:
        # Activation is confirmed; only now is the old tree disposable.
        if previous is not None:
            shutil.rmtree(previous, ignore_errors=True)
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
pytest tests/unit/test_reference_bundle.py -v
```

Expected: PASS. If `test_a_matching_digest_passes_silently` looks odd, simplify it to
compute the digest inline — the unused literal is dead weight, delete it.

- [ ] **Step 5: Lint and type-check**

```bash
ruff check vntyper/scripts/reference_bundle.py && mypy vntyper/scripts/reference_bundle.py
```

Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add vntyper/scripts/reference_bundle.py tests/unit/test_reference_bundle.py
git commit -m "feat(references): transactional bundle install with constrained extraction

Downloads and extraction now happen in a staging directory beside the target and
are activated by rename, so an interruption cannot leave a half-populated tree
that the next run treats as complete. Archive members with absolute paths, ..
components, device nodes or escaping links are rejected by name.

Refs #217

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 4: Derivations and `--from-source`  *(PR-1)*

**Files:**
- Modify: `vntyper/scripts/install_references.py`
- Modify: `vntyper/scripts/install_references_config.json`
- Modify: `vntyper/scripts/cli_parser.py:279-316` (add `--from-source`)
- Test: `tests/unit/test_install_references_derivations.py`

**Interfaces:**
- Consumes: `reference_bundle.verify_sha256` (Task 3).
- Produces:
  - `derive_region_fasta(source_fasta: Path, region: str, destination: Path, samtools: str) -> Path`
  - `run_derivations(install_config: dict, output_dir: Path, config: dict) -> None`
  - a `--from-source` flag on `install-references`

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_install_references_derivations.py`:

```python
"""The MUC1 region FASTAs are derived, not shipped - and the derivation is verified.

`samtools faidx chr1.hg19.fa chr1:155158000-155163000` reproduces the FASTA this
repository used to track, byte for byte (md5 c9737129069d4855b433b178ebb21e1c). That
is why the file could be removed from git: the bundle build asserts the derivation
against a digest committed in the release spec, so a silent change in a UCSC
chromosome file turns the build red instead of publishing different sequence under an
unchanged name.
"""

import subprocess

import pytest

from vntyper.scripts import install_references

pytestmark = pytest.mark.unit


class TestDeriveRegionFasta:
    def test_the_command_names_the_region_and_redirects_to_the_destination(self, tmp_path, monkeypatch):
        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess, "run",
            lambda cmd, **kw: calls.append(cmd) or subprocess.CompletedProcess(cmd, 0, "", ""),
        )

        destination = tmp_path / "muc1_region_hg38.fa"
        install_references.derive_region_fasta(
            source_fasta=tmp_path / "chr1.hg38.fa",
            region="chr1:155184000-155194000",
            destination=destination,
            samtools="samtools",
        )

        assert "faidx" in calls[0]
        assert "chr1:155184000-155194000" in calls[0]
        assert str(destination) in calls[0]

    def test_paths_with_spaces_are_quoted_into_single_operands(self, tmp_path, monkeypatch):
        """run_command runs this under bash, so an unquoted path with a space splits."""
        import shlex

        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess, "run",
            lambda cmd, **kw: calls.append(cmd) or subprocess.CompletedProcess(cmd, 0, "", ""),
        )

        source = tmp_path / "my refs" / "chr1.hg38.fa"
        install_references.derive_region_fasta(source, "chr1:1-2", tmp_path / "o.fa", "samtools")

        assert str(source) in shlex.split(calls[0])

    def test_a_failed_derivation_raises_rather_than_producing_an_empty_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            install_references.subprocess, "run",
            lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "no such region"),
        )
        with pytest.raises(RuntimeError, match="faidx|deriv"):
            install_references.derive_region_fasta(tmp_path / "c.fa", "chr1:1-2", tmp_path / "o.fa", "samtools")
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
pytest tests/unit/test_install_references_derivations.py -v
```

Expected: FAIL, `AttributeError: module ... has no attribute 'derive_region_fasta'`.

- [ ] **Step 3: Add the derivation config**

In `vntyper/scripts/install_references_config.json`, add a top-level `derivations` block
after `own_repository_references`:

```json
  "derivations": [
    {
      "kind": "shark", "assembly": "hg19",
      "output": "muc1_region_hg19.fa",
      "from": "hg19", "region": "chr1:155158000-155163000",
      "tool": "samtools faidx",
      "expected_sha256": "<computed in Step 3 and committed here>"
    },
    {
      "kind": "shark", "assembly": "hg38",
      "output": "muc1_region_hg38.fa",
      "from": "hg38", "region": "chr1:155184000-155194000",
      "tool": "samtools faidx",
      "expected_sha256": "<computed in Step 3 and committed here>"
    },
    {
      "kind": "literal", "config_key": "muc1_reference_vntr",
      "output": "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
      "from_seeds": ["MUC1_motifs_Rev_com.fa", "filter_config.json"],
      "tool": "generate_vntr_reference",
      "expected_sha256": "<computed in Step 3 and committed here>"
    }
  ],
  "common_references": [
    {"config_key": "muc1_motifs_rev_com",        "installed_path": "MUC1_motifs_Rev_com.fa"},
    {"config_key": "code_adVNTR_RUs",            "installed_path": "code-adVNTR_RUs.fa"},
    {"config_key": "advntr_reference_vntr_hg19", "installed_path": "vntr_db_advntr/hg19_muc1.db"},
    {"config_key": "advntr_reference_vntr_hg38", "installed_path": "vntr_db_advntr/hg38_muc1.db"}
  ]
```

**This is the one schema the whole milestone depends on. Fix it here and do not vary it.**
Every genome entry is named by its physical id and carries `kind`, `installed_path`,
`asset`/`asset_sha256` (bundle path) and `source_url`/`source_sha256` (`--from-source`
path). Every derivation carries either `kind: "shark"` plus `assembly` — so its config key
comes from the registry — or `kind: "literal"` plus an explicit `config_key`. Every common
reference carries `config_key` plus `installed_path`. Task 10's `canonical_reference_keys`
reads exactly these fields and nothing else, and a schema-validation test rejects an entry
missing any of them.

Compute the digests once with the local references and paste them in:

```bash
export PATH=~/miniforge3/envs/vntyper/bin:$PATH
samtools faidx reference/alignment/chr1.hg19.fa chr1:155158000-155163000 | sha256sum
samtools faidx reference/alignment/chr1.hg38.fa chr1:155184000-155194000 | sha256sum
```

- [ ] **Step 4: Write the implementation**

In `vntyper/scripts/install_references.py`, add the import and the two functions after
`process_own_repository_references`:

**Import only `reference_bundle`.** `command_builders` (11 commits/yr, and it drags in
`samtools_command_fragments`) and `utils` (8 commits/yr, imports `command_builders`) would
both join the base-image content-hash set — see the Global Constraints. Inline the helper
and shell out directly instead.

```python
import shlex
import subprocess

from vntyper.scripts.reference_bundle import verify_sha256


def _quote(value: object) -> str:
    """Quote one operand for a shell command.

    A three-line copy of `command_builders.quote_path`, deliberately duplicated: importing
    that module would put it and its own dependencies into the base-image content hash,
    where every future edit costs a 25-120 minute rebuild. `tests/unit/test_shell_quoting.py`
    pins that both behave identically.
    """
    return shlex.quote(str(value))


def derive_region_fasta(source_fasta: Path, region: str, destination: Path, samtools: str) -> Path:
    """Cut a MUC1 region out of a chromosome FASTA.

    This is how both SHARK references are produced. The hg19 one reproduces the FASTA
    this repository used to track byte-for-byte, which is what makes deriving rather
    than shipping them safe.

    Args:
        source_fasta: Indexed chromosome FASTA to cut from.
        region: `chr1:start-end`, in the source's own contig naming.
        destination: File to write.
        samtools: samtools command prefix from the install config.

    Returns:
        Path: `destination`.

    Raises:
        RuntimeError: If samtools fails.
    """
    command = f"{samtools} faidx {_quote(source_fasta)} {_quote(region)} > {_quote(destination)}"
    completed = subprocess.run(command, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if completed.returncode != 0:
        message = f"faidx derivation failed for {destination.name}: {completed.stderr.strip()}"
        logger.error(message)
        raise RuntimeError(message)
    return destination


def run_derivations(install_config: dict[str, Any], output_dir: Path, references: dict[str, Path]) -> None:
    """Run every configured derivation and verify it against its committed digest.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree being populated.
        references: Physical id to installed chromosome FASTA, for `from` lookups.

    Raises:
        RuntimeError: If a source reference is missing.
        ValueError: If a derived file does not match `expected_sha256`.
    """
    samtools = install_config.get("samtools_path", "samtools")
    for spec in install_config.get("derivations", []):
        source = references.get(spec["from"])
        if source is None or not source.exists():
            raise RuntimeError(f"cannot derive {spec['output']}: {spec['from']} is not installed")
        destination = output_dir / spec["output"]
        derive_region_fasta(source, spec["region"], destination, samtools)
        verify_sha256(destination, spec["expected_sha256"])
        logger.info(f"Derived {spec['output']} from {source.name} at {spec['region']}")
```

Add `"samtools_path": "samtools"` next to `"bwa_path": "bwa"` in the install config.

- [ ] **Step 4b: Implement `install_from_source` and route `main` to it**

Threading a flag through signatures is not enough — `main` (`:726-773`) must branch. Add:

```python
def install_from_source(
    install_config: dict[str, Any],
    output_dir: Path,
    references: list[str],
    aligners: dict[str, dict[str, Any]],
    index_threads: int,
    release_spec: dict[str, Any] | None = None,
) -> dict[str, Path]:
    """Build every requested reference from its upstream source.

    This is what the bundle build workflow runs, so it is the definition of what a
    bundle contains - there is no second derivation implementation to drift from it.

    When `release_spec` is given (release builds), every URL and digest comes from the
    spec and the downloaded bytes are verified **before** decompression or indexing.
    Without it, the URLs and `source_sha256` values in install_references_config.json
    are used. The spec is what pins Ensembl to an explicit release rather than the
    mutable `current` path the shipped config still names.

    Returns:
        dict[str, Path]: physical id -> installed chromosome FASTA, for `run_derivations`.
    """
```

It must: download each selected physical reference, `verify_sha256` it, decompress,
index with the enabled aligners, install the common seeds, then call `run_derivations`
with the returned mapping. `main` selects it when `from_source` is true and
`install_from_bundle` (Task 11) otherwise.

Add `--release-spec PATH` to the `install-references` parser, documented as "release
builds only: take every source URL and digest from this file". Test that supplying it
overrides a config URL and that a digest mismatch fails before any decompression.

- [ ] **Step 5: Add the `--from-source` flag**

In `vntyper/scripts/cli_parser.py`, after the `--references` argument of `parser_install`
(`:308-316`):

```python
    parser_install.add_argument(
        "--from-source",
        action="store_true",
        help=(
            "Build references from their upstream sources instead of downloading the "
            "published bundle. Slower (downloads and BWA-indexes six chromosome FASTAs) "
            "but needs no access to the reference release. This is the path the bundle "
            "build workflow itself runs."
        ),
    )
```

Thread it through `cli_handlers.handle_install_references` (`:139-160`) and
`install_references.main`'s signature as `from_source: bool = False`.

- [ ] **Step 6: Run the tests**

```bash
pytest tests/unit/test_install_references_derivations.py tests/unit/test_install_references.py \
       tests/unit/test_cli_parser.py -v
```

Expected: all PASS.

- [ ] **Step 7: Verify the derivation end-to-end against the real references**

```bash
export PATH=~/miniforge3/envs/vntyper/bin:$PATH
vntyper install-references -d /tmp/m5_src --from-source --references hg19 hg38 --threads 8
md5sum /tmp/m5_src/muc1_region_hg19.fa reference/muc1_region_hg19.fa
```

Expected: identical md5 (`c9737129069d4855b433b178ebb21e1c`), and
`/tmp/m5_src/muc1_region_hg38.fa` exists at 10,228 bytes.

- [ ] **Step 8: Commit**

```bash
git add vntyper/scripts/install_references.py vntyper/scripts/install_references_config.json \
        vntyper/scripts/cli_parser.py vntyper/scripts/cli_handlers.py \
        tests/unit/test_install_references_derivations.py
git commit -m "feat(install-references): derive the MUC1 region FASTAs, add --from-source

Both SHARK references are cuts of a chromosome FASTA, verified against a digest
committed in the install config. The hg19 one reproduces the tracked file byte for
byte, which is what lets the tracked copies be removed later. --from-source is the
path the bundle build workflow runs, so the derivation logic stays exercised.

Refs #152 #217

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 5: The reusable bundle-build workflow  *(PR-1 — last task before merge)*

**Files:**
- Create: `.github/workflows/build-reference-bundles.yml`
- Test: manual dispatch (Task 6)

**Interfaces:**
- Consumes: `--from-source` from Task 4.
- Produces: a draft GitHub release in the **calling** repository carrying per-assembly
  `.tar.gz` assets, `SHA256SUMS`, `release-manifest.json`, `verification-report.json`,
  `BUILD_INFO.json`.

- [ ] **Step 1: Write the workflow**

Create `.github/workflows/build-reference-bundles.yml`:

```yaml
# Builds a versioned reference bundle and publishes it as a DRAFT release in the
# calling repository (berntpopp/vntyper-data).
#
# It lives here, not in the data repo, for the reason phentrieve does the same: the
# builder is software and belongs with the software. The data repo pins this file at a
# full commit SHA so a bundle is always reproducible from a known builder.
#
# Two checkouts are required. workflow_call does NOT check out the caller, and the
# seeds and the release spec live there while the builder lives here.
name: Build reference bundles

on:
  workflow_call:
    inputs:
      release_tag:
        description: Release tag declared by a committed JSON spec in the caller.
        required: true
        type: string
      release_spec_path:
        description: Path to that spec, relative to the caller's repository root.
        required: true
        type: string
      source_commit:
        description: The VNtyper commit this builder was pinned at.
        required: true
        type: string

permissions:
  contents: write

jobs:
  build:
    runs-on: ubuntu-24.04
    timeout-minutes: 240
    steps:
      - name: Check out the data repository (seeds and release spec)
        uses: actions/checkout@v5
        with:
          path: data

      - name: Check out VNtyper at the pinned builder commit
        uses: actions/checkout@v5
        with:
          repository: hassansaei/VNtyper
          ref: ${{ inputs.source_commit }}
          path: vntyper

      - name: Record both resolved SHAs
        id: shas
        run: |
          set -euo pipefail
          echo "data=$(git -C data rev-parse HEAD)"     >> "$GITHUB_OUTPUT"
          echo "builder=$(git -C vntyper rev-parse HEAD)" >> "$GITHUB_OUTPUT"

      - uses: conda-incubator/setup-miniconda@v3
        with:
          miniforge-version: latest
          activate-environment: vntyper
          environment-file: vntyper/conda/environment_vntyper.yml

      - name: Install VNtyper
        shell: bash -el {0}
        run: pip install --no-deps -e vntyper

      - name: Stage the seeds into the reference tree
        shell: bash -el {0}
        run: |
          set -euo pipefail
          mkdir -p refs
          cp data/seeds/MUC1_motifs_Rev_com.fa data/seeds/code-adVNTR_RUs.fa refs/
          cp data/seeds/vntr_db_advntr.zip refs/
          python vntyper/vntyper/scripts/verify_seeds.py \
            --spec "data/${{ inputs.release_spec_path }}" --seeds data/seeds

      - name: Build every reference from source
        shell: bash -el {0}
        run: |
          set -euo pipefail
          vntyper install-references \
            --output-dir refs --from-source \
            --release-spec "data/${{ inputs.release_spec_path }}" \
            --references hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl \
            --threads 4

      - name: Assemble assets, checksums and manifests
        shell: bash -el {0}
        run: |
          set -euo pipefail
          python vntyper/scripts/bundle_release.py \
            --refs refs --spec "data/${{ inputs.release_spec_path }}" \
            --tag "${{ inputs.release_tag }}" --out dist \
            --data-sha "${{ steps.shas.outputs.data }}" \
            --builder-sha "${{ steps.shas.outputs.builder }}"
          ( cd dist && sha256sum *.tar.gz > SHA256SUMS )

      - name: Publish a DRAFT release
        env:
          GH_TOKEN: ${{ github.token }}
        run: |
          set -euo pipefail
          gh release create "${{ inputs.release_tag }}" dist/* \
            --repo "${{ github.repository }}" \
            --draft --title "VNtyper references ${{ inputs.release_tag }}" \
            --notes "Builder: hassansaei/VNtyper@${{ steps.shas.outputs.builder }}
          Spec: ${{ inputs.release_spec_path }} @ ${{ steps.shas.outputs.data }}

          Verify before publishing: check verification-report.json, then
          \`sha256sum --check SHA256SUMS\`."
```

- [ ] **Step 2: Write the two helper scripts the workflow calls**

Create `vntyper/scripts/verify_seeds.py` — reads the spec's `seeds` block and calls
`reference_bundle.verify_sha256` on each; exits non-zero on the first mismatch.

Create `scripts/bundle_release.py` — walks `refs/`, groups files into the seven assets
named in the spec (§4.1 of the design), writes `release-manifest.json` (path, size,
sha256, source URL, derivation command per entry), `BUILD_INFO.json`
(`bwa --version`, `samtools --version`, runner, both SHAs) and
`verification-report.json` (every digest checked and its verdict), then tars each asset.

Both are ordinary Python and get unit tests in `tests/unit/test_bundle_release.py`
covering: asset grouping, manifest completeness (every file in `refs/` appears exactly
once), and a spec entry with no matching file failing loudly.

- [ ] **Step 3: Validate the workflow syntax**

```bash
python -c "import yaml,sys; yaml.safe_load(open('.github/workflows/build-reference-bundles.yml'))"
pytest tests/unit/test_workflow_consistency.py -v
make ci-local
```

`yaml.safe_load` proves only that it parses — it validates no Actions semantics. AGENTS.md
requires `make ci-local` after any workflow edit; `make check-all` is
`format-check lint type-check-all test-unit` (`Makefile:409`) and runs neither
`lint-actions` nor anything Docker.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/build-reference-bundles.yml vntyper/scripts/verify_seeds.py \
        scripts/bundle_release.py tests/unit/test_bundle_release.py
git commit -m "ci: reusable workflow that builds and drafts a reference bundle release

Lives here rather than in the data repo because the builder is software; the data
repo pins it at a full SHA. Checks out both repositories explicitly, since
workflow_call does not check out the caller and the seeds live there.

Refs #217

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

# INTERLUDE — publish `refs-v1`  (Task 6)

PR-1 is merged. Create the data repo, pin it at PR-1's **merge commit**, build, verify,
publish.

---

### Task 6: Create `berntpopp/vntyper-data` and publish `refs-v1`  *(between the PRs)*

**Files:** (in the new repository, not this one)
- Create: `LICENSE`, `ATTRIBUTION.md`, `README.md`, `.gitignore`, `seeds/`, `releases/refs-v1.json`, `.github/workflows/release-data.yml`

**Interfaces:**
- Consumes: `build-reference-bundles.yml` from Task 5, at this branch's head SHA.
- Produces: the published `refs-v1` release and the per-asset SHA-256 digests that Task 8 commits.

- [ ] **Step 1: Merge PR-1, then capture the merge commit**

```bash
gh pr merge --repo hassansaei/VNtyper <PR-1 number> --merge
git fetch origin main
BUILDER_SHA=$(git rev-parse origin/main)
echo "$BUILDER_SHA"
```

This merged SHA is what the data repo pins. It must be a **merge commit on `main`**, not a
branch head: the release spec records it as the provenance of every published byte, and it
has to stay reachable and reviewed forever.

Confirm PR-1 is actually in it before proceeding:

```bash
git show "$BUILDER_SHA:vntyper/scripts/reference_bundle.py" >/dev/null && echo "builder present"
git show "$BUILDER_SHA:.github/workflows/build-reference-bundles.yml" >/dev/null && echo "workflow present"
```

**Published releases are immutable — there is no edit path.** If the builder ever changes
after publication, the answer is a new tag, never a re-cut `refs-v1`: commit
`releases/refs-v2.json` pinned at the new merge SHA, add `refs-v2` to `release-data.yml`'s
`options`, build and verify it, update `release_tag` and every `asset_sha256` in
`install_references_config.json`, and only then proceed. Task 13's Step 1 enforces this by
refusing to delete anything the published bundle does not already carry.

- [x] **Step 2: Create the repository and copy the seeds — ALREADY DONE (2026-08-11)**

`berntpopp/vntyper-data` exists, is public, MIT-licensed, and already carries `LICENSE`,
`ATTRIBUTION.md`, `README.md`, `.gitignore` and all seven seeds. The seeds were taken from
`origin/main` and verified byte-identical; PR-1 does not touch `reference/`, so they are
also the bytes at `BUILDER_SHA`. **Re-verify before writing the spec** rather than assuming:

```bash
git clone https://github.com/berntpopp/vntyper-data /tmp/vntyper-data
cd /tmp/vntyper-data
for f in seeds/*; do
  n=$(basename "$f")
  a=$(git -C ~/development/VNtyper show "$BUILDER_SHA:reference/$n" | sha256sum | cut -d' ' -f1)
  b=$(sha256sum "$f" | cut -d' ' -f1)
  [ "$a" = "$b" ] && echo "OK   $n" || echo "DIFF $n - re-materialise from BUILDER_SHA"
done
```

The recorded digests, for `releases/refs-v1.json`'s `seeds` block:

```
7e6589f2388f3a08da6fb3bffa1fe22f10a5e03ec618b883fa6f07bcf1cb3e47  MUC1_motifs_Rev_com.fa
b570810bd7d43d9cb8f7e64aa080d71a13722071c9a2b4b453d9d9a4465e5603  MUC1_motifs_Rev_com.fa.fai
c21d631cf894e388c8cb76d7bbd2a51ebf4a27cd1a9158f50971cd831d8aa26e  code-adVNTR_RUs.fa
13efec9be6844aa9e3dc4547b529909e588f5895f015e82ee242cd9e6d8503b5  code-adVNTR_RUs.fa.fai
90a619f6aa2ee7d038b6d8703a5736d92fd483e8b4bfad4a5ad07480bf8f7ff1  vntr_db_advntr.zip
d2190ed78695efe9b1b8105c97479391b81129cf641410dfb88feb1c1ffea085  filter_config.json
6dcf12357cb648cbf66def71c9abde1e76b5bd6de62aa27a3757055052645c33  generate_vntr_reference.py
```
```

Those digests go into the spec's `seeds` block in the next step.

- [ ] **Step 3: Write the release spec**

`releases/refs-v1.json`, following `berntpopp/phentrieve-data/releases/*.json`. Pin
Ensembl to an **explicit release**, never `current` — `install_references_config.json:72,78`
currently track `current`, which is why no existing base image records which Ensembl
release it contains. Fill `sha256` for each source by downloading it once and hashing it.

```json
{
  "release_tag": "refs-v1",
  "source_commit": "<BUILDER_SHA from step 1>",
  "bwa_version": "<bwa 2>&1 | head -3>",
  "samtools_version": "<samtools --version | head -1>",
  "sources": {
    "hg19":         {"url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz", "sha256": "..."},
    "hg38":         {"url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz", "sha256": "..."},
    "GRCh37":       {"url": "https://ftp.ncbi.nlm.nih.gov/.../GCF_000001405.25_GRCh37.p13/.../chr1.fna.gz", "sha256": "..."},
    "GRCh38":       {"url": "https://ftp.ncbi.nlm.nih.gov/.../GCF_000001405.40_GRCh38.p14/.../chr1.fna.gz", "sha256": "..."},
    "hg19_ensembl": {"url": "https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.1.fa.gz", "sha256": "..."},
    "hg38_ensembl": {"url": "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz", "sha256": "..."}
  },
  "seeds": {
    "MUC1_motifs_Rev_com.fa": {"sha256": "..."},
    "code-adVNTR_RUs.fa":     {"sha256": "..."},
    "vntr_db_advntr.zip":     {"sha256": "..."},
    "filter_config.json":     {"sha256": "..."}
  },
  "derivations": [
    {"output": "muc1_region_hg19.fa", "command": "samtools faidx chr1.hg19.fa chr1:155158000-155163000", "expected_sha256": "..."},
    {"output": "muc1_region_hg38.fa", "command": "samtools faidx chr1.hg38.fa chr1:155184000-155194000", "expected_sha256": "..."},
    {"output": "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa", "command": "python generate_vntr_reference.py", "expected_sha256": "..."}
  ]
}
```

- [ ] **Step 4: Write the dispatch workflow** (LICENSE, ATTRIBUTION.md, README.md and
      .gitignore are already written and pushed — only `release-data.yml` remains, because
      it is the one file that must pin `BUILDER_SHA`)

`.github/workflows/release-data.yml`:

```yaml
name: Build VNtyper reference release

on:
  workflow_dispatch:
    inputs:
      release_tag:
        description: Release tag declared by a committed JSON spec.
        required: true
        type: choice
        options:
          - refs-v1

permissions:
  contents: write

jobs:
  build:
    uses: hassansaei/VNtyper/.github/workflows/build-reference-bundles.yml@<BUILDER_SHA>
    with:
      release_tag: ${{ inputs.release_tag }}
      release_spec_path: releases/${{ inputs.release_tag }}.json
      source_commit: <BUILDER_SHA>
```

Then update `README.md`'s **Status** section — it currently says `refs-v1` is not yet
published and is blocked on the builder landing in VNtyper. That stops being true here.

- [ ] **Step 5: Dispatch, verify, publish**

```bash
git add -A && git commit -m "feat: refs-v1 release specification" && git push
gh workflow run release-data.yml -f release_tag=refs-v1 --repo berntpopp/vntyper-data
gh run watch --repo berntpopp/vntyper-data
```

When the draft appears:

```bash
gh release download refs-v1 --repo berntpopp/vntyper-data --dir /tmp/refs-v1
cd /tmp/refs-v1 && sha256sum --check SHA256SUMS
python -m json.tool verification-report.json | head -50
```

Expected: every checksum OK; the report shows all three derivations matching their
`expected_sha256`. **Then** publish:

```bash
gh release edit refs-v1 --repo berntpopp/vntyper-data --draft=false
gh release view refs-v1 --repo berntpopp/vntyper-data \
  --json assets --jq '.assets[] | "\(.name)"'
```

- [ ] **Step 6: Record the per-asset digests for Task 8**

```bash
cd /tmp/refs-v1 && sha256sum *.tar.gz | tee /tmp/refs-v1-asset-digests.txt
```

---

---

# PR 2 — the consumers  (`fix/milestone-5-references-configuration`)

Tasks 1, 2, 7-17. The bundle is published; tracked reference data is deleted in Task 13.
→ base rebuild #2.

---

### Task 1: The reference key resolver  *(PR-2)*

**Files:**
- Modify: `vntyper/scripts/reference_registry.py` (append after `resolve_chromosome_name`, before the validation section at `:496`)
- Test: `tests/unit/test_reference_keys.py`

**Interfaces:**
- Consumes: `get_coordinate_system`, `get_reference_source`, `normalize_assembly_name`, `ASSEMBLY_METADATA` — all already in this module.
- Produces:
  - `REFERENCE_KINDS: dict[str, dict[str, str]]`
  - `physical_reference_id(assembly: str) -> str`
  - `ucsc_family(assembly: str) -> str`
  - `reference_keys(kind: str, assembly: str) -> tuple[str, ...]`

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_reference_keys.py`:

```python
"""The single mapping from an assembly label to the config keys that name its files.

Eight labels denote six physical files: `GRCh37`/`hg19_ncbi` are one NCBI file and
`GRCh38`/`hg38_ncbi` are another. Keying on the label cannot work, because
`install_references_config.json` only ever produces a `GRCh37` and a `GRCh38` entry -
a label-keyed writer would never emit `bwa_reference_hg38_ncbi`, and that run would
silently fall back to the UCSC FASTA.
"""

import pytest

from vntyper.scripts.reference_registry import (
    list_assemblies,
    physical_reference_id,
    reference_keys,
    ucsc_family,
)

pytestmark = pytest.mark.unit

PHYSICAL_ID = {
    "hg19": "hg19",
    "hg38": "hg38",
    "GRCh37": "GRCh37",
    "hg19_ncbi": "GRCh37",
    "GRCh38": "GRCh38",
    "hg38_ncbi": "GRCh38",
    "hg19_ensembl": "hg19_ensembl",
    "hg38_ensembl": "hg38_ensembl",
}

UCSC_FAMILY = {
    "hg19": "hg19", "GRCh37": "hg19", "hg19_ncbi": "hg19", "hg19_ensembl": "hg19",
    "hg38": "hg38", "GRCh38": "hg38", "hg38_ncbi": "hg38", "hg38_ensembl": "hg38",
}


class TestPhysicalIdentity:
    @pytest.mark.parametrize("label,expected", sorted(PHYSICAL_ID.items()))
    def test_each_label_maps_to_its_physical_file(self, label, expected):
        assert physical_reference_id(label) == expected

    def test_the_two_ncbi_aliases_collapse_to_one_file(self):
        assert physical_reference_id("GRCh38") == physical_reference_id("hg38_ncbi")
        assert physical_reference_id("GRCh37") == physical_reference_id("hg19_ncbi")

    def test_every_registry_label_is_covered(self):
        """A label added to the registry without a physical id is a silent wrong file."""
        assert set(list_assemblies()) == set(PHYSICAL_ID)

    def test_there_are_exactly_six_physical_files(self):
        assert len(set(PHYSICAL_ID.values())) == 6


class TestUcscFamily:
    @pytest.mark.parametrize("label,expected", sorted(UCSC_FAMILY.items()))
    def test_family_is_the_ucsc_name_of_the_coordinate_system(self, label, expected):
        assert ucsc_family(label) == expected


class TestPhysicalKeyedKinds:
    """bwa and cram: contig naming differs by source, so the source is load-bearing."""

    @pytest.mark.parametrize("kind,prefix", [("bwa", "bwa_reference"), ("cram", "cram_reference")])
    def test_exact_key_first_then_ucsc_family(self, kind, prefix):
        assert reference_keys(kind, "hg38_ensembl") == (
            f"{prefix}_hg38_ensembl",
            f"{prefix}_hg38",
        )

    @pytest.mark.parametrize("kind,prefix", [("bwa", "bwa_reference"), ("cram", "cram_reference")])
    def test_the_existing_label_override_contract_survives(self, kind, prefix):
        """test_reference_resolution.py pins that a config may specialise hg19_ncbi."""
        assert reference_keys(kind, "hg19_ncbi")[0] == f"{prefix}_hg19_ncbi"

    @pytest.mark.parametrize("label", ["GRCh38", "hg38_ncbi"])
    def test_both_ncbi_labels_reach_the_same_physical_key(self, label):
        assert "bwa_reference_GRCh38" in reference_keys("bwa", label)

    def test_a_label_specific_key_is_offered_before_the_physical_one(self):
        assert reference_keys("bwa", "hg38_ncbi") == (
            "bwa_reference_hg38_ncbi", "bwa_reference_GRCh38", "bwa_reference_hg38",
        )

    def test_a_ucsc_label_yields_a_single_key_not_a_duplicate_pair(self):
        assert reference_keys("bwa", "hg38") == ("bwa_reference_hg38",)


class TestCoordinateKeyedKinds:
    """advntr and shark: only two databases and two regions exist, by coordinate system."""

    @pytest.mark.parametrize(
        "label,expected",
        [(label, UCSC_FAMILY[label]) for label in sorted(UCSC_FAMILY)],
    )
    def test_advntr_key_follows_the_coordinate_system(self, label, expected):
        assert reference_keys("advntr", label) == (f"advntr_reference_vntr_{expected}",)

    @pytest.mark.parametrize(
        "label,expected",
        [(label, UCSC_FAMILY[label]) for label in sorted(UCSC_FAMILY)],
    )
    def test_shark_key_follows_the_coordinate_system(self, label, expected):
        assert reference_keys("shark", label) == (f"muc1_region_fasta_{expected}",)

    @pytest.mark.parametrize("label", ["hg38_ncbi", "hg38_ensembl"])
    def test_the_four_labels_that_used_to_fall_through_now_reach_grch38(self, label):
        """pipeline.py's raw dict sent these to the hg19 database. Regression for F2."""
        assert reference_keys("advntr", label) == ("advntr_reference_vntr_hg38",)


class TestRejection:
    def test_an_unknown_kind_is_rejected_by_name(self):
        with pytest.raises(ValueError, match="minimap2"):
            reference_keys("minimap2", "hg19")

    def test_an_unknown_assembly_is_rejected_with_the_supported_list(self):
        with pytest.raises(ValueError, match="hg38_ensembl"):
            reference_keys("bwa", "b37")
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
export PATH=~/miniforge3/envs/vntyper/bin:$PATH
pytest tests/unit/test_reference_keys.py -v
```

Expected: collection error, `ImportError: cannot import name 'physical_reference_id'`.

- [ ] **Step 3: Write the implementation**

Append to `vntyper/scripts/reference_registry.py`, immediately before the
`# Validation Functions` banner at `:496`:

```python
# =============================================================================
# Reference Key Resolution
# =============================================================================

# How each kind of reference file is keyed in `config["reference_data"]`.
#
# `physical`  - the file differs by reference source, because contig naming does
#               (`chr1` vs `1` vs `NC_000001.11`). Eight labels, six files.
# `coordinate_system` - only two files exist, one per coordinate system: the adVNTR
#               databases `hg19_muc1.db`/`hg38_muc1.db` and the two MUC1 region
#               FASTAs. Source naming is irrelevant to both.
REFERENCE_KINDS: dict[str, dict[str, str]] = {
    "bwa": {"prefix": "bwa_reference", "keyed_by": "physical"},
    "cram": {"prefix": "cram_reference", "keyed_by": "physical"},
    "advntr": {"prefix": "advntr_reference_vntr", "keyed_by": "coordinate_system"},
    "shark": {"prefix": "muc1_region_fasta", "keyed_by": "coordinate_system"},
}

# The physical file an (coordinate system, source) pair names. These ids are also the
# entry names in install_references_config.json, which is what lets the writer and the
# readers agree without either of them owning the name.
_PHYSICAL_IDS: dict[tuple[str, str], str] = {
    ("GRCh37", "ucsc"): "hg19",
    ("GRCh38", "ucsc"): "hg38",
    ("GRCh37", "ncbi"): "GRCh37",
    ("GRCh38", "ncbi"): "GRCh38",
    ("GRCh37", "ensembl"): "hg19_ensembl",
    ("GRCh38", "ensembl"): "hg38_ensembl",
}

_UCSC_BY_COORDINATE_SYSTEM = {"GRCh37": "hg19", "GRCh38": "hg38"}


def physical_reference_id(assembly_name: str) -> str:
    """Collapse an assembly label onto the physical file it names.

    `GRCh38` and `hg38_ncbi` are two labels for one NCBI file, as are `GRCh37` and
    `hg19_ncbi`. Keying on the label instead of the file means the writer emits
    `bwa_reference_GRCh38` while a `hg38_ncbi` run looks for
    `bwa_reference_hg38_ncbi`, misses, and silently uses UCSC sequence.

    Args:
        assembly_name (str): Assembly label (canonical or alias).

    Returns:
        str: One of `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`, `hg38_ensembl`.

    Raises:
        ValueError: If the assembly is unknown.

    Examples:
        >>> physical_reference_id("hg38_ncbi")
        'GRCh38'
        >>> physical_reference_id("hg38_ensembl")
        'hg38_ensembl'
    """
    identity = (get_coordinate_system(assembly_name), get_reference_source(assembly_name))
    physical = _PHYSICAL_IDS.get(identity)
    if physical is None:  # pragma: no cover - unreachable while validate_registry passes
        raise ValueError(f"No physical reference id for {identity}")
    return physical


def ucsc_family(assembly_name: str) -> str:
    """Return the UCSC-named assembly sharing this label's coordinate system.

    Args:
        assembly_name (str): Assembly label (canonical or alias).

    Returns:
        str: `hg19` or `hg38`.

    Raises:
        ValueError: If the assembly is unknown.

    Examples:
        >>> ucsc_family("hg38_ensembl")
        'hg38'
    """
    return _UCSC_BY_COORDINATE_SYSTEM[get_coordinate_system(assembly_name)]


def reference_keys(kind: str, assembly_name: str) -> tuple[str, ...]:
    """Config keys naming this reference, most specific first.

    Physical-keyed kinds return up to three keys - the label key, the physical key, then
    the UCSC-family key - de-duplicated with order preserved, so a UCSC label collapses to
    one. The label tier preserves an existing tested capability: a replacement config may
    specialise or disable one accepted label. Coordinate-keyed kinds return a single key,
    whose name is unchanged from before this milestone.

    Callers must resolve by **membership**, not truthiness: a key that is present with
    value `None` is an authoritative "disabled" and must not fall through. See
    `reference_resolution.resolve_from_mapping`.

    Args:
        kind (str): One of `REFERENCE_KINDS`.
        assembly_name (str): Assembly label (canonical or alias).

    Returns:
        tuple[str, ...]: One or two `config["reference_data"]` keys.

    Raises:
        ValueError: If the kind or the assembly is unknown.

    Examples:
        >>> reference_keys("bwa", "hg38_ensembl")
        ('bwa_reference_hg38_ensembl', 'bwa_reference_hg38')
        >>> reference_keys("bwa", "hg38_ncbi")
        ('bwa_reference_hg38_ncbi', 'bwa_reference_GRCh38', 'bwa_reference_hg38')
        >>> reference_keys("advntr", "hg38_ncbi")
        ('advntr_reference_vntr_hg38',)
    """
    spec = REFERENCE_KINDS.get(kind)
    if spec is None:
        supported = ", ".join(sorted(REFERENCE_KINDS))
        raise ValueError(f"Unknown reference kind '{kind}'. Supported kinds: {supported}")

    prefix = spec["prefix"]
    family_key = f"{prefix}_{ucsc_family(assembly_name)}"
    if spec["keyed_by"] == "coordinate_system":
        return (family_key,)

    label = normalize_assembly_name(assembly_name, warn_deprecated=False)
    ordered = (
        f"{prefix}_{label}",
        f"{prefix}_{physical_reference_id(assembly_name)}",
        family_key,
    )
    return tuple(dict.fromkeys(ordered))
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
pytest tests/unit/test_reference_keys.py -v
```

Expected: all PASS.

- [ ] **Step 5: Confirm the registry still validates and nothing else moved**

```bash
pytest tests/unit/test_reference_registry.py tests/unit/test_reference_resolution.py -v
ruff check vntyper/scripts/reference_registry.py && mypy vntyper/scripts/reference_registry.py
```

Expected: all PASS, no lint or type errors.

- [ ] **Step 6: Commit**

```bash
git add vntyper/scripts/reference_registry.py tests/unit/test_reference_keys.py
git commit -m "feat(registry): resolve reference config keys from physical identity

Eight assembly labels denote six physical files - GRCh37/hg19_ncbi are one NCBI
file and GRCh38/hg38_ncbi another - so a label-keyed lookup can never resolve the
_ncbi aliases and silently falls back to UCSC sequence. Key on
(coordinate_system, reference_source) instead; the six resulting ids are already
the entry names in install_references_config.json.

Refs #163

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 2: Membership-based resolution helper  *(PR-2)*

**Files:**
- Modify: `vntyper/scripts/reference_resolution.py`
- Test: `tests/unit/test_reference_resolution.py` (extend)

**Interfaces:**
- Consumes: `reference_keys` from Task 1.
- Produces:
  - `ResolvedReference` — frozen dataclass with fields `key: str`, `value: str | None`, `is_fallback: bool`
  - `resolve_from_mapping(kind: str, assembly: str, mapping: dict) -> ResolvedReference | None`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_reference_resolution.py`:

```python
class TestResolveFromMapping:
    """Membership, not truthiness - the rule `configured_reference_candidates` already uses.

    A key present with value None is a deliberate "disabled". Falling through it would
    silently re-enable a reference an operator switched off.
    """

    def test_the_physical_key_wins_when_present(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        mapping = {"bwa_reference_hg38_ensembl": "/refs/ensembl.fa", "bwa_reference_hg38": "/refs/ucsc.fa"}
        resolved = resolve_from_mapping("bwa", "hg38_ensembl", mapping)
        assert resolved.key == "bwa_reference_hg38_ensembl"
        assert resolved.value == "/refs/ensembl.fa"
        assert resolved.is_fallback is False

    def test_the_family_key_is_used_and_flagged_when_the_physical_key_is_absent(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        resolved = resolve_from_mapping("bwa", "hg38_ensembl", {"bwa_reference_hg38": "/refs/ucsc.fa"})
        assert resolved.key == "bwa_reference_hg38"
        assert resolved.value == "/refs/ucsc.fa"
        assert resolved.is_fallback is True

    def test_a_present_null_is_authoritative_and_does_not_fall_through(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        mapping = {"bwa_reference_hg38_ensembl": None, "bwa_reference_hg38": "/refs/ucsc.fa"}
        resolved = resolve_from_mapping("bwa", "hg38_ensembl", mapping)
        assert resolved.key == "bwa_reference_hg38_ensembl"
        assert resolved.value is None
        assert resolved.is_fallback is False

    def test_no_key_present_at_all_returns_none(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        assert resolve_from_mapping("bwa", "hg38_ensembl", {}) is None

    def test_both_ncbi_labels_resolve_the_same_entry(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        mapping = {"bwa_reference_GRCh38": "/refs/ncbi.fna"}
        for label in ("GRCh38", "hg38_ncbi"):
            resolved = resolve_from_mapping("bwa", label, mapping)
            assert resolved.value == "/refs/ncbi.fna", label
            assert resolved.is_fallback is False, label
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
pytest tests/unit/test_reference_resolution.py::TestResolveFromMapping -v
```

Expected: FAIL, `ImportError: cannot import name 'resolve_from_mapping'`.

- [ ] **Step 3: Write the implementation**

In `vntyper/scripts/reference_resolution.py`, extend the imports and add the helper above
`configured_reference_candidates`:

```python
from dataclasses import dataclass

from vntyper.scripts.reference_registry import (
    get_coordinate_system,
    get_reference_source,
    list_assemblies,
    reference_keys,
)


@dataclass(frozen=True)
class ResolvedReference:
    """Which config key supplied a reference, and whether it was a fallback.

    Attributes:
        key: The `reference_data` key that was present.
        value: Its value. `None` means present-and-disabled, not missing.
        is_fallback: True when the UCSC-family key stood in for an absent physical
            key, which means the run uses UCSC sequence for a non-UCSC request.
    """

    key: str
    value: str | None
    is_fallback: bool


def resolve_from_mapping(kind: str, reference_assembly: str, mapping: dict) -> ResolvedReference | None:
    """Resolve a reference from a config mapping by key membership.

    Args:
        kind: One of `reference_registry.REFERENCE_KINDS`.
        reference_assembly: Supported assembly label.
        mapping: A dict of config keys to values, e.g. `config["reference_data"]`.

    Returns:
        The first key that is *present* in the mapping, or None when none are.
        Presence wins over truthiness so an explicit null stays authoritative.

    Raises:
        ValueError: If the kind or the assembly is unknown.
    """
    keys = reference_keys(kind, reference_assembly)
    for index, key in enumerate(keys):
        if key in mapping:
            return ResolvedReference(key=key, value=mapping[key], is_fallback=index > 0)
    return None
```

- [ ] **Step 4: Delegate `configured_reference_candidates` to it**

Replace the body of `configured_reference_candidates` (`:36-52`) with:

```python
    reference_data = config.get("reference_data", {})
    values: list[tuple[str, Any]] = []
    for source, kind in (("config_cram_reference", "cram"), ("config_bwa_reference", "bwa")):
        resolved = resolve_from_mapping(kind, reference_assembly, reference_data)
        values.append((source, resolved.value if resolved is not None else None))
    return tuple(values)
```

The `get_reference_source` / `list_assemblies` imports become unused here — remove them if
`ruff` flags them, but keep `get_coordinate_system` if other functions still use it.

**Behaviour note.** This is *not* a pure refactor once Task 1's label tier is in place: it
strictly widens what resolves. `test_label_specific_reference_keys_override_the_family_fallback_including_null`
must still pass unchanged — it is the reason the label tier exists. Run it explicitly and do
not edit it.

- [ ] **Step 5: Run the full existing suite for this module**

```bash
pytest tests/unit/test_reference_resolution.py tests/unit/test_alignment_preflight_reference.py \
       tests/unit/test_reference_binding.py tests/unit/test_reference_cache_binding.py -v
```

Expected: all PASS **without editing any existing test**. If
`test_label_specific_reference_keys_override_the_family_fallback_including_null` goes red,
the label tier in Task 1 is wrong — fix Task 1, never the test.

- [ ] **Step 6: Commit**

```bash
git add vntyper/scripts/reference_resolution.py tests/unit/test_reference_resolution.py
git commit -m "refactor(references): one membership-based resolver for all reference kinds

configured_reference_candidates computed the UCSC family inline; it now delegates
to reference_keys. Behaviour is unchanged - the existing tests pinning that a
present-but-null key is authoritative stay green - but the rule now has one
implementation that the remaining readers can share.

Refs #163

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 7: Route the BWA reader and the logging-safety reader

**Files:**
- Modify: `vntyper/scripts/cli_handlers.py:261-275`
- Modify: `vntyper/scripts/cli_logging_safety.py:43-51`
- Test: `tests/unit/test_cli_handlers_bwa_selection.py` (create), `tests/unit/test_cli_logging_safety.py` (extend)

**Interfaces:**
- Consumes: `resolve_from_mapping`, `ResolvedReference` (Task 2).
- Produces, both in `cli_handlers`:
  - `select_bwa_reference(config: dict, reference_assembly: str, *, required: bool = True) -> str | None`.
    `required=True` (the run path) raises `ValueError` when nothing resolves;
    `required=False` (the safety guard) returns `None`. An **unknown assembly** raises in
    both cases and must propagate. A custom exception class would have been the obvious
    way to separate those, but `AGENTS.md:140` forbids one — hence the keyword.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_cli_handlers_bwa_selection.py`:

```python
"""BWA reference selection, and the safety guard that shares it.

Two readers used to collapse every assembly onto bwa_reference_hg19/hg38 independently:
cli_handlers for the run, and cli_logging_safety._selected_bwa_reference for the check
that stops --log-file naming an operator input. Because the guard runs before
setup_logging opens the log in append mode, a guard looking at the wrong file means
logging can append into a reference FASTA.
"""

import pytest

from vntyper.scripts.cli_handlers import select_bwa_reference

pytestmark = pytest.mark.unit


def config_with(**reference_data):
    return {"reference_data": dict(reference_data)}


class TestSelection:
    def test_the_physical_key_is_preferred(self):
        cfg = config_with(bwa_reference_hg38_ensembl="/refs/ens.fa", bwa_reference_hg38="/refs/ucsc.fa")
        assert select_bwa_reference(cfg, "hg38_ensembl") == "/refs/ens.fa"

    @pytest.mark.parametrize("label", ["GRCh38", "hg38_ncbi"])
    def test_both_ncbi_labels_select_the_ncbi_reference(self, label):
        cfg = config_with(bwa_reference_GRCh38="/refs/ncbi.fna", bwa_reference_hg38="/refs/ucsc.fa")
        assert select_bwa_reference(cfg, label) == "/refs/ncbi.fna"

    def test_the_ucsc_family_fallback_warns_and_names_both_keys(self, caplog):
        cfg = config_with(bwa_reference_hg38="/refs/ucsc.fa")
        with caplog.at_level("WARNING"):
            assert select_bwa_reference(cfg, "hg38_ensembl") == "/refs/ucsc.fa"
        message = " ".join(record.getMessage() for record in caplog.records)
        assert "hg38_ensembl" in message and "bwa_reference_hg38" in message and "ucsc" in message.lower()

    def test_an_exact_null_is_authoritative_and_fails_closed(self):
        cfg = config_with(bwa_reference_hg38_ensembl=None, bwa_reference_hg38="/refs/ucsc.fa")
        with pytest.raises(ValueError, match="hg38_ensembl"):
            select_bwa_reference(cfg, "hg38_ensembl")

    def test_nothing_configured_fails_closed_naming_every_key_tried(self):
        with pytest.raises(ValueError) as excinfo:
            select_bwa_reference(config_with(), "hg38_ensembl")
        assert "bwa_reference_hg38_ensembl" in str(excinfo.value)
        assert "bwa_reference_hg38" in str(excinfo.value)
```

Append to `tests/unit/test_cli_logging_safety.py`:

```python
class TestTheGuardUsesTheSameReferenceTheRunWill:
    """A guard that inspects a different file than the run uses is not a guard."""

    @pytest.mark.parametrize(
        "label,key",
        [
            ("hg38_ensembl", "bwa_reference_hg38_ensembl"),
            ("hg38_ncbi", "bwa_reference_GRCh38"),
            ("GRCh37", "bwa_reference_GRCh37"),
            ("hg19_ensembl", "bwa_reference_hg19_ensembl"),
        ],
    )
    def test_the_exact_reference_cannot_be_used_as_a_log_file(self, tmp_path, label, key):
        from vntyper.scripts.cli_logging_safety import validate_pipeline_log_destination

        reference = tmp_path / "exact.fa"
        reference.write_text(">chr1\nACGT\n")
        config = {"reference_data": {key: str(reference), "bwa_reference_hg38": str(tmp_path / "ucsc.fa")}}
        args = argparse.Namespace(
            fastq1="r1.fq", fastq2="r2.fq", bam=None, cram=None, bed_file=None,
            reference_fasta=None, reference_assembly=label, log_file=str(reference),
        )
        with pytest.raises(ValueError, match="exact.fa"):
            validate_pipeline_log_destination(str(reference), args, config)

    @pytest.mark.parametrize("suffix", [".amb", ".ann", ".bwt", ".pac", ".sa"])
    def test_a_bwa_sidecar_of_the_exact_reference_cannot_be_used_either(self, tmp_path, suffix):
        from vntyper.scripts.cli_logging_safety import validate_pipeline_log_destination

        reference = tmp_path / "exact.fa"
        reference.write_text(">chr1\nACGT\n")
        sidecar = tmp_path / f"exact.fa{suffix}"
        sidecar.write_bytes(b"\x00")
        config = {"reference_data": {"bwa_reference_hg38_ensembl": str(reference)}}
        args = argparse.Namespace(
            fastq1="r1.fq", fastq2="r2.fq", bam=None, cram=None, bed_file=None,
            reference_fasta=None, reference_assembly="hg38_ensembl", log_file=str(sidecar),
        )
        with pytest.raises(ValueError):
            validate_pipeline_log_destination(str(sidecar), args, config)
```

**The guard's real entry point is `validate_pipeline_log_destination(log_file, args, config)`
(`cli_logging_safety.py:104-108`)** — use that signature, not an invented one, and
`import argparse` at the top of the test module.

- [ ] **Step 2: Run the tests to verify they fail**

```bash
pytest tests/unit/test_cli_handlers_bwa_selection.py tests/unit/test_cli_logging_safety.py -v
```

Expected: `ImportError: cannot import name 'select_bwa_reference'`, and the sidecar and
exact-key cases failing because the guard inspects the UCSC file.

- [ ] **Step 3: Implement `select_bwa_reference`**

In `vntyper/scripts/cli_handlers.py`, replace `:261-275` with a call to a new module-level
function:

First add the imports `cli_handlers.py` does not have (`:22-34` imports none of these):

```python
from vntyper.scripts.reference_registry import get_reference_source, reference_keys
from vntyper.scripts.reference_resolution import resolve_from_mapping
```

```python
def select_bwa_reference(config: dict[str, Any], reference_assembly: str, *, required: bool = True) -> str | None:
    """Resolve the BWA reference for an assembly, failing closed.

    Both this and `cli_logging_safety` must agree, or the guard that refuses to let
    `--log-file` name an operator input inspects a different file from the one the run
    opens for writing.

    Args:
        config: Pipeline configuration.
        reference_assembly: Supported assembly label.

    Returns:
        str: Path to the reference FASTA.

    Raises:
        ValueError: If no configured key resolves and `required` is True.
        ValueError: If `reference_assembly` is not a supported label. This must NOT be
            swallowed by callers - an unknown assembly is a configuration error, not a
            missing file.
    """
    resolved = resolve_from_mapping("bwa", reference_assembly, config.get("reference_data", {}))
    keys = ", ".join(reference_keys("bwa", reference_assembly))
    if resolved is None or not resolved.value:
        if not required:
            return None
        message = (
            f"No BWA reference configured for --reference-assembly {reference_assembly!r}. "
            f"Tried: {keys}. Run `vntyper install-references` or set one of those keys."
        )
        logger.error(message)
        raise ValueError(message)
    if resolved.is_fallback:
        # Name the effective source plainly. Deriving it from the key suffix would be
        # fragile; the fallback key is always the UCSC family key, so the source is ucsc.
        logger.warning(
            f"--reference-assembly {reference_assembly!r} has no {keys.split(', ')[0]!r} entry; "
            f"falling back to {resolved.key!r}. This run therefore uses 'ucsc' sequence, "
            f"not {get_reference_source(reference_assembly)!r}. "
            f"Run `vntyper install-references` to install the requested reference."
        )
    return resolved.value
```

Call it from `handle_run` where the old block was, and pass `bwa_reference` on unchanged.

- [ ] **Step 4: Route the safety guard through it**

In `vntyper/scripts/cli_logging_safety.py`, replace `_selected_bwa_reference`'s body
(`:43-51`) with:

```python
def _selected_bwa_reference(args: argparse.Namespace, config: dict[str, Any]) -> Path | None:
    """The reference the run will actually use, resolved exactly as `cli_handlers` does.

    Returns None when nothing resolves; the caller's other checks still apply and the
    run itself will fail closed later with a clearer message.
    """
    from vntyper.scripts.cli_handlers import select_bwa_reference

    assembly = args.reference_assembly or config.get("default_values", {}).get("reference_assembly", "hg19")
    # required=False turns "nothing configured" into None. An unknown assembly still
    # raises and MUST propagate: swallowing it fails the guard open, and the guard runs
    # before setup_logging opens the log file in append mode.
    selected = select_bwa_reference(config, assembly, required=False)
    return Path(selected) if selected else None
```

- [ ] **Step 5: Run the tests to verify they pass**

```bash
pytest tests/unit/test_cli_handlers_bwa_selection.py tests/unit/test_cli_logging_safety.py \
       tests/unit/test_cli_handlers.py -v
```

Expected: all PASS.

- [ ] **Step 6: Record the effective reference in the run summary**

Adding keys to a dict does not put them in the report — `generate_report.py` and its
templates select known fields explicitly, and nothing consumes these names today. Three
places must change together, and all three belong in this task's commit:

1. `summary.py:26-42`, where the top-level mapping is built — add
   `reference_assembly_requested`, `reference_key_used`, `reference_path` and
   `reference_source_effective`. `select_bwa_reference` returns only a path, so have it
   also return the `ResolvedReference` (or expose a sibling that does) rather than
   re-deriving the key at the call site.
2. `generate_report.py` — pass them into the template context.
3. The report template — render them, in the block that already shows run parameters.

Then extend `tests/unit/test_generate_report.py` with a case asserting all four appear
when the fallback was taken, and a case asserting the requested and effective sources are
shown as different.

- [ ] **Step 7: Commit**

```bash
git add vntyper/scripts/cli_handlers.py vntyper/scripts/cli_logging_safety.py \
        vntyper/scripts/pipeline.py tests/unit/test_cli_handlers_bwa_selection.py \
        tests/unit/test_cli_logging_safety.py tests/unit/test_generate_report.py
git commit -m "fix(references): one BWA reference resolver for the run and the safety guard

cli_logging_safety._selected_bwa_reference collapsed every assembly onto
bwa_reference_hg19/hg38 independently of cli_handlers, so with an exact
assembly-specific key configured the guard inspected the UCSC file while the run
used the exact one - and --log-file could then append into that reference or its
BWA sidecars, because the guard runs before setup_logging opens the file. Both now
call select_bwa_reference, which fails closed and warns when a UCSC-family
fallback degrades the requested source.

Refs #163

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 8: Route the adVNTR reader

**Files:**
- Modify: `vntyper/scripts/pipeline.py:484-505`
- Test: `tests/unit/test_pipeline_advntr_reference.py` (create)

**Interfaces:**
- Consumes: `resolve_from_mapping` (Task 2).
- Produces: no new public names; `pipeline.py` behaviour only.

- [ ] **Step 1: Write the failing test**

```python
"""Which adVNTR database a run uses, for every assembly label.

pipeline.py used a four-entry dict with `.get(label, "hg19")`. The four
source-qualified labels were absent, so hg38_ncbi and hg38_ensembl silently loaded the
GRCh37 database - the wrong coordinate system, with no warning. Genotype-affecting.
"""

import pytest

from vntyper.scripts.reference_registry import list_assemblies
from vntyper.scripts.reference_resolution import resolve_from_mapping

pytestmark = pytest.mark.unit

EXPECTED = {
    "hg19": "hg19", "GRCh37": "hg19", "hg19_ncbi": "hg19", "hg19_ensembl": "hg19",
    "hg38": "hg38", "GRCh38": "hg38", "hg38_ncbi": "hg38", "hg38_ensembl": "hg38",
}


@pytest.mark.parametrize("label", sorted(EXPECTED))
def test_every_label_reaches_its_own_coordinate_systems_database(label):
    mapping = {
        "advntr_reference_vntr_hg19": "/refs/vntr_db_advntr/hg19_muc1.db",
        "advntr_reference_vntr_hg38": "/refs/vntr_db_advntr/hg38_muc1.db",
    }
    resolved = resolve_from_mapping("advntr", label, mapping)
    assert resolved.value.endswith(f"{EXPECTED[label]}_muc1.db")


@pytest.mark.parametrize("label", ["hg38_ncbi", "hg38_ensembl"])
def test_the_two_labels_that_used_to_load_the_wrong_database(label):
    mapping = {"advntr_reference_vntr_hg19": "/hg19.db", "advntr_reference_vntr_hg38": "/hg38.db"}
    assert resolve_from_mapping("advntr", label, mapping).value == "/hg38.db"


def test_no_label_is_left_unmapped():
    assert set(list_assemblies()) == set(EXPECTED)
```

**These tests pass without touching `pipeline.py`, which makes them worthless on their
own.** They exercise `resolve_from_mapping`, not the production branch at
`pipeline.py:480-505`. Extract the decision into a named function and test *that*:

```python
def select_advntr_reference(config: dict, reference_assembly: str) -> str | None:
    """The adVNTR database for an assembly, by coordinate system.

    Args:
        config: Pipeline configuration.
        reference_assembly: Supported assembly label.

    Returns:
        str | None: Database path, or None when no key is configured.
    """
    resolved = resolve_from_mapping("advntr", reference_assembly, config.get("reference_data", {}))
    return resolved.value if resolved is not None else None
```

and add, in the same test module:

```python
@pytest.mark.parametrize("label", sorted(EXPECTED))
def test_the_pipeline_itself_selects_the_right_database(label):
    from vntyper.scripts.pipeline import select_advntr_reference

    config = {"reference_data": {
        "advntr_reference_vntr_hg19": "/refs/hg19_muc1.db",
        "advntr_reference_vntr_hg38": "/refs/hg38_muc1.db",
    }}
    assert select_advntr_reference(config, label).endswith(f"{EXPECTED[label]}_muc1.db")
```

- [ ] **Step 2: Run the tests and watch the pipeline one fail**

```bash
pytest tests/unit/test_pipeline_advntr_reference.py -v
```

Expected: the `resolve_from_mapping` cases PASS (they pin the resolver); the
`select_advntr_reference` case FAILS with `ImportError`.

- [ ] **Step 3: Replace the raw dict**

In `vntyper/scripts/pipeline.py`, add `select_advntr_reference` at module scope (body
above) and replace `:485-493` with:

```python
                advntr_reference = select_advntr_reference(config, reference_assembly)
```

Leave the explicit-`advntr_reference` branch at `:494-501` and the fail-closed check at
`:503-505` unchanged.

- [ ] **Step 4: Prove the old behaviour is gone**

```bash
grep -n 'ref_map' vntyper/scripts/pipeline.py
pytest tests/unit/ -k advntr -v
```

Expected: no `ref_map` remains; adVNTR tests PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/pipeline.py tests/unit/test_pipeline_advntr_reference.py
git commit -m "fix(advntr): select the database by coordinate system, not a partial dict

The four-entry map defaulted unknown labels to hg19, so --reference-assembly
hg38_ncbi and hg38_ensembl silently loaded the GRCh37 database. Route through the
registry resolver, which covers all eight labels.

Fixes the adVNTR half of #163

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 9: SHARK assembly-keyed selection

**Files:**
- Modify: `vntyper/modules/shark/shark_filtering.py`
- Modify: `vntyper/modules/shark/shark_config.json`
- Test: `tests/unit/test_shark_filtering.py` (rewrite two classes, add two)

**Interfaces:**
- Consumes: `resolve_from_mapping` (Task 2).
- Produces: `select_muc1_region_fasta(config: dict, main_config: dict, reference_assembly: str) -> str`

- [ ] **Step 1: Rewrite the two characterisation classes as selection tests**

In `tests/unit/test_shark_filtering.py`, update `SHARK_CONFIG` at `:31`:

```python
SHARK_CONFIG = {
    "shark_settings": {
        "muc1_region_fasta_hg19": "reference/muc1_region_hg19.fa",
        "muc1_region_fasta_hg38": "reference/muc1_region_hg38.fa",
    }
}
```

Replace `TestReferenceAssemblyIsAccceptedAndIgnored` (`:179`) and
`TestNonHg19AssemblyLogsAWarning` (`:236`) with:

```python
class TestReferenceAssemblySelectsTheRegion:
    """#152. SHARK filters on exact k-mers, and 40.6% of the hg38 MUC1 region's
    canonical 17-mers are absent from the hg19 region - so an hg38 run against the
    hg19 reference cannot see them. Measured across seven cohort samples: the hg38
    reference retains 3.2-34.7% more reads. The parameter now selects.
    """

    @pytest.mark.parametrize(
        "assembly,expected",
        [
            ("hg19", "muc1_region_hg19.fa"), ("GRCh37", "muc1_region_hg19.fa"),
            ("hg19_ncbi", "muc1_region_hg19.fa"), ("hg19_ensembl", "muc1_region_hg19.fa"),
            ("hg38", "muc1_region_hg38.fa"), ("GRCh38", "muc1_region_hg38.fa"),
            ("hg38_ncbi", "muc1_region_hg38.fa"), ("hg38_ensembl", "muc1_region_hg38.fa"),
        ],
    )
    def test_each_assembly_gets_its_coordinate_systems_region(self, tmp_path, captured_command, assembly, expected):
        filter_with(tmp_path, reference_assembly=assembly)
        assert expected in captured_command[0]["command"]

    def test_no_warning_is_logged_for_a_supported_assembly(self, tmp_path, captured_command, caplog):
        with caplog.at_level(logging.WARNING):
            filter_with(tmp_path, reference_assembly="hg38")
        assert not [r for r in caplog.records if "reference_assembly" in r.getMessage()]


class TestTheInstalledTreeWins:
    """shark_config.json is a separate file from config.json, so install-references
    cannot update it. Without reference_data as the first layer, `install-references
    --output-dir X` would leave SHARK pointing at a CWD-relative path - #163's exact
    bug, reintroduced for SHARK.
    """

    def test_reference_data_overrides_the_shipped_shark_config(self, tmp_path, captured_command):
        main = {**MAIN_CONFIG, "reference_data": {"muc1_region_fasta_hg38": "/custom/refs/muc1_region_hg38.fa"}}
        filter_with(tmp_path, main_config=main, reference_assembly="hg38")
        assert "/custom/refs/muc1_region_hg38.fa" in captured_command[0]["command"]

    def test_a_present_null_in_reference_data_fails_closed(self, tmp_path):
        main = {**MAIN_CONFIG, "reference_data": {"muc1_region_fasta_hg38": None}}
        with pytest.raises(ValueError, match="hg38"):
            filter_with(tmp_path, main_config=main, reference_assembly="hg38")

    def test_a_structurally_legacy_config_still_works(self, tmp_path, captured_command):
        legacy = {"shark_settings": {"muc1_region_fasta": "/legacy/region.fa"}}
        filter_with(tmp_path, config=legacy, reference_assembly="hg19")
        assert "/legacy/region.fa" in captured_command[0]["command"]

    def test_an_incomplete_keyed_config_does_not_masquerade_as_legacy(self, tmp_path):
        """Only hg19 configured, plus a stray flat key: an hg38 run must not take it."""
        partial = {"shark_settings": {"muc1_region_fasta_hg19": "/a.fa", "muc1_region_fasta": "/flat.fa"}}
        with pytest.raises(ValueError, match="hg38"):
            filter_with(tmp_path, config=partial, reference_assembly="hg38")


# NOTE: this class does NOT belong in tests/unit/test_shark_filtering.py - `:28` applies
# `pytestmark = pytest.mark.unit` to the whole module, so adding an integration marker
# here gives it both and `pytest -m unit` still selects it. Put it in
# tests/integration/test_shark_case_sensitivity.py with only the integration marker.
class TestSharkIsCaseInsensitive:
    """Verified against shark 1.2.0 (h077b44d_5): its to_int table maps a/c/g/t and
    A/C/G/T identically, which is why the soft masking in the derived references is
    harmless. A future base rebuild could change that, and ordinary cohort reads need
    not contain a lowercase-only discriminating k-mer to catch it.
    """

    def test_a_lowercase_only_match_is_retained_exactly_as_its_uppercase_copy(self, tmp_path):
        import shutil
        import subprocess

        shark = shutil.which("shark")
        if shark is None:
            pytest.skip("shark binary not on PATH")

        # Every 17-mer overlapping the read is lowercase in `lower`, uppercase in `upper`.
        sequence = "acgt" * 40
        lower = tmp_path / "lower.fa"
        upper = tmp_path / "upper.fa"
        lower.write_text(f">r\n{sequence}\n")
        upper.write_text(f">r\n{sequence.upper()}\n")

        read = sequence[20:100].upper()          # reads are always uppercase
        quality = "I" * len(read)
        for name in ("r1", "r2"):
            (tmp_path / f"{name}.fastq").write_text(f"@read1\n{read}\n+\n{quality}\n")

        def retained(reference):
            out1 = tmp_path / f"{reference.stem}_1.fq"
            out2 = tmp_path / f"{reference.stem}_2.fq"
            subprocess.run(
                [shark, "-r", str(reference), "-1", str(tmp_path / "r1.fastq"),
                 "-2", str(tmp_path / "r2.fastq"), "-o", str(out1), "-p", str(out2)],
                check=True, capture_output=True,
            )
            return out1.read_text()

        assert retained(lower) == retained(upper), (
            "SHARK has become case-sensitive; the derived references are soft-masked, "
            "so both derivations now need an upper-casing step and new expected digests"
        )
```

Complete the last test's body when implementing — build the FASTQ with
`Bio.SeqIO`-free plain text (four lines per read: `@r1`, the sequence, `+`, `I` * len).

- [ ] **Step 2: Run the tests to verify they fail**

```bash
pytest tests/unit/test_shark_filtering.py -v
```

Expected: the selection tests FAIL — `shark_filtering.py:76` reads a single flat key that
`SHARK_CONFIG` no longer has, so `ValueError: muc1_region_fasta not defined`.

- [ ] **Step 3: Make the shipped config assembly-keyed**

`vntyper/modules/shark/shark_config.json`:

```json
{
    "shark_settings": {
      "muc1_region_fasta_hg19": "reference/muc1_region_hg19.fa",
      "muc1_region_fasta_hg38": "reference/muc1_region_hg38.fa"
    }
}
```

- [ ] **Step 4: Implement the layered selection**

In `vntyper/modules/shark/shark_filtering.py`, delete the warning at `:64-70` and the
docstring paragraph at `:50-54`, then replace `:76-78` with a call to:

```python
def select_muc1_region_fasta(config: dict, main_config: dict, reference_assembly: str) -> str:
    """Pick the MUC1 region FASTA for an assembly, most authoritative source first.

    1. ``main_config["reference_data"]`` - what `install-references` writes, so a
       custom ``--output-dir`` is honoured. `shark_config.json` is a separate file
       that `--config-path` never touches, so without this layer an installed tree
       would be ignored.
    2. ``config["shark_settings"]`` keyed by assembly - the shipped default.
    3. the legacy flat ``muc1_region_fasta`` key, **only** when the config carries no
       keyed entry at all. A partly-populated keyed config is incomplete, not legacy.

    Resolution is by key *membership*: a key present with value ``None`` is a
    deliberate "disabled" and does not fall through.

    Args:
        config: The shark_config dictionary.
        main_config: The main configuration dictionary.
        reference_assembly: Supported assembly label.

    Returns:
        str: Path to the region FASTA.

    Raises:
        ValueError: If nothing resolves, naming the assembly and every key tried.
    """
    from vntyper.scripts.reference_registry import reference_keys
    from vntyper.scripts.reference_resolution import resolve_from_mapping

    (key,) = reference_keys("shark", reference_assembly)
    settings = config.get("shark_settings", {})

    resolved = resolve_from_mapping("shark", reference_assembly, main_config.get("reference_data", {}))
    if resolved is not None:
        if resolved.value:
            return resolved.value
        raise ValueError(f"reference_data[{resolved.key!r}] is null; SHARK is disabled for {reference_assembly!r}")

    if key in settings:
        if settings[key]:
            return settings[key]
        raise ValueError(f"shark_settings[{key!r}] is null; SHARK is disabled for {reference_assembly!r}")

    keyed = [name for name in settings if name.startswith("muc1_region_fasta_")]
    if not keyed and settings.get("muc1_region_fasta"):
        return settings["muc1_region_fasta"]

    raise ValueError(
        f"No SHARK MUC1 region FASTA for reference_assembly {reference_assembly!r}. "
        f"Tried reference_data[{key!r}] and shark_settings[{key!r}]."
    )
```

Update `run_shark_filter`'s signature docs: `reference_assembly` now *selects*.

- [ ] **Step 5: Run the tests to verify they pass**

```bash
pytest tests/unit/test_shark_filtering.py -v
pytest tests/unit/ -k shark -v
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add vntyper/modules/shark/shark_filtering.py vntyper/modules/shark/shark_config.json \
        tests/unit/test_shark_filtering.py
git commit -m "fix(shark): select the MUC1 region FASTA by coordinate system

40.6% of the hg38 MUC1 region's canonical 17-mers are absent from the hg19 region,
so --reference-assembly hg38 --extra-modules shark filtered against a Bloom filter
that could not see them. Measured on seven cohort samples: the hg38 reference
retains 3.2-34.7% more reads. reference_data is the first lookup layer because
shark_config.json is a separate file that install-references cannot update.

Fixes #152

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 10: The writer emits canonical keys, and the bundle is pinned

**Files:**
- Modify: `vntyper/scripts/install_references.py` (`update_config`, `main`)
- Modify: `vntyper/scripts/install_references_config.json`
- Modify: `vntyper/config.json`
- Test: `tests/unit/test_install_references_config_writing.py`, `tests/unit/test_reference_key_agreement.py`

**Interfaces:**
- Consumes: `reference_keys` (Task 1), `reference_bundle` (Task 3), `refs-v1` digests (Task 6).
- Produces: `canonical_reference_keys(install_config: dict, output_dir: Path) -> dict[str, Path]`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_install_references_config_writing.py`:

```python
"""What `install-references --config-path` writes into config.json.

It used to write only ucsc_*/ncbi_*/ensembl_*/vntyper_*/own_repo_* keys, which nothing
in vntyper/ reads, while the seven keys the pipeline does read kept their shipped
relative paths. After an install into a custom --output-dir the run then died at
pipeline.py:154-156. Two of the old keys also named files that cannot be used: the
.gz rather than the extracted FASTA, and the .zip rather than the two .db files.
"""

import json

import pytest

from vntyper.scripts.install_references import canonical_reference_keys

pytestmark = pytest.mark.unit


class TestCanonicalKeys:
    def test_a_genome_entry_writes_the_key_the_pipeline_reads(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert "bwa_reference_hg38" in keys
        assert "ucsc_hg38" not in keys

    def test_the_written_path_is_the_extracted_fasta_not_the_gz(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert str(keys["bwa_reference_hg38"]).endswith("chr1.hg38.fa")

    def test_the_advntr_keys_name_the_databases_not_the_zip(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert str(keys["advntr_reference_vntr_hg19"]).endswith("vntr_db_advntr/hg19_muc1.db")
        assert str(keys["advntr_reference_vntr_hg38"]).endswith("vntr_db_advntr/hg38_muc1.db")

    def test_both_shark_regions_are_written(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert set(keys) >= {"muc1_region_fasta_hg19", "muc1_region_fasta_hg38"}

    def test_every_written_key_is_one_the_registry_knows(self, tmp_path, install_config):
        from vntyper.scripts.reference_registry import REFERENCE_KINDS, list_assemblies, reference_keys

        known = {k for kind in REFERENCE_KINDS for a in list_assemblies() for k in reference_keys(kind, a)}
        known |= {"muc1_reference_vntr", "code_adVNTR_RUs", "muc1_motifs_rev_com"}
        assert set(canonical_reference_keys(install_config, tmp_path)) <= known


class TestUpdateConfigIsAtomic:
    def test_a_write_failure_leaves_the_previous_config_intact(self, tmp_path, monkeypatch):
        from vntyper.scripts import install_references

        config_path = tmp_path / "config.json"
        config_path.write_text(json.dumps({"reference_data": {"bwa_reference_hg19": "old"}}))
        monkeypatch.setattr(install_references.json, "dump", lambda *a, **k: (_ for _ in ()).throw(OSError("disk full")))
        with pytest.raises((OSError, SystemExit)):
            install_references.update_config(config_path, {"bwa_reference_hg19": "new"})
        assert json.loads(config_path.read_text())["reference_data"]["bwa_reference_hg19"] == "old"
```

Add an `install_config` fixture to `tests/unit/conftest.py` that loads the real
`vntyper/scripts/install_references_config.json` — the tests must exercise the shipped
config, not a hand-written stub, or a drift between them goes unnoticed.

Create `tests/unit/test_reference_key_agreement.py`:

```python
"""The writer's output must be resolvable by every reader, for every label.

A test that only compares expected key *tuples* would have passed a design that keys
on the assembly label - which can never resolve hg19_ncbi/hg38_ncbi, because those are
aliases of GRCh37/GRCh38 and the writer only ever emits the alias-free name.
"""

import pytest

from vntyper.scripts.install_references import canonical_reference_keys
from vntyper.scripts.reference_registry import list_assemblies
from vntyper.scripts.reference_resolution import resolve_from_mapping

pytestmark = pytest.mark.unit

KINDS = ["bwa", "advntr", "shark"]


@pytest.mark.parametrize("label", sorted(list_assemblies()))
@pytest.mark.parametrize("kind", KINDS)
def test_a_full_installation_resolves_without_a_fallback(tmp_path, install_config, label, kind):
    written = {key: str(path) for key, path in canonical_reference_keys(install_config, tmp_path).items()}
    resolved = resolve_from_mapping(kind, label, written)
    assert resolved is not None, f"{kind}/{label} resolves to nothing after a full install"
    assert resolved.value, f"{kind}/{label} resolved to an empty value"
    assert not resolved.is_fallback, (
        f"{kind}/{label} fell back to {resolved.key} after a full install - "
        "a complete installation must never degrade the requested reference source"
    )
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
pytest tests/unit/test_install_references_config_writing.py tests/unit/test_reference_key_agreement.py -v
```

Expected: `ImportError: cannot import name 'canonical_reference_keys'`.

- [ ] **Step 3: Rewrite `install_references_config.json`**

Give every genome entry its physical id as its name, a `kind`, an `installed_path` naming
the **extracted** file, the bundle `asset` and its `asset_sha256` from
`/tmp/refs-v1-asset-digests.txt`, and the `source_url`/`source_sha256` pair for
`--from-source`. Add the top-level bundle pointer:

```json
  "bundle": {
    "repository": "berntpopp/vntyper-data",
    "release_tag": "refs-v1",
    "common_asset": "vntyper-references-refs-v1-muc1.tar.gz",
    "common_asset_sha256": "..."
  },
```

`common_asset` is **not** selectable — it carries the MUC1 FASTAs and both adVNTR
databases and installs on every run. Today `own_repository_references` is unconditional
(`:771`) but `vntyper_references` is filtered by `--references` (`:670`), so
`--references hg19` — what `scheduled-tests.yml:103-110` runs — installs no adVNTR
database at all. Making the common asset implicit fixes that.

- [ ] **Step 4: Implement `canonical_reference_keys` and use it in `main`**

Add the module-scope imports `install_references.py` lacks — `os` (for the atomic
replace) and the two registry names:

```python
import os

from vntyper.scripts.reference_registry import REFERENCE_KINDS, reference_keys
```

```python
def canonical_reference_keys(install_config: dict[str, Any], output_dir: Path) -> dict[str, Path]:
    """Map every installed reference onto the config key the pipeline reads.

    Reads exactly the schema fixed in Task 4 and nothing else. Genome keys are derived
    from the registry rather than written by hand, so the writer and the readers cannot
    drift: neither of them owns the name.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Directory the references were installed into.

    Returns:
        dict[str, Path]: Absolute paths, keyed by `reference_data` key.

    Raises:
        KeyError: If an entry is missing a required schema field. Failing here is
            correct - a silently skipped entry is a silently missing reference.
    """
    prefix = REFERENCE_KINDS["bwa"]["prefix"]
    written: dict[str, Path] = {}

    for section in ("ucsc_references", "ncbi_references", "ensembl_references"):
        for physical_id, entry in install_config.get(section, {}).items():
            written[f"{prefix}_{physical_id}"] = (output_dir / entry["installed_path"]).resolve()

    for entry in install_config.get("common_references", []):
        written[entry["config_key"]] = (output_dir / entry["installed_path"]).resolve()

    for spec in install_config.get("derivations", []):
        if spec["kind"] == "shark":
            (key,) = reference_keys("shark", spec["assembly"])
        else:
            key = spec["config_key"]
        written[key] = (output_dir / spec["output"]).resolve()

    # Only name files that are actually there. A partial install (`--references hg38`)
    # must not write a key pointing at an hg19 FASTA nobody downloaded - that is the
    # same class of defect as #163 itself, just from the other direction.
    return {key: path for key, path in written.items() if path.exists()}
```

Add a schema-validation test: every genome entry has `kind` and `installed_path`, every
derivation has `kind` plus either `assembly` or `config_key`, every common reference has
both its fields — and a fixture with a field removed raises rather than silently dropping
the reference.

Replace `main`'s `:781-820` block with `update_config(config_path, canonical_reference_keys(install_config, output_dir))`.

Make `update_config` atomic: write to `config_path.with_suffix(".json.tmp")`, then
`os.replace`.

- [ ] **Step 5: Ship all six physical keys in `vntyper/config.json`**

```json
  "reference_data": {
    "muc1_reference_vntr": "reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
    "advntr_reference_vntr_hg19": "reference/vntr_db_advntr/hg19_muc1.db",
    "advntr_reference_vntr_hg38": "reference/vntr_db_advntr/hg38_muc1.db",
    "muc1_region_fasta_hg19": "reference/muc1_region_hg19.fa",
    "muc1_region_fasta_hg38": "reference/muc1_region_hg38.fa",
    "bwa_reference_hg19": "reference/alignment/chr1.hg19.fa",
    "bwa_reference_hg38": "reference/alignment/chr1.hg38.fa",
    "bwa_reference_GRCh37": "reference/alignment/chr1.GRCh37.fna",
    "bwa_reference_GRCh38": "reference/alignment/chr1.GRCh38.fna",
    "bwa_reference_hg19_ensembl": "reference/alignment/chr1.hg19_ensembl.fa",
    "bwa_reference_hg38_ensembl": "reference/alignment/chr1.hg38_ensembl.fa",
    "cram_reference_hg19": null,
    "cram_reference_hg38": null,
    "code_adVNTR_RUs": "reference/code-adVNTR_RUs.fa",
    "muc1_motifs_rev_com": "reference/MUC1_motifs_Rev_com.fa"
  },
```

- [ ] **Step 6: Run the tests**

```bash
pytest tests/unit/test_install_references_config_writing.py tests/unit/test_reference_key_agreement.py \
       tests/unit/test_install_references.py -v
```

Expected: all PASS, including all 24 agreement cases with `is_fallback` False.

**Ordering note.** Step 7 runs a *default* install, which is the bundle path — implemented
in Task 11. Either run Task 11 before this step, or run Step 7 with `--from-source`. The
plan assumes the former; if you are executing strictly in order, use `--from-source` here
and repeat the check after Task 11.

**Migrate the source readers in the same commit.** `process_ucsc_references` reads `url`
and `target_path` (`install_references.py:389-400`), but the rewritten schema uses
`source_url` and `installed_path`. Update those readers, or the source path silently reads
nothing. A schema test that loads the shipped config and asserts every reader's required
fields are present catches this.

- [ ] **Step 7: Reproduce the original #163 bug report and confirm it is fixed**

```bash
export PATH=~/miniforge3/envs/vntyper/bin:$PATH
cp vntyper/config.json /tmp/m5_cfg.json
vntyper --config-path=/tmp/m5_cfg.json install-references --output-dir=/tmp/m5_refs \
  --from-source --references GRCh38
python -c "
import json; d=json.load(open('/tmp/m5_cfg.json'))['reference_data']
import pathlib
for k,v in sorted(d.items()):
    print(('OK ' if v and pathlib.Path(v).exists() else '   '), k, '=', v)"
```

Expected: `bwa_reference_GRCh38` names an existing extracted `.fna`, both
`advntr_reference_vntr_*` name existing `.db` files, both `muc1_region_fasta_*` exist. No
key names a `.gz` or a `.zip`.

- [ ] **Step 8: Commit**

```bash
git add vntyper/scripts/install_references.py vntyper/scripts/install_references_config.json \
        vntyper/config.json tests/unit/test_install_references_config_writing.py \
        tests/unit/test_reference_key_agreement.py tests/unit/conftest.py
git commit -m "fix(install-references): write the keys the pipeline actually reads

--config-path wrote ucsc_*/ncbi_*/own_repo_* keys that nothing in vntyper/ reads,
while the pipeline kept its shipped relative paths - so an install into a custom
--output-dir left the run pointing at files that were never created. Two of the old
keys also named unusable files: the .gz rather than the extracted FASTA, and the
.zip rather than the two .db files.

Keys are now derived from the registry, so the writer and the readers cannot drift.
config.json ships all six physical BWA keys, so a complete installation never takes
the source-degrading UCSC-family fallback. The MUC1/adVNTR common asset installs
unconditionally, fixing a pre-existing gap where --references hg19 installed no
adVNTR database at all.

Fixes #163

Co-Authored-By: ElenaPianfetti <REPLACE@users.noreply.github.com>
Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 11: Bundle fetch in `install-references`

**Files:**
- Modify: `vntyper/scripts/install_references.py` (`main`)
- Test: `tests/unit/test_reference_bundle_install.py`

**Interfaces:**
- Consumes: `reference_bundle` (Task 3), the bundle block (Task 10).
- Produces: `install_from_bundle(install_config: dict, output_dir: Path, references: list[str]) -> None`

- [ ] **Step 1: Write the failing test**

Cover, with a fake release served from `tmp_path` and `download_file` monkeypatched:
happy path (assets fetched, verified, extracted, activated); `asset_sha256` mismatch →
`ValueError` naming the asset and both digests; a `BUILD_INFO.json` whose `bwa_version`
differs from local → warning logged and `bwa index` re-run; a 404 → error naming the
repository, tag and asset; `--references hg19` → the common asset is still installed;
interruption mid-extract → previous tree intact.

- [ ] **Step 2: Run it and watch it fail**

```bash
pytest tests/unit/test_reference_bundle_install.py -v
```

- [ ] **Step 3: Implement `install_from_bundle`**

Fetch `https://github.com/{repository}/releases/download/{release_tag}/{asset}` for each
selected physical id plus the common asset, `verify_sha256` against the committed
`asset_sha256`, `safe_extract` into a `staged_install` staging directory, then activate.
Route `main` to it unless `from_source` is set.

**`release-manifest.json` and `BUILD_INFO.json` must be inside every tarball, not fetched
separately.** A separately downloaded metadata file has no committed digest to check it
against, so trusting it would reintroduce exactly the trust hole the committed
`asset_sha256` closes: an attacker who can replace a release asset can replace a loose
JSON file too. Task 5's `bundle_release.py` therefore writes both into each archive, and
they are covered by that archive's `asset_sha256`. Only after the archive verifies does
the installer read them — comparing `BUILD_INFO.json`'s `bwa_version` with local
`bwa 2>&1 | head -3` and re-indexing locally with a warning on mismatch, and checking every
extracted file against the manifest.

Update Task 5's asset list accordingly: the release still publishes a top-level
`SHA256SUMS` and `verification-report.json` for the human reviewing the draft, but the
installer never depends on them.

- [ ] **Step 4: Run the tests, then a real install**

```bash
pytest tests/unit/test_reference_bundle_install.py -v
rm -rf /tmp/m5_bundle && vntyper install-references -d /tmp/m5_bundle --references hg38
ls /tmp/m5_bundle/alignment/ /tmp/m5_bundle/vntr_db_advntr/ /tmp/m5_bundle/*.fa
```

Expected: `chr1.hg38.fa` plus its five BWA index files, both `.db` files, both
`muc1_region_*.fa`, all three MUC1 motif FASTAs.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/install_references.py tests/unit/test_reference_bundle_install.py
git commit -m "feat(install-references): install from the pinned reference bundle

References now come from a versioned, checksummed release in berntpopp/vntyper-data
rather than from six third-party hosts at build time. Each asset's SHA-256 is
committed here rather than taken from the release's own SHA256SUMS, which is
co-hosted with the assets it describes and so cannot be the root of trust - and
because this file is a base-image content-hash input, changing any reference byte
now necessarily changes the base tag.

Refs #217

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 12: `bcrypt` and `shlex.quote` — the parked base-rebuild batch

**Files:**
- Modify: `docker/requirements-web.txt:16`, `pyproject.toml` (`web` and `dev` extras)
- Modify: `tests/unit/test_version_consistency.py` (the `UNDECLARED_IMPORT_ALLOWANCES` entry and its rationale comment)
- Modify: `vntyper/scripts/install_references.py:239-291`
- Modify: `tests/unit/test_shell_quoting.py` (docstring + three cases)

- [ ] **Step 1: Read the version the published base actually resolves**

```bash
docker run --rm --entrypoint bash ghcr.io/hassansaei/vntyper-base:latest \
  -c 'conda run -n vntyper pip show bcrypt | head -2'
```

Do **not** use the local version (5.0.0 here) — the pin must not change runtime behaviour.

- [ ] **Step 2: Edit requirements, extras and the allowance in ONE commit**

`docker/requirements-web.txt`: replace `passlib[bcrypt]==1.7.4` with `bcrypt==<that version>`.
Mirror in `pyproject.toml`'s `web` and `dev` extras.
Delete the `"bcrypt": "passlib"` entry from `UNDECLARED_IMPORT_ALLOWANCES` **and** the
bcrypt rationale comment immediately above it. **Locate both by symbol, not by line
number** — the allowance is near `:337-350`, and `:304-313` is a different test's body
(the required-image-binaries assertion), which must not be touched:

```bash
grep -n 'UNDECLARED_IMPORT_ALLOWANCES' -A 20 tests/unit/test_version_consistency.py | grep -n bcrypt
```

Add a positive assertion in the same commit: `bcrypt` now appears in
`web_requirements()`, and `passlib` does not.

`test_undeclared_import_allowances_are_all_still_needed` fires on "module now declared
directly", so splitting these across commits leaves the tree red at the intermediate commit.

- [ ] **Step 3: Restore the `shlex.quote` fix**

In `install_references.execute_aligner_index` (`:239-291`), quote `ref_path`, `index_dir`
and `index_base` with `command_builders.quote_path` before `.format(**params)`.

- [ ] **Step 4: Invert the three characterisation cases**

`tests/unit/test_shell_quoting.py` — rewrite the module docstring paragraph that says the
fifth site is "still unquoted, and the tests below characterise it that way on purpose",
and invert the three `execute_aligner_index` cases to assert `shlex.split()` yields one
operand per path, matching the four already-quoted sites.

- [ ] **Step 5: Run the affected tiers**

```bash
pytest tests/unit/test_version_consistency.py tests/unit/test_shell_quoting.py -v
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add docker/requirements-web.txt pyproject.toml tests/unit/test_version_consistency.py \
        vntyper/scripts/install_references.py tests/unit/test_shell_quoting.py
git commit -m "fix(docker): declare bcrypt directly and quote the aligner index command

docker/app/utils.py imports bcrypt at :7 and uses it at :59 and :83, but
requirements-web.txt declared only passlib[bcrypt]==1.7.4 - and nothing under
docker/app/ imports passlib, so bcrypt survived solely on an extra of an
otherwise-unused distribution. That exact area already produced a silent total
failure once: passlib 1.7.4's backend probe versus bcrypt >= 4.1 meant no cohort
passphrase ever worked (#179).

execute_aligner_index interpolated paths into a shell=True command unquoted. Both
were parked on a base-image rebuild, which this milestone provides.

Fixes #193

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 13: Delete the tracked reference data AND rewrite the Docker stage

> **These two were separate tasks and must not be.** `Dockerfile.base:85-99` copies only
> `install_references.py` into `/opt/ir` and runs it directly, so the moment Tasks 3-4 add
> imports of `reference_bundle`, `reference_registry`, `command_builders` and `utils`, that
> stage cannot import them and the base build fails at its most expensive step. Deleting
> the tracked data in one commit and fixing the Dockerfile in the next leaves an
> intermediate commit whose base build is broken twice over. One task, one commit.

**Files:**
- Delete: `reference/*.fa`, `reference/*.fa.fai`, `reference/vntr_db_advntr.zip`, `reference/md5_checksums.txt`, `reference/generate_vntr_reference.py`, `reference/filter_config.json`
- Modify: `.gitignore`, `.dockerignore:89-100`, `pyproject.toml:166-206`, `SPEC.md:23-29`, `AGENTS.md` trap 10

- [ ] **Step 1: Prove the published bundle carries every file about to be deleted**

```bash
for f in $(git ls-files reference/ | grep -E '\.(fa|fai|zip)$'); do
  n=$(basename "$f")
  test -f "/tmp/m5_bundle/$n" -o -f "/tmp/m5_bundle/alignment/$n" \
    && echo "OK   $n" || echo "MISS $n"
done
md5sum reference/muc1_region_hg19.fa /tmp/m5_bundle/muc1_region_hg19.fa
```

Expected: every file `OK`, and identical md5 for the SHARK reference. **Do not proceed on
any `MISS`** — R4 in the spec is about exactly this ordering.

- [ ] **Step 1b: Make the installer importable in the `refs` stage**

Tasks 3-4 kept `install_references.py` and `reference_bundle.py` free of `vntyper` imports
except each other, so only **two** Python files need copying — and only `reference_bundle.py`
is new to the content-hash set. Copying `command_builders` and `utils` instead would have
added two modules with 11 and 8 commits per year respectively to that set, making every
edit to either cost a 25-120 minute base rebuild.

```dockerfile
COPY vntyper/__init__.py                            /opt/ir/vntyper/__init__.py
COPY vntyper/scripts/__init__.py                    /opt/ir/vntyper/scripts/__init__.py
COPY vntyper/scripts/install_references.py          /opt/ir/vntyper/scripts/
COPY vntyper/scripts/install_references_config.json /opt/ir/vntyper/scripts/
COPY vntyper/scripts/reference_bundle.py            /opt/ir/vntyper/scripts/
ENV PYTHONPATH=/opt/ir
RUN conda run -n vntyper python -m vntyper.scripts.install_references --output-dir /opt/refs ...
```

`reference_registry` is deliberately absent: `canonical_reference_keys` imports it
function-locally, and Docker passes no `--config-path` (`Dockerfile:50-56`), so that code
never runs in the image.

Guard the property so it cannot silently rot:

```python
def test_the_docker_installer_modules_import_nothing_they_cannot_reach():
    """A module-scope vntyper import the refs stage cannot resolve breaks the base build.

    It breaks it at the most expensive step, and no unit test would have caught it -
    only a 25-120 minute Docker run. Hence this AST check.
    """
    import ast
    import pathlib
    import re

    copied = set(
        re.findall(r"vntyper/scripts/(\w+)\.py", pathlib.Path("docker/Dockerfile.base").read_text())
    )
    for name in ("install_references", "reference_bundle"):
        tree = ast.parse(pathlib.Path(f"vntyper/scripts/{name}.py").read_text())
        imported = {
            node.module.split(".")[-1]
            for node in tree.body                      # module scope only - nested imports are fine
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("vntyper.")
        }
        assert imported <= copied, f"{name}.py imports {sorted(imported - copied)}, which the refs stage lacks"
```

- [ ] **Step 2: Delete, and update the ignore files**

```bash
git rm reference/*.fa reference/*.fa.fai reference/vntr_db_advntr.zip reference/md5_checksums.txt
git rm reference/generate_vntr_reference.py reference/filter_config.json
```

`.gitignore`: replace the `/reference/chr1*`, `/reference/*/` block with a comment saying
reference data now comes from `berntpopp/vntyper-data`, and ignore the whole populated tree
except the retained `README.md`, `pseudonymize.py` and `pseudonymize_config.json`.

`.dockerignore:89-100`: its comments declare tracked reference files as image inputs.
Rewrite them to say the image fetches the bundle.

`pyproject.toml:166-206`: Ruff excludes `reference/` *because* it was a base-hash input.
That reason has expired — remove the exclusion, or restate it truthfully.

`SPEC.md:23-29`: its proof command runs against a deleted file. Repoint it at a bundle
install or drop it.

`AGENTS.md` trap 10: `reference/**` is no longer a base input. Say so, and say the bundle
pin in `install_references_config.json` is what triggers a rebuild instead.

- [ ] **Step 3: Verify nothing still reads a deleted path**

```bash
grep -rn "muc1_region_hg19.fa\|vntr_db_advntr.zip\|generate_vntr_reference\|md5_checksums" \
  --include='*.py' --include='*.json' --include='*.yml' --include='*.md' \
  vntyper/ tests/ docker/ scripts/ .github/ Makefile | grep -v docs/superpowers/
```

Expected: only runtime path *strings* in `config.json` and Docker guards, never a
repository path.

- [ ] **Step 4: Run the full suite**

```bash
make format && make test-unit
```

Expected: green. Any failure here names a consumer the migration missed.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "chore(references): move reference data to berntpopp/vntyper-data

reference/ keeps README.md and pseudonymize* (pseudonymize.py and pseudonymize_config.json). The MUC1 region FASTAs are
derivations verified against digests in the release spec; the three non-derivable
seeds live in the data repo's seeds/ alongside generate_vntr_reference.py and its
filter config. reference/** leaves the base-image content-hash set - the bundle pin
in install_references_config.json, already in that set, is what triggers a rebuild
now.

This also removes the raw.githubusercontent.com/.../main/... coupling that PR #164
tripped over: reference files come from a release, not a mutable branch.

Refs #217

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 14: Content hash, workflows and the image tier

> Runs in the **same commit** as Task 13 — see the note there. Split into its own task only
> because it has its own verification loop.

**Files:**
- Modify: `docker/Dockerfile.base:72-105`
- Modify: `.github/workflows/docker-base.yml:66`, `docker-build.yml:96`, `docker-build.yml:168`
- Modify: `Makefile:511-513`, `tests/unit/test_workflow_consistency.py:82-85`
- Modify: `tests/docker/test_image_structure.py:56-57`, `tests/docker/image_probe.py`

- [ ] **Step 1: Update the workflow-consistency test first**

Remove `"reference/**"` from the expected path list at
`tests/unit/test_workflow_consistency.py:82-85` and run it — it must now FAIL, because the
workflows still list it. That failure is the checklist.

```bash
pytest tests/unit/test_workflow_consistency.py -v
```

- [ ] **Step 2: Remove `reference/**` and ADD the newly-copied installer modules**

`reference/**` leaves the set, and `reference_bundle.py` **enters** it — otherwise the base
could be built against a stale copy of a module the installer imports, which is precisely
the failure the content hash exists to prevent. It is one file, not four, because Tasks 3-4
kept the installer self-contained.

Edit `docker-base.yml:66`, `docker-build.yml:96`, `docker-build.yml:168` and
`Makefile:511-513` (`BASE_INPUTS`) identically: remove `reference/**`, add exactly
`vntyper/scripts/reference_bundle.py`. That is the only module the refs stage newly copies.
Update the expected list in `tests/unit/test_workflow_consistency.py:82-85` to match.
Re-run the test — it must now PASS.

- [ ] **Step 3: Rewrite the `refs` stage**

In `docker/Dockerfile.base`, delete `COPY reference/ /opt/refs/` and its comment block
(`:78-83`), and change the install step to fetch the bundle. Extend the guards:

```dockerfile
    test -d /opt/refs/vntr_db_advntr || { echo "adVNTR database missing"; exit 1; } && \
    test -f /opt/refs/muc1_region_hg19.fa || { echo "SHARK hg19 region missing"; exit 1; } && \
    test -f /opt/refs/muc1_region_hg38.fa || { echo "SHARK hg38 region missing"; exit 1; } && \
```

- [ ] **Step 4: Extend the image tests**

`tests/docker/test_image_structure.py`: add `muc1_region_hg38.fa` to the probed set and add
minimum sizes for the four new BWA references. `tests/docker/image_probe.py:81-105` already
discovers path-valued main and SHARK config entries; extend its BWA-sidecar check beyond
hg19/hg38 to all six physical keys. **The relative-path assertion the spec assumed does not exist.** `test_workdir_is_prefix`
(`:243-245`) has that docstring but asserts only
`image_metadata["Config"]["WorkingDir"] == "/opt/vntyper"`, and
`test_every_declared_reference_exists` (`:345-355`) checks resolved existence, not
relativity. Add the missing assertion, since the whole "layout did not move" argument
rests on it:

```python
@pytest.mark.parametrize("key", sorted(MINIMUM_SIZES))
def test_declared_reference_paths_stay_relative(probe: dict[str, Any], key: str) -> None:
    """They resolve against WORKDIR; an absolute path would bind the image to a host."""
    assert not Path(probe["references"][key]["declared"]).is_absolute()
```

- [ ] **Step 5: Build the base locally and run the image tier**

```bash
make docker-build-base
make docker-build DOCKER_BASE_IMAGE=vntyper-base:local
pytest tests/docker/ -v
```

Expected: the base build is dramatically faster than the 25-120 min it used to take, and
the image tier passes.

- [ ] **Step 6: Commit**

```bash
git add docker/Dockerfile.base .github/workflows/ Makefile \
        tests/unit/test_workflow_consistency.py tests/docker/
git commit -m "build(docker): fetch the reference bundle instead of building it in-image

The base build downloaded six chr1 genomes from UCSC/NCBI/Ensembl and BWA-indexed
them - the single most expensive step, and the main reason the image is split in
two. It now fetches a pinned, checksummed bundle. reference/** leaves the
content-hash set at all three call sites, the Makefile mirror and the consistency
test; the bundle pin in install_references_config.json, already hashed, is what
triggers a rebuild now. The in-image layout is unchanged.

Refs #217

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 15: CLI/server assembly parity

**Files:**
- Modify: `docker/app/main.py:58-63`
- Test: `tests/unit/web/test_assembly_parity.py` (create)

- [ ] **Step 1: Write the failing test**

```python
"""The CLI and the server must accept the same assemblies.

`vntyper online --reference-assembly hg38_ensembl` subsets and uploads the BAM, then
gets a 422 because the server enum has only four values - reported to the user as a
generic submission failure, after the upload.
"""

import pytest

from vntyper.scripts.reference_registry import list_assemblies

pytestmark = pytest.mark.unit


def test_the_server_enum_matches_the_registry():
    from docker.app.main import ReferenceAssembly

    assert {member.value for member in ReferenceAssembly} == set(list_assemblies())
```

- [ ] **Step 2: Run it and watch it fail**

Expected: the sets differ by the four source-qualified labels.

- [ ] **Step 3: Add the four missing members**

```python
class ReferenceAssembly(str, Enum):
    HG19 = "hg19"
    HG38 = "hg38"
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
    HG19_NCBI = "hg19_ncbi"
    HG38_NCBI = "hg38_ncbi"
    HG19_ENSEMBL = "hg19_ensembl"
    HG38_ENSEMBL = "hg38_ensembl"
```

- [ ] **Step 4: Run the web tier**

```bash
pytest tests/unit/web/ -v
```

- [ ] **Step 5: Commit**

```bash
git add docker/app/main.py tests/unit/web/test_assembly_parity.py
git commit -m "fix(web): accept every assembly the CLI accepts

vntyper online accepted all eight registry labels and submitted them verbatim, but
the server enum had four - so the four source-qualified labels uploaded the BAM and
then got a 422, surfaced as a generic submission failure. A parity test against
list_assemblies() stops the two drifting again.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 16: Documentation

**Files:** `docs/cli/install-references.md`, `docs/getting-started/reference-setup.md`, `docs/user-guide/reference-assemblies.md`, `docs/user-guide/configuration.md`, `docs/pipeline/optional-modules.md`, `docs/development/golden-cohort-gate.md`, `README.md`, `reference/README.md`

- [ ] **Step 1: Rewrite the reference-setup and install-references pages**

Document: the bundle is the default; `--from-source` exists and what it costs; the trust
model (digests committed here, `SHA256SUMS` corroborates); `berntpopp/vntyper-data` and its
release lifecycle; that `--references hg19` still installs the common MUC1/adVNTR assets.

- [ ] **Step 2: Rewrite the reference-assemblies page**

Add the physical identity table from the top of this plan. State that `GRCh38` and
`hg38_ncbi` name one file. State that a config missing a physical key falls back to the
UCSC family with a warning, and that a complete installation never does.

- [ ] **Step 3: Correct the SHARK documentation**

`docs/pipeline/optional-modules.md` currently states that assembly does not select a region
(added by `39fe0eb` for #187). That is now wrong. Replace it with the selection behaviour
and the measurement: 40.6 % of the hg38 region's canonical 17-mers are absent from the hg19
region; the hg38 reference retains 3.2–34.7 % more reads across the seven cohort samples.

- [ ] **Step 4: Add the golden-gate prerequisite**

`docs/development/golden-cohort-gate.md:897-905` assumes both trees carry tracked reference
data. Add a first step: install and verify the published bundle into the shared reference
tree.

- [ ] **Step 5: Trim `reference/README.md`**

It documents `generate_vntr_reference.py` and the adVNTR database curation, both of which
now live in the data repo. Keep the `pseudonymize.py` sections; point the rest at
`berntpopp/vntyper-data`.

- [ ] **Step 6: Build the docs and commit**

```bash
mkdocs build --strict
git add docs/ README.md reference/README.md
git commit -m "docs: reference bundles, assembly identity, and SHARK region selection

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KH1Z9wjJ2dZspU68VmXpTR"
```

---

### Task 17: Golden-cohort gate and release

**Files:** `scripts/golden_cohort/matrix.py:108-109`; then the release cycle.

- [ ] **Step 1: Actually generate cases for the two uncovered labels**

`hg19_ncbi` and `hg38_ncbi` are the labels the physical-identity change affects most, and
the gate does not cover them. **Adding them to `ASSEMBLIES` (`matrix.py:109`) achieves
nothing**: `rg -n '\bASSEMBLIES\b' scripts/golden_cohort/` shows it is referenced only by
its own definition and a docstring. Cases are derived from fixture filenames and
directories at `:230-283`, and the counts are pinned at `:172-190`.

So: generate explicit alias cases that reuse the GRCh37/GRCh38 NCBI BAMs under the
`hg19_ncbi`/`hg38_ncbi` labels, then update the pinned case counts, the matrix tests, and
`docs/development/golden-cohort-gate.md`. Verify with:

```bash
python -c "from scripts.golden_cohort import matrix; print(len(matrix.build_matrix()))"
```

Expected: the documented total rises by the number of alias cases added.

- [ ] **Step 2: Install the bundle into both gate trees, then run the gate**

```bash
export PATH=~/miniforge3/envs/vntyper/bin:$PATH
vntyper install-references -d reference --references hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl
python scripts/golden_cohort/run_gate.py   # ~14 min on 32 cores
python -m json.tool results/golden_cohort/result.json | head -40
```

Genotype-affecting changes to expect and confirm as intended: the adVNTR database for
`hg38_ncbi`/`hg38_ensembl` (previously the GRCh37 database), and SHARK's region for every
GRCh38 run. Generate the CRAM fixtures rather than passing `--allow-matrix-drift`.

- [ ] **Step 3: Full verification**

```bash
make check-all        # format-check, lint, type-check-all, test-unit
make ci-local         # workflow linting - check-all does NOT run lint-actions
make ci-local-docker  # the image jobs - check-all does NOT run these either
pytest tests/docker/ -v
```

All four are required: this milestone changes workflows *and* the Dockerfile, and
`make check-all` covers neither.

- [ ] **Step 4: Open the PR with the gate result**

```bash
gh pr create --repo hassansaei/VNtyper --base main \
  --head fix/milestone-5-references-configuration \
  --title "fix: milestone 5 — references and configuration (#163 #152 #193 #217)" \
  --body-file docs/superpowers/specs/2026-08-11-milestone-5-references-and-configuration-design.md
```

(PR-1 was opened and merged after Task 5, with the title
`feat: reference bundle builder and --from-source (#217)` and a body pointing at §4.6.)

Post the golden-gate result as a PR comment. Note that `Docker Build` takes 39–69 min and
builds from a GitHub clone rather than the PR branch, so `ci-tests.yml` → `ci-success` is
the gate to watch.

- [ ] **Step 5: Merge, then hand the release to an operator**

Bump `vntyper/version.py`, `CITATION.cff` and `docs/about/changelog.md` to the same patch
version — all three must agree (AGENTS.md trap 12). Merge. Delete the branch.

**Stop there. Do not tag.** `AGENTS.md` "Never" list: *"Never create, move, or push a
release tag as an agent."* And `.github/workflows/publish-pypi.yml` has **no `push.tags`
trigger** — production starts only from an authenticated `repository_dispatch` of type
`vntyper_release` whose `client_payload.tag` names an **already existing** reviewed tag. A
tag push publishes nothing; an earlier belief that it published to PyPI immediately is out
of date, superseded by the milestone-6 release controller.

The runbook to hand over:

1. Confirm the ten exact-full-SHA checks are green on the merge commit.
2. *Operator* creates the annotated tag on that commit.
3. *Operator* sends the authenticated `vntyper_release` dispatch naming that tag.

Inspect without writing at any point:

```bash
gh workflow run publish-pypi.yml -f tag=vX.Y.Z
```

- [ ] **Step 6: Bump the data repo pin and close the issues**

```bash
MERGE_SHA=$(git rev-parse origin/main)
# in /tmp/vntyper-data: update release-data.yml's `uses:` and `source_commit` to $MERGE_SHA
```

Close #163, #152, #193 and #217. Close PR #164 with a comment thanking
`@ElenaPianfetti`, pointing at the superseding PR, and noting the commits that carry
`Co-Authored-By`.

---

## Self-Review

**Spec coverage.** Every §9 deliverable maps to a task: registry resolver → Task 1;
five readers → Tasks 2, 7, 8, 9; writer and config → Task 10; bundle install → Tasks 3, 11;
derivations and `--from-source` → Task 4; build workflow and data repo → Tasks 5, 6;
deletion and migration checklist → Task 13; Docker and content hash → Task 14; #193 →
Task 12; server parity → Task 15; docs → Task 16; gate and release → Task 17.

**Known gaps, deliberately left to implementation.** Three items are specified by
behaviour and test rather than by finished code, because their content depends on values
that only exist after Pass 1 runs: `scripts/bundle_release.py` and
`vntyper/scripts/verify_seeds.py` (Task 5 Step 2), and the digests in `refs-v1.json`
(Task 6 Step 3) and `install_references_config.json` (Task 10 Step 3). Each states its
inputs, outputs and tests. The SHARK case-insensitivity test body (Task 9 Step 1) is left
one step short — build the FASTQ as plain text and assert equal retained counts.

**Type consistency.** `reference_keys(kind, assembly) -> tuple[str, ...]` is used
identically in Tasks 1, 2, 9, 10. `resolve_from_mapping(kind, assembly, mapping) ->
ResolvedReference | None` with fields `key`/`value`/`is_fallback` is used identically in
Tasks 2, 7, 8, 9, 10. `select_bwa_reference(config, assembly) -> str` is defined in Task 7
and consumed by `cli_logging_safety` in the same task. `canonical_reference_keys(install_config,
output_dir) -> dict[str, Path]` is defined in Task 10 and consumed by
`test_reference_key_agreement` in the same task.
