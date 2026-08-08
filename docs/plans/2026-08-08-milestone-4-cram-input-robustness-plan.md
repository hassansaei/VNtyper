# Milestone 4 — CRAM and Input Robustness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps
> use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give VNtyper an input contract, so a CRAM or BAM run either succeeds or fails at
submission naming exactly what is missing, and no input is ever silently discarded.

**Architecture:** One preflight, run after `validate_bam_file` and before any stage,
resolves and *proves* everything a run needs (index, reference, read layout) and returns a
frozen `AlignmentPlan`. The decisions live in three pure modules that need no filesystem;
only the samtools calls stay behind I/O. `process_bam_to_fastq` consumes the plan instead
of rediscovering its contents mid-run.

**Tech Stack:** Python 3.10+ (floor — CI matrix 3.10–3.13), samtools 1.20 / htslib 1.23,
pytest, ruff, mypy, Docker + Celery for the web layer.

**Spec:** [`2026-08-08-milestone-4-cram-input-robustness-spec.md`](2026-08-08-milestone-4-cram-input-robustness-spec.md).
Read §3 (measured evidence) before writing any CRAM code — several intuitions about
htslib are wrong and the measurements say which.

## Global Constraints

Copied verbatim from AGENTS.md. Every task's requirements implicitly include these.

- **Run pytest from the repo root.** `tests/parametrization.py` opens
  `tests/test_data_config.json` by relative path at collection time.
- **Activate the `vntyper` conda env and put its bin FIRST on `PATH`.** The base env and
  `~/.local/bin` shadow it, and `make check-all` will otherwise report on the wrong tools.
- **Every new test file needs `pytestmark = pytest.mark.unit` after the imports.** CI runs
  `pytest -m unit`; an unmarked file silently never runs.
  `tests/unit/test_marker_hygiene.py` fails the build naming the offending file.
- **Ruff only. Line length 120, double quotes, `target-version = "py310"`.** Do not
  reformat to 88.
- **Code must run on Python 3.10.** Your local interpreter is newer. PEP 604 (`X | None`)
  and PEP 585 (`list[str]`) are available; `match`, `Self` and PEP 695 generics are not.
- **Google-style docstrings** (`Args:` / `Returns:` / `Raises:`) on public functions.
- **`logger = logging.getLogger(__name__)` after the imports.** Never `logging.info(...)`,
  never `basicConfig`, never a per-module handler. f-strings in log calls are the style.
- **No custom exception classes.** `logger.error(msg)` then
  `raise ValueError(msg)` / `RuntimeError(msg)` with the same message.
- **Never lower the coverage floor** (`fail_under` in `pyproject.toml`). Add the test.
  `make patch-coverage` must be ≥ 80% on changed lines.
- **`--config-path` replaces the whole config, it does not merge** (trap 2). Every read of
  a new config key uses `.get` with the shipped default; a `KeyError` must not abort a run
  over a threshold.
- **`run_command` uses `shell=True` deliberately** (trap 9). Quoting can only happen at
  construction time, in `command_builders.py`, via `quote_path`.
- **Tool invocations are not quoted; paths, regions and sample names are.**
  `config["tools"]` holds command *prefixes* (`advntr` is `"mamba run -n envadvntr advntr"`).
- **Keep files under ~650 LOC.** `fastq_bam_processing.py` is already 612. New logic goes
  in new modules; extract the region you touch rather than growing the file.
- **Touch a file, add tests for it** — including for functions you only modify.
- **If you touch `.github/workflows/`, `make check-all` is not enough — run `make ci-local`.**
- **Never claim tests pass without showing the command output.**

**Standard verification block** — referenced by every task as "run the standard block":

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate vntyper
export PATH="$CONDA_PREFIX/bin:$PATH"
cd "$(git rev-parse --show-toplevel)"
make format && make test-unit
```

---

## File Structure

| File | Status | Responsibility |
| --- | --- | --- |
| `vntyper/scripts/alignment_contract.py` | create | The `AlignmentPlan` dataclass; what each format requires; the exact text of every failure message. Pure. |
| `vntyper/scripts/reference_resolution.py` | create | Ordered reference candidates; which candidate wins given probe outcomes; which header contigs a FASTA fails to cover. Pure. |
| `vntyper/scripts/read_layout.py` | create | The `paired`/`single`/`mixed`/`empty` verdict; which FASTQs feed downstream and which are unconsumed. Pure. |
| `vntyper/scripts/alignment_preflight.py` | create | The only new module that shells out: resolve or build the index, probe-decode each candidate, count layout. Returns `AlignmentPlan`. |
| `vntyper/scripts/alignment_index.py` | modify | Add `resolve_cram_index`. `resolve_bam_index` is **unchanged** — the offset extractor parses BAI directly. |
| `vntyper/scripts/command_builders.py` | modify | Thread the slice and index (P1, P2); make slice indexing optional (P3); add the indexed unplaced-fetch builder (P4); accept a real `cram_ref_option`. |
| `vntyper/scripts/fastq_bam_processing.py` | modify | Consume `AlignmentPlan`. Net LOC must not increase — the index-resolution block moves out. |
| `vntyper/scripts/chromosome_utils.py` | modify | `detect_naming_convention` scores against matching contigs, threshold from config. |
| `vntyper/scripts/pipeline.py` | modify | Call the preflight; route the layout; stop binding FASTQs to `_`. |
| `vntyper/scripts/alignment_processing.py` | modify | `align_and_sort_fastq` accepts a single FASTQ. |
| `vntyper/scripts/kestrel_genotyping.py` | modify | `construct_kestrel_command` accepts a single FASTQ. |
| `vntyper/scripts/cli_parser.py` | modify | `--reference-fasta`. |
| `vntyper/scripts/cli_handlers.py` | modify | Thread `--reference-fasta` into `run_pipeline`. |
| `vntyper/config.json` | modify | The §6 config surface. |
| `docker/entrypoint.sh` | modify | `--no-capture-output` ×4. |
| `scripts/make_cram_fixtures.py` | modify | Declared fixtures by default, `--all`, plus a reference-dependent fixture. |

---

## Wave 0 — #213

Lands first and alone. Every later failure is illegible until it does.

### Task 1: `conda run --no-capture-output`

**Files:**
- Modify: `docker/entrypoint.sh:158,163,168,186`
- Test: `tests/docker/test_entrypoint_streaming.py` (create)

**Interfaces:**
- Consumes: nothing.
- Produces: nothing importable. The behaviour change is that the container streams logs.

- [ ] **Step 1: Write the failing test**

`tests/docker/test_entrypoint_streaming.py`:

```python
"""The image must stream logs while a command runs, not buffer them until it exits.

`conda run` without `--no-capture-output` buffers all child stdout and stderr until the
child exits, so a running pipeline produces no output at all in `docker logs` (#213).
The rest of the Docker suite invokes the container through `bash -c` and passes the flag
itself, so it exercises the correct form while production used the buffering one -- this
test runs the image **through its entrypoint**, which is the path users take.
"""

import subprocess

import pytest

pytestmark = [pytest.mark.docker, pytest.mark.slow]

ENTRYPOINT_CONDA_RUN_SITES = 4


def test_entrypoint_passes_no_capture_output(repo_root):
    """Every `conda run` in the entrypoint streams rather than buffers."""
    entrypoint = (repo_root / "docker" / "entrypoint.sh").read_text()
    conda_run_lines = [ln.strip() for ln in entrypoint.splitlines() if "conda run" in ln and "exec" in ln]
    assert len(conda_run_lines) == ENTRYPOINT_CONDA_RUN_SITES, (
        f"expected {ENTRYPOINT_CONDA_RUN_SITES} exec'd `conda run` sites, found {len(conda_run_lines)}: "
        f"{conda_run_lines}"
    )
    missing = [ln for ln in conda_run_lines if "--no-capture-output" not in ln]
    assert not missing, f"these `conda run` invocations buffer their child's output: {missing}"


def test_entrypoint_streams_before_exit(docker_image):
    """A line reaches the caller while the child is still running."""
    proc = subprocess.Popen(
        [
            "docker", "run", "--rm", "--entrypoint", "/opt/vntyper/docker/entrypoint.sh",
            docker_image, "vntyper", "--version",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    out, _ = proc.communicate(timeout=120)
    assert proc.returncode == 0, out
    assert out.strip(), "the entrypoint produced no output at all"
```

- [ ] **Step 2: Run it and watch it fail**

```bash
pytest tests/docker/test_entrypoint_streaming.py::test_entrypoint_passes_no_capture_output -v
```

Expected: FAIL — `these conda run invocations buffer their child's output: [...]` listing
all four lines.

- [ ] **Step 3: Add the flag at all four sites**

In `docker/entrypoint.sh`, at lines 158, 163, 168 and 186, change

```bash
exec conda run -n "${CONDA_ENV}" "$@"
```

to

```bash
exec conda run --no-capture-output -n "${CONDA_ENV}" "$@"
```

and the same for the `celery ... beat` invocation at `:168` and the `uvicorn` invocation
at `:186`. Add this comment above the first one:

```bash
    # --no-capture-output is load-bearing, not cosmetic: without it conda buffers all
    # child stdout and stderr until the child exits, so a running pipeline emits nothing
    # at all in `docker logs` and every "it hangs" report arrives with no stage log
    # attached (#213). For the long-running services below it also restores live
    # streaming, which matters for the same reason.
```

- [ ] **Step 4: Run it and watch it pass**

```bash
pytest tests/docker/test_entrypoint_streaming.py::test_entrypoint_passes_no_capture_output -v
```

Expected: PASS.

- [ ] **Step 5: Check the shell lints**

```bash
make lint-docker
```

Expected: no new findings.

- [ ] **Step 6: Commit**

```bash
git add docker/entrypoint.sh tests/docker/test_entrypoint_streaming.py
git commit -m "fix(docker): stream pipeline logs by passing conda run --no-capture-output

Closes #213"
```

---

## Wave 1 — the keystone (#225 and the preflight seam)

Everything in Wave 2 consumes `AlignmentPlan`. This wave must land before any of it
starts. It is four tasks; do them in order.

### Task 2: `alignment_contract.py` — the plan and the messages

**Files:**
- Create: `vntyper/scripts/alignment_contract.py`
- Test: `tests/unit/test_alignment_contract.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `AlignmentPlan` frozen dataclass with fields
    `input_path: str`, `file_format: str`, `index_path: str | None`,
    `reference_path: str | None`, `reference_source: str | None`,
    `uncovered_contigs: tuple[str, ...]`, `layout: str`.
  - `index_candidate_names(in_path: str | Path, file_format: str) -> tuple[str, str]`
  - `missing_index_message(in_path: str, file_format: str, tried: Sequence[str]) -> str`
  - `unresolvable_reference_message(in_path: str, contig: str, m5: str | None,
    attempts: Sequence[tuple[str, str, str]]) -> str`
    where each attempt is `(source, path, reason)`.
  - Constants `FORMAT_BAM = "bam"`, `FORMAT_CRAM = "cram"`.

- [ ] **Step 1: Write the failing test**

`tests/unit/test_alignment_contract.py`:

```python
"""The alignment contract: what a format requires, and how a missing piece is named."""

import pytest

from vntyper.scripts.alignment_contract import (
    FORMAT_BAM,
    FORMAT_CRAM,
    AlignmentPlan,
    index_candidate_names,
    missing_index_message,
    unresolvable_reference_message,
)

pytestmark = pytest.mark.unit


class TestIndexCandidateNames:
    def test_bam_carries_bai_under_both_spellings(self):
        assert index_candidate_names("/d/sample.bam", FORMAT_BAM) == ("/d/sample.bam.bai", "/d/sample.bai")

    def test_cram_carries_crai_under_both_spellings(self):
        assert index_candidate_names("/d/sample.cram", FORMAT_CRAM) == ("/d/sample.cram.crai", "/d/sample.crai")

    def test_an_unknown_format_is_rejected_rather_than_guessed(self):
        with pytest.raises(ValueError, match="unknown alignment format"):
            index_candidate_names("/d/sample.sam", "sam")


class TestMissingIndexMessage:
    def test_it_names_the_file_the_format_and_every_name_tried(self):
        message = missing_index_message("/d/s.cram", FORMAT_CRAM, ["/d/s.cram.crai", "/d/s.crai"])
        assert "/d/s.cram" in message
        assert "/d/s.cram.crai" in message and "/d/s.crai" in message
        assert "index" in message.lower()

    def test_it_says_what_to_do_about_it(self):
        assert "samtools index" in missing_index_message("/d/s.cram", FORMAT_CRAM, ["/d/s.cram.crai"])


class TestUnresolvableReferenceMessage:
    def test_it_names_the_contig_the_digest_and_every_candidate_with_its_reason(self):
        message = unresolvable_reference_message(
            "/d/s.cram",
            contig="chr1",
            m5="4047df6392200a575d41cc627247f231",
            attempts=[
                ("--reference-fasta", "", "not supplied"),
                ("config cram_reference_hg38", "/r/g.fa", "file does not exist"),
                ("config bwa_reference_hg38", "/r/chr1.hg38.fa", "probe exited 1: MD5 checksum reference mismatch"),
            ],
        )
        assert "chr1" in message
        assert "4047df6392200a575d41cc627247f231" in message
        assert "--reference-fasta" in message
        assert "MD5 checksum reference mismatch" in message
        assert "/r/chr1.hg38.fa" in message

    def test_a_cram_with_no_m5_still_produces_a_usable_message(self):
        message = unresolvable_reference_message("/d/s.cram", contig="1", m5=None, attempts=[("cli", "", "not supplied")])
        assert "1" in message
        assert "no M5" in message


class TestAlignmentPlan:
    def test_it_is_frozen_so_no_stage_can_mutate_the_contract_mid_run(self):
        plan = AlignmentPlan(
            input_path="/d/s.cram",
            file_format=FORMAT_CRAM,
            index_path="/d/s.cram.crai",
            reference_path="/r/g.fa",
            reference_source="--reference-fasta",
            uncovered_contigs=(),
            layout="paired",
        )
        with pytest.raises(AttributeError):
            plan.reference_path = "/r/other.fa"  # type: ignore[misc]

    def test_cram_ref_option_is_the_flag_fragment_the_builders_interpolate(self):
        plan = AlignmentPlan(
            input_path="/d/s.cram",
            file_format=FORMAT_CRAM,
            index_path=None,
            reference_path="/r/genome ref.fa",
            reference_source="cli",
            uncovered_contigs=(),
            layout="paired",
        )
        assert plan.cram_ref_option == "-T '/r/genome ref.fa'"

    def test_no_reference_yields_an_empty_fragment_so_no_ref_crams_are_unchanged(self):
        plan = AlignmentPlan(
            input_path="/d/s.cram",
            file_format=FORMAT_CRAM,
            index_path=None,
            reference_path=None,
            reference_source=None,
            uncovered_contigs=(),
            layout="paired",
        )
        assert plan.cram_ref_option == ""
```

- [ ] **Step 2: Run it and watch it fail**

```bash
pytest tests/unit/test_alignment_contract.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'vntyper.scripts.alignment_contract'`.

- [ ] **Step 3: Write the module**

`vntyper/scripts/alignment_contract.py`:

```python
"""
vntyper/scripts/alignment_contract.py

Module Purpose:
---------------
What a run needs from its input alignment, decided without touching the filesystem.

VNtyper used to discover an alignment's requirements -- its index, its reference -- inside
stage code at the moment each was needed, and delegate resolution of the second to
htslib's ambient environment. That is late (the run has already been accepted and
started), opaque (htslib reports what it could not fetch, not what the operator should
supply) and unbounded in time (a network ``REF_PATH`` blocks with no timeout, which is
why #178 is a hang rather than an error). This module holds the half of that decision
which is pure: the names an index can carry, and the exact wording of every failure.
:mod:`vntyper.scripts.alignment_preflight` holds the half that runs samtools.

:class:`AlignmentPlan` is the frozen result. Stages consume it instead of re-deriving its
contents, so there is no second place where a run can decide it needs something.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.command_builders import quote_path

logger = logging.getLogger(__name__)

#: The two alignment containers `pipeline.py` accepts. FASTQ has no header and no index
#: and never reaches this module.
FORMAT_BAM = "bam"
FORMAT_CRAM = "cram"

#: Index suffix per container. BAM's BAI is deliberately the only BAM entry -- see
#: `alignment_index.resolve_bam_index` for why a CSI must not be resolved here.
_INDEX_SUFFIX = {FORMAT_BAM: "bai", FORMAT_CRAM: "crai"}


@dataclass(frozen=True)
class AlignmentPlan:
    """Everything a run needs from its input, resolved and proven before any stage runs.

    Frozen on purpose: the contract is settled once, at submission, and a stage that
    could mutate it would reintroduce the second decision point this module exists to
    remove.

    Attributes:
        input_path (str): The alignment the run was given.
        file_format (str): ``"bam"`` or ``"cram"``.
        index_path (str | None): The index to use, existing or freshly built. None only
            for a format/mode that provably needs none.
        reference_path (str | None): The reference proven to decode this CRAM, or None
            for a BAM and for a reference-free (``no_ref=1``) CRAM.
        reference_source (str | None): Which candidate won, for the log and the report.
        uncovered_contigs (tuple[str, ...]): Header contigs the winning reference does
            not contain. Empty is the good case; non-empty is a warning, because
            ``samtools view -P`` chases mates onto other contigs and will abort the run
            if it reaches one of these (spec §3.7).
        layout (str): The read layout verdict from :mod:`vntyper.scripts.read_layout`.
    """

    input_path: str
    file_format: str
    index_path: str | None
    reference_path: str | None
    reference_source: str | None
    uncovered_contigs: tuple[str, ...]
    layout: str

    @property
    def cram_ref_option(self) -> str:
        """The ``-T <ref>`` fragment the command builders interpolate verbatim.

        Returns:
            str: ``-T <quoted path>`` when a reference was proven, and ``""`` otherwise
            -- which is what a reference-free CRAM needs and what a BAM ignores.

        Note:
            Quoted here rather than at the call site because the builders interpolate
            this fragment **unquoted** by contract (it is a flag fragment, not a path),
            and a reference path with a space in it would otherwise split into two
            arguments.
        """
        if not self.reference_path:
            return ""
        return f"-T {quote_path(self.reference_path)}"


def index_candidate_names(in_path: str | Path, file_format: str) -> tuple[str, str]:
    """Both names an index for this alignment can carry, in resolution order.

    Args:
        in_path (str | Path): The alignment.
        file_format (str): ``"bam"`` or ``"cram"``.

    Returns:
        tuple[str, str]: ``<file>.<suffix>`` then ``<stem>.<suffix>``. Both spellings
        occur in the wild and the upload endpoint accepts both (#210).

    Raises:
        ValueError: If the format is neither BAM nor CRAM. Guessing a suffix for an
            unknown container would build an index nothing can read.
    """
    suffix = _INDEX_SUFFIX.get(file_format)
    if suffix is None:
        message = f"unknown alignment format: {file_format!r}; expected {FORMAT_BAM!r} or {FORMAT_CRAM!r}"
        logger.error(message)
        raise ValueError(message)
    path = Path(in_path)
    return (f"{path}.{suffix}", str(path.with_suffix(f".{suffix}")))


def missing_index_message(in_path: str, file_format: str, tried: Sequence[str]) -> str:
    """The message a user gets when an index is absent and could not be built.

    Args:
        in_path (str): The alignment.
        file_format (str): ``"bam"`` or ``"cram"``.
        tried (Sequence[str]): Every name that was looked for.

    Returns:
        str: A single line naming the file, the names tried and the command that fixes it.
    """
    names = ", ".join(tried)
    return (
        f"No index found for {file_format.upper()} {in_path}. Looked for: {names}. "
        f"The region slice is random alignment retrieval and htslib requires an index for it. "
        f"Build one with `samtools index {in_path}`, or re-run with a writable output directory "
        f"so VNtyper can build one itself."
    )


def unresolvable_reference_message(
    in_path: str,
    contig: str,
    m5: str | None,
    attempts: Sequence[tuple[str, str, str]],
) -> str:
    """The message a user gets when no candidate reference decodes the CRAM.

    Args:
        in_path (str): The CRAM.
        contig (str): The contig whose reference could not be resolved.
        m5 (str | None): That contig's ``M5`` digest from the header, or None when the
            header carries no ``M5`` (in which case no digest-based lookup can ever work
            and the message says so).
        attempts (Sequence[tuple[str, str, str]]): ``(source, path, reason)`` per
            candidate, in the order they were tried.

    Returns:
        str: A multi-line message naming the contig, its digest and every candidate with
        the reason it failed.
    """
    digest = m5 if m5 else "no M5 in the header, so no digest lookup can succeed"
    lines = [
        f"Could not resolve a reference that decodes {in_path}.",
        f"  contig: {contig}",
        f"  M5:     {digest}",
        "  candidates tried, in order:",
    ]
    for source, path, reason in attempts:
        shown = path if path else "(not set)"
        lines.append(f"    - {source}: {shown} -- {reason}")
    lines.append("  Supply one with --reference-fasta <genome.fa>.")
    return "\n".join(lines)
```

- [ ] **Step 4: Run it and watch it pass**

```bash
pytest tests/unit/test_alignment_contract.py -v
```

Expected: PASS, 10 tests.

- [ ] **Step 5: Confirm branch coverage on the new module**

```bash
pytest tests/unit/test_alignment_contract.py --cov=vntyper.scripts.alignment_contract \
  --cov-branch --cov-report=term-missing -q
```

Expected: 100%. If a line is uncovered, add the test rather than deleting the line.

- [ ] **Step 6: Commit**

```bash
git add vntyper/scripts/alignment_contract.py tests/unit/test_alignment_contract.py
git commit -m "feat(preflight): add the alignment contract and its failure messages

Refs #225 #209"
```

### Task 3: `resolve_cram_index`

**Files:**
- Modify: `vntyper/scripts/alignment_index.py`
- Test: `tests/unit/test_alignment_index.py`

**Interfaces:**
- Consumes: `index_candidate_names` from Task 2.
- Produces: `resolve_cram_index(in_cram: str | Path) -> str | None`.

**Do not touch `resolve_bam_index`.** Its BAI-only behaviour is deliberate and is pinned
by `test_a_csi_index_is_ignored_and_a_bai_is_built_instead`: the only consumer of its
return value is `extract_unmapped_from_offset.get_last_chunk_end`, which rejects any file
whose first four bytes are not `BAI\x01`. Widening it is a change to that reader first.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_alignment_index.py`:

```python
class TestResolveCramIndex:
    """CRAM index resolution: both spellings, and nothing else.

    Unlike the BAM resolver this has no downstream parser to satisfy -- the CRAM path
    uses samtools' own indexed retrieval rather than the BAI offset extractor -- so the
    only question is which of the two names exists.
    """

    def test_it_finds_the_file_dot_crai_spelling(self, tmp_path):
        cram = tmp_path / "sample.cram"
        cram.write_bytes(b"CRAM")
        index = tmp_path / "sample.cram.crai"
        index.write_bytes(b"")
        assert resolve_cram_index(cram) == str(index)

    def test_it_finds_the_stem_dot_crai_spelling(self, tmp_path):
        cram = tmp_path / "sample.cram"
        cram.write_bytes(b"CRAM")
        index = tmp_path / "sample.crai"
        index.write_bytes(b"")
        assert resolve_cram_index(cram) == str(index)

    def test_the_file_spelling_wins_when_both_exist(self, tmp_path):
        cram = tmp_path / "sample.cram"
        cram.write_bytes(b"CRAM")
        (tmp_path / "sample.cram.crai").write_bytes(b"")
        (tmp_path / "sample.crai").write_bytes(b"")
        assert resolve_cram_index(cram) == str(tmp_path / "sample.cram.crai")

    def test_it_returns_none_when_there_is_no_index(self, tmp_path):
        cram = tmp_path / "sample.cram"
        cram.write_bytes(b"CRAM")
        assert resolve_cram_index(cram) is None

    def test_a_bai_beside_a_cram_is_not_mistaken_for_an_index(self, tmp_path):
        cram = tmp_path / "sample.cram"
        cram.write_bytes(b"CRAM")
        (tmp_path / "sample.cram.bai").write_bytes(b"")
        assert resolve_cram_index(cram) is None
```

Add `resolve_cram_index` to the module's import at the top of the test file.

- [ ] **Step 2: Run it and watch it fail**

```bash
pytest tests/unit/test_alignment_index.py -k Cram -v
```

Expected: FAIL — `ImportError: cannot import name 'resolve_cram_index'`.

- [ ] **Step 3: Implement it**

Append to `vntyper/scripts/alignment_index.py`:

```python
def resolve_cram_index(in_cram: str | Path) -> str | None:
    """Find an existing **CRAI**, under either of the two names it can carry.

    ``<file>.crai`` is tried first, then ``<stem>.crai`` -- the same two spellings
    :func:`resolve_bam_index` handles for BAI, and the same two the upload endpoint and
    the Celery worker already accept (`docker/app/tasks.py`).

    This is a separate function rather than a widening of :func:`resolve_bam_index`, and
    the asymmetry is deliberate. That function is BAI-only because the value it returns is
    parsed by ``extract_unmapped_from_offset.get_last_chunk_end``, which walks the BAI
    container itself. The CRAM path has no such reader: it hands the index to samtools,
    which resolves it. So the two have different constraints and must not share an
    implementation that satisfies only one of them.

    Args:
        in_cram (str | Path): The CRAM whose index is wanted.

    Returns:
        str | None: The existing CRAI, or None when neither spelling is present. A
        ``.bai`` beside a CRAM is not an index for it and is not returned.
    """
    for candidate in index_candidate_names(in_cram, FORMAT_CRAM):
        if Path(candidate).exists():
            return candidate
    return None
```

and add to the imports at the top of the file:

```python
from vntyper.scripts.alignment_contract import FORMAT_CRAM, index_candidate_names
```

- [ ] **Step 4: Run the whole file and watch it pass**

```bash
pytest tests/unit/test_alignment_index.py -v
```

Expected: PASS, including every pre-existing `resolve_bam_index` test unchanged.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/alignment_index.py tests/unit/test_alignment_index.py
git commit -m "feat(cram): resolve an existing .crai under both of its names

Refs #225"
```

### Task 4: `alignment_preflight.py` — index resolution before the slice

**Files:**
- Create: `vntyper/scripts/alignment_preflight.py`
- Test: `tests/unit/test_alignment_preflight.py`
- Modify: `vntyper/config.json`

**Interfaces:**
- Consumes: `AlignmentPlan`, `index_candidate_names`, `missing_index_message` (Task 2);
  `resolve_bam_index`, `resolve_cram_index` (Task 3);
  `build_samtools_index_command` (existing).
- Produces:
  `ensure_index(in_path: str, file_format: str, output_dir: str | Path, output_name: str,
  config: dict, threads: int) -> str` — returns the index path, building one into
  `output_dir` when absent. Raises `ValueError` with `missing_index_message` when it
  cannot.

This task delivers index resolution only. The reference probe (Task 8) and the layout
count (Task 12) are added to this module by later tasks; each brings its own tests.

- [ ] **Step 1: Add the config surface**

In `vntyper/config.json`, add a top-level `"cram"` block and the two reference keys:

```jsonc
  "cram": {
    "allow_ambient_reference_resolution": false,
    "unmapped_scan": "indexed",
    "reference_probe_region": null
  },
  "read_layout": {
    "mixed_tolerance": 0.0
  },
```

and inside the existing `"reference_data"` object:

```jsonc
    "cram_reference_hg19": null,
    "cram_reference_hg38": null,
```

- [ ] **Step 2: Write the failing test**

`tests/unit/test_alignment_preflight.py`:

```python
"""Index resolution happens before the slice that needs it, and never beside the input.

The region slice is random alignment retrieval, so htslib requires an index for it. The
pipeline used to resolve one seven lines *after* that call for BAM and not at all for
CRAM (#225), so an alignment arriving without its index died inside
`process_bam_to_fastq` rather than being rejected at submission.
"""

from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_contract import FORMAT_BAM, FORMAT_CRAM
from vntyper.scripts.alignment_preflight import ensure_index

pytestmark = pytest.mark.unit

CONFIG = {"tools": {"samtools": "samtools"}}


class TestEnsureIndexReusesWhatIsThere:
    def test_an_existing_bai_is_reused_and_nothing_is_built(self, tmp_path):
        bam = tmp_path / "s.bam"
        bam.write_bytes(b"BAM")
        (tmp_path / "s.bam.bai").write_bytes(b"")
        with patch("vntyper.scripts.alignment_preflight.run_command") as run:
            resolved = ensure_index(str(bam), FORMAT_BAM, tmp_path, "output", CONFIG, threads=4)
        assert resolved == str(tmp_path / "s.bam.bai")
        run.assert_not_called()

    def test_an_existing_crai_is_reused_and_nothing_is_built(self, tmp_path):
        cram = tmp_path / "s.cram"
        cram.write_bytes(b"CRAM")
        (tmp_path / "s.crai").write_bytes(b"")
        with patch("vntyper.scripts.alignment_preflight.run_command") as run:
            resolved = ensure_index(str(cram), FORMAT_CRAM, tmp_path, "output", CONFIG, threads=4)
        assert resolved == str(tmp_path / "s.crai")
        run.assert_not_called()


class TestEnsureIndexBuildsIntoTheOutputDirectory:
    def test_a_missing_bai_is_built_into_output_never_beside_the_input(self, tmp_path):
        inputs = tmp_path / "in"
        inputs.mkdir()
        out = tmp_path / "out"
        out.mkdir()
        bam = inputs / "s.bam"
        bam.write_bytes(b"BAM")
        with patch("vntyper.scripts.alignment_preflight.run_command", return_value=True) as run:
            resolved = ensure_index(str(bam), FORMAT_BAM, out, "output", CONFIG, threads=4)
        assert resolved == str(out / "output_input.bam.bai")
        command = run.call_args[0][0]
        assert "-o" in command and str(out) in command
        assert str(inputs) not in command.split("-o")[1].split()[0]

    def test_a_missing_crai_is_built_into_output_too(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        cram = tmp_path / "s.cram"
        cram.write_bytes(b"CRAM")
        with patch("vntyper.scripts.alignment_preflight.run_command", return_value=True) as run:
            resolved = ensure_index(str(cram), FORMAT_CRAM, out, "output", CONFIG, threads=4)
        assert resolved == str(out / "output_input.cram.crai")
        assert "samtools index" in run.call_args[0][0]

    def test_the_build_is_threaded(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        cram = tmp_path / "s.cram"
        cram.write_bytes(b"CRAM")
        with patch("vntyper.scripts.alignment_preflight.run_command", return_value=True) as run:
            ensure_index(str(cram), FORMAT_CRAM, out, "output", CONFIG, threads=8)
        assert "-@ 8" in run.call_args[0][0]


class TestEnsureIndexFailsWithAMessageThatNamesWhatIsMissing:
    def test_a_failed_build_names_the_file_and_both_index_spellings(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        cram = tmp_path / "s.cram"
        cram.write_bytes(b"CRAM")
        with patch("vntyper.scripts.alignment_preflight.run_command", return_value=False):
            with pytest.raises(ValueError) as excinfo:
                ensure_index(str(cram), FORMAT_CRAM, out, "output", CONFIG, threads=4)
        message = str(excinfo.value)
        assert str(cram) in message
        assert "s.cram.crai" in message and "s.crai" in message
        assert "samtools index" in message
```

- [ ] **Step 3: Run it and watch it fail**

```bash
pytest tests/unit/test_alignment_preflight.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'vntyper.scripts.alignment_preflight'`.

- [ ] **Step 4: Write the module**

`vntyper/scripts/alignment_preflight.py`:

```python
"""
vntyper/scripts/alignment_preflight.py

Module Purpose:
---------------
Prove, before any stage runs, that the run has everything its input requires.

This is the only new module in the input-contract work that shells out.
:mod:`vntyper.scripts.alignment_contract`,
:mod:`vntyper.scripts.reference_resolution` and :mod:`vntyper.scripts.read_layout` hold
the decisions; this holds the samtools calls that turn a question about a file into an
answer.

The ordering is the point. The region slice is random alignment retrieval and htslib
requires an index for it, but index handling used to sit *below* that call for BAM and
was absent entirely for CRAM (#225). Everything here runs before the first stage
command is built, so a missing piece is a submission-time message rather than a stage
log an hour into a cohort run.
"""

from __future__ import annotations

import logging
from pathlib import Path

from vntyper.scripts.alignment_contract import (
    FORMAT_BAM,
    FORMAT_CRAM,
    index_candidate_names,
    missing_index_message,
)
from vntyper.scripts.alignment_index import resolve_bam_index, resolve_cram_index
from vntyper.scripts.command_builders import build_samtools_index_command
from vntyper.scripts.utils import run_command

logger = logging.getLogger(__name__)

_RESOLVERS = {FORMAT_BAM: resolve_bam_index, FORMAT_CRAM: resolve_cram_index}


def ensure_index(
    in_path: str,
    file_format: str,
    output_dir: str | Path,
    output_name: str,
    config: dict,
    threads: int,
) -> str:
    """Return an index for this alignment, building one if there is none.

    An existing index is reused under either of the two names it can carry. A missing one
    is built into ``output_dir`` and never beside the input: the input directory holds
    patient data and is routinely mounted read-only (#162, #210).

    Args:
        in_path (str): The alignment.
        file_format (str): ``"bam"`` or ``"cram"``.
        output_dir (str | Path): The run's output directory -- where a built index goes.
        output_name (str): Base name for produced files.
        config (dict): Configuration; ``config["tools"]["samtools"]`` is read.
        threads (int): Thread count for the build.

    Returns:
        str: Path to a usable index.

    Raises:
        ValueError: If the format is unknown, or if no index exists and one could not be
            built. The message names the alignment, both index spellings and the command
            that fixes it.
    """
    resolver = _RESOLVERS.get(file_format)
    if resolver is None:
        message = f"unknown alignment format: {file_format!r}"
        logger.error(message)
        raise ValueError(message)

    existing = resolver(in_path)
    if existing is not None:
        logger.info(f"Reusing existing {file_format.upper()} index: {existing}")
        return existing

    suffix = Path(in_path).suffix
    built = str(Path(output_dir) / f"{output_name}_input{suffix}.{'bai' if file_format == FORMAT_BAM else 'crai'}")
    command = build_samtools_index_command(
        samtools_path=config["tools"]["samtools"],
        bam_file=in_path,
        output_bai=built,
        threads=threads,
    )
    log_file = str(Path(output_dir) / f"{output_name}_preflight_index.log")
    logger.info(f"No index beside the input; building one into the output directory: {command}")
    if not run_command(str(command), log_file, critical=False):
        message = missing_index_message(in_path, file_format, list(index_candidate_names(in_path, file_format)))
        logger.error(message)
        raise ValueError(message)
    return built
```

- [ ] **Step 5: Run it and watch it pass**

```bash
pytest tests/unit/test_alignment_preflight.py -v
```

Expected: PASS, 6 tests. `test_the_build_is_threaded` will still FAIL — that is correct;
Task 5 adds `threads` to `build_samtools_index_command`. Mark it `xfail` only if you are
doing Task 5 next in the same session; otherwise do Task 5 first and return here.

- [ ] **Step 6: Commit**

```bash
git add vntyper/scripts/alignment_preflight.py tests/unit/test_alignment_preflight.py vntyper/config.json
git commit -m "feat(preflight): resolve or build the alignment index before the slice

Closes #225"
```

### Task 5: command-builder performance and shape changes (P1, P2, P3, P4)

**Files:**
- Modify: `vntyper/scripts/command_builders.py`
- Test: `tests/unit/test_command_builders.py`

**Interfaces:**
- Produces:
  - `build_samtools_index_command(*, samtools_path, bam_file, output_bai=None, threads=1) -> str`
  - `build_samtools_slice_command(*, samtools_path, in_bam, output_bam, region=None,
    bed_file=None, cram_ref_option="", threads=1, index_output=True) -> str`
  - `build_cram_unmapped_indexed_command(*, samtools_path, in_bam, unmapped_bam, threads,
    cram_ref_option="") -> str`
  - `build_cram_reference_probe_command(*, samtools_path, in_bam, region=None,
    bed_file=None, cram_ref_option="", threads=1) -> str`

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_command_builders.py`:

```python
class TestSliceThreadingAndIndexing:
    """P1 and P3 from the spec: the slice was single-threaded, and in the default
    (non-fast) mode it indexed a file the merge then overwrote."""

    def test_the_slice_is_threaded(self):
        command = build_samtools_slice_command(
            samtools_path="samtools", in_bam="/i/s.bam", output_bam="/o/s.bam",
            region="chr1:1-2", threads=8,
        )
        assert "-@ 8" in command

    def test_indexing_the_slice_can_be_skipped(self):
        command = build_samtools_slice_command(
            samtools_path="samtools", in_bam="/i/s.bam", output_bam="/o/s.bam",
            region="chr1:1-2", threads=8, index_output=False,
        )
        assert "index" not in command
        assert "&&" not in command

    def test_indexing_the_slice_is_still_the_default(self):
        command = build_samtools_slice_command(
            samtools_path="samtools", in_bam="/i/s.bam", output_bam="/o/s.bam", region="chr1:1-2",
        )
        assert "samtools index" in command


class TestIndexThreading:
    """P2: `samtools index` ran single-threaded everywhere."""

    def test_the_index_is_threaded(self):
        assert "-@ 4" in build_samtools_index_command(samtools_path="samtools", bam_file="/o/s.bam", threads=4)

    def test_one_thread_emits_no_flag_so_existing_command_strings_are_unchanged(self):
        assert build_samtools_index_command(samtools_path="samtools", bam_file="/o/s.bam") == "samtools index /o/s.bam"


class TestIndexedUnmappedFetch:
    """P4: the whole-file scan is O(file); the indexed unplaced fetch is O(index)."""

    def test_it_fetches_the_unplaced_reads_through_the_index(self):
        command = build_cram_unmapped_indexed_command(
            samtools_path="samtools", in_bam="/i/s.cram", unmapped_bam="/o/u.bam",
            threads=8, cram_ref_option="-T /r/g.fa",
        )
        assert "'*'" in command
        assert "-f 12" in command
        assert "-T /r/g.fa" in command
        assert "|" not in command, "the indexed form is one command; a pipe means the stream form leaked in"

    def test_it_carries_no_reference_fragment_for_a_reference_free_cram(self):
        command = build_cram_unmapped_indexed_command(
            samtools_path="samtools", in_bam="/i/s.cram", unmapped_bam="/o/u.bam", threads=1,
        )
        assert "-T" not in command


class TestReferenceProbe:
    """The probe must have the same shape as the slice it authorises (spec section 3.7):
    a -P-less probe passes on a chr1-only reference the real slice then rejects."""

    def test_the_probe_carries_the_same_fetch_pairs_flag_as_the_slice(self):
        probe = build_cram_reference_probe_command(
            samtools_path="samtools", in_bam="/i/s.cram", region="chr1:1-2", cram_ref_option="-T /r/g.fa",
        )
        slice_command = build_samtools_slice_command(
            samtools_path="samtools", in_bam="/i/s.cram", output_bam="/o/s.bam",
            region="chr1:1-2", cram_ref_option="-T /r/g.fa",
        )
        assert "-P" in probe, "a probe without -P authorises a reference the slice will reject"
        assert "-P" in slice_command

    def test_the_probe_counts_and_writes_nothing(self):
        probe = build_cram_reference_probe_command(
            samtools_path="samtools", in_bam="/i/s.cram", region="chr1:1-2", cram_ref_option="-T /r/g.fa",
        )
        assert "-c" in probe
        assert "-o" not in probe
        assert "-b" not in probe

    def test_the_probe_accepts_a_bed_file_the_way_the_slice_does(self):
        probe = build_cram_reference_probe_command(
            samtools_path="samtools", in_bam="/i/s.cram", bed_file="/o/r.bed", cram_ref_option="",
        )
        assert "-L /o/r.bed" in probe
```

- [ ] **Step 2: Run them and watch them fail**

```bash
pytest tests/unit/test_command_builders.py -k "Threading or Indexed or Probe" -v
```

Expected: FAIL — `TypeError: ... got an unexpected keyword argument 'threads'` and
`ImportError` for the two new builders.

- [ ] **Step 3: Implement**

In `vntyper/scripts/command_builders.py`:

Add a helper near `quote_path`:

```python
def _thread_flag(threads: int) -> str:
    """The ``-@ N`` fragment, or empty for a single thread.

    Kept empty at ``threads=1`` so that command strings for existing single-threaded
    call sites are byte-identical to the ones this module produced before threading was
    added -- the test suite compares several of them literally.

    Args:
        threads (int): Thread count.

    Returns:
        str: ``"-@ N "`` with a trailing space, or ``""``.
    """
    return f"-@ {quote_path(threads)} " if threads and int(threads) > 1 else ""
```

Change `build_samtools_index_command` to take `threads: int = 1` and emit
`_thread_flag(threads)` after `index`. Change `build_samtools_slice_command` to take
`threads: int = 1` and `index_output: bool = True`; emit `_thread_flag(threads)` after
`view -P -b`, and return without the `&& samtools index ...` half when
`index_output` is False. Document why:

```python
        index_output (bool): Whether to index the slice. Default True. The merge in
            non-fast mode ``os.replace``s the merged BAM onto this exact path and then
            re-indexes it, so indexing here first is work that is thrown away on every
            default run. Fast mode keeps it: there is no merge, and ``pipeline.py``
            consumes ``output_sliced.bam`` with its index.
```

Add the two new builders:

```python
def build_cram_unmapped_indexed_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    unmapped_bam: str | Path,
    threads: int,
    cram_ref_option: str = "",
) -> str:
    """Fetch the unplaced unmapped reads through the CRAM index.

    :func:`build_cram_unmapped_filter_command` decodes and re-encodes the **entire**
    file as SAM text to find flag-12 reads. The BAM path never did that -- it reads the
    BAI's last chunk offset and seeks straight to the unplaced block -- and the
    asymmetry is why a CRAM run looks like a hang where the equivalent BAM run does not
    (#178). ``'*'`` is htslib's region syntax for the unplaced reads, so this is the
    CRAM equivalent of that seek.

    Measured on a fixture with 200 000 placed and 4 000 unplaced reads: both forms
    recover the same 4 000 reads.

    Args:
        samtools_path (str): samtools invocation from config.
        in_bam (str | Path): The input CRAM. Must be indexed -- the preflight guarantees it.
        unmapped_bam (str | Path): Where the unmapped reads are written.
        threads (int): Thread count.
        cram_ref_option (str): Pre-formatted flag fragment, interpolated verbatim.

    Returns:
        str: One samtools invocation. No pipe, so no ``pipefail`` prefix is needed.
    """
    return (
        f"{samtools_path} view -b -f 12 {cram_ref_option} {_thread_flag(threads)}"
        f"{quote_path(in_bam)} '*' -o {quote_path(unmapped_bam)}"
    )


def build_cram_reference_probe_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    region: str | None = None,
    bed_file: str | Path | None = None,
    cram_ref_option: str = "",
    threads: int = 1,
) -> str:
    """Decode the run's own target region and count the records, writing nothing.

    This authorises a candidate reference, and it therefore has **the same shape as the
    slice it authorises**: the same ``-P``, the same region or BED, the same reference
    fragment. That is not a stylistic choice. ``-P`` (``--fetch-pairs``) chases mates
    onto other contigs, so a probe without it succeeds against a chr1-only reference
    that the real slice then fails on -- measured, spec section 3.7. A probe that can
    pass where the real command fails is a fix that masks a failure mode instead of
    removing one.

    Args:
        samtools_path (str): samtools invocation from config.
        in_bam (str | Path): The CRAM being probed.
        region (str | None): Region string, as for the slice.
        bed_file (str | Path | None): BED file, as for the slice. Wins if both are given.
        cram_ref_option (str): The candidate's ``-T`` fragment, interpolated verbatim.
        threads (int): Thread count.

    Returns:
        str: ``samtools view -c -P ...``.

    Raises:
        ValueError: If neither a region nor a BED file is given, matching the slice.
    """
    if bed_file is not None:
        target = f"-L {quote_path(bed_file)}"
    elif region is not None:
        target = quote_path(region)
    else:
        message = "build_cram_reference_probe_command needs either a region or a bed_file"
        logger.error(message)
        raise ValueError(message)

    return f"{samtools_path} view -c -P {cram_ref_option} {_thread_flag(threads)}{quote_path(in_bam)} {target}"
```

- [ ] **Step 4: Run the whole builder suite**

```bash
pytest tests/unit/test_command_builders.py -v
```

Expected: PASS. If a pre-existing test compares a command string literally and now
fails on whitespace, fix the **test's** expectation only where `_thread_flag` legitimately
changed the string; do not add compensating spaces to the builder.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/command_builders.py tests/unit/test_command_builders.py
git commit -m "perf(samtools): thread the slice and index, and stop indexing a slice the merge overwrites

Adds the indexed unplaced-read fetch and the reference probe, both of which
the CRAM contract needs.

Refs #178 #209 #225"
```

### Task 6: wire the plan into `process_bam_to_fastq`

**Files:**
- Modify: `vntyper/scripts/fastq_bam_processing.py:118-236`
- Modify: `vntyper/scripts/pipeline.py:225-335`
- Test: `tests/unit/test_fastq_bam_processing.py`, `tests/unit/test_input_tree_is_never_written.py`

**Interfaces:**
- Consumes: `AlignmentPlan` (Task 2), `ensure_index` (Task 4), the builders (Task 5).
- Produces: `process_bam_to_fastq(..., plan: AlignmentPlan)` — the `file_format` argument
  is removed; the plan carries it.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_fastq_bam_processing.py`:

```python
class TestTheSliceIsAuthorisedBeforeItRuns:
    def test_the_index_is_resolved_before_the_slice_command_is_built(self, tmp_path):
        """#225: the slice is random retrieval and needs an index; resolution used to sit
        seven lines below it for BAM and be absent for CRAM."""
        calls: list[str] = []

        def record(command, *_args, **_kwargs):
            calls.append(str(command))
            return True

        plan = AlignmentPlan(
            input_path=str(tmp_path / "s.cram"),
            file_format=FORMAT_CRAM,
            index_path=str(tmp_path / "s.cram.crai"),
            reference_path=str(tmp_path / "g.fa"),
            reference_source="--reference-fasta",
            uncovered_contigs=(),
            layout="paired",
        )
        with patch("vntyper.scripts.fastq_bam_processing.run_command", side_effect=record):
            process_bam_to_fastq(
                in_bam=plan.input_path, output=tmp_path, output_name="output", threads=4,
                config=CONFIG, bed_file=None, plan=plan, fast_mode=True,
            )
        slice_calls = [c for c in calls if " view " in c and " -b " in c]
        assert slice_calls, "no slice command was built"
        assert "-T" in slice_calls[0], "the slice ran without the proven reference"

    def test_the_cram_slice_carries_the_reference_from_the_plan(self, tmp_path):
        """#209: cram_ref_option was the empty string unconditionally."""
        calls: list[str] = []
        plan = AlignmentPlan(
            input_path=str(tmp_path / "s.cram"), file_format=FORMAT_CRAM,
            index_path=str(tmp_path / "s.cram.crai"), reference_path="/r/genome.fa",
            reference_source="config cram_reference_hg38", uncovered_contigs=(), layout="paired",
        )
        with patch("vntyper.scripts.fastq_bam_processing.run_command",
                   side_effect=lambda c, *a, **k: (calls.append(str(c)), True)[1]):
            process_bam_to_fastq(
                in_bam=plan.input_path, output=tmp_path, output_name="output", threads=4,
                config=CONFIG, bed_file=None, plan=plan, fast_mode=True,
            )
        assert all("-T /r/genome.fa" in c for c in calls if " view " in c and "-c -P" not in c)
```

- [ ] **Step 2: Run it and watch it fail**

```bash
pytest tests/unit/test_fastq_bam_processing.py -k "authorised or reference_from_the_plan" -v
```

Expected: FAIL — `TypeError: process_bam_to_fastq() got an unexpected keyword argument 'plan'`.

- [ ] **Step 3: Change the signature and delete the inline index block**

In `fastq_bam_processing.py`:

- Replace the `file_format="bam"` parameter with `plan: AlignmentPlan`.
- Replace `cram_ref_option = ""` (line 133) with `cram_ref_option = plan.cram_ref_option`.
- Delete the `bam_bai = resolve_bam_index(in_bam)` block (lines 167-177) — the preflight
  did that. Use `plan.index_path` in the call to `extract_unmapped_reads_from_offset`.
- Replace every `file_format` reference with `plan.file_format`.
- In the CRAM branch, select the scan by config:

```python
            scan = config.get("cram", {}).get("unmapped_scan", "indexed")
            if scan == "indexed":
                command_filter = build_cram_unmapped_indexed_command(
                    samtools_path=samtools_path, in_bam=in_bam, unmapped_bam=unmapped_bam,
                    threads=threads, cram_ref_option=cram_ref_option,
                )
            else:
                command_filter = build_cram_unmapped_filter_command(
                    samtools_path=samtools_path, in_bam=in_bam, unmapped_bam=unmapped_bam,
                    threads=threads, cram_ref_option=cram_ref_option,
                )
```

- Pass `threads=threads` and `index_output=fast_mode` to `build_samtools_slice_command`.

In `pipeline.py`, build the plan before the region is resolved and pass it to both call
sites. Remove `file_format="cram"` from the CRAM call.

- [ ] **Step 4: Run the affected suites**

```bash
pytest tests/unit/test_fastq_bam_processing.py tests/unit/test_pipeline.py \
       tests/unit/test_input_tree_is_never_written.py -v
```

Expected: PASS.

- [ ] **Step 5: Check the file did not grow**

```bash
wc -l vntyper/scripts/fastq_bam_processing.py
```

Expected: ≤ 612. The deleted index block should more than pay for the plan plumbing. If
it grew, extract the CRAM branch into a helper rather than accepting it.

- [ ] **Step 6: Run the standard block and commit**

```bash
git add vntyper/scripts/fastq_bam_processing.py vntyper/scripts/pipeline.py tests/unit/
git commit -m "refactor(pipeline): consume the alignment plan instead of rediscovering it

Refs #225 #209 #178"
```

---

## Wave 2 — parallel worktrees

These four tasks share no files. Each gets its own worktree, created with
`superpowers:using-git-worktrees`, branched from the tip of Wave 1.

**#209 and #178 are one task (Task 7), not two.** Both edit the CRAM branch of
`process_bam_to_fastq` and both edit `command_builders.py`; running them as separate
agents produces a merge conflict in the file the whole milestone turns on.

### Task 7: the reference contract and the unreachable hang (#209 + #178)

**Files:**
- Create: `vntyper/scripts/reference_resolution.py`
- Modify: `vntyper/scripts/alignment_preflight.py`, `vntyper/scripts/cli_parser.py`,
  `vntyper/scripts/cli_handlers.py`, `vntyper/scripts/pipeline.py`
- Test: `tests/unit/test_reference_resolution.py`,
  `tests/unit/test_alignment_preflight.py`, `tests/unit/test_ref_path_is_pinned.py`

**Interfaces:**
- Consumes: `AlignmentPlan`, `unresolvable_reference_message` (Task 2);
  `build_cram_reference_probe_command` (Task 5).
- Produces:
  - `reference_candidates(cli_reference: str | None, config: dict, assembly: str,
    header_ur: str | None) -> list[tuple[str, str | None]]` — `(source, path)` in order.
  - `uncovered_contigs(header_contigs: Sequence[str], fasta_contigs: Sequence[str]) -> tuple[str, ...]`
  - `resolve_reference(...) -> tuple[str | None, str | None, tuple[str, ...]]` in
    `alignment_preflight` — `(path, source, uncovered)`.
  - `pin_reference_resolution(config: dict) -> None` in `alignment_preflight`.

- [ ] **Step 1: Write the failing tests for candidate ordering**

`tests/unit/test_reference_resolution.py`:

```python
"""Which reference a CRAM run tries, in which order, and what it does not cover.

The winner is decided by probing (see alignment_preflight); this module only decides
what gets probed and in what order.
"""

import pytest

from vntyper.scripts.reference_resolution import reference_candidates, uncovered_contigs

pytestmark = pytest.mark.unit

CONFIG = {
    "reference_data": {
        "cram_reference_hg38": "/r/genome.hg38.fa",
        "bwa_reference_hg38": "reference/alignment/chr1.hg38.fa",
    }
}


class TestCandidateOrder:
    def test_the_cli_flag_is_tried_first(self):
        candidates = reference_candidates("/u/mine.fa", CONFIG, "hg38", header_ur=None)
        assert candidates[0] == ("--reference-fasta", "/u/mine.fa")

    def test_the_config_cram_reference_is_tried_before_the_shipped_chr1(self):
        sources = [source for source, _ in reference_candidates(None, CONFIG, "hg38", header_ur=None)]
        assert sources.index("config cram_reference_hg38") < sources.index("config bwa_reference_hg38")

    def test_the_header_ur_is_tried_last(self):
        candidates = reference_candidates(None, CONFIG, "hg38", header_ur="/lab/genome.fa")
        assert candidates[-1] == ("header UR:", "/lab/genome.fa")

    def test_unset_candidates_are_still_listed_so_the_failure_message_can_name_them(self):
        candidates = reference_candidates(None, {"reference_data": {}}, "hg38", header_ur=None)
        assert ("--reference-fasta", None) in candidates

    def test_an_unknown_assembly_does_not_raise_because_a_config_may_be_replaced_wholesale(self):
        assert reference_candidates(None, {}, "hg99", header_ur=None)


class TestUncoveredContigs:
    def test_a_full_reference_covers_everything(self):
        assert uncovered_contigs(["chr1", "chr2"], ["chr1", "chr2", "chr3"]) == ()

    def test_a_chr1_only_reference_reports_what_it_misses(self):
        assert uncovered_contigs(["chr1", "chr2", "chrX"], ["chr1"]) == ("chr2", "chrX")

    def test_it_is_order_stable_so_the_warning_text_is_deterministic(self):
        assert uncovered_contigs(["chrX", "chr2"], ["chr1"]) == ("chrX", "chr2")
```

- [ ] **Step 2: Run and watch fail; implement `reference_resolution.py`**

```bash
pytest tests/unit/test_reference_resolution.py -v
```

Expected: FAIL, `ModuleNotFoundError`. Then write the module with `reference_candidates`
returning, in order: `("--reference-fasta", cli_reference)`,
`(f"config cram_reference_{assembly}", ...)`, `(f"config bwa_reference_{assembly}", ...)`,
`("header UR:", header_ur)`. Read every config key with `.get` and shipped defaults
(trap 2). `uncovered_contigs` returns header contigs absent from the FASTA, in header
order.

- [ ] **Step 3: Write the failing test for the probe**

Append to `tests/unit/test_alignment_preflight.py`:

```python
class TestReferenceIsProvenNotInferred:
    def test_the_first_candidate_that_decodes_wins(self, tmp_path):
        outcomes = {"/r/a.fa": False, "/r/b.fa": True}

        def probe(command, *_a, **_k):
            return next((ok for path, ok in outcomes.items() if path in command), False)

        with patch("vntyper.scripts.alignment_preflight.run_command", side_effect=probe):
            path, source, _ = resolve_reference(
                in_path="/i/s.cram", candidates=[("cli", "/r/a.fa"), ("config", "/r/b.fa")],
                region="chr1:1-2", bed_file=None, config=CONFIG, threads=1,
                output_dir=tmp_path, output_name="output", header_contigs=["chr1"],
            )
        assert (path, source) == ("/r/b.fa", "config")

    def test_a_cram_that_needs_no_reference_resolves_to_none_and_still_runs(self, tmp_path):
        with patch("vntyper.scripts.alignment_preflight.run_command", return_value=True):
            path, source, _ = resolve_reference(
                in_path="/i/s.cram", candidates=[], region="chr1:1-2", bed_file=None,
                config=CONFIG, threads=1, output_dir=tmp_path, output_name="output",
                header_contigs=["chr1"],
            )
        assert path is None and source is None

    def test_no_candidate_decoding_raises_a_message_naming_every_one(self, tmp_path):
        with patch("vntyper.scripts.alignment_preflight.run_command", return_value=False):
            with pytest.raises(ValueError) as excinfo:
                resolve_reference(
                    in_path="/i/s.cram", candidates=[("--reference-fasta", None), ("config", "/r/b.fa")],
                    region="chr1:1-2", bed_file=None, config=CONFIG, threads=1,
                    output_dir=tmp_path, output_name="output", header_contigs=["chr1"],
                )
        message = str(excinfo.value)
        assert "--reference-fasta" in message and "/r/b.fa" in message

    def test_the_probe_uses_the_run_s_own_region_not_a_hardcoded_one(self, tmp_path):
        seen: list[str] = []
        with patch("vntyper.scripts.alignment_preflight.run_command",
                   side_effect=lambda c, *a, **k: (seen.append(str(c)), True)[1]):
            resolve_reference(
                in_path="/i/s.cram", candidates=[("cli", "/r/a.fa")], region="chr1:155184000-155194000",
                bed_file=None, config=CONFIG, threads=1, output_dir=tmp_path,
                output_name="output", header_contigs=["chr1"],
            )
        assert "chr1:155184000-155194000" in seen[0]
        assert "-P" in seen[0]
```

- [ ] **Step 4: Implement `resolve_reference`**

Try each candidate whose path is non-None and exists; run
`build_cram_reference_probe_command` through `run_command(..., critical=False)`; the
first that returns True wins. Before returning, compute `uncovered_contigs` from the
winner's `.fai` (or by running `samtools faidx`-free parsing of the FASTA headers) and
log a warning naming them when non-empty. When no candidate wins **and** the no-reference
probe also fails, raise `ValueError(unresolvable_reference_message(...))`. When the
candidate list is empty, probe once with no `-T`; success means a reference-free CRAM.

- [ ] **Step 5: Write the failing test for the pinned `REF_PATH`**

`tests/unit/test_ref_path_is_pinned.py`:

```python
"""No samtools call VNtyper makes may block on a network reference fetch.

Measured on the pinned samtools 1.20 / htslib 1.23: with REF_PATH pointed at an endpoint
that accepts TCP and never responds, `samtools view -h` produced no output, no error and
did not exit -- exit 124 at a 25 s timeout. That is #178's reported shape verbatim. A
local-only REF_PATH fails in milliseconds instead, naming the reference.

Pinning is belt and braces on top of always passing -T: -T bypasses REF_PATH entirely
(measured: 0.00 s against the same unresponsive endpoint), but a call site this milestone
did not think of would otherwise still be able to reach the network.
"""

import os

import pytest

from vntyper.scripts.alignment_preflight import pin_reference_resolution

pytestmark = pytest.mark.unit


def test_an_unset_ref_path_is_pinned_to_a_local_only_value(monkeypatch):
    monkeypatch.delenv("REF_PATH", raising=False)
    pin_reference_resolution({})
    assert "http://" not in os.environ["REF_PATH"]
    assert "https://" not in os.environ["REF_PATH"]


def test_an_operator_s_network_ref_path_is_overridden_by_default(monkeypatch):
    monkeypatch.setenv("REF_PATH", "http://www.ebi.ac.uk/ena/cram/md5/%s")
    pin_reference_resolution({})
    assert "ebi.ac.uk" not in os.environ["REF_PATH"]


def test_it_can_be_opted_out_of_for_a_site_with_its_own_refget_server(monkeypatch):
    monkeypatch.setenv("REF_PATH", "http://refget.internal/%s")
    pin_reference_resolution({"cram": {"allow_ambient_reference_resolution": True}})
    assert os.environ["REF_PATH"] == "http://refget.internal/%s"


def test_a_config_replaced_wholesale_does_not_raise(monkeypatch):
    monkeypatch.delenv("REF_PATH", raising=False)
    pin_reference_resolution({"tools": {}})
    assert os.environ["REF_PATH"]
```

- [ ] **Step 6: Implement `pin_reference_resolution` and call it from the preflight**

- [ ] **Step 7: Add `--reference-fasta`**

`cli_parser.py`, in the `pipeline` subparser beside `--cram`:

```python
    parser_pipeline.add_argument(
        "--reference-fasta",
        type=str,
        help=(
            "Reference FASTA used to decode a CRAM. Only needed when the shipped chr1 "
            "reference cannot decode it -- VNtyper probes the candidates it knows about "
            "first and names this flag if none of them works."
        ),
    )
```

Thread it through `cli_handlers.handle_pipeline` into `run_pipeline`.

- [ ] **Step 8: Run everything and commit**

```bash
pytest tests/unit/test_reference_resolution.py tests/unit/test_alignment_preflight.py \
       tests/unit/test_ref_path_is_pinned.py -v
make format && make test-unit
git add -A
git commit -m "feat(cram): prove the reference at submission and pin REF_PATH

Closes #209
Closes #178"
```

### Task 8: `detect_naming_convention` (#165)

**Files:**
- Modify: `vntyper/scripts/chromosome_utils.py:137-197`, `vntyper/config.json`
- Test: `tests/unit/test_chromosome_utils.py`

**Interfaces:** unchanged signature —
`detect_naming_convention(contig_names: list[str], config: dict | None = None) -> str`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_chromosome_utils.py`:

```python
class TestDecoyAndAltContigsDoNotPoisonTheDenominator:
    """#165: contigs that match no convention were counted in `total`, so a header whose
    every primary contig agreed could still return 'unknown'."""

    def test_the_issue_s_own_example_resolves_to_ucsc(self):
        primaries = [f"chr{n}" for n in list(range(1, 23)) + ["X", "Y", "M"]]
        decoys = [f"chr{i}_gl{i:06d}_random" for i in range(68)]
        contigs = primaries + decoys
        assert len(contigs) == 93
        assert detect_naming_convention(contigs) == "ucsc"

    def test_a_header_that_is_genuinely_ambiguous_is_still_unknown(self):
        assert detect_naming_convention(["chr1", "2", "NC_000003.12"]) == "unknown"

    def test_a_header_with_no_recognisable_contig_is_unknown_not_a_zero_division(self):
        assert detect_naming_convention(["scaffold_1", "scaffold_2"]) == "unknown"

    def test_the_threshold_comes_from_config_not_from_a_literal(self):
        contigs = ["chr1", "chr2", "1", "2", "X"]
        assert detect_naming_convention(contigs, {"assembly_detection": {"naming_convention_threshold": 0.9}}) == "unknown"
```

- [ ] **Step 2: Run and watch fail**

```bash
pytest tests/unit/test_chromosome_utils.py -k Decoy -v
```

Expected: FAIL — `assert 'unknown' == 'ucsc'`.

- [ ] **Step 3: Fix the denominator**

Replace `total = len(contig_names)` with the count of contigs that matched *some*
convention, and return `"unknown"` when that count is zero. Read the threshold with
`config.get("assembly_detection", {}).get("naming_convention_threshold", 0.5)`. Keep the
existing warning log and extend it to report how many contigs were classified.

- [ ] **Step 4: Run the whole file plus the guard's suite**

```bash
pytest tests/unit/test_chromosome_utils.py tests/unit/test_assembly_guard.py \
       tests/unit/test_region_utils.py tests/unit/test_builders.py -v
```

Expected: PASS. `test_builders.py:263-268` asserts a specific `unknown` result — check
whether that header is *genuinely* ambiguous or was only ambiguous because of the bug. If
the latter, change the test deliberately and say so in the commit message.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/chromosome_utils.py vntyper/config.json tests/unit/
git commit -m "fix(chromosome): score naming convention over classified contigs only

Closes #165"
```

### Task 9: single-end support (#161)

**Files:**
- Create: `vntyper/scripts/read_layout.py`
- Modify: `vntyper/scripts/alignment_preflight.py`, `vntyper/scripts/pipeline.py`,
  `vntyper/scripts/alignment_processing.py`, `vntyper/scripts/kestrel_genotyping.py`,
  `vntyper/scripts/command_builders.py` (fastp and bwa single-end forms)
- Test: `tests/unit/test_read_layout.py`, plus additions to the four modified suites

**Interfaces:**
- Produces:
  - `classify_layout(paired: int, total: int, mixed_tolerance: float) -> str` returning
    `"paired"`, `"single"`, `"mixed"` or `"empty"`.
  - `route_fastqs(layout: str, r1: str, r2: str, other: str, single: str,
    sizes: Mapping[str, int]) -> tuple[tuple[str, ...], tuple[str, ...]]` returning
    `(consumed, unconsumed_nonempty)`.

- [ ] **Step 1: Write the failing test**

`tests/unit/test_read_layout.py`:

```python
"""Which reads a run consumes, and the guarantee that it consumes all of them.

`samtools fastq -0` collects every read with READ1 and READ2 both unset -- which is every
single-end read -- and both call sites in pipeline.py used to bind that path to `_`. A
single-end BAM therefore had 100% of its reads written to output_other.fastq.gz and then
discarded, and the run continued on two empty FASTQs (#161).
"""

import pytest

from vntyper.scripts.read_layout import classify_layout, route_fastqs

pytestmark = pytest.mark.unit


class TestClassifyLayout:
    def test_all_reads_flagged_paired_is_paired(self):
        assert classify_layout(paired=1000, total=1000, mixed_tolerance=0.0) == "paired"

    def test_no_read_flagged_paired_is_single(self):
        assert classify_layout(paired=0, total=1000, mixed_tolerance=0.0) == "single"

    def test_a_mixture_is_named_rather_than_coerced(self):
        assert classify_layout(paired=600, total=1000, mixed_tolerance=0.0) == "mixed"

    def test_a_tolerance_lets_a_few_stray_reads_be_absorbed(self):
        assert classify_layout(paired=999, total=1000, mixed_tolerance=0.01) == "paired"
        assert classify_layout(paired=1, total=1000, mixed_tolerance=0.01) == "single"

    def test_no_reads_at_all_is_empty_not_single(self):
        assert classify_layout(paired=0, total=0, mixed_tolerance=0.0) == "empty"


class TestRouteFastqs:
    def test_paired_consumes_r1_and_r2(self):
        consumed, unconsumed = route_fastqs(
            "paired", "r1.gz", "r2.gz", "o.gz", "s.gz",
            sizes={"r1.gz": 900, "r2.gz": 900, "o.gz": 0, "s.gz": 0},
        )
        assert consumed == ("r1.gz", "r2.gz")
        assert unconsumed == ()

    def test_single_consumes_the_other_file_where_samtools_puts_unpaired_reads(self):
        consumed, unconsumed = route_fastqs(
            "single", "r1.gz", "r2.gz", "o.gz", "s.gz",
            sizes={"r1.gz": 0, "r2.gz": 0, "o.gz": 1800, "s.gz": 0},
        )
        assert consumed == ("o.gz",)
        assert unconsumed == ()

    def test_a_non_empty_file_nothing_consumes_is_reported_never_dropped(self):
        _, unconsumed = route_fastqs(
            "paired", "r1.gz", "r2.gz", "o.gz", "s.gz",
            sizes={"r1.gz": 900, "r2.gz": 900, "o.gz": 4200, "s.gz": 0},
        )
        assert unconsumed == ("o.gz",)

    def test_empty_consumes_nothing_and_reports_nothing(self):
        assert route_fastqs("empty", "r1.gz", "r2.gz", "o.gz", "s.gz",
                            sizes={"r1.gz": 0, "r2.gz": 0, "o.gz": 0, "s.gz": 0}) == ((), ())
```

- [ ] **Step 2: Run and watch fail, then implement `read_layout.py`**

- [ ] **Step 3: Write the failing tests for the single-end command forms**

Append to `tests/unit/test_command_builders.py` and
`tests/unit/test_kestrel_genotyping.py`:

```python
def test_fastp_accepts_a_single_input(self):
    command = build_fastp_command(
        fastp_path="fastp", threads=4, fastq_1="/o/single.fastq.gz", fastq_2=None,
        output="/o", output_name="output", compression_level=4, qualified_quality_phred=20,
        dup_calc_accuracy=3, length_required=50, disable_adapter_trimming=False, deduplication=False,
    )
    assert "--in2" not in command and "--out2" not in command
    assert "--in1 /o/single.fastq.gz" in command


def test_bwa_aligns_a_single_fastq(self):
    command = build_bwa_align_sort_command(
        bwa_path="bwa", samtools_path="samtools", threads=4, reference="/r/chr1.fa",
        fastq1="/o/single.fastq.gz", fastq2=None, sorted_bam="/o/out.bam",
    )
    assert command.count("fastq.gz") == 1


def test_kestrel_accepts_a_single_fastq(self):
    command = construct_kestrel_command(
        kmer_size=20, kestrel_path="/k/kestrel.jar", reference_vntr="/r/vntr.fa",
        output_dir="/o", fastq_1="/o/single.fastq.gz", fastq_2=None, vcf_out="/o/out.vcf",
        java_path="java", java_memory="12g", max_align_states=30, max_hap_states=30,
        log_level="INFO", sample_name="s",
    )
    assert "/o/single.fastq.gz" in command
    assert command.count("fastq.gz") == 1


def test_kestrel_still_rejects_no_fastq_at_all(self):
    with pytest.raises(ValueError, match="FASTQ input"):
        construct_kestrel_command(
            kmer_size=20, kestrel_path="/k/kestrel.jar", reference_vntr="/r/vntr.fa",
            output_dir="/o", fastq_1=None, fastq_2=None, vcf_out="/o/out.vcf",
            java_path="java", java_memory="12g", max_align_states=30, max_hap_states=30,
            log_level="INFO", sample_name="s",
        )
```

- [ ] **Step 4: Implement.** Change the three builders and
  `construct_kestrel_command`'s guard from `if not fastq_1 or not fastq_2` to
  `if not fastq_1`. Change `align_and_sort_fastq`'s `fastq2` to `Path | None`.

- [ ] **Step 5: Wire the routing into `pipeline.py`.** Replace both
  `fastq1, fastq2, _, _ = process_bam_to_fastq(...)` bindings with the four-value form
  and a `route_fastqs` call; raise `ValueError` naming the file and its size when
  `unconsumed_nonempty` is non-empty; log the layout.

- [ ] **Step 6: Derive single-end fixtures**

Add to `scripts/make_cram_fixtures.py` a `--single-end` mode producing a BAM whose reads
carry no `0x1` flag, derived from a cohort BAM with
`samtools view -h <in> | awk 'BEGIN{OFS="\t"} /^@/{print;next} {$2=and($2,compl(0xC1)); print}'`
— or, more robustly, with pysam clearing `is_paired`, `is_read1`, `is_read2`,
`mate_is_unmapped` and `is_proper_pair`. Write it beside the CRAM fixtures and register
it in `tests/test_data_config.json`.

- [ ] **Step 7: Run everything and commit**

```bash
pytest tests/unit/test_read_layout.py tests/unit/test_command_builders.py \
       tests/unit/test_kestrel_genotyping.py tests/unit/test_alignment_processing.py \
       tests/unit/test_pipeline.py -v
make format && make test-unit
git add -A
git commit -m "feat(input): support single-end alignments and never discard a read set

Closes #161"
```

### Task 10: fixture deriver (`--all` and the reference-dependent fixture)

**Files:**
- Modify: `scripts/make_cram_fixtures.py`
- Test: `tests/unit/test_make_cram_fixtures.py`

- [ ] **Step 1:** Write a failing test asserting that with no `--all`, the deriver
  selects only samples named in `tests/test_data_config.json`, and that `--all` selects
  every discovered sample.
- [ ] **Step 2:** Run it, watch it fail.
- [ ] **Step 3:** Add the `--all` flag; default the sample list to the declared fixtures.
- [ ] **Step 4:** Add `build_reference_dependent_fixture` producing a
  reference-compressed CRAM plus a copy of its reference, so #209's path is testable.
- [ ] **Step 5:** `pytest tests/unit/test_make_cram_fixtures.py -v`, then commit.

---

## Wave 3 — integration

### Task 11: golden-cohort gate and the performance measurement

- [ ] **Step 1: Baseline before the change**

```bash
git stash && time python scripts/golden_cohort_gate.py --json /tmp/baseline.json ; git stash pop
```

- [ ] **Step 2: Derive the CRAM fixtures**

```bash
make cram-fixtures
```

Do **not** pass `--allow-matrix-drift`.

- [ ] **Step 3: Run the gate on the branch**

```bash
time python scripts/golden_cohort_gate.py --json /tmp/after.json
```

Expected: ~14 min on 32 cores. Compare genotypes sample by sample; any difference is a
blocker, not a note.

- [ ] **Step 4: Prove the two CRAM scans agree**

```bash
VNTYPER_CRAM_UNMAPPED_SCAN=stream python scripts/golden_cohort_gate.py --json /tmp/stream.json
python - <<'PY'
import json
a = json.load(open("/tmp/after.json")); b = json.load(open("/tmp/stream.json"))
assert a["samples"] == b["samples"], "indexed and stream scans disagree"
print("indexed == stream across", len(a["samples"]), "samples")
PY
```

- [ ] **Step 5: Record both wall-clocks in the spec's §5b Gate paragraph and commit.**

### Task 12: full gate and PR

- [ ] **Step 1:** `make check-all`
- [ ] **Step 2:** `make patch-coverage` — must be ≥ 80%
- [ ] **Step 3:** `make ci-local` **only if** `.github/workflows/` was touched
- [ ] **Step 4:** Run the Codex adversarial gate on the final diff (see spec §10).
- [ ] **Step 5:** Bump `vntyper/version.py`, `CITATION.cff` and
  `docs/about/changelog.md` to 2.0.10 — patch only.
- [ ] **Step 6:** Open one PR against `main` with `Closes #213`, `Closes #225`,
  `Closes #209`, `Closes #178`, `Closes #165`, `Closes #161` in the body.

---

## Self-review

**Spec coverage.** §4.1 seam → Tasks 2, 4, 6. §4.2 reference contract → Task 7. §4.3(a)
`REF_PATH` → Task 7 step 5. §4.3(b) bounded scan → Task 5 + Task 6 step 3. §4.4 no silent
discard → Task 9. §5b P1–P3 → Task 5; P4 → Task 5 + Task 6. §5 #213 → Task 1. §5 #225 →
Tasks 3, 4, 6. §5 #165 → Task 8. §5 #161 → Task 9. Fixture deriver → Task 10. §7 A-PERF-1
and A-178-2 → Task 11. §7 A-ALL-1 → Task 12.

**Placeholders.** Task 7 steps 4 and 6, and Task 9 steps 2, 4 and 5, describe the
implementation in prose rather than showing it. That is deliberate: those are the three
places where the shape depends on what Wave 1 actually landed, and the tests above them
are complete and executable, so the implementer has an exact target. Task 10's steps are
prose for the same reason — the deriver's internals are not yet read.

**Type consistency.** `AlignmentPlan` field names are identical in Tasks 2, 4, 6 and 7.
`cram_ref_option` is a property in Task 2 and consumed as one in Task 6.
`index_candidate_names` returns a 2-tuple in Task 2 and is unpacked as an iterable in
Tasks 3 and 4. `classify_layout`/`route_fastqs` names match between Task 9's test and its
interface block.
