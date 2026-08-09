# Milestone 4 — CRAM and Input Robustness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps
> use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give VNtyper an input contract, so a BAM or CRAM run proves every prerequisite
before any processing stage does work, naming exactly what is missing, and no input is
ever silently discarded. The produced-read layout is necessarily decided after lossless
conversion and before any FASTQ consumer, because the four output files are its evidence.

**Architecture:** An owned alignment-preflight boundary runs first in `run_pipeline`. It
opens the regular alignment and materialises a run-local view of that exact inode (an
atomic procfd view, with verified same-filesystem hardlink fallback) before quickcheck,
header, target or coverage-region reads. An operator BED is likewise copied byte-exactly
into a pipeline-owned file before proof. The boundary then builds the co-located fresh
index, chooses the unmapped scan from `samtools idxstats`, and proves the reference for
every consumer shape. The generated BAI/CRAI is opened before publication; post-preflight
random-access consumers use its retained procfd or verified-hardlink path explicitly.
The boundary returns the target and frozen `AlignmentPlan` whose alignment and index
bindings remain alive through coverage and are released by `AlignmentPlan.close()`. Read layout is
**not** part of that plan — it is decided after conversion, from the FASTQs that were
actually written. Decisions live in pure modules; only the preflight shells out.

**Tech Stack:** Python 3.10+ (floor; CI matrix 3.10–3.13), samtools 1.20 / htslib 1.23,
pytest, ruff, mypy, Docker + Celery for the web layer.

**Spec:** [`2026-08-08-milestone-4-cram-input-robustness-spec.md`](2026-08-08-milestone-4-cram-input-robustness-spec.md).

**Status:** H1/H2/H3 and their final follow-ups are integrated and independently reviewed.
The final-candidate A-PERF-1 rerun is recorded; the combined no-HIGH review remains unchecked.

> **Read §3 of the spec before writing any CRAM code.** Six intuitions about htslib in
> this area are wrong and the measurements say which. In particular: an index in the
> output directory is invisible (§3.9); `-P` and `-c` cannot be combined (§3.10); `-P`
> fetches cross-contig mates (§3.7); a resolvable `UR:` silently rescues a candidate that
> should fail (§3.10); `'*'` loses *placed* read-unmapped records and `idxstats` detects
> them (§3.13), but zero placed records alone do not prove literal-`'*'` completeness
> (§3.11).

## Global Constraints

From AGENTS.md. Every task's requirements implicitly include these.

- **Run pytest from the repo root.** `tests/parametrization.py` opens
  `tests/test_data_config.json` by relative path at collection time.
- **Activate the `vntyper` conda env and put its bin FIRST on `PATH`.** The base env and
  `~/.local/bin` shadow it, and `make check-all` otherwise reports on the wrong tools.
- **Every new test file needs `pytestmark = pytest.mark.unit`** after the imports, or CI
  silently never runs it. `tests/unit/test_marker_hygiene.py` enforces this.
- **A test that must import `app.tasks` or `app.main` goes under `tests/unit/web/`** —
  that package's `conftest.py` is what puts `docker/` on `sys.path`.
- **Ruff only. Line length 120, double quotes, `target-version = "py310"`.**
- **Code must run on Python 3.10.** PEP 604 (`X | None`) and PEP 585 (`list[str]`) are
  available; `match`, `Self` and PEP 695 generics are not.
- **Google-style docstrings** (`Args:` / `Returns:` / `Raises:`) on public functions.
- **`logger = logging.getLogger(__name__)` after the imports.** Never `logging.info(...)`,
  never `basicConfig`. f-strings in log calls are the established style.
- **No custom exception classes.** `logger.error(msg)` then `raise ValueError(msg)` /
  `RuntimeError(msg)` with the same message.
- **Never lower the coverage floor.** `make patch-coverage` must be ≥ 80% on changed lines.
- **`--config-path` replaces the whole config** (trap 2). Every new key is read with
  `.get` and the shipped default.
- **`run_command` runs one string under bash with `shell=True`** (trap 9). Quoting happens
  only in `command_builders.py`, via `quote_path`. Tool invocations are *not* quoted
  (`config["tools"]` holds command prefixes); paths, regions and sample names *are*.
- **Keep files under ~650 LOC.** Re-measured final implementation values are
  `fastq_bam_processing.py` **625**, `pipeline.py` **665**,
  `kestrel_genotyping.py` **865**, and `alignment_preflight.py` **644** lines. Further
  edits to the two over-limit files extract the touched responsibility rather than grow it.
- **Touch a file, add tests for it**, including for functions you only modify.
- **If you touch `.github/workflows/`, run `make ci-local`,** not just `make check-all`.
- **Never claim tests pass without showing the command output.**

**Standard verification block** — every task ends with this:

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
| `vntyper/scripts/alignment_contract.py` | create | `AlignmentPlan`, including lifetime-owned alignment/index bindings and `stable_index_path`; index names per format; failure messages and the `preflight_error.json` payload. Contract module, not pure. |
| `vntyper/scripts/alignment_binding.py` | create | Opened-alignment binding; procfd view, verified hardlink fallback, generated-index ownership, exact cleanup and FD lifetime. |
| `vntyper/scripts/alignment_index_binding.py` | create | Open the fresh temporary BAI/CRAI before publication and retain its exact inode through procfd or verified-hardlink view. |
| `vntyper/scripts/alignment_consumer_commands.py` | create | Build post-preflight slice/unmapped commands from the complete proven plan. Pure. |
| `vntyper/scripts/samtools_command_fragments.py` | create | Shared quoted thread/reference/explicit-index command fragments. Pure. |
| `vntyper/scripts/reference_resolution.py` | create | Validated ordered reference candidates; known or unavailable FASTA contig coverage. Pure. |
| `vntyper/scripts/idxstats_parsing.py` | create | Parse an `idxstats` table and decide the scan. Fails closed. Pure. |
| `vntyper/scripts/read_layout.py` | create | The layout verdict from FASTQ record counts, and which files are consumed vs stranded. Pure. |
| `vntyper/scripts/alignment_preflight.py` | create | Orchestrate the view, fresh index, scan choice and reference proof; pin `REF_PATH`. |
| `vntyper/scripts/preflight_command_io.py` | create | The one captured-command executor: atomic logs plus optional process-group deadline, kill and reap. |
| `vntyper/scripts/kestrel_command.py` | create | `construct_kestrel_command`, extracted without behavioural change from `kestrel_genotyping.py` and re-exported for existing callers; final source is 865 lines. |
| `vntyper/scripts/alignment_index.py` | modify | Add `resolve_any_index`; `resolve_bam_index` unchanged. |
| `vntyper/scripts/command_builders.py` | modify | Thread slice + index; optional slice indexing; indexed unplaced fetch/count; reference probe; quoted `reference_path` and explicit retained `index_path` on slice, indexed filter and **depth**. |
| `vntyper/scripts/fastq_bam_processing.py` | modify | Consume `AlignmentPlan`. Must end **below** 649 lines. |
| `vntyper/scripts/chromosome_utils.py` | modify | Classified-only denominator, strict majority, config threshold. |
| `vntyper/scripts/pipeline.py` | modify | Call the preflight; route the layout; give coverage the reference and the view. |
| `vntyper/scripts/pipeline_cleanup.py` | create | Close plans without replacing an active primary pipeline outcome. |
| `vntyper/scripts/pipeline_coverage.py` | create | The final coverage consumer; closes the plan after it returns. |
| `vntyper/scripts/pipeline_inputs.py` | create | Input selection, output/archive ownership and trailing-separator normalization. |
| `vntyper/scripts/cli_logging_safety.py` | create | Exclusive regular-file/alias checks before a pipeline log is opened. |
| `vntyper/scripts/archive_safety.py` | create | Shared CLI/web descriptor-anchored archive creation and stale-archive cleanup. |
| `docker/app/archive_delivery.py` | create | Descriptor-bound individual/cohort downloads and worker-owned cohort-input snapshots. |
| `vntyper/scripts/alignment_processing.py` | modify | `align_and_sort_fastq` accepts a single FASTQ. |
| `vntyper/scripts/cli_parser.py`, `cli_handlers.py` | modify | `--reference-fasta`; allow `--fastq1` without `--fastq2`. |
| `docker/app/tasks.py`, `docker/app/main.py` | modify | Read `preflight_error.json`; surface its message in job status. |
| `vntyper/config.json` | modify | The spec §6 surface. |
| `scripts/make_cram_fixtures.py` | modify | Declared fixtures by default, `--all`, reference-dependent and single-end fixtures. |

---

## Wave 0 — #213 (COMPLETE)

Landed in `18ae985` and `949c9de`. Six `conda run` sites: four in `docker/entrypoint.sh`,
two in `docker/app/tasks.py` (`build_vntyper_command` and the inline cohort launcher).
`tests/unit/web/test_conda_run_streams_output.py` scans both files as source text.
`make check-all` green, 2846 unit tests.

---

## Wave 1 — the keystone

Everything in Wave 2 consumes what this wave builds. Do these in order; do not start
Wave 2 until `make check-all` is green on the tip of Wave 1.

### Task 1: `alignment_contract.py`

**Files:** create `vntyper/scripts/alignment_contract.py`,
`tests/unit/test_alignment_contract.py`

**Interfaces produced:**
- `FORMAT_BAM = "bam"`, `FORMAT_CRAM = "cram"`
- `INDEX_SUFFIXES: dict[str, tuple[str, ...]]` — `bam` → `("bai", "csi")`,
  `cram` → `("crai", "csi")`. The CRAM CSI support is based on a direct measurement
  with the pinned samtools 1.20/htslib stack; CRAI remains first.
- `index_candidate_names(in_path, file_format, *, bai_only: bool = False) -> tuple[str, ...]`
  — `<file>.<suffix>` then `<stem>.<suffix>` for each suffix in order.
- `AlignmentPlan` frozen dataclass: `input_path: str`, `view_path: str`,
  `file_format: str`, `index_path: str`, `reference_path: str | None`,
  `reference_source: str`, `uncovered_contigs: tuple[str, ...]`, `unmapped_scan: str`,
  `binding: AlignmentBinding | None`; `close() -> None` removes the owned view then
  releases the descriptor exactly once.
- `missing_index_message(in_path, file_format, tried) -> str` — final round-3 wording
  names failure to build the fresh run-local index; it does not tell users to create a
  supplied index the implementation deliberately will not trust.
- `unresolvable_reference_message(in_path, contig, m5, attempts) -> str`
- `preflight_error_payload(code: str, message: str, attempts) -> dict` — what
  `preflight_error.json` holds.

**`AlignmentPlan` has no `layout` field.** Layout cannot be known before conversion
(spec §4.4); putting it here is the circular dependency Codex round 1 found.

- [x] **Step 1: write the failing tests.** Cover: both index spellings per suffix; CSI
      included for BAM and excluded under `bai_only=True`; CRAM never offers `.bai`;
      an unknown format raises `ValueError` matching `"unknown alignment format"`;
      `missing_index_message` contains the path, every supplied name inspected and the
      `samtools index` diagnostic command;
      `unresolvable_reference_message` contains the contig, the M5 (or the words
      `no M5` when it is None), every candidate source, every candidate path and every
      reason; the dataclass is frozen; and:

```python
def test_the_plan_exposes_a_reference_path_not_a_preformatted_shell_fragment():
    plan = _plan(reference_path="/r/genome ref.fa")
    assert plan.reference_path == "/r/genome ref.fa"
    assert not hasattr(plan, "cram_ref_option")


def test_the_error_payload_carries_no_absolute_worker_paths_beyond_the_input():
    payload = preflight_error_payload("reference_unresolved", "msg", [("cli", None, "not supplied")])
    assert set(payload) == {"code", "message", "candidates"}
```

- [x] **Step 2:** `pytest tests/unit/test_alignment_contract.py -v` → FAIL
      (`ModuleNotFoundError`).
- [x] **Step 3:** implement. `AlignmentPlan` carries only `reference_path`; command
      builders receive that path and own quoting and option construction (§5 #209).
- [x] **Step 4:** `pytest tests/unit/test_alignment_contract.py -v` → PASS.
- [x] **Step 5:** `pytest ... --cov=vntyper.scripts.alignment_contract --cov-branch
      --cov-report=term-missing` → 100%. Add tests, never delete lines, to reach it.
- [x] **Step 6:** commit — `feat(preflight): add the alignment contract and its messages`.

### Task 2: `idxstats_parsing.py` — fail-closed scan selection

**Files:** create `vntyper/scripts/idxstats_parsing.py`,
`tests/unit/test_idxstats_parsing.py`

**Interfaces produced:**
- `SCAN_INDEXED = "indexed"`, `SCAN_STREAM = "stream"`
- `parse_idxstats(text: str) -> tuple[int, int] | None` — `(placed_unmapped, unplaced)`,
  or `None` when the table is not well-formed.
- `choose_scan(configured: str, idxstats_text: str | None, exit_ok: bool, *,
  indexed_count_text: str | None = None, indexed_count_exit_ok: bool = False) ->
  tuple[str, str]` — `(scan, reason)`. Indexed is authorised only when the table is valid,
  its placed sum is zero and the exact literal-`'*'` flag-4 count equals its terminal-`*`
  count. Forced indexed raises when valid evidence proves that consumer incomplete.

Spec §4.5: **anything not well-formed selects `stream`.** An unparsable table must never
read as "column 4 summed to zero".

- [x] **Step 1: write the failing tests.**

```python
GOOD = "chr1\t20000\t600\t50\n*\t0\t0\t80\n"
CLEAN = "chr1\t20000\t600\t0\n*\t0\t0\t80\n"

class TestParse:
    def test_it_returns_placed_and_unplaced_counts(self):
        assert parse_idxstats(GOOD) == (50, 80)

    def test_a_clean_file_has_no_placed_unmapped_reads(self):
        assert parse_idxstats(CLEAN) == (0, 80)

    @pytest.mark.parametrize("bad", [
        "", "   \n", "chr1\t20000\t600\n", "chr1\t20000\t600\t50\textra\n",
        "chr1\t20000\tsix\t50\n*\t0\t0\t80\n", "chr1\t20000\t600\t-1\n*\t0\t0\t80\n",
        "chr1\t20000\t600\t50\n", "*\t0\t0\t1\n*\t0\t0\t2\n",
    ])
    def test_anything_malformed_is_rejected_rather_than_guessed(self, bad):
        assert parse_idxstats(bad) is None

class TestChooseScan:
    def test_auto_picks_indexed_only_when_both_counts_prove_nothing_would_be_lost(self):
        assert choose_scan(
            "auto", CLEAN, exit_ok=True,
            indexed_count_text="80\n", indexed_count_exit_ok=True,
        )[0] == SCAN_INDEXED

    def test_auto_falls_back_to_stream_when_placed_unmapped_reads_exist(self):
        scan, reason = choose_scan("auto", GOOD, exit_ok=True)
        assert scan == SCAN_STREAM and "50" in reason

    def test_auto_falls_back_to_stream_on_a_malformed_table(self):
        assert choose_scan("auto", "garbage", exit_ok=True)[0] == SCAN_STREAM

    def test_auto_falls_back_to_stream_when_idxstats_itself_failed(self):
        assert choose_scan("auto", CLEAN, exit_ok=False)[0] == SCAN_STREAM

    def test_forcing_indexed_where_reads_would_be_lost_raises_rather_than_dropping(self):
        with pytest.raises(ValueError, match="50"):
            choose_scan("indexed", GOOD, exit_ok=True)

    def test_zero_placed_but_incomplete_literal_star_fetch_is_not_authorised(self):
        scan, reason = choose_scan(
            "auto", CLEAN, exit_ok=True,
            indexed_count_text="20\n", indexed_count_exit_ok=True,
        )
        assert scan == SCAN_STREAM and "20" in reason and "80" in reason

    def test_stream_is_always_allowed_because_it_is_never_lossy(self):
        assert choose_scan("stream", GOOD, exit_ok=True)[0] == SCAN_STREAM
```

- [x] **Step 2:** run → FAIL. **Step 3:** implement. **Step 4:** run → PASS.
- [x] **Step 5:** 100% branch coverage on the module. **Step 6:** commit.

### Task 3: `resolve_any_index`

**Files:** modify `vntyper/scripts/alignment_index.py`, `tests/unit/test_alignment_index.py`

**Interfaces produced:** `resolve_any_index(in_path, file_format) -> str | None` — the
first existing name from `index_candidate_names(...)`. The final round-3 design uses this
only to enumerate/protect supplied patient artifacts; §3.15 forbids consuming any supplied
index because a valid wrong-sample index can return an empty slice with exit 0.

**Do not touch `resolve_bam_index`.** It remains a BAI-only legacy
preflight/protection contract, and
`test_a_csi_index_is_ignored_and_a_bai_is_built_instead` pins that. Production indexed
recovery no longer consumes its result through the optional offset parser (§3.20);
widening preflight index policy is a separate change.

- [x] **Step 1: write the failing tests** — `.cram.crai` and `.crai` both found, file
      spelling wins; CRAM CSI is found only after both CRAI spellings, while a `.bai`
      beside a CRAM is not an index for it; for BAM a `.csi` is found by
      `resolve_any_index` but still **not** by `resolve_bam_index`; `None` when nothing
      exists.
- [x] **Step 2:** run → FAIL (`ImportError`). **Step 3:** implement. **Step 4:** run the
      whole file → PASS, every pre-existing test unchanged. **Step 5:** commit.

### Task 4: command-builder changes (P1–P4 and the new forms)

**Files:** modify `vntyper/scripts/command_builders.py`,
`tests/unit/test_command_builders.py`; create `tests/unit/test_nonfast_slice_union.py`

**Interfaces produced** — note every one takes `reference_path`, **not** a preformatted
fragment, and quotes it here (spec §5 #209):

```python
_thread_flag(threads: int) -> str                      # "-@ N " or "" at threads<=1
_reference_flag(reference_path: str | Path | None) -> str   # "-T <quoted> " or ""

build_samtools_index_command(*, samtools_path, bam_file, output_bai=None, threads=1) -> str
build_samtools_slice_command(*, samtools_path, in_bam, output_bam, region=None,
                             bed_file=None, reference_path=None, index_path=None, threads=1,
                             index_output=True, exclude_unmapped=False) -> str
build_cram_unmapped_indexed_command(*, samtools_path, in_bam, unmapped_bam, threads,
                                    reference_path=None, index_path=None) -> str
build_cram_unmapped_filter_command(*, samtools_path, in_bam, unmapped_bam, threads,
                                   reference_path=None) -> str     # signature change
build_cram_reference_probe_command(*, samtools_path, in_bam, region=None, bed_file=None,
                                   reference_path=None, threads=1) -> str
build_samtools_idxstats_command(*, samtools_path, in_bam, threads=1) -> str
build_samtools_unmapped_indexed_count_command(*, samtools_path, in_bam, threads=1) -> str
build_samtools_depth_command(*, samtools_path, threads, region, bam_file,
                             coverage_output, reference_path=None,
                             index_path=None) -> str   # signature change
```

- [x] **Step 1: write the failing tests.** The four that matter most:

```python
def test_the_probe_never_combines_P_with_c_because_samtools_rejects_that():
    probe = build_cram_reference_probe_command(
        samtools_path="samtools", in_bam="/o/view.cram", region="chr1:1-2", reference_path="/r/g.fa")
    assert "-P" in probe
    assert " -c" not in probe, "samtools: The options -P and -c cannot be combined"
    assert "-o /dev/null" in probe


def test_the_probe_has_the_same_shape_as_the_slice_it_authorises():
    kwargs = dict(samtools_path="samtools", in_bam="/o/view.cram", reference_path="/r/g.fa")
    probe = build_cram_reference_probe_command(bed_file="/o/r.bed", **kwargs)
    sliced = build_samtools_slice_command(output_bam="/o/s.bam", bed_file="/o/r.bed", **kwargs)
    for flag in ("-P", "-L /o/r.bed", "-T /r/g.fa"):
        assert flag in probe and flag in sliced


def test_depth_carries_the_reference_so_coverage_does_not_die_on_a_cram():
    command = build_samtools_depth_command(
        samtools_path="samtools", threads=4, region="chr1:1-2",
        bam_file="/o/view.cram", coverage_output="/o/d.txt", reference_path="/r/g.fa")
    assert "--reference /r/g.fa" in command


def test_a_reference_path_with_a_space_is_quoted_not_split():
    command = build_samtools_slice_command(
        samtools_path="samtools", in_bam="/o/v.cram", output_bam="/o/s.bam",
        region="chr1:1-2", reference_path="/r/my genome.fa")
    assert "-T '/r/my genome.fa'" in command
```

Plus: slice threaded (`-@ 8`); `index_output=False` emits no `&&` and no `index`;
`index_output=True` still the default; index threaded; `threads=1` emits no `-@` so
existing literal comparisons hold; the indexed fetch contains `'*'`, `-f 4` and **no**
pipe; `build_samtools_idxstats_command` carries no `-T` (spec §3.13: idxstats needs no
reference). `exclude_unmapped=True` emits `-F 4`; the default remains byte-identical for
fast mode and direct slice consumers.

- [x] **Step 2:** run → FAIL. **Step 3:** implement.
- [x] **Step 4:** `pytest tests/unit/test_command_builders.py tests/unit/test_nonfast_slice_union.py -v` → PASS. Where a
      pre-existing test compares a command string literally, fix the *expectation* only
      where `_thread_flag`/`_reference_flag` legitimately changed it; never pad the
      builder to match a stale string.
- [x] **Step 5:** commit — `perf(samtools): thread the slice and index, and give every CRAM
      command its reference`.

### Task 5: `alignment_preflight.py`

**Files:** create `vntyper/scripts/alignment_preflight.py`,
`vntyper/scripts/alignment_binding.py`,
`vntyper/scripts/reference_resolution.py`, `vntyper/scripts/reference_uri_policy.py`,
`vntyper/scripts/preflight_input_io.py`, `vntyper/scripts/alignment_preflight_logs.py`,
`vntyper/scripts/reference_resolution_environment.py`,
focused alignment-preflight unit modules, `tests/unit/test_reference_resolution.py`,
`tests/unit/test_reference_uri_policy.py`, `tests/unit/test_ref_path_is_pinned.py`, and
`tests/unit/test_preflight_input_reads.py`;
modify `vntyper/scripts/pipeline_alignment.py`, `vntyper/config.json`

**Interfaces produced:**
- `build_alignment_view(in_path, output_dir, output_name, file_format, config, threads,
  *, bai_only=False) -> tuple[str, str, AlignmentBinding]` → `(view_path, index_path,
  binding)`; it opens the input, publishes a procfd view or verified same-filesystem
  hardlink fallback, and returns the still-open binding to the plan.
- `choose_unmapped_scan(view_path, config, threads, output_dir, output_name, *, file_format) -> str`
- `resolve_reference(view_path, candidates, region, bed_file, config, threads, output_dir, output_name, header_contigs, m5, *, coverage_region=None) -> tuple[str | None, str, tuple[str, ...]]`
- `pin_reference_resolution(config) -> str | None` and
  `restore_reference_resolution(previous: str | None) -> None`
- `run_preflight(...) -> AlignmentPlan`
- `capture_command(command, log_file, cwd=None, *, timeout_seconds=None) -> tuple[bool, str]`
  — the stdout-capturing seam `idxstats` needs and the process-group deadline reference
  probes require; `run_command` returns only a bool and streams stdout to a log.
- `ordered_reference_candidates(config, reference_assembly, reference_fasta)` and
  `uncovered_reference_contigs(header_contigs, reference_contigs)` — pure candidate-policy
  and reference-coverage decisions consumed by the I/O preflight.
- `header_reference_scheme(value)`, `ref_path_remote_scheme(value)`,
  `allow_ambient_reference_resolution(config)`, `first_remote_header_reference(header)`
  and `enforce_header_reference_policy(header, *, allow_ambient)` — pure context-specific
  URI and strict-boolean policy used by local `REF_PATH` validation and the owned
  post-header boundary.

**Config** (spec §6) — add `bam.unmapped_scan` (`"auto"`),
`cram.allow_ambient_reference_resolution` (false), `cram.local_ref_path`
(`"%2s/%2s/%s"`), `cram.reference_probe_timeout_seconds` (120), `cram.unmapped_scan` (`"auto"`),
`cram.reference_candidate_order`, `reference_data.cram_reference_hg19`/`_hg38` (null),
`utils.preflight_text_max_bytes` (1048576),
`assembly_detection.naming_convention_threshold` (0.5) and `primary_contig_patterns`.
There is **no** `read_layout` block.

The default `allow_ambient_reference_resolution=false` carries the no-network-block
guarantee. `true` is an informed operator opt-in that may block later samtools stages;
only reference probes remain bounded by `reference_probe_timeout_seconds`. Task 5 does
not add a general stage timeout (spec §9).

- [x] **Step 1: write the failing tests for the view** (spec §4.5):

```python
class TestAlignmentView:
    def test_a_supplied_index_is_never_trusted_and_a_fresh_one_is_built(self, tmp_path): ...
    def test_the_index_is_built_beside_the_view_never_beside_the_input(self, tmp_path): ...
    def test_a_valid_wrong_sample_index_cannot_authorise_an_empty_slice(self, tmp_path): ...
    def test_the_input_directory_is_untouched_even_when_it_is_read_only(self, tmp_path): ...
    def test_the_procfd_view_identifies_the_opened_input_after_path_replacement(self, tmp_path): ...
    def test_procfs_failure_uses_only_a_verified_same_filesystem_hardlink(self, tmp_path): ...
    def test_close_removes_only_the_exact_owned_view_before_releasing_the_descriptor(self, tmp_path): ...
```

- [x] **Step 2:** run → FAIL. **Step 3:** implement `AlignmentBinding` and
      `build_alignment_view`: open the regular input first, atomically publish the procfd
      view, then fall back only to a same-filesystem hardlink whose identity matches the
      opened descriptor. Do not restore the superseded pathname-reuse design.
- [x] **Step 4:** write and run the failing tests for `choose_unmapped_scan` for both BAM
      and CRAM (delegates to
      `idxstats_parsing.choose_scan`, feeding it `capture_command`'s output and exit) and
      for `resolve_reference`:

```python
def test_explicit_candidates_are_probed_before_the_no_reference_one(self): ...
def test_the_no_reference_candidate_is_last_and_is_recorded_as_htslib_resolved(self):
    """spec 4.2: a no--T probe cannot distinguish 'reference-free' from 'UR: resolved it',
    so the source string must not claim the former."""
def test_the_first_candidate_that_decodes_wins(self): ...
def test_no_candidate_decoding_raises_naming_every_one_with_its_reason(self): ...
def test_the_probe_uses_the_runs_own_region_and_bed_not_a_hardcoded_one(self): ...
def test_the_winner_also_decodes_the_independent_coverage_region(self): ...
def test_stream_mode_additionally_proves_a_complete_whole_file_decode(self): ...
def test_fast_mode_skips_scan_policy_idxstats_and_whole_stream_probe(self): ...
def test_every_reference_probe_has_a_bounded_process_group_deadline(self): ...
def test_a_reference_not_covering_every_header_contig_logs_a_warning(self, caplog): ...
```

- [x] **Step 5:** write and run the failing tests for CRAM-only `REF_PATH` in
      `tests/unit/test_ref_path_is_pinned.py`: an unset `REF_PATH` is pinned to a
      local-only value; an operator's `http://` value is overridden by default; the
      override is skipped when `allow_ambient_reference_resolution` is true; and
      **`restore_reference_resolution` puts back the previous value, including "unset", in
      nested CRAM scopes and after `BaseException`** (spec §4.5 — `run_pipeline` is called
      in-process by tests). Because `os.environ` is process-global, hold one re-entrant
      lock across the complete pin/restore lifetime: same-thread nesting remains LIFO,
      overlapping CRAM threads wait, and validation failure releases the lock. Pin the warning contract:
      the opt-in may block later stages even though every reference probe remains bounded.
      Assert BAM and FASTQ neither validate nor mutate this policy.
      Add pure URI-policy tests for local/relative/`file://` values, anchored remote
      schemes with and without `//`, every duplicate `UR` field, colon-separated
      `REF_PATH`, strict boolean waiver validation and path-free contig/scheme rejection.
      In the owned boundary, test a remote CRAM header fails under
      `REFERENCE_POLICY_FAILURE` immediately after header read, before
      assembly/target/`capture_command`, and persists that exact curated diagnostic;
      only the actual boolean `true` remains the waiver. The real A-178-5 integration
      rewrites the purpose CRAM with a remote-first/local-last duplicate localhost `UR`,
      bypasses localhost proxies in both environment-variable spellings, waits for the
      listener backlog and asserts no connection, captured command or stage artifact;
      retain the ambient A-178-1 blackhole test unchanged.
- [x] **Step 6:** implement, then run the focused URI, REF_PATH, boundary and real CRAM
      integration suites → PASS. Pre-probe reference/BED/FAI reads reject FIFOs without
      waiting and apply the configured text-byte bound; do not describe `O_NONBLOCK` as a
      timeout for a stalled regular file on FUSE or NFS (spec §4.5, A-INPUT-1).
- [x] **Step 7:** commit — `feat(preflight): prove the index, the scan and the reference
      before any stage runs`.

### Task 6: wire the plan into the pipeline

**Files:** modify `vntyper/scripts/fastq_bam_processing.py`, `vntyper/scripts/pipeline.py`;
tests in `tests/unit/test_fastq_bam_command_wiring.py`, `tests/unit/test_pipeline.py`,
`tests/unit/test_input_tree_is_never_written.py` and
`tests/integration/test_read_only_alignment_preflight.py`

- [x] **Step 1: write the failing tests** — the slice runs against `plan.view_path`, not
      `in_bam`; the slice carries `-T` when the plan has a reference; the CRAM unmapped
      command is chosen by `plan.unmapped_scan`; `index_output=False` in non-fast mode
      (P3) and `True` in fast mode; the non-fast slice passes
      `exclude_unmapped=True` so complete recovery is merged as a disjoint union, while
      fast mode passes `False`; **coverage runs against the view with the plan's
      reference and the exact independently preflighted coverage region**; the input tree
      is byte-identical afterwards for both formats. A real
      nonempty all-unplaced BAM selects indexed and preserves its exact five-QNAME
      multiset (§3.20).
- [x] **Step 2:** run → FAIL. **Step 3:** implement:
      replace the `file_format="bam"` parameter with `plan: AlignmentPlan`; delete the
      inline `resolve_bam_index` block (lines ~167-177); route every post-preflight
      random-access builder through the plan helpers and `plan.stable_index_path`, while
      retaining `plan.index_path` only as the owned public/provenance pathname;
      replace `cram_ref_option = ""` with the plan's reference; branch the CRAM scan on
      `plan.unmapped_scan`. In `pipeline.py`, invoke the owned input-alignment preflight
      boundary before header/region resolution; the boundary resolves and materialises the
      exact BED internally, calls `run_preflight`, and returns both BED and plan. Pass the
      plan to both call sites, and change
      `input_bam = Path(cram)` at `:468` to the view plus the plan's reference. Route
      cleanup through `close_alignment_plan(plan, preserve_primary=...)`, so a cleanup
      failure is observable on a successful run but never replaces a primary failure.
      Use `calculate_alignment_coverage(...)` as the final consumer: it consumes the
      plan's view/reference and calls `plan.close()` only after coverage returns.
      Pass `exclude_unmapped=not fast_mode` to the slice builder. For indexed BAM and
      CRAM plans, call `build_cram_unmapped_indexed_command` (BAM has
      `reference_path=None`); remove only the production import/call of the optional BAI
      offset extractor, leaving its legacy module and parser tests intact.
      Mutation-test the lifetime by atomically replacing `plan.index_path` with a valid
      wrong-sample BAI after preflight: the bound `-X` slice and depth results must remain
      29,736 records and 5,001 rows, while removing `-X` reproduces the silent zero-record
      result.
- [x] **Step 4:** run the three suites → PASS.
- [x] **Step 5:** `wc -l vntyper/scripts/fastq_bam_processing.py` → **must be < 649**.
      If not, extract the CRAM branch into a helper rather than accepting it.
- [x] **Step 6:** standard block, then commit.

---

## Wave 2 — parallel worktrees

Branch each from the tip of Wave 1, created with `superpowers:using-git-worktrees`.
**Tasks 7 and 8 are one worktree**: both edit the CRAM branch and `command_builders.py`.

### Task 7: finish #209 + #178 (CLI surface, web transport, fixtures)

- [x] **Step 1:** add `--reference-fasta` to `cli_parser.py` beside `--cram`; thread it
      through `cli_handlers.handle_pipeline` into `run_pipeline`. Test that it reaches the
      preflight's candidate list first.
- [x] **Step 2:** write the failing test for `preflight_error.json`: the preflight writes
      it into the output directory *before* raising, and its payload has exactly
      `code`, `message`, `candidates`.
- [x] **Step 3:** write the failing tests in `tests/unit/web/` for the transport — the
      worker reads the file when the pipeline exits non-zero and stores `code`/`message`
      on the job hash; the status endpoint returns that message when present and its
      existing generic text otherwise; **and the response contains no absolute worker
      path** (the endpoint is unauthenticated and its genericity is a deliberate no-leak
      property).
- [x] **Step 4:** implement; run → PASS.
- [x] **Step 5:** build the reference-dependent CRAM fixture (Task 10) and add the
      integration tests for A-209-1/2/3. **The fixture's `UR:` target must be renamed in
      the test**, or a resolvable `UR:` silently rescues the candidate that should fail
      (spec §3.10).
- [x] **Step 6:** standard block; commit `Closes #209` and `Closes #178`.

### Task 8: `detect_naming_convention` (#165)

**Files:** modify `vntyper/scripts/chromosome_utils.py`, `vntyper/scripts/region_utils.py`,
`vntyper/config.json`; test `tests/unit/test_chromosome_utils.py`,
`tests/unit/test_region_utils.py`

- [x] **Step 1: write the failing tests.**

```python
def test_the_issues_own_93_contig_header_resolves_to_ucsc(self):
    primaries = [f"chr{n}" for n in list(range(1, 23)) + ["X", "Y", "M"]]
    decoys = [f"chr{i}_gl{i:06d}_random" for i in range(68)]
    contigs = primaries + decoys
    assert len(contigs) == 93
    assert detect_naming_convention(contigs) == "ucsc"

def test_a_two_way_fifty_fifty_split_is_unknown_not_whichever_is_checked_first(self):
    """The old code tested `>= 0.5` and returned the first passing convention, so an
    exact tie silently became UCSC. Narrowing the denominator makes ties far more
    reachable, so the tie rule has to change in the same commit."""
    assert detect_naming_convention(["chr1", "chr2", "1", "2"]) == "unknown"

def test_a_header_with_no_classifiable_contig_is_unknown_and_does_not_divide_by_zero(self):
    assert detect_naming_convention(["scaffold_1", "scaffold_2"]) == "unknown"

def test_the_threshold_comes_from_config_not_from_a_literal(self):
    contigs = ["chr1", "chr2", "chr3", "1"]
    assert detect_naming_convention(
        contigs, {"assembly_detection": {"naming_convention_threshold": 0.9}}) == "unknown"
```

Also pin the public resolver and target boundary: 50/50 and zero-classifiable headers
raise a `ValueError` naming the unresolved convention; an unambiguous header still
resolves; mutating the naming policy and replacing a file at the same path are both
observed on the next call.

- [x] **Step 2:** run → FAIL. **Step 3:** implement: denominator = contigs matching *some*
      convention; return `unknown` when that is zero; a convention wins on a **strict
      majority with no tie**; threshold from config with the 0.5 default. The public
      resolver raises before target construction when the result is `unknown`, the legacy
      wrapper preserves that refusal, and the path-only chromosome cache is removed.
- [x] **Step 4:** run `test_chromosome_utils.py`, `test_assembly_guard.py`,
      `test_region_utils.py`, `test_builders.py`. `test_builders.py:263-268` asserts an
      `unknown` result — decide whether that header is *genuinely* ambiguous or was only
      ambiguous because of the bug, and if the latter, change it deliberately and say so
      in the commit message.
- [x] **Step 5:** standard block; commit `Closes #165`.

### Task 9: single-end support (#161)

**Files:** create `vntyper/scripts/read_layout.py`, `vntyper/scripts/kestrel_command.py`;
modify `pipeline.py`, `alignment_processing.py`, `kestrel_genotyping.py`,
`command_builders.py`, `cli_handlers.py`

**Interfaces produced:**
- `classify_layout(r1: int, r2: int, other: int, single: int) -> str` →
  `"paired"` | `"single"` | `"mixed"` | `"empty"`
- `route_fastqs(layout, paths, counts) -> tuple[tuple[str, ...], tuple[str, ...]]` →
  `(consumed, stranded_nonempty)`

- [x] **Step 1: write the failing tests for `read_layout`.**

```python
def test_r1_and_r2_populated_is_paired(self):
    assert classify_layout(r1=900, r2=900, other=0, single=0) == "paired"

def test_everything_in_other_is_single_because_samtools_puts_unpaired_reads_there(self):
    """`samtools fastq -0` collects reads with READ1 and READ2 both unset -- every
    single-end read -- and both pipeline.py call sites bound it to `_` (#161)."""
    assert classify_layout(r1=0, r2=0, other=1800, single=0) == "single"

def test_unequal_r1_and_r2_is_mixed_because_that_is_a_truncated_conversion(self):
    assert classify_layout(r1=900, r2=880, other=0, single=0) == "mixed"

def test_both_populated_is_mixed_and_is_rejected_rather_than_coerced(self):
    assert classify_layout(r1=900, r2=900, other=4200, single=0) == "mixed"

def test_nothing_at_all_is_empty_not_single(self):
    assert classify_layout(r1=0, r2=0, other=0, single=0) == "empty"

def test_a_stranded_non_empty_file_is_reported_never_dropped(self):
    _, stranded = route_fastqs("paired", PATHS, {"r1": 900, "r2": 900, "other": 4200, "single": 0})
    assert stranded == ("other.gz",)
```

- [x] **Step 2:** run → FAIL. **Step 3:** implement `read_layout.py`.
- [x] **Step 4:** extract `construct_kestrel_command` into `vntyper/scripts/kestrel_command.py`
      without behavioural change, re-export it from `kestrel_genotyping` so existing
      imports keep working, and confirm the final module is 865 lines. The extraction is
      the bounded responsibility split required by AGENTS.md rule 3, not a Kestrel-design
      change.
- [x] **Step 5:** write failing tests for the single-input command forms — fastp with no
      `--in2`/`--out2`; bwa with one FASTQ; Kestrel with one FASTQ; Kestrel still raising
      when `fastq_1` is absent. Change the guard from `if not fastq_1 or not fastq_2` to
      `if not fastq_1`, and `fastq2` to `Path | None` in `align_and_sort_fastq`.
- [x] **Step 6:** wire routing into `pipeline.py`: replace both
      `fastq1, fastq2, _, _ = process_bam_to_fastq(...)` bindings with the four-value form
      plus `classify_layout`/`route_fastqs`; raise naming the file and its record count
      when anything is stranded; log the layout.
**Fast-mode scope:** keep the existing explicit `--fast-mode` behaviour: it skips
unmapped-read recovery and retains the historical target slice. This is not a
`read_layout` tolerance or configuration escape hatch; normal-mode conversion must still
route every produced FASTQ or fail naming it.
- [x] **Step 7:** allow `--fastq1` without `--fastq2` in `cli_handlers.py:227`.
- [x] **Step 8:** derive the single-end fixture with pysam, clearing `is_paired`,
      `is_read1`, `is_read2`, `is_proper_pair` and `mate_is_unmapped`; register it in
      `tests/test_data_config.json`; add the integration test for A-161-1.
- [x] **Step 9:** standard block; commit `Closes #161`.

### Task 10: fixture deriver and the golden-cohort matrix

- [x] **Step 1:** failing test — with no `--all`, the deriver selects only the samples
      declared in `tests/test_data_config.json`; with `--all`, every discovered sample.
- [x] **Step 2:** implement the flag and the declared-only default.
- [x] **Step 3:** add `build_reference_dependent_fixture` (its own docstring already
      describes it) producing a reference-compressed CRAM plus a copy of its reference —
      #209's path cannot be tested without one.
- [x] **Step 4:** add a CRAM fixture containing **placed** read-unmapped records, for A-SCAN-1.
- [x] **Step 5:** extend `scripts/golden_cohort/matrix.py` so the CRAM cases cover both
      scan modes; the current matrix declares only two CRAM cases (`matrix.py:126`,
      `:143`), so deriving more fixtures does not by itself run them.
- [x] **Step 6:** commit.

---

## Wave 3 — integration

### Task 11: golden cohort and the performance measurement

- [x] **Step 1:** baseline — three runs on `main`, alternating with the branch, on an idle
      host, through the harness's `run` subcommand. Time the BAM cases that complete on
      both revisions; separately record every case the new no-discard rule refuses so an
      early failure is never compared with completed genotyping. Each declared
      mixed-layout refusal must contain the stable causal routing diagnostic; an unrelated
      exit 1 is not evidence.
- [x] **Step 2:** `make cram-fixtures` (never `--allow-matrix-drift`).
- [x] **Step 3:** three runs on the branch. Report median and range per arm. A regression
      is called only when the slower arm's *best* run is worse than the faster arm's *worst*
      (spec §5b).
- [x] **Step 4:** prove A-178-2 — when preflight authorises both strategies, indexed and
      stream must produce the **same read set**, not merely the same genotype. When
      placed-unmapped evidence or the exact literal-`'*'`/terminal-`*` count comparison
      rejects indexed, preserve that rejection and record the stream read set plus a raw
      diagnostic of the reads `'*'` would lose. The indexed result must also contain the
      exact stable guard diagnostic and its declared causal count (placed-unmapped counts
      are 11,571 for `7a61` and 329 for `b178`); a different preflight exit is not evidence
      for A-178-2:

```bash
python scripts/golden_cohort_gate.py run \
  --side after --tree "$PWD" --data-dir "$PWD/tests/data" \
  --run-root /tmp/m4-final-cram-evidence \
  --marker vntyper.scripts.alignment_contract --expect-marker present \
  --no-probes --no-cohort \
  --case indexed_safe_indexed_cram \
  --case indexed_safe_stream_cram \
  --case 7a61_hg38_ensembl_indexed_cram \
  --case 7a61_hg38_ensembl_stream_cram \
  --case b178_hg19_indexed_cram \
  --case b178_hg19_stream_cram
```

The harness must exit zero. Each result must record
`effective_unmapped_scan == observed_unmapped_scan`. The indexed-safe pair must both
succeed with the same declared read set. Each lossy forced-indexed case must exit 1 with
its exact declared guard count, no observed extraction command and no unmapped BAM; its
stream partner must record the complete expected read set and the declared raw loss.

Plus, per CRAM sample, compare `samtools view -c` on the unmapped BAM and a sorted
read-name digest of it.

- [x] **Step 5:** record both wall-clocks and the read-set result in spec §5b and §10.

### Task 12: full gate, review, PR, release

- [x] **Step 1:** `make check-all`
- [x] **Step 2:** `make patch-coverage` — must be ≥ 80%
- [x] **Step 3:** `make ci-local` **only if** `.github/workflows/` was touched — not applicable; the final diff does not touch workflow files.
- [ ] **Step 4:** Codex adversarial gate on the final **diff** (spec §10, concluding round). Loop
      until no HIGH survives; record the verdict.
- [ ] **Step 5:** `superpowers:requesting-code-review`, then address findings.
- [x] **Step 6:** bump `vntyper/version.py`, `CITATION.cff` and
      `docs/about/changelog.md` to **2.0.10** — patch only, all three files (trap 12).
- [ ] **Step 7:** one PR against `main` whose body carries `Closes #213`, `Closes #225`,
      `Closes #209`, `Closes #178`, `Closes #165`, `Closes #161`.
- [ ] **Step 8:** after merge, release with
      `gh release create v2.0.10 --target main` — **never push a `v*` tag** (it publishes
      to PyPI immediately and irreversibly).

### Final H1/H2/H3 reconciliation: binding, ownership and archival

**Files:** modify `vntyper/scripts/alignment_binding.py`, `alignment_index_binding.py`,
`alignment_consumer_commands.py`, `samtools_command_fragments.py`, `alignment_contract.py`,
`pipeline.py`, `pipeline_cleanup.py`, `pipeline_coverage.py`, `pipeline_inputs.py`,
`cli_logging_safety.py`, `archive_safety.py`, `docker/app/tasks.py`; tests in
`tests/unit/test_alignment_binding_lifecycle.py`, `test_pipeline_coverage.py`,
`test_pipeline_archive_destination_ownership.py`, `test_cli_logging_safety.py`,
`test_archive_safety.py`, `tests/unit/web/test_response_hygiene.py`, and web archive-task tests.

**Interfaces and lifetime:** `AlignmentPlan.binding: AlignmentBinding | None`,
`AlignmentPlan.stable_index_path: str` and
`AlignmentPlan.close() -> None`; `build_alignment_view(in_path, ...) -> (view_path,
index_path, binding)`; `close_alignment_plan(plan, *, preserve_primary: bool) -> None`; and
`calculate_alignment_coverage(plan=..., ...) -> str`. The binding opens the regular input,
installs a procfd view or verified hardlink fallback, opens the fresh temporary index before
publication, and lives until coverage has consumed both exact inodes. `close()` deletes
only exact owned fallback entries before closing the index and alignment FDs.

- [x] **H1:** test path replacement after input open, generated-index replacement after
      preflight, procfs fallbacks, exact owned-entry removal, cleanup failure with and
      without a primary failure, nested `REF_PATH` restoration and overlapping-thread
      serialization; implement both binding lifetimes and primary-outcome preservation.
- [x] **H2:** test exclusive regular-file ownership for outputs and logs, alias rejection,
      patient-tree scope, and `archive_base_name()` trailing-separator normalization;
      implement the shared pipeline/CLI ownership checks before every writer opens a path.
- [x] **H3:** test reports-before-archive, view removal before traversal, descriptor
      anchoring with `O_NOFOLLOW`, symlink/hardlink/special-entry refusal, CLI/web shared
      writer behaviour, descriptor-bound delivery/cohort snapshots, and web
      rollback/quarantine failure handling; implement archive creation only after report
      completion. Final H3 publication and descriptor-snapshot re-reviews are clean.
- [x] **A-PERF-1:** alternate three 18-case runs per arm on clean `ddf49a1` and
      `388f157`; record 90.00/91.13/90.03 s versus 87.83/87.49/87.52 s, exact loads and
      the non-overlapping no-regression verdict in spec §5b.
- [ ] **Close-out:** run and record the combined final-diff adversarial review with no
      HIGH findings.

---

## Self-review

**Spec coverage.** §1 web transport → Task 7 steps 2-4. §4.1a view → Task 5 step 3.
§4.2 reference contract → Tasks 4, 5, 7. §4.3(a) `REF_PATH` → Task 5 step 5.
§4.3(b) scan → Task 2 + Task 6. §4.4 layout → Task 9. §4.5 contracts → Tasks 2, 5.
§5 #213 → Wave 0 (done). §5 #225 → Tasks 1, 3, 5, 6. §5 #209 → Tasks 4, 5, 7.
§5 #178 → Tasks 2, 4, 5, 6, 11. §5 #165 → Task 8. §5 #161 → Task 9.
§5b P1-P3 → Task 4; P4 → Tasks 2, 4, 6. §6 config → Task 5 step 1. §7 → Tasks 11, 12.
§4.1a's binding/lifetime and §4.5 cleanup/`REF_PATH` semantics → Final H1. §4.6 output,
log and archive ownership → Final H2/H3. A-VIEW-1/2 → Final H1; A-OWN-1 and
A-ARCHIVE-1/2 → Final H2/H3; the historical `3083ed9` timing and final-rerun placeholder
→ Final close-out.

**Placeholder scan.** Tasks 5-10 give complete test code for the decisions and prose for
the wiring. That is deliberate and bounded: the wiring's exact shape depends on what Wave 1
landed, and every prose step names the file, the line region and the property to assert, so
there is no step whose target a fresh implementer has to invent.

**Type consistency.** `AlignmentPlan` has no `layout` field in Task 1 and none is read in
Tasks 5, 6 or 9. `reference_path` (a path) is the builder parameter everywhere in Task 4;
No preformatted reference fragment exists in `AlignmentPlan`; command builders receive
`reference_path` and own quoting. `choose_scan` returns `(scan, reason)` in Task 2 and is
unpacked as such in Task 5.
`classify_layout` takes four counts in Task 9's tests and its interface block.
Task 5 uses the singular `resolve_any_index` for general supplied-artifact metadata and the
deliberate BAI-only `resolve_bam_index` for non-fast BAM. Neither result is consumed: every
run builds a fresh BAI or CRAI beside its view.
